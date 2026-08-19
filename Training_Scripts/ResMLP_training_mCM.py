#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D masked reconstruction training script with block-wise softmax.

Notes
-----
- If the HDF5 rows are already unpadded, they are used directly.
- If the HDF5 rows are padded to lmid_max per metabolite-experiment block,
  the generator trims them back to the unpadded representation first.
- Since block-wise softmax guarantees per-block sum = 1, a block-sum error
  metric is included and should stay very close to 0.
- For final reconstruction at inference time, you can still do:
      recon = target * obs_mask + pred * (1.0 - obs_mask)

Important masking behavior
--------------------------
NEXP_KEEP is now the maximum number of experiments to keep.

For example:
    NEXP_KEEP = 3

means each augmented sample randomly keeps:
    1, 2, or 3 experiments

It no longer means always keep exactly 3 experiments.
"""

from __future__ import annotations

# ============================================================
# USER SETTINGS
# ============================================================
NISO = 13

# Maximum number of experiments to keep.
# If NEXP_KEEP = 3, each augmented sample randomly keeps 1, 2, or 3 experiments.
NEXP_KEEP = 3

NANOR = -1.0

H5_PATH = "" #CHANGE THIS
LENGTHS_CSV = "Trained_Models/ResMLP_inpainting_mCM/met_lengths_probaliy_mCM.csv"

Nmask_train = 900
Nmask_val   = 100
BATCH_SIZE  = 2048

# Dense model size
HIDDEN_DIM   = 2048
N_RES_BLOCKS = 4
DROPOUT      = 0.10

# Fine-tuning optimizer settings.
# The pretrained neural-network weights are loaded from CKPT_PATH_old, but a
# weights-only file does not contain Adam's momentum or variance variables.
# Therefore, Adam starts with a new optimizer state at this low learning rate.
LEARNING_RATE = 3e-5
ADAM_CLIPNORM = 1.0

INITIAL_EPOCH = 0
EPOCHS        = 9

MODEL_OUT = f"CCM_mask_{Nmask_train}_pyk_mseHOLE_blockSoftmax_1D_randExp1to{NEXP_KEEP}.keras"
JSON_OUT  = f"CCM_mask_{Nmask_train}_pyk_mseHOLE_blockSoftmax_1D_randExp1to{NEXP_KEEP}.json"
CKPT_PATH = f"CCM_mask_{Nmask_train}_mseHOLE_pyk_blockSoftmax_1D_randExp1to{NEXP_KEEP}.weights.h5"
CKPT_PATH_old = f"CCM_mask_1000_mseHOLE_blockSoftmax_1D_randExp1to3.weights.h5"

# Separate crash-recovery checkpoint directory.
# BackupAndRestore stores the model weights, Adam optimizer variables, and
# training progress here. One recovery checkpoint is overwritten by the next.
BACKUP_DIR = (
    f"CCM_mask_{Nmask_train}_pyk_mseHOLE_blockSoftmax_1D_"
    f"randExp1to{NEXP_KEEP}_training_backup"
)

# Save the recovery state every this many training batches.
# With a very long epoch, batch-level backup avoids losing the entire epoch.
BACKUP_SAVE_FREQ = 5000

# Optional fit-time pipeline parallelism
FIT_WORKERS = 4
FIT_MAX_QUEUE_SIZE = 16

# Reproducibility
GLOBAL_SEED = 0
TRAIN_MASK_SEED = 54321
VAL_MASK_SEED   = 12345

# ============================================================
# Reproducibility and multi-GPU setup
# ============================================================
from numpy.random import seed
seed(GLOBAL_SEED)

import os
import numpy as np
import h5py
import pandas as pd
from typing import List

import tensorflow as tf
tf.random.set_seed(GLOBAL_SEED)
tf.keras.utils.set_random_seed(GLOBAL_SEED)

print("TensorFlow version:", tf.__version__)
print("TensorFlow op determinism is OFF")
print("Fixed validation masking is preserved by VAL_MASK_SEED")

gpus = tf.config.list_physical_devices("GPU")
print("Physical GPUs visible to TF:", gpus)
print("Number of visible GPUs:", len(gpus))

for gpu in gpus:
    try:
        tf.config.experimental.set_memory_growth(gpu, True)
    except Exception as e:
        print(f"Could not set memory growth for {gpu}: {e}")

if len(gpus) > 1:
    strategy = tf.distribute.MirroredStrategy()
else:
    strategy = tf.distribute.get_strategy()

print("Num replicas in sync:", strategy.num_replicas_in_sync)

from tensorflow import keras
import tensorflow.keras.backend as K
from tensorflow.keras.callbacks import (
    BackupAndRestore,
    ModelCheckpoint,
    ReduceLROnPlateau,
)


# ============================================================
# Length loader and padded-vector extractor
# ============================================================
def load_lmid_per_met(format_a_csv_path: str) -> np.ndarray:
    df = pd.read_csv(format_a_csv_path)

    required = {"met_idx", "true_length"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {sorted(required)}")

    df = df.sort_values("met_idx")
    lmid_per_met = df["true_length"].to_numpy(dtype=np.int64)

    if lmid_per_met.min() < 1:
        raise ValueError("true_length must be >= 1")

    return lmid_per_met


def build_keep_idx_from_padded(
    nmet: int,
    niso: int,
    lmid_max: int,
    lmid_per_met: np.ndarray,
) -> np.ndarray:
    """
    Build indices that remove per-block padding and recover the unpadded vector.

    Expected padded layout:
      for each metabolite k:
        for each experiment e:
          block of length lmid_max, of which only first lmid_per_met[k] are real
    """
    lmid_per_met = np.asarray(lmid_per_met, dtype=np.int64)

    keep: List[int] = []
    stride_met = niso * lmid_max

    for k in range(nmet):
        lk = int(lmid_per_met[k])
        base_k = k * stride_met

        for e in range(niso):
            start = base_k + e * lmid_max
            keep.extend(range(start, start + lk))

    return np.asarray(keep, dtype=np.int64)


# ============================================================
# Block metadata for unpadded representation
# ============================================================
def build_block_metadata(lmid_per_met: np.ndarray, niso: int):
    """
    Build metadata for the unpadded vector representation.

    Unpadded vector order:
      metabolite 0, exp 0, exp 1, ..., exp NISO-1,
      metabolite 1, exp 0, exp 1, ..., exp NISO-1,
      ...
    """
    lmid_per_met = np.asarray(lmid_per_met, dtype=np.int64)
    nmet = int(lmid_per_met.size)

    offsets = np.zeros(nmet, dtype=np.int64)
    if nmet > 1:
        offsets[1:] = np.cumsum(niso * lmid_per_met[:-1], dtype=np.int64)

    exp_ids = np.arange(niso, dtype=np.int64)[None, :]

    block_starts = (
        offsets[:, None] + exp_ids * lmid_per_met[:, None]
    ).reshape(-1)

    block_lengths = np.repeat(lmid_per_met, niso)

    nblocks = int(block_starts.size)
    expected_len = int(np.sum(block_lengths))
    max_lmid = int(block_lengths.max())

    local_idx = np.arange(max_lmid, dtype=np.int64)[None, :]
    block_valid_mask = local_idx < block_lengths[:, None]

    # Valid positions in the padded block-logit representation.
    # Flatten order is block-major, then local position.
    valid_padded_idx = np.flatnonzero(
        block_valid_mask.reshape(-1)
    ).astype(np.int32)

    # Segment ids for the unpadded output vector.
    segment_ids = np.repeat(
        np.arange(nblocks, dtype=np.int32),
        block_lengths.astype(np.int32),
    )

    return {
        "block_starts": block_starts.astype(np.int64),
        "block_lengths": block_lengths.astype(np.int64),
        "nblocks": nblocks,
        "expected_len": expected_len,
        "max_lmid": max_lmid,
        "block_valid_mask": block_valid_mask.astype(bool),
        "valid_padded_idx": valid_padded_idx,
        "segment_ids": segment_ids,
    }


# ============================================================
# Custom block-wise softmax layer
# ============================================================
@keras.utils.register_keras_serializable(package="Custom")
class MaskedBlockSoftmax(keras.layers.Layer):
    """
    Converts flat logits of shape (batch, nblocks * max_lmid)
    into block-wise softmax probabilities over valid entries only,
    then gathers only the valid entries back into the unpadded vector
    of shape (batch, expected_len).
    """

    def __init__(
        self,
        nblocks,
        max_lmid,
        block_valid_mask,
        valid_padded_idx,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.nblocks = int(nblocks)
        self.max_lmid = int(max_lmid)
        self.block_valid_mask = np.asarray(block_valid_mask, dtype=bool)
        self.valid_padded_idx = np.asarray(valid_padded_idx, dtype=np.int32)

    def call(self, logits):
        x = tf.reshape(logits, (-1, self.nblocks, self.max_lmid))

        valid_mask = tf.constant(self.block_valid_mask, dtype=tf.bool)
        neg_large = tf.constant(-1e9, dtype=logits.dtype)

        x = tf.where(valid_mask[None, :, :], x, neg_large)
        probs = tf.nn.softmax(x, axis=-1)

        probs_flat = tf.reshape(probs, (-1, self.nblocks * self.max_lmid))
        gather_idx = tf.constant(self.valid_padded_idx, dtype=tf.int32)

        out = tf.gather(probs_flat, gather_idx, axis=1)
        return out

    def get_config(self):
        config = super().get_config()
        config.update({
            "nblocks": self.nblocks,
            "max_lmid": self.max_lmid,
            "block_valid_mask": self.block_valid_mask.tolist(),
            "valid_padded_idx": self.valid_padded_idx.tolist(),
        })
        return config


# ============================================================
# Vectorized generator with deterministic validation masks
#
# y_true payload:
#   channel 0: target vector
#   channel 1: hole mask, 1 only on missing valid pixels
# ============================================================
class CreateAugment1DBlockSoftmax(keras.utils.Sequence):
    def __init__(
        self,
        h5,
        split: str,
        Nmask: int,
        lmid_per_met: np.ndarray,
        block_starts: np.ndarray,
        block_lengths: np.ndarray,
        batch_size: int = 512,
        shuffle: bool = True,
        nanor: float = -1.0,
        niso: int = 13,
        nexp_keep: int = 3,
        train_mask_seed: int = 54321,
        val_mask_seed: int = 12345,
    ):
        self.h5 = h5
        self.split = split
        self.Nmask = int(Nmask)
        self.batch_size = int(batch_size)
        self.shuffle = bool(shuffle)

        self.nanor = float(nanor)
        self.niso = int(niso)

        # This is now the maximum number of experiments to keep.
        # The actual number kept per sample is randomly chosen from:
        # 1, 2, ..., self.nexp_keep
        self.nexp_keep = int(nexp_keep)

        if self.nexp_keep < 1 or self.nexp_keep > self.niso:
            raise ValueError(
                f"nexp_keep must satisfy 1 <= nexp_keep <= niso, "
                f"got {self.nexp_keep} vs {self.niso}"
            )

        self.train_mask_seed = int(train_mask_seed)
        self.val_mask_seed = int(val_mask_seed)

        self.is_validation = split != "training"
        self.base_seed = self.val_mask_seed if self.is_validation else self.train_mask_seed

        self.lmid_per_met = np.asarray(lmid_per_met, dtype=np.int64)
        self.nmet = int(self.lmid_per_met.size)

        self.block_starts = np.asarray(block_starts, dtype=np.int64)
        self.block_lengths = np.asarray(block_lengths, dtype=np.int64)

        self.nblocks = int(self.block_starts.size)
        self.expected_len = int(np.sum(self.block_lengths))
        self.max_lmid = int(self.block_lengths.max())

        x0 = np.asarray(self.h5[self.split][0]).reshape(-1)
        self.raw_len = int(x0.size)

        self.keep_idx = None

        if self.raw_len == self.expected_len:
            print(f"[{split}] Detected UNPADDED vectors, length={self.raw_len}")
        else:
            denom = self.nmet * self.niso

            if self.raw_len % denom != 0:
                raise ValueError(
                    f"[{split}] raw_len={self.raw_len} not divisible by "
                    f"nmet*niso={denom}, cannot infer lmid_max"
                )

            lmid_max = self.raw_len // denom

            print(
                f"[{split}] Detected PADDED vectors, "
                f"raw_len={self.raw_len}, inferred lmid_max={lmid_max}"
            )

            self.keep_idx = build_keep_idx_from_padded(
                self.nmet,
                self.niso,
                lmid_max,
                self.lmid_per_met,
            )

            if self.keep_idx.size != self.expected_len:
                raise RuntimeError(
                    f"[{split}] keep_idx size {self.keep_idx.size} "
                    f"!= expected_len {self.expected_len}"
                )

        local_idx = np.arange(self.max_lmid, dtype=np.int64)[None, :]
        self.block_idx_mat = self.block_starts[:, None] + local_idx
        self.block_valid_mat = local_idx < self.block_lengths[:, None]

        self.on_epoch_end()

    def __len__(self):
        return int(
            np.ceil(
                len(self.h5[self.split]) * self.Nmask / float(self.batch_size)
            )
        )

    def on_epoch_end(self):
        self.indexes = np.arange(len(self.h5[self.split]) * self.Nmask)

        if self.shuffle:
            np.random.shuffle(self.indexes)

    def __getitem__(self, idx):
        indexes = self.indexes[
            idx * self.batch_size:(idx + 1) * self.batch_size
        ]

        return self.__data_generation(indexes)

    def _fetch_batch_unpadded(self, base_idxs: np.ndarray) -> np.ndarray:
        base_idxs = np.asarray(base_idxs, dtype=np.int64)

        uniq, inv = np.unique(base_idxs, return_inverse=True)

        raw = np.asarray(self.h5[self.split][uniq], dtype=np.float32)
        raw = raw.reshape(raw.shape[0], -1)
        raw = raw[inv]

        if self.keep_idx is not None:
            raw = raw[:, self.keep_idx]

        if raw.shape[1] != self.expected_len:
            raise ValueError(
                f"Expected {self.expected_len} columns, got {raw.shape[1]}"
            )

        return raw

    def _make_deterministic_rng_seeds(
        self,
        batch_augmented_indices: np.ndarray,
    ) -> np.ndarray:
        batch_augmented_indices = np.asarray(
            batch_augmented_indices,
            dtype=np.int64,
        )

        return (
            self.base_seed + batch_augmented_indices * 1000003
        ).astype(np.int64)

    def __data_generation(self, idxs):
        bs = len(idxs)
        idxs = np.asarray(idxs, dtype=np.int64)

        # ----------------------------------------------------
        # 1) Vectorized HDF5 read, converted to unpadded 1D
        # ----------------------------------------------------
        base_idxs = idxs // self.Nmask
        x_batch = self._fetch_batch_unpadded(base_idxs)

        # ----------------------------------------------------
        # 2) Generate masks
        #
        # Validation:
        #   deterministic per augmented index
        #
        # Training:
        #   random every time
        #
        # Experiment selection:
        #   randomly keeps 1, 2, ..., self.nexp_keep experiments
        # ----------------------------------------------------
        if self.is_validation:
            sample_seeds = self._make_deterministic_rng_seeds(idxs)

            rnd_exp = np.empty((bs, self.niso), dtype=np.float32)
            keepmet_mask = np.empty((bs, self.nmet), dtype=bool)
            forced_cols = np.empty((bs,), dtype=np.int64)
            nexp_keep_each = np.empty((bs,), dtype=np.int64)

            for i in range(bs):
                rng = np.random.default_rng(int(sample_seeds[i]))

                rnd_exp[i] = rng.random(self.niso, dtype=np.float32)
                keepmet_mask[i] = rng.random(self.nmet) >= 0.5
                forced_cols[i] = rng.integers(0, self.nmet)

                # Random number of experiments to keep:
                # 1, 2, ..., self.nexp_keep
                nexp_keep_each[i] = rng.integers(
                    1,
                    self.nexp_keep + 1,
                )

        else:
            rnd_exp = np.random.rand(bs, self.niso).astype(np.float32)
            keepmet_mask = np.random.rand(bs, self.nmet) >= 0.5
            forced_cols = np.random.randint(0, self.nmet, size=bs)

            # Random number of experiments to keep:
            # 1, 2, ..., self.nexp_keep
            nexp_keep_each = np.random.randint(
                1,
                self.nexp_keep + 1,
                size=bs,
            )

        # Rank experiments using random scores.
        # Smaller values are selected first.
        exp_order = np.argsort(rnd_exp, axis=1)

        keepexp_mask = np.zeros((bs, self.niso), dtype=bool)

        for i in range(bs):
            n_keep_i = int(nexp_keep_each[i])
            selected_exps_i = exp_order[i, :n_keep_i]
            keepexp_mask[i, selected_exps_i] = True

        # Guarantee that every sample keeps at least one metabolite.
        empty_rows = ~keepmet_mask.any(axis=1)

        if np.any(empty_rows):
            keepmet_mask[empty_rows, forced_cols[empty_rows]] = True

        # ----------------------------------------------------
        # 3) Select kept values blockwise
        #
        # block_keep shape:
        #   (batch, nmet, niso) reshaped to (batch, nblocks)
        #
        # A block is observed only if:
        #   metabolite is kept AND experiment is kept
        # ----------------------------------------------------
        block_keep = (
            keepmet_mask[:, :, None] & keepexp_mask[:, None, :]
        ).reshape(bs, self.nblocks)

        selected = block_keep[:, :, None] & self.block_valid_mat[None, :, :]

        row_idx, block_idx, local_idx = np.nonzero(selected)
        vec_idx = self.block_idx_mat[block_idx, local_idx]

        # ----------------------------------------------------
        # 4) Build masked inputs and observation mask
        # ----------------------------------------------------
        masked_input = np.full(
            (bs, self.expected_len),
            self.nanor,
            dtype=np.float32,
        )

        obs_mask = np.zeros(
            (bs, self.expected_len),
            dtype=np.float32,
        )

        masked_input[row_idx, vec_idx] = x_batch[row_idx, vec_idx]
        obs_mask[row_idx, vec_idx] = 1.0

        # ----------------------------------------------------
        # 5) Hole mask, only missing valid pixels
        # ----------------------------------------------------
        hole_mask = 1.0 - obs_mask

        # ----------------------------------------------------
        # 6) y payload
        #
        # channel 0 = target
        # channel 1 = hole mask
        # ----------------------------------------------------
        y_payload = np.concatenate(
            [
                x_batch[..., None],
                hole_mask[..., None],
            ],
            axis=-1,
        ).astype(np.float32)

        return [masked_input, obs_mask], y_payload


# ============================================================
# Build metadata
# ============================================================
lmid_per_met = load_lmid_per_met(LENGTHS_CSV)
meta = build_block_metadata(lmid_per_met, niso=NISO)

block_starts = meta["block_starts"]
block_lengths = meta["block_lengths"]
nblocks = meta["nblocks"]
expected_len = meta["expected_len"]
max_lmid = meta["max_lmid"]
block_valid_mask = meta["block_valid_mask"]
valid_padded_idx = meta["valid_padded_idx"]
segment_ids = meta["segment_ids"]

print(f"nmet={len(lmid_per_met)}")
print(f"nblocks={nblocks}")
print(f"expected unpadded length={expected_len}")
print(f"max_lmid={max_lmid}")
print(f"padded block-logit length={nblocks * max_lmid}")

SEGMENT_IDS_TF = tf.constant(segment_ids, dtype=tf.int32)


# ============================================================
# Data and generators
# ============================================================
print("Opening HDF5:", H5_PATH)
h5 = h5py.File(H5_PATH, "r")

traingen = CreateAugment1DBlockSoftmax(
    h5=h5,
    split="training",
    Nmask=Nmask_train,
    lmid_per_met=lmid_per_met,
    block_starts=block_starts,
    block_lengths=block_lengths,
    batch_size=BATCH_SIZE,
    shuffle=True,
    nanor=NANOR,
    niso=NISO,
    nexp_keep=NEXP_KEEP,
    train_mask_seed=TRAIN_MASK_SEED,
    val_mask_seed=VAL_MASK_SEED,
)

valgen = CreateAugment1DBlockSoftmax(
    h5=h5,
    split="validation",
    Nmask=Nmask_val,
    lmid_per_met=lmid_per_met,
    block_starts=block_starts,
    block_lengths=block_lengths,
    batch_size=BATCH_SIZE,
    shuffle=False,
    nanor=NANOR,
    niso=NISO,
    nexp_keep=NEXP_KEEP,
    train_mask_seed=TRAIN_MASK_SEED,
    val_mask_seed=VAL_MASK_SEED,
)

print(
    f"raw_len(train)={traingen.raw_len}, expected_len={traingen.expected_len}, "
    f"NISO={NISO}, NEXP_KEEP_MAX={NEXP_KEEP}, VAL_MASK_SEED={VAL_MASK_SEED}"
)


# ============================================================
# Losses and metrics
#
# y_true payload:
#   y_true[..., 0:1] = target
#   y_true[..., 1:2] = hole mask
# ============================================================
def split_target_and_hole(y_true):
    target = y_true[..., 0:1]
    hole_mask = y_true[..., 1:2]
    return target, hole_mask


@keras.utils.register_keras_serializable(package="Custom")
def masked_mse_loss(y_true, y_pred):
    target, hole_mask = split_target_and_hole(y_true)

    sq_err = tf.square(y_pred - target) * hole_mask

    denom = tf.reduce_sum(hole_mask, axis=[1, 2]) + K.epsilon()
    numer = tf.reduce_sum(sq_err, axis=[1, 2])

    return numer / denom


@keras.utils.register_keras_serializable(package="Custom")
def masked_mae_metric(y_true, y_pred):
    target, hole_mask = split_target_and_hole(y_true)

    abs_err = tf.abs(y_pred - target) * hole_mask

    denom = tf.reduce_sum(hole_mask, axis=[1, 2]) + K.epsilon()
    numer = tf.reduce_sum(abs_err, axis=[1, 2])

    return numer / denom


@keras.utils.register_keras_serializable(package="Custom")
def hole_fraction_metric(y_true, y_pred):
    _, hole_mask = split_target_and_hole(y_true)
    return tf.reduce_mean(hole_mask)


@keras.utils.register_keras_serializable(package="Custom")
def block_sum_abs_err_metric(y_true, y_pred):
    """
    Mean absolute deviation of each predicted block sum from 1.
    With correct block-softmax this should be very close to 0.
    """
    pred = tf.squeeze(y_pred, axis=-1)

    pred_t = tf.transpose(pred, perm=[1, 0])

    block_sums_t = tf.math.unsorted_segment_sum(
        pred_t,
        SEGMENT_IDS_TF,
        nblocks,
    )

    block_sums = tf.transpose(block_sums_t, perm=[1, 0])

    return tf.reduce_mean(tf.abs(block_sums - 1.0))


# ============================================================
# Model definition
# ============================================================
class BlockSoftmax1DModel:
    def __init__(
        self,
        input_dim: int,
        nblocks: int,
        max_lmid: int,
        block_valid_mask: np.ndarray,
        valid_padded_idx: np.ndarray,
        hidden_dim: int = 2048,
        n_res_blocks: int = 4,
        dropout: float = 0.10,
    ):
        self.input_dim = int(input_dim)
        self.nblocks = int(nblocks)
        self.max_lmid = int(max_lmid)

        self.block_valid_mask = np.asarray(block_valid_mask, dtype=bool)
        self.valid_padded_idx = np.asarray(valid_padded_idx, dtype=np.int32)

        self.hidden_dim = int(hidden_dim)
        self.n_res_blocks = int(n_res_blocks)
        self.dropout = float(dropout)

    def residual_block(self, x, hidden_dim, dropout, name_prefix):
        skip = x

        y = keras.layers.Dense(
            hidden_dim * 2,
            name=f"{name_prefix}_dense1",
        )(x)

        y = keras.layers.Activation(
            "gelu",
            name=f"{name_prefix}_gelu1",
        )(y)

        y = keras.layers.Dropout(
            dropout,
            name=f"{name_prefix}_drop1",
        )(y)

        y = keras.layers.Dense(
            hidden_dim,
            name=f"{name_prefix}_dense2",
        )(y)

        y = keras.layers.Dropout(
            dropout,
            name=f"{name_prefix}_drop2",
        )(y)

        x = keras.layers.Add(name=f"{name_prefix}_add")([skip, y])

        x = keras.layers.LayerNormalization(
            name=f"{name_prefix}_ln",
        )(x)

        return x

    def prepare_model(self):
        masked_input = keras.layers.Input(
            shape=(self.input_dim,),
            name="masked_input",
        )

        obs_mask = keras.layers.Input(
            shape=(self.input_dim,),
            name="obs_mask",
        )

        # Missing entries become 0 in the value channel.
        # obs_mask carries missingness explicitly.
        masked_values_only = keras.layers.Multiply(
            name="apply_obs_mask",
        )([masked_input, obs_mask])

        x = keras.layers.Concatenate(
            name="concat_input",
        )([masked_values_only, obs_mask])

        x = keras.layers.LayerNormalization(name="input_ln")(x)

        x = keras.layers.Dense(
            self.hidden_dim,
            name="stem_dense",
        )(x)

        x = keras.layers.Activation(
            "gelu",
            name="stem_gelu",
        )(x)

        x = keras.layers.Dropout(
            self.dropout,
            name="stem_drop",
        )(x)

        x = keras.layers.LayerNormalization(name="stem_ln")(x)

        for i in range(self.n_res_blocks):
            x = self.residual_block(
                x,
                hidden_dim=self.hidden_dim,
                dropout=self.dropout,
                name_prefix=f"resblk{i + 1}",
            )

        logits = keras.layers.Dense(
            self.nblocks * self.max_lmid,
            name="block_logits",
        )(x)

        outputs = MaskedBlockSoftmax(
            nblocks=self.nblocks,
            max_lmid=self.max_lmid,
            block_valid_mask=self.block_valid_mask,
            valid_padded_idx=self.valid_padded_idx,
            name="block_softmax_output",
        )(logits)

        outputs = keras.layers.Reshape(
            (self.input_dim, 1),
            name="output_reshape",
        )(outputs)

        return keras.models.Model(
            inputs=[masked_input, obs_mask],
            outputs=outputs,
            name="BlockSoftmax1DImputer",
        )


# ============================================================
# Build, compile, load, and train
# ============================================================
def directory_has_contents(path: str) -> bool:
    """Return True when a recovery directory exists and is not empty."""
    if not os.path.isdir(path):
        return False

    try:
        with os.scandir(path) as entries:
            return any(True for _ in entries)
    except OSError as exc:
        print(f"Could not inspect backup directory {path}: {exc}")
        return False


# Check before loading pretrained weights. If a recovery backup exists,
# BackupAndRestore will restore it when model.fit() starts. Loading the old
# pretrained weights first would be unnecessary and potentially confusing.
backup_exists = directory_has_contents(BACKUP_DIR)

with strategy.scope():
    model = BlockSoftmax1DModel(
        input_dim=expected_len,
        nblocks=nblocks,
        max_lmid=max_lmid,
        block_valid_mask=block_valid_mask,
        valid_padded_idx=valid_padded_idx,
        hidden_dim=HIDDEN_DIM,
        n_res_blocks=N_RES_BLOCKS,
        dropout=DROPOUT,
    ).prepare_model()

    model.compile(
        optimizer=keras.optimizers.Adam(
            learning_rate=LEARNING_RATE,
            clipnorm=ADAM_CLIPNORM,
        ),
        loss=masked_mse_loss,
        metrics=[
            masked_mae_metric,
            hole_fraction_metric,
            block_sum_abs_err_metric,
        ],
    )

model.summary()

print(f"Configured Adam learning rate: {LEARNING_RATE:.2e}")
print(f"Configured Adam clipnorm: {ADAM_CLIPNORM}")

if backup_exists:
    print("=" * 72)
    print("INTERRUPTED-TRAINING BACKUP FOUND")
    print(f"Recovery directory: {BACKUP_DIR}")
    print(
        "When model.fit() starts, BackupAndRestore will restore the model "
        "weights, Adam optimizer state, and saved training progress."
    )
    print("The original pretrained weights will not be loaded separately.")
    print("=" * 72)

elif os.path.exists(CKPT_PATH_old):
    try:
        model.load_weights(CKPT_PATH_old)
        print(f"Loaded pretrained weights from: {CKPT_PATH_old}")
        print(
            "Adam has a new optimizer state and will fine-tune with "
            f"learning_rate={LEARNING_RATE:.2e}."
        )
    except Exception as e:
        raise RuntimeError(
            f"Could not load pretrained checkpoint {CKPT_PATH_old}. "
            "Training was stopped to avoid running 900 masks per sample from "
            "random initialization with the fine-tuning learning rate."
        ) from e
else:
    raise FileNotFoundError(
        f"No recovery backup was found in {BACKUP_DIR}, and the pretrained "
        f"weights file was not found: {CKPT_PATH_old}"
    )

with open(JSON_OUT, "w") as f:
    f.write(model.to_json())

print(f"Saved model architecture JSON to: {JSON_OUT}")

# Best-validation weights only. This file does not contain Adam state.
best_weights_checkpoint = ModelCheckpoint(
    filepath=CKPT_PATH,
    monitor="val_loss",
    verbose=1,
    save_best_only=True,
    save_weights_only=True,
    mode="min",
)

# Crash recovery. This stores the model, Adam optimizer variables, epoch, and
# batch progress in BACKUP_DIR. The directory remains after an interruption and
# is removed automatically after model.fit() completes successfully.
backup_restore = BackupAndRestore(
    backup_dir=BACKUP_DIR,
    save_freq=BACKUP_SAVE_FREQ,
    delete_checkpoint=True,
)

# Reduce the already-small fine-tuning learning rate when validation loss stops
# improving. The optimizer's current learning rate is included in the recovery
# checkpoint because it is part of the optimizer state.
reduce_lr = ReduceLROnPlateau(
    monitor="val_loss",
    factor=0.5,
    patience=1,
    min_delta=1e-7,
    min_lr=1e-7,
    verbose=1,
    mode="min",
)

callbacks = [
    backup_restore,
    best_weights_checkpoint,
    reduce_lr,
]

print(f"Best-validation weights path: {CKPT_PATH}")
print(f"Crash-recovery directory: {BACKUP_DIR}")
print(
    f"Model and Adam recovery state will be saved every "
    f"{BACKUP_SAVE_FREQ:,} training batches."
)
print(
    "If training is interrupted, rerun this same script without deleting "
    "the recovery directory."
)

try:
    history = model.fit(
        traingen,
        validation_data=valgen,
        initial_epoch=INITIAL_EPOCH,
        epochs=EPOCHS,
        workers=FIT_WORKERS,
        max_queue_size=FIT_MAX_QUEUE_SIZE,
        use_multiprocessing=False,
        callbacks=callbacks,
    )

    # Save the full model using the weights that achieved the best validation
    # loss, rather than automatically saving the final epoch's weights.
    if os.path.exists(CKPT_PATH):
        model.load_weights(CKPT_PATH)
        print(f"Reloaded best validation weights from: {CKPT_PATH}")
    else:
        print(
            f"WARNING: Best-validation checkpoint was not found at "
            f"{CKPT_PATH}. Saving the current model weights instead."
        )

    model.save(MODEL_OUT)

    if os.path.exists(MODEL_OUT):
        print(f"Saved full model to: {MODEL_OUT}")
    else:
        print(f"WARNING: save did not create file at: {MODEL_OUT}")

finally:
    h5.close()
    print("Closed HDF5 file.")

