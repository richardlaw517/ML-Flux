#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

# ============================================================
# REPRODUCIBILITY
# ============================================================
from numpy.random import seed as np_seed
np_seed(0)

import tensorflow as tf
tf.random.set_seed(0)
tf.keras.utils.set_random_seed(0)

# ============================================================
# IMPORTS
# ============================================================
import os
from typing import Optional

import numpy as np
import pandas as pd
from tensorflow import keras
from keras.models import model_from_json



# ============================================================
# CSV / LENGTH HELPERS
# ============================================================
def load_lmid_per_met(csv_path: str) -> np.ndarray:
    df = pd.read_csv(csv_path)
    required = {"met_idx", "true_length"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {sorted(required)}")

    df = df.sort_values("met_idx").reset_index(drop=True)
    lmid_per_met = df["true_length"].to_numpy(dtype=np.int64)

    if lmid_per_met.min() < 1:
        raise ValueError("true_length values must be >= 1")

    return lmid_per_met


def build_keep_idx(nmet: int, niso: int, lmid_max: int, lmid_per_met: np.ndarray) -> np.ndarray:
    lmid_per_met = np.asarray(lmid_per_met, dtype=np.int64)

    if lmid_per_met.shape[0] != nmet:
        raise ValueError(f"lmid_per_met must have length nmet={nmet}")

    if lmid_per_met.max() > lmid_max:
        raise ValueError(f"Found true_length={lmid_per_met.max()} > lmid_max={lmid_max}")

    keep = []
    stride_met = niso * lmid_max

    for k in range(nmet):
        lk = int(lmid_per_met[k])
        base_k = k * stride_met
        for e in range(niso):
            start = base_k + e * lmid_max
            keep.extend(range(start, start + lk))

    return np.asarray(keep, dtype=np.int64)


def infer_storage_layout_from_raw_len(raw_len: int, lmid_per_met: np.ndarray, niso: int):
    expected_len = int(niso * np.sum(lmid_per_met))
    nmet = int(lmid_per_met.size)

    if raw_len == expected_len:
        keep_idx = None
        layout_name = "UNPADDED"
    else:
        denom = nmet * niso
        if raw_len % denom != 0:
            raise ValueError(
                f"Cannot infer padded layout from raw_len={raw_len}. "
                f"expected_len={expected_len}, nmet={nmet}, niso={niso}"
            )
        lmid_max = raw_len // denom
        keep_idx = build_keep_idx(nmet, niso, lmid_max, lmid_per_met)
        layout_name = f"PADDED(lmid_max={lmid_max})"

    return expected_len, keep_idx, layout_name


def to_unpadded_batch(x_raw: np.ndarray, keep_idx: Optional[np.ndarray]) -> np.ndarray:
    x_raw = np.asarray(x_raw, dtype=np.float32)
    if keep_idx is None:
        return x_raw.astype(np.float32, copy=False)
    return x_raw[:, keep_idx].astype(np.float32, copy=False)


def repad_to_fixed_length(
    unpadded_batch: np.ndarray,
    keep_idx_padded: np.ndarray,
    padded_len: int
) -> np.ndarray:
    unpadded_batch = np.asarray(unpadded_batch, dtype=np.float32)
    out = np.zeros((unpadded_batch.shape[0], padded_len), dtype=np.float32)
    out[:, keep_idx_padded] = unpadded_batch
    return out


# ============================================================
# BLOCK METADATA FOR CUSTOM SOFTMAX LAYER
# ============================================================
def build_block_metadata(lmid_per_met: np.ndarray, niso: int):
    lmid_per_met = np.asarray(lmid_per_met, dtype=np.int64)
    nmet = int(lmid_per_met.size)

    offsets = np.zeros(nmet, dtype=np.int64)
    if nmet > 1:
        offsets[1:] = np.cumsum(niso * lmid_per_met[:-1], dtype=np.int64)

    exp_ids = np.arange(niso, dtype=np.int64)[None, :]
    block_starts = (offsets[:, None] + exp_ids * lmid_per_met[:, None]).reshape(-1)
    block_lengths = np.repeat(lmid_per_met, niso)

    nblocks = int(block_starts.size)
    expected_len = int(np.sum(block_lengths))
    max_lmid = int(block_lengths.max())

    local_idx = np.arange(max_lmid, dtype=np.int64)[None, :]
    block_valid_mask = local_idx < block_lengths[:, None]

    valid_padded_idx = np.flatnonzero(block_valid_mask.reshape(-1)).astype(np.int32)
    segment_ids = np.repeat(np.arange(nblocks, dtype=np.int32), block_lengths.astype(np.int32))

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
# CUSTOM BLOCK-WISE SOFTMAX LAYER
# ============================================================
@keras.utils.register_keras_serializable(package="Custom")
class MaskedBlockSoftmax(keras.layers.Layer):
    def __init__(self, nblocks, max_lmid, block_valid_mask, valid_padded_idx, **kwargs):
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
        config.update(
            {
                "nblocks": self.nblocks,
                "max_lmid": self.max_lmid,
                "block_valid_mask": self.block_valid_mask.tolist(),
                "valid_padded_idx": self.valid_padded_idx.tolist(),
            }
        )
        return config


# ============================================================
# MODEL LOADER
# ============================================================
def load_blocksoftmax_1d_model(json_path: str, weights_path: str):
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Could not find json file: {json_path}")
    if not os.path.exists(weights_path):
        raise FileNotFoundError(f"Could not find weights file: {weights_path}")

    with open(json_path, "r", encoding="utf-8") as f:
        model_json = f.read()

    model = keras.models.model_from_json(
        model_json,
        custom_objects={"MaskedBlockSoftmax": MaskedBlockSoftmax},
    )
    model.load_weights(weights_path)


    print(f"Loaded model JSON   : {json_path}")
    print(f"Loaded model weights: {weights_path}")
    return model


# ============================================================
# MAIN INFERENCE FUNCTION
# ============================================================
def inpaintLabelsFromScratched_mCM_ANN(
    scratchedLabelSet=None,
    modelType: str = "mCM",
    enforce_observed_pixels: bool = True,
    lengths_csv_path: str = "Trained_Models/ResMLP_Label_Inpainting/met_lengths_probality_mCM.csv",
    save_debug_files: bool = True,
    nanor: float = -1.0,
    niso: int = 13,
    pred_batch_size: int = 2048,
    return_padded_10608: bool = True,
    model_json_path=r"Trained_Models/ResMLP_Label_Inpainting/ResMLP_inpainting_mCM.json",
    model_weights_path=r"Trained_Models/ResMLP_Label_Inpainting/ResMLP_inpainting_mCM.h5"
  ):
    """
    Notes
    -----
    For modelType == "MammalianCCM_Packed80" in this updated version:
    - 'Packed80' is treated as a compatibility alias for the new 1D model
    - model input is [masked_input, obs_mask]
    - no gamma-related processing is used
    - no inference-time clipping or normalization is applied
    """

    if modelType != "mCM":
        raise ValueError(
            "This updated inference file only implements the 1D block-softmax path "
            "through modelType='MammalianCCM_Packed80'."
        )

    if scratchedLabelSet is None:
        raise ValueError("scratchedLabelSet must be provided.")

    # --------------------------------------------------------
    # Metadata from CSV
    # --------------------------------------------------------
    lmid_per_met = load_lmid_per_met(lengths_csv_path)
    nmet = int(lmid_per_met.size)
    meta = build_block_metadata(lmid_per_met, niso=niso)
    expected_len = int(meta["expected_len"])

    # This is for compatibility with your padded .data workflow
    lmid_max_output = 12
    padded_len_10608 = int(nmet * niso * lmid_max_output)

    keep_idx_10608 = build_keep_idx(
        nmet=nmet,
        niso=niso,
        lmid_max=lmid_max_output,
        lmid_per_met=lmid_per_met,
    )

    if keep_idx_10608.size != expected_len:
        raise RuntimeError(
            f"keep_idx_10608 size {keep_idx_10608.size} does not match expected_len {expected_len}"
        )

    # --------------------------------------------------------
    # Load model
    # --------------------------------------------------------
    model = load_blocksoftmax_1d_model(model_json_path, model_weights_path)

    print("=" * 80)
    print("Using 1D block-softmax inference path through modelType='MammalianCCM_Packed80'")
    print(f"lengths_csv_path         : {lengths_csv_path}")
    print(f"model_json_path          : {model_json_path}")
    print(f"model_weights_path       : {model_weights_path}")
    print(f"nmet                     : {nmet}")
    print(f"niso                     : {niso}")
    print(f"expected_len             : {expected_len}")
    print(f"return_padded_10608      : {return_padded_10608}")
    print(f"enforce_observed_pixels  : {enforce_observed_pixels}")
    print(f"nanor                    : {nanor}")
    print("=" * 80)

    # --------------------------------------------------------
    # Reshape input
    # --------------------------------------------------------
    scratchedLabelSet = np.asarray(scratchedLabelSet, dtype=np.float32)
    if scratchedLabelSet.ndim == 1:
        scratchedLabelSet = scratchedLabelSet[None, :]
    scratchedLabelSet = scratchedLabelSet.reshape(scratchedLabelSet.shape[0], -1)

    raw_len = int(scratchedLabelSet.shape[1])

    expected_len_check, keep_idx_in, layout_name = infer_storage_layout_from_raw_len(
        raw_len=raw_len,
        lmid_per_met=lmid_per_met,
        niso=niso,
    )
    if expected_len_check != expected_len:
        raise RuntimeError(
            f"Expected length mismatch: {expected_len_check} vs {expected_len}"
        )

    print("scratchedLabelSet reshaped:", scratchedLabelSet.shape)
    print("Input vector layout       :", layout_name)

    # --------------------------------------------------------
    # Convert to unpadded vector if needed
    # --------------------------------------------------------
    x_unpadded = to_unpadded_batch(scratchedLabelSet, keep_idx_in)

    if x_unpadded.shape[1] != expected_len:
        raise ValueError(
            f"After unpadding, expected {expected_len} columns but got {x_unpadded.shape[1]}"
        )

    obs_mask = (x_unpadded != nanor).astype(np.float32)
    true_raw_for_blend = np.where(obs_mask > 0.5, x_unpadded, 0.0).astype(np.float32)

    masked_input = x_unpadded.astype(np.float32, copy=False)

    # --------------------------------------------------------
    # Predict
    # --------------------------------------------------------
    y_pred = model.predict(
        [masked_input, obs_mask],
        batch_size=pred_batch_size,
        verbose=1,
    )

    y_pred = np.asarray(y_pred, dtype=np.float32)

    if y_pred.ndim == 3 and y_pred.shape[-1] == 1:
        pred_vec = np.squeeze(y_pred, axis=-1)
    elif y_pred.ndim == 2:
        pred_vec = y_pred
    else:
        raise ValueError(f"Unexpected model output shape: {y_pred.shape}")

    if pred_vec.shape[1] != expected_len:
        raise ValueError(
            f"Predicted vector has length {pred_vec.shape[1]}, expected {expected_len}"
        )

    pred_vec = pred_vec.astype(np.float32, copy=False)

    # --------------------------------------------------------
    # Blend back observed pixels
    # --------------------------------------------------------
    if enforce_observed_pixels:
        pred_vec = (
            true_raw_for_blend + (1.0 - obs_mask) * pred_vec
        ).astype(np.float32, copy=False)

    # --------------------------------------------------------
    # Diagnostics
    # --------------------------------------------------------
    unpadded_sums = np.sum(pred_vec, axis=1)
    print("Per-sample predicted sums (unpadded):")
    print(unpadded_sums)

    print(
        "Unpadded sum stats: min/mean/max =",
        float(unpadded_sums.min()),
        float(unpadded_sums.mean()),
        float(unpadded_sums.max()),
    )

    if save_debug_files:
        np.savetxt("pred_unpadded_raw.txt", pred_vec, fmt="%.9g")
        np.savetxt("pred_sums_unpadded_raw.txt", unpadded_sums, fmt="%.8g")

    # --------------------------------------------------------
    # Return padded 10608 if requested
    # --------------------------------------------------------
    if return_padded_10608:
        pred_padded = repad_to_fixed_length(
            unpadded_batch=pred_vec,
            keep_idx_padded=keep_idx_10608,
            padded_len=padded_len_10608,
        )

        padded_sums = np.sum(pred_padded, axis=1)
        print("Per-sample predicted sums (padded 10608):")
        print(padded_sums)

        print(
            "Padded 10608 sum stats: min/mean/max =",
            float(padded_sums.min()),
            float(padded_sums.mean()),
            float(padded_sums.max()),
        )

        if save_debug_files:
            np.savetxt("pred_padded_10608.txt", pred_padded, fmt="%.9g")
            np.savetxt("pred_sums_padded_10608.txt", padded_sums, fmt="%.8g")

        return pred_padded

    return pred_vec
