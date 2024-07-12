# Fix random seed number for reproducibility
from numpy.random import seed
from joblib import Parallel, delayed
seed(0)
import tensorflow as tf
tf.random.set_seed(0)

import tensorflow as tf
from tensorflow import keras

import os
import numpy as np
from keras.callbacks import ModelCheckpoint
from sklearn.model_selection import train_test_split
from keras.models import model_from_json

from pconv_layer import PConv2D

## For more information into formulation: https://www.youtube.com/watch?v=AZr64OxshLo
## Metric

def dice_coef(y_true, y_pred):
    y_true_f = keras.backend.flatten(y_true)
    y_pred_f = keras.backend.flatten(y_pred)
    intersection = keras.backend.sum(y_true_f * y_pred_f)
    return (2. * intersection) / (keras.backend.sum(y_true_f + y_pred_f))

# class InpaintingModel:
#   '''
#   Build UNET like model for image inpaining task.
#   '''
#   def prepare_model(self, input_size=(80,80,1)):
#     input_image = keras.layers.Input(input_size)
#     input_mask = keras.layers.Input(input_size, name='encoder_input')

#     conv1, mask1, conv2, mask2 = self.__encoder_layer(80, input_image, input_mask, ['conv1', 'conv2'])
#     conv3, mask3, conv4, mask4 = self.__encoder_layer(160, conv2, mask2, ['conv3', 'conv4'])
#     conv5, mask5, conv6, mask6 = self.__encoder_layer(320, conv4, mask4, ['conv5', 'conv6'])
#     conv7, mask7, conv8, mask8 = self.__encoder_layer(640, conv6, mask6, ['conv7', 'encoder_output'])

#     conv9, mask9, conv10, mask10 = self.__decoder_layer(640, 320, conv8, mask8, conv7, mask7, ['conv9', 'conv10'])
#     conv11, mask11, conv12, mask12 = self.__decoder_layer(320, 160, conv10, mask10, conv5, mask5, ['conv11', 'conv12'])
#     conv13, mask13, conv14, mask14 = self.__decoder_layer(160, 80, conv12, mask12, conv3, mask3, ['conv13', 'conv14'])
#     conv15, mask15, conv16, mask16 = self.__decoder_layer(80, 1, conv14, mask14, conv1, mask1, ['conv15', 'decoder_output'])

#     outputs = keras.layers.Conv2D(1, (3, 3), activation='sigmoid', padding='same')(conv16)

#     #conv5, mask5, conv6, mask6 = self.__decoder_layer(160, 80, conv4, mask4, conv3, mask3, ['conv5', 'conv6'])
#     #conv7, mask7, conv8, mask8 = self.__decoder_layer(80, 1, conv6, mask6, conv1, mask1, ['conv7', 'decoder_output'])

#     #outputs = keras.layers.Conv2D(1, (3, 3), activation='sigmoid', padding='same')(conv8)

#     return keras.models.Model(inputs=[input_image, input_mask], outputs=[outputs])

#   def __encoder_layer(self, filters, in_layer, in_mask, names):
#     conv1, mask1 = PConv2D(filters, (3,3), strides=1, padding='same', name=names[0])([in_layer, in_mask])
#     conv1 = keras.activations.relu(conv1)

#     conv2, mask2 = PConv2D(filters, (3,3), strides=2, padding='same', name=names[1])([conv1, mask1])
#     conv2 = keras.layers.BatchNormalization()(conv2, training=True)
#     conv2 = keras.activations.relu(conv2)

#     return conv1, mask1, conv2, mask2

#   def __decoder_layer(self, filter1, filter2, in_img, in_mask, share_img, share_mask, names):
#     up_img = keras.layers.UpSampling2D(size=(2,2))(in_img)
#     up_mask = keras.layers.UpSampling2D(size=(2,2))(in_mask)
#     concat_img = keras.layers.Concatenate(axis=3)([share_img, up_img])
#     concat_mask = keras.layers.Concatenate(axis=3)([share_mask, up_mask])

#     conv1, mask1 = PConv2D(filter1, (3,3), padding='same', name=names[0])([concat_img, concat_mask])
#     conv1 = keras.activations.relu(conv1)

#     conv2, mask2 = PConv2D(filter2, (3,3), padding='same', name=names[1])([conv1, mask1])
#     conv2 = keras.layers.BatchNormalization()(conv2)
#     conv2 = keras.activations.relu(conv2)

#     return conv1, mask1, conv2, mask2

# #model = InpaintingModel().prepare_model()


def inpaintLabelsFromScratched(scratchedLabelSet,modelType):
    # Prepare function calls based on which model is used
    if modelType == 'Simple':
        modelArchitectureFile = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_Simple.json'
        modelWeightsFiles = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_Simple.h5'
        newShape1 = (1,12,10,1)
        expandedShape = [1,16,16,1]
        finalLabelVectorLength = 120
        indexOfLabels = (slice(None),slice(2,14),slice(3,13))
    elif modelType == 'UpperGly':
        modelArchitectureFile = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_UpperGly.json'
        modelWeightsFiles = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_UpperGly.h5'
        newShape1 = (1,5,7,1)
        expandedShape = [1,8,8,1]
        finalLabelVectorLength = 35
        indexOfLabels = (slice(None),slice(2,7),slice(1,8))
    elif modelType == 'Gly13C2H':
        modelArchitectureFile = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_Glycolysis.json'
        modelWeightsFiles = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_Glycolysis.h5'
        newShape1 = (1,16,12,1)
        expandedShape = [1,16,16,1]
        finalLabelVectorLength = 192
        indexOfLabels = (slice(None),slice(0,16),slice(2,14))
    elif modelType == 'GlyPPP':
        modelArchitectureFile = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_GlyPPP.json'
        modelWeightsFiles = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_GlyPPP.h5'
        newShape1 = (1,72,48,1)
        expandedShape = [1,80,80,1]
        finalLabelVectorLength = 3456
        indexOfLabels = (slice(None),slice(4,76),slice(16,64))
    elif modelType == 'MammalianCCM':
        modelArchitectureFile = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_CCM.json'
        modelWeightsFiles = 'Trained_Models/PCNN_Label_Inpainting/PCONV_Inpainting_CCM.h5'
        newShape1 = (64,75,1)
        expandedShape = [1,80,80,1]
        finalLabelVectorLength = 4800
        indexOfLabels = (slice(None),slice(8,72),slice(3,78))
    else: raise ValueError("Unexpected model name. Accepted models are 'Simple','UpperGly','Gly13C2H','GlyPPP', or 'MammalianCCM'")
    
    # Generate the mask from the scratchedLabelSet
    mask = np.where(scratchedLabelSet == -1,1,0)

    # Load and compile models
    json_file = open(modelArchitectureFile, 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json, custom_objects = {"PConv2D": PConv2D})
    model.load_weights(modelWeightsFiles)
    model.compile(optimizer='adam', loss='mean_absolute_error', metrics=[dice_coef],run_eagerly=True)

    # Reshape scratched labels (to account for single or multiple datasets)
    scratchedLabelSet = np.reshape(scratchedLabelSet,(scratchedLabelSet.shape[0],-1))
    print(scratchedLabelSet.shape)
    scratchedLabelSet = np.transpose(scratchedLabelSet)
    print(scratchedLabelSet.shape)

    # Reshape and generate all expanded images
    all_scratchedLabelSet = np.empty((scratchedLabelSet.shape[0],expandedShape[1],expandedShape[2],expandedShape[3]))
    all_masks = np.empty((scratchedLabelSet.shape[0],expandedShape[1],expandedShape[2],expandedShape[3]))
    for i in range(scratchedLabelSet.shape[0]):
        singleScratchedLabelSet = scratchedLabelSet[i,:]

        # Generate the mask from the scratchedLabelSet
        mask = np.where(singleScratchedLabelSet == -1,1,0)
        mask = 1-mask

        # Reshape labeling and masks to fit model image shape
        singleScratchedLabelSet = singleScratchedLabelSet.reshape(newShape1)
        mask = mask.reshape(newShape1)
        masked_images_expand = np.ones(expandedShape)
        masks_expand = np.ones(expandedShape)
        masked_images_expand[indexOfLabels] = singleScratchedLabelSet
        masks_expand[indexOfLabels] = mask
        
        all_scratchedLabelSet[i,:] = masked_images_expand
        all_masks[i,:] = masks_expand

        print("Rehaped "+str(i+1)+" labelings out of "+str(scratchedLabelSet.shape[0]))

    # Inpaint (predict) labeling and reshape to proper output size
    y_pred = model.predict([all_scratchedLabelSet, all_masks])

    # Extract meaningful labeling data from each image
    all_fullLabelSet = np.empty(scratchedLabelSet.shape)
    for i in range(y_pred.shape[0]):
        single_y_pred = np.reshape(y_pred[i,:,:,:],expandedShape)
        fullLabelSet = single_y_pred[indexOfLabels]
        fullLabelSet = fullLabelSet.reshape(1,finalLabelVectorLength)
        all_fullLabelSet[i,:] = fullLabelSet
    
    return  all_fullLabelSet
