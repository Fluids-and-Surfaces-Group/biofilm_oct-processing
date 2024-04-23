# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:52:42 2023

@author: Cornelius
"""

"""Rescale measurements from 14 micron to 12 micron pixels. The rescaled scan is 
thresholded to minimise the change in biovolume."""

import numpy as np
from skimage import io
from skimage.transform import resize
from skimage.util import img_as_ubyte

import os

def main():
    
    folder = "D:\\Iteration3\\processed"
    saveFolder = "D:\\Iteration3\\rescaled"
    filenrs = np.arange(504,576)
    
    pixelSizeRef = 12.0
    pixelSizeWrong = 14.41
    
    
    for filenr in filenrs:
        filepath = os.path.join(folder,f'{filenr:04}'+'_processed.tif')
        print("Resizing " + str(filenr))
        image = io.imread(filepath)
        
        
        # Calculate shape after rescaling
        imShape = image.shape
        target_shape = (imShape[0]*pixelSizeWrong/pixelSizeRef, imShape[1], imShape[2])
                        #*pixelSizeWrong/pixelSizeRef,imShape[2])
        
        # Rescale image. This introduces grey values
        rescaledImage = resize(image, target_shape, order=1,preserve_range=True)
        
        # Calculate biovolume in source scan, one dimension with wrong scale
        volRef = np.sum(image)*pixelSizeWrong*pixelSizeRef/255
        
        # Calculate resulting biovolume depending on threshold value
        vols = np.zeros(255)
        error = np.zeros(255)
        for thresh in range (255):
            vols[thresh] = len(rescaledImage[rescaledImage>thresh])*pixelSizeRef**2
            error[thresh] = np.abs(volRef-vols[thresh])/volRef
        
        # Set threshold to minimise change in biovolume and re-binarise image
        threshChosen = np.argmin(error)
        print("Chosen threshold: " + str(threshChosen) + ", Remaining relative error: " + str(error[threshChosen]))
        finalImage = np.where(rescaledImage>thresh,255,0)
            
        savepath = os.path.join(saveFolder, f'{filenr:04}'+'_processed.tif')
        io.imsave(savepath,img_as_ubyte(finalImage), check_contrast=False)

if __name__ == "__main__":
    main()
