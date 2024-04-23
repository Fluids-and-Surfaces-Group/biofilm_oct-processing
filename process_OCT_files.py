# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:01:19 2023

@author: Cornelius
"""
from OCT_utils import find_substratum, level_substratum, threshold_OCT, remove_outliers, trimTop
from skimage.util import img_as_ubyte
from skimage.transform import resize
from skimage import io
import numpy as np
import os

folder = ("C:\\Users\\Cornelius\\OneDrive - KTH\\Dokument\\Promotion\\"
          "Measurements\\KIT\\OCT\\Iteration2\\rotated\\")
folder = ("D:\\Iteration2\\rotated\\")
saveFolder = ("D:\\Iteration2\\processed\\")

folder = "E:\Iteration3\\rotated\\"
saveFolder = "E:\Iteration3\\processed"

for filenr in range(491,577):
    print("Processing file number " + str(filenr))
    filename = os.path.join(folder,'Wittig_'+f'{filenr:04}'+'_Mode3D.tif')
    image = io.imread(filename)
    image = image.transpose((2,1,0))
    image = np.flip(image,axis=(0,1,2))
    
    print("Locating substratum")
    subPos = find_substratum(image,filter_radius=11)
    print("Levelling image")
    levelledImage = level_substratum(image,subPos)
    image = None
    
    "Do not do the following. Rescaling here messes up the thresholding and strongly increases the apparent biomass."
    # if filenr >= 491 and filenr <= 576:
    #     print("Rescaling image")
    #     imShape = levelledImage.shape
    #     target_shape = (imShape[0]*14/12, imShape[1]*14/12,imShape[2])
          
    #     rescaledImage = resize(levelledImage, target_shape, order=1,preserve_range=True)
    #     levelledImage = rescaledImage
      
    
    print("Thresholding")
    binImg, thresh = threshold_OCT(levelledImage)
    levelledImage = None
    
    print("Denoising")
    radius = 2
    denoisedImage = remove_outliers(img_as_ubyte(binImg),radius, 1)
    binImg = None
    
    print("Trimming")
    dZ = 2.1e-6
    maximumHeight = 5e-4
    trimmedImage, biofilmTop = trimTop(denoisedImage,maximumHeight,dZ)
    denoisedImage = None
    
        
    print("Saving")
    savePath = os.path.join(saveFolder,f'{filenr:04}'+'_processed.tif')
    io.imsave(savePath,img_as_ubyte(trimmedImage))
    trimmedImage = None
