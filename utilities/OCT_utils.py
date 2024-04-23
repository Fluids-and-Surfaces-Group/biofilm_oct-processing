# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:02:44 2023
Contains functions for processing of OCT data.

@author: Cornelius
"""

import csv
import numpy as np
from scipy.signal import medfilt2d
from skimage.measure import block_reduce
from skimage.filters import rank
from skimage.morphology import disk, ball
from matplotlib import pyplot as plt
from PIL import Image


def find_substratum(image,filter_radius=11):
    "Find substratum by detecting the maximum intensity in each A-scan"
    substratumPos = np.argmax(image,axis=2)
    
    # Use median filter to remove holes in the substratum
    smoothSub = medfilt2d(substratumPos,kernel_size=filter_radius) 
    smoothSub = smoothSub.astype(np.uint16)
    return smoothSub


def level_substratum(image, substratumPos):
    "Build new image from content of image above substratumPos"
    "XXX Take care of areas outside the FOV as well as walls"
    levelledImage = np.zeros(np.shape(image),dtype=np.uint8)
    zDim = np.shape(image)[2]
    for x in range(np.shape(image)[0]):
        for y in range(np.shape(image)[1]):
            # Read cutoff point from substratum location
            cutoff = substratumPos[x,y]
            # Move 6 pixels to take width of substratum into account
            cutoff = cutoff + 6
            # Bound highest cutoff to zDim
            if cutoff > zDim-1:
                    cutoff = zDim-1
            # Enter values into array        
            levelledImage[x,y,:zDim-cutoff] = image[x,y,cutoff:]
    
    "Cut off upper half to remove volume that never contains data in this campaign."
    levelledImage = levelledImage[:,:,:zDim//2]
    return levelledImage
 
def threshold_OCT(rawImage):
    """Thresholding of OCT signals based on the shape of a typical OCT scan. 
    The signal is typically represented by a bump in the right side of a Gaussian 
    noise distribution. This method finds the first inflection point to the right
    of the mode of the histogram (excluding zero values). The value given as personal
    experience is required to reliably match the manually determined thresholds."""
    # Calculate histogram
    hist=np.histogram(rawImage,255,range=[0,255])
    
    # Calculate second derivative
    histDiff = np.diff(hist[0],n=2,prepend=0)
    
    # Find maximum location of histogram, excluding zero value from rotating/expanding
    intPeak = np.argmax(hist[0][1:])
    
    # Find maximum in histDiff to the right of intPeak 
    personalExperience = 3
    thresh = np.argmax(histDiff[intPeak:]) + intPeak + personalExperience
    
    # Mask image
    binImage = rawImage>thresh
    
    return binImage, thresh
 
def remove_outliers(image, radius=2, threshold=50):
    """"Removes outliers similar to the remove outliers function in FIJI. A 
    rank median filter is used to remove bright points that would disappear using
    the filter without filling in gaps."""
    # From https://forum.image.sc/t/using-remove-outliers-filter-in-pyimagej/55460/2
    footprint_function = disk if image.ndim == 2 else ball
    footprint = footprint_function(radius=radius)
    median_filtered = rank.median(image, footprint)
    outliers = (
        (image > median_filtered + threshold)
        | (image < median_filtered - threshold)
    )
    output = np.where(outliers, median_filtered, image)
    return output

def trimTop(image,maximumHeight,dZ):
    """Remove part of image above the first layer that does not contain any biofilm.
    If the biofilm appears to be larger than maximumHeight, cut there instead."""
    # Find first empty x/y-plane above the floor of the channel
    
    try:
        biofilmTop = np.nonzero(image.sum(axis=(0,1))==0)[0][0]
        if biofilmTop > maximumHeight//dZ:   
            raise Exception('Apparent height larger than expected')
            print('Apparent height larger than expected')
    except:
        # Variant 1: use slice of minimum signal
        # biofilmTop = np.argmin(image.sum(axis=(0,1)))
        
        # Variant 2: cut at maximum height above the substratum
        biofilmTop = int(maximumHeight//dZ)
        print("No empty slice found! Trimming at " + str(maximumHeight*1e3) + "mm.")
    # Remove part of image above top of biofilm
    trimmedImage = image[:,:,:biofilmTop]

    return trimmedImage, biofilmTop


def calculateHeightMap(image, dX, dY, dZ, downsampleFactor):
    "Calculates maps of biofilm-bulk-interface and biofilm thickness"
    
    # Read image, scale to [0,1]
    #image = io.imread(filename)//255

    # Reorder axes to x/y/z
    #image = image.transpose((1,2,0))
    
    # Flip vertical axis
    #image = np.flip(image,2)
    
    # Calculate height map in m
    thickness = np.sum(image,2)*dZ
    hMap = last_nonzero(image,2,0)*dZ 
    
    
    # Coarsen data
    if downsampleFactor != 1:
        thickness_coarse = block_reduce(thickness, block_size=(downsampleFactor,downsampleFactor), func=np.mean, cval=0)
        hMap_coarse = block_reduce(hMap, block_size=(downsampleFactor,downsampleFactor), func=np.mean, cval=0)
    else:
        thickness_coarse = thickness
        hMap_coarse = hMap
    
    # Generate X/Y grid
    X,Y = np.meshgrid(np.arange(np.shape(hMap_coarse)[1]),np.arange(np.shape(hMap_coarse)[0]))
    X = X*dX*downsampleFactor
    Y = Y*dY*downsampleFactor
    
    return X,Y,thickness_coarse,hMap_coarse

def plotHeightMap(X,Y,hMap,*args):
    "Plot the height map supplied in hMap. *args may be used to pass vmax to the colorbar."
    # args: vmax, label
    
    # Plot result
    # Define figure and font sizes
    SMALL_SIZE = 17
    MEDIUM_SIZE = 26
    BIGGER_SIZE = 30
    
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
    plt.rc('figure', figsize=(16,8))  # Default size of figure
    
    plt.figure()
    if args[0]:
        # Adjust colorbar limit
        im = plt.imshow(hMap, extent=(np.min(X)*1e3,np.max(X)*1e3,np.min(Y)*1e3,np.max(Y)*1e3), vmin=0, vmax=args[0])
        #im.cmap.set_under('m')
    else:
        plt.imshow(hMap, extent=(np.min(X)*1e3,np.max(X)*1e3,np.min(Y)*1e3,np.max(Y)*1e3))
    # Ajust labels as needed
    plt.colorbar(label=r"Biofilm height [$\mu m$]")
    plt.ylabel("Cross-stream position [$mm$]")
    plt.xlabel("Streamwise position in [$mm$]")
    plt.savefig(args[1], bbox_inches='tight')
    plt.close()

def plotThicknessOverHeight(X,Y,thickness,hMap,savepath):
    "Plot the biofilm thickness divided by the biofilm height"
    # args: vmax, label
    
    # Plot result
    # Define figure and font sizes
    SMALL_SIZE = 17
    MEDIUM_SIZE = 26
    BIGGER_SIZE = 30
    
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
    plt.rc('figure', figsize=(16,8))  # Default size of figure
    
    plt.figure()
    plt.imshow(thickness/hMap, extent=(np.min(X)*1e3,np.max(X)*1e3,np.min(Y)*1e3,np.max(Y)*1e3))
    plt.clim(0,1)
    # Ajust labels as needed
    plt.colorbar(label="Local solidity [-]")
    plt.ylabel("Cross-stream position [$mm$]")
    plt.xlabel("Streamwise position [$mm$]")
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()

def saveMapCSV(X,Y,hMap,filename):
    "Saves a given height map to a csv file."
    saveData = np.stack((np.ravel(X), np.ravel(Y), np.ravel(hMap)),axis=1)
    saveFile = open(filename, 'w+', newline ='') 
    with saveFile:     
        write = csv.writer(saveFile) 
        write.writerows(saveData) 
    
def first_nonzero(arr, axis, invalid_val=-1):
    """Find the first nonzero element of an array along axis. If none are found,
    invalid-val is returned. For OCT scans, this should typically be replaced by
    regionprops.bbox from skimage.measure. That function is much faster"""
    # from: https://stackoverflow.com/questions/47269390/how-to-find-first-non-zero-value-in-every-column-of-a-numpy-array
    mask = arr!=0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val) 

def last_nonzero(arr, axis, invalid_val=-1):
    """Find the last nonzero element of an array along axis. If none are found,
    invalid-val is returned. For finding bounding boxes in OCT scans, this should typically be replaced by
    regionprops.bbox from skimage.measure. That function is much faster"""
    # from: https://stackoverflow.com/questions/47269390/how-to-find-first-non-zero-value-in-every-column-of-a-numpy-array
    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def createTimelapse(images,saveName,fps):
    "Creates a gif out of the images supplied as a list of filenames."
    
    frames =  []
    for image in images:
        newFrame = Image.open(image)
        frames.append(newFrame)
    
    duration = 1/fps*1000
    frames[0].save(saveName, format='GIF', append_images=frames[1:], save_all=True,
                   duration=duration, loop=0)

