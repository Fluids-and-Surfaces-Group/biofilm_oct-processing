# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:10:27 2024

@author: Cornelius Wittig
"""
import numpy as np
from scipy.optimize import curve_fit
import warnings
from skimage.measure import regionprops

def fft_xcorr2D(image1, image2):
    """Calculate 2D cross correlation using fft / Wiener chinchin theorem"""
    # Subtract mean
    image1 = image1 - image1.mean()
    image2 = image2 - image2.mean()
    
    # Cross-correlation theorem
    xcorr = np.fft.irfft2(np.fft.rfft2(image1)*np.conj(np.fft.rfft2(image2)))

    # Return correlation of image1 and image2 
    return xcorr

def calculateCoverage(image, useMIP=False):
    """Calculate coverage of substratum at bottom layer"""
    
    coverage = np.count_nonzero(image[:,:,0])/np.size(image[:,:,0])
    
    if useMIP:
        "Alternative definition of coverage as portion of area that has biofilm above it"
        MIP = np.max(image,axis=2)
        coverage = np.count_nonzero(MIP)/np.size(MIP)
    return coverage

def filmConcentrationByHeight(image):
    """Calculate "coverage" of measurement area at each height"""
   
    concentration = np.zeros(np.shape(image)[2])
    for height in range(np.shape(image)[2]):
        concentration[height] = np.count_nonzero(image[:,:,height])/np.size(image[:,:,height])
    return concentration

def avgRoughnessHeight(hMap,dA,A_t):
    """Calculates the average roughness height of hMap using pixel area dA, and total area A_t."""
    h_avg = np.sum(hMap*dA)/A_t
    
    return h_avg

def averageHeightDeviation(hMap,h_avg,dA,A_t):
    """Calculate average height deviation of hMap using average roughness height 
    h_avg, pixel area dA, and total area A_t."""
    Ra = np.sum(np.abs(hMap-h_avg)*dA)/A_t
    
    return Ra

def RMSroughnessHeight(hMap,h_avg,dA,A_t):
    """Calculate RMS roughness height of hMap using average roughness height 
    h_avg, pixel area dA, and total area A_t."""
    h_rms = np.sqrt(np.sum((hMap-h_avg)**2*dA)/A_t)
    
    return h_rms

def effectiveSlope(hMap, dA, A_t, axis=0):
    "Calculate the effective slope of hMap in direction axis using pixel area dA and total area A_t."
    slope = np.diff(hMap,axis=axis)
    ES = 1/A_t * np.sum(np.abs(slope)*dA)
    
    return ES

def calculate_porosity(hMap,dY,k_c,A_t):
    "Calculate porosity according to Jouybari et al. (2021)"
    
    
    
def average_1D_autocorr(data,nFFT=1024,axis=0):
    "Calculate average of autocorrelations. Takes a 2D array as input."
    
    # Subtract mean from series
    data = data - np.mean(data, axis=axis, keepdims=True)
    
    # Normalise by standard deviation along axis
    # This often raises a divide by zero warning, because std ist zero if a 
    # row/column does not contain biofilm. Not scaling in this case is valid,
    # as the values stay at constant zero.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in divide")
        data = data/np.std(data,axis=axis,keepdims=True)
    
    
    # Calculate FFT
    dataFT = np.fft.rfft(data,n=nFFT,axis=axis)
    
    # Calculate autocorrelation via FFT
    autocorrs = np.fft.irfft(dataFT*np.conjugate(dataFT),n=nFFT,axis=axis)
    
    # Normalise by length of signal
    autocorrs /= autocorrs.shape[axis]

    # Average along the other axis
    mean1DAutocorr = np.mean(autocorrs,axis=1-axis)
    
    return mean1DAutocorr


def funcExp(t, a, tau, c):
    """Function for fitting autocorrelation decay"""
    return a * np.exp(-t / tau) + c

def corrLen(correlationData, nFFT, dX):
    """ Calculate correlation length, here defined as the distance [mm] at which the 
   decay of the autocorrelation reaches 1/e"""
    
    xData = np.arange(0,nFFT/2)*dX
    
    # Catch nan errors
    if np.isnan(correlationData).any():
        corrLen = np.nan
    else:
        aFit, tauFit, cFit = curve_fit(funcExp,xData,correlationData)[0]
        corrLen = tauFit # in mm
    
    return corrLen
    
def median_height(hMap):
    "Calculate median height of existing biofilm"
    
    biofilmHeights = hMap[hMap!=0]
    medianHeight = np.median(biofilmHeights)
    
    return medianHeight
    
def average_1D_premultiplied_spectrum(data,dX,nFFT=1024,axis=0):
    "Calculates the 1D premultiplied wavenumber spectrum"
    
    # Using nFFT instead to force common dimensionality. Zero-padding up to 1024 in most cases.
    # nx = np.shape(data)[axis]
        
    # Wave numbers in 1/m 
    kx = np.arange(nFFT/2+1)/(nFFT*dX) 
    
    # Subtract mean
    deviation = data - np.mean(data, axis=axis, keepdims=True);
    
    
    # Computing premultiplied PSD using FFT
    X = np.fft.rfft(deviation,n=nFFT,axis=axis)
    unscaledS = np.abs(X)**2/nFFT*dX
    
    # Average along the other axis
    meanUnscaledS = np.mean(unscaledS,axis=1-axis)
    
    meanS = meanUnscaledS * kx
    # Remove mirrored part, not needed because rFFT already does this
    # meanS = meanS[nFFT//2:]
    # kx = kx[nFFT//2:]
    
    return meanS, kx
  
def mean_run_lengths_area(label_img,dX,dY):
    """Calculates the mean run lengths of the biofilm colonies at the substratum 
    Also returns the mean area footprint of the colonies
    label_img takes a labelled image stack. Labelling can be done via cc3d.connected_components or skimage.measure.label"""
    
    #Find connected regions / colonies / structures at substratum level
    #labels = cc3d.connected_components(data[:,:,2], connectivity=8)
    
    # Consider one layer near the bottom for the footprint. 
    # Offset of two pixels is chosen to avoid influence of minor substratum roughness in certain channels.
    labels = label_img[:,:,2]
    regions = regionprops(labels)
    
    nColonies = np.unique(labels).size - 1 
    
    runLengthsX = np.zeros(nColonies)
    runLengthsY = np.zeros(nColonies)
    areas = np.zeros(nColonies)
    
    for i in range(nColonies):  
        
        # Find extreme points of each colony
        
        "My own implementation. Slow"
        # # i+1 to skip background element
        # startX = np.nanmin(first_nonzero(labels==i+1,axis=1,invalid_val=np.nan))
        # stopX = np.nanmax(last_nonzero(labels==i+1,axis=1,invalid_val=np.nan))
        
        # startY = np.nanmin(first_nonzero(labels==i+1,axis=0,invalid_val=np.nan))
        # stopY = np.nanmax(last_nonzero(labels==i+1,axis=0,invalid_val=np.nan))
        
        # # calculate maximum extent in X and Y, +1 is needed, see single pixel elements
        # runLengthsX[i] = (stopX-startX+1)
        # runLengthsY[i] = (stopY-startY+1)
        
        # areas[i] = np.sum(labels==i+1)
        
        "Regionprops"
        runLengthsX[i] = (regions[i].bbox[3] - regions[i].bbox[1] + 1)  #  + 1 to allow for one px long structures
        runLengthsY[i] = (regions[i].bbox[2] - regions[i].bbox[0] + 1)  
        
        # Contact patch sizes
        areas[i] = regions[i].area
        
        # regionprops.bbox contains (y_min, x_min, z_min, y_max, x_max, z_max) in physical coordinates

    # runLengthsX = runLengthsX.astype(np.int32)
    # runLengthsY = runLengthsY.astype(np.int32)
    # areas = areas.astype(np.int32)
    
    # Average run lengths and area over all colonies, mm or mm^2
    runLengthsX = runLengthsX*dX
    runLengthsY = runLengthsY*dY
    meanFootprint = np.mean(areas)*dX*dY
    
    # Expected values
    #expRLX,_ = np.histogram(runLengthsX,bins=np.max(runLengthsX),density=True)
    #RLXBins = np.arange(1,np.max(runLengthsX)+1)
    
    
    # plt.figure()
    # plt.bar(RLXBins,expRLX*RLXBins,align="center")
    # plt.show()
    
    # bins = np.logspace(0,3,99)
    # expAreas,_ = np.histogram(areas,bins=bins,density=True)
    # areasBins = np.arange(1,np.max(areas)+1)
    
    # plt.figure()
    # plt.loglog(bins[:-1], expAreas,'o')
    
    return runLengthsX, runLengthsY, meanFootprint, areas


def find_centroid_lines(label_img,dX,dY,dZ):
    "Finds the lines representing the centroids of all microcolonies"
    # For each colony, find the centroid in each slice.
    # Takes labelled image stack as input.    
    
    # Create empty list of lists for the centroid lines
    centroidLines = [ [] for _ in range(np.max(label_img)) ]
    
    # Iterate through all slices
    for sliceNr in range(label_img.shape[2]):
        # For each unique value larger than zero in the layer, calculate the centroid
        # Append to the respective centroid line
        imgSlice = label_img[:,:,sliceNr]
        regions = regionprops(imgSlice)
        
        # Iterate through all colonies present in this slice. Drop background (0) from unique
        for colonyInd, colony in enumerate(np.unique(imgSlice)[1:]):
            centroid = regions[colonyInd].centroid
            centroid = np.flip(centroid) # Switch to (x,y) from (y,x)
            centroid = centroid 
            centroidLines[colony-1].append(centroid) # -1 to get back to zero-indexed values
        
        
    # Convert list of list of tuples to list of arrays
    for elemInd, elem in enumerate(centroidLines):
        centroidLines[elemInd] = np.array(elem)
        centroidLines[elemInd] = centroidLines[elemInd] * dX # Convert index position to mm
    return centroidLines

def find_colony_bbox(label_img,dX,dY,dZ):
    """Calculates colony heights, i.e. highest position reached by a colony, and 
    colony thickness, i.e. the run length of the colony in the vertical direction.
    Also returns the width and length of each colony."""
    
    # regionprops.bbox contains (y_min, x_min, z_min, y_max, x_max, z_max) in physical coordinates
    
    # Calculate regionprops
    regions = regionprops(label_img)
    
    colonyHeights = np.zeros(len(regions)) 
    colonyThickness = np.zeros(len(regions)) 
    colonyLengths = np.zeros(len(regions)) 
    colonyWidths = np.zeros(len(regions)) 
    
    for colony in range(0,len(regions)):
        colonyHeights[colony] = regions[colony].bbox[5] 
        colonyThickness[colony] = regions[colony].bbox[5] - regions[colony].bbox[2] + 1 #  + 1 to allow for one px long structures 
        colonyLengths[colony] = regions[colony].bbox[4] - regions[colony].bbox[1] + 1  #  + 1 to allow for one px long structures
        colonyWidths[colony] = regions[colony].bbox[3] - regions[colony].bbox[0] + 1  #  + 1 to allow for one px long structures
       
    colonyHeights = colonyHeights * dZ
    colonyThickness = colonyThickness * dZ
    colonyLengths = colonyLengths * dX
    colonyWidths = colonyWidths * dY
    
    return colonyHeights, colonyThickness, colonyLengths, colonyWidths   

def calculate_leaning_rates_and_angles(centroidLines,dZ):
    "Calculates the leaning rates and mean leaning angles of the centroid lines away from the vertical"
    
    leaningRates = [] # Leaning rates of all colonies
    meanLR = np.ones((len(centroidLines),2))*np.nan # Average leaning rate of each colony
    for lineInd, line in enumerate(centroidLines):
        
        # If two vertical slices exist, calculate the angle
        if line.shape[0] > 1:
            leaningRate = np.diff(line, axis=0) / dZ # Calculate derivative of centroid position with regards to the height coordinate
            meanLR[lineInd] = np.mean(leaningRate,axis=0)
        else:
            leaningRate = []
        
        leaningRates.append(leaningRate)
        
    meanLeaningAngles = np.arctan(meanLR)/np.pi*180
    
    return leaningRates, meanLeaningAngles

def calculate_slenderness(label_img, dX, dZ, angles):
    """Calculates the slenderness of the microcolonies. Here, the slenderness is
    atypically defined as the ratio of the length of the centroid line of the 
    microcolony along the streamwise and vertical direction, and the streamwise 
    run length of the contact patch."""
    
    regions = regionprops(label_img)
    regions_bottom = regionprops(label_img[:,:,0])
    
    centroidLineLen = np.ones(len(regions)) * np.nan
    colonyLength = np.ones(len(regions)) * np.nan
    
    for colony in range(0,len(regions_bottom)):
        
        regionInd = regions_bottom[colony].label-1 # Indices of all colonies that touch the bottom
        
        # Length of centroid line
        colonyThickness = (regions[regionInd].bbox[5] - regions[regionInd].bbox[2]) * dZ # Vertical extent of the colony
        centroidLineLen[regionInd] = colonyThickness/np.cos(angles[regionInd][0])
        
        # Run length at substratum
        colonyLength[regionInd] = (regions_bottom[colony].bbox[3] - regions_bottom[colony].bbox[1]) * dX # Streamwise extent of the colony
    
    slenderness = np.abs(centroidLineLen) / colonyLength
    return slenderness

    
def calculate_mean_meanCS(label_img,dX):
    """Calculates the mean cross section of each microcolony and averages all of them."""
    
    nLayers = label_img.shape[2] 
    allAreas = np.ones((np.max(label_img),nLayers)) * np.nan
    
    for layer in range(nLayers):
        regions = regionprops(label_img[:,:,layer])
        
        layerAreas = [region.area for region in regions]   
        layerLabels = [region.label-1 for region in regions]
        allAreas[layerLabels,layer] = layerAreas 
        
    meanAreas = np.nanmean(allAreas,axis=1)
    
    mean_meanCS = np.nanmean(meanAreas)
    
    return mean_meanCS
    
