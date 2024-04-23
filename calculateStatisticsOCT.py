# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:05:57 2023
Calculate statistics of individual biofilms. Show time evolution for each channel.

@author: Cornelius
"""
import os
import numpy as np
import utilities.biofilm_stats as bs
from utilities.OCT_utils import calculateHeightMap
from scipy.stats import skew, kurtosis
from skimage import io
from skimage.measure import label
from skimage.morphology import remove_small_objects

import pandas as pd


def calculate_statistics(folder,filenr,skipFiles,dX,dY,dZ,nFFT):
    " Calculate biofilm statistics for one file and return a dictionary"
    
    if filenr in skipFiles:
        # Skip faulty measurement and fill with nan. 
        biovol = np.nan
        biovolHeight = np.nan
        coverage = np.nan
        coverageMIP = np.nan
        filmConc = np.nan*np.ones((238,))
        h_avg = np.nan
        Ra = np.nan
        h_rms = np.nan
        ESx = np.nan
        ESy = np.nan
        skewness = np.nan
        kurt = np.nan
        meanAutocorrX = np.nan*np.ones((nFFT//2,))
        meanAutocorrY = np.nan*np.ones((nFFT//2,))
        corrLenX = np.nan
        corrLenY = np.nan
        medianHeight = np.nan
        meanKSpecX = np.nan*np.ones((nFFT//2+1,))
        meanKSpecY = np.nan*np.ones((nFFT//2+1,))
        meanRLX = np.nan
        meanRLY = np.nan
        runLengthsX = np.nan
        runLengthsY = np.nan
        meanFootprint = np.nan
        colonyHeights = np.nan
        colonyThickness = np.nan
        colonyLengths = np.nan
        colonyWidths = np.nan
        areas = np.nan
        mean_meanCS = np.nan
        centroidLines = np.nan
        leaningRates = np.nan
        leaningAngles = np.nan
        slenderness = np.nan
        
        
    else:    
        print("Generating stats for file number " + str(filenr))
        filename = os.path.join(folder,f'{filenr:04}'+'_processed.tif')
        
        # Area element for integration, in mm^2
        dA = dX * dY
        
        # Read image and convert true value to 1
        image = io.imread(filename)//255
        
        # Calculate biofilm height and thickness
        X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
        
        # Identify individual structures
        label_img = label(image)
        
        # Total measurement area
        A_t = np.size(hMap)*dA
           
        # Total biofilm volume in cubic millimetres
        biovol = np.sum(image)*dX*dY*dZ
        
        # Biovolume based on integrating the height map
        biovolHeight = np.sum(image)*dX*dY*dZ

        # Surface coverage
        coverage = bs.calculateCoverage(image,useMIP=False)
        coverageMIP = bs.calculateCoverage(image,useMIP=True)
        
        # Biofilm concentration by layer
        filmConcTemp = bs.filmConcentrationByHeight(image)
        
        # right-pad all entries to 238 entries, i.e. the maximum half-channel height
        filmConc = np.pad(filmConcTemp,[0,238-len(filmConcTemp)],constant_values=np.nan)
        
        # Calculate roughness parameters
        h_avg = bs.avgRoughnessHeight(hMap,dA,A_t)
        Ra = bs.averageHeightDeviation(hMap,h_avg,dA,A_t)
        h_rms = bs.RMSroughnessHeight(hMap,h_avg,dA,A_t)
        
        # Calculate effective slope; X-Axis is index 1
        ESx = bs.effectiveSlope(hMap, dA, A_t, 1)
        ESy = bs.effectiveSlope(hMap, dA, A_t, 0)
        skewness = skew(np.ravel(hMap))
        kurt = kurtosis(np.ravel(hMap))
        
        # Calculate autocorrelation along each axis; X-Axis is index 1
        meanAutocorrX = bs.average_1D_autocorr(hMap,nFFT=nFFT,axis=1)[:nFFT//2]
        meanAutocorrY = bs.average_1D_autocorr(hMap,nFFT=nFFT,axis=0)[:nFFT//2]
        
        # Calculate correlation length of current directional autocorrelations. 
        #corrLen: length at which decay reaches 1/e
        corrLenX = bs.corrLen(meanAutocorrX,nFFT,dX)
        corrLenY = bs.corrLen(meanAutocorrY,nFFT,dY)
        
        # Calculate median height
        medianHeight = bs.median_height(hMap)
        
        # Calculate average premultiplied wavenumber spectra in both directions
        meanKSpecX, kx = bs.average_1D_premultiplied_spectrum(hMap, dX, nFFT=nFFT, axis=1)
        meanKSpecY, ky = bs.average_1D_premultiplied_spectrum(hMap, dX, nFFT=nFFT, axis=0)
                
        # Average run lengths
        tmpRunLengthArea = bs.mean_run_lengths_area(label_img, dX, dY)
        
        meanRLX = np.mean(tmpRunLengthArea[0])
        meanRLY = np.mean(tmpRunLengthArea[1])
        
        meanFootprint = tmpRunLengthArea[2]
        
        runLengthsX = tmpRunLengthArea[0]
        runLengthsY = tmpRunLengthArea[1]
        areas = tmpRunLengthArea[3]
        
        # Colony heights and thicknesses
        colonyHeights, colonyThickness, colonyLengths, colonyWidths = bs.find_colony_bbox(label_img,dX,dY,dZ)
        
        
        # Mean of the average cross-section
        mean_meanCS = bs.calculate_mean_meanCS(label_img, dX)
        
        # For colony morphology analysis, drop tiny colonies
        cleanImage = remove_small_objects(image.astype(bool),min_size=3*3*3)
        cleanLabel_img = label(cleanImage)
        
        centroidLines = bs.find_centroid_lines(cleanLabel_img,dX,dY,dZ)
        
        leaningRates, leaningAngles = bs.calculate_leaning_rates_and_angles(centroidLines,dZ)
        
        slenderness = bs.calculate_slenderness(cleanLabel_img, dX, dZ, leaningAngles)
        
    channelDict = { 'Coverage': coverage, 
                    'CoverageMIP': coverageMIP, 
                    'Biovolume': biovol,
                    'BiovolHeight': biovolHeight,
                    'medianHeight': medianHeight, 
                    'Skewness': skewness,
                    'Kurtosis': kurt,
                    'FilmConcentration': [filmConc],
                    'RoughnessHeight': h_avg,
                    'Ra': Ra,
                    'h_rms': h_rms, 
                    'ES_x': ESx, 
                    'ES_y': ESy,
                    'meanAutocorrX': [meanAutocorrX],
                    'meanAutocorrY': [meanAutocorrY],
                    'corrLenX': corrLenX,
                    'corrLenY': corrLenY,
                    'meanKSpecX': [meanKSpecX],
                    'meanKSpecY': [meanKSpecY],
                    'meanRLX': meanRLX,
                    'meanRLY': meanRLY,
                    'meanFootprint': meanFootprint,
                    'runLengthsX': [runLengthsX],
                    'runLengthsY': [runLengthsY],
                    'areas': [areas],
                    'colonyHeights': [colonyHeights],
                    'colonyThickness': [colonyThickness],
                    'colonyLengths': [colonyLengths],
                    'colonyWidths': [colonyWidths],
                    'mean_meanCS': mean_meanCS,
                    'centroidLines': [centroidLines],
                    'leaningRates': [leaningRates],
                    'meanLeaningAngles': [leaningAngles],
                    'slenderness': [slenderness],
                             
                    }
        
    return channelDict  

# =============================================================================
# Start of the evaluation script
# =============================================================================
def main():
    # Decide which dataset to evaluate. iteration 2 is Re=100, iteration 3 is Re=300
    iteration = 2
    
    # Number of channels in the grid
    nChannels = 12
    
    # Lateral pixel size in mm
    dX = dY = 12e-3
    # Vertical pixel size in mm
    dZ = 2.1e-3 
    
    # Size of nFFT in autocorrelations
    nFFT = 1024
    
    
    
    # Channel heights in mm
    heights = np.array([1., np.sqrt(2), 2.])
    
    
    if iteration==2:
        # Source and target folders
        folder = "D:\\Iteration2\\registered"
        saveFolder = "D:\\Iteration2\\results"
        
        # Define measurement numbers for channels
        h1 = [177, 180, 183, 186]
        hs2 =  [175, 178, 181, 184]
        h2 = [176, 179, 182, 185]
        
        # Width of the channels in m
        width = 20e-3
        
        # Flow rate
        Q = 1e-6
        
        # Bulk velocity
        u_bulk = Q/(heights*20*1e-6)
        
        # Hydraulic diameter
        D_h = 2*(heights*20*1e-6)/((heights+20)*1e-3)
        
        # Viscosity at 24 C
        mu = 0.9096e-3
        
        # Density of water
        rho = 1000
        
        # Diffusivity, order of magnitude
        diffusivity = 1e-9
        
        # Calculate Reynolds number and estimate a Péclet number
        Re = u_bulk*D_h/mu*rho
        Pe = heights*1e-3*u_bulk/diffusivity
        
    # Define channels, iteration 3
    elif iteration == 3:
        folder = "D:\\Iteration3\\registered"
        saveFolder = "E:\\Iteration3\\results"
    
        h1 = [elem + 292 for elem in h1]
        hs2 = [elem + 292 for elem in hs2]
        h2 = [elem + 292 for elem in h2]
        
        # Flow rate
        Q = 2.63e-6
        mu = 0.8502e-3  # Viscosity at 27 C
    
        # Bulk velocity
        u_bulk = Q/(heights*20*1e-6)
    
        # Dimensionless numbers
        Re = u_bulk*D_h/mu*rho
        Pe = heights*1e-3*u_bulk/diffusivity
         
    else:
        print("Iteration not defined")
        
        
    "Estimate wall shear stress:"
    # Using just channel height
    tau_est_0 = 8 * mu * Q /(width*(heights*1e-3)**2) # from Darcy-Weisbach equation with f_d = 64/Re for laminar flow.
    # Using hydraulic diameter
    tau_est_1 = 8 * mu * (Q /(width*heights*1e-3)) / D_h # from Darcy-Weisbach equation with f_d = 64/Re for laminar flow.
    # From plane Poiseuille flow
    tau_est_2 = 6 * mu * (Q /(width*heights*1e-3)) / (heights*1e-3)
    
    # Using plane Poiseuille estimate
    tau_est = tau_est_2
    
    # Files to be skipped due to fibres / air bubbles etc. trapped in the FOV
    skipFiles = [195, 207, 219, 231, 243, 255, 267, 279, 291, 303, 315, 327, 339, 351, 363, 375, 387, 399, 411, 423, 435, 447, 459, # Bubbles
                 475,487, # Fibre and bubbles
                 614, 626, 638, 650, 662, 674, 686, 698, 710, 722, # Huge structure appears and got stuck
                 ]
    # Initialise dataFrame
    allData = pd.DataFrame(columns=['Channel','Height', 'tau_w','Re', 'Pe', 'Day', 
                                    'Coverage','Biovolume', 'medianHeight', 'Skewness', 
                                    'Kurtosis','FilmConcentration', 'RoughnessHeight', 
                                    'Ra', 'h_rms', 'ES_x', 'ES_y','meanAutocorrX',
                                    'meanAutocorrY','meanRLX','meanRLY','meanFootprint','runLengthsX','runLengthsY','areas'])
    
    # Definitions of calculated parameters:
    # biovol            # Biovolume in mm^3
    # coverage          # Substratum coverage
    # filmConc          # Biofilm concentration by layer
    # h_avg_list        # Average roughness height
    # Ra                # average height deviation
    # h_rms             # rms roughness height
    # ESx               # Effective slope in streamwise direction
    # ESy               # Effective slope in cross-channel direction
    # skewness          # Skewness of height map
    # kurt              # Kurtosis of height map
    # meanAutocorrX     # Average autocorrelation along X
    # meanAutocorrY     # Average autocorrelation along Y
    # corrLenX          # Correlation length along X in mm
    # corrLenY          # Correlation length along Y in mm
    # medianHeight      # Median height of existing biofilm in mm
    # meanKSpecX        # Average premultiplied wavenumber spectrum in X
    # meanKSpecY        # Average premultiplied wavenumber spectrum in Y
    # meanRLX           # Mean run length in streamwise direction in mm
    # meanRLY           # Mean run length in streamwise direction in mm
    # meanFootprint     # Mean footprint area of microcolonies in mm^2
    # meanArea          # Mean of the mean cross-sectional area of the microcolonies in mm^2
    # areas             # Footprint areas of microcolonies in px
    # centroidLines     # Centroid lines of microcolonies
    # leaningRates      # Leaning rates of centroid lines
    # meanLeaningAngles # Mean leaning angles of microcolonies in degrees
    # slenderness       # Slenderness of microcolonies   
    # colonyHeights       # Total heights of individual colonies in mm
    # colonyThickness     # Total vertical extent of individual colonies in mm
    # colonyLengths       # Total lengths of individual colonies in mm
    # colonyWidths        # Total widths of individual colonies in mm                 
        
    
    # Group channel heights    
    channelHeights = [h1, hs2, h2]
        
    filenumbers = []
    days = []
    channels = []
    heightIDs = []
    
    # Generate set of filenumbers to evaluate
    for heightID, channelHeight in enumerate(channelHeights):
            for channel in channelHeight:
                for timeStep in range(24):
                    # Get file number for current measurement
                    filenumber = channel + timeStep*nChannels
                    
                    day = timeStep / 2
                    
                    if filenumber >= 491:
                        # Skip measurement 491 and an additional set of measurements that was taken to diagnose mechanical issues
                        filenumber = filenumber + 13
                        
                    # Catch unrecorded files, including the jump in index
                    if filenumber <= 731:
                        filenumbers.append(filenumber)
                        days.append(day)
                        channels.append(channel)
                        heightIDs.append(heightID)
                        
    # Calculate statistics for the filenumbers                   
    for fileID, filenr in enumerate(filenumbers):
        
        # Get index for flow parameters of current channel
        filesPerHeight = len(channelHeights[0])*24
        heightID = fileID//(filesPerHeight)
    
        # Calculate biofilm statistics
        channelDataDict = calculate_statistics(folder,filenr,skipFiles,dX,dY,dZ,nFFT)  
        
        # Convert dicts to dataframe
        channelData = pd.DataFrame(channelDataDict)
        
        # Add current scan to full dataset
        allData = pd.concat([allData, channelData], ignore_index=True)
    
    # add metadata to dataframe
    allData['Channel'] = channels
    allData["Height"] = heights[heightIDs]*1e-3
    allData['tau_w'] = tau_est[heightIDs]
    allData['Re'] = Re[heightIDs] 
    allData['Pe'] = Pe[heightIDs]
    allData['Day'] = days
    
    # Save dataframe as pickle file 
    print("Saving DataFrame...")
    savePathPickle =  os.path.join(saveFolder,"stats_DataFrame_new.pkl")
    allData.to_pickle(savePathPickle)

if __name__ == "__main__":
    main()
