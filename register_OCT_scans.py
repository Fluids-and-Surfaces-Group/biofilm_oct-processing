# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 11:06:35 2023

@author: Cornelius
"""

"Image registration"
from skimage.registration import phase_cross_correlation
#from skimage.registration._phase_cross_correlation import _upsampled_dft
from skimage.util import img_as_ubyte
#from scipy.ndimage import fourier_shift
from skimage import io
import matplotlib.pyplot as plt
import numpy as np
import os

def register_OCT_scans(inputFolder,inputFiles,targetFolder,xRange,yRange,testLayer=3):#,firstFilenr,lastFilenr):
    # Takes filenames for a timeseries of OCT scans and uses cross-correlation along one plane to align 
    # the targetScan onto the refScan. Only translation is considered.
    
    shifts = []
    for i, file in enumerate(inputFiles):
        if i < len(inputFiles)-1:
            # Extracting file number from name. Needs to be changed depending on naming convention.
            #filenr = int(file[1:4])
            
            
            # Load reference image
            inputPath = os.path.join(inputFolder,file)
            print(inputPath)
            refData = io.imread(inputPath)
            refImage = refData[:,:,testLayer]
            refShape = refImage.shape
            
            # Load the following image as targetImage
            targetPath = os.path.join(inputFolder,inputFiles[i+1])#f'{filenr+12:04}'+'_processed.tif')
            print(targetPath)
            targetData = io.imread(targetPath)
            targetImage = targetData[:,:,testLayer]
            targetShape = targetImage.shape 
        
            # Check if image sizes are equal. If not, zero-pad smaller image / dimension on the right / bottom side in order not to influence the offset
            #print(refShape)
            #print(targetShape)
            if refShape[0]<targetShape[0]:
                refImage = np.pad(refImage,((0,targetShape[0]-refShape[0]),(0,0)))
            elif refShape[0]>targetShape[0]:
                targetImage = np.pad(targetImage,((0,refShape[0]-targetShape[0]),(0,0)))
                
            if refShape[1]<targetShape[1]:
                refImage = np.pad(refImage,((0,0),(0,targetShape[1]-refShape[1])))
            elif refShape[1]>targetShape[1]:
                targetImage = np.pad(targetImage,((0,0),(0,refShape[1]-targetShape[1])))
            #print(refImage.shape)
            #print(targetImage.shape)
            
            # Calculate pixel-accurate displacement
            shift, error, diffphase = phase_cross_correlation(refImage, targetImage)
            print(f'Detected pixel offset (y, x): {shift}, error {error}, dphase {diffphase}')
    
            shifts.append(shift)
                
    shifts = np.array(shifts,dtype=int)
    # Do not include! shifts = np.fliplr(shifts) # change to [xShift, yShift]
    totalShifts = np.cumsum(shifts,axis=0)
    
    
    # Find maximum displacement in positive and negative direction for both axes
    maxPosShift = [0, 0]
    if np.max(totalShifts,axis=0)[0]>0:
        maxPosShift[0] = np.max(totalShifts,axis=0)[0]
    if np.max(totalShifts,axis=0)[1]>0:
        maxPosShift[1] = np.max(totalShifts,axis=0)[1]
        
    maxNegShift = [0, 0]
    if np.min(totalShifts,axis=0)[0]<0:
        maxNegShift[0] = np.abs(np.min(totalShifts,axis=0)[0])
    if np.min(totalShifts,axis=0)[1]<0:
        maxNegShift[1] = np.abs(np.min(totalShifts,axis=0)[1])
    
    
    for i, file in enumerate(inputFiles):
        if i == 0:
            # Only trim original image, no shift
            xShift = yShift = 0    
        if i >= 1:
            # Translate all frames starting from the second time step
            xShift = totalShifts[i-1,0]
            yShift = totalShifts[i-1,1]
        
        
        inputPath = os.path.join(inputFolder,file)
        print(inputPath + "xShift: " + str(xShift) + "yShift: " + str(yShift))
        
        sourceData = io.imread(inputPath)
         

        # Roll data into correct position
        translatedData = np.roll(sourceData,xShift,axis=0)
        translatedData = np.roll(translatedData,yShift,axis=1)
        
        # Dimensions of data that is available in all time steps
        sourceDims = np.shape(sourceData)
        targetDimX = sourceDims[0]-maxPosShift[0]-maxNegShift[0]
        targetDimY = sourceDims[1]-maxPosShift[1]-maxNegShift[1]
        

        # Only keep the part that is available in all time steps
        translatedData = translatedData[maxPosShift[0]:maxPosShift[0]+targetDimX,
                                        maxPosShift[1]:maxPosShift[1]+targetDimY,
                                        :]
            
        # If necessary, trim images to relevant area¨
        channelNr = np.mod(i,12)
        translatedData = translatedData[yRange[0]:yRange[1],
                                        xRange[0]:xRange[1],:]
        
        
        print("Saving " + file)
        savePath = os.path.join(targetFolder,file)
        io.imsave(savePath,img_as_ubyte(translatedData))
            
            
    
    return

    

# Select dataset, iteration 2 is Re=100, iteration 3 is Re=300
iteration = 2

if iteration == 2:
    firstFilenr = 175
    lastFilenr = 462
    inputFolder = r"D:\Iteration2\processed"
    targetFolder = r"D:\Iteration2\registered"

    
    # Define ROI for each channel, manually decided to exclude regions with reflections or strong noise
    xRange = [[190, 932], [171,936],[58,941],[50,908],[287,912],[46,934],
              [46,931],   [167,712],[48,936],[50,941],[95,923], [40,915]]
              
              
    yRange = [[45, 1004], [15,886], [0,805], [0,800], [0,803],  [0,811],
              [0,810],    [4,1006], [5,717], [12,883],[0,800],  [6,803]]
              
      

elif iteration == 3:
    firstFilenr = 467
    lastFilenr = 731
    inputFolder = r"D:\Iteration3\processed_rescaled_merged\\"
    targetFolder = r"D:\Iteration3\registered\\"
    
    # Define ROI for each channel, manually decided to exclude regions with reflections or strong noise
    xRange = [[57,943], [193,926],[50,936],[38,921],[296,917],[147,934],
              [120,919],[76,894], [55,937],[42,860],[160,937],[26,836]]
              
    yRange = [[94,1007],[0,749],  [0,757], [0,736], [0,727],  [0,726],
              [0,755],  [0,769],  [0,734], [0,839], [0,746],  [0,751]]
        
else:
    print("No experiment of this name exists!")

# Set FOV to full to detect correct ranges. Comment out after.
# xRange=yRange = [[0,-1]]*12

for channel in range(12):
    # channel=1
    inputFiles = []
    for filenr in np.arange(firstFilenr+channel,lastFilenr+channel+1,12):
        if filenr >= 491:
           # Skip 491 and the 2 pm measurement
           filenr = filenr + 13
        if filenr <=731:          
            inputFiles.append(f'{filenr:04}'+'_processed.tif')
    
    register_OCT_scans(inputFolder, inputFiles, targetFolder, 
                       xRange[channel], yRange[channel])
  # fig = plt.figure(figsize=(8, 3))
  # ax1 = plt.subplot(1, 4, 1)
  # ax2 = plt.subplot(1, 4, 2, sharex=ax1, sharey=ax1)
  # ax3 = plt.subplot(1, 4, 3)
  # ax4 = plt.subplot(1, 4, 4)
  
  # ax1.imshow(image1, cmap='gray')
  # ax1.set_axis_off()
  # ax1.set_title('Reference image')
  
  # ax2.imshow(image2, cmap='gray')
  # ax2.set_axis_off()
  # ax2.set_title('Offset image')
  
  # # Show the output of a cross-correlation to show what the algorithm is
  # # doing behind the scenes
  # image_product = np.fft.fft2(image1) * np.fft.fft2(image2).conj()
  # cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
  # ax3.imshow(cc_image.real)
  # ax3.set_axis_off()
  # ax3.set_title("Cross-correlation")
  
  # ax4.imshow(refImage-targetImage)
  # ax4.set_axis_off()
  # ax4.set_title('Difference of images')
  
  # plt.show()
  
  # input("Press any key...")
  
