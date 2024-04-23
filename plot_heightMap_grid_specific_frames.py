# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 11:18:33 2024

@author: Cornelius
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as ag1
from matplotlib.patches import Rectangle
import os
from utilities.OCT_utils import calculateHeightMap
from skimage import io

condition = 1

# Define input and output folders
if condition in (1,2,4):
    sourceFolder = "D:\\Iteration2\\registered"  
elif condition in (3,5,6):
    sourceFolder = r"D:\Iteration3\registered" 
else:
    print("Condition does not exist")
    
figureFolder = r"C:\Users\Cornelius\OneDriveKTH\Dokument\Promotion\python\OCT\plot_heightMap_grid_specific_frames"

"Define measurements to be added to grid"
if condition == 1:
    imageNrs = np.linspace(188,308,6).astype(np.int32) # Case 1 Channel 2 H2
elif condition == 2:
    imageNrs = np.linspace(187,307,6).astype(np.int32) # case 2 Channel 1 HS2
elif condition == 3:
    imageNrs = np.linspace(486,606,6).astype(np.int32) # Case 3 Channel 8 H2
elif condition == 4:
    imageNrs = np.linspace(189,309,6).astype(np.int32) # case 4 Channel 3 H1
elif condition == 5:
    imageNrs = np.linspace(479,599,6).astype(np.int32) # case 5 Channel 1 HS2       
elif condition == 6:
    imageNrs = np.linspace(484,604,6).astype(np.int32) # case 6 Channel 6 H1
        
"Define savepath"
savepath = os.path.join(figureFolder, "timeSeries_case" + str(condition) +".eps")


"Skip 491 and the 2pm measurement"
if imageNrs[0]>400:
    imageNrs[1:] = imageNrs[1:]+13

# Axial resolution in m
dX = 12e-6
dY = 12e-6
dZ = 2.1e-6 

# Maximum height value for colorbar
hMax = 300

# Number of images
nImages = len(imageNrs) 

# Number of columns and rows
nCols = 3
nRows = int(np.ceil(nImages/nCols))

MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
plt.rc('figure', figsize=(12,9))  # Default size of figure
plt.rcParams['figure.dpi'] = 300

"Set up figure"
fig1 = plt.figure()
grid = ag1.AxesGrid(fig1, 111,  # similar to subplot(142)
                    nrows_ncols=(nRows, nCols),
                    axes_pad=(0.2,0.3),
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_pad=0.2,
                    )



for i,imageNr in enumerate(imageNrs):

    print("Generating map for file number " + str(imageNr))
    filename = os.path.join(sourceFolder,f'{imageNr:04}'+'_processed.tif')
    
    # Load image and calculate height map
    image = io.imread(filename)//255
    X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)

    
    # Get dimensions of height map
    Xmax = hMap.shape[0]
    Ymax = hMap.shape[1]
    
    # Plot current map
    im = grid[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                          np.min(Y)*1e3,np.max(Y)*1e3),
                      vmax=hMax, interpolation='none')
    
    
    grid[i].set_title("Day " + str(i+1))
    
    # For the last row, add xlabels
    if i > nCols * (nRows-1) - 1:
        grid[i].set_xlabel("x [$mm$]")
    
    # For the first column, add ylabels
    if np.mod(i,nCols) == 0:
        grid[i].set_ylabel("y [$mm$]")
    
    
    # Add rectangles around streamer regions in case 1, comment out for figure in appendix
    if imageNrs[0] == 188:
        ax = plt.gca()
        X1start = 3.8
        Y1start = 7.8
        X2start = 1.2
        Y2start = 5
        Xlen = 3.6
        Ylen = 1
        grid[i].add_patch(Rectangle((X1start, Y1start), Xlen, Ylen,edgecolor='r',
                      alpha=1, facecolor='none'))
        grid[i].add_patch(Rectangle((X2start, Y2start), Xlen, Ylen,edgecolor='tab:orange',
                      alpha=1, facecolor='none'))
        
for i in range(np.mod(nRows*nCols,nImages)):
    grid[-(i+1)].remove()

# Draw ticks and colorbar
plt.xticks()
grid[0].set_xticks(np.arange(np.min(X)*1e3, np.max(X)*1e3, 2.5))
grid[0].set_yticks(np.arange(np.min(Y)*1e3, np.max(Y)*1e3, 3.0))
cax = grid.cbar_axes[0].colorbar(im)    
cax.set_label(r"h [$\mu m$]")

plt.savefig(savepath, bbox_inches='tight')
