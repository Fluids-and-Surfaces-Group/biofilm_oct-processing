# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:12:09 2024

@author: Cornelius
"""


import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as ag1
import os
from OCT_utils import calculateHeightMap
from skimage import io
from matplotlib.patches import Rectangle



sourceFolder = "D:\\Iteration2\\registered" 
figureFolder = r"C:\Users\Cornelius\OneDriveKTH\Dokument\Promotion\python\OCT\plot_streamer_development"

# Define scans to be included
imageNrs = np.linspace(188,308,6).astype(np.int32) # Case 1 Channel 2 H2

"Define savepath"
savepath = os.path.join(figureFolder, "streamer_development_case1.eps")


# Axial resolution in m
dX = 12e-6
dY = 12e-6
dZ = 2.1e-6 

# Maximum height value
hMax = 300

# Number of images
nImages = len(imageNrs) 

# Number of columns
nCols = 2
nRows = int(np.ceil(nImages))

MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
plt.rc('figure', figsize=(9,12))  # Default size of figure
plt.rcParams['figure.dpi'] = 200

"Set up figure"
fig1 = plt.figure()
grid = ag1.AxesGrid(fig1, 111,  # similar to subplot(142)
                    nrows_ncols=(nRows, nCols),
                    axes_pad=(0.2,0.1),
                    share_all=False,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_pad=0.2,
                    cbar_size="1%",
                    )



for i,imageNr in enumerate(imageNrs):

    print("Generating map for file number " + str(imageNr))
    filename = os.path.join(sourceFolder,f'{imageNr:04}'+'_processed.tif')
    
    image = io.imread(filename)//255
    X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)

    # Select region of interest containing streamers
    X = X[0:800, 0:765]
    Y = Y[0:800, 0:765]
    hMap = hMap[0:800,0:765]*1e6 # Convert to microns
    
    "Position of streamers in mm"
    X1 = np.array([3.8,7.0])
    Y1 = np.array([7.0, 7.8])
    
    X2 = np.array([1.2, 4.8])
    Y2 = np.array([4.15, 5.15])
    
    Xlen = np.max([np.diff(X1), np.diff(X2)])
    Ylen = np.max([np.diff(Y1), np.diff(Y2)])
    
    X1[1] = X1[0] + Xlen 
    X2[1] = X2[0] + Xlen 
    
    Y1[1] = Y1[0] + Ylen
    Y2[1] = Y2[0] + Ylen

    "Corresponding indices"
    Xi1Start = (X1[0]//(dX*1e3)).astype(np.int32)
    Xi2Start = (X2[0]//(dX*1e3)).astype(np.int32)
   
    Yi1Start = 800-(Y1[1]//(dY*1e3)).astype(np.int32)
    Yi2Start = 800-(Y2[1]//(dY*1e3)).astype(np.int32)

    XiLen = (Xlen //(dX*1e3)).astype(np.int32)
    YiLen = (Ylen //(dY*1e3)).astype(np.int32)
    
    Xi1 = [Xi1Start, Xi1Start + XiLen]
    Xi2 = [Xi2Start, Xi2Start + XiLen]
    
    Yi1 = [Yi1Start, Yi1Start + YiLen]
    Yi2 = [Yi2Start, Yi2Start + YiLen]
    
    # Xi1 = (X1 //(dX*1e3)).astype(np.int32)
    # Yi1 = np.sort(800-(Y1 //(dY*1e3)).astype(np.int32))
    
    # Xi2 = (X2 //(dX*1e3)).astype(np.int32)
    # Yi2 = np.sort(800-(Y2 //(dY*1e3)).astype(np.int32))
    
    # Xlen = np.max([Xi1.max()-Xi1.min(), Xi2.max()-Xi2.min()])
    # Ylen = np.max([Yi1.max()-Yi1.min(), Yi2.max()-Yi2.min()])
    
    "Adjust lengths to actual sample"
    
    
    
    "Isolate two streamers and the respective ranges"
    imSlice1 = hMap[Yi1[0]:Yi1[1],Xi1[0]:Xi1[1]]
    xTickPos1 = np.linspace(0,imSlice1.shape[1]-1,num=5)
    xrange1 = np.round(np.linspace(X1[0],X1[1],num=5),2)
    yTickPos1 = np.linspace(0,imSlice1.shape[0]-1,num=3)
    yrange1 = np.linspace(Y1[0],Y1[1],num=3)

    imSlice2 = hMap[Yi2[0]:Yi2[1],Xi2[0]:Xi2[1]]
    xTickPos2 = np.linspace(0,imSlice2.shape[1]-1,num=5)
    xrange2 = np.round(np.linspace(X2[0],X2[1],num=5),2)
    yTickPos2 = np.linspace(0,imSlice2.shape[0]-1,num=3)
    yrange2 = np.linspace(Y2[0],Y2[1],num=3)

    # Plot regions of interest
    im = grid[2*i].imshow(imSlice1, vmax=hMax, interpolation='none')
    im2 = grid[2*i+1].imshow(imSlice2, vmax=hMax, interpolation='none')
    
    
    # Add titles
    if i==0:
        grid[0].set_title(r"Region 1")
        grid[1].set_title(r"Region 2")
        
    # Remove ticks
    if i == len(imageNrs)-1:
        grid[2*i].set_xticks([],[])
        grid[2*i+1].set_xticks([],[])
    
    
    "Add rectangles"
    X1start = 10
    Y1start = 30
    X1len = 80
    Y1len = 30
    
    X2start = 120
    Y2start = 30
    X2len = 70
    Y2len = 25
    
    grid[2*i].add_patch(Rectangle((X1start, Y1start), X1len, Y1len,edgecolor='r',
                  alpha=1, facecolor='none'))
    if i >= 1:
        grid[2*i+1].add_patch(Rectangle((X2start, Y2start), X2len, Y2len,edgecolor='r',
                      alpha=1, facecolor='none'))
    "Remove yticks for clarity"
    grid[2*i].set_yticks([],[])
    
  
cax = grid.cbar_axes[0].colorbar(im)    
cax.set_label(r"h [$\mu m$]")

plt.savefig(savepath, bbox_inches='tight')
