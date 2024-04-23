# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:56:17 2023

@author: Cornelius
"""

"Generate a 4x12 grid of height maps. One channel height vs time."


from OCT_utils import calculateHeightMap
import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as ag1
from skimage import io



plt.ioff()

iteration = 2

# Define channel heights
folder = "D:\\Iteration2\\registered"
saveFolder = "D:\\Iteration2\\results"

h1 = [177, 180, 183, 186]
hs2 =  [175, 178, 181, 184]
h2 = [176, 179, 182, 185]

# For SMILS poster
#h1 = [177, 178, 179]

# xRangeH1 = [[75, 945], [91, 947], [100, 943], [61, 928]]
# yRangeH1 = [[1, -1],   [1, -1],   [1, -1],    [1, -1]]

# xRangeHs2 = [[95, 953],  [60, 920], [100, 938],[70, 945] ]
# yRangeHs2 = [[45, 1004], [1, -1],   [1, -1],   [9, 1009]]

# xRangeH2 = [[60, 940],  [58, 925], [50, 913], [63, 932]]
# yRangeH2 = [[12, 1007], [1, -1],   [1, -1],   [1, -1]]
# Define channel heights, iteration 3
if iteration == 3:
    folder = "E:\\Iteration3\\processed"
    saveFolder = "E:\\Iteration3\\results"

    h1 = [elem + 292 for elem in h1]
    hs2 = [elem + 292 for elem in hs2]
    h2 = [elem + 292 for elem in h2]

    # xRangeH1 = [[100, 975], [120, 970], [111, 965], [93, 958] ]
    # yRangeH1 = [[90, 1008], [105, -1], [94, 1015], [84, 1003]]

    # xRangeHs2 = [[56, 871],  [270, 970], [310, 946], [250, 960]]
    # yRangeHs2 = [[97, 1000], [110, -1], [87, -1], [56, -1]]
    
    # xRangeH2 = [[68, 933],  [125, 980], [85, 967], [85, 953]]
    # yRangeH2 = [[66, 1009], [100, -1], [90, 1009], [91, 1010]]
# Axial resolution in m
dZ = 2.1e-6 

# Maximum height value
hMax = 500

# Number of channels
nChannels = 12 

# Use every skipNr'th time step
skipNr = 2

# Number of pages per series
pages = 1

heights = np.array([1., np.sqrt(2), 2.])

SMALL_SIZE = 17
MEDIUM_SIZE = 21
BIGGER_SIZE = 23

plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
plt.rc('figure', figsize=(8.27, 3*11.69)) # A4 size in inches
plt.rcParams ['figure.dpi'] = 600

plt.rcParams['axes.grid'] = False # Turn off axes

# 1 mm channels:

for page in range(pages):
    fig1 = plt.figure()
    grid1 = ag1.AxesGrid(fig1, 111,  # similar to subplot(142)
                        nrows_ncols=((24//skipNr)//pages, 4),
                        axes_pad=(0.3,0.1),
                        share_all=True,
                        label_mode="L",
                        cbar_location="bottom",
                        cbar_mode="single",
                        cbar_pad=0.4
                        )
    i=0
    for timeStep in range(0,24//pages,skipNr):
        for heightInd, channel in enumerate(h1):
            filenr = channel + timeStep*nChannels + page * 288 // pages
            if filenr >= 491:
                # Skip 491 and the 2 pm measurement
                filenr = filenr + 13
            # Lateral pixel size in m
            if filenr >= 491 and filenr <= 575:
                dX = dY = 14e-6
            else:
                dX = dY = 12e-6
            print("Generating map for file number " + str(filenr))
            filename = os.path.join(folder,f'{filenr:04}'+'_processed.tif')
            
            if filenr <= 731:
                image = io.imread(filename)//255
                X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
                # hMap = hMap[yRangeH1[0][0]:yRangeH1[0][1],xRangeH1[0][0]:xRangeH1[0][1]]
                
                if i<3:
                    im = grid1[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                                         np.max(Y)*1e3,np.min(Y)*1e3),
                                         vmax=hMax)
                    grid1[i].set_title("H = "+ f'{heights[heightInd]:1.2}'+ " mm")
                    firstFrame = False
                else: 
                    im = grid1[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                                         np.max(Y)*1e3,np.min(Y)*1e3),
                                         vmax=hMax)
                grid1[i].axis('off')
            i = i + 1
    
    #plt.title("1 mm channels, first half")
    #plt.xlabel("Cross-stream position in mm")
    #plt.ylabel("Streamwise position in mm")
    cax = grid1.cbar_axes[0].colorbar(im)    
    cax.set_label(r"Biofilm height [$\mu m$]")
    savepath = os.path.join(saveFolder, "grid_1mm_channels_part" + str(page) + ".png")
    plt.savefig(savepath, bbox_inches='tight')
plt.close(fig1)   
    

# 1.41 mm channels:

for page in range(pages):
    fig2 = plt.figure()
    grid2 = ag1.AxesGrid(fig2, 111,  # similar to subplot(142)
                        nrows_ncols=((24//skipNr)//pages, 4),
                        axes_pad=(0.3,0.1),
                        share_all=True,
                        label_mode="L",
                        cbar_location="bottom",
                        cbar_mode="single",
                        cbar_pad=0.4
                        )
    i=0
    for timeStep in range(0,24//pages,skipNr):
        for channel in hs2:
            filenr = channel + timeStep*nChannels + page * 288 // pages
            if filenr >= 491:
                # Skip 491 and 2 pm series
                filenr = filenr + 13
            # Lateral pixel size in m
            if filenr >= 491 and filenr <= 575:
                dX = dY = 14e-6
            else:
                dX = dY = 12e-6
            print("Generating map for file number " + str(filenr))
            filename = os.path.join(folder,f'{filenr:04}'+'_processed.tif')
            if filenr <= 731:  
                image = io.imread(filename)//255
                X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
                # hMap = hMap[yRangeHs2[0][0]:yRangeHs2[0][1],xRangeHs2[0][0]:xRangeHs2[0][1]]
                im = grid2[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                                         np.max(Y)*1e3,np.min(Y)*1e3),
                                       vmax=hMax)
            i = i + 1
    
    #plt.title("1 mm channels, first half")
    #plt.xlabel("Cross-stream position in mm")
    #plt.ylabel("Streamwise position in mm")
    grid2.cbar_axes[0].colorbar(im)
    savepath = os.path.join(saveFolder, "grid_1.4mm_channels_part" + str(page) + ".png")
    plt.savefig(savepath, bbox_inches='tight')
plt.close(fig2)   
   

# 2 mm channels:
for page in range(pages):
    fig3 = plt.figure()
    grid3 = ag1.AxesGrid(fig3, 111,  # similar to subplot(142)
                        nrows_ncols=((24//skipNr)//pages, 4),
                        axes_pad=0.1,
                        share_all=True,
                        label_mode="L",
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_pad=0.4
                        )
    i=0
    for timeStep in range(0,24//pages,skipNr):
        " Need to fix missing measurements on this level!!!!!!"
        for channel in h2:
            filenr = channel + timeStep*nChannels + page * 288 // pages
            if filenr >= 491:
                # Skip 491 and 2 pm series
                filenr = filenr + 13
            # Lateral pixel size in m
            if filenr >= 491 and filenr <= 575:
                dX = dY = 14e-6
            else:
                dX = dY = 12e-6
            print("Generating map for file number " + str(filenr))
            filename = os.path.join(folder,f'{filenr:04}'+'_processed.tif')
            if filenr <= 731:
                image = io.imread(filename)//255
                X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
                # hMap = hMap[yRangeH2[0][0]:yRangeH2[0][1],xRangeH2[0][0]:xRangeH2[0][1]]
                im = grid3[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                                         np.max(Y)*1e3,np.min(Y)*1e3),
                                       vmax=hMax)
            i = i + 1
    
    #plt.title("3 mm channels, first half")
    #plt.xlabel("Cross-stream position in mm")
    #plt.ylabel("Streamwise position in mm")
    grid3.cbar_axes[0].colorbar(im)
    savepath = os.path.join(saveFolder, "grid_2mm_channels_part" + str(page) + ".png")
    plt.savefig(savepath, bbox_inches='tight')
plt.close(fig3)   
plt.ion()  
