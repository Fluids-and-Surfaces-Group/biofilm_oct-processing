# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 11:18:33 2024

@author: Cornelius
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpl_toolkits.axes_grid1 as ag1
import os
from utilities.OCT_utils import calculateHeightMap
from skimage import io



def main(imageNrs):
    "Define heightMaps to be added to grid"
    
    # Define source folder depending on iteration
    if np.max(imageNrs) < 465 :
        sourceFolder = "D:\\Iteration2\\registered" 
    else:
        sourceFolder = "E:\\Iteration3\\processed"
    
    # Define output folder
    figureFolder = r"C:\Users\Cornelius\OneDriveKTH\Dokument\Promotion\python\OCT\plot_two_maps_equal_FOV"
    
    # Axial resolution in m
    dX = 12e-6
    dY = 12e-6
    dZ = 2.1e-6 
    
    # Number of images
    nImages = len(imageNrs) 
    
    # Number of columns and rows
    nCols = 2
    nRows = int(np.ceil(nImages/nCols))
    
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16
    
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
    plt.rc('figure', figsize=(12,9))  # Default size of figure
    plt.rcParams['figure.dpi'] = 200
    
    
    "----------------------------Plot hMaps---------------------------------------"
    # Maximum height value for colormap
    hMax = 300
    
    
    "Set up figure"
    fig1 = plt.figure()
    grid = ag1.AxesGrid(fig1, 111,  # similar to subplot(142)
                        nrows_ncols=(nRows, nCols),
                        axes_pad=(0.4,0.5),
                        share_all=True,
                        label_mode="L",
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_pad=0.2,
                        )
    
    
    
    for i,imageNr in enumerate(imageNrs):
    
        print("Generating map for file number " + str(imageNr))
        filename = os.path.join(sourceFolder,f'{imageNr:04}'+'_processed.tif')
        
        # Load scan and generate height map
        image = io.imread(filename)//255
        X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
        
        
        "Trim both figures to one size"
        X = X[0:800, 0:765]
        Y = Y[0:800, 0:765]
        hMap = hMap[0:800,0:765]
        print(hMap.shape)
        
        # Plot height maps
        im = grid[i].imshow(hMap*1e6, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                              np.min(Y)*1e3,np.max(Y)*1e3),
                            vmax=hMax, interpolation='none')
        # Add titles
        if i==0:
            grid[i].set_title(r"Case 1, $\tau_w=0.068$ Pa") 
        else:
            grid[i].set_title(r"Case 4, $\tau_w=0.27$ Pa")
        "For the last row, add xlabels"
        if i > nCols * (nRows-1) - 1:
            grid[i].set_xlabel("x [$mm$]")
        
        "For the first column, add ylabels"
        if np.mod(i,nCols) == 0:
            grid[i].set_ylabel("y [$mm$]")
    
    # Add hightlights for detail slices    
    grid[0].hlines(6.6,6,7.2,color='r')
    grid[1].hlines(5.4,6,7.2,color='r')
    
    grid[0].set_xticks([0,3,6,9])
    grid[0].set_yticks([0,3,6,9])
    
    # Add colorbar
    cax = grid.cbar_axes[0].colorbar(im)    
    cax.set_label(r"h [$\mu m$]")
    
    savepath = os.path.join(figureFolder, "hMaps_comparison.eps")
    plt.savefig(savepath, bbox_inches='tight')
    
    "----------------------------Plot hOverThick----------------------------------"
    
    # Maximum height value
    hMax = 1
    
    
    "Set up figure"
    fig2 = plt.figure()
    grid2 = ag1.AxesGrid(fig2, 111,  # similar to subplot(142)
                        nrows_ncols=(nRows, nCols),
                        axes_pad=(0.4,0.5),
                        share_all=True,
                        label_mode="L",
                        cbar_location="right",
                        cbar_mode="single",
                        cbar_pad=0.2,
                        )
    
    
    
    for i,imageNr in enumerate(imageNrs):
    
        print("Generating map for file number " + str(imageNr))
        filename = os.path.join(sourceFolder,f'{imageNr:04}'+'_processed.tif')
        
        # Load scan and calculate height map
        image = io.imread(filename)//255
        X, Y, thickness, hMap = calculateHeightMap(image, dX, dY, dZ, 1)
        
        # Calculate solodity measure thickness over height
        ToH = thickness/hMap
        
        "Trim both figures to one size"
        X = X[0:800, 0:765]
        Y = Y[0:800, 0:765]
        ToH = ToH[0:800, 0:765]
        
        print(ToH.shape)
        
        # Plot map using white as bad-color to avoid issues with transparency in .eps files
        cmap = mpl.colormaps.get_cmap('viridis')  
        cmap.set_bad(color='white')    
        im = grid2[i].imshow(ToH, extent=(np.min(X)*1e3,np.max(X)*1e3,
                                              np.min(Y)*1e3,np.max(Y)*1e3),
                            vmin = 0, vmax=1, cmap=cmap, interpolation='none')
        
        # Add titles
        if i==0:
            grid2[i].set_title(r"Case 1, $\tau_w=0.068$ Pa")
    
        else:
            grid2[i].set_title(r"Case 4, $\tau_w=0.27$ Pa")
            
        "For the last row, add xlabels"
        if i > nCols * (nRows-1) - 1:
            grid2[i].set_xlabel("x [$mm$]")
        
        "For the first column, add ylabels"
        if np.mod(i,nCols) == 0:
            grid2[i].set_ylabel("y [$mm$]")
        
     
    # Set nice ticks   
    grid2[0].set_xticks([0,3,6,9])
    grid2[0].set_yticks([0,3,6,9])
    
    # Add colorbar
    cax = grid2.cbar_axes[0].colorbar(im)    
    cax.set_label(r"T/h [-]")
    
    savepath = os.path.join(figureFolder, "ToH_comparison.eps")
    plt.savefig(savepath, bbox_inches='tight')
    
    "--------------------------Plot typical structure-----------------------------"
    
    
    "Set up figure"
    fig3 = plt.figure()
    grid3 = ag1.AxesGrid(fig3, 111,  # similar to subplot(142)
                        nrows_ncols=(nRows, nCols),
                        axes_pad=(0.4,0.5),
                        share_all=True,
                        label_mode="L",
                        )
    
    
    for i,imageNr in enumerate(imageNrs):
        
        print("Generating map for file number " + str(imageNr))
        filename = os.path.join(sourceFolder,f'{imageNr:04}'+'_processed.tif')
        
        # Load scan
        image = io.imread(filename)//255
        
        # For each scan, select the relevant vertical slice containing the typical structures
        if i==0:
            imSlice = image[247,500:600,:].T
            xTickPos = np.linspace(0,imSlice.shape[1]-1,num=5)
            xrange = np.linspace(X[0,500]*1e3,X[0,600]*1e3,num=5)
            yTickPos = np.linspace(0,(imSlice.shape[0]),num=3)
            yrange = np.int32(np.round(np.linspace(0,(imSlice.shape[0])*dZ*1e6,num=3)))
            
            # Fill right plot with black background to ensure equal plot size
            grid3[1].imshow(np.zeros(imSlice.shape),cmap="binary_r",aspect=dZ/dX, interpolation='none')
        else:
            imSlice = image[343,500:600,:].T
    
        # Plot structures in correct aspect ratio
        grid3[i].imshow(imSlice,origin='lower',cmap="binary_r",aspect=dZ/dX, interpolation='none')
        
        if i==0:
            grid3[i].set_title(r"Case 1, $\tau_w=0.068$ Pa")
    
        else:
            grid3[i].set_title(r"Case 4, $\tau_w=0.27$ Pa")
            
        "For the last row, add xlabels"
        if i > nCols * (nRows-1) - 1:
            grid3[i].set_xlabel("x [$mm$]")
        
        "For the first column, add ylabels"
        if np.mod(i,nCols) == 0:
            grid3[i].set_ylabel("y [$\mu m$]")
        
    
    grid3[0].set_xticks(xTickPos,xrange)
    grid3[0].set_yticks(yTickPos,yrange)
    
    savepath = os.path.join(figureFolder, "typicalStructures.eps")
    plt.savefig(savepath, bbox_inches='tight')
 
    
imageNrs = [332,333]       
if __name__ == "__main__":
    main(imageNrs)
 
