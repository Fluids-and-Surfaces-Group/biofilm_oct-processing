# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:49:39 2024

@author: Cornelius
"""

"Plot the equivalent biofilm height on one day in all conditions. Use the replicates to create errorbars"


import pandas as pd
import seaborn as sns
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle

# =============================================================================
# Load dataframe
# =============================================================================
# Input folders containing the stored dataframes
folder1 = "D:\\Iteration2\\results"
folder2 = "E:\\Iteration3\\results"

# Name of the dataframe containing the biofilm statistics
pklFile = "stats_DataFrame_new.pkl"

# Destination folder for output figures 
figureFolder = r"C:\Users\Cornelius\OneDriveKTH\Dokument\Promotion\python\OCT\plot_biovolume_errorbars"


filepath1 = os.path.join(folder1,pklFile)
filepath2 = os.path.join(folder2,pklFile)

print("Loading data...")
# Load low Re and high Re cases
with (open(filepath1, "rb")) as openfile:
     dF1 = pickle.load(openfile)
     
with (open(filepath2, "rb")) as openfile:
     dF2 = pickle.load(openfile)
     

# Concatenate into one DataFrame
dF = dF2
dF = pd.concat([dF1, dF2], ignore_index=True)


# Means along same condition
# Need to average all unique combinations of day&tau_w
dFmeans = dF.groupby(['Day','tau_w'], as_index=False).mean(numeric_only=True)

print("""Dropping Channel 471 from mean calculations. Very spotty/whispy tall structures are visible that have a completely different morphology. """)
      
dFDrop = dF.drop(dF[dF["Channel"]==471].index)

"To drop the flow loop with much higher growth (183-186), uncomment the following block."
# print("""Dropping Channel 186 from mean calculations. Much larger biovolume. """)
# dFDrop = dFDrop.drop(dFDrop[dFDrop["Channel"]==186].index)
# print("""Dropping Channel 184 from mean calculations. Much larger biovolume. """)
# dFDrop = dFDrop.drop(dFDrop[dFDrop["Channel"]==184].index)
# print("""Dropping Channel 185 from mean calculations. Much larger biovolume. """)
# dFDrop = dFDrop.drop(dFDrop[dFDrop["Channel"]==185].index)
# print("""Dropping Channel 183 from mean calculations. Much larger biovolume. """)
# dFDrop = dFDrop.drop(dFDrop[dFDrop["Channel"]==183].index)

print("Normalising biovolume with FOV")
dFDrop["meanHeight"] = dFDrop["Biovolume"] /dFDrop["FOVArea"]*1000 # mean biofilm height in microns 


# =============================================================================
# Set plot parameters
# =============================================================================
tauVals = np.unique(dFDrop["tau_w"])

# Define colors for tau color grading
cmap = mpl.colormaps.get_cmap('plasma')
cSteps = np.linspace(0,255,num=6,dtype=np.int16)
tauColors = cmap(cSteps)

MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('axes', titlesize=BIGGER_SIZE)    # title fontsize
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels    
plt.rc('figure', figsize=(12,9))  # Default size of figure
plt.rcParams['figure.dpi'] = 100


# Define labels
ticklabels = []
colorlabels = []
for i,tauVal in enumerate(tauVals):
    colorlabels.append(r"$\tau_w =$ " + f'{tauVal:0.2}' + "Pa")
    ticklabels.append(i+1)

# =============================================================================
# # Create figure for mean height
# =============================================================================
for day in np.unique(dFDrop["Day"]):
    plt.figure()
    # sns.scatterplot(dFDrop[dFDrop["Day"]==day],x="tau_w",y="meanHeight",palette=cmap)#,label=colorlabels)#,hue="tau_w")#,label=colorlabels)
    sns.barplot(dFDrop[dFDrop["Day"]==day],x="tau_w",y="meanHeight",palette=tauColors,label=colorlabels)#,hue="tau_w")#,label=colorlabels)
    
    plt.xticks(ticks = range(6), labels=ticklabels)
    plt.xlabel(r"Case")
    plt.ylabel(r"$\overline{h}$ [$\mu$m]")
    
    savepath = os.path.join(figureFolder,"meanHeight_barplot_day" + str(day) + ".eps")
    plt.savefig(savepath, bbox_inches="tight")


# =============================================================================
# # Create figure for coverage
# =============================================================================
    
    plt.figure()
    # sns.scatterplot(dFDrop[dFDrop["Day"]==day],x="tau_w",y="Coverage",hue="tau_w")#,palette=tauColors)#,label=colorlabels)#,label=colorlabels)
    sns.barplot(dFDrop[dFDrop["Day"]==day],x="tau_w",y="Coverage",palette=tauColors,label=colorlabels)#,hue="tau_w")#,label=colorlabels)

    # plt.legend()
    plt.xticks(ticks = range(6), labels=ticklabels)
    plt.xlabel(r"Case")
    plt.ylabel(r"SC [-]")
    # plt.xlim([0,1.05*np.max(tauVals)])
    savepath = os.path.join(figureFolder,"coverage_barplot_day" + str(day) + ".eps")
    plt.savefig(savepath, bbox_inches="tight")
