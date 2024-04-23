# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:10:06 2024

@author: Cornelius
"""

"Plot figures for KIT paper"


import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import ast
import pickle

# Converter function to restore numpy arrays from csv format
def from_np_array(array_string):
    array_string = array_string.replace('nan','np.nan')
    array_string = ','.join(array_string.replace('[ ', '[').split())
    return np.array(ast.literal_eval(array_string))

def linFun(x, a): # Linear function without intercept
    return a*x
    
def powFun(x, a, k): # Power function for power law
    return a*(x**(k))


def powFun2D(XY,a1,k1,k2): # Power funcion for 2D power law
    x,y = XY
    return a1 * (x**k1 * y**k2)


def powLinFun2D(XY,a1,k,c): # Power funcion in x and linear in y, with offset, for 2D power law
    x,y = XY
    return a1 * (x * y**k) + c

def powLinFun(XY,a1,k,): # Power funcion in x and linear in y, without offset
    x,y = XY
    return a1 * x * y**k

def powFun3D(XYZ,a1,k1,k2,k3): # Power funcion for 3D power law
    x,y,z = XYZ
    return a1 * (x**k1 * y**k2 * z**k3)

# =============================================================================
# Load data
# =============================================================================
def main():
    # Define input folders and files
    folder1 = "D:\\Iteration2\\results"
    folder2 = "E:\\Iteration3\\results"
    
    pklFile = "stats_DataFrame_new.pkl"
    
    filepath1 = os.path.join(folder1,pklFile)
    filepath2 = os.path.join(folder2,pklFile)
    
    # Define output folder
    figureFolder = r"C:\Users\Cornelius\OneDriveKTH\Dokument\Promotion\python\OCT\plot_figures_KIT_paper"
    
    
    print("Loading data...")
    # Load low Re and high Re cases
    with (open(filepath1, "rb")) as openfile:
         dF1 = pickle.load(openfile)
         
    with (open(filepath2, "rb")) as openfile:
         dF2 = pickle.load(openfile)
         
    
    # Concatenate into one DataFrame
    dF = dF2
    dF = pd.concat([dF1, dF2], ignore_index=True)
    
    # =============================================================================
    # Adjust dataframe
    # =============================================================================
        
    print("""Dropping Channel 471 from mean calculations. Very spotty/whispy tall 
          structures are visible that have a completely different morphology. """)
    dFDrop = dF.drop(dF[dF["Channel"]==471].index)
    
    print("Normalising biovolume with FOV")
    dFDrop["meanHeight"] = dFDrop["Biovolume"] /dFDrop["FOVArea"]*1000 # mean biofilm height in microns 
    
    # Average over replicates
    dFDropMeans = dFDrop.groupby(['Day','tau_w'], as_index=False).mean(numeric_only=True)
    
    
    # =============================================================================
    # Plot results
    # =============================================================================
    
    # Define colors for tau color grading
    cmap = mpl.colormaps.get_cmap('plasma')
    cSteps = np.linspace(0,255,num=6,dtype=np.int16)
    tauColors = cmap(cSteps)
    
    cSteps= np.linspace(0,255,num=16,dtype=np.int16)
    # dayColors = cmap(cSteps)
    
    # Biovolume over time
    fParsTau = []
    plt.figure()
    for cStep, tauVal in enumerate(dFDropMeans["tau_w"].drop_duplicates().sort_values()):
        # Get time series for case with tau=tauVal
        singleSeries = dFDropMeans[(dFDropMeans["Day"]<=7) & (dFDropMeans["tau_w"]==tauVal)]    
        
        # Linear least-squares fit for each condition
        fPars, pcov = curve_fit(linFun,singleSeries["Day"],singleSeries["meanHeight"])
        fParsTau.append(fPars)
        
        # Plot data
        plt.plot(singleSeries["Day"],singleSeries["meanHeight"],"s-",color=tauColors[cStep],label=r"$\tau_w = $" + f'{tauVal:0.2}' + "Pa")
    
        # Plot linear fit
        plt.plot(singleSeries["Day"],linFun(singleSeries["Day"],fPars[0]),"--",color=tauColors[cStep])        
    
    # Save growth rates for later use
    fParsTau = np.vstack(fParsTau)
    
    # Get axis limits
    ylim = plt.ylim()
        
    plt.xlabel("Day")
    plt.ylabel(r"$\langle\overline{T}\rangle $ [$\mu m$]")
    plt.legend()
    
    savepathUnscaledBiovol = os.path.join(figureFolder,"meanHeightOverTimeFit.eps")
    plt.savefig(savepathUnscaledBiovol, bbox_inches="tight")
    plt.show()
    
    "-----------------------------------------------------------------------------"
    # Biovolume over tau, use every third time step only, starting at 0.5 days
    plt.figure()
    for cStep, dayVal in enumerate(dFDropMeans["Day"].drop_duplicates().sort_values()[1::3]):
        if cStep < len(tauColors) and dayVal <= 7:
            # Get time series for case with tau=tauVal
            singleSeries = dFDropMeans[(dFDropMeans["Day"]<=7) & (dFDropMeans["Day"]==dayVal)]
            
            # Plot data        
            plt.plot(singleSeries["tau_w"],singleSeries["meanHeight"],"s-",color=tauColors[cStep],label=r"$t = $" + str(dayVal) + "d")
        
    
    plt.ylim(ylim) # Setting ylim equal to previous figure
    plt.xlabel(r"$\tau_w$ [Pa]")
    plt.ylabel(r"$\langle\overline{T}\rangle $ [$\mu m$]")
    plt.legend()
    
    savepathUnscaledBiovol = os.path.join(figureFolder,"meanHeightOverTauFit.eps")
    plt.savefig(savepathUnscaledBiovol, bbox_inches="tight")
    plt.show()
    
    "-----------------------------------------------------------------------------"
    # coverage over time
    plt.figure()
    for cStep, tauVal in enumerate(dFDropMeans["tau_w"].drop_duplicates().sort_values()):
        # Get time series for case with tau=tauVal
        singleSeries = dFDropMeans[(dFDropMeans["Day"]<=7) & (dFDropMeans["tau_w"]==tauVal)]
    
        plt.plot(singleSeries["Day"],singleSeries["Coverage"],"s-",color=tauColors[cStep],label=r"$\tau_w = $" + f'{tauVal:0.2}' + "Pa")
        
    plt.xlabel("Day")
    plt.ylabel(r"$\langle SC \rangle $[-]")
    plt.legend()
    plt.ylim((-0.01, 0.35)) 
    savepathUnscaled_h_rms = os.path.join(figureFolder,"coverage_unscaled.eps")
    plt.savefig(savepathUnscaled_h_rms, bbox_inches="tight")
    plt.show()
    
    "-----------------------------------------------------------------------------"
    
    # coverage over tau, use every third time step only, starting at 0.5 days
    plt.figure()
    for cStep, dayVal in enumerate(dFDropMeans["Day"].drop_duplicates().sort_values()[1::3]):
        if cStep < len(tauColors) and dayVal <= 7:
            # Get time series for case with tau=tauVal
            singleSeries = dFDropMeans[(dFDropMeans["Day"]<=7) & (dFDropMeans["Day"]==dayVal)]
    
            plt.plot(singleSeries["tau_w"],singleSeries["Coverage"],"s-",color=tauColors[cStep],label=r"$t = $" + f'{dayVal}' + "d")
    
    plt.xlabel(r"$\tau_w$ [Pa]")
    plt.ylabel(r"$\langle SC \rangle $[-]")
    plt.legend()
    plt.ylim((-0.01, 0.35)) # set ylim equal to figure above
    savepathUnscaledBiovol = os.path.join(figureFolder,"CoverageOverTau_unscaled.eps")
    plt.savefig(savepathUnscaledBiovol, bbox_inches="tight")
    plt.show()
    
    "-----------------------------------------------------------------------------"
    """Calculate the growth rate sigma for each tau condition, then fit the resulting 
    growth rates to sigma(tau) = a * tau ^ k_tau.
    """
    
    # Get unique values of tau_w    
    tauVals = dFDropMeans["tau_w"].drop_duplicates().sort_values()
    
    # Compare r-squared value for power-law fit and 1/tau scaling
    fParsSigmaPLaw, pcoSigmaPLaw = curve_fit(powFun,tauVals,fParsTau[:,0])
    r2Fit = r2_score(fParsTau[:,0],powFun(tauVals,fParsSigmaPLaw[0],fParsSigmaPLaw[1]))
    r2Estimate = r2_score(fParsTau[:,0],powFun(tauVals,0.25,-1))
    print("Fit result for growth rate over tau, r-squared = " + str(r2Fit))
    print("1/4 tau^-1 estimate for growth rate over tau, r-squared = " + str(r2Estimate))
    
    
    
    # Plot growth rates and 0.25*1/tau_w scaling
    plt.figure()
    plt.plot(tauVals,fParsTau[:,0],'-s',color=tauColors[0])
    plt.plot(np.linspace(tauVals.min(),tauVals.max()),powFun(np.linspace(tauVals.min(),tauVals.max()),0.25,-1),color=tauColors[3])
    
    plt.xlabel(r"$\tau_w [Pa]$")
    plt.ylabel(r"$\sigma$ [$\mu$m/day]")
    plt.legend(["Growth rates",r"$\tau_w^{-1}$-scaling"])
    
    savepath = os.path.join(figureFolder,"growthRates.eps")
    plt.tight_layout()
    plt.savefig(savepath)
    

if __name__ == "__main__":
    main()
