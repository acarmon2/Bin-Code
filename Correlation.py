#!/usr/bin/env python
# -*- coding: utf-8 -*-

#########################################################################################################
# functions to calculate the cross correlation in time domain, in Spectral.csv has the power response
# of spectral parameters. 
# Pandas is the library to read .csv data (comma separate dataset, not point and comma), each row has the
# specral response of input pulse with different lambda separation
# matplotlib contain the basis to plot
import numpy as np
import pandas as pd
from scipy import special
import matplotlib.pyplot as plt
import time
from Bin import * 

##########################################################################################################
# Global calculations to extract the values of response for FBG
Spectral0 = pd.read_csv("Spectral.csv")
# convert dataframe to numpy array, first column is time, second delta lambda = 0, third delta lambda = 1, ...
Specter = Spectral0.values
Ns = np.shape(Specter)[0]                                                                               # Number of data in symbol section                                                                                                  # Number of data for delay sections
DeltaM = np.shape(Specter)[1]

##########################################################################################################
# Main functions for correlation part, Ns corresponds to number of points to power section spectral
# differences is contained in Specter matrix (numpy array) and NBi is a vector with the distribution of blank 
# spaces, Mi contains L numbers and NBi contains L-1 numbers, but add one more to easy implementation
def CrossCorr(C0, M1, NB0, MB1):
    global Specter
    C1 = M1[::-1]                                                                                       # The detector is a correlational but with inverter order
    NB1 = MB1[::-1]                                                                                     # Inverter the order of spaces too
    L0 = np.shape(C0)[0]                                                                                # Length of codes
    NB0 = np.append(NB0, [0])                                                                           # No blank spaces adjacent to final Bragg grating
    NB1 = np.append(NB1, [0])
    L1 = np.sum(NB1)                                                                                    # Total blank spaces
    Npoints = np.shape(Specter)[0] 
    Dim0 = (L0*Npoints) + L1
    Dim1 = (2*Dim0)
    b = 0
    VectorCC = np.zeros(Dim1, dtype = float)                                                            # Cross Correlation vector
    # Run until the L0-1 because the final FBG does not have delay associated
    for i in range(0, L0):
        Vector0 = np.zeros(Dim0, dtype = float)                                                         # Spectral response for individual symbol over all the grid
        a = 0
        for j in range(0, L0):
            if((C0[i] >= 0) and (C1[j] >= 0)):                                                          # Case of not insert blank spaces like Kim publication
                delt = np.abs(C0[i] - C1[j])                                                            # Spectral distance
                Vector0[a:(a+Npoints)] = Vector0[a:(a+Npoints)] + np.transpose(Specter[:,[delt+1]])     # insert the power response of FBG
                a = a + Npoints + NB1[j]                                                                # Insert Blank spaces
            else:                                                                                       # Kim case
                a = a + Npoints + NB1[j]                                                              
        # Add the propagation of symbol to the crosscorrelation function
        VectorCC[b:(b+Dim0)] = VectorCC[b:(b+Dim0)] + Vector0[:]                                        # Sum of the power (Vector of crosscorelation)
        b = b + Npoints + NB0[i]
    # return the cross correlation deleting the final space ot last one FBG
    return VectorCC

##########################################################################################################
# Function Stats to produce the final probability error based in Q function. Gen is the generator function
# q0 correspond to the maximum number of users. BTest is the default profile of blank spaces between 
# adjacent FBG
def Stats(Gen, q0, BTest):
    L0 = np.shape(Gen)[0]
    Codes = UCodes(Gen, q0)
    # Values to normalize the crosscorrelation values
    ACorrMAx = np.zeros(q0, dtype = float)
    # Mean power for cross correlation detector (var + average^2)
    PowM = np.zeros(q0, dtype = float)
    for i in range(0, q0):
        CrossVal = CrossCorr(Codes[i], Codes[i], BTest[i], BTest[i])
        ACorrMAx[i] = np.amax(CrossVal)
        CrossVal = CrossVal/ACorrMAx[i]
        PowM[i] = np.var(CrossVal) + ((np.average(CrossVal))**2)
    # Matrix of variance values
    Sigma0 = np.zeros((q0, q0), dtype = float)
    for i in range(0, q0):
        for j in range(0, q0):
            Corr0 = CrossCorr(Codes[i], Codes[j], BTest[i], BTest[j])
            NCorr0 = Corr0/(np.sqrt(ACorrMAx[i]*ACorrMAx[j]))
            Sigma0[i][j] = np.var(NCorr0) + ((np.average(NCorr0))**2)
        # Eliminate variance of the auto correlation
        Sigma0[i][i] = 0.0
    # Average over all combinations
    S0 = np.mean(Sigma0)
    # get the vector of SIR to produce the final BER with the q-func
    # The value 0.2711155 correspond to calibration with the Fathallah simulation
    SIR = np.zeros((q0-1), dtype = float)
    BER = np.zeros((q0-1), dtype = float)
    for i in range(0, (q0-1)):
        SIR[i] = (0.2711155)/((i+1)*S0)
        BER[i] = 0.5*special.erfc(np.sqrt((SIR[i])/2))
    return BER