#!/usr/bin/env python
# -*- coding: utf-8 -*-

#########################################################################################################
# functions to calculate the cross correlation in time domain, in Spectral.csv has the power response
# of correlational parameters. 
# pandas is the library ro read .csv data (comma separate dataset, not point and comma), each row has the
# specral response of input pulse with different lambda separation
# matplotlib contain the basis to plot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from Bin import * 

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
    # CCVect = VectorCC[0:(Dim1-NB1[L0-1])]
    return VectorCC

##########################################################################################################
# Function Stats to produce the final probability error based in Q function. Gen is the generator function
# q0 correspond to the maximum number of users. BTest is the default profile of blank spaces between 
# adjacent FBG
def Stats(Gen, q0, BTest):
    Codes = UCodes(Gen, q0)
    # Values to normalize the crosscorrelation values
    ACorrMAx = np.zeros(q0, dtype = float)
    for i in range(0, q0):
        ACorrMAx[i] = np.amax(CrossCorr(Codes[i], Codes[i], BTest, BTest))
    print ACorrMAx
    # Matrix of variance values
    Sigma0 = np.zeros((q0, q0), dtype = float)
    for i in range(0, q0):
        for j in range(0, q0):
            Corr0 = CrossCorr(Codes[i], Codes[j], BTest, BTest)
            NCorr0 = Corr0/(np.sqrt(ACorrMAx[i]*ACorrMAx[j]))
            Sigma0[i][j] = np.var(NCorr0)
        # Eliminate variance of the auto correlation
        Sigma0[i][i] = 0.0
    # Average over all combinations
    print Sigma0
    S0 = np.mean(Sigma0)
    return S0

##########################################################################################################
# Extract the values of response for FBG
Spectral0 = pd.read_csv("Spectral.csv")
# convert dataframe to numpy array, first column is time, second delta lambda = 0, third delta lambda = 1, ...
Specter = Spectral0.values
Ns = np.shape(Specter)[0]                                                                               # Number of data in symbol section
Nb = 50                                                                                                 # Number of data for delay sections
DeltaM = np.shape(Specter)[1]

# Fathallah case, generator used by original article with FBG parameters of Spectral.csv
Gen0 = np.array([13,9,12,11,14,10,16,20,17,18,15,19], dtype = np.int_)
q0 = 29
Codes = UCodes(Gen0, q0)
L0 = np.shape(Codes)[1]
A0 = Codes[0]
A1 = Codes[1]
print A0
print A1
#B0 = np.array([220, 230, 150, 258, 340, 456, 100, 300, 119, 234, 76])
B0 = np.zeros((L0 - 1), dtype = int) + Nb

S0 = Stats(Gen0, q0, B0)
print S0

CC00 = CrossCorr(A0, A0, B0, B0)
CC2 = np.amax(CC00)
CC00 = CC00/CC2
X0 = np.arange(np.shape(CC00)[0])

#plt.plot(X0, CC00)
#plt.show()