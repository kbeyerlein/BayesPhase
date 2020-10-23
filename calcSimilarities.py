#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:20:14 2020

@author: beyerlein
"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

# path to modules used to calculate similarity
# modPath = "/gpfs/cfel/cxi/scratch/user/gelisio/shared/calcSimilarities/"
modPath = "/Users/beyerlein/Work/LaserSiHighPressurePhases/Analysis/"
sys.path.append(modPath)

import ImageAnalysis as iA
import powderAnalysis as pA
import diffraction as diff
import dataLoaders as dL
import similarity as sim

#print ('Number of arguments:', len(sys.argv), 'arguments.')
#print ('Argument List:', str(sys.argv))

if (len(sys.argv) != 4 ):
    print ("Usage: ")
    print ("  python calcSimilarities.py <imageFile> <phaseFile> <outputFile>")
    sys.exit(0)

imgFile = str(sys.argv[1])
phaseFile = str(sys.argv[2])
outputFile = str(sys.argv[3])

# Load image info from .Q file
x,y, z, imgI = iA.getSpotImage(imgFile, "Q")
print ("Loaded Image: ", imgFile)

# Remove spots overlapping with Si-I
l1 = len(x)
x,y, z, imgI = iA.removeSiSpots(imgFile, x,y, z, imgI,  convertQ=False)
imgQ = np.array([x,y,z]).T
imgQr = np.linalg.norm(imgQ, axis=1)
    # remove hkl =000
i = np.where(imgQr>0)
imgQ, imgQr, imgI = [imgQ[i], imgQr[i], imgI[i]]
l2 = len(imgQr)
print ("Spots in image: ", l1)
print ("Spots Remaining after removing Si-I: ", l2 )
print ("Fraction of Spots Remaining: ", l2/l1)


# Load phase info from hkl files
Si_Qr = {}
Si_Q = {}
Si_hkl = {}
Si_SF = {}
phaseNames = []

with open(phaseFile, 'r') as f:
    for line in f.readlines():
        l = line.strip().split(' ')
        name = l[0]
        phaseNames.append(name)
        hklFile = l[1]
        Si_hkl[name],Si_Q[name],Si_Qr[name],Si_SF[name] = dL.loadHKLH5File(hklFile, path="")
        # remove 0,0,0 spot
        i = np.where(Si_Qr[name]>0)
        Si_hkl[name], Si_Q[name], Si_Qr[name], Si_SF[name] = [Si_hkl[name][i], Si_Q[name][i], Si_Qr[name][i], Si_SF[name][i]]
        print ("Loaded phase: ", name)
       
NPhases = len(phaseNames)


# calculate Qr Similarity Matrix
d_sigma = 0.01
th_sigma = 0.07
QrSimilarity, bestQrHKLIdx = sim.makeQrSimilaritySpotPhaseMatrix(imgQr, Si_Qr, phaseNames, sigma = d_sigma)
print ("Calculated Qr Spot Phase Similarity using sigma = ", d_sigma)

# calculate Vector Q similarity Matrix
twoPtSimilarity = sim.make2PtSimilaritySpotPhaseMatrix(imgQ, Si_Q, Si_Qr, phaseNames, d_sigma = d_sigma, th_sigma = th_sigma)
print ("Calculated Q Spot Phase Similarity using sigmas = ", d_sigma, th_sigma)


np.savez(outputFile, phaseNames = phaseNames, 
         d_sigma = d_sigma, th_sigma = th_sigma, 
         QrSimilarity = QrSimilarity, bestQrHKLIdx = bestQrHKLIdx, 
         twoPtSimilarity = twoPtSimilarity)
print ("Wrote similarities to file: ", outputFile)

