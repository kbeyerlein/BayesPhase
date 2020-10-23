#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 09:49:26 2018

@author: beyerlein
"""
import numpy as np
import matplotlib.pyplot as plt
import os

import dataLoaders as loader

L = 51095.755  # detector distance in pixels (found by Luca)
L = 51390  # value refined by hand 51240 
V = 300000      # accelerating voltage
lam = np.sqrt(150.4/V)


def getDetDist(imgName):
    if (imgName[0]=='c'):
        return L/3.985  # correction param found by hand: 3.965
    else:
        return L

def pixToQ(r, detDist):
    th = np.arctan2(r,detDist)/2.0
    return 2*np.sin(th)/lam

def genSpotPowder(imgName, dq, filetype="hit"):
    if (filetype == "hit"):
        r,phi,x,y,I,circ = loader.loadHitFile(imgName)
    else:
        print ("Filetype not supported.")
        return 0
    dd= getDetDist(imgName)
    q = pixToQ(r, dd)
    qmax = pixToQ(np.max(r)+10.0, dd)
    
    nbins = np.int(qmax/dq)
    
    Is,bins = np.histogram(q, bins=nbins, weights=I, range=(0,qmax))
    
    return bins[0:-1], Is

def getSiQs():
    Si_hkl, Si_d, Si_twoTheta = loader.loadHKLFile("Si_hkl.lst", path = os.getenv('SI_PATH', './'))
    return 1.0/Si_d

def getPhasePowder(filename):
    if ("Listing" in filename):
        temp, h, k, l, d, dstar, \
        temp, temp, temp, temp, I = loader.loadListingFile(filename)
    elif("SF" in filename):
        x,x,x,d,x,x,x,x,I,x = loader.loadHKLFile(filename)
        dstar = 1/d
    else:
        print ("Cannot load filetype.")
        return 0
    
    return dstar, I

def removeSiPeaks(q,I, dPeak):
    Si_qs = getSiQs()
    dq = q[1]-q[0]
    n = len(I)
    J = np.copy(I)
    for Si_q in Si_qs:
        qPMin = Si_q-dPeak
        qPMax = Si_q+dPeak
        iPMin = np.int((qPMin-q.min())/dq)
        iPMax = np.int((qPMax-q.min())/dq)
        if (iPMin<n):
            if (iPMax>(n-1)):
                J[iPMin:n-1]= 0
            else:
                J[iPMin:iPMax]= 0
    return J

def getMaxVal(d):
    maxes = [np.max(x) for x in d.values()]
    return np.max(maxes)
def getMinVal(d):
    mins = [np.min(x) for x in d.values()]
    return np.min(mins)

def sumPatts(qs, Is, dq):
    maxq = getMaxVal(qs)
    minq = getMinVal(qs)
    nq = int((maxq-minq)/dq) +1
    newQs = np.linspace(minq, maxq, nq)
    newIs = np.zeros(len(newQs))
    
    for x in qs:
        for i,I in enumerate(Is[x]):
            if (I != 0):
                j = int((qs[x][i]-minq)/dq)
                newIs[j]+=I
    return newQs, newIs