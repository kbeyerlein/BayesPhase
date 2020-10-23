#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:31:58 2018

@author: beyerlein
"""

import numpy as np
import matplotlib.pyplot as plt

import dataLoaders as loader
import powderAnalysis as pA



def getSpotImage(imgName, filetype="hit", path=""):
    if (filetype == "hit"):
        r,phi,x,y,I,circ = loader.loadHitFile(imgName, path=path)
        return x,y,I
    elif (filetype =="spot.Q"):
        i,qx,qy,qz,qr,qphi,size,delta,I,dia,circ,mate = loader.loadSpotQFile(imgName, path=path)
        return qx/(2*np.pi),qy/(2*np.pi),qz/(2*np.pi),I
    elif (filetype =="Q"):
        i,qx,qy,qz,qr,qphi,size,delta,I,dia,circ,mate = loader.loadQFile(imgName, path=path)
        return qx/(2*np.pi),qy/(2*np.pi),qz/(2*np.pi),I
    else:
        print ("Filetype not supported.")
        return 0
    # Convert to qx, qy
    return x,y,I

def removeSiSpots(imgName, x,y,z, I,dq=0.01, convertQ = True):
    r = np.sqrt(x**2+y**2+z**2)
    if (convertQ):
        dd = pA.getDetDist(imgName)
        qs = pA.pixToQ(r,dd)
    else:
        qs = r
    
    Si_qs = pA.getSiQs()
    
    diffs = [np.absolute(Si_qs-q) for q in qs]
    
    newX=[]
    newY=[]
    newZ=[]
    newI=[]
    for i,d in enumerate(diffs):
        if (np.min(d)>dq):
            newX.append(x[i])
            newY.append(y[i])
            newZ.append(z[i])
            newI.append(I[i])
            
    return np.array(newX), np.array(newY), np.array(newZ), np.array(newI)