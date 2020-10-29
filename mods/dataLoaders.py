#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 09:41:55 2018

@author: beyerlein
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py

def loadFile(f, nSkip, nSkipEnd = 0):
    
    with open(f, "r") as dataFile:
        data = dataFile.readlines()
        # labels = data[0].split()
        data = [map(float, \
                    data[i].replace("[","").replace("]","").replace("%","").replace("\t-\n","\t0\n").split()\
                    ) for i in range(nSkip,len(data)-nSkipEnd)]

    return np.array(list(zip(*data)))

def composeFilename(path, fileName, suffix):
    if (len(path)>0):
        if (path[-1]!="/"):
            path = path + "/"
        fileName = path + fileName
    if (len(suffix)>0):
        if (fileName[-len(suffix):] != suffix):
            fileName = fileName+suffix
    return fileName
    

def loadHKLFile(file, path=""):
    filename = composeFilename(path, file, "")
    return loadFile(filename,1)

def loadHitFile(imgName, path = "../Luca-analysis/FirstAnalysis/"):
    filename = path+imgName+".hit"
    return loadFile(filename,3)

def loadSpotQFile(imgName, path = "../Luca-analysis/Analysis-Oct2018/"):
    filename = path+"/"+imgName+"/"+imgName+".spot.Q"
    return loadFile(filename, 5)

def loadQFile(imgName, path = "../Luca-analysis/Analysis-Oct2018/"):
    filename = composeFilename(path, imgName, ".Q")
    return loadFile(filename, 1)

def loadListingFile(filename):
    if ("bt8" in filename):
        return loadFile(filename, 63, 1)
    else:
        return loadFile(filename, 78, 1)

def loadCell(file,path="../Ken Beyerlein diffraction analysis/New Si-phases/"):
    filename = path+file
    params = loadFile(filename,1)
    cell = {}
    cell['a'] = params[0][0]
    cell['b'] = params[1][0]
    cell['c'] = params[2][0]
    cell['al'] = params[3][0]
    cell['be'] = params[4][0]
    cell['ga'] = params[5][0]
    cell['sg'] = int(params[6][0])
    return cell

def loadHKLH5File(file, path = "../hklLists/"):
    f = h5py.File( composeFilename(path,file, ".HKL.h5" ))
    K = f['Q'][()] / (2*np.pi)
    Kr = f['Qr'][()] / (2*np.pi)
    hkl = f['hkl'][()]
    SF = f['SF'][()]
    f.close()
    return hkl, K, Kr, SF
    
