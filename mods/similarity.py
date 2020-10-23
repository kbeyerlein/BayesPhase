#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:36:02 2020

@author: beyerlein
"""

import numpy as np

def similarity(eps, sigma=0.01):
    #eps = (qp-qi)/qi
    return np.exp(-eps**2/(2*sigma**2))

### MAKING SIMILARITY MATRICES

# Calculate similarity of image spot Qr with best matching phaseQrs to get indication of which spots can match a given phase.
def calcSimilarity_imageQr_phase(imageQr, phaseQr, sigma=0.01):
    bestHKLIdx = np.zeros(len(imageQr)).astype(int)
    values = np.zeros(len(imageQr))
    for i,q in enumerate (imageQr):
        dels = np.absolute(phaseQr-q)
        m = np.min(dels)
        # only calculate similarity for best matches
        if (m<5*sigma):
            bestHKLIdx[i] = np.where(dels==m)[0][0].astype(int)
            values[i]=similarity(dels[bestHKLIdx[i]])
    return values, bestHKLIdx

def makeQrSimilaritySpotPhaseMatrix(imageQr, PhaseQrs, phaseNames, sigma = 0.01):
    matrix = []
    bestIdxs = []
    for phase in phaseNames:
        scores, bestHKLIdx = calcSimilarity_imageQr_phase(imageQr,PhaseQrs[phase], sigma = sigma)
        matrix.append(scores)
        bestIdxs.append(bestHKLIdx)
    return np.array(matrix), np.array(bestIdxs).astype(int)

def make1darray2d(a):
    return a.reshape((1,len(a)))

def calcVectorDistanceMatrix(a,b):
    # Try to handle case to convert 1d arrays into 2d
    if (a.ndim == 1):
        a = make1darray2d(a)
    if (b.ndim == 1):
        b = make1darray2d(b)
        
    dist = np.zeros((len(a), len(b)))
    for j,q in enumerate(b):
        d = a-q
        dist[:,j]=np.linalg.norm(d, axis=1)
    return dist

def calcVectorCosAngleMatrix(a,b):
    # Try to handle case to convert 1d arrays into 2d
    if (a.ndim == 1):
        a = make1darray2d(a)
    if (b.ndim == 1):
        b = make1darray2d(b)
        
    aDotb = np.matmul(a,b.T)
    aMag = np.linalg.norm(a, axis=1)
    bMag = np.linalg.norm(b, axis=1)
    
    ab = np.outer(aMag,bMag)
    return aDotb/ab

# Calculate similarity of image spot Qr with best matching phaseQrs to get indication of which spots can match a given phase.
def calcSimilarity_imageQr_phaseQr(imageQr, phaseQr, sigma=0.01):
    values = np.zeros((len(imageQr), len(phaseQr)))
    for i,q in enumerate (imageQr):
        dels = np.absolute(phaseQr-q)
        # only calculate similarity for best matches
        idxs = np.where(dels<5*sigma)
        values[i][idxs]= similarity(dels[idxs[0]])
    return values

# imgQ : Nx3 2d numpy array of spot q vectors
# phaseQ: Nx3 2d numpy array of phase q vectors
# phaseQr: N 1d numpy array of phase q vectors
def calc2PtSimilarity_imageQ_phaseQ(imgQ, phaseQ, phaseQr, d_sigma=0.01, th_sigma=0.17, goodieThreshold=0.75):
    
    #myList = []
    finalscores =[]
    
    # Get out of here, nothing to do
    if (len(imgQ) ==0):
        return np.array([]), np.array([])
    
    # Calculate all distances and angles between observed spots
    # twoPtDist: Nspots x Nspots array, containing distances, symmetric
    twoPtDist = calcVectorDistanceMatrix(imgQ, imgQ)

    cosAng = calcVectorCosAngleMatrix(imgQ, imgQ)
    

    # find best matching hkls for phase
    # hkl_scores: Nspots x Nhkls array containing similarities
    imgQr = np.linalg.norm(imgQ, axis=1)
    hkl_scores = calcSimilarity_imageQr_phaseQr(imgQr, phaseQr, sigma = d_sigma)

    #Calculate distance between potential hkls
    # totHKLD: For each spot, j, calc distance between good hkl, k, and good hkls,idxs2, for other spots, l. 
    for j in range(len(imgQ)):
        idxs = np.where(hkl_scores[j]>goodieThreshold)

        totHKLD=[]
        totHKLCosAng=[]
        for k in idxs[0]:
            hklD=[]
            hklCosAng=[]
            for l in range(len(imgQ)):
                if (j==l):
                    hklD.append(np.zeros(len(idxs[0])))
                    hklCosAng.append(np.ones(len(idxs[0])))
                else:
                    idxs2 = np.where(hkl_scores[l]>goodieThreshold)
                    
                    if (len(idxs2[0])>0):

                        hklD.append(calcVectorDistanceMatrix(phaseQ[idxs2], phaseQ[k]))

                        hklCosAng.append(calcVectorCosAngleMatrix(phaseQ[idxs2], phaseQ[k]))
                    else:
                        hklD.append([])
                        hklCosAng.append([])

            totHKLD.append(hklD)
            totHKLCosAng.append(hklCosAng)

        #Calculate distance simlarities
        twoPtScores=[]
        for k,hklidx in enumerate(idxs[0]):
            delList=[]
            angDelList=[]
            for l,dSpot in enumerate(twoPtDist[j]):
                if (dSpot>0):
                    
                    hklD = totHKLD[k][l]
                    
                    if (len(hklD)>0):
                        temp = np.zeros(len(hklD))
                        # Compare distance between hkls and image spots
                        dels = np.absolute(np.array(hklD)-dSpot)
                        idxs2 = np.where(dels<5*d_sigma)

                        # Compare difference between hklCosAngles and imageCosAngles
                        hklCosAng = totHKLCosAng[k][l]
                        angDels = np.absolute(np.array(hklCosAng)- cosAng[j][l])

                        if (len(idxs2[0])>0):

                            temp[idxs2[0]] = similarity(dels[idxs2], sigma=d_sigma)*similarity(angDels[idxs2], sigma=th_sigma)

                            if (np.max(temp) > 0):
                                # collect best matching hkl info
                                bestmatch=np.where(temp==np.max(temp))[0]
                                delList.append((np.array(hklD)-dSpot)[bestmatch[0]])
                                angDelList.append((np.array(hklCosAng)- cosAng[j][l])[bestmatch[0]])
                            else:
                                # for indices that do not have a reasonable match, use -99 value
                                delList.append(-99)
                                angDelList.append(-99)
                        else:
                            delList.append(-99)
                            angDelList.append(-99)
                    else:
                        delList.append(-99)
                        angDelList.append(-99)
                    
                else:
                    delList.append(-99)
                    angDelList.append(-99)
                    
            delList = np.array(delList)
            angDelList = np.array(angDelList)
            
            if (len(np.where(delList>-99.0)[0])>0):
                #calculate average distance and angle from set if more than one best matching hkl
                avgDel=np.average(delList[np.where(delList>-99.)])
                avgAngDel=np.average(angDelList[np.where(angDelList>-99.)])
                
                # calculate similarity score as product of average spot distance similarity and abvg angle distance similarity
                twoPtScores.append(similarity(avgDel, sigma=d_sigma)*similarity(avgAngDel, sigma=th_sigma))
            else:
                twoPtScores.append(0.0)
            
        # Take best score from hkl list
        if (len(twoPtScores)>0):
            finalscores.append(np.max(np.array(twoPtScores)))
        else:
            finalscores.append(0.0)
        #myList.append(np.unique(np.where(np.array(twoPtScores)>(sim_min**2))[1]))
    
    # output spotwise score and array index of goodies
    # best match index currently not supported.
    return np.array(finalscores), np.array(-99*np.ones(len(imgQ)))

def make2PtSimilaritySpotPhaseMatrix(imageQ, phaseQs, phaseQrs, phaseNames, d_sigma = 0.01, th_sigma = 0.17):
    matrix = []
    for phase in phaseNames:
        if (phase == "Si-I"):
            matrix.append(np.zeros(len(imageQ)))
        else:
            scores, idx = calc2PtSimilarity_imageQ_phaseQ(imageQ, phaseQs[phase], phaseQrs[phase], d_sigma = d_sigma, th_sigma = th_sigma)
            matrix.append(scores)
            
    return np.array(matrix)

### MANIPULATING SIMILARITY MATRICES

def getMaxSimilarity(similarityMatrix, thresh = 0.5):
    i = np.where(similarityMatrix>thresh)
    return np.max(similarityMatrix[i])

def getMinSimilarity(similarityMatrix, thresh = 0.5):
    i = np.where(similarityMatrix>thresh)
    return np.min(similarityMatrix[i])

def getAvgSimilarity(similarityMatrix, thresh = 0.5):
    i = np.where(similarityMatrix>thresh)
    if (len(i[0]) == 0):
        return np.nan
    else:
        return np.average(similarityMatrix[i])

def getStdDevSimilarity(similarityMatrix, thresh = 0.5):
    i = np.where(similarityMatrix>thresh)
    if (len(i[0]) == 0):
        return 0
    else:
        return np.std(similarityMatrix[i])

def getPhaseAvgSimilarity(similarityMatrix, thresh = 0.5):
    return np.array([getAvgSimilarity(x,thresh=thresh) for x in similarityMatrix])

def getPhaseStdDevSimilarity(similarityMatrix, thresh = 0.5):
    return np.array([getStdDevSimilarity(x,thresh=thresh) for x in similarityMatrix])


def findGoodies(similarityMatrix, goodThresh= 0.75):
    return np.where(similarityMatrix > goodThresh)

def findBesties(similarityMatrix, goodThresh = 0.75 ):
    return [np.argmax(scores) if np.max(scores) > goodThresh else -1 for scores in similarityMatrix.T ]
