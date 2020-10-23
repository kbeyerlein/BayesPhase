#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 16:03:45 2018

@author: beyerlein
"""
import numpy as np
import dataLoaders as dl

def convertEtoLam(E):
    return 12398.42/E

# Create Unit cell vectors from unit cell constants
def get_abcVectors(a,b,c,al,be,ga):
    al_r  = np.radians(al)
    be_r = np.radians(be)
    ga_r = np.radians(ga)
    A = [a,0,0]
    B = [b*np.cos(ga_r), b*np.sin(ga_r), 0]
    temp = (np.cos(al_r)-np.cos(be_r)*np.cos(ga_r))/np.sin(ga_r)
    C = [c*np.cos(be_r), c*temp, c*np.sqrt(1.0-np.cos(be_r)**2 - temp**2)]
    return A,B,C

# Create reciprocal lattice vectors from unit cell vectors
def get_recipVectors (a,b,c):
    V = np.dot(a,np.cross(b,c))
    aStar = np.cross(b,c)/V
    bStar = np.cross(c,a)/V
    cStar = np.cross(a,b)/V
    return aStar,bStar,cStar

def getBraggTheta(h,k,l,aStar,bStar,cStar, lam):
    dStar = np.linalg.norm(h*aStar+k*bStar+l*cStar)
    return (np.arcsin(dStar*lam/2.0))


def SG1(h,k,l):
    return True

def SG14(h,k,l):
    if (k==0 and l%2 ==1):
        return False
    if (h==0 and k%2 ==1 and l == 0):
        return False
    return True


def SG64(h,k,l):
    if ((h+k)%2 ==1):
        return False
    if (h==0 and (k%2 ==1)):
        return False
    if (k==0 and (h%2 ==1 or l%2 ==1)):
        return False
    if (l==0 and (h%2 ==1 or k%2 ==1)):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%2 ==1):
        return False
    return True

def SG72(h,k,l):
    if ((h+k+l)%2 ==1):
        return False
    if (h==0 and (k%2 ==1 or l%2 ==1)):
        return False
    if (k==0 and (h%2 ==1 or l%2 ==1)):
        return False
    if (l==0 and ((h+k)%2 ==1)):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%2 ==1):
        return False
    return True

def SG74(h,k,l):
    if ((h+k+l)%2 ==1):
        return False
    if (h==0 and (k+l)%2 ==1):
        return False
    if (k==0 and (h+l)%2 ==1):
        return False
    if (l==0 and (h%2 ==1 or k%2 ==1)):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%2 ==1):
        return False
    return True

def SG88(h,k,l):
    if ((h+k+l)%2 ==1):
        return False
    if (h==0 and (k+l)%2 ==1):
        return False
    if (k==0 and (h+l)%2 ==1):
        return False
    if (l==0 and (h%2 ==1 or k%2 ==1)):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%4 > 0):
        return False
    return True

def SG92(h,k,l):
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%4 > 0):
        return False
    return True

def SG93(h,k,l):
    if (h==0 and k==0 and l%2 ==1):
        return False
    return True

def SG114(h,k,l):
    if ( h==k and l%2 == 1):
        return False
    if ( h%2==1 and k==0 and l==0):
        return False
    if ( h==0 and k%2==1 and l==0):
        return False
    return True

def SG140(h,k,l):
    if ((h+k+l)%2 == 1):
        return False
    if (h==0 and (k%2 ==1 or l%2 ==1)):
        return False
    if (k==0 and (h%2 ==1 or l%2 ==1)):
        return False
    if (h==k and l%2 ==1 ):
        return False
    return True

def SG141(h,k,l):
    if ((h+k+l)%2 == 1):
        return False
    if (h==0 and (k+l)%2 ==1):
        return False
    if (k==0 and (h+l)%2 ==1):
        return False
    if (l==0 and (h%2 ==1 or k%2 ==1)):
        return False
    if (h==k and l!=0 and (2*h+l)%4>0):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%4 > 0):
        return False
    if (h==k and l==0 and h%2==1):
        return False
    return True

def SG194(h,k,l):
    if (h==k and l%2 ==1):
        return False
    if (h==0 and k==0 and l%2==1):
        return False
    return True

def SG206(h,k,l):
    if ((h+k+l)%2 == 1):
        return False
    if (h==0 and (k%2 ==1 or l%2 ==1)):
        return False
    if (k==0 and (h%2 ==1 or l%2 ==1)):
        return False
    if (l==0 and (h%2 ==1 or k%2 ==1)):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%2 ==1):
        return False
    return True

def SG225(h,k,l):
    if ((h+k)%2 ==1 or (h+l)%2 ==1 or (k+l)%2==1):
        return False
    if (h==0 and (k%2==1 or l%2 ==1)):
        return False
    if (k==0 and (h%2==1 or l%2 ==1)):
        return False
    if (l==0 and (h%2==1 or k%2 ==1)):
        return False
    if (h==k and (h+l)%2==1):
        return False
    if (h==l and (h+k)%2==1):
        return False
    if (k==l and (h+k)%2==1):
        return False
    if (k==0 and l==0 and h%2 ==1):
        return False
    if (h==0 and l==0 and k%2 ==1):
        return False
    if (h==0 and k==0 and l%2 ==1):
        return False
    if (h==0 and k==l and k%2==1):
        return False
    if (k==0 and h==l and h%2==1):
        return False
    if (l==0 and h==k and h%2==1):
        return False
    return True

def SG227(h,k,l):
    if ((h+k)%2 ==1 or (h+l)%2 ==1 or (k+l)%2==1):
        return False
    if (h==0 and (k%2==1 or l%2 ==1 or (k+l)%4 > 0)):
        return False
    if (k==0 and (h%2==1 or l%2 ==1 or (h+l)%4 > 0)):
        return False
    if (l==0 and (h%2==1 or k%2 ==1 or (h+k)%4 > 0)):
        return False
    if (h==k and (h+l)%2==1):
        return False
    if (h==l and (h+k)%2==1):
        return False
    if (k==l and (h+k)%2==1):
        return False
    if (k==0 and l==0 and h%4>0):
        return False
    if (h==0 and l==0 and k%4>0):
        return False
    if (h==0 and k==0 and l%4>0):
        return False
    if (h==0 and k==l and k%2==1):
        return False
    if (k==0 and h==l and h%2==1):
        return False
    if (l==0 and h==k and h%2==1):
        return False
    return True
    
    


def noSG(h,k,l):
    print ("Invalid Space Group")
    return False

def CheckHKL(sg, h,k,l):
    switcher = {
            1:SG1,
            14:SG14,
            64:SG64,
            72:SG72,
            74:SG74,
            88:SG88,
            92:SG92,
            93:SG93,
            96:SG92,
            114:SG114,
            140:SG140,
            141:SG141,
            148:SG1,
            191:SG1,
            194:SG194,
            206:SG206,
            225:SG225,
            227:SG227
    }
    selRules = switcher.get(sg, noSG)
    return selRules(h,k,l)

def MakeSpotList(SGn, a, b,c, qmax):
    q= []
    hkl = []
    aStar,bStar,cStar = get_recipVectors(a,b,c)
    hmax = np.int(qmax/np.linalg.norm(aStar)) + 2
    kmax = np.int(qmax/np.linalg.norm(bStar)) + 2
    lmax = np.int(qmax/np.linalg.norm(cStar)) + 2
    for h in np.arange(hmax,-hmax, -1):
        for k in np.arange(kmax,-kmax, -1):
            for l in np.arange(lmax,-lmax, -1):
                if (CheckHKL(SGn, h,k,l)):
                    tempq = np.linalg.norm(h*aStar+k*bStar+l*cStar)
                    if (tempq<=qmax and tempq>0):
                        q.append(h*aStar+k*bStar+l*cStar)
                        hkl.append(np.str(h)+np.str(k)+np.str(l))
    return hkl, q