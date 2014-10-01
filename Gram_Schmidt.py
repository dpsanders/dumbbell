# -*- coding: utf-8 -*-
"""
Created on Thu May 22 15:33:46 2014

@author: rafael
"""

import numpy as np
from numpy import array

la = np.linalg


def proy(uu, vv):
    l = len(uu)
    if uu.any() == 0.:
        return np.zeros(l)
    
    else:
        coef = 1.*np.dot(uu,vv)/(la.norm(uu)**2)
        return coef*uu


def GS(vecs):
    N = len(vecs) #número de vectores
    dim = len(vecs[0]) #dimensiones de los vectores
    ovecs = []
    onvecs = []
    
    for k in range(N):
        proyections = np.zeros(dim)
        for u in ovecs:
            proyections += proy(u, vecs[k])
        #print proyections
        ovecs.append( vecs[k] - proyections )
        #print vecs[k]
        #print vecs[k] - proyections
        onvecs.append(ovecs[k]/(la.norm(ovecs[k])))
    return array(ovecs), array(onvecs)  # el primer resultado son los vectores ortogonales y el segundo los ortoNORMales
