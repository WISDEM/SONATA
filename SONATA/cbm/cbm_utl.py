#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:27:46 2019

@author: gu32kij
"""
# Third party modules
import numpy as np


def trsf_sixbysix(M, T):
    """
    Transform six-by-six compliance/stiffness matrix. 
    change of reference frame in engineering (or Voigt) notation.
    
    Parameters
    ----------
    M : np.ndarray
        6x6 Siffness or Mass Matrix
    T : np.ndarray
        Transformation Matrix
        
    Returns
    ----------
    res : np.ndarray
        Transformed 6x6 matrix
    """

    TS_1 = np.dot(np.dot(T.T, M[0:3, 0:3]), T)
    TS_2 = np.dot(np.dot(T.T, M[3:6, 0:3]), T)
    TS_3 = np.dot(np.dot(T.T, M[0:3, 3:6]), T)
    TS_4 = np.dot(np.dot(T.T, M[3:6, 3:6]), T)

    tmp_1 = np.vstack((TS_1, TS_2))
    tmp_2 = np.vstack((TS_3, TS_4))
    res = np.hstack((tmp_1, tmp_2))
    return res
