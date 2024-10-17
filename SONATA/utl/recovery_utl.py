#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13, 2020

@author: RF
"""

import numpy as np



def vabs_recovery_utl(F, M):
    """

    """

    loads = {
        "u": np.array([[0, 1, 1, 1], [1, 1, 1, 1]]),  # displacement dof
        "Ci1": np.array([[0, 1, 0, 0], [1, 1, 0, 0]]),  # direction cosine matrix
        "Ci2": np.array([[0, 0, 1, 0], [1, 0, 1, 0]]),  # direction cosine matrix
        "Ci3": np.array([[0, 0, 0, 1], [1, 0, 0, 1]]),  # direction cosine matrix
        "F": np.array(
            [[0.0, F[0], F[1], F[2]],
             [1.0, F[0], F[1], F[2]]]),  # forces, N (F1: axial force; F2,F3: sectional transverse shear forces)
        "M": np.array(
            [[0.0, M[0], M[1], M[2]],
             [1.0, M[0], M[1], M[2]]]), # moments, Nm (M1: torsion; M2: flap; M3: lag)
        "f": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # distributed forces, N/m (including both applied forces and inertial forces)
        "df": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # first derivative
        "ddf": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # second derivative
        "dddf": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # third derivative
        "m": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # distributed moments, Nm/m (including both applied forces and inertial moments)
        "dm": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # first derivative
        "ddm": np.array([[0, 0, 0, 0], [1, 0, 0, 0]]),  # second derivative
        "dddm": np.array([[0, 0, 0, 0], [1, 0, 0, 0]])}  # third derivative

    return loads

def anba_recovery_utl(F, M):
    """
    Uses loads given in sonata coord system input and converts them to anba loads
    """

    loads = {
        "F": np.array(
            [[0, F[0], F[1], F[2]],
             [1, F[0], F[1], F[2]]]),
        # forces, N (F1: shear force in x-direction; F2: shear force in y -direction; F3: axial force)
        "M": np.array(
            [[0, M[0], M[1], M[2]],
             [1, M[0], M[1], M[2]]])}  # moments, Nm (M1: bending moment around x; M2: bending moment around y; M3: torsional moment)

    return loads