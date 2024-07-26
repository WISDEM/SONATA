#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:22:29 2019

@author: gu32kij
"""
# Third party modules
import numpy as np
from OCC.Core.Geom import Geom_Plane
from OCC.Core.gp import gp_Dir, gp_Pnt
from scipy.interpolate import interp1d

# First party modules
from SONATA.cbm.topo.BSplineLst_utils import intersect_BSplineLst_with_plane
from SONATA.utl.blade_utl import check_uniformity


class interpBSplineLst(object):
    """
    
    
    """

    def __init__(self, BSplineLst, xgrid, xvalues):
        """
        Parameters
        --------
        BSplineLst : 
        xgrid : array
        xvalues : array
        """

        self.BSplineLst = BSplineLst
        if check_uniformity(xgrid, xvalues) == False:
            print("WARNING:\t The reference axis is not uniformly defined along x!")
        self._f_xint = interp1d(xgrid, xvalues, bounds_error=False, fill_value="extrapolate")

    def interpolate(self, grid):
        """
        interpolates the BSplineLst at a certain radial station by intersecting
        with a plane.
        
        Parameters
        --------
        grid : float
            nondimensional x coordinate in the blade ref. frame
        
        Retruns
        --------
        res, resCoords : (array,array)
            tuple of arrays. res contains the coordinates, while resCoords 
            contains the BSpline List parameter coordinates i and u.
            
        """

        if not isinstance(grid, np.ndarray):
            grid = np.atleast_1d(grid)
        res = []
        resCoords = []
        for gp in grid:
            plane = Geom_Plane(gp_Pnt(float(self._f_xint(gp)), 0, 0), gp_Dir(1, 0, 0))
            IntCoords, IntPnts = intersect_BSplineLst_with_plane(self.BSplineLst, plane)
            if len(IntPnts) > 1:
                print("WARNING:\t More than one intersection point was found!")

            try:
                res.append(IntPnts[0].Coord())
                resCoords.append(IntCoords[0])
            except:
                print("No Intersection Found")
                raise Exception

        return (np.asarray(res), np.asarray(resCoords))

if __name__ == "__main__":
    pass
