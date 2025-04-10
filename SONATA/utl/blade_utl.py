#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:20:51 2019

@author: gu32kij
"""
# Third party modules
import numpy as np
from scipy.interpolate import interp1d

# PythonOCC Libraries
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.TopoDS import TopoDS_Vertex, TopoDS_Wire
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.Approx import Approx_ChordLength

# First party modules
from SONATA.cbm.topo.utils import (Array_to_PntLst, PntLst_to_npArray,
                                   lin_pln_intersect,)


class assert_isdone(object):
    '''
    raises an assertion error when IsDone() returns false, with the error
    specified in error_statement
    -> this is from the pythonocc-utils utility-may not use it?
    '''
    def __init__(self, to_check, error_statement):
        self.to_check = to_check
        self.error_statement = error_statement

    def __enter__(self, ):
        if self.to_check.IsDone():
            pass
        else:
            raise AssertionError(self.error_statement)

    def __exit__(self, type, value, traceback):
        pass

def interp_loads(loads, grid_loc):
    """
    Interpolates the loads at the given radial station (grid location)
    
    Parameters
    ----------
    loads : dict
        dictionary of the following keys and values, (default=None)
        for detailed information see the VABSConfig documentation or the 
        VABS user manual
        F : nparray([[grid, F1, F2, F3]]) 
        M : nparray([[grid, M1, M2, M3]]) 
        f : nparray([[grid, f1, f2, f2]])
        df : nparray([[grid, f1', f2', f3']])
        dm :  nparray([[grid, m1', m2', m3']])
        ddf : nparray([[grid, f1'', f2'', f3'']])
        ddm : nparray([[grid, m1'', m2'', m3'']])
        
    grid_loc : float
        location of interpolation
        
    Returns
    ----------
    sectional_load : dict
    dictionary of the following keys and values, (default=None)
        for detailed information see the VABSConfig documentation or the 
        VABS user manual
        F : nparray([F1, F2, F3]) 
        M : nparray([M1, M2, M3]) 
        f : nparray([f1, f2, f2])
        df : nparray([f1', f2', f3'])
        dm :  nparray([m1', m2', m3'])
        ddf : nparray([f1'', f2'', f3''])
        ddm : nparray([m1'', m2'', m3''])
        
    """

    d = {}
    for k, item in loads.items():
        fit = interp1d(item[:, 0], item[:, 1:], axis=0)
        d[k] = fit(grid_loc)
    return d


def interp_airfoil_position(airfoil_position, airfoils, grid_loc):
    """
    
    
    Parameters
    ----------
    airfoil_position: tuple
        ([0.3, 1.0], ['n0012', 'n0012'])
    airfoils: list
        [Airfoil: n0012, Airfoil: naca23012]
    """

    # TBD: Extrapolation
    if grid_loc in airfoil_position[0]:
        afname = airfoil_position[1][airfoil_position[0].index(grid_loc)]
        return next((x for x in airfoils if x.name == afname), None)

    # find closest value:
    min_idx = np.argmin([abs(x - grid_loc) for x in airfoil_position[0]])
    min_val = airfoil_position[0][min_idx]

    if grid_loc > min_val:
        iv_idx = (min_idx, min_idx + 1)
    else:
        iv_idx = (min_idx - 1, min_idx)

    iv_val = tuple(airfoil_position[0][iv_idx[0] : iv_idx[1] + 1])
    iv_af = tuple(airfoil_position[1][iv_idx[0] : iv_idx[1] + 1])
    k = (grid_loc - iv_val[0]) / (iv_val[1] - iv_val[0])

    # select af from airfoils
    af1 = next((x for x in airfoils if x.name == iv_af[0]), None)
    af2 = next((x for x in airfoils if x.name == iv_af[1]), None)

    if af1 == af2:
        return af1

    # return transformed airfoil
    return af1.transformed(af2, k)


def make_loft(elements, solid=False, ruled=False, tolerance=1e-6, 
              continuity=GeomAbs_C2, max_degree=8, check_compatibility=True, **kwargs):
    """
    A set of sections that are used to generate a surface with the 
        BRepOffsetAPI_ThruSections function from OCC
        

    Parameters
    ----------
    elements : list
        list of OCC.TopoDS_Wire or TopoDS_Vertex
        A set of sections that are used to generate a surface with the 
        BRepOffsetAPI_ThruSections function from OCC.
    solid : bool, optional
        solid or surface. The default is False.
    ruled : bool, optional
        linear are nurbs surface. The default is False.
    tolerance : float, optional
        tolerance. The default is 1e-6.
    continuity : int, optional
        DESCRIPTION. The default is 4 (GeomAbs_C2).
    max_degree : int, optional
             The order of the fitted NURBS surface. The default is 8.
    check_compatibility : TYPE, optional
        DESCRIPTION. The default is True.
    **kwargs : TYPE
        DESCRIPTION.


    Raises
    ------
    TypeError
        DESCRIPTION.

    Returns
    -------
    loft : TopoDS_Shape
        surface of the ThruSections Loft

    """
    
    generator = BRepOffsetAPI_ThruSections(solid, ruled, tolerance)
    generator.SetMaxDegree(max_degree)
    generator.SetParType(Approx_ChordLength)
    for i in elements:
        if isinstance(i, TopoDS_Wire):
            generator.AddWire(i)
        elif isinstance(i, TopoDS_Vertex):
            generator.AddVertex(i)
        else:
            raise TypeError("elements is a list of TopoDS_Wire or TopoDS_Vertex, found a %s fool" % i.__class__)

    generator.CheckCompatibility(check_compatibility)
    generator.SetContinuity(continuity)
    generator.Build()
    
    with assert_isdone(generator, 'failed lofting'):
        loft = generator.Shape() 
        return loft
        

def check_uniformity(grid, values, tol=1e-6):
    """
    Checks the uniformity of the values along the grid by calculating the 
    gradient and checking if it's constant with respect to a giving tolererance
    
    Parameters
    ----------
    grid : array
    values : array
    tol : float, optional
    
    Returns
    ----------
    bool
    """
    grad = np.gradient(values, grid)
    mean = np.mean(grad)
    return all(mean - tol < x < mean + tol for x in grad)


def array_pln_intersect(array, ax2):
    """
    intersects an array of connecting points with the yz plane of the 
    ax2 coordinate system.
    
    Parameters:
        ax2 : gp_Ax2
            right handed coordinate system
        array : 
    
    """
    factors = []
    coords = []
    for i in range(len(array) - 1):
        coord = []
        factor = []
        for p1, p2 in zip(array[i], array[i + 1]):
            # print(p1,p2)
            pnt, lmb = lin_pln_intersect(ax2.XDirection().Coord(), ax2.Location().Coord(), p1, p2)
            coord.append(pnt)
            factor.append(lmb)
        coords.append(np.asarray(coord))
        factors.append(np.asarray(factor))

    # ==== extrapolate ===
    coords = np.asarray(coords)
    factors = np.asarray(factors)

    result = []
    # iterate first over points j and than over airfoils i
    for j, lmbs in enumerate(factors.swapaxes(0, 1)):
        found = False
        for i, l in enumerate(lmbs):
            if 0 <= l <= 1:
                # print(i,j,l,res2[i,j])
                result.append(coords[i, j])
                found = True
                break

        if found == False:
            # print(j, lmbs, lmbs[-1]>0,lmbs[0]<0)
            # try last
            if lmbs[-1] > 0:
                i = len(lmbs) - 1
                result.append(coords[i, j])
            # try first
            elif lmbs[0] < 0:
                i = 0
                result.append(coords[i, j])
            else:
                result.append(np.array([np.nan, np.nan, np.nan]))

    return np.asarray(result)
