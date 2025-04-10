#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:10:35 2019

@author: gu32kij
"""
# Third party modules
from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt


class Component(object):
    """Describes a Sonata Component Object
    
    Attributes
    ----------
    name : str
        name of the component
        
    cosy : gp_Ax2
        Describes a right-handed coordinate system in 3D space. It is part of 
        the parent class 'Component'. It is an instance from the opencascade 
        gp_Ax2 class
        
    """

    __slots__ = ("name", "Ax2")

    def __init__(self, name="NONAME", *args, **kwargs):
        self.name = name
        self.Ax2 = gp_Ax2(*args)
        # self.units = {'mass': 'g', 'length': 'mm', 'force' : 'N'}
        # TODO: How to define Units in the global SONATA context

    def __repr__(self):
        """__repr__ is the built-in function used to compute the "official" 
        string reputation of an object, """
        return "Component: " + self.name

if __name__ == "__main__":
    C = Component(gp_Pnt(0, 10, 0), gp_Dir(1, 0, 0), name="TestComponent")
