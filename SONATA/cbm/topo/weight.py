# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 12:49:54 2017

@author: TPflumm
"""
# Core Library modules
import os

# Third party modules
import numpy as np
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,)
from OCC.Core.Geom2d import Geom2d_BezierCurve, Geom2d_Circle
from OCC.Core.gp import (gp_Ax2, gp_Ax2d, gp_Circ, gp_Circ2d, gp_Dir,
                         gp_Dir2d, gp_Pnt, gp_Pnt2d, gp_Vec2d,)
from OCC.Display.SimpleGui import init_display

# First party modules
from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d

if __name__ == "__main__":
    os.chdir("../../..")



def WeightStyle1_BezierPoint(P1, P2, R, W):
    Vec12 = gp_Vec2d(P1, P2)
    n12 = Vec12.GetNormal()
    n12.Normalize()
    P_tmp = P1.Translated(Vec12.Multiplied(0.5))
    P12 = P_tmp.Translated(n12.Multiplied(W * R))
    return P12


class Weight(object):
    def __init__(self, ID, Pnt2d, Diameter, MatID, weight_style=0):
        self.ID = ID
        self.Pnt2d = Pnt2d
        # self.Y = YPos
        self.D = Diameter
        self.MatID = MatID
        self.style = weight_style

        # Circ = gp_Circ2d(Ax2d,R)
        # Center = gp_Pnt2d(self.X,self.Y)
        Ax2d = gp_Ax2d(self.Pnt2d, gp_Dir2d(1, 0))
        self.Curve = Geom2d_Circle(Ax2d, self.D / 2)

    @property
    def wire(self):
        return self.__build_wire()

    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.ID, self.X, self.Y, self.D, self.MatID, self.style)

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.ID, self.X, self.Y, self.D, self.MatID, self.style) = state
        Center = gp_Pnt2d(self.X, self.Y)
        Ax2d = gp_Ax2d(Center, gp_Dir2d(1, 0))
        self.Curve = Geom2d_Circle(Ax2d, self.D / 2)

    def __build_wire(self):

        Center = gp_Pnt(self.Pnt2d.X(), self.Pnt2d.Y(), 0)
        Ax2 = gp_Ax2(Center, gp_Dir(0, 0, 1))
        Circ = gp_Circ(Ax2, self.D / 2)
        aedge = BRepBuilderAPI_MakeEdge(Circ).Edge()
        wire = BRepBuilderAPI_MakeWire(aedge)
        return wire.Wire()


# ======================================================
#       MAIN
# ======================================================
if __name__ == "__main__":
    # ==========================
    # DISPLAY CONFIG:
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.Context.SetDeviationAngle(0.000001)  # 0.001 default. Be careful to scale it to the problem.
    display.Context.SetDeviationCoefficient(0.000001)  # 0.001 default. Be careful to scale it to the problem.
    # show_coordinate_system(display) #CREATE AXIS SYSTEM for Visualization

    BW = Weight(1, 0.1, 0.1, 0.25, 1, 0)
    display.DisplayShape(BW.Curve, color="BLACK")

    display.set_bg_gradient_color(20, 6, 111, 200, 200, 200)

    display.View_Top()
    display.FitAll()

    start_display()
