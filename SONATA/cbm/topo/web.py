# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:28:42 2017

@author: TPflumm
"""
from OCC.Core.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.Core.gp import gp_Pnt2d

from SONATA.cbm.topo.utils import point2d_list_to_TColgp_Array1OfPnt2d
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import ParaLst_from_BSplineLst, BSplineLst_from_ParaLst



#
# from OCC.Display.SimpleGui import init_display
# display, start_display, add_menu, add_function_to_menu = init_display()
# from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge2d

class Web(object):
      #Build Webs:
        #TODO: CHECK IF WEB DEFINITION INTERSECT EACH OTHER
        #TODO: SORT WEBS BY POS1 VALUES:
    def __init__(self, ID, Pos1, Pos2, curvature=0, SegmentLst=None):
        """

        ID:         web id (?)
        Pos1:       Pos of start point of web
        Pos2:       Pos of end point of web
        curvature:  value that defines the chordwise displacement of the ellipsis, i.e. the curvature of a web
        SegmentLst: List of segments

        """

        self.ID = ID
        self.Pos1 = Pos1
        self.Pos2 = Pos2
        
        if curvature is None:
            self.curvature = 0
        else: self.curvature = curvature
        
        self.Pos1_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos1, SegmentLst)    # p1 - start point at arc
        self.Pos2_Pnt2d = SegmentLst[0].get_Pnt2d(self.Pos2, SegmentLst)    # p2 - end point at arc

        BSplineLst_straight = [Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.Pos1_Pnt2d, self.Pos2_Pnt2d])).Curve()]


        if self.curvature == 0:  # if no curvature is assigned or curvature is assigned to be zero
            self.BSplineLst = BSplineLst_straight

        self.wr_nodes =[] 
        self.wl_nodes = []
        self.wr_cells =[] 
        self.wl_cells = []
        
        #self.BackBSplineLst = trim_BSplineLst(BSplineLst, S1, S2, start, end):
        #Intersect Segment0_Boundary_BSplineLst with self.BSpline
        #[self.IntPnts,self.IntPnts_Pnt2d] = intersect_BSplineLst_with_BSpline(self.Segment0_Boundary_BSplineLst,self.BSpline_Line)
        
    def __getstate__(self):
        """Return state values to be pickled."""
        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        self.Pos1_coordinates = (self.Pos1_Pnt2d.X(),self.Pos1_Pnt2d.Y())
        self.Pos2_coordinates = (self.Pos2_Pnt2d.X(),self.Pos2_Pnt2d.Y())
        return (self.ID, self.Pos1, self.Pos2, self.Pos1_coordinates, self.Pos2_coordinates, self.Para_BSplineLst, self.wr_nodes, self.wl_nodes)   
    
    
    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (self.ID, self.Pos1, self.Pos2, self.Pos1_coordinates, self.Pos2_coordinates, self.Para_BSplineLst, self.wr_nodes, self.wl_nodes)  = state
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)
        self.Pos1_Pnt2d = gp_Pnt2d(self.Pos1_coordinates[0],self.Pos1_coordinates[1])
        self.Pos2_Pnt2d = gp_Pnt2d(self.Pos2_coordinates[0],self.Pos2_coordinates[1])
        