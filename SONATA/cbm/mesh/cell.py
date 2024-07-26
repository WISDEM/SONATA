# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 13:39:55 2017

@author: TPflumm
"""
# Core Library modules
import operator

# Third party modules
import numpy as np
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeWire,)
from OCC.Core.Geom2dAPI import (Geom2dAPI_PointsToBSpline,
                                Geom2dAPI_ProjectPointOnCurve,)
from OCC.Core.gp import gp_Pnt2d, gp_Vec2d

# First party modules
from SONATA.cbm.topo.utils import (PolygonArea, calc_angle_between,
                                   point2d_list_to_TColgp_Array1OfPnt2d,)


class Cell(object):
    __slots__ = ("id", "nodes", "theta_1", "theta_3", "MatID", "structured", "interior_nodes", "strain", "strainM", "stress", "stressM", "sf", "failure_mode")
    class_counter = 1

    def __init__(self, nodeLst):  # int
        self.id = self.__class__.class_counter
        self.__class__.class_counter += 1
        self.nodes = nodeLst  # [node,node,node,nodes]      !!!counterclockwise direction!!!
        # self.face  = []                     #[rear,top,front,bottom]        !!!counterclockwise direction!!!
        # self.wire  = self.build_wire()      #TopoDs_wire
        # self.wire  = None     #TopoDs_wire
        # self.neighbours = []                #-[Cell_ID,CELL_ID... ]
        self.theta_1 = [0] * 9  # Ply coordinate system is formed by rotating the global coordinate system in the right-hand sense about the amount 0<Theta_1<260.
        # Theta_1[0:9] is a list storing nine real numbers for the layer plane angles at the nodes of ths element. For simplification, if the
        # ply orinetation can be considered as uniform this element. Theta_1[0] stores the layer plane angles and Theta_1[1] = 540, and all the
        # rest can be zeroes or other real numbers because they do not enter the calculation. If the elements'''
        self.theta_3 = None  # The Ply coordiate system is rotated about y3 in the right hand sense by the amount -90<Theta_3<90 to for the material system.

        self.MatID = None  # material id,
        self.structured = True
        self.interior_nodes = []
        # Element quality critiria
        # self.minimum_edge = None

    #        self.maximum_edge = None
    #        self.shape_quality = None
    #        self.minimum_jacobinan = None
    # AREA RATIO to neighbors

    # THE AVERAGE OF 3D strain and Stress Results at Gaussian Points within each element.
    #        self.strain = Strain()      #[psilon11,2epsilon12,2epsilon13,epsilon22,2epsilon24,epsilon33
    #        self.strainM = Strain()   #[epsilon11,2epsilon12,2epsilon13,epsilon22,2epsilon24,epsil$b_\text{nodes}$on33]M
    #        self.stress =  Stress()     #sigma11,sigma12,sigma13,sigma22,sigma23,sigma33]
    #        self.stressM = Stress()   #[sigma11,sigma12,sigma13,sigma22,sigma23,sigma33]M

    @property
    def theta_11(self):
        return self.theta_1[0]

    @property
    def area(self):
        return self.calc_area()

    @property
    def orientation(self):
        return self.calc_orientation()

    def __repr__(self):
        """
        tells Python how to represent an the Cell object (when using a print 
        statement) for a general purposes we use  def __repr__(self):
            
        Returns: String
        """
        STR = ""
        STR += str("Cell %s w. nodes:\t" % (self.id))
        for n in self.nodes:
            STR += str("%i, " % (n.id))

        return STR

    def __getstate__(self):
        """Return state values to be pickled."""
        return (self.id, self.nodes, self.theta_3, self.MatID, self.theta_1, self.structured, self.interior_nodes)

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.id, self.nodes, self.theta_3, self.MatID, self.theta_1, self.structured, self.interior_nodes = state
        # self.wire = self.build_wire()

    def calc_theta_1(self):
        """This method calculates the theta_1 vector. theta_1[0] represents the 
        ply coordinate system, which is formed by roating the global coordinate 
        system in the right-hand sense about x1 by the amount theta_1[0] 
        (theta_11). Afterwards the ply coordinate system us ritated avizt y3 in 
        the right-hand sens by the amount of Theta_3 to form the material 
        system. 
        For a detailed description see docs/man/VABS-Manual.pdf Figure 4.
        
        theta_11 is calculated as the angle between the x-axis and the Vector 
        from Node 1 to Node 2.
        
        Returns:
           None, but stores the theta_1 definition      
        """
        theta_1 = [0] * 9
        if self.structured:
            v0 = gp_Vec2d(gp_Pnt2d(0, 0), gp_Pnt2d(1, 0))
            v1 = gp_Vec2d(self.nodes[1].Pnt2d, self.nodes[2].Pnt2d)
            try:
                theta_11 = (v0.Angle(v1)) * 180 / np.pi
            except:
                print("WARNING: Vector with Null Magnitude at cell", self)
                theta_11 = 0
            # print 'v0 Magnitude:',v0.Magnitude(), 'v1 Magnitude:',v1.Magnitude(),
            if theta_11 < 0:
                theta_11 = 360 + theta_11
            theta_1[0] = theta_11
            theta_1[1] = 540
        else:
            theta_1[0] = 0
            theta_1[1] = 540
        self.theta_1 = theta_1
        return None

    def calc_area(self):
        """Calculates and returns the surface area of the cell"""
        corners = []
        for node in self.nodes:
            corners.append(node.coordinates)
        return PolygonArea(corners)

    def calc_orientation(self):
        """Calculates the orientation of the cell.
        
        Returns:
            - True: if counterclockwise 
            - False: else
        """
        corners = []
        for node in self.nodes:
            corners.append(node.coordinates)
        if PolygonArea(corners) > 0:
            return True
        else:
            return False

    def invert_nodes(self):
        self.nodes.reverse()

    def cell_node_distance(self, node):
        P_distances = []

        for i in range(0, len(self.nodes) - 1):
            spline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.nodes[i].Pnt2d, self.nodes[i + 1].Pnt2d])).Curve()
            projection = Geom2dAPI_ProjectPointOnCurve(node.Pnt2d, spline)

            for j in range(1, projection.NbPoints() + 1):
                P_distances.append(projection.Distance(j))

        return min(P_distances or [10e6])
    
    def split_quads(self):
        """method that splits quad cells into triangles and returns the list of
        cells [originalcell, newcell]"""

        if len(self.nodes) == 3:
            return [self]
        elif len(self.nodes) == 4:
            newcell = Cell([self.nodes[0], self.nodes[2], self.nodes[3]])
            newcell.theta_1 = self.theta_1
            newcell.theta_3 = self.theta_3
            newcell.MatID = self.MatID
            self.nodes = [self.nodes[0], self.nodes[1], self.nodes[2]]
            return [self, newcell]
        else:
            return []

    def closest_cell_edge(self, node):
        P_distances = []
        for i in range(0, len(self.nodes) - 1):
            spline = Geom2dAPI_PointsToBSpline(point2d_list_to_TColgp_Array1OfPnt2d([self.nodes[i].Pnt2d, self.nodes[i + 1].Pnt2d])).Curve()
            projection = Geom2dAPI_ProjectPointOnCurve(node.Pnt2d, spline)
            for j in range(1, projection.NbPoints() + 1):
                P_distances.append([projection.Distance(j), i])

        min_index, min_value = min(enumerate(P_distances), key=operator.itemgetter(1))
        return min_value[1]
