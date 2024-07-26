# Third party modules
import numpy as np
from OCC.Core.gp import gp_Pnt2d
from scipy.spatial import distance
import matplotlib.pyplot as plt

# First party modules
from SONATA.cbm.mesh.mesh_byprojection import \
    mesh_by_projecting_nodes_on_BSplineLst
from SONATA.cbm.mesh.mesh_improvements import (
    modify_sharp_corners, second_stage_improvements,)
from SONATA.cbm.mesh.mesh_utils import (
    equidistant_nodes_on_BSplineLst,
    grab_nodes_on_BSplineLst, merge_nodes_if_too_close,
    remove_duplicates_from_list_preserving_order,)
from SONATA.cbm.topo.BSplineLst_utils import (BSplineLst_from_dct,
                                              ProjectPointOnBSplineLst,
                                              copy_BSpline,
                                              discretize_BSplineLst,
                                              findPnt_on_BSplineLst,
                                              trim_BSplineLst,)
from SONATA.cbm.topo.cutoff import cutoff_layer
from SONATA.cbm.topo.layer_utils import get_layer
from SONATA.cbm.topo.offset import shp_parallel_offset
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import (BSplineLst_from_ParaLst,
                                                      ParaLst_from_BSplineLst,)
from SONATA.cbm.topo.utils import isclose
from SONATA.cbm.topo.wire_utils import build_wire_from_BSplineLst


class Layer(object):
    """ 
    The layer object is constructed from multiple BSplineCurveSegments. It is the basis for all future operations. 
    The object can be constructed from either a discrete formulation of point tables or from an existing TopoDS_Wire.
    """

    def __init__(self, ID, Boundary_BSplineLst, globalStart, globalEnd, thickness, Orientation=0, MatID=1, **kwargs):
        self.ID = ID  # First single Digit: SegmentNb, Last 3 Digits: LayerNb; e.g.: 1029, Segment1, Layer29
        self.Boundary_BSplineLst = Boundary_BSplineLst  # List of Geom2d_BSplineCurve, Geom_BSplineCurve
        self.S1 = globalStart  # Starting Point in S coordinates
        self.S2 = globalEnd  # End Point in S coordinates
        self.thickness = thickness  # in mm
        self.Orientation = Orientation  # in deg
        self.MatID = MatID  # type: int
        self.cells = []  # type: List
        self.ivLst = []  # type: List
        self.cumA_ivLst = []  # type: List
        self.cumB_ivLst = []  # type: List
        self.inverse_ivLst = []  # type: List
        self.a_nodes = []  # type: List
        self.b_nodes = []  # type: List

        # KWARGS:
        if kwargs.get("name") == None:
            self.name = "DEFAULT"
        else:
            self.name = kwargs.get("name")

        if (kwargs.get("cutoff_style") == None) or (type(kwargs.get("cutoff_style")) is not int):  # cutoff_style (step, linear, smooth_bezier)
            self.cutoff_style = 2
        else:
            self.cutoff_style = kwargs.get("cutoff_style")

        if (kwargs.get("join_style") == None) or (type(kwargs.get("join_style")) is not int):  # offset algorithm join_style = 1#( 1:round,2:mitre,3:bevels)
            self.join_style = 1
        else:
            self.join_style = kwargs.get("join_style")

    @property
    def StartPoint(self):  # gp_Pnt2d
        return self.BSplineLst[0].StartPoint()

    @property
    def EndPoint(self):  # gp_Pnt2d
        return self.BSplineLst[-1].EndPoint()

    @property
    def a_BSplineLst(self):  # gp_Pnt2d
        return self.BSplineLst

    @property
    def b_BSplineLst(self):  # gp_Pnt2d
        return self.Boundary_BSplineLst

    @property
    def IsClosed(self):
        return self.a_BSplineLst[0].StartPoint().IsEqual(self.a_BSplineLst[-1].EndPoint(), 1e-5)

    def __str__(self):
        # we can tell Python how to prepresent an object of our class (when using a print statement) for general purposes use  __repr__(self):
        return str(
            "LayerID: \tStart[-]: \tEnd[-]: \tthickness[-]: \tOrientation[deg]: \tMatID \tName:\t\n"
            "%s, \t%s, \t%s, \t%s, \t\t%s, \t\t%s, \t%s, " % (self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.name)
        )

    def __getstate__(self):
        """Return state values to be pickled."""
        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        self.Para_Boundary_BSplineLst = ParaLst_from_BSplineLst(self.Boundary_BSplineLst)
        return (self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.Para_Boundary_BSplineLst, self.Para_BSplineLst)

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        self.ID, self.S1, self.S2, self.thickness, self.Orientation, self.MatID, self.Para_Boundary_BSplineLst, self.Para_BSplineLst = state
        self.Boundary_BSplineLst = BSplineLst_from_ParaLst(self.Para_Boundary_BSplineLst)
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)
        self.build_wire()

    def build_wire(self):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)

    def display_bsplinelst(self, bsplinelst, color = 'red'):
        for bspline in bsplinelst:
            u_min, u_max = bspline.FirstParameter(), bspline.LastParameter()
            # Extract points for plotting
            num_points = 100  # Number of points to plot
            u_values = [u_min + (u_max - u_min) * i / (num_points - 1) for i in range(num_points)]
            x_values = [bspline.Value(u).X() for u in u_values]
            y_values = [bspline.Value(u).Y() for u in u_values]
        
            plt.plot(x_values, y_values, color = color)

    # Function to check if two line segments (p1, q1) and (p2, q2) intersect
    def do_intersect(self, p1, q1, p2, q2):
        def orientation(p, q, r):
            val = (float(q[1] - p[1]) * (r[0] - q[0])) - (float(q[0] - p[0]) * (r[1] - q[1]))
            if val > 0:
                return 1  # Clockwise
            elif val < 0:
                return 2  # Counterclockwise
            else:
                return 0  # Collinear

        def on_segment(p, q, r):
            if min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1]):
                return True
            return False

        o1 = orientation(p1, q1, p2)
        o2 = orientation(p1, q1, q2)
        o3 = orientation(p2, q2, p1)
        o4 = orientation(p2, q2, q1)

        # General case
        if o1 != o2 and o3 != o4:
            return True

        # Special cases
        if o1 == 0 and on_segment(p1, q1, p2):
            return True
        if o2 == 0 and on_segment(p1, q1, q2):
            return True
        if o3 == 0 and on_segment(p2, q2, p1):
            return True
        if o4 == 0 and on_segment(p2, q2, q1):
            return True

        return False

    # Function to check if a shape intersects itself
    def shape_intersects_itself(self, coords):
        num_points = len(coords)
        for i in range(num_points - 1):  # Loop to (num_points - 1) to avoid out-of-bounds access
            for j in range(i + 2, num_points - (1 if i == 0 else 0)):  # Adjusting the range to avoid out-of-bounds
                if self.do_intersect(coords[i], coords[i + 1], coords[j], coords[(j + 1) % num_points]):
                    return True
        return False 

    def build_layer(self, l0=1):
        npArray = discretize_BSplineLst(self.Boundary_BSplineLst, 1.2e-6 * l0)
        self.offlinepts = shp_parallel_offset(npArray, self.thickness, self.join_style)
        if self.shape_intersects_itself(self.offlinepts[1:-2]):
            print("WARNING: THERE IS AN INTERSECTION IN THE STRUCTURE")
        OffsetBSplineLst = BSplineLst_from_dct(self.offlinepts, angular_deflection=15, tol_interp=1e-8 * l0, cutoff_style = 2)
        OffsetBSplineLst = cutoff_layer(self.Boundary_BSplineLst, OffsetBSplineLst, self.S1, self.S2, self.cutoff_style)
        self.BSplineLst = OffsetBSplineLst
        self.display_bsplinelst(self.BSplineLst)

    def determine_a_nodes(self, SegmentLst, global_minLen, display=None):
        """ """
        unmeshed_ids = []
        for seg in SegmentLst:
            if seg.LayerLst:
                unmeshed_ids.append(int(seg.LayerLst[-1].ID + 1))

        new_a_nodes = []
        for iv_counter, iv in enumerate(self.inverse_ivLst):
            if int(iv[2]) in unmeshed_ids:  # if
                # print iv, "equidistand nodes on BsplineLst of LayerLst entry"
                eq_nodes = []
                BSplineLst = self.a_BSplineLst
                iv_BSplineLst = trim_BSplineLst(BSplineLst, iv[0], iv[1], self.S1, self.S2)
                if iv_counter == 0 and len(self.inverse_ivLst) > 1:  # first but not last
                    IncStart = True
                    IncEnd = False

                elif iv_counter == 0 and len(self.inverse_ivLst) == 1:  # first and last
                    IncStart = True
                    IncEnd = True

                elif iv_counter == len(self.inverse_ivLst) - 1 and len(self.inverse_ivLst) > 1:  # last but not first
                    # print iv, 'Last but not first'
                    if iv_counter == len(self.inverse_ivLst) - 1 and iv[1] == 1 and self.inverse_ivLst[0][0] == 0:
                        IncStart = False
                        IncEnd = False
                    else:
                        IncStart = False
                        IncEnd = True

                else:
                    IncStart = False
                    IncEnd = False

                eq_nodes = equidistant_nodes_on_BSplineLst(iv_BSplineLst, True, IncStart, IncEnd, minLen=global_minLen, LayerID=self.ID)
                new_a_nodes.extend(eq_nodes)

            else:
                # only use once for each layer!
                # print iv, "use nodes b_nodes of layer", int(iv[2])
                tmp_layer = get_layer(int(iv[2]), SegmentLst)
                # iv_BSplineLst = trim_BSplineLst(tmp_layer.b_BSplineLst,iv[0],iv[1],tmp_layer.S1,tmp_layer.S2)
                iv_BSplineLst = trim_BSplineLst(self.a_BSplineLst, iv[0], iv[1], self.S1, self.S2)
                # iv_BSplineLst = self.a_BSplineLst
                isClosed = iv_BSplineLst[0].StartPoint().IsEqual(iv_BSplineLst[-1].EndPoint(), 1e-5)

                if isClosed:
                    tmp_nodes = tmp_layer.b_nodes
                else:
                    tmp_nodes = [tmp_layer.a_nodes[0]] + tmp_layer.b_nodes + [tmp_layer.a_nodes[-1]]

                disco_nodes = grab_nodes_on_BSplineLst(tmp_nodes, iv_BSplineLst)
                new_a_nodes.extend(disco_nodes)

        self.a_nodes = remove_duplicates_from_list_preserving_order(new_a_nodes)
        self.a_nodes = merge_nodes_if_too_close(self.a_nodes, self.a_BSplineLst, global_minLen, 0.01)

    def mesh_layer(self, SegmentLst, global_minLen, proj_tol_1=9e-2, proj_tol_2=4e-1, crit_angle_1=110, alpha_crit_2=60, growing_factor=1.8, shrinking_factor=0.01, display=None, l0=None):
        """
        The mesh layer function discretizes the layer, which is composed of a 
        a_BsplineLst and a b_BsplineLst. Between the a_BsplineLst and the 
        b_BsplineLst the cells are created. First nodes on the a_BSplineLst are
        determined with the determine_a_nodes procedure. Subsequently the 
        a_nodes are projected onto the b_BsplineLst. The following functions 
        (modify_cornerstyle_one, modify_sharp_corners and
        second_stage_improvements) try to improve the quality of the mesh.
        everything is stored in the layer.cells and is returned
        
        Parameters
        --------
        SegmentLst:               The overall list of Segments within the 
                                segemet. This list is needed to determine
                                the a_nodes
        global_minLen:          
        proj_tol_1 = 5e-2:      tolerance value to determine a distance, 
                                in which the resulting projection point 
                                has to be. 
                                (mesh_by_projecting_nodes_on_BSplineLst)
        proj_tol_2 = 2e-1:      tolerance value to determine a distance, 
                                in which the resulting projection point 
                                has to be. (modify_sharp_corners)
        crit_angle_1 = 115:     is the critical angle to determine a corner 
                                if 2 projection points are found.    
        alpha_crit_2 = 60:      is the critical angle to refine  a corner 
        growing_factor = 1.8:   critical growing factor of cell before 
                                splitting 
        shrinking_factor = 0.10:  critical shrinking factor for cells 
                                before merging nodes
        
        Returns
        -------
        self.cells: (list of cells) 
        """
        self.determine_a_nodes(SegmentLst, global_minLen, display)
        self.a_nodes, self.b_nodes, self.cells = mesh_by_projecting_nodes_on_BSplineLst(
            self.a_BSplineLst, self.a_nodes, self.b_BSplineLst, self.thickness, proj_tol_1, crit_angle_1, LayerID=self.ID, refL=l0, display=display
        )
        # enhanced_cells = modify_cornerstyle_one(cells,self.b_BSplineLst)
        self.cells, nb_nodes = modify_sharp_corners(self.cells, self.b_BSplineLst, global_minLen, self.thickness, self.ID, proj_tol_2, alpha_crit_2, display=display)
        self.b_nodes.extend(nb_nodes)
        try:
            self.cells, nb_nodes = second_stage_improvements(self.cells, self.b_BSplineLst, global_minLen, self.ID, growing_factor, shrinking_factor, display=display)
            self.b_nodes.extend(nb_nodes)
        except:
            pass

        # self.b_nodes = sorted(self.b_nodes, key=lambda Node: (Node.parameters[1],Node.parameters[2]))

        for c in self.cells:
            c.calc_theta_1()
            c.theta_3 = self.Orientation
            c.MatID = int(self.MatID)
            c.structured = True

        return self.cells

    def set_layer_origin(self):
        """
        this procedure reorders the self.BSplineLst to an origin if the layer 
        is closed. The Origin is detected by searching for an orthogonal 
        projection of the StartPoint of the self.Boundary_BSplineLst. If no 
        projection is found it takes the closest neighbor of the discrete 
        offlinepts (Offset Line Points).
                
        """
        # Determine Origin as point
        if self.IsClosed:
            org = self.Boundary_BSplineLst[0].StartPoint()
            proj = ProjectPointOnBSplineLst(self.BSplineLst, org, 3 * self.thickness)

            if len(proj) > 0:
                OriPnt = proj[0]
            elif len(proj) == 0:
                distarr = distance.cdist(self.offlinepts, np.asarray([org.Coord()]))
                OriPnt = gp_Pnt2d(self.offlinepts[distarr.argmin()][0], self.offlinepts[distarr.argmin()][1])

        else:
            print("Only apply member layer.set_layer_origin method when layer is Closed!")

        # print(OriPnt.Coord())
        # Reorder Sequence of BSplines of BSplinesLst
        OriPara = findPnt_on_BSplineLst(OriPnt, self.BSplineLst)
        OBSplineLst = []
        CorrectOrigin = False

        for i, item in enumerate(self.BSplineLst):
            if i == OriPara[0]:
                First = item.FirstParameter()
                Last = item.LastParameter()

                if isclose(OriPara[1], First) == True:
                    OBSplineLst.append(item)
                    CorrectOrigin = True

                elif isclose(OriPara[1], Last) == True:
                    CorrectOrigin = False
                    BSplineCurve2 = item

                else:
                    CorrectOrigin = False
                    BSplineCurve1 = copy_BSpline(item)
                    BSplineCurve1.Segment(OriPara[1], Last)
                    BSplineCurve2 = copy_BSpline(item)
                    BSplineCurve2.Segment(First, OriPara[1])
                    OBSplineLst.append(BSplineCurve1)

            elif i > OriPara[0]:
                OBSplineLst.append(item)
            else:
                None

        for i, item in enumerate(self.BSplineLst):
            if i < OriPara[0]:
                OBSplineLst.append(item)
            else:
                None

        if CorrectOrigin == False:
            OBSplineLst.append(BSplineCurve2)
        else:
            None

        #        if BSplineLst_Orientation(OBSplineLst,11) == True:
        #            OBSplineLst = reverse_BSplineLst(OBSplineLst)

        self.BSplineLst = OBSplineLst
        return None


# execute the following code if this file is executed as __main__
if __name__ == "__main__":
    L1 = Layer()
    pass
