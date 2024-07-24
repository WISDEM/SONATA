# Third party modules
import numpy as np
# First party modules
from SONATA.cbm.mesh.mesh_core import gen_core_cells
from SONATA.cbm.mesh.mesh_utils import (grab_nodes_of_cells_on_BSplineLst,
    grab_nodes_on_BSplineLst, remove_dublicate_nodes,
    remove_duplicates_from_list_preserving_order,)
from SONATA.cbm.topo.BSplineLst_utils import (copy_BSplineLst,
                                              get_BSplineLst_Pnt2d,
                                              reverse_BSplineLst,
                                              trim_BSplineLst,)
from SONATA.cbm.topo.layer import Layer
from SONATA.cbm.topo.layer_utils import get_layer
from SONATA.cbm.topo.para_Geom2d_BsplineCurve import (BSplineLst_from_ParaLst,
                                                      ParaLst_from_BSplineLst,)
from SONATA.cbm.topo.projection import (
    chop_interval_from_layup,insert_interval_in_layup, inverse_relevant_cummulated_layup_boundaries, relevant_cummulated_layup_boundaries,
    sort_layup_projection,)
from SONATA.cbm.topo.wire_utils import build_wire_from_BSplineLst


class Segment(object):
    """ 
    The Segment object is constructed from multiple Layers obejcts. 
    Each Segment has one Segment boundary.
    """

    def __init__(self, ID=0, **kwargs):  # gets called whenever we create a new instance of a class
        """ Initialize with BSplineLst:             Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=True, Boundary = BSplineLst)
            Initialize with airfoil database:       Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, airfoil = 'naca23012')
            Initialize from file:                   Segment(ID, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i], OCC=False, filename = 'naca23012.dat')   
        #empty initialization with no Boundary: Segment(item, Layup = Configuration.Layup[1], CoreMaterial = Configuration.SEG_CoreMaterial[i])"""

        self.Segment0 = None
        self.WebLst = None
        self.BSplineLst = None
        self.ID = ID
        self.Layup = kwargs.get("Layup")
        self.CoreMaterial = kwargs.get("CoreMaterial")
        self.OCC = kwargs.get("OCC")
        self.Theta = kwargs.get("Theta")
        self.scale_factor = kwargs.get("scale_factor")

        self.LayerLst = []
        self.cells = []  # list of both layer and c_cells
        self.l_cells = []  # list of layer cells
        self.c_cells = []  # list of core cells
        self.boundary_ivLst = np.array([[0.0, 1.0, 0.0]])
        self.inv_cumivLst = []
        self.final_Boundary_ivLst = []
        self.Projection = relevant_cummulated_layup_boundaries(self.Layup)
        self.wire = None

        if self.OCC == True:
            self.BSplineLst = kwargs.get("Boundary")

        elif self.OCC == False:
            if kwargs.get("airfoil") != None:
                BSplineLst_tmp = self.BSplineLst_from_airfoil_database(kwargs.get("airfoil"), 30, self.scale_factor)


    def __repr__(self):
        return "{}: {}".format(self.__class__.__name__, self.ID)

    def __getstate__(self):
        """Return state values to be pickled."""
        self.Para_BSplineLst = ParaLst_from_BSplineLst(self.BSplineLst)
        return (
            self.ID,
            self.Layup,
            self.CoreMaterial,
            self.OCC,
            self.Theta,
            self.scale_factor,
            self.Projection,
            self.LayerLst,
            self.Para_BSplineLst,
            self.cells,
            self.l_cells,
            self.c_cells,
            self.boundary_ivLst,
            self.inv_cumivLst,
            self.final_Boundary_ivLst,
        )

    def __setstate__(self, state):
        """Restore state from the unpickled state values."""
        (
            self.ID,
            self.Layup,
            self.CoreMaterial,
            self.OCC,
            self.Theta,
            self.scale_factor,
            self.Projection,
            self.LayerLst,
            self.Para_BSplineLst,
            self.cells,
            self.l_cells,
            self.c_cells,
            self.boundary_ivLst,
            self.inv_cumivLst,
            self.final_Boundary_ivLst,
        ) = state
        self.BSplineLst = BSplineLst_from_ParaLst(self.Para_BSplineLst)
        self.build_wire()

    def ivLst_to_BSplineLst(self, ivLst):
        """The member function ivLst_to_BSplineLst generates the 
        BSplineLst from the InvervalLst definitions. It loops through all 
        intervals,trims them accordingly and assembles them into the 
        iv_BSplineLst, which is returned.
        
        Args:   self: the attributes of the Segment class
                iVLst: the iVLst is a entry of the Layup Projection. 
                It has the following form: array([[ 0.2  ,  0.3  ,  6.   ],
                                                  [ 0.3  ,  0.532,  8.   ]])
        returns: iv_BSplineLst: (list of BSplines)
        """
        if self.ID == 0:
            self.Segment0 = self

        iv_BSplineLst = []
        for iv in ivLst:
            if int(iv[2]) == 0:
                BSplineLst = self.BSplineLst
                start = 0.0
                end = 1.0

            elif int(iv[2]) < 0:
                WebID = -int(iv[2]) - 1
                if self.ID == WebID + 1:  # BACK
                    BSplineLst = reverse_BSplineLst(copy_BSplineLst(self.WebLst[WebID].BSplineLst))
                    start = self.WebLst[WebID].Pos2
                    end = self.WebLst[WebID].Pos1

                else:  # FRONT
                    BSplineLst = self.WebLst[WebID].BSplineLst
                    start = self.WebLst[WebID].Pos1
                    end = self.WebLst[WebID].Pos2

            elif self.ID * 1000 < int(iv[2]) < self.ID * 1000 + 1000:
                lid = iv[2] - (self.ID * 1000)
                layer = self.LayerLst[int(lid) - 1]
                BSplineLst = layer.BSplineLst
                start = layer.S1
                end = layer.S2

            else:
                layer = self.Segment0.LayerLst[int(iv[2]) - 1]
                BSplineLst = layer.BSplineLst
                start = layer.S1
                end = layer.S2

            S1 = iv[0]
            S2 = iv[1]
            iv_BSplineLst.extend(trim_BSplineLst(BSplineLst, S1, S2, start, end))
        return iv_BSplineLst

    def get_BsplineLst_plus(self, lid, SegmentLst, WebLst, layer_attr="a_BSplineLst"):
        BSplineLst = []
        start = None
        end = None

        if lid > 0:
            layer = get_layer(lid, SegmentLst)
            BSplineLst = getattr(layer, layer_attr)
            start = layer.S1
            end = layer.S2

        return (BSplineLst, start, end)

    def get_Pnt2d(self, S, SegmentLst, WebLst=None):
        """returns a Pnt2d for the coresponding layer number and the coordinate S"""
        a = self.final_Boundary_ivLst
        # select interval were S can be found

        for it, e in enumerate(a):
            if e[0] < S <= e[1]:
                lid = e[2]
                break

        (BSplineLst, start, end) = self.get_BsplineLst_plus(lid, SegmentLst, WebLst)
        return get_BSplineLst_Pnt2d(BSplineLst, S, start, end)
    def build_layers(self, WebLst=None, Segment0=None, display=None, l0=None, **kwargs):
        """The build_layers member function of the class Segment generates all Layer objects and it's associated wires
        and return the relevant_boundary_BSplineLst"""
        # ivLst is a list of the each layer start, end, and id
        # cum_ivLst is the cumulative list of ivLsts
        cum_ivLst = self.boundary_ivLst

        if self.Layup.size != 0:
            for i in range(1, len(self.Layup) + 1):
                print("STATUS:\t Building Segment %d, Layer: %d" % (self.ID, i))
                begin = float(self.Layup[i - 1][0])
                end = float(self.Layup[i - 1][1])
                # Deleting all ivs from the ivLst inside of the current layer's begin and end and storing
                # the deleted ivs in ivLst to use for the reference layer.
                ivLst = chop_interval_from_layup(cum_ivLst, begin, end)
                ivLst = sort_layup_projection([ivLst])[0]
                # Creating reference layer from the chopped ivLst
                relevant_boundary_BSplineLst = self.ivLst_to_BSplineLst(ivLst)

                # CREATE LAYER Object
                lid = int((self.ID * 1000) + i)
                tmp_Layer = Layer(
                    lid,
                    relevant_boundary_BSplineLst,
                    self.Layup[i - 1][0],
                    self.Layup[i - 1][1],
                    self.Layup[i - 1][2],
                    self.Layup[i - 1][3],
                    self.Layup[i - 1][4],
                    cutoff_style=2,
                    join_style=1,
                    name="test",
                )

                tmp_Layer.build_layer(l0)
                tmp_Layer.ivLst = ivLst

                if tmp_Layer.IsClosed:
                    # tmp_Layer.BSplineLst = set_BSplineLst_to_Origin(tmp_Layer.BSplineLst, self.Theta)
                    tmp_Layer.set_layer_origin()

                # tmp_Layer.cumB_ivLst = cummulated_layup_boundaries(self.Layup)[i-1]
                tmp_Layer.cumB_ivLst = cum_ivLst
                cum_ivLst = insert_interval_in_layup(cum_ivLst, begin, end, value=self.ID * 1000 + i)
                cum_ivLst = sort_layup_projection([cum_ivLst])[0]
                # tmp_Layer.cumA_ivLst = cummulated_layup_boundaries(self.Layup)[i]
                tmp_Layer.cumA_ivLst = cum_ivLst
                tmp_Layer.inverse_ivLst = inverse_relevant_cummulated_layup_boundaries(self.Layup)[i - 1]

                tmp_Layer.build_wire()

                self.LayerLst.append(tmp_Layer)

            return relevant_boundary_BSplineLst

    def mesh_layers(self, SegmentLst, global_minLen, WebLst=None, display=None, l0=None):
        """
        More Commenting!!!!
        """
        self.cells = []  # list of both layer and c_cells
        self.l_cells = []  # list of layer cells
        self.c_cells = []  # list of core cells

        if self.LayerLst:
            np.set_printoptions(suppress=True)
            # initialize inv_ivLst
            self.inv_cumivLst = np.array([[0, 1, self.LayerLst[-1].ID + 1]])

            if self.ID == 0:  # concatenate ivLsts of the previous segments!
                ivCollector = self.inv_cumivLst
                for seg in SegmentLst[1:]:
                    if seg.ID == 1:
                        tmp_ivLst = chop_interval_from_layup(seg.inv_cumivLst, WebLst[0].Pos1, WebLst[0].Pos2)
                        for iv in tmp_ivLst:
                            ivCollector = insert_interval_in_layup(ivCollector, iv[0], iv[1], value=iv[2])

                    elif seg.ID == len(WebLst) + 1:
                        tmp_ivLst = chop_interval_from_layup(seg.inv_cumivLst, WebLst[-1].Pos2, WebLst[-1].Pos1)
                        for iv in tmp_ivLst:
                            ivCollector = insert_interval_in_layup(ivCollector, iv[0], iv[1], value=iv[2])
                    else:
                        # print seg.inv_cumivLst, WebLst[seg.ID-2].Pos2,WebLst[seg.ID-1].Pos2
                        tmp_ivLst = chop_interval_from_layup(seg.inv_cumivLst, WebLst[seg.ID - 1].Pos1, WebLst[seg.ID - 2].Pos1)
                        for iv in tmp_ivLst:
                            ivCollector = insert_interval_in_layup(ivCollector, iv[0], iv[1], value=iv[2])
                        tmp_ivLst = chop_interval_from_layup(seg.inv_cumivLst, WebLst[seg.ID - 2].Pos2, WebLst[seg.ID - 1].Pos2)
                        for iv in tmp_ivLst:
                            ivCollector = insert_interval_in_layup(ivCollector, iv[0], iv[1], value=iv[2])
                    # print sort_layup_projection([ivCollector])[0]

                ivCollector = sort_layup_projection([ivCollector])[0]
                self.inv_cumivLst = ivCollector

            for i, layer in enumerate(reversed(self.LayerLst)):
                print("STATUS:\t Meshing Segment %s, Layer %s" % (self.ID, len(self.LayerLst) - i))

                layer.inverse_ivLst = chop_interval_from_layup(self.inv_cumivLst, layer.S1, layer.S2)
                layer.inverse_ivLst = sort_layup_projection([layer.inverse_ivLst])[0]
                layer.mesh_layer(SegmentLst, global_minLen, display=display, l0=1.5 * l0)
                self.inv_cumivLst = insert_interval_in_layup(self.inv_cumivLst, layer.S1, layer.S2, value=layer.ID)
                self.l_cells.extend(layer.cells)

            self.cells.extend(self.l_cells)
            return self.l_cells
        else:
            return []

    def mesh_core(self, SegmentLst, WebLst, core_cell_area, display=None):
        if self.ID == 0 and len(SegmentLst) > 1:
            pass

        elif self.CoreMaterial == 0:
            pass

        else:
            print("STATUS:\t Meshing Segment %s, Core" % (self.ID))
            # print self.final_Boundary_ivLst
            core_a_nodes = []
            for iv in self.final_Boundary_ivLst:
                (BSplineLst, start, end) = self.get_BsplineLst_plus(int(iv[2]), SegmentLst, WebLst, layer_attr="a_BSplineLst")
                iv_BSplineLst = trim_BSplineLst(BSplineLst, iv[0], iv[1], start, end)

                # print iv, "use nodes a_nodes of layer", int(iv[2])
                tmp_layer = get_layer(iv[2], SegmentLst)
                disco_nodes = grab_nodes_on_BSplineLst(tmp_layer.a_nodes, iv_BSplineLst)
                core_a_nodes.extend(disco_nodes)

            core_a_nodes = remove_dublicate_nodes(core_a_nodes)
            core_a_nodes = remove_duplicates_from_list_preserving_order(core_a_nodes)
            [c_cells, c_nodes] = gen_core_cells(core_a_nodes, core_cell_area)

            for c in c_cells:
                c.structured = False
                c.theta_3 = 0
                c.MatID = int(self.CoreMaterial)
                c.calc_theta_1()
                if c.area < 1e-10:
                    print(c.nodes)
            #                    display.DisplayShape(c.nodes[0].Pnt2d)
            #                    display.DisplayShape(c.nodes[1].Pnt2d, color='RED')
            #                    display.DisplayShape(c.nodes[2].Pnt2d, color='GREEN')
            #                    display.DisplayShape(c.wire,color='WHITE')

            self.c_cells.extend(c_cells)
            self.cells.extend(self.c_cells)
        return self.c_cells

    def determine_final_boundary(self, WebLst=None, Segment0=None):
        """The member function determin_final_boundary2 generates the 
        BSplineLst that encloses all Layers of the Segement. This final 
        boundary is needed for the generation of the subordinate segments, 
        where the final boundary BSplineLst is intersected with the Webs.
        
        Args:   self: only the attributes of the Segment class itself.
        
        returns: None, but assignes the final_Boundary_BSplineLst class 
            attribute
        """
        cum_ivLst = self.boundary_ivLst
        if self.Layup.size != 0:
            for i in range(1, len(self.Layup) + 1):
                cum_ivLst = insert_interval_in_layup(cum_ivLst, float(self.Layup[i - 1][0]), float(self.Layup[i - 1][1]), value=self.ID * 1000 + i)

        self.final_Boundary_ivLst = sort_layup_projection([cum_ivLst])[0]
        self.final_Boundary_BSplineLst = self.ivLst_to_BSplineLst(self.final_Boundary_ivLst)
        return None

    def build_wire(self):  # Builds TopoDS_Wire from connecting BSplineSegments and returns it
        self.wire = build_wire_from_BSplineLst(self.BSplineLst)

    def build_segment_boundary_from_WebLst(self, WebLst, Segment0):
        print("STATUS:\t Building Segment Boundaries %s" % (self.ID))
        i = self.ID - 1

        if self.ID == 0:
            self.boundary_ivLst = None

        if self.ID == 1:
            # CREATE SEGMENT BOUNDARY 1
            ivLst = chop_interval_from_layup(Segment0.final_Boundary_ivLst, WebLst[i].Pos1, WebLst[i - 1].Pos2)
            ivLst = insert_interval_in_layup(ivLst, WebLst[i].Pos2, WebLst[i].Pos1, value=-WebLst[i].ID)
            self.boundary_ivLst = sort_layup_projection([ivLst])[0]

        elif self.ID == len(WebLst) + 1:
            # CREATE LAST BOUNDARY
            ivLst = chop_interval_from_layup(Segment0.final_Boundary_ivLst, WebLst[i - 1].Pos2, WebLst[i - 1].Pos1)
            ivLst = insert_interval_in_layup(ivLst, WebLst[i - 1].Pos1, WebLst[i - 1].Pos2, value=-WebLst[i - 1].ID)
            self.boundary_ivLst = sort_layup_projection([ivLst])[0]

        else:
            ivLst1 = chop_interval_from_layup(Segment0.final_Boundary_ivLst, WebLst[i].Pos1, WebLst[i - 1].Pos1)
            ivLst1 = insert_interval_in_layup(ivLst1, WebLst[i - 1].Pos1, WebLst[i - 1].Pos2, value=-WebLst[i - 1].ID)
            ivLst2 = chop_interval_from_layup(Segment0.final_Boundary_ivLst, WebLst[i - 1].Pos2, WebLst[i].Pos2)
            ivLst2 = insert_interval_in_layup(ivLst2, WebLst[i].Pos2, WebLst[i].Pos1, value=-WebLst[i].ID)
            ivLst = np.vstack((ivLst1, ivLst2))
            self.boundary_ivLst = sort_layup_projection([ivLst])[0]
        self.BSplineLst = self.ivLst_to_BSplineLst(self.boundary_ivLst)
        self.build_wire()

        return None
