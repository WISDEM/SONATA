# -*- coding: utf-8 -*-
"""Defines the Crosssectional Beam Model (CBM) class
Created on Wed Jan 03 13:56:37 2018
@author: TPflumm

https://numpydoc.readthedocs.io/en/latest/format.html
 """

# Core Library modules
import os
import copy
import math
from collections import OrderedDict

# Basic PYTHON Modules:
import pickle as pkl
from datetime import datetime

# Third party modules
import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.gp import gp_Ax2

# First party modules
from SONATA.cbm.cbm_utl import trsf_sixbysix
from SONATA.cbm.classBeamSectionalProps import BeamSectionalProps
from SONATA.cbm.classCBMConfig import CBMConfig
from SONATA.cbm.display.display_mesh import plot_cells
from SONATA.cbm.mesh.cell import Cell
from SONATA.cbm.mesh.consolidate_mesh import consolidate_mesh_on_web
from SONATA.cbm.mesh.mesh_utils import (grab_nodes_of_cells_on_BSplineLst,
                                        sort_and_reassignID,)
from SONATA.cbm.mesh.node import Node
from SONATA.cbm.topo.BSplineLst_utils import (get_BSplineLst_length,)
from SONATA.cbm.topo.segment import Segment
from SONATA.cbm.topo.utils import getID
from SONATA.cbm.topo.web import Web
from SONATA.cbm.topo.weight import Weight
from SONATA.vabs.classStrain import Strain
from SONATA.vabs.classStress import Stress
from SONATA.cbm.topo.projection import (
    chop_interval_from_layup, sort_layup_projection,)
# from SONATA.vabs.classVABSConfig import VABSConfig
from SONATA.classMaterial import IsotropicMaterial, OrthotropicMaterial

from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.Core.gp import gp_Pnt2d

try:
    import dolfin as do
    from SONATA.anbax.anbax_utl import build_dolfin_mesh, anbax_recovery, ComputeShearCenter, ComputeTensionCenter, ComputeMassCenter
    import sys
    from anba4.anbax import anbax



except:
    print("dolfin and anbax could not be imported!")
    pass


class CBM(object):
    """ 
    This Class includes the SONATA Dicipline Module for Structural 
    Composite Beam Modelling (CBM). 
    Design Variables are passed in form of the configuration Object or a 
    configuration file.          

    Attributes
    ----------
    config : configuration
        Pointer to the Configuration object.
        
    materials: list
        List of Materials object instances.
    
    SegmentLst: list
        list of Segment object instances
    
    refL : float, default: 1.0
        reference length to account for different dimensions and sizes in the 
        cross-section. This length is approximately the circumference of the 
        outer curve.  

    Methods
    -------     
    cbm_stpexport_topo(export_filename=None)
        exports all Layer wires and the Segment0 Boundary Wire as .step
       
    __cbm_generate_SegmentLst(**kwargs)
        psydo private method of the cbm class to generate the list of 
        Segments in the instance. 
        
    cbm_gen_topo(**kwargs)
        generates the topology.
         
    cbm_gen_mesh(**kwargs)
        generates the dicretization of topology 
        
    cbm_review_mesh()
        prints a summary of the mesh properties to the screen 
        
    cbm_run_vabs()
        runs the solver VABS (Variational Asymptotic Beam Sectional Analysis)
    
    cbm_run_anbax()
        runs the solver anbax from macro morandini
        
    cbm_post_2dmesh(attribute='MatID', title='NOTITLE', **kw)
        displays the mesh with specified attributes with matplotlib
      
    cbm_post_3dtopo()
        displays the topology with the pythonocc 3D viewer
    
    cbm_post_3dmesh()
        displays the 2D mesh with the pythonocc 3D viewer
        
    cbm_set_DymoreMK(x_offset=0)
        Converts the Units of CBM to DYMORE/PYMORE/MARC units and returns the 
        array of the beamproperties with Massterms(6), Stiffness(21), 
        damping(1) and curvilinear coordinate(1)
    
    
    Notes
    ----------
    1 
        For computational efficiency it is sometimes not suitable to recalculate 
        the topology or the crosssection every iteration, maybe design flags 
        to account for that.
    2 
        make it possible to construct an instance by passing the topology 
        and/or mesh


    Examples
    --------
    >>> config = CBMConfig(fname)
    >>> job = CBM(config)
    >>> job.cbm_gen_topo()
    >>> job.cbm_gen_mesh()
    >>> job.cbm_review_mesh()
    >>> job.cbm_run_vabs()
    >>> job.cbm_post_2dmesh(title='Hello World!')

    """

    # __slots__ = ('config' , 'materials' , 'SegmentLst' , 'WebLst' , 'BW' , 'mesh', 'BeamProperties', 'display', )
    def __init__(self, Configuration, materials=None, **kwargs):
        """
        Initialize attributes.

        Parameters
        ----------
        Configuration : <Configuration>
            Pointer to the <Configuration> object.
        materials : dict(id, Materials)
        """

        self.config = Configuration
        if isinstance(materials, dict):
            self.materials = materials
        else:
            print("Materials not in dictionary format. Check yaml input file.")

        self.name = "cbm_noname"
        if kwargs.get("name"):
            self.name = kwargs.get("name")

        self.Ax2 = gp_Ax2()
        self.SegmentLst = []
        self.WebLst = []
        self.BW = None

        self.mesh = []
        self.BeamProperties = None
        self.display = None

        self.refL = 1.0

        self.startTime = datetime.now()
        self.exportLst = []  # list that contains all objects to be exported as step
        self.surface3d = None  # TODO: Remove definition and set it up in classBlade
        self.Blade = None  # TODO: Remove definition and set it up in classBlade

        # wire = kwargs.get('wire') #in the blade reference frame!
        self.Ax2 = kwargs.get("Ax2")
        self.BoundaryBSplineLst = kwargs.get("BSplineLst")
        self.Theta = 0
        self.refL = get_BSplineLst_length(self.BoundaryBSplineLst)

        self.cutoff_style = kwargs.get("cutoff_style")

    def _cbm_generate_SegmentLst(self, **kwargs):
        """
        psydo private method of the cbm class to generate the list of 
        Segments in the instance. 
        
        """
        self.SegmentLst = []  # List of Segment Objects

        # TODO cleanup this mess!
        for k, seg in self.config.segments.items():
            if k == 0:

                self.SegmentLst.append(Segment(k, **seg, Theta=self.Theta, OCC=True, Boundary=self.BoundaryBSplineLst))

            else:
                self.SegmentLst.append(Segment(k, **seg, Theta=self.Theta))

        sorted(self.SegmentLst, key=getID)
        self.refL = get_BSplineLst_length(self.SegmentLst[0].BSplineLst)
        return None
    
    def find_bspline_ends(self, bspline):
        start_point = gp_Pnt2d()
        end_point = gp_Pnt2d()
        bspline.D0(bspline.FirstParameter(), start_point)
        bspline.D0(bspline.LastParameter(), end_point)
        return start_point, end_point

    def check_bspline_intersections(self, Boundary_BSplineLst):
        # Checking for bspline intersections in a bspline list not including start and end points
        start_tol = 1e-5
        intersection_points = []
        intersected = False
        for i in range(len(Boundary_BSplineLst)):
            for j in range(len(Boundary_BSplineLst)):
                if i != j:
                    start1, end1 = self.find_bspline_ends(Boundary_BSplineLst[i])
                    start2, end2 = self.find_bspline_ends(Boundary_BSplineLst[j])
                    intersector = Geom2dAPI_InterCurveCurve(Boundary_BSplineLst[i],Boundary_BSplineLst[j])
                    if intersector.NbPoints() > 0:
                        for k in range(1, intersector.NbPoints()+1):
                            if not intersector.Point(k).IsEqual(start1, start_tol) and not intersector.Point(k).IsEqual(end1, start_tol) and not intersector.Point(k).IsEqual(start2, start_tol) and not intersector.Point(k).IsEqual(end2, start_tol):
                                intersected = True
                                intersection_points.append(intersector.Point(k))
        return [intersected, intersection_points]
    
    def display_bsplinelst(self, bsplinelst, color = 'blue'):
        for bspline in bsplinelst:
            u_min, u_max = bspline.FirstParameter(), bspline.LastParameter()
            # Extract points for plotting
            num_points = 100  # Number of points to plot
            u_values = [u_min + (u_max - u_min) * i / (num_points - 1) for i in range(num_points)]
            x_values = [bspline.Value(u).X() for u in u_values]
            y_values = [bspline.Value(u).Y() for u in u_values]
            plt.plot(x_values, y_values, color = color)

    def check_for_bspline_intersections(self, segment):
        # Getting all boundary Bsplines
        ivLst = chop_interval_from_layup(segment.boundary_ivLst, 0, 1)
        ivLst = sort_layup_projection([ivLst])[0]
        # Creating reference layer from the chopped ivLst
        Boundary_BSplineLst = segment.ivLst_to_BSplineLst(ivLst)
        [intersected, intersection_pnt] = self.check_bspline_intersections(Boundary_BSplineLst)
        if intersected:
            print("WARNING: There is an intersection in the structure.")
            plt.figure()
            self.display_bsplinelst(self.SegmentLst[0].BSplineLst, 'black')
            self.display_bsplinelst(Boundary_BSplineLst, 'blue')
            for points in intersection_pnt:
                plt.plot(points.X(), points.Y(), 'x', color = 'red', linewidth = 4, markersize = 10)

    def cbm_gen_topo(self, **kwargs):
        """
        CBM Method that generates the topology. It starts by generating the 
        list of Segments. It continous to gen all layers for Segment 0. 
        Subsequently the webs are defined and afterwards the layers of the 
        remaining Segments are build. The Balance Weight is defined at the end
        of this method.        
        """
        # Generate SegmentLst from config:
        self.SegmentLst = []
        self._cbm_generate_SegmentLst(**kwargs)
        # Build Segment 0:
        self.SegmentLst[0].build_wire()
        self.SegmentLst[0].build_layers(l0=self.refL, cutoff_style = self.cutoff_style, **kwargs)
        self.SegmentLst[0].determine_final_boundary()
        self.check_for_bspline_intersections(self.SegmentLst[0])

        # Build Webs:
        self.WebLst = []
        if len(self.config.webs) > 0:
            for k, w in self.config.webs.items(): 
                print('STATUS:\t Building Web %s' %(k+1))
                self.WebLst.append(Web(k, w['Pos1'], w['Pos2'], w['curvature'], self.SegmentLst))
            sorted(self.SegmentLst, key=getID)  
            
        #Build remaining Segments:
        if len(self.config.webs) > 0:
            for i, seg in enumerate(self.SegmentLst[1:], start=1):
                seg.Segment0 = self.SegmentLst[0]
                seg.WebLst = self.WebLst
                seg.build_segment_boundary_from_WebLst(self.WebLst, self.SegmentLst[0])
                seg.build_layers(self.WebLst, self.SegmentLst[0], l0=self.refL, cutoff_style = self.cutoff_style)
                seg.determine_final_boundary(self.WebLst, self.SegmentLst[0])
                self.check_for_bspline_intersections(seg)
                seg.build_wire()

        self.BW = None
        # Balance Weight:
        if self.config.setup["BalanceWeight"] == True:
            # print('STATUS:\t Building Balance Weight')
            # self.BW = Weight(0, self.config.bw['XPos'], self.config.bw['YPos'], self.config.bw['Diameter'], self.config.bw['Material'])
            p = self.SegmentLst[0].det_weight_Pnt2d(self.config.bw["s"], self.config.bw["t"])
            self.BW = Weight(0, p, self.config.bw["Diameter"], self.config.bw["Material"])
            # print(p.Coord())
        return None

    def cbm_gen_mesh(self, **kwargs):
        """
        CBM Method that generates the dicretization of topology and stores the 
        cells and nodes in both the <Layer> instances, the <Segment> instances 
        and the attribute self.mesh that is a list of <Cell> instances


        Parameters:
        ----------
        split_quads : bool, optional
            This option can be passed as keyword argument and splits the quads 
            (4 node cells) in mesh into two cells of 3 nodes      
        
        
        Notes:
        ----------  
        More option keyword arguments shall be possible in the future      
        
        Examples:
        ----------  
        >>> job.cbm_gen_mesh(splitquads=True)
        
        """

        split_quads = False
        if "split_quads" in kwargs:
            if type(kwargs["split_quads"]) == bool:
                split_quads = kwargs["split_quads"]
            else:
                print("split_quads must provide a boolean value")

        self.mesh = []
        Node.class_counter = 1
        Cell.class_counter = 1
        # meshing parameters:
        Resolution = self.config.setup["mesh_resolution"]  # Nb of Points on Segment0
        global_minLen = round(self.refL / Resolution, 5)

        core_cell_area = 1.0 * global_minLen ** 2
        bw_cell_area = 0.7 * global_minLen ** 2
        web_consolidate_tol = 0.5 * global_minLen

        # ===================MESH SEGMENT
        for j, seg in enumerate(reversed(self.SegmentLst)):
            self.mesh.extend(seg.mesh_layers(self.SegmentLst, global_minLen, self.WebLst, display=self.display, l0=self.refL))

        # ===================MESH CORE
        if self.config.flags["mesh_core"]:
            for j, seg in enumerate(reversed(self.SegmentLst)):
                # if seg.ID == 1:
                # core_cell_area = 1.6*global_minLen**2
                # print(core_cell_area)
                self.mesh.extend(seg.mesh_core(self.SegmentLst, self.WebLst, core_cell_area, display=self.display))



        # ===================consolidate mesh on web interface
        for web in self.WebLst:
            #print web.ID,  'Left:', SegmentLst[web.ID].ID, 'Right:', SegmentLst[web.ID+1].ID,
            print('STATUS:\t Consolidate Mesh on Web Interface', web.ID)  
            (web.wl_nodes, web.wl_cells) = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID].cells, web.BSplineLst)            
            (web.wr_nodes, web.wr_cells) = grab_nodes_of_cells_on_BSplineLst(self.SegmentLst[web.ID+1].cells, web.BSplineLst)

            if not web.wl_nodes or not web.wl_cells or not web.wr_nodes or not web.wr_cells:  # in case there was no mesh in a segment
                print('STATUS:\t No mesh on Web Interface ' + str(web.ID) + ' to be consolodated.')

                # This message gets printed in cases where the web doesn't get
                # meshed because there is no isotropic layer.
                # However, may also get printed in other cases.
                # There are no other warnings printed about not having
                # an isotropic core on webs.
                print('STATUS:\t If web does not appear, ensure that webs have isotropic core layer.')
            else:
                newcells = consolidate_mesh_on_web(web, web_consolidate_tol, self.display)
                self.mesh.extend(newcells)


        print("STATUS:\t Splitting Quads into Trias")
        tmp = []
        for c in self.mesh:
            tmp.extend(c.split_quads())
        self.mesh = tmp
        
        # invert nodes list of all cell to make sure they are counterclockwise for vabs in the right coordinate system!
        for c in self.mesh:
            if c.orientation == False:
                c.invert_nodes()
        (self.mesh, nodes) = sort_and_reassignID(self.mesh)
        
        # Mesh cleanup
        # 1. Identify all nodes with exact repeated coordinates
        # (tolerance 1e-14)
        # 2. Keep only the lower number node and replace all cells node entries
        # as needed
        # Redo the sort and reassignID call
        
        # Find number of nodes
        n_nodes = 0
        
        for cell_i in self.mesh:
            n_nodes = np.maximum(n_nodes, np.max([n.id for n in cell_i.nodes]))

        node_coords = np.full((n_nodes+1, 2), np.nan)
        node_list = (n_nodes+1)* [None]
        
        for ind,cell_i in enumerate(self.mesh):
            for n in cell_i.nodes:
                node_coords[n.id] = [n.Pnt2d.X(), n.Pnt2d.Y()]
                node_list[n.id] = n
        
        same_coords = (n_nodes+1)* [[]]
        
        for ind in range(node_coords.shape[0]):
            
            same_coords[ind] = np.where(np.linalg.norm(node_coords
                                       - node_coords[ind], axis=1) < 1e-14)[0]
        
        reduced_sets = [group for group in same_coords if group.shape[0]>1]
        
        node_sets = len(reduced_sets)*[None]
        set_inds = 0
                
        for curr in reduced_sets:
            
            added = False
            
            for i in range(set_inds):
                if np.intersect1d(curr, node_sets[i]).shape[0] > 0:
                    added=True
                    node_sets[i] = np.unique(np.hstack((node_sets[i], curr)))
            
            if not added:
                node_sets[set_inds] = curr
                set_inds += 1
                
        node_sets = node_sets[:set_inds]
        # WARNING: It could be possible that a node is in multiple node sets.
        # Using a large enough tolerance should limit this risk

        for cell in self.mesh:
            for i,n in enumerate(cell.nodes):
                
                for check_set in node_sets:
                    if n.id in check_set:
                        cell.nodes[i] = node_list[check_set[0]]
        
        # remove any cell that has repeated nodes
        remove_inds = []
        
        for ind,cell in enumerate(self.mesh):
            nodes = [n.id for n in cell.nodes]
            
            _,counts = np.unique(nodes, return_counts=True)
            if counts.max() > 1:
                remove_inds += [ind]
        
        if len(remove_inds) > 0:
            print("Removing Cells with repeated nodes.")
            
            self.mesh = [cell
                         for i, cell in enumerate(self.mesh)
                         if i not in remove_inds]
        
        (self.mesh, nodes) = sort_and_reassignID(self.mesh)
        
        return None

    def cbm_custom_mesh(self, nodes, cells, materials, split_quads=True,
                        theta_11=None, theta_3=None):
        """
        Give a custom mesh to the section model.

        Parameters
        ----------
        nodes : (N, 2) numpy.ndarray
            Coordinates of each node. First column is x, second is y.
        cells : (M, 4) numpy.ndarray
            List of nodes for each element.
            Element orientation is set based on the vector between nodes
            indexed 1 and 2.
        materials : length N list
            Material for each cell.
        split_quads : bool, optional
            Flag for if quad elements should be split into triangles after
            reading the custom mesh.
        theta_3 : float, optional
            Value for fiber orientation angle to be passed down into SONATA
            and ANBA. If None, then zero is passed down. Units are degrees.
            The default value is None.

        Returns
        -------
        None.
        
        Notes
        -----
        
        Cannot currently set the ply orientation except by ordering the
        cell nodes appropriately for the cell orientation.

        """
        
        if theta_3 is None:
            theta_3 = 0
        
        # Generate node list, do not need to explicitly save since each cell
        # saves its own nodes.
        node_list = nodes.shape[0] * [None]
        
        for ind in range(nodes.shape[0]):
            
            node_list[ind] = Node(gp_Pnt2d(nodes[ind, 0], nodes[ind, 1]),
                                  ['Custom Layer', ind, 0])
            
        self.mesh = cells.shape[0] * [None]
        
        for ind in range(cells.shape[0]):
            
            c = Cell([node_list[nind] for nind in cells[ind]])
            
            if c.orientation == False:
                c.invert_nodes()
            
            c.calc_theta_1()
            c.theta_3 = theta_3
            c.MatID = int(materials[ind])
            c.structured = True
            
            if theta_11 is not None:
                c.theta_1[0] = theta_11[ind]
            
            self.mesh[ind] = c
        
        if split_quads:
            print("STATUS:\t Splitting Quads into Trias")
            tmp = []
            for c in self.mesh:
                tmp.extend(c.split_quads())
            self.mesh = tmp
        
        return None

    def cbm_run_anbax(self):
        """interface method to run the solver anbax from marco.morandini 
        
        Notes
        ----------
        To be defined.

        """



        self.mesh, nodes = sort_and_reassignID(self.mesh)

        # # plot before conversion
        # x_coord_sonata = np.zeros(len(nodes))
        # y_coord_sonata = np.zeros(len(nodes))
        # for i in range(len(nodes)):
        #     x_coord_sonata[i] = nodes[i].coordinates[0]  # x1
        #     y_coord_sonata[i] = nodes[i].coordinates[1]  # x2
        #
        # plt.plot(x_coord_sonata, y_coord_sonata)



        try:
            (mesh, matLibrary, materials, plane_orientations, fiber_orientations, maxE) = build_dolfin_mesh(self.mesh, nodes, self.materials)
        except:
            print('\n')
            print('==========================================\n\n')
            print('Error, Anba4 wrapper called, likely ')
            print('Anba4 _or_ Dolfin are not installed\n\n')
            print('==========================================\n\n')


        #TBD: pass it to anbax and run it!
        anba = anbax(mesh, 1, matLibrary, materials, plane_orientations, fiber_orientations, maxE)
        tmp_TS = anba.compute().getValues(range(6),range(6))    # get stiffness matrix
        tmp_MM = anba.inertia().getValues(range(6),range(6))    # get mass matrix

        # Define transformation T (from ANBA to SONATA/VABS coordinates)
        B = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        T = np.dot(np.identity(3), np.linalg.inv(B))

        self.BeamProperties = BeamSectionalProps()
        self.BeamProperties.TS = trsf_sixbysix(tmp_TS, T)
        self.BeamProperties.MM = trsf_sixbysix(tmp_MM, T)

        # self.BeamProperties.Xm = np.array(ComputeMassCenter(self.BeamProperties.MM))  # mass center - is already allocated from mass matrix
        self.BeamProperties.Xt = np.array(ComputeTensionCenter(self.BeamProperties.TS)) # tension center
        self.BeamProperties.Xs = np.array(ComputeShearCenter(self.BeamProperties.TS))   # shear center


        # --- Stress & Strain recovery --- #
        if  self.config.anbax_cfg.recover_flag == True:
            print("STATUS:\t Running ANBAX Stress & Strain Recovery:")
            [tmp_StressF_tran, tmp_StressF_M_tran, tmp_StrainF_tran, tmp_StrainF_M_tran] = \
                anbax_recovery(anba, len(self.mesh), self.config.anbax_cfg.F.tolist(), self.config.anbax_cfg.M.tolist(), self.config.anbax_cfg.voigt_convention, T)

            # ASSIGN stresses and strains to mesh elements:
            for i,c in enumerate(self.mesh):
                #                  [s_11[i],                   s_12[i],                   s_13[i],                   s_22[i],                   s_23[i],                   s_33[i]])
                c.stress =  Stress([tmp_StressF_tran[i,0,0],   tmp_StressF_tran[i,0,1],   tmp_StressF_tran[i,0,2],   tmp_StressF_tran[i,1,1],   tmp_StressF_tran[i,1,2],   tmp_StressF_tran[i,2,2]])
                c.stressM = Stress([tmp_StressF_M_tran[i,0,0], tmp_StressF_M_tran[i,0,1], tmp_StressF_M_tran[i,0,2], tmp_StressF_M_tran[i,1,1], tmp_StressF_M_tran[i,1,2], tmp_StressF_M_tran[i,2,2]])
                #                  [e_11[i],                   e_12[i],                   e_13[i],                   e_22[i],                   e_23[i],                   e_33[i]])
                c.strain =  Strain([tmp_StrainF_tran[i,0,0],   tmp_StrainF_tran[i,0,1],   tmp_StrainF_tran[i,0,2],   tmp_StrainF_tran[i,1,1],   tmp_StrainF_tran[i,1,2],   tmp_StrainF_tran[i,2,2]])
                c.strainM = Strain([tmp_StrainF_M_tran[i,0,0], tmp_StrainF_M_tran[i,0,1], tmp_StrainF_M_tran[i,0,2], tmp_StrainF_M_tran[i,1,1], tmp_StrainF_M_tran[i,1,2], tmp_StrainF_M_tran[i,2,2]])


        return


    def cbm_run_viscoelastic(self, test_elastic=True, test_tau0=True):
        """
        Calculation of viscoelastic 6x6 matrices at the section.

        Returns
        -------
        None.

        """

        self.mesh, nodes = sort_and_reassignID(self.mesh)

        try:
            (mesh, matLibrary, materials, plane_orientations,
             fiber_orientations, maxE) = build_dolfin_mesh(self.mesh,
                                                       nodes, self.materials)
        except:
            print('\n')
            print('==========================================\n\n')
            print('Error, Anba4 wrapper called, likely ')
            print('Anba4 _or_ Dolfin are not installed\n\n')
            print('==========================================\n\n')


        # Call ANBAX with baseline properties
        anba = anbax(mesh, 1, matLibrary, materials, plane_orientations,
                     fiber_orientations, maxE)
        
        tmp_TS = anba.compute().getValues(range(6),range(6))    # get stiffness matrix
        tmp_MM = anba.inertia().getValues(range(6),range(6))    # get mass matrix

        # Define transformation T (from ANBA to SONATA/VABS coordinates)
        B = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        T = np.dot(np.identity(3), np.linalg.inv(B))

        self.BeamProperties = BeamSectionalProps()
        self.BeamProperties.TS = trsf_sixbysix(tmp_TS, T)
        self.BeamProperties.MM = trsf_sixbysix(tmp_MM, T)

        # self.BeamProperties.Xm = np.array(ComputeMassCenter(self.BeamProperties.MM))  # mass center - is already allocated from mass matrix
        self.BeamProperties.Xt = np.array(ComputeTensionCenter(self.BeamProperties.TS)) # tension center
        self.BeamProperties.Xs = np.array(ComputeShearCenter(self.BeamProperties.TS))   # shear center

        # Recover the mapping from sectional Forces and Moments to Strain
        # in a non-invasive way from ANBA
        
        # 1. Initialize memory on each element to store the mapping
        # ASSIGN stresses and strains to mesh elements:
        for i,c in enumerate(self.mesh):
            
            c.fm_to_strain = np.zeros((6,6))
        
        # 2. Call ANBA looping over unit forces/moments
        for i in range(6):
            
            F = [0.0, 0.0, 0.0]
            M = [0.0, 0.0, 0.0]
            
            if i < 3:
                F[i] = 1.0
            else:
                M[i - 3] = 1.0
            
            # This ends up being potentially excessively slow since
            # it does a new calculation for each stress and strain field (4),
            # but only need the global coordinate strain field.
            [tmp_StressF_tran, tmp_StressF_M_tran, tmp_StrainF_tran, tmp_StrainF_M_tran] = \
                anbax_recovery(anba, len(self.mesh), F, M,
                               self.config.anbax_cfg.voigt_convention, T)
        
        
            # 3. Store Strain results in each case / element
            for j,c in enumerate(self.mesh):
            
                # This creates a Strain class object that clearly identifies
                # elasticity strain tensor components (epsilon) versus
                # engineering shear strain components (gamma).
                # While is is probably unneccesary for computation, it is done for
                # code clarity.
                curr_strain = Strain([tmp_StrainF_tran[j,0,0],
                                      tmp_StrainF_tran[j,0,1],
                                      tmp_StrainF_tran[j,0,2],
                                      tmp_StrainF_tran[j,1,1],
                                      tmp_StrainF_tran[j,1,2],
                                      tmp_StrainF_tran[j,2,2]])

                c.fm_to_strain[:, i] = np.array([curr_strain.epsilon11,
                                                 curr_strain.epsilon22,
                                                 curr_strain.epsilon33,
                                                 curr_strain.gamma23,
                                                 curr_strain.gamma13,
                                                 curr_strain.gamma12])

        # 4. Create a material dictionary for each time scale.
        #    The will loop over those time scales
        time_scale_list = []
        
        for MatID in self.materials:
            time_scale_list += self.materials[MatID]\
                                    .viscoelastic['time_scales_v'].tolist()
        
        time_scale_list = np.sort(np.unique(time_scale_list)).tolist()
        
        time_scale_mat_dicts = len(time_scale_list) * [None]
                
        for i, tau in enumerate(time_scale_list):
            
            curr_materials = OrderedDict()
            
            for MatID in self.materials:
                
                mat = self.materials[MatID]
                
                if hasattr(mat, 'viscoelastic'):
                    found_time_scale = tau in mat.viscoelastic['time_scales_v'].tolist()
                else:
                    found_time_scale = False
                
                if not found_time_scale:
                    # Set material as isotropic with no stiffness
                    # nu value is irrelevant since E=0.0
                    curr_dict = {'nu' : 0.0,
                                 'E'  : 0.0}
                    
                    curr_materials[MatID] = IsotropicMaterial(ID=MatID,
                                                                **curr_dict)
                else:
                    # find index of time scale
                    time_scale_ind = np.argmax(tau 
                                          == mat.viscoelastic['time_scales_v'])
                
                if found_time_scale and mat.orth == 0:
                    
                    curr_dict = {'E' : mat.viscoelastic['E_v'][time_scale_ind],
                                 'nu' : mat.nu,
                                 'name' : mat.name + ' time scale : {}'.format(tau)}

                    curr_materials[MatID] = IsotropicMaterial(ID=MatID,
                                                                **curr_dict)
                elif found_time_scale and mat.orth == 1:
                    
                    curr_dict = {'E_1' : mat.viscoelastic['E_1_v'][time_scale_ind],
                                 'E_2' : mat.viscoelastic['E_2_v'][time_scale_ind],
                                 'E_3' : mat.viscoelastic['E_3_v'][time_scale_ind],
                                 'G_12' : mat.viscoelastic['G_12_v'][time_scale_ind],
                                 'G_13' : mat.viscoelastic['G_13_v'][time_scale_ind],
                                 'G_23' : mat.viscoelastic['G_23_v'][time_scale_ind],
                                 'nu' : mat.nu.tolist(),
                                 'name' : mat.name + ' time scale : {}'.format(tau)}

                    curr_materials[MatID] = OrthotropicMaterial(ID=MatID,
                                                                flag_mat=False,
                                                                **curr_dict)
            time_scale_mat_dicts[i] = curr_materials

        # 5. calculate integrated element contributions
        # This mapping is just the partial product that maps force/moments
        # back to forces/moments (or time derivatives to be integrated)
        
        viscoelastic_6x6 = len(time_scale_list) * [None]
        
        for tau_ind, tau in enumerate(time_scale_list):
            
            material_dict = time_scale_mat_dicts[tau_ind]
            
            force_to_forcedot = np.zeros((6,6))
    
            ze = np.zeros((6,6))
            ze[0, -1] = 1.0 # shear stress x=2 -> shear force x
            ze[1, -2] = 1.0 # shear stress y=3 -> shear force y
            ze[2, 0] = 1.0 # axial stress -> axial force
    
            for i,c in enumerate(self.mesh):
    
                cxy = c.center
    
                ze[3, 0] = cxy[1] # axial stress -> moment around x
                ze[4, 0] = -cxy[0] # axial stress -> moment around y
    
                ze[-1, -2] = cxy[0] # shear zy (13) stress -> moment around z/torsion
                ze[-1, -1] = -cxy[1] # shear zx (12) stress -> moment around z/torsion
    
                De = material_dict[c.MatID].rotated_constitutive_tensor(
                                                    c.theta_1[0], c.theta_3)
    
                force_to_forcedot += c.area * (ze @ De @ c.fm_to_strain)
                
                
            viscoelastic_6x6[tau_ind] = trsf_sixbysix(force_to_forcedot @ tmp_TS, T)

        if test_elastic:
            
            force_to_forcedot = np.zeros((6,6))
    
            ze = np.zeros((6,6))
            ze[0, -1] = 1.0 # shear stress x=2 -> shear force x
            ze[1, -2] = 1.0 # shear stress y=3 -> shear force y
            ze[2, 0] = 1.0 # axial stress -> axial force
    
            for i,c in enumerate(self.mesh):
    
                cxy = c.center
    
                ze[3, 0] = cxy[1] # axial stress -> moment around x
                ze[4, 0] = -cxy[0] # axial stress -> moment around y
    
                ze[-1, -2] = cxy[0] # shear zy (13) stress -> moment around z/torsion
                ze[-1, -1] = -cxy[1] # shear zx (12) stress -> moment around z/torsion
    
                De = self.materials[c.MatID].rotated_constitutive_tensor(
                                                         c.theta_1[0], c.theta_3)
    
                force_to_forcedot += c.area * (ze @ De @ c.fm_to_strain)
            
            error = np.linalg.norm(force_to_forcedot - np.eye(6))
            
            print('Error in recovering elastic 6x6 with mappings is: {:.3e}'
                  .format(error))

        if test_tau0:
            # Test that all of the viscoelastic matrices add up to the 
            # default calculated matrix
            
            sum_6x6 = np.zeros((6,6))
            
            for mat6x6 in viscoelastic_6x6:
                sum_6x6 += mat6x6
            
            error = np.linalg.norm(self.BeamProperties.TS - sum_6x6)
            
            mag = np.linalg.norm(self.BeamProperties.TS)
            
            print('Adding all time scale 6x6 v. baseline error: '
                  + 'absolute: {:.3e}, relative: {:.3e}'
                  .format(error, error/mag))
            
            print('This error should be small if the (sum of the elastic'
                  + ' modulus and shear modulus at all time scales) equals'
                  + ' the reference elastic and shear moduli.')

        self.BeamProperties.TSv = viscoelastic_6x6
        self.BeamProperties.tau = time_scale_list
        
        return viscoelastic_6x6


    def cbm_exp_stress_strain_map(self, station_ind, station_pos, **kwargs):
        """
        Calculation and save of the map from the 6x6 applied forces to the
        internal stress and strain at each element
        
        Parameters
        ----------
        station_ind : int
            Index of station along the blade. Used in the filename of the
            output.
        station_pos : float
            Fraction of the blade length. Included in output for verification
            when the outputs are used.
        output_folder : str, optional
            Folder to output mapping files to.
        **kwargs : TYPE
            Ignored.

        Returns
        -------
        None.
        
        Notes
        -----
        Stress and strain maps are for the SONATA force coordinates to the
        local (material coordinate direction) stress and strain.
        
        Outputs are (6,6,n_elem). First dimension is stress or strain in order
        [11, 22, 33, 23, 13, 12]. Second dimension is the internal forces then
        moments. Third dimension is the elements.
        
        Outputs are of the form 
        'blade_station{station_ind}_stress_strain_map.npz'
        
        Outputs are engineering shear strain (not elasticity tensor components)

        """

        self.mesh, nodes = sort_and_reassignID(self.mesh)

        try:
            (mesh, matLibrary, materials, plane_orientations,
             fiber_orientations, maxE) = build_dolfin_mesh(self.mesh,
                                                       nodes, self.materials)
        except:
            print('\n')
            print('==========================================\n\n')
            print('Error, Anba4 wrapper called, likely ')
            print('Anba4 _or_ Dolfin are not installed\n\n')
            print('==========================================\n\n')


        # Call ANBAX with baseline properties
        anba = anbax(mesh, 1, matLibrary, materials, plane_orientations,
                     fiber_orientations, maxE)
        
        tmp_TS = anba.compute().getValues(range(6),range(6))    # get stiffness matrix
        tmp_MM = anba.inertia().getValues(range(6),range(6))    # get mass matrix

        # Define transformation T (from ANBA to SONATA/VABS coordinates)
        B = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        T = np.dot(np.identity(3), np.linalg.inv(B))

        self.BeamProperties = BeamSectionalProps()
        self.BeamProperties.TS = trsf_sixbysix(tmp_TS, T)
        self.BeamProperties.MM = trsf_sixbysix(tmp_MM, T)

        # self.BeamProperties.Xm = np.array(ComputeMassCenter(self.BeamProperties.MM))  # mass center - is already allocated from mass matrix
        self.BeamProperties.Xt = np.array(ComputeTensionCenter(self.BeamProperties.TS)) # tension center
        self.BeamProperties.Xs = np.array(ComputeShearCenter(self.BeamProperties.TS))   # shear center

        # Recover the mapping from sectional Forces and Moments to Strain
        # in a non-invasive way from ANBA
        
            
        fc_to_strain_m = np.zeros((6,6,len(self.mesh)))
        fc_to_stress_m = np.zeros((6,6,len(self.mesh)))
        elem_areas = np.zeros((len(self.mesh)))
        elem_materials = np.zeros((len(self.mesh)))
        elem_cxy = np.zeros((len(self.mesh), 2))
        
        # 1. Initialize memory on each element to store the mapping
        # ASSIGN stresses and strains to mesh elements:
        for i,c in enumerate(self.mesh):
            
            elem_areas[i] = c.area
            elem_materials[i] = c.MatID
            elem_cxy[i, :] = c.center

        # 2. Call ANBA looping over unit forces/moments
        
        # Mapping should match `beam_struct_eval`
        external_to_internal_ind = [2, 0, 1]
        for i in range(6):
            
            F = [0.0, 0.0, 0.0]
            M = [0.0, 0.0, 0.0]
            
            if i < 3:
                F[external_to_internal_ind[i]] = 1.0
            else:
                M[external_to_internal_ind[i - 3]] = 1.0
            
            # This ends up being potentially excessively slow since
            # it does a new calculation for each stress and strain field (4),
            # but only need the material coordinate strain field.
            [tmp_StressF_tran, tmp_StressF_M_tran, tmp_StrainF_tran, tmp_StrainF_M_tran] = \
                anbax_recovery(anba, len(self.mesh), F, M,
                               self.config.anbax_cfg.voigt_convention, T)
        
        
            # 3. Store Strain results in each case / element
            for j,c in enumerate(self.mesh):
            
                # This creates a Strain class object that clearly identifies
                # elasticity strain tensor components (epsilon) versus
                # engineering shear strain components (gamma).
                # While is is probably unneccesary for computation, it is done for
                # code clarity.
                curr_strain = Strain([tmp_StrainF_M_tran[j,0,0],
                                      tmp_StrainF_M_tran[j,0,1],
                                      tmp_StrainF_M_tran[j,0,2],
                                      tmp_StrainF_M_tran[j,1,1],
                                      tmp_StrainF_M_tran[j,1,2],
                                      tmp_StrainF_M_tran[j,2,2]])
                
                
                curr_stress = Stress([tmp_StressF_M_tran[j,0,0],
                                      tmp_StressF_M_tran[j,0,1],
                                      tmp_StressF_M_tran[j,0,2],
                                      tmp_StressF_M_tran[j,1,1],
                                      tmp_StressF_M_tran[j,1,2],
                                      tmp_StressF_M_tran[j,2,2]])

                fc_to_strain_m[:, i, j] = np.array([curr_strain.epsilon11,
                                                    curr_strain.epsilon22,
                                                    curr_strain.epsilon33,
                                                    curr_strain.gamma23,
                                                    curr_strain.gamma13,
                                                    curr_strain.gamma12])

                fc_to_stress_m[:, i, j] = np.array([curr_stress.sigma11,
                                                    curr_stress.sigma22,
                                                    curr_stress.sigma33,
                                                    curr_stress.sigma23,
                                                    curr_stress.sigma13,
                                                    curr_stress.sigma12])
                
        # Save the material names as strings to be output
        material_names = (np.max([mat for mat in self.materials]) + 1)*['None']
        
        for mat in self.materials:
            material_names[mat] = self.materials[mat].name
            
        
        if 'output_folder' in kwargs.keys():
            folder = kwargs['output_folder']
        else:
            folder = 'stress-map'
        
        # Output File with Data
        fname = os.path.join(folder,
             'blade_station{:04d}_stress_strain_map.npz'.format(station_ind))
        
        os.makedirs(folder, exist_ok=True)
        
        np.savez(fname,
                 station_pos=station_pos,
                 station_ind=station_ind,
                 fc_to_strain_m=fc_to_strain_m,
                 fc_to_stress_m=fc_to_stress_m,
                 elem_areas=elem_areas,
                 elem_materials=elem_materials,
                 material_names=np.asarray(material_names),
                 elem_cxy=elem_cxy,
                 allow_pickle=False)

        return

    def cbm_exp_BeamDyn_beamprops(self, Theta=0, solver="vabs"):
        """ 
        Converts the Beam Properties of CBM to the correct coordinate System of
        BeamDyn and returns the 6x6 Stiffness matrix, the 6x6 MassMatrix.
        
        The geometry of the blade is defined by key-point coordinates and initial
        twist angles (in units of degree) in the blade local coordinate system
        (IEC standard blade system where Zr is along blade axis from root to
        tip, Xr directs normally toward the suction side, and Yr directs 
        normally toward the trailing edge).
        https://openfast.readthedocs.io/en/master/source/user/beamdyn/input_files.html
        
        Parameters
        ----------
        Theta: float, optional
            is the angle of rotation of the coordinate system in "radians"
        solver: str, optional
        
        Returns
        ----------
            tuple of arrays
            (6x6 StiffnessMatrix, 6x6MassMatrix)
            
            
        Notes:
        ----------
        - Following the station location parameter η, there are two 
        6×6 matrices providing the structural and inertial properties for this
        cross-section. First is the stiffness matrix and then the mass matrix. 
        We note that these matrices are defined in a local coordinate system 
        along the blade axis with Zl directing toward the unit tangent vector 
        of the blade reference axis.
        - Does this create an oblique cross-section!?
        
        
        """
        if solver == "vabs" or solver == "anbax":
            if Theta != 0:
                tmp_bp = self.BeamProperties.rotate(Theta)
            else:
                tmp_bp = self.BeamProperties

        else:
            print("Check solver for BeamDyn Beam Property input.")

        tmp_bp = copy.deepcopy(tmp_bp)

        # transform to BeamDYN Coordinates
        B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
        T = np.dot(np.identity(3), np.linalg.inv(B))

        tmp_TS = trsf_sixbysix(tmp_bp.TS, T)
        tmp_MM = trsf_sixbysix(tmp_bp.MM, T)
        return (tmp_TS, tmp_MM)



    def cbm_exp_dymore_beamprops(self, eta, Theta=0, solver="vabs", units={"mass": "kg", "length": "m", "force": "N"}):
        """
        Converts the Units of CBM to DYMORE/PYMORE/MARC units and returns the 
        array of the beamproperties with Massterms(6), Stiffness(21), 
        damping(1) and curvilinear coordinate(1)

        Parameters
        ----------
        
        eta : float, 
            is the beam curvilinear coordinate of the beam from 0 to 1. 
        
        Theta: float
            is the angle of rotation of the coordinate system in "radians"

        Returns
        ----------
        arr : ndarray
            [Massterms(6) (m00, mEta2, mEta3, m33, m23, m22) 
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]
            
            
        Notes
        ----------
        - Unit Convertion takes sooo much time. Commented out for now!
        
        """
        if solver == "vabs" or solver == "anbax":
            if Theta != 0:
                tmp_bp = self.BeamProperties.rotate(Theta)
            else:
                tmp_bp = self.BeamProperties

        else:
            print("Check solver for Dymore Beam Property input.")


        MM = tmp_bp.MM
        MASS = np.array([MM[0, 0], MM[2, 3], MM[0, 4], MM[5, 5], MM[4, 5], MM[4, 4]])
        STIFF = tmp_bp.TS[np.tril_indices(6)[1], np.tril_indices(6)[0]]
        mu = 0.0
        return np.hstack((MASS, STIFF, mu, eta))


    def cbm_post_2dmesh(self, attribute="MatID", title="NOTITLE", **kw):
        """
        CBM Postprocessing method that displays the mesh with matplotlib.
        
        Parameters
        ----------
        attribute : string, optional
            Uses the string to look for the cell attributes. 
            The default attribute is MatID. Possible other attributes can be 
            fiber orientation (theta_3) or strains and stresses. 
            If BeamProperties are already calculated by VABS or something
            similar, elastic-axis, center-of-gravity... are displayed.
        title : string, optional
            Title to be placed over the plot.
        **kw : keyword arguments, optional
            are passed to the lower "plot_cells" function. Such options are: 
            VABSProperties=None, title='None', plotTheta11=False, 
            plotDisplacement=False, savepath
            
        Returns
        ----------
        (fig, ax) : tuple
            figure and axis handler are returned by the method
        
        Examples
        ----------
        >>> job.cbm_post_2dmesh(title='Hello World!', attribute='theta_3', plotTheta11=True)

        """
        mesh, nodes = sort_and_reassignID(self.mesh)
        fig, ax = plot_cells(self.mesh, nodes, attribute, self.materials, self.BeamProperties, title, **kw)
        return fig, ax

#%%############################################################################
#                           M    A    I    N                                  #
###############################################################################
if __name__ == "__main__":
    plt.close("all")
    fname = "jobs/debug/issue20/sec_config.yml"
    fname = "jobs/VariSpeed/uh60a_cbm_advanced/sec_config_R2000.yml"
    # fname = 'jobs/AREA/R250/sec_config.yml'
    # fname = 'jobs/PBortolotti/sec_config.yml'
    config = CBMConfig(fname)

    job = CBM(config)

    job.cbm_gen_topo()
    job.cbm_gen_mesh(split_quads=True)

    job.cbm_review_mesh()
    job.cbm_run_vabs(rm_vabfiles=False)
    # AnbaBeamProperties = job.cbm_run_anbax()

    # job.cbm_post_2dmesh(title='Hello World!')
    job.cbm_post_3dtopo()
#    job.config.vabs_cfg.recover_flag = 1
#    job.config.vabs_cfg.M = [0,2000e4,0]
