# -*- coding: utf-8 -*-
"""Defines the Crosssectional Beam Model (CBM) class
Created on Wed Jan 03 13:56:37 2018
@author: TPflumm

https://numpydoc.readthedocs.io/en/latest/format.html
 """

# Core Library modules
import copy
import math
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
# from SONATA.vabs.classVABSConfig import VABSConfig

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
        self.SegmentLst[0].build_layers(l0=self.refL, **kwargs)
        self.SegmentLst[0].determine_final_boundary()

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
                seg.build_layers(self.WebLst, self.SegmentLst[0], l0=self.refL)
                seg.determine_final_boundary(self.WebLst, self.SegmentLst[0])
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
            # mesh,nodes = sort_and_reassignID(mesh)

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
                print('STATUS:\t No mesh on Web Interface ' + str(web.ID) + ' to be consolodated')
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
            print('Error, Anba4 wrapper called, but ')
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
