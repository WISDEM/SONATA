#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:38:19 2018

@author: Tobias Pflumm
"""
# Third party modules
import matplotlib.pyplot as plt
import numpy as np
import yaml
# from jsonschema import validate
from OCC.Core.gp import (gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pln,
                         gp_Pnt, gp_Pnt2d, gp_Trsf, gp_Vec,)
from scipy.interpolate import interp1d

# First party modules
from SONATA.cbm.classCBM import CBM
from SONATA.cbm.classCBMConfig import CBMConfig
from SONATA.cbm.display.display_utils import (display_Ax2,
                                              display_cbm_SegmentLst,
                                              display_config,)

from SONATA.cbm.topo.BSplineLst_utils import (BSplineLst_from_dct,
                                              set_BSplineLst_to_Origin2,)
from SONATA.cbm.topo.to3d import bsplinelst_to3d, pnt_to3d, vec_to3d
from SONATA.cbm.topo.utils import (Array_to_PntLst, PntLst_to_npArray,)
from SONATA.cbm.topo.wire_utils import (equidistant_Points_on_wire,)
from SONATA.classAirfoil import Airfoil
from SONATA.classComponent import Component
from SONATA.classMaterial import read_materials
from SONATA.utl.blade_utl import (array_pln_intersect, check_uniformity,
                                  interp_airfoil_position, interp_loads,
                                  make_loft,)
from SONATA.utl.converter_WT import converter_WT
from SONATA.utl.interpBSplineLst import interpBSplineLst
from SONATA.utl.plot import plot_beam_properties
from SONATA.utl.trsf import trsf_af_to_blfr, trsf_blfr_to_cbm
from SONATA.anbax.classANBAXConfig import ANBAXConfig



# from SONATA.airconics_blade_cad.blade_cst import blade_cst
# import SONATA.airconics_blade_cad.airconics.liftingSurface as liftingSurface

def rotate(xo, yo, xp, yp, angle):
    ## Rotate a point clockwise by a given angle around a given origin.
    # angle *= -1.
    qx = xo + np.cos(angle) * (xp - xo) - np.sin(angle) * (yp - yo)
    qy = yo + np.sin(angle) * (xp - xo) + np.cos(angle) * (yp - yo)
    return qx, qy

class Blade(Component):
    """
    SONATA Blade component object.
    
    Attributes
    ----------                      
    coordinates :  ndarray
        Describes the axis LE coordinates in meters along the span.
        nparray([[grid, x, y, z]]).
        The grid represents the nondimensional x position along the Blade from 
        0 to 1
    
    chord : ndarray
        Describes the blades chord lenght in meters in spanwise direction. 
        nparray([[grid, chord]]) 

    twist : ndarray
        Describes the blades twist angles in !radians! in spanwise direction. 
        nparray([[grid, twist]]) 

    pitch_axis : ndarray
        Describes the blades pitch-axis location in 1/chord lengths from the 
        leading edge. nparray([[grid, pitch_axis]]) 

    airfoils : ndarray
        array of grid location and airfoil instance 
        nparray([[grid, airfoil instance]],dtype = object)
        
    sections : ndarray
        array of CBM cross-sections 
        nparray([[grid, CBM instance]],dtype = object)
        
    beam_properties : ndarray
        array of grid location and VABSSectionalProp instance
        nparray([[grid, beam_properties]],dtype = object)
              
        
    Methods
    -------
    blade_matrix : ndarray
        Summons all the blades global properties in one array
        nparray([[grid, x, y, z, chord, twist, pitch_axis,....]])
    
    
    Notes
    --------
    Units: meter (m), Newton (N), kilogramm (kg), degree (deg), Kelvin (K),



    See Also
    --------
    Component,
    

    ToDo
    -----
    - Include the possibity to rotate the beam_properties non-twisted frame. 
        Default is the twisted frame
    -
    

    Examples
    --------
    Initialize Blade Instance:
    
    >>> job = Blade(name='UH-60A_adv')
    
    >>> job.read_yaml(yml.get('components').get('blade'), airfoils, materials)
    
    >>> job.blade_gen_section()
    >>> job.blade_run_anbax()
    >>> job.blade_plot_sections()
    >>> job.blade_post_3dtopo(flag_lft = True, flag_topo = True)

    """

    __slots__ = (
        "blade_ref_axis",
        "chord",
        "twist",
        "curvature",
        "pitch_axis",
        "airfoils",
        "sections",
        "beam_properties",
        "beam_ref_axis",
        "f_chord",
        "f_twist",
        "materials",
        "blade_ref_axis_BSplineLst",
        "f_blade_ref_axis",
        "beam_ref_axis_BSplineLst",
        "f_beam_ref_axis",
        "f_pa",
        "f_curvature_k1",
        "anba_beam_properties",
        "wopwop_bsplinelst",
        "wopwop_pnts",
        "wopwop_vecs",
        "display",
        "start_display",
        "add_menu",
        "add_function_to_menu",
        "yml",
        "loft"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.beam_properties = None
        self.loft=None
        
        if 'filename' in kwargs:
            filename = kwargs.get('filename')
            with open(filename, 'r') as myfile:
                inputs  = myfile.read()
                yml = yaml.load(inputs, Loader = yaml.FullLoader)
                self.yml = yml

            
            airfoils = [Airfoil(af) for af in yml.get('airfoils')]
            self.materials = read_materials(yml.get('materials'))
            
            self.read_yaml(yml.get('components').get('blade'), airfoils, **kwargs)

            
#    def __repr__(self):
#        """__repr__ is the built-in function used to compute the "official" 
#        string reputation of an object, """
#        return 'Blade: '+ str(self.name)
    
    def _read_ref_axes(self, yml_ra, flag_ref_axes_wt=False, c2_axis=False, tmp_chord = [], tmp_pa = [], tmp_tw=[]):
        """
        reads and determines interpolates function for the reference axis of 
        the blade
        
        Parameters
        ----------
        yml_ra : dict
            yaml style dict data describes the referenceaxis with non-dim
            grid stations and x,y,z values

        Returns
        -------
        BSplineLst : list of OCC.GeomBSplines
            DESCRIPTION.
        f_ra : function
            BSplineLst interpolation function
        tmp_ra : np.ndarray
            DESCRIPTION.

        """
        tmp_ra = {}

        if flag_ref_axes_wt:
            # adapt reference axis provided in yaml file to match with SONATA (equiv. rotorcraft) format
            # x_SONATA equiv. to z_wind
            # y_SONATA equiv. to -y_wind
            # z_SONATA equiv. to x_wind
            tmp_ra['x'] = np.asarray((yml_ra.get('z').get('grid'), yml_ra.get('z').get('values'))).T
            tmp_ra['y'] = np.asarray((yml_ra.get('y').get('grid'), np.negative(yml_ra.get('y').get('values')))).T
            tmp_ra['z'] = np.asarray((yml_ra.get('x').get('grid'), yml_ra.get('x').get('values'))).T
        else:
            tmp_ra['x'] = np.asarray((yml_ra.get('x').get('grid'),yml_ra.get('x').get('values'))).T
            tmp_ra['y'] = np.asarray((yml_ra.get('y').get('grid'),yml_ra.get('y').get('values'))).T
            tmp_ra['z'] = np.asarray((yml_ra.get('z').get('grid'),yml_ra.get('z').get('values'))).T
        
        f_ref_axis_x = interp1d(tmp_ra['x'][:,0], tmp_ra['x'][:,1], bounds_error=False, fill_value='extrapolate')
        f_ref_axis_y = interp1d(tmp_ra['y'][:,0], tmp_ra['y'][:,1], bounds_error=False, fill_value='extrapolate')
        f_ref_axis_z = interp1d(tmp_ra['z'][:,0], tmp_ra['z'][:,1], bounds_error=False, fill_value='extrapolate')
        
        x_blra = np.unique(np.sort(np.hstack((tmp_ra['x'][:,0], tmp_ra['y'][:,0], tmp_ra['z'][:,0]))))
        tmp_ra = np.vstack((x_blra, f_ref_axis_x(x_blra), f_ref_axis_y(x_blra), f_ref_axis_z(x_blra))).T
        
        if c2_axis:
            f_chord = interp1d(tmp_chord[:,0], tmp_chord[:,1], bounds_error=False, fill_value='extrapolate')
            chord = f_chord(x_blra)
            f_pitch_axis = interp1d(tmp_pa[:,0], tmp_pa[:,1], bounds_error=False, fill_value='extrapolate')
            pitch_axis = f_pitch_axis(x_blra)
            f_tw = interp1d(tmp_tw[:,0], tmp_tw[:,1], bounds_error=False, fill_value='extrapolate')
            twist_rad = f_tw(x_blra)
            # Get the absolute offset between mid chord and pitch axis (rotation center)
            ch_offset = chord * (0.5 - pitch_axis)
            # Rotate it by the twist
            z , y = rotate(0., 0., 0., ch_offset, -twist_rad)
            tmp_ra[:,2] -= y
            tmp_ra[:,3] += z

        if check_uniformity(tmp_ra[:, 0], tmp_ra[:, 1]) == False:
            print("WARNING:\t The blade beference axis is not uniformly defined along x")

        # print(tmp_ra[:,1:])
        BSplineLst = BSplineLst_from_dct(tmp_ra[:, 1:], angular_deflection=5, twoD=False)
        f_ra = interpBSplineLst(BSplineLst, tmp_ra[:, 0], tmp_ra[:, 1])
        return (BSplineLst, f_ra, tmp_ra)

    def _get_local_Ax2(self, x):
        """
        

        Parameters
        ----------
        x : float
            non-dimensional grid location

        Returns
        -------
        local_Ax2 : OCC.gp_Ax2
            return the gp_AX2 coordinatesystem

        """
        # interpolate blade_ref_axis
        res, resCoords = self.f_beam_ref_axis.interpolate(x)
        # print(res)
        p = gp_Pnt()
        vx = gp_Vec()
        v2 = gp_Vec()
        
        #determine local the local cbm coordinate system Ax2
        self.beam_ref_axis_BSplineLst[int(resCoords[0,0])].D2(resCoords[0,1],p,vx,v2)
        vz = gp_Vec(-vx.Z(),0,vx.X()).Normalized()
        tmp_Ax2 = gp_Ax2(p, gp_Dir(vz), gp_Dir(vx))
        local_Ax2 = tmp_Ax2.Rotated(gp_Ax1(p, gp_Dir(vx)), float(self.f_twist(x)))
        return local_Ax2

    def _interpolate_cbm_boundary(self, x, fs=1.1, nPoints=4000):
        """
        interpolates a cbm boundary BSplineLst from the blade definition at a 
        certain grid station. Following the procedure: 
        Determine all important neighboring airfoil positions
        discretize all airfoils equidistantly with the same number of Points. 
        Use these Points to performe a plane_line_intersection with the local 
        coordinate system Ax2. Find the correct intersection and extrapolate if
        necessary over the blade boundaries.
        Transfer the points to the local cbm frame and performe a BSpline 
        interpolation.        

        Parameters
        -------
        x : float
            nondimensional grid location
        
        Returns
        -------
        BoundaryBSplineLst : BSplineLst
            of the Boundary for the CBM Crosssection in the cbm frame

        ToDo
        -------
        - Use equidistant_Points_on_BSplineLst instead of equidistant_Points_on_wire 
            to capture corners
            
        """
        ax2 = self._get_local_Ax2(x)

        a = float(self.f_chord(x)) * float(self.f_pa(x))
        b = float(self.f_chord(x)) * (1 - float(self.f_pa(x)))
        beta = self.Ax2.Angle(self._get_local_Ax2(x))
        x0 = x - (np.sin(beta) * a * fs / self.f_blade_ref_axis.interpolate(1.0)[0][0, 0])
        x1 = x + (np.sin(beta) * b * fs / self.f_blade_ref_axis.interpolate(1.0)[0][0, 0])

        # select all airfoil in the interval between x0 < x1 and their closest neighbors
        idx0 = np.searchsorted(self.airfoils[:, 0], x0, side="left") - 1
        idx1 = np.searchsorted(self.airfoils[:, 0], x1, side="right")

        if idx0 < 0:
            idx0 = 0

        afs = self.airfoils[idx0 : idx1 + 1]

        # transform airfoils from nondimensional coordinates to coordinates
        afs = self.airfoils[idx0 : idx1 + 1]
        wireframe = []
        tes = []
        for item in afs:
            xi = item[0]
            af = item[1]
            (wire, te_pnt) = af.trsf_to_blfr(self.f_blade_ref_axis.interpolate(xi)[0][0], float(self.f_pa(xi)), float(self.f_chord(xi)), float(self.f_twist(xi)))
            wireframe.append(wire)
            tes.append(te_pnt)

        if len(wireframe) > 1:
            tmp = []
            for w in wireframe:
                PntLst = equidistant_Points_on_wire(w, nPoints)
                tmp.append(PntLst_to_npArray(PntLst))
            array = np.asarray(tmp)
            te_array = np.expand_dims(PntLst_to_npArray(tes), axis=1)
            result = array_pln_intersect(array, ax2)
            te_res = array_pln_intersect(te_array, ax2)

        else:
            w = wireframe[0]
            PntLst = equidistant_Points_on_wire(w, nPoints)
            result = PntLst_to_npArray(PntLst)
            te_res = PntLst_to_npArray(tes)

        trsf = trsf_blfr_to_cbm(self.Ax2, ax2)
        PntLst = Array_to_PntLst(result)
        te_pnt = Array_to_PntLst(te_res)[0]
        # Going from blade frame to cbm frame
        PntLst = [p.Transformed(trsf) for p in PntLst]
        te_pnt = te_pnt.Transformed(trsf)

        array = PntLst_to_npArray(PntLst)
        # array = np.flipud(array)
        # print(array)
        BSplineLst = BSplineLst_from_dct(array[:, 0:2], angular_deflection=30, tol_interp=1e-6)
        BoundaryBSplineLst = set_BSplineLst_to_Origin2(BSplineLst, gp_Pnt2d(te_pnt.Coord()[0], te_pnt.Coord()[1]))

        return BoundaryBSplineLst

    def read_yaml(self, yml, airfoils, stations=None, npts=11, wt_flag=False, **kwargs):
        """
        reads the Beam or Blade dictionary
        generates the blade matrix and airfoil to represent all given 
        information at every grid point by interpolating the input data 
        and assign them to the class attribute twist, choord, coordinates
        and airfoil_positions with the first column representing the 
        non-dimensional radial location

        Parameters
        ----------
        airfoils : list
            Is the database of airfoils
        
        """
        self.name = self.yml.get('name')
        print('STATUS:\t Reading YAML Dictionary for Beam/Blade: %s' % (self.name))
        c2_axis = kwargs.get('flags',{}).get('c2_axis')
        #Read chord, twist and nondim. pitch axis location and create interpolation
        tmp_chord = np.asarray((yml.get('outer_shape_bem').get('chord').get('grid'),yml.get('outer_shape_bem').get('chord').get('values'))).T
        tmp_tw = np.asarray((yml.get('outer_shape_bem').get('twist').get('grid'),yml.get('outer_shape_bem').get('twist').get('values'))).T
        tmp_pa = np.asarray((yml.get('outer_shape_bem').get('pitch_axis').get('grid'),yml.get('outer_shape_bem').get('pitch_axis').get('values'))).T
        
        #Read blade & beam reference axis and create BSplineLst & interpolation instance
        (self.blade_ref_axis_BSplineLst, self.f_blade_ref_axis, tmp_blra) = self._read_ref_axes(yml.get('outer_shape_bem').get('reference_axis'), flag_ref_axes_wt=kwargs.get('flags', {}).get('flag_ref_axes_wt'), c2_axis=c2_axis, tmp_chord = tmp_chord, tmp_pa = tmp_pa, tmp_tw=tmp_tw)

        if not yml.get('outer_shape_bem').get('beam_reference_axis'):
            #  In case beam reference axis is not defined in yaml file, use identical coordinates for beam reference and reference axis
            (self.beam_ref_axis_BSplineLst, self.f_beam_ref_axis, tmp_bera) = self._read_ref_axes(yml.get('outer_shape_bem').get('reference_axis'), flag_ref_axes_wt=kwargs.get('flags', {}).get('flag_ref_axes_wt'), c2_axis=c2_axis, tmp_chord = tmp_chord, tmp_pa = tmp_pa, tmp_tw=tmp_tw)
        else:
            (self.beam_ref_axis_BSplineLst, self.f_beam_ref_axis, tmp_bera) = self._read_ref_axes(yml.get('outer_shape_bem').get('beam_reference_axis'), flag_ref_axes_wt=kwargs.get('flags', {}).get('flag_ref_axes_wt'))
        
        if c2_axis:
            tmp_pa[:,1]=0.5
        self.f_chord = interp1d(tmp_chord[:,0], tmp_chord[:,1], bounds_error=False, fill_value='extrapolate')
        self.f_twist = interp1d(tmp_tw[:,0], tmp_tw[:,1], bounds_error=False, fill_value='extrapolate')
        self.f_pa = interp1d(tmp_pa[:,0], tmp_pa[:,1], bounds_error=False, fill_value='extrapolate')
        
        #Read airfoil information 
        airfoil_position = (yml.get('outer_shape_bem').get('airfoil_position').get('grid'),yml.get('outer_shape_bem').get('airfoil_position').get('labels'))
        tmp = []
        for an in airfoil_position[1]:
            tmp.append(next((x for x in airfoils if x.name == an), None).id)
        arr = np.asarray([airfoil_position[0],tmp]).T

        #Read CBM Positions
        if kwargs.get('flags',{}).get('flag_wt_ontology'):
            if stations is not None:
                cs_pos = stations
            else:
                cs_pos = np.linspace(0.0, 1.0, npts)
        else:
            if stations is None:
                cs_pos = np.asarray([cs.get('position') for cs in yml.get('internal_structure_2d_fem').get('sections')])
            else:
                cs_pos = stations
            
        x = np.unique(np.sort(np.hstack((tmp_chord[:,0], tmp_tw[:,0],
                                         tmp_blra[:,0], tmp_bera[:,0],
                                         tmp_pa[:,0], arr[:,0], cs_pos))))

        self.airfoils = np.asarray([[x, interp_airfoil_position(airfoil_position, airfoils, x)] for x in x])
        self.blade_ref_axis = np.hstack((np.expand_dims(x, axis=1), self.f_blade_ref_axis.interpolate(x)[0]))
        self.beam_ref_axis = np.hstack((np.expand_dims(x, axis=1), self.f_beam_ref_axis.interpolate(x)[0]))
        self.chord = np.vstack((x, self.f_chord(x))).T
        self.twist = np.vstack((x, self.f_twist(x))).T
        self.pitch_axis = np.vstack((x, self.f_pa(x))).T
        self.f_curvature_k1 = interp1d(x, np.gradient(self.twist[:,1],self.beam_ref_axis[:,1]))  # determine twist per unit length, i.e. the twist gradient at a respective location



        #Generate CBMConfigs
        if kwargs.get('flags',{}).get('flag_wt_ontology'):
            cbmconfigs = converter_WT(self, cs_pos, yml, self.materials, mesh_resolution = kwargs.get('flags').get('mesh_resolution'))
            
        else:
            lst = [[cs.get("position"), CBMConfig(cs, self.materials)] for cs in yml.get("internal_structure_2d_fem").get("sections")]
            cbmconfigs = np.asarray(lst)



        #Generate CBMs
        tmp = []
        for x, cfg in cbmconfigs:
            print(self.name, x)
            # get local beam coordinate system, and local cbm_boundary
            tmp_Ax2 = self._get_local_Ax2(x)
            tmp_blra = self.f_beam_ref_axis.interpolate(x)[0][0]
            BoundaryBSplineLst = self._interpolate_cbm_boundary(x)
            cs_name = self.name + '_section_R'+ ("%.3f" % x).replace('.','')
            tmp.append([x, CBM(cfg, materials=self.materials, name=cs_name, Ax2=tmp_Ax2, BSplineLst=BoundaryBSplineLst)])
        self.sections = np.asarray(tmp)

        return None

    @property
    def blade_matrix(self):
        """
         getter method for the property blade_matrix to retrive the full
        information set of the class in one reduced array

        Returns
        -------
        np.ndarray
            blade matrix of bl_ra, chord, twist, pa, 

        """
        return np.column_stack((self.blade_ref_axis, self.chord[:, 1], self.twist[:, 1], self.pitch_axis[:, 1]))


    def blade_gen_section(self, topo_flag=True, mesh_flag=True, **kwargs):
        """
        generates and meshes all cross-sections of the blade

        Parameters
        ----------
        topo_flag : bool, optional
            If this flag is true the topology of each cross-section is 
            generated. The default is True.
        mesh_flag : bool, optional
            IF this flag is set true, the discretization of each cross-section 
            is generated if a topology is generated beforehand. 
            The default is True.
        **kwargs : TYPE
            keyword arguments can be passed down to the cbm_gen_mesh function

        Returns
        -------
        None.

        """
        for (x, cs) in self.sections:
            if topo_flag:
                print("STATUS:\t Building Section at grid location %s" % x)
                cs.cbm_gen_topo()
            if mesh_flag:
                print("STATUS:\t Meshing Section at grid location %s" % x)
                cs.cbm_gen_mesh(**kwargs)
        return None

    def blade_run_anbax(self, loads=None, **kwargs):
        """
        runs anbax for every section

        Parameters
        ----------
        loads : dict, optional
            dictionary of the following keys and values, (default=None)
            F : nparray([[grid, F1, F2, F3]])
            M : nparray([[grid, M1, M2, M3]])

        """

        ac = ANBAXConfig()
        lst = []
        for (x, cs) in self.sections:
            if loads:
                ac.recover_flag = 1
                load = interp_loads(loads, x)
                for k,v in load.items():
                    setattr(ac,k,v)

            cs.config.anbax_cfg = ac
            print("STATUS:\t Running ANBAX at grid location %s" % (x))
            cs.cbm_run_anbax(**kwargs)
            lst.append([x, cs.BeamProperties])
        # self.anba_beam_properties = np.asarray(lst)
        self.beam_properties = np.asarray(lst)
        return None      
        
        
     
    def blade_exp_beam_props(self, cosy='local', style='DYMORE', eta_offset=0, solver='vabs', filename = None):
        """
        Exports the beam_properties in the 
        
        Parameters
        ----------
        cosy : str, optional
            either 'global' for the global beam coordinate system or 
            'local' for a coordinate system that is always pointing with 
            the chord-line (in the twisted frame)
        
        style : str, optional
            select the style you want the beam_properties to be exported
            'DYMORE' will return an array of the following form:
            [[Massterms(6) (m00, mEta2, mEta3, m33, m23, m22)
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]]
            ...
            
        eta_offset : float, optional
            if the beam eta coordinates from start to end of the beam doesn't 
            coincide with the global coorinate system of the blade. The unit
            is in nondimensional r coordinates (x/Radius)
            
        solver : str, optional
            solver : if multiple or other solvers than vabs were applied, use 
            this option
        
        filename : str, optional
            if the user wants to write the output to a file. 
            
        Returns
        ----------
        arr : ndarray
            an array that reprensents the beam properties for the 
        """

        lst = []
        for cs in self.sections:
            # collect data for each section
            R = self.blade_ref_axis[-1, 1]
            # eta = -eta_offset/(1-eta_offset) + (1/(1-eta_offset))*cs[0]
            eta = (cs[0] * R) - (eta_offset * R)
            if style == "DYMORE":
                lst.append(cs[1].cbm_exp_dymore_beamprops(eta=eta, solver=solver))

            elif style == "BeamDyn":
                lst.append(cs[1].cbm_exp_BeamDyn_beamprops(eta=eta, solver=solver))

            elif style == "CAMRADII":
                pass

            elif style == "CPLambda":
                pass

        arr = np.asarray(lst)

        return arr

    def blade_plot_attributes(self):
        """
        plot the coordinates, chord, twist and pitch axis location of the blade

        Returns
        -------
        None.

        """
        fig, ax = plt.subplots(3, 2)
        fig.suptitle(self.name, fontsize=16)
        fig.subplots_adjust(wspace=0.25, hspace=0.25)

        ax[0][0].plot(self.blade_ref_axis[:, 0], self.blade_ref_axis[:, 1], "k.-")
        ax[0][0].set_ylabel("x-coordinate [m]")

        ax[1][0].plot(self.blade_ref_axis[:, 0], self.blade_ref_axis[:, 2], "k.-")
        ax[1][0].set_ylabel("y-coordinate [m]")

        ax[2][0].plot(self.blade_ref_axis[:, 0], self.blade_ref_axis[:, 3], "k.-")
        ax[2][0].set_ylabel("z-coordinate [m]")

        ax[0][1].plot(self.chord[:, 0], self.chord[:, 1], "k.-")
        ax[0][1].set_ylabel("chord [m]")

        ax[1][1].plot(self.twist[:, 0], self.twist[:, 1], "k.-")
        ax[1][1].set_ylabel("twist [rad]")

        ax[2][1].plot(self.pitch_axis[:, 0], self.pitch_axis[:, 1], "k.-")
        ax[2][1].set_ylabel("pitch axis location [1/chord]")

        #        ax3d = fig.add_subplot(326, projection='3d')
        #        for bm, af in zip(self.blade_matrix, self.airfoil):
        #            tmp_shape = af.coordinates[:,0].shape
        #            arr = af.coordinates*bm[4]
        #            ax3d.plot(np.ones(tmp_shape)*bm[1],arr[:,0],arr[:,1])
        plt.show()

    def blade_plot_beam_props(self, **kwargs):
        """
        plots the beam properties of the blade

        Parameters
        ----------
        **kwargs : TYPE
            keyword arguments can be passed down to the plot such as 
            sigma=None, ref=None, x_offset = 0, description = True

        Returns
        -------
        None.

        """
        plot_beam_properties(self.blade_exp_beam_props(), **kwargs)

    def blade_plot_sections(self, **kwargs):
        """
        plots the different sections of the blade
        """      
        for (x,cs) in self.sections:
            print('STATUS:\t Plotting section at grid location %s' % x)
            string = 'Blade: '+ str(self.name) + '; Section %.3f: ' % x
            cs.cbm_post_2dmesh(title=string, section = str(x), **kwargs)
        return None    
    
    def blade_post_3dtopo(self, flag_wf=True, flag_lft=False, flag_topo=False, flag_mesh=False, flag_wopwop=False):
        """
        generates the wireframe and the loft surface of the blade

        Returns
        ----------
        loft : OCC.TopoDS_surface
            the 3D surface of the blade
        
        wireframe : list
            list of every airfoil_wire scaled and rotated at every grid point
            
            
        ToDo
        ----------

        """
        (self.display, self.start_display, self.add_menu, self.add_function_to_menu) = display_config(cs_size=0.5, DeviationAngle=1e-4, DeviationCoefficient=1e-4)

        if flag_wf:
            wireframe = []

            # visualize blade and beam reference axis
            for s in self.blade_ref_axis_BSplineLst:
                self.display.DisplayShape(s, color="RED")

            for s in self.beam_ref_axis_BSplineLst:
                self.display.DisplayShape(s, color="GREEN")

            # airfoil wireframe
            for bm, afl in zip(self.blade_matrix, self.airfoils[:, 1]):
                (wire, te_pnt) = afl.trsf_to_blfr(bm[1:4], bm[6], bm[4], bm[5])
                wireframe.append(wire)
                self.display.DisplayShape(wire, color='BLACK')


        if flag_lft:
            # # step/iges file export
            # from jobs.RFeil.utls.import_export_step_files import STEPExporter
            # AP214_stepExporter = STEPExporter('loft_AP214.step', schema='AP214CD')  # init for writing step file; alternatively: schema='AP203'

            for i in range(len(wireframe)-1):
                # loft = make_loft(wireframe[i:i+2], ruled=True, tolerance=1e-2, continuity=1, check_compatibility=True)
                loft = make_loft(wireframe[i:i+2], ruled=True, tolerance=1e-6, continuity=1, check_compatibility=True)
                self.display.DisplayShape(loft, transparency=0.5, update=True)
                # if self.loft is not None:
                #     self.display.DisplayShape(self.loft, transparency=0.2, update=True, color="GREEN")
            #     AP214_stepExporter.add_shape(loft)  # add each lofted shape to the AP203_stepExporter component to generate full blade
            # AP214_stepExporter.write_file()  # write step file


        if flag_topo:
            for (x, cs) in self.sections:
                # display sections
                display_Ax2(self.display, cs.Ax2, length=0.2)
                display_cbm_SegmentLst(self.display, cs.SegmentLst, self.Ax2, cs.Ax2)

        if flag_wopwop:
            for bspl in self.wopwop_bsplinelst:
                for s in bspl:
                    self.display.DisplayShape(s, color="GREEN")

            for i, cs in enumerate(self.wopwop_pnts):
                for j, p1 in enumerate(cs):
                    v2 = self.wopwop_vecs[i][j]
                    v1 = gp_Vec(p1.XYZ())
                    v2.Normalize()
                    v2.Multiply(0.1)
                    v3 = v1.Added(v2)
                    p2 = gp_Pnt(v3.XYZ())

        self.display.View_Iso()
        self.display.FitAll()
        self.start_display()
        
    def blade_exp_beam_props(self, cosy='local', style='DYMORE', eta_offset=0, solver='anbax', filename = None):
        """
        Exports the beam_properties in the 
        
        Parameters
        ----------
        cosy : str, optional
            either 'global' for the global beam coordinate system or 
            'local' for a coordinate system that is always pointing with 
            the chord-line (in the twisted frame)
        
        style : str, optional
            select the style you want the beam_properties to be exported
            'DYMORE' will return an array of the following form:
            [[Massterms(6) (m00, mEta2, mEta3, m33, m23, m22)
            Stiffness(21) (k11, k12, k22, k13, k23, k33,... k16, k26, ...k66)
            Viscous Damping(1) mu, Curvilinear coordinate(1) eta]]
            ...
            
        eta_offset : float, optional
            if the beam eta coordinates from start to end of the beam doesn't 
            coincide with the global coorinate system of the blade. The unit
            is in nondimensional r coordinates (x/Radius)
            
        solver : str, optional
            solver : if multiple or other solvers than vabs were applied, use 
            this option
        
        filename : str, optional
            if the user wants to write the output to a file. 
            
        Returns
        ----------
        arr : ndarray
            an array that reprensents the beam properties for the 
        """

        lst = []
        for cs in self.sections:
            # collect data for each section
            R = self.blade_ref_axis[-1, 1]
            # eta = -eta_offset/(1-eta_offset) + (1/(1-eta_offset))*cs[0]
            eta = (cs[0] * R) - (eta_offset * R)
            if style == "DYMORE":
                lst.append(cs[1].cbm_exp_dymore_beamprops(eta=eta, solver=solver))

            elif style == "BeamDyn":
                lst.append(cs[1].cbm_exp_BeamDyn_beamprops(eta=eta, solver=solver))

            elif style == "CAMRADII":
                pass

            elif style == "CPLambda":
                pass

        arr = np.asarray(lst)

        return arr







# ====== M A I N ==============
if __name__ == "__main__":
    plt.close("all")

    #% ====== WindTurbine ==============
    with open("../jobs/PBortolotti/IEAonshoreWT.yaml", "r") as myfile:
        inputs = myfile.read()
    with open("../jobs/PBortolotti/IEAontology_schema.yaml", "r") as myfile:
        schema = myfile.read()
    # validate(yaml.load(inputs), yaml.load(schema))
    yml = yaml.load(inputs)

    airfoils = [Airfoil(af) for af in yml.get("airfoils")]
    materials = read_materials(yml.get("materials"))

    job = Blade(name="IEAonshoreWT")
    job.read_yaml(yml.get("components").get("blade"), airfoils, materials, wt_flag=True)