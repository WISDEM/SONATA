import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval

import matplotlib
import matplotlib.pyplot as plt

import sys
import pytest

sys.path.append(os.path.join(os.path.dirname( os.path.realpath(__file__)),
                             '..'))

import utils

def test_6x6_iea22mw():
    
    original_backend = matplotlib.get_backend()
    matplotlib.use('Agg')
    
    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = 'iea_22_280_rwt.yaml'
    job_name = 'IEA22'
    filename_str = run_dir + job_str
    
    # ===== Define flags ===== #
    flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
    flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
    
    # --- plotting flags ---
    # Define mesh resolution, i.e. the number of points along the profile that is used for out-to-inboard meshing of a 2D blade cross section
    mesh_resolution = 400
    # For plots within blade_plot_sections
    attribute_str           = 'MatID'  # default: 'MatID' (theta_3 - fiber orientation angle)
                                                # others:  'theta_3' - fiber orientation angle
                                                #          'stress.sigma11' (use sigma_ij to address specific component)
                                                #          'stressM.sigma11'
                                                #          'strain.epsilon11' (use epsilon_ij to address specific component)
                                                #          'strainM.epsilon11'
    
    # 2D cross sectional plots (blade_plot_sections)
    flag_plotTheta11        = False      # plane orientation angle
    flag_recovery           = False
    flag_plotDisplacement   = True     # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots
    # 3D plots (blade_post_3dtopo)
    flag_wf                 = True      # plot wire-frame
    flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
    flag_topo               = True      # plot mesh topology
    c2_axis                 = False
    flag_DeamDyn_def_transform = True               # transform from SONATA to BeamDyn coordinate system
    flag_write_BeamDyn = True                       # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
    flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'     # applied only when exported to BeamDyn files
    
    # Shape of corners
    choose_cutoff = 2    # 0 step, 2 round
    
    
    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology,
                  "flag_ref_axes_wt": flag_ref_axes_wt,
                  "attribute_str": attribute_str,
                  "flag_plotDisplacement": flag_plotDisplacement,
                  "flag_plotTheta11": flag_plotTheta11,
                  "flag_wf": flag_wf,
                  "flag_lft": flag_lft,
                  "flag_topo": flag_topo,
                  "mesh_resolution": mesh_resolution,
                  "flag_recovery": flag_recovery,
                  "c2_axis": c2_axis}
    
    
    # ===== User defined radial statiosns ===== #
    # Define the radial stations for cross sectional analysis (only used for
    # flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
    radial_stations = [0.000000000000000e+00, 6.896551724137931e-02,
                       1.034482758620690e-01, 1.379310344827586e-01,
                       1.724137931034483e-01, 2.068965517241379e-01, 
                       2.413793103448276e-01, 2.758620689655172e-01,
                       3.103448275862069e-01, 3.448275862068966e-01,
                       3.793103448275862e-01, 4.137931034482759e-01,
                       4.482758620689655e-01, 
                       4.827586206896551e-01, 5.172413793103449e-01,
                       5.517241379310345e-01, 5.862068965517241e-01,
                       6.206896551724138e-01, 6.551724137931034e-01,
                       6.896551724137931e-01, 
                       7.241379310344828e-01, 7.586206896551724e-01,
                       7.931034482758621e-01, 8.275862068965517e-01,
                       8.620689655172413e-01, 8.965517241379310e-01,
                       9.310344827586207e-01, 
                       9.655172413793103e-01, 1.000000000000000e+00]
    
    # ===== Execute SONATA Blade Component Object ===== #
    # name          - job name of current task
    # filename      - string combining the defined folder directory and the job
    #                 name
    # flags         - communicates flag dictionary (defined above)
    # stations      - input of radial stations for cross sectional analysis
    job = Blade(name=job_name, filename=filename_str, flags=flags_dict,
                cutoff_style = choose_cutoff, stations=radial_stations)
    
    # ===== Build & mesh segments ===== #
    job.blade_gen_section(topo_flag=True, mesh_flag=True)
    
    # Define flags
    flag_csv_export = False

    # Update flags dictionary
    flags_dict['flag_csv_export'] = flag_csv_export
    flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
    flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
    flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
    Loads_dict = {"Forces":[1.,1.,1.],"Moments":[1.,1.,1.]}
    
    # Set damping for BeamDyn input file
    delta = np.array([0.03, 0.03, 0.06787])
    zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.)
    omega = np.array([0.508286, 0.694685, 4.084712])*2*np.pi
    
    mu1 = 2*zeta[0]/omega[0]
    mu2 = 2*zeta[1]/omega[1]
    mu3 = 2*zeta[2]/omega[2]
    mu = np.array([mu1, mu2, mu3, mu2, mu1, mu3])
    beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir,
                     job_str, mu)
    
    plt.close('all')
    matplotlib.use(original_backend)
    
    reference_file = 'ref_iea22mw_bd_blade.dat'
    test_file = 'IEA_22MW_BeamDyn_Blade.dat'
    
    ref_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            reference_file)
    
    test_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             test_file)
    
    utils.compare_bd_blade(ref_path, test_path, tolerance=1e-9)
        

if __name__ == "__main__":
    pytest.main(["-s", "test_iea22mw.py"])
