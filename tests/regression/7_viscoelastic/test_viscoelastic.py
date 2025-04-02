import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval, strain_energy_eval


import matplotlib
import matplotlib.pyplot as plt

import sys
import pytest

sys.path.append(os.path.join(os.path.dirname( os.path.realpath(__file__)),
                             '..'))

import utils

def test_one_viscoelastic():
    """
    One term test of viscoelastic SONATA output. Uses orthotropic material, but
    inputs are equivalent to isotropic.

    Returns
    -------
    None.

    """
    
    original_backend = matplotlib.get_backend()
    matplotlib.use('Agg')
    
    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = '7_one_term.yaml'
    job_name = 'Box_Beam'
    filename_str = run_dir + job_str
    
    # ===== Define flags ===== #
    flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
    flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
    
    flag_viscoelastic = True # flag to run viscoelastic 6x6 matrices calculation.
    
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
    flag_recovery           = False     # Set to True to Plot stresses/strains
    flag_plotDisplacement   = True     # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots
    
    # 3D plots (blade_post_3dtopo)
    flag_wf                 = True      # plot wire-frame
    flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
    flag_topo               = True      # plot mesh topology
    c2_axis                 = False
    flag_DeamDyn_def_transform = True               # transform from SONATA to BeamDyn coordinate system
    flag_write_BeamDyn = True                       # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
    flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'     # applied only when exported to BeamDyn files
    flag_OpenTurbine_transform = True
    flag_write_OpenTurbine = True
    
    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
                  "attribute_str": attribute_str,
                  "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
                  "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
                  "flag_recovery": flag_recovery, "c2_axis": c2_axis}
    
    
    # ===== User defined radial stations ===== #
    # Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
    # radial_stations =  [0., 1.]
    radial_stations = [0.0]
    # radial_stations = [.7]
    
    # ===== Execute SONATA Blade Component Object ===== #
    # name          - job name of current task
    # filename      - string combining the defined folder directory and the job name
    # flags         - communicates flag dictionary (defined above)
    # stations      - input of radial stations for cross sectional analysis
    # stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
    job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations)  # initialize job with respective yaml input file
    
    # ===== Build & mesh segments ===== #
    job.blade_gen_section(topo_flag=True, mesh_flag = True)
    
    
    # ===== Recovery Analysis + BeamDyn Outputs ===== #
    
    # # Define flags
    flag_csv_export = False                         # export csv files with structural data
    # Update flags dictionary
    flags_dict['flag_csv_export'] = flag_csv_export
    flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
    flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
    flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
    flags_dict['viscoelastic'] = flag_viscoelastic
    flags_dict['flag_OpenTurbine_transform'] = flag_OpenTurbine_transform
    flags_dict['flag_write_OpenTurbine'] = flag_write_OpenTurbine
    
    Loads_dict = {"Forces":[0.0,0.0,0.0],"Moments":[0.0,1.0e3,0.0]}
    
    # Set damping for BeamDyn input file
    mu = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
    # job.blade_run_viscoelastic()
    
    beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir,
                     job_str, mu)

    plt.close('all')
    matplotlib.use(original_backend)
    
    elastic_file = 'Box_Beam1_BeamDyn_Blade.dat'
    visco_file = 'Box_Beam1_BeamDyn_Blade_Viscoelastic.dat'
    
    elastic_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            elastic_file)

    visco_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             visco_file)
    
    
    mass, stiff = utils.load_bd_blade(elastic_path)
    
    prony_stiff = utils.load_bd_visco(visco_path)
    
    # Because input has elastic properties at time scales that sum to the 
    # total used in warping field, 6x6 matrices should add to match.
    tot_prony = prony_stiff[0][0] + prony_stiff[0][1]
        
    assert np.allclose(tot_prony[:3], stiff[0][:3], atol=0.00001), \
        "First 3 rows of viscoelastic do not add to match the elastic."
    
    assert np.allclose(tot_prony[3:], stiff[0][3:],
                       atol=0.0001*tot_prony[3:].max()), \
        "Last 3 rows of viscoelastic do not add to match the elastic."
    
    # Check ratios of terms match that of inputs
    mask = np.abs(tot_prony) > 1.0
    
    assert np.allclose((prony_stiff[0][0] / tot_prony)[mask], 1./21.), \
        "Ratio of viscoelastic terms does not match ratio of inputs."
    
def test_two_viscoelastic():
    """
    Two term test of viscoelastic SONATA output. Uses isotropic material.

    Verifies against reference file that gave correct damping in GEBT for
    properties used as inputs here.

    Returns
    -------
    None.

    """
    
    original_backend = matplotlib.get_backend()
    matplotlib.use('Agg')
    
    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = '7_two_term.yaml'
    job_name = 'Box_Beam'
    filename_str = run_dir + job_str
    
    # ===== Define flags ===== #
    flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
    flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
    
    flag_viscoelastic = True # flag to run viscoelastic 6x6 matrices calculation.
    
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
    flag_recovery           = False     # Set to True to Plot stresses/strains
    flag_plotDisplacement   = True     # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots
    
    # 3D plots (blade_post_3dtopo)
    flag_wf                 = True      # plot wire-frame
    flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
    flag_topo               = True      # plot mesh topology
    c2_axis                 = False
    flag_DeamDyn_def_transform = True               # transform from SONATA to BeamDyn coordinate system
    flag_write_BeamDyn = True                       # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
    flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'     # applied only when exported to BeamDyn files
    flag_OpenTurbine_transform = True
    flag_write_OpenTurbine = True
    
    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
                  "attribute_str": attribute_str,
                  "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
                  "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
                  "flag_recovery": flag_recovery, "c2_axis": c2_axis}
    
    
    # ===== User defined radial stations ===== #
    # Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
    # radial_stations =  [0., 1.]
    radial_stations = [0.0]
    # radial_stations = [.7]
    
    # ===== Execute SONATA Blade Component Object ===== #
    # name          - job name of current task
    # filename      - string combining the defined folder directory and the job name
    # flags         - communicates flag dictionary (defined above)
    # stations      - input of radial stations for cross sectional analysis
    # stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
    job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations)  # initialize job with respective yaml input file
    
    # ===== Build & mesh segments ===== #
    job.blade_gen_section(topo_flag=True, mesh_flag = True)
    
    
    # ===== Recovery Analysis + BeamDyn Outputs ===== #
    
    # # Define flags
    flag_csv_export = False                         # export csv files with structural data
    # Update flags dictionary
    flags_dict['flag_csv_export'] = flag_csv_export
    flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
    flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
    flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
    flags_dict['viscoelastic'] = flag_viscoelastic
    flags_dict['flag_OpenTurbine_transform'] = flag_OpenTurbine_transform
    flags_dict['flag_write_OpenTurbine'] = flag_write_OpenTurbine
    
    Loads_dict = {"Forces":[0.0,0.0,0.0],"Moments":[0.0,1.0e3,0.0]}
    
    # Set damping for BeamDyn input file
    mu = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
    # job.blade_run_viscoelastic()
    
    beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir,
                     job_str, mu)

    plt.close('all')
    matplotlib.use(original_backend)
    
    ref_file = 'ref_two_term_bd_blade_visc.dat'
    visco_file = 'Box_Beam2_BeamDyn_Blade_Viscoelastic.dat'
    
    ref_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            ref_file)

    visco_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             visco_file)
    
    ref_prony_stiff = utils.load_bd_visco(ref_path)
    
    prony_stiff = utils.load_bd_visco(visco_path)
    
    for i in range(3):
        
        print("term {:} error: {:}".format(i,
            np.abs(ref_prony_stiff[0][i] - prony_stiff[0][i]).max()))
        
        assert np.allclose(prony_stiff[0][i], ref_prony_stiff[0][i],
                           atol=1e-3), \
            "Viscoelastic 6x6 doesn't match reference."
            

if __name__ == "__main__":
    pytest.main(["-s", "test_viscoelastic.py"])
