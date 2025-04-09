"""
Regression pytest for Beam Example 0
"""

import os
import numpy as np
from SONATA.classBlade import Blade

def test_6x6_beam0():
    
    # ===== Provide Path Directory & Yaml Filename ===== #
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = '0_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91.yaml'  
    # note: for better meshing convergence, units specified in yaml are in 'mm' instead of 'm'
    job_name = 'box_beam_SmithChopra91'
    
    filename_str = run_dir + job_str
    
    
    # ===== Define flags ===== #
    flag_wt_ontology        = False
    flag_ref_axes_wt        = False
    attribute_str           = 'MatID'
    flag_plotTheta11        = False     # plane orientation angle
    flag_recovery           = False
    flag_plotDisplacement   = False     # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots
    flag_wf                 = True      # plot wire-frame
    flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
    flag_topo               = True      # plot mesh topology

    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology,
                  "flag_ref_axes_wt": flag_ref_axes_wt,
                  "attribute_str": attribute_str,
                  "flag_plotDisplacement": flag_plotDisplacement,
                  "flag_plotTheta11": flag_plotTheta11,
                  "flag_wf": flag_wf,
                  "flag_lft": flag_lft,
                  "flag_topo": flag_topo,
                  "flag_recovery": flag_recovery}
    
    radial_stations = [0.0, 1.0]
    
    # Create job structure
    job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations)
    job.blade_gen_section(topo_flag=True, mesh_flag = True)
    
    job.blade_run_anbax()
    
    ###########################################################################
    ######## Checks on if answer looks consistent with previous runs ##########
    ###########################################################################
    
    # Reference 6x6 timoshenko stiffness matrix
    ref_TS = np.array([[ 5.50900295e+06,  3.59669546e+01, -4.99861142e+00,
                            -6.08011530e+06,  2.70016342e+03, -4.99226065e+02],
                       [ 3.59669547e+01,  4.30331012e+05, -1.63748365e+01,
                            4.76228312e+01,  3.01674716e+06, -6.92650427e+02],
                       [-4.99861142e+00, -1.63748365e+01,  1.83506118e+05,
                            2.22593606e+02,  9.99057189e+01,  3.13632190e+06],
                       [-6.08011530e+06,  4.76228312e+01,  2.22593606e+02,
                            5.24670335e+07, -4.13539061e+03,  7.97500454e+02],
                       [ 2.70016342e+03,  3.01674716e+06,  9.99057189e+01,
                            -4.13539061e+03,  1.72356571e+08, -3.47095240e+03],
                       [-4.99226065e+02, -6.92650427e+02,  3.13632190e+06,
                            7.97500454e+02, -3.47095240e+03,  4.29057741e+08]])
    
    # Reference mass matrix
    ref_MM = np.array([[55.35444883, 0., 0., 0., 0.0006865, 0.00167523],
                       [0., 55.35444883, 0., -0.0006865, 0., 0.],
                       [0., 0., 55.35444883, -0.00167523, 0., 0.],
                       [0., -0.0006865, -0.00167523, 6096.39902105, 0., 0.],
                       [0.0006865, 0., 0., 0., 1757.09202585, -0.01111826],
                       [0.00167523, 0., 0., 0., -0.01111826, 4339.3069952]])
    
    for i in range(job.beam_properties.shape[0]):
        
        # 6x6 timoshenko stiffness matrix
        assert np.allclose(job.beam_properties[i, 1].TS, ref_TS), \
            "Stiffness matrix does not match."
        
        # 6x6 mass matrix
        assert np.allclose(job.beam_properties[i, 1].MM, ref_MM), \
            "Mass matrix does not match."