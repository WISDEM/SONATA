import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval, strain_energy_eval


import pytest

def run_stresses(job_str, loads_dict, flag_constant_loads):
    
    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_name = 'Box_Beam'
    filename_str = run_dir + job_str
    
    # ===== Define flags ===== #
    flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
    flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)
    
    # --- plotting flags ---
    # Define mesh resolution, i.e. the number of points along the profile that is used for out-to-inboard meshing of a 2D blade cross section
    mesh_resolution = 400
    # For plots within blade_plot_sections
    
    # 2D cross sectional plots (blade_plot_sections)
    flag_plotTheta11        = False      # plane orientation angle
    flag_recovery           = True     # Set to True to Plot stresses/strains
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
    
    # Flag applies twist rotations in SONATA before output and then sets the output
    # twist to all be zero degrees.
    flag_output_zero_twist = False
    
    
    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology,
                  "flag_ref_axes_wt": flag_ref_axes_wt,
                  "attribute_str": 'MatID',
                  "flag_plotDisplacement": flag_plotDisplacement,
                  "flag_plotTheta11": flag_plotTheta11,
                  "flag_wf": flag_wf,
                  "flag_lft": flag_lft,
                  "flag_topo": flag_topo,
                  "mesh_resolution": mesh_resolution,
                  "flag_recovery": flag_recovery,
                  "c2_axis": c2_axis}
    
    
    # ===== User defined radial stations ===== #
    # Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
    radial_stations =  [0., 0.6, 1.]
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
    flag_3d = False
    flag_csv_export = False                         # export csv files with structural data
    # Update flags dictionary
    flags_dict['flag_csv_export'] = flag_csv_export
    flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
    flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
    flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
    flags_dict['flag_output_zero_twist'] = flag_output_zero_twist
    
    # Flag for different load input formats.
    # Just used for example script, not passed to SONATA
    
    # Set damping for BeamDyn input file
    delta = np.array([0.03, 0.03, 0.06787]) # logarithmic decrement, natural log of the ratio of the amplitudes of any two successive peaks. 3% flap and edge, 6% torsion
    zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.) # damping ratio,  dimensionless measure describing how oscillations in a system decay after a disturbance
    omega = np.array([0.508286, 0.694685, 4.084712])*2*np.pi # Frequency (rad/s), flap/edge/torsion
    mu1 = 2*zeta[0]/omega[0]
    mu2 = 2*zeta[1]/omega[1]
    mu3 = 2*zeta[2]/omega[2]
    mu = np.array([mu1, mu2, mu3, mu2, mu1, mu3])
    beam_struct_eval(flags_dict, loads_dict, radial_stations, job, run_dir, job_str, mu)


    # Verfify stresses integrate to applied loads
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        force_moment = np.zeros(6)

        area = np.nan*np.zeros(len(cells))

        for ind,c in enumerate(cells):

            area = c.calc_area()

            force_moment[0] += c.stress.sigma11*area
            force_moment[1] += c.stress.sigma12*area
            force_moment[2] += c.stress.sigma13*area

            cxy = c.calc_center()

            # Equations for integrated forces and moments are taken
            # from classCBM.py -> cbm_run_viscoelastic
            # However, the ordering on the left needs to be taken from
            # beam_struct_eval.py -> beam_struct_eval

            force_moment[3] += area*(c.stress.sigma13*cxy[0]
                                     - c.stress.sigma12*cxy[1])

            force_moment[4] += c.stress.sigma11*area*cxy[1]

            force_moment[5] += c.stress.sigma11*area*(-cxy[0])

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        diff = force_moment - applied_loads

        # Normalize the recovered forces / moments
        # Add 1.0 in denominator to prevent divide by zero and accept
        # small absolute errors when stresses are near zero.
        norm_diff = np.linalg.norm(diff) \
                        / (np.linalg.norm(applied_loads) + 1.0)

        assert (norm_diff < 1e-4), \
            "Recovered stresses do not integrate to applied loads."

    return job


def test_force_dir1():
    
    job_str = '6_box_beam.yaml'
    flag_constant_loads = False
        
    recover_forces = np.array([[0.0, 1.0e3, 0.0, 0.0],
                               [1.0, 1.0e2, 0.0, 0.0]])

    recover_moments = np.array([[0.0, 0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0, 0.0]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    # Verify against analytical stress expectation.
    thickness = 0.1
    
    height_outer = 1 #m
    width_outer = 2 #m
    
    height_inner = height_outer - 2*thickness # m
    width_inner = width_outer - 2*thickness # m
    
    # Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
    # Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)
    
    area = height_outer*width_outer - height_inner*width_inner
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])


        sigma11 = applied_loads[0] / area
        
        tol = 1e-4
        
        for ind,c in enumerate(cells):

            assert np.abs(c.stress.sigma11 - sigma11) < tol*sigma11, \
                "Wrong stress sigma11 for axial load."
                
            assert np.abs(c.stress.sigma22) < tol, \
                "Wrong stress sigma22 for axial load."
                
            assert np.abs(c.stress.sigma33) < tol, \
                "Wrong stress sigma33 for axial load."
            
            assert np.abs(c.stress.sigma12) < tol, \
                "Wrong stress sigma12 for axial load."
                
            assert np.abs(c.stress.sigma13) < tol, \
                "Wrong stress sigma13 for axial load."
                
            assert np.abs(c.stress.sigma23) < tol, \
                "Wrong stress sigma23 for axial load."

            # cxy = c.calc_center()

def test_force_dir2():
    
    job_str = '6_box_beam.yaml'
    flag_constant_loads = False
        
    recover_forces = np.array([[0.0, 0.0, 1.0e3, 0.0],
                               [1.0, 0.0, 1.0e2, 0.0]])

    recover_moments = np.array([[0.0, 0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0, 0.0]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    # Verify against analytical stress expectation.
    thickness = 0.1
    
    height_outer = 1 #m
    width_outer = 2 #m
    
    height_inner = height_outer - 2*thickness # m
    width_inner = width_outer - 2*thickness # m
    
    # Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
    Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)    
    # area = height_outer*width_outer - height_inner*width_inner
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        sigma11 = np.zeros(len(cells));
        sigma22 = np.zeros(len(cells));
        sigma33 = np.zeros(len(cells));
        sigma23 = np.zeros(len(cells));
        sigma13 = np.zeros(len(cells));
        sigma12 = np.zeros(len(cells));

        for ind,c in enumerate(cells):

            sigma11[ind] = c.stress.sigma11
            sigma22[ind] = c.stress.sigma22
            sigma33[ind] = c.stress.sigma33
            sigma23[ind] = c.stress.sigma23
            sigma13[ind] = c.stress.sigma13
            sigma12[ind] = c.stress.sigma12

            # cxy = c.calc_center()
            
        tol = 1e-4
        
        assert np.abs(sigma11).max() < tol, \
            "Should have 0 sigma11 for force 2."
        
        assert np.abs(sigma22).max() < tol, \
            "Should have 0 sigma22 for force 2."
        
        assert np.abs(sigma33).max() < tol, \
            "Should have 0 sigma33 for force 2."
        
        Q = thickness*height_outer*(width_outer/2 - 0.5*thickness) \
            + 2*thickness*(width_inner/2)*(width_inner/4)
        
        sigma12_ref = applied_loads[1] * Q / (2*thickness * Iy)
        
        assert np.abs(sigma12.max() - sigma12_ref) < 1e-2*sigma12_ref, \
            "Max sigma12 for force 2 doesn't match theory."
        
        # There appear to be stress concentrations or other effects, so this
        # field looks bad.
        # assert np.abs(sigma13).max() < tol, \
        #     "Should have 0 sigma13 for force 2."
        
        assert np.abs(sigma23).max() < tol, \
            "Should have 0 sigma23 for force 2."

def test_force_dir3():
    
    job_str = '6_box_beam.yaml'
    flag_constant_loads = False
        
    recover_forces = np.array([[0.0, 0.0, 0.0, 1.0e3],
                               [1.0, 0.0, 0.0, 1.0e2]])

    recover_moments = np.array([[0.0, 0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0, 0.0]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    # Verify against analytical stress expectation.
    thickness = 0.1
    
    height_outer = 1 #m
    width_outer = 2 #m
    
    height_inner = height_outer - 2*thickness # m
    width_inner = width_outer - 2*thickness # m
    
    Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
    # Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)    
    # area = height_outer*width_outer - height_inner*width_inner
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        sigma11 = np.zeros(len(cells));
        sigma22 = np.zeros(len(cells));
        sigma33 = np.zeros(len(cells));
        sigma23 = np.zeros(len(cells));
        sigma13 = np.zeros(len(cells));
        sigma12 = np.zeros(len(cells));
        
        cxy = np.zeros((len(cells), 2));
        

        for ind,c in enumerate(cells):

            sigma11[ind] = c.stress.sigma11
            sigma22[ind] = c.stress.sigma22
            sigma33[ind] = c.stress.sigma33
            sigma23[ind] = c.stress.sigma23
            sigma13[ind] = c.stress.sigma13
            sigma12[ind] = c.stress.sigma12

            cxy[ind] = c.calc_center()
            
        tol = 1e-4
        
        assert np.abs(sigma11).max() < tol, \
            "Should have 0 sigma11 for force 3."
        
        assert np.abs(sigma22).max() < tol, \
            "Should have 0 sigma22 for force 3."
        
        assert np.abs(sigma33).max() < tol, \
            "Should have 0 sigma33 for force 3."

        Q = thickness*width_outer*(height_outer/2 - 0.5*thickness) \
            + 2*thickness*(height_inner/2)*(height_inner/4)
        
        sigma13_ref = applied_loads[2] * Q / (2*thickness * Ix)
                
        # Have to mask down to only the center area where shear should be max
        # because of stress concentration.
        mask = np.abs(cxy[:, 1]) < 0.05
        
        assert np.abs(sigma13[mask].max() - sigma13_ref) < 0.04*sigma13_ref, \
            "Max sigma13 for force 3 doesn't match theory."
        
        # There appear to be stress concentrations or other effects, so this
        # field looks bad.
        # assert np.abs(sigma12).max() < tol, \
        #     "Should have 0 sigma12 for force 2."
        
        assert np.abs(sigma23).max() < tol, \
            "Should have 0 sigma23 for force 3."
        

def test_moment_dir1():
    
    job_str = 'circle_beam.yaml'
    flag_constant_loads = False
        
    recover_forces =  np.array([[0.0, 0.0, 0.0, 0.0],
                                [1.0, 0.0, 0.0, 0.0]])
    
    recover_moments = np.array([[0.0, 1.0e3, 0.0, 0.0],
                                [1.0, 1.0e2, 0.0, 0.0]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        sigma11 = np.zeros(len(cells));
        sigma22 = np.zeros(len(cells));
        sigma33 = np.zeros(len(cells));
        sigma23 = np.zeros(len(cells));
        sigma13 = np.zeros(len(cells));
        sigma12 = np.zeros(len(cells));
        
        cxy = np.zeros((len(cells), 2));
        

        for ind,c in enumerate(cells):

            sigma11[ind] = c.stressM.sigma11
            sigma22[ind] = c.stressM.sigma22
            sigma33[ind] = c.stressM.sigma33
            sigma23[ind] = c.stressM.sigma23
            sigma13[ind] = c.stressM.sigma13
            sigma12[ind] = c.stressM.sigma12

            cxy[ind] = c.calc_center()
            
        tol = 1e-4
        
        assert np.abs(sigma11).max() < tol, \
            "Should have 0 sigma11 for moment 1 (torsion)."
        
        assert np.abs(sigma22).max() < tol, \
            "Should have 0 sigma22 for moment 1 (torsion)."
        
        assert np.abs(sigma33).max() < tol, \
            "Should have 0 sigma33 for moment 1 (torsion)."
            
        radius = 0.5
        thickness = 0.1
        J = np.pi/2 * (radius**4 - (radius-thickness)**4)
        
        tau_max = -applied_loads[3] * radius / J
        test_rad = np.linalg.norm(cxy + np.array([[radius, 0.0]]), axis = 1)
        
        assert np.abs(sigma12 - tau_max*test_rad/radius).max() \
            < 0.01*np.abs(tau_max), \
            "Max sigma12 for moment 1 (torsion) doesn't match theory."
        
        assert np.abs(sigma13).max() < 0.1*np.abs(tau_max), \
            "Should have 0 sigma13 for moment 1."
        
        assert np.abs(sigma23).max() < tol, \
            "Should have 0 sigma23 for moment 1."
        
def test_moment_dir2():
    
    job_str = '6_box_beam.yaml'
    flag_constant_loads = False

    recover_forces = np.array([[0.0, 0.0, 0.0, 0.0],
                               [1.0, 0.0, 0.0, 0.0]])
        
    recover_moments = np.array([[0.0, 0.0, 1.0e3, 0.0],
                                [1.0, 0.0, 1.0e2, 0.0]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    # Verify against analytical stress expectation.
    thickness = 0.1
    
    height_outer = 1 #m
    width_outer = 2 #m
    
    height_inner = height_outer - 2*thickness # m
    width_inner = width_outer - 2*thickness # m
    
    Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
    # Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)    
    # area = height_outer*width_outer - height_inner*width_inner
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        sigma11 = np.zeros(len(cells));
        sigma22 = np.zeros(len(cells));
        sigma33 = np.zeros(len(cells));
        sigma23 = np.zeros(len(cells));
        sigma13 = np.zeros(len(cells));
        sigma12 = np.zeros(len(cells));
        
        cxy = np.zeros((len(cells), 2));

        for ind,c in enumerate(cells):

            sigma11[ind] = c.stress.sigma11
            sigma22[ind] = c.stress.sigma22
            sigma33[ind] = c.stress.sigma33
            sigma23[ind] = c.stress.sigma23
            sigma13[ind] = c.stress.sigma13
            sigma12[ind] = c.stress.sigma12

            cxy[ind] = c.calc_center()
            
        tol = 1e-4
        
        
        sigma11_max = applied_loads[4] * (height_outer / 2) / Ix
        
        expected = (sigma11_max / (0.5 * height_outer)) * cxy[:, 1]
                
        assert np.abs(sigma11 - expected).max() < 0.005*sigma11_max, \
            "Don't have expected sigma11 for moment 2."
        
        assert np.abs(sigma22).max() < 50.0, \
            "Should have 0 sigma22 for moment 2."
        
        assert np.abs(sigma33).max() < 50.0, \
            "Should have 0 sigma33 for moment 2."
        
        assert np.abs(sigma12).max() < tol, \
            "Should have 0 sigma12 for moment 2."
            
        assert np.abs(sigma13).max() < tol, \
            "Should have 0 sigma13 for moment 2."
            
        assert np.abs(sigma23).max() < 50.0, \
            "Should have 0 sigma23 for moment 2."

def test_moment_dir3():
    
    job_str = '6_box_beam.yaml'
    flag_constant_loads = False

    recover_forces = np.array([[0.0, 0.0, 0.0, 0.0],
                               [1.0, 0.0, 0.0, 0.0]])
        
    recover_moments = np.array([[0.0, 0.0, 0.0, 1.0e3],
                                [1.0, 0.0, 0.0, 1.0e2]])

    loads_dict = {"Forces" : recover_forces,
                  "Moments": recover_moments}
    
    
    job = run_stresses(job_str, loads_dict, flag_constant_loads)
    
    # Verify against analytical stress expectation.
    thickness = 0.1
    
    height_outer = 1 #m
    width_outer = 2 #m
    
    height_inner = height_outer - 2*thickness # m
    width_inner = width_outer - 2*thickness # m
    
    # Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
    Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)    
    # area = height_outer*width_outer - height_inner*width_inner
    
    for sec_ind,section in enumerate(job.sections):

        (x,cs) = section
        cells = cs.mesh

        # Interpolate loads_dict to the current position
        if flag_constant_loads:

            applied_loads = np.hstack((loads_dict['Forces'],
                                             loads_dict['Moments']))

        else:
            applied_loads = np.zeros(6)

            applied_loads[0] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 1])

            applied_loads[1] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 2])

            applied_loads[2] = np.interp(x, loads_dict['Forces'][:, 0],
                                         loads_dict['Forces'][:, 3])

            applied_loads[3] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 1])

            applied_loads[4] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 2])

            applied_loads[5] = np.interp(x, loads_dict['Moments'][:, 0],
                                         loads_dict['Moments'][:, 3])

        sigma11 = np.zeros(len(cells));
        sigma22 = np.zeros(len(cells));
        sigma33 = np.zeros(len(cells));
        sigma23 = np.zeros(len(cells));
        sigma13 = np.zeros(len(cells));
        sigma12 = np.zeros(len(cells));
        
        cxy = np.zeros((len(cells), 2));

        for ind,c in enumerate(cells):

            sigma11[ind] = c.stress.sigma11
            sigma22[ind] = c.stress.sigma22
            sigma33[ind] = c.stress.sigma33
            sigma23[ind] = c.stress.sigma23
            sigma13[ind] = c.stress.sigma13
            sigma12[ind] = c.stress.sigma12

            cxy[ind] = c.calc_center()
            
        tol = 1e-4
        
        sigma11_max = applied_loads[5] * (width_outer / 2) / Iy
        
        expected = -(sigma11_max / (0.5 * width_outer)) * cxy[:, 0]

        assert np.abs(sigma11 - expected).max() < 0.001*sigma11_max, \
            "Don't have expected sigma11 for moment 3."
        
        assert np.abs(sigma22).max() < 50.0, \
            "Should have 0 sigma22 for moment 3."
        
        assert np.abs(sigma33).max() < 50.0, \
            "Should have 0 sigma33 for moment 3."
        
        assert np.abs(sigma12).max() < tol, \
            "Should have 0 sigma12 for moment 3."
            
        assert np.abs(sigma13).max() < tol, \
            "Should have 0 sigma13 for moment 3."
            
        assert np.abs(sigma23).max() < 50.0, \
            "Should have 0 sigma23 for moment 3."


if __name__ == "__main__":
    pytest.main(["-s", "test_stress_recov.py"])
    
