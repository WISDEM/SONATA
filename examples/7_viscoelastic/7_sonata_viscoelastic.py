import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval, strain_energy_eval


# Path to yaml file
run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
job_str = '7_hollow_rect.yaml'
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
radial_stations = np.linspace(0, 1, 21).tolist()
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
flags_dict['viscoelastic'] = flag_viscoelastic
flags_dict['flag_OpenTurbine_transform'] = flag_OpenTurbine_transform
flags_dict['flag_write_OpenTurbine'] = flag_write_OpenTurbine

# Flag for different load input formats.
# Just used for example script, not passed to SONATA
flag_constant_loads = False

if flag_constant_loads:
    # forces, N (F1: axial force
    #            F2: x-direction shear force
    #            F3: y-direction shear force)
    # moments, Nm (M1: torsional moment,
    #              M2: bending moment about x, (axis parallel to chord)
    #              M3: bending moment around y)
    Loads_dict = {"Forces":[0.0,0.0,0.0],"Moments":[0.0,1.0e3,0.0]}
else:

    # Forces and moments have a first column of the station (normalized length)
    # The next three columns are force/moment values at the given stations.
    # See above for the description of what the columns are.
    # Linear interpolation is used between stations.
    # Set forces or moments at the 0.0 station to have analytical stress/strain
    # output to compare to.

    recover_Forces = np.array([[0.0, 0.0, 0.0, 0.0],
                               [1.0, 0.0, 0.0, 0.0]])

    recover_Moments = np.array([[0.0, 1.0e3, 0.0, 0.0],
                                [1.0, 0.0, 0.0, 0.0]])

    Loads_dict = {"Forces" : recover_Forces,
                  "Moments": recover_Moments}

# Set damping for BeamDyn input file
mu = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# job.blade_run_viscoelastic()

beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir,
                 job_str, mu)

# ===== PLOTS ===== #
# job.blade_plot_attributes()
# job.blade_plot_beam_props()

# saves figures in folder_str/figures if savepath is provided:
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11,
                        plotDisplacement=flag_plotDisplacement,
                        savepath=run_dir)
if flag_3d:
    job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'],
                          flag_lft=flags_dict['flag_lft'],
                          flag_topo=flags_dict['flag_topo'])

# ===== Expected Damping Behavior ===== #

assert len(job.materials) == 1, 'Following calculations assume single material.'


assert job.materials[1].viscoelastic['time_scales_v'].shape[0] == 2,\
    'Following calculations assume single viscoelastic decaying time scale.'

# Frequency to evaluate viscoelastic properties at. 
# If checking a modal damping factor, this should be that mode's natural freq.
# Note this needs to be the actual response frequency and that the viscoelastic
# material adds nonlinearity and thus may not be simple to evaluate
freq = 1.0534000 # Hz
omega  = 2 * np.pi * freq

if job.materials[1].orth == 0:
    
    E_i = job.materials[1].viscoelastic['E_v'][0]
    tau_i = job.materials[1].viscoelastic['time_scales_v'][0]

    E_inf = job.materials[1].viscoelastic['E_v'][1]
    
    storage_mod = E_inf + (omega**2 * tau_i**2 * E_i)/(omega**2 * tau_i**2 + 1)
    loss_mod = (omega * tau_i * E_i)/(omega**2 * tau_i**2 + 1)
    tan_delta = loss_mod / storage_mod
    zeta = tan_delta / (2.0)
    
    
    print('\nStorage Modulus at {:.3f} Hz: {:.3e}'.format(freq, storage_mod))
    print('For linear eigenanalysis, this should match elastic value of {:.3e}'
          .format(job.materials[1].E))
    
    print('zeta [fraction critical damping] : {:.4e}'.format(zeta))
    
elif job.materials[1].orth == 1:
    
    mat_prop_dirs = ['E_1_v', 'E_2_v', 'E_3_v', 'G_12_v', 'G_13_v', 'G_23_v']
    
    ref_values = job.materials[1].E.tolist() + job.materials[1].G.tolist()
    
    tau_i = job.materials[1].viscoelastic['time_scales_v'][0]
        
    for ind, prop_key in enumerate(mat_prop_dirs):
        
        E_i = job.materials[1].viscoelastic[prop_key][0]
        E_inf = job.materials[1].viscoelastic[prop_key][1]
        
        storage_mod = E_inf + (omega**2 * tau_i**2 * E_i)/(omega**2 * tau_i**2 + 1)
        loss_mod = (omega * tau_i * E_i)/(omega**2 * tau_i**2 + 1)
        tan_delta = loss_mod / storage_mod
        zeta = tan_delta / (2.0)
        
        
        print('\nFor property: {:s}'.format(prop_key))
        print('Storage Modulus at {:.3f} Hz: {:.3e}'.format(freq, storage_mod))
        print('For linear eigenanalysis, this should match elastic value of {:.3e}'
              .format(ref_values[ind]))
        
        print('zeta [fraction critical damping] : {:.4e}'.format(zeta))

    print('\nIf all zeta values match, that should be modal zeta for a mode at' 
          + ' that frequency.')

print('\nIf (all) storage modulus match the elastic reference values, then'
      + ' mode/modal frequency is expected to be the same between'
      + ' viscoelastic and elastic models at the frequency checked here.')