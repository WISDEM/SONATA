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

# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
              "attribute_str": attribute_str,
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
              "flag_recovery": flag_recovery, "c2_axis": c2_axis}


# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
radial_stations =  [0., 0.25, 0.5, 0.75, 1.]
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
flag_3d = False
flag_csv_export = False                         # export csv files with structural data
# Update flags dictionary
flags_dict['flag_csv_export'] = flag_csv_export
flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert

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

job.blade_run_viscoelastic()

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
