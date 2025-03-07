import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval, strain_energy_eval


# Path to yaml file
run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
job_str = '8_box_beam.yaml'
job_name = 'Box_Beam'
filename_str = run_dir + job_str

# ===== Define flags ===== #
flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)

# --- plotting flags ---
# Define mesh resolution, i.e. the number of points along the profile that is used for out-to-inboard meshing of a 2D blade cross section
mesh_resolution = 400
# For plots within blade_plot_sections
attribute_str           = 'stress.sigma11'  # default: 'MatID' (theta_3 - fiber orientation angle)
                                            # others:  'theta_3' - fiber orientation angle
                                            #          'stress.sigma11' (use sigma_ij to address specific component)
                                            #          'stressM.sigma11'
                                            #          'strain.epsilon11' (use epsilon_ij to address specific component)
                                            #          'strainM.epsilon11'

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
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
              "attribute_str": attribute_str,
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
              "flag_recovery": flag_recovery, "c2_axis": c2_axis}


# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
radial_stations =  [0.0, 1.0]
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
flag_constant_loads = True

if flag_constant_loads:
    # forces, N (F1: axial force
    #            F2: x-direction shear force
    #            F3: y-direction shear force)
    # moments, Nm (M1: torsional moment,
    #              M2: bending moment about x, (axis parallel to chord)
    #              M3: bending moment around y)
    Loads_dict = {"Forces":[1.0e4,0.0,0.0],"Moments":[0.0,0.0,0.0]}
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
delta = np.array([0.03, 0.03, 0.06787]) # logarithmic decrement, natural log of the ratio of the amplitudes of any two successive peaks. 3% flap and edge, 6% torsion
zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.) # damping ratio,  dimensionless measure describing how oscillations in a system decay after a disturbance
omega = np.array([0.508286, 0.694685, 4.084712])*2*np.pi # Frequency (rad/s), flap/edge/torsion
mu1 = 2*zeta[0]/omega[0]
mu2 = 2*zeta[1]/omega[1]
mu3 = 2*zeta[2]/omega[2]
mu = np.array([mu1, mu2, mu3, mu2, mu1, mu3])
beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir, job_str, mu)

# ===== PLOTS ===== #
# job.blade_plot_attributes()
# job.blade_plot_beam_props()

# saves figures in folder_str/figures if savepath is provided:
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=flag_plotDisplacement, savepath=run_dir)
if flag_3d:
    job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'], flag_lft=flags_dict['flag_lft'], flag_topo=flags_dict['flag_topo'])

# ===== Strain Energy Recovery ===== #

if flag_recovery:
    # This example only has one material, but in other cases can use 'MatID'
    # input to filter by a specific material.
    total_energy, directional_energy, strain_all, energy_all \
        = strain_energy_eval(job, MatID=None)

    print('\n\nFraction of Strain Energy in Components:')
    component_order = ['11', '22', '33', '23', '13', '12']

    for ind,component in enumerate(component_order):
        print('Strain Energy {} : {:.4f}'.format(component,
                                         directional_energy[ind]/total_energy))

# ===== Analytical Calculations for Stress Verification ===== #

if flag_constant_loads:
    reference_moments = Loads_dict['Moments']
    reference_forces = Loads_dict['Forces']

else:
    # Take forces/moments at the root.
    reference_moments = Loads_dict['Moments'][0, 1:]
    reference_forces = Loads_dict['Forces'][0, 1:]

print('\n\nAnalytical Stresses for Rectangular Input:')

thickness = 0.1

height_outer = 1 #m
width_outer = 2 #m


height_inner = height_outer - 2*thickness # m
width_inner = width_outer - 2*thickness # m

Ix = (1/12)*((height_outer**3)*width_outer - (height_inner**3)*width_inner)
Iy = (1/12)*((width_outer**3)*height_outer - (width_inner**3)*height_inner)

area = height_outer*width_outer - height_inner*width_inner

print('\nAverage sigma11 for just Fx is: {:.2f}'.format(
    reference_forces[0]/area))


# ===== Checking Stress-Strain Relations at Elements for Orientation ===== #

# 1. Pick an element in the center of each side
# 2. Recover the stress and strain at that element
# 3. Transform that stress-strain to the expected fiber orientation
# 4. Check if the elastic modulus along that axis is appropriate.

section = job.sections[0]
(x,cs) = section
cells = cs.mesh

desired_coords_list = [np.array([0.0, 0.45]),
                       np.array([0.0, -0.45]),
                       np.array([0.95, 0.0])]
cell_list = len(desired_coords_list) * [None]

for coord_ind, desired_coords in enumerate(desired_coords_list):
    dist = np.zeros(len(cells))
    theta_1 = np.zeros(len(cells))
    theta_3 = np.zeros(len(cells))
    
    for ind,c in enumerate(cells):
        dist[ind] = np.linalg.norm(c.calc_center() - desired_coords)
        theta_1[ind] = c.theta_1[0]
        theta_3[ind] = c.theta_3
    
    cell_ind = np.argmin(dist)
    cell_list[coord_ind] = cells[cell_ind]

c = cell_list[0]
c.stress.tensor
c.strain.tensor

print("Theta_1: top={:.3f}, bottom={:.3f}".format(cell_list[0].theta_1[0],
                                                  cell_list[1].theta_1[0]))

print("Theta_3: top={:.3f}, bottom={:.3f}".format(cell_list[0].theta_3,
                                                  cell_list[1].theta_3))

print("Stress Tensor (global coords) - top: ")
print(cell_list[0].stress.tensor)

print("Stress Tensor (global coords) - bottom: ")
print(cell_list[1].stress.tensor)

print("Stress Tensor (global coords) - leading edge: ")
print(cell_list[2].stress.tensor)

print("")
print("Stress Tensor (local coords) - top: ")
print(cell_list[0].stressM.tensor)

print("Stress Tensor (local coords) - bottom: ")
print(cell_list[1].stressM.tensor)

print("Stress Tensor (local coords) - leading edge: ")
print(cell_list[2].stressM.tensor)