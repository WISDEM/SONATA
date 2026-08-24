import os
import time
from SONATA.classBlade import Blade


#------------------------------------------------------------------------------------
# Box beam analysis with SONATA
#------------------------------------------------------------------------------------
"""Box beam analysis with SONATA.

Problem Description:
====================

This script performs a structural analysis of a composite box beam using SONATA.
The analysis is based on the Smith-Chopra 1991 reference case and follows the
Popescu/Hodges 2000 paper specifications.

Problem Setup:
==============
- Rectangular box beam geometry with dimensions 24.2062 mm × 13.6398 mm (width × height)
- Beam length: 762 mm
- Composite layup: 6-ply antisymmetric laminate [30°/0°/30°/0°/30°/0°]
- Each ply thickness: 0.127 mm
- Material: AS4/3501-6 graphite/epoxy composite
- Analysis stations: root (r/R=0.0) and tip (r/R=1.0)

Analysis Objectives:
=====================
- Generate finite element mesh for the composite cross-section
- Compute sectional properties using b3_secfem (Asymptotic Numerical Beam Analysis eXtended)
- Evaluate structural response and material distribution
- Validate against reference solutions from literature

Reference:
=============
- Smith, E.C., Chopra, I. (1991). "Formulation and evaluation of an analytical model for
composite box beams." Journal of the American Helicopter Society, 36(3), 23-35.
- Popescu, B., Hodges, D.H. (2000). "On asymptotically correct Timoshenko-like anisotropic
beam theory." International Journal of Solids and Structures, 37(3), 535-558.

Units:
======
All dimensions in mm for better meshing convergence
"""

start_time = time.time()
print('Current working directory is:', os.getcwd())

#-----------------------------------------
# File paths
#-----------------------------------------
run_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep
# Note: for better meshing convergence, units specified in yaml are in 'mm' instead of 'm'
job_str = '0_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91.yaml'
job_name = 'box_beam_SmithChopra91'
filename_str = run_dir + job_str

#-----------------------------------------
# Analysis configuration
#-----------------------------------------
# General flags
flag_3d = False
flag_wt_ontology = False
flag_ref_axes_wt = False

# Plotting configuration
attribute_str = 'MatID'
flag_plotTheta11 = False        # plane orientation angle
flag_recovery = False
flag_plotDisplacement = False   # Needs recovery flag to be activated, shows displacements from loadings in cross sectional plots

# 3D visualization flags
flag_wf = True                  # plot wire-frame
flag_lft = True                 # plot lofted shape of blade surface
flag_topo = True                # plot mesh topology

# Shape configuration of corners
choose_cutoff = 2               # 0: step, 2: round corners

# Create configuration dictionary
flags_dict = {
    "flag_wt_ontology": flag_wt_ontology,
    "flag_ref_axes_wt": flag_ref_axes_wt,
    "attribute_str": attribute_str,
    "flag_plotDisplacement": flag_plotDisplacement,
    "flag_plotTheta11": flag_plotTheta11,
    "flag_wf": flag_wf,
    "flag_lft": flag_lft,
    "flag_topo": flag_topo,
    "flag_recovery": flag_recovery
}

#-----------------------------------------
# Execute SONATA analysis
#-----------------------------------------
# Define analysis stations
radial_stations = [0.0, 1.0]

print(f"Initializing blade analysis: {job_name}")
job = Blade(
    name=job_name,
    filename=filename_str,
    flags=flags_dict,
    stations=radial_stations
)

print("Generating sections and mesh...")
job.blade_gen_section(topo_flag=True, mesh_flag=True)

print("Running b3_secfem analysis...")
job.blade_run_b3_secfem()

#-----------------------------------------
# Generate plots and summary
#-----------------------------------------
print("Generating plots for box beam analysis...")
# Plot sections
job.blade_plot_sections(
    attribute=attribute_str,
    plotTheta11=flag_plotTheta11,
    plotDisplacement=flag_plotDisplacement
)

# Plot 3D topology
if flag_3d:
    job.blade_post_3dtopo(
        flag_wf=flags_dict['flag_wf'],
        flag_lft=flags_dict['flag_lft'],
        flag_topo=flags_dict['flag_topo']
    )

# Summary of analysis
elapsed_time = time.time() - start_time
print(f"\n{'='*60}")
print("ANALYSIS SUMMARY")
print(f"{'='*60}")
print(f"Job name: {job_name}")
print(f"Blade file: {filename_str}")
print(f"Number of radial stations: {len(radial_stations)}")
print(f"Station locations: {radial_stations}")
print("Analysis type: Box beam with antisymmetric layup")
print(f"Mesh generation: {'Enabled' if True else 'Disabled'}")
print(f"Topology analysis: {'Enabled' if True else 'Disabled'}")
print(f"3D visualization: {'Enabled' if flag_3d else 'Disabled'}")
print(f"Displacement plots: {'Enabled' if flag_plotDisplacement else 'Disabled'}")
print(f"Theta11 plots: {'Enabled' if flag_plotTheta11 else 'Disabled'}")
print(f"Recovery analysis: {'Enabled' if flag_recovery else 'Disabled'}")
print(f"Computational time: {elapsed_time:.2f} seconds")
print(f"{'='*60}")
print("Analysis completed successfully!")
