# -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os

import time
from SONATA.classBlade import Blade


# ==============
# Main
# ==============

start_time = time.time()

print('Current working directory is:', os.getcwd())


# ===== Provide Path Directory & Yaml Filename ===== #
run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
job_str = '0_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91.yaml'  # note: for better meshing convergence, units specified in yaml are in 'mm' instead of 'm'
job_name = 'box_beam_SmithChopra91'

filename_str = run_dir + job_str


# ===== Define flags ===== #
flag_3d                 = False
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
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
              "attribute_str": attribute_str,
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo,
              "flag_recovery": flag_recovery}

radial_stations = [0.0, 1.0]

# Create job structure
job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations)
job.blade_gen_section(topo_flag=True, mesh_flag = True)

job.blade_run_anbax()

# ===== PLOTS ===== #
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=flag_plotDisplacement) #, savepath=folder_str)
if flag_3d:
    job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'], flag_lft=flags_dict['flag_lft'], flag_topo=flags_dict['flag_topo'])

print("--- Computational time: %s seconds ---" % (time.time() - start_time))

# EOF
