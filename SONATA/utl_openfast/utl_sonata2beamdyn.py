# -*- coding: utf-8 -*-
"""
Created on Thursday Oct 10 10:52:28 2019

@author: Roland Feil
"""
# ============================================= #
"""
Convert SONATA structural results to OpenFAST input files

Functions:
    write_beamdyn_axis  - writes '*_BeamDyn.dat' file
    write_beamdyn_prop  - writes '*_BeamDyn_Blade.dat' file

Inputs:
    wt_name         - name of currently investigated concept; i.e. wind turbine 
    byml            - yml.get('components').get('blade')
    radial_stations - radial stations along the blade span for inertia and stiffness characteristics
    beam_stiff      - matrix containing beam stiffness values
    beam_inertia    - matrix containing beam inertia values

Outputs: 
    - None -

"""
# ============================================= #

import os
import numpy as np
from scipy.interpolate import interp1d

from SONATA.cbm.cbm_utl import trsf_sixbysix

if __name__ == '__main__':
    os.chdir('../..')

# --- Convert from SONATA/VABS coordinates to BeamDyn coordinates ---#
def convert_structdef_SONATA_to_beamdyn(cs_pos, SONATA_beam_prop):
    """
    Convert structural characteristics from SONATA definition to BeamDyn definition
    
    Inputs:
        cs_pos              - array of radial stations along the span
        SONATA_beam_prop    - data struct containing the direct results from VABS in SONATA/VABS definition; equiv to job.beam_properties
    
    Outputs:
        BeamDyn_beam_prop   - converted data struct in BeamDyn definition
    """

    # initially used matrices and arrays
    beam_stiff_init = np.zeros([len(cs_pos), 6, 6])
    beam_inertia_init = np.zeros([len(cs_pos), 6, 6])
    beam_mass_center_init = np.zeros([len(cs_pos), 2])
    beam_neutral_axes_init = np.zeros([len(cs_pos), 2])
    beam_geometric_center_init = np.zeros([len(cs_pos), 2])
    beam_shear_center_init = np.zeros([len(cs_pos), 2])

    BeamDyn_beam_prop = {}
    BeamDyn_beam_prop['beam_section_mass'] = np.zeros([len(cs_pos), 1])
    BeamDyn_beam_prop['beam_stiff'] = np.zeros([len(cs_pos), 6, 6])
    BeamDyn_beam_prop['beam_inertia'] = np.zeros([len(cs_pos), 6, 6])
    BeamDyn_beam_prop['beam_mass_center'] = np.zeros([len(cs_pos), 2])
    BeamDyn_beam_prop['beam_neutral_axes'] = np.zeros([len(cs_pos), 2])
    BeamDyn_beam_prop['beam_geometric_center'] = np.zeros([len(cs_pos), 2])
    BeamDyn_beam_prop['beam_shear_center'] = np.zeros([len(cs_pos), 2])

    # --------------------------------------- #
    # retrieve & allocate VABS results
    for i in range(len(SONATA_beam_prop)):
        if SONATA_beam_prop[i, 1] is not None:
            BeamDyn_beam_prop['beam_section_mass'][i] = SONATA_beam_prop[i, 1].m00  # mass per unit span (absolute - no transform needed)
            
            # use init for the following in order to transform afterwards
            beam_mass_center_init[i] = np.array(SONATA_beam_prop[i, 1].Xm[:])  # Center of gravity (mass center)
            beam_neutral_axes_init[i] = np.array(SONATA_beam_prop[i, 1].Xt[:])  # Neutral axes (tension center)
            beam_geometric_center_init[i] = np.array(SONATA_beam_prop[i, 1].Xg[:])  # Geometric center
            beam_shear_center_init[i] = np.array(SONATA_beam_prop[i, 1].Xs[:])  # Generalized Shear Center of the Cross Section
            for j in range(6):
                beam_stiff_init[i, j, :] = np.array(SONATA_beam_prop[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
                beam_inertia_init[i, j, :] = np.array(SONATA_beam_prop[i, 1].MM[j, :])  # receive 6x6 mass matrix
                # beam_inertia_init[i, j, :] = np.array(SONATA_beam_prop[i, 1].MMatMC[j, :])  # receive 6x6 mass matrix at mass center
        else:  # If no solution from VABS available
            print('Radial station ' + str(cs_pos[i]) + ' did not run successfully in VABS. Check yaml input or change location!')
            # ToDo: instead of break -> skip or interpolate in between working radial stations

    # --------------------------------------- #
    #  rotate VABS results from SONATA/VABS def to BeamDyn def coordinate system
    
    # B = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # transformation matrix
    B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])  # NEW transformation matrix
    T = np.dot(np.identity(3), np.linalg.inv(B))

    for i in range(len(cs_pos)):
        BeamDyn_beam_prop['beam_stiff'][i, :, :] = trsf_sixbysix(beam_stiff_init[i, :, :], T)
        BeamDyn_beam_prop['beam_inertia'][i, :, :] = trsf_sixbysix(beam_inertia_init[i, :, :], T)
        BeamDyn_beam_prop['beam_mass_center'][i, :] = [-beam_mass_center_init[i, 0], beam_mass_center_init[i, 1]]
        BeamDyn_beam_prop['beam_neutral_axes'][i, :] = [-beam_neutral_axes_init[i, 0], beam_neutral_axes_init[i, 1]]
        BeamDyn_beam_prop['beam_geometric_center'][i, :] = [-beam_geometric_center_init[i, 0], beam_geometric_center_init[i, 1]]
        BeamDyn_beam_prop['beam_shear_center'][i, :] = [-beam_shear_center_init[i, 0], beam_shear_center_init[i, 1]]

    print('STATUS:\t Structural characteristics of VABS converted from SONATA/VABS to BeamDyn coordinate system definition!')

    return BeamDyn_beam_prop






# --- Write BeamDyn file with blade reference line locations ---#
def write_beamdyn_axis(folder, flags_dict, wt_name, ra, twist):

    n_pts = 50
    grid = np.linspace(0, 1, n_pts)

    f_interp = interp1d(ra[:,0], ra[:,3])
    kp_xr = f_interp(grid)
    f_interp = interp1d(ra[:,0], -ra[:,2])
    kp_yr = f_interp(grid)
    f_interp = interp1d(ra[:,0], ra[:,1])
    kp_zr = f_interp(grid)
    f_interp = interp1d(twist[:,0], np.rad2deg(twist[:,1]))
    twist_interp = f_interp(grid)

    data = np.vstack((kp_xr, kp_yr, kp_zr, twist_interp)).T

    # file = open(folder + '00_analysis/analysis/' + wt_name + '_BeamDyn.dat', 'w')
    file = open(os.path.join(folder , wt_name + '_BeamDyn.dat'), 'w')
    file.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
    file.write('%s blade\n' % (wt_name))
    file.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
    file.write('True          Echo            - Echo input data to "<RootName>.ech" (flag)\n')
    file.write('True          QuasiStaticInit - Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n')
    file.write(' 0            rhoinf          - Numerical damping parameter for generalized-alpha integrator\n')
    file.write(' 2            quadrature      - Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n')
    file.write('"DEFAULT"     refine          - Refinement factor for trapezoidal quadrature (-). DEFAULT = 1 [used only when quadrature=2]\n')
    file.write('"DEFAULT"     n_fact          - Factorization frequency (-). DEFAULT = 5\n')
    file.write('"DEFAULT"     DTBeam          - Time step size (s).\n')
    file.write('"DEFAULT"     load_retries    - Number of factored load retries before quitting the aimulation\n')
    file.write('"DEFAULT"     NRMax           - Max number of iterations in Newton-Ralphson algorithm (-). DEFAULT = 10\n')
    file.write('"DEFAULT"     stop_tol        - Tolerance for stopping criterion (-)\n')
    file.write('"DEFAULT"     tngt_stf_fd     - Flag to use finite differenced tangent stiffness matrix (-)\n')
    file.write('"DEFAULT"     tngt_stf_comp   - Flag to compare analytical finite differenced tangent stiffness matrix  (-)\n')
    file.write('"DEFAULT"     tngt_stf_pert   - perturbation size for finite differencing (-)\n')
    file.write('"DEFAULT"     tngt_stf_difftol- Maximum allowable relative difference between analytical and fd tangent stiffness (-)\n')
    file.write('True          RotStates       - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n')
    file.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
    file.write('          1   member_total    - Total number of members (-)\n')
    file.write('         %u   kp_total        - Total number of key points (-) [must be at least 3]\n' % (n_pts))
    file.write('     1     %u                 - Member number; Number of key points in this member\n' % (n_pts))
    file.write('\t\t kp_xr \t\t\t kp_yr \t\t\t kp_zr \t\t initial_twist\n')
    file.write('\t\t  (m)  \t\t\t  (m)  \t\t\t  (m)  \t\t   (deg)\n')
    if flags_dict['flag_write_BeamDyn_unit_convert'] == 'mm_to_m':  # convert units from mm (yaml input) to m
        for i in range(n_pts):
            file.write('\t %.5e \t %.5e \t %.5e \t %.5e \n' % (1.e-3*data[i, 0], 1.e-3*data[i, 1], 1.e-3*data[i, 2], 1.e-3*data[i, 3]))
        print('STATUS: converted from mm to m for export to BeamDyn.dat file')
    else:  # no unit or scaling conversion
        for i in range(n_pts):
            file.write('\t %.5e \t %.5e \t %.5e \t %.5e \n' % (data[i, 0], data[i, 1], data[i, 2], data[i, 3]))

    file.write('---------------------- MESH PARAMETER ------------------------------------------\n')
    file.write('          10   order_elem     - Order of interpolation (basis) function (-)\n')
    file.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
    file.write('"%s"    BldFile - Name of file containing properties for blade (quoted string)\n' % (wt_name + '_BeamDyn_Blade.dat'))
    file.write('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------\n')
    file.write('False         UsePitchAct - Whether a pitch actuator should be used (flag)\n')
    file.write('        200   PitchJ      - Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]\n')
    file.write('      2E+07   PitchK      - Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]\n')
    file.write('     500000   PitchC      - Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]\n')
    file.write('---------------------- OUTPUTS -------------------------------------------------\n')
    file.write('False          SumPrint       - Print summary data to "<RootName>.sum" (flag)\n')
    file.write('"ES10.3E2"    OutFmt         - Format used for text tabular output, excluding the time channel.\n')
    file.write('          1   NNodeOuts      - Number of nodes to output to file [0 - 9] (-)\n')
    file.write('          1,          2,          3,          4,          5,          6    OutNd          - Nodes whose values will be output  (-)\n')
    file.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')
    file.write('"RootFxr, RootFyr, RootFzr"\n')
    file.write('"RootMxr, RootMyr, RootMzr"\n')
    # file.write('"N1Fxl, N1Fyl, N1Fzl"\n')
    # file.write('"N1Mxl, N1Myl, N1Mzl"\n')
    file.write('"TipTDxr, TipTDyr, TipTDzr"\n')
    # file.write('"TipRDxr, TipRDyr, TipRDzr"\n')
    file.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')
    file.write('---------------------------------------------------------------------------------------\n')

    file.close()

    print('Finished writing BeamDyn File')

    return None

# --- Write BeamDyn_Blade file with blade properties ---#
def write_beamdyn_prop(folder, flags_dict, wt_name, radial_stations,
                       beam_stiff, beam_inertia, mu, format_name='BeamDyn'):
    n_pts = len(radial_stations)

    # file = open(folder + '00_analysis/analysis/' + wt_name + '_BeamDyn_Blade.dat', 'w')
    file = open(os.path.join(folder , wt_name + '_' + format_name 
                             + '_Blade.dat'), 'w')
    
    if not (format_name == 'BeamDyn'):
        file.write(' ------- NOT A BeamDyn Input File. Format is: {:}'
                   .format(format_name) + ' --------------------------\n')
        
    
    file.write(' ------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
    file.write(' Test Format 1\n')
    file.write(' ---------------------- BLADE PARAMETERS --------------------------------------\n')
    file.write('%u   station_total    - Number of blade input stations (-)\n' % (n_pts))
    file.write(' 1   damp_type        - Damping type: 0: no damping; 1: damped\n')
    file.write('  ---------------------- DAMPING COEFFICIENT------------------------------------\n')
    file.write('   mu1        mu2        mu3        mu4        mu5        mu6\n')
    file.write('   (-)        (-)        (-)        (-)        (-)        (-)\n')
    file.write('\t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\n' % (mu[0], mu[1], mu[2], mu[3], mu[4], mu[5])) 
    file.write(' ---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')

    if flags_dict['flag_write_BeamDyn_unit_convert'] == 'mm_to_m':  # convert units from mm (yaml input) to m
        beam_stiff = convert_stiff_matrix(beam_stiff, flags_dict)
        beam_inertia = convert_inertia_matrix(beam_inertia, flags_dict)
        print('STATUS: converted from mm to m for export to BeamDyn_Blade.dat file')

    for i in range(n_pts):
        file.write('\t %.6f \n' % (radial_stations[i]))
        # write stiffness matrices
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_stiff[i, j, 0], beam_stiff[i, j, 1], beam_stiff[i, j, 2], beam_stiff[i, j, 3], beam_stiff[i, j, 4],
            beam_stiff[i, j, 5]))
        file.write('\n')

        # write inertia properties
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_inertia[i, j, 0], beam_inertia[i, j, 1], beam_inertia[i, j, 2], beam_inertia[i, j, 3],
            beam_inertia[i, j, 4], beam_inertia[i, j, 5]))
        file.write('\n')
        # ToDO: check correct translation of stiffness and mass matrices from VABS and anbax !!!
    file.close()

    print('Finished writing BeamDyn_Blade File')

    return None

# --- Write BeamDyn_Blade file with blade properties ---#
def write_beamdyn_viscoelastic(folder, flags_dict, wt_name, radial_stations,
                               time_scales, beam_viscoelastic,
                               format_name='BeamDyn'):
    n_pts = len(radial_stations)

    output_name = os.path.join(folder , wt_name + '_' + format_name
                               + '_Blade_Viscoelastic.dat')
    
    # file = open(folder + '00_analysis/analysis/' + wt_name + '_BeamDyn_Blade.dat', 'w')
    file = open(output_name, 'w')
    
    if not (format_name == 'BeamDyn'):
        file.write(' ------- NOT A BeamDyn Input File. Format is: {:}'
                   .format(format_name) + ' --------------------------\n')
        
    file.write(' ------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
    file.write(' Test Format 1\n')
    file.write(' ---------------------- BLADE PARAMETERS --------------------------------------\n')
    file.write('{:d}   station_total        - Number of blade input stations (-)\n'
               .format(n_pts))
    file.write('{:d}   time_scales_total    - Number of blade input stations (-)\n'
               .format(len(time_scales)))
    file.write('  ---------------------- Time Scales (s) ------------------------------------\n')
    file.write(''.join(str(t)+'   ' for t in time_scales)[:-3] + '\n')
    file.write(' ---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')

    assert not (flags_dict['flag_write_BeamDyn_unit_convert'] == 'mm_to_m'), \
        'Viscoelastic beam property output does not support unit conversion.'

    for i in range(n_pts):
        file.write('\t %.6f \n' % (radial_stations[i]))
        # loop over time scales
        for k in range(len(time_scales)):
            
            curr_stiff = beam_viscoelastic[i, k, :, :]
        
            # write stiffness matrices
            for j in range(6):
                file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
                curr_stiff[j, 0], curr_stiff[j, 1], curr_stiff[j, 2],
                curr_stiff[j, 3], curr_stiff[j, 4], curr_stiff[j, 5]))
            file.write('\n')

        # ToDO: check correct translation of stiffness and mass matrices from VABS and anbax !!!
    file.close()

    print('STATUS:\t Finished writing BeamDyn_Blade_Viscoelastic file.')
    
    return None


def convert_stiff_matrix(beam_stiff, flags_dict):
    if flags_dict['flag_write_BeamDyn_unit_convert'] == 'mm_to_m':
        for j in range(6):  # row
            for k in range(6):  # column
                if j <= 2 and k <= 2:  # axial stiffness, EA
                    beam_stiff[:, j,k] = beam_stiff[:, j,k]  # N (left top part of matrix)
                elif j <= 2 and k > 2:  # extension & trans-shear to EI $ GJ coupling (right top part of matrix)
                    beam_stiff[:, j,k] = beam_stiff[:, j,k] * 1e-3  # Nmm to Nm
                elif j > 2 and k <= 2:  # extension & trans-shear to EI $ GJ coupling (left bottom part of matrix)
                    beam_stiff[:, j, k] = beam_stiff[:, j, k] * 1e-3  # Nmm to Nm
                elif j > 2 and k > 2 :  # EI and GJ coupling terms
                    beam_stiff[:, j,k] = beam_stiff[:, j,k] * 1e-6  # Nmm2 to Nm2

    return beam_stiff

def convert_inertia_matrix(beam_inertia, flags_dict):
    if flags_dict['flag_write_BeamDyn_unit_convert'] == 'mm_to_m':
        for j in range(6):  # row
            for k in range(6):  # column
                if (j == 0 and k == 0) or (j == 1 and k == 1) or (j == 2 and k == 2):  # mass per unit length
                    beam_inertia[:, j,k] = beam_inertia[:, j,k] * 1e3  # kg/mm -> kg/m
                elif (j == 3 and k == 3) or (j == 4 and k == 4) or (j == 5 and k == 5):  # mass moments of inertia
                    beam_inertia[:, j, k] = beam_inertia[:, j, k] * 1e-6  # kg mm2 -> kg m2
                elif (j == 3 and k == 4) or (j == 4 and k == 3):  # mass moments of inertia
                    beam_inertia[:, j, k] = beam_inertia[:, j, k] * 1e-6  # kg mm2 -> kg m2

    return beam_inertia

# ==============
# Main
# ==============

if __name__ == '__main__':
    import yaml
    from SONATA.classAirfoil import Airfoil
    #from SONATA.classBlade import Blade

    from SONATA.classAirfoil import Airfoil
    from SONATA.classMaterial import read_materials

    # provide primary path; used for providing yaml input file as well as output directory
    folder = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
    filename = (folder + 'IEAonshoreWT_BAR_005a.yaml')
    with open(filename, 'r') as myfile:
        inputs = myfile.read()
        yml = yaml.load(inputs, Loader=yaml.FullLoader)

    airfoils = [Airfoil(af) for af in yml.get('airfoils')]
    materials = read_materials(yml.get('materials'))
    # Blade.read_IEA37(yml.get('components').get('blade'), airfoils, **kwargs)
    wt_name = yml.get('name')

    # radial_stations = [0.0, 0.5, 1.0] # must include 0 and 1!
    radial_stations = np.linspace(0.0, 1.0, 11)

    beam_stiff = np.zeros([len('radial_stations'), 6, 6])
    beam_inertia = np.zeros([len('radial_stations'), 6, 6])

    write_beamdyn_axis(folder, wt_name, yml.get('components').get('blade'))
    write_beamdyn_prop(folder, wt_name, radial_stations, beam_stiff, beam_inertia)


