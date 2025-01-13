#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Nov 20 13:33:33 2019

@author: Roland Feil
"""


from builtins import len, range
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

from SONATA.cbm.cbm_utl import trsf_sixbysix
from SONATA.utl_openfast.utl_sonata2beamdyn import convert_structdef_SONATA_to_beamdyn, \
    write_beamdyn_axis, write_beamdyn_prop, write_beamdyn_viscoelastic

# from SONATA.utl.analytical_rectangle.utls_analytical_rectangle import utls_analytical_rectangle



def beam_struct_eval(flags_dict, loads_dict, cs_pos, job, folder_str, job_str, mu):

    """
    Analyse, transform, evaluate and plot structural results from VABS and/or ANBAX

    Functions:
    beam_struct_eval                    - parent function (retrieve & transform data, structure of data evaluation)
    plot_beam_props_6by6                - function to plot the 6x6 stiffness and mass matrices
    plot_beam_mass_distribution         - function to plot the beam mass per unit length distribution
    plot_vabs_anbax                     - function to plot the 6x6 stiffness and mass matrices from both VABS and ANBAX for code-to-code verification
    vabs_export_beam_struct_properties  - csv export of structural beam properties
    anbax_export_beam_struct_properties - csv export of structural beam properties

    Inputs:
    flags_dict          - dictionary containing relevant flags
    loads_dict          - dictionary containing the applied loads for recovery analysis
    cs_pos              - radial station of blade cross sections
    job                 - contains the whole blade data (yaml file content, wires, mesh, etc.)
    folder_str          - name of operating folder
    job_str             - name of job that is currently under investigation 

    Outputs: 
    - None -


    Dictionary Examples:
    flag_recovery           = True
    flags_dict = {"flag_recovery": flag_recovery}

    Forces =  np.array([0, 0, 0])  # forces, N (F1: axial force; F2,F3: sectional transverse shear forces)
    Moments =  np.array([0, 5000000, 0])  # moments, Nm
    loads_dict = {"Forces": Forces, "Moments": Moments}
    
    Notes
    -----
    
    Flag options for `flags_dict` include `'viscoelastic' : True` to evaluate
    viscoelastic 6x6 matrices and save them in an additional file. This flag
    is optional and is assumed `False` by default.

    """
    
    if 'viscoelastic' not in flags_dict.keys():
        flags_dict['viscoelastic'] = False

    # --- ANBAX --- #
    # --------------------------------------- #
    # clear var before initializing anbax
    job.beam_properties = None


    # --------------------------------------- #
    # --- ANBAX --- #
    if flags_dict['flag_recovery'] == True and not flags_dict['viscoelastic']:

        if np.asarray(loads_dict['Forces']).shape == (3,):
            # Assume the format is provided for just a uniform load at all
            # sections

            # forces, N (F1: shear force in x-direction;
            # F2: shear force in y -direction;
            # F3: axial force)
            Forces = [float(loads_dict["Forces"][1]),
                      float(loads_dict["Forces"][2]),
                      float(loads_dict["Forces"][0])]

            # moments, Nm (M1: bending moment around x;
            #              M2: bending moment around y;
            #              M3: torsional moment)
            Moments = [float(loads_dict["Moments"][1]),
                       float(loads_dict["Moments"][2]),
                       float(loads_dict["Moments"][0])]

            loads = {# sonata coord system input converted to anbax coordinates
                "F":    np.array([[0] + Forces,
                                  [1] + Forces]),
                "M":    np.array([[0] + Moments,
                                  [1] + Moments]),
                }

        else:
            # Input has first column of station (normalized). Other columns
            # are same as above, but have indices starting from 1.
            loads = {
                "F" : loads_dict["Forces"][:, [0, 2, 3, 1]],
                "M" : loads_dict["Moments"][:, [0, 2, 3, 1]]
                    }

        job.blade_run_anbax(loads)  # run anbax
    elif not flags_dict['viscoelastic']:
        job.blade_run_anbax()  # run anbax
    else:
        # flags_dict['viscoelastic'] == True
        
        if flags_dict['flag_recovery']:
            print('Recovery of stress/strain not supported with viscoelastic'
                  + ' material simulations.')
        
        job.blade_run_viscoelastic()

    # init used matrices and arrays
    anbax_beam_stiff_init = np.zeros([len(cs_pos), 6, 6])
    anbax_beam_inertia_init = np.zeros([len(cs_pos), 6, 6])
    anbax_beam_stiff = np.zeros([len(cs_pos), 6, 6])
    anbax_beam_inertia = np.zeros([len(cs_pos), 6, 6])
    anbax_beam_section_mass = np.zeros([len(cs_pos), 1])
    
    if flags_dict['viscoelastic']:
        anbax_beam_viscoelastic \
            = np.zeros((len(cs_pos), len(job.beam_properties[0][1].tau), 6, 6))

    # --------------------------------------- #
    # retrieve & allocate ANBAX results
    for i in range(len(job.beam_properties)):
        anbax_beam_section_mass[i] = job.beam_properties[i, 1].m00  # receive mass per unit span
        for j in range(6):
            anbax_beam_stiff_init[i, j, :] = np.array(job.beam_properties[i, 1].TS[j, :])  # receive 6x6 timoshenko stiffness matrix
            anbax_beam_inertia_init[i, j, :] = np.array(job.beam_properties[i, 1].MM[j, :])  # receive 6x6 mass matrix

        if flags_dict['viscoelastic']:
            for k in range(len(job.beam_properties[0][1].tau)):
                anbax_beam_viscoelastic[i, k, :, :] = job.beam_properties[i, 1].TSv[k]
    
    # --------------------------------------- #
    #  rotate anbax results from SONATA/VABS def to BeamDyn def coordinate system (for flag_DeamDyn_def_transform = True)
    if flags_dict['flag_DeamDyn_def_transform']:
        print('STATUS:\t Transform to BeamDyn coordinates')
        # B = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # transformation matrix
        B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])  # transformation matrix
        T = np.dot(np.identity(3), np.linalg.inv(B))
        for n_sec in range(len(cs_pos)):
            anbax_beam_stiff[n_sec, :, :] = trsf_sixbysix(anbax_beam_stiff_init[n_sec, :, :], T)
            anbax_beam_inertia[n_sec, :, :] = trsf_sixbysix(anbax_beam_inertia_init[n_sec, :, :], T)
            
            if flags_dict['viscoelastic']:
                for k in range(len(job.beam_properties[0][1].tau)):
                    anbax_beam_viscoelastic[i, k, :, :] = trsf_sixbysix(
                                        anbax_beam_viscoelastic[i, k, :, :], T)
                    
            
        str_ext = '_BeamDyn_def'
        coordsys = 'BeamDyn'

        print('STATUS:\t Structural characteristics of ANBAX converted from SONATA/VABS to BeamDyn coordinate system definition!')
    else:
        anbax_beam_stiff = anbax_beam_stiff_init
        anbax_beam_inertia = anbax_beam_inertia_init
        str_ext = ''
        coordsys = 'VABS/SONATA'


    # --------------------------------------- #
    # Export beam structural properties to csv file
    if flags_dict['flag_csv_export']:
        print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
        anbax_export_beam_struct_properties(folder_str, job_str, cs_pos, coordsys=coordsys, solver='anbax', beam_stiff=anbax_beam_stiff,
                                        beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)



    # --------------------------------------- #
    # Export beam structural properties to csv file
    # if flags_dict['flag_csv_export']:
    #     print('STATUS:\t Export csv files with structural blade characeristics from ANBAX to: ' + folder_str + 'csv_export/')
    #     export_beam_struct_properties(folder_str, job_str, cs_pos, solver='anbax', beam_stiff=anbax_beam_stiff, beam_inertia=anbax_beam_inertia, beam_mass_per_length=anbax_beam_section_mass)

    # ToDo: also export BeamDyn files for results from anbax as soon as the verification is completed
    # --------------------------------------- #
    # write BeamDyn input files
    np.savetxt('anbax_BAR00.txt', np.array([cs_pos, anbax_beam_stiff[:, 3, 3], anbax_beam_stiff[:, 4, 4], anbax_beam_stiff[:, 5, 5], anbax_beam_stiff[:, 2, 2], anbax_beam_inertia[:, 0, 0]]).T)
    if flags_dict['flag_write_BeamDyn'] & flags_dict['flag_DeamDyn_def_transform']:
        print('STATUS:\t Write BeamDyn input files')
        refine = int(30/len(cs_pos))  # initiate node refinement parameter
        write_beamdyn_axis(folder_str, flags_dict, job.yml.get('name'), job.blade_ref_axis, job.twist)
        write_beamdyn_prop(folder_str, flags_dict, job.yml.get('name'), cs_pos, anbax_beam_stiff, anbax_beam_inertia, mu)

        if flags_dict['viscoelastic']:
            
            print('STATUS:\t Writing viscoelastic BeamDyn input file.')
            write_beamdyn_viscoelastic(folder_str, flags_dict,
                                       job.yml.get('name'), cs_pos,
                                       job.beam_properties[0][1].tau,
                                       anbax_beam_viscoelastic)



# ============================================= #
def plot_beam_props_6by6(cs_pos, data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    # fig = plt.figure(tight_layout=True, figsize=(9, 6), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title)
    for i in range(len(data[0, :, 0])):
        for j in range(len(data[0, 0, :])):
            if j >= i:
                ax = fig.add_subplot(len(data[0, :, 0]), len(data[0, 0, :]), k)
                ax.plot(cs_pos, data[:, i, j], '-k')
                plt.ylim(1.1 * min(-1, min(data[:, i, j])), 1.1 * max(1, max(data[:, i, j])))
                ax.set_xlabel('r/R')
                if fig_title == 'Mass matrix':
                    ax.set_title('$m_{%i %i}$' % ((i + 1), (j + 1)))
                elif fig_title == 'Stiffness matrix':
                    ax.set_title('$k_{%i %i}$' % ((i + 1), (j + 1)))
                ax.grid(True)
            k = k + 1
    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None


def plot_beam_axes(cs_pos, vabs_beam_mass_center, vabs_beam_neutral_axes, vabs_beam_geometric_center,
               vabs_beam_shear_center, save_path):

    fig = plt.figure(tight_layout=True, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.axes(projection='3d')
    ax.plot3D(cs_pos, vabs_beam_mass_center[:,0], vabs_beam_mass_center[:,1], label='Mass center')
    ax.plot3D(cs_pos, vabs_beam_neutral_axes[:,0], vabs_beam_neutral_axes[:,1], label='Neutral axes')
    ax.plot3D(cs_pos, vabs_beam_geometric_center[:,0], vabs_beam_geometric_center[:,1], label='Geometric center')
    ax.plot3D(cs_pos, vabs_beam_shear_center[:,0], vabs_beam_shear_center[:,1], label='Shear Center')
    ax.set_xlabel('r/R')
    ax.set_ylabel('chordwise location, m')
    ax.set_zlabel('thickness location, m')
    ax.legend()

    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None

def plot_beam_mass_distribution(cs_pos, vabs_beam_mass_distribution, save_path):

    fig = plt.figure(tight_layout=True, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(cs_pos, vabs_beam_mass_distribution[:,0])
    plt.xlabel('r/R')
    plt.ylabel('Mass per unit length, N/m')
    # plt.ylim([0, 1000])

    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None




# ============================================= #
def plot_vabs_anbax(cs_pos, vabs_data, anbax_data, fig_title, save_path):
    # plots 6x6 matrix
    k = 1
    fig = plt.figure(tight_layout=True, figsize=(14, 10), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title)
    for i in range(len(vabs_data[0, :, 0])):
        for j in range(len(vabs_data[0, 0, :])):
            if j >= i:
                ax = fig.add_subplot(len(vabs_data[0, :, 0]), len(vabs_data[0, 0, :]), k)
                ax.plot(cs_pos, vabs_data[:, i, j], '--k')
                ax.plot(cs_pos, anbax_data[:, i, j], ':r')
                plt.ylim(1.1 * min(-1, min(vabs_data[:, i, j]), min(anbax_data[:, i, j])), 1.1 * max(1, max(vabs_data[:, i, j]), max(anbax_data[:, i, j])))
                ax.set_xlabel('r/R')
                if fig_title == 'Mass matrix':
                    ax.set_title('$m_{%i %i}$' % ((i + 1), (j + 1)))
                elif fig_title == 'Stiffness matrix':
                    ax.set_title('$k_{%i %i}$' % ((i + 1), (j + 1)))
                ax.grid(True)
            k = k + 1
    plt.show()
    fig.savefig(''.join(save_path), dpi=300)

    return None



# ============================================= #
def vabs_export_beam_struct_properties(folder_str, job_str, radial_stations, coordsys, solver, beam_stiff, beam_inertia, beam_mass_per_length,
                                  beam_mass_center, beam_neutral_axes, beam_geometric_center, beam_shear_center):

    if solver=='vabs':
        export_name_general = 'vabs_beam_properties_general.csv'
        export_name_stiff = 'vabs_beam_properties_stiff_matrices.csv'
        export_name_mass = 'vabs_beam_properties_mass_matrices.csv'
    else:
        print('Define correct solver name (vabs or anbax) when calling export_beam_struct_properties')

    # -------------------------------------------------- #
    # Export mass per unit length for the defined radial stations
    if os.path.isdir(folder_str + 'csv_export/') == False:
        os.mkdir(folder_str + 'csv_export/')
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_general]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if coordsys == 'BeamDyn':
            beam_prop_writer.writerow(['Coordinate system:', 'Beamdyn coordinates'])
        elif coordsys == 'VABS/SONATA':
            beam_prop_writer.writerow(['Coordinate system:', 'VABS/SONATA coordinates'])
        else:
            beam_prop_writer.writerow(['Coordinate system:', 'to be verified'])

        beam_prop_writer.writerow(['section in r/R', 'Mass per unit length [kg/m]',
                                      'Mass center (chordwise), m', 'Mass center (thickness), m',
                                      'Neutral axes (chordwise), m', 'Neutral axes (thickness), m',
                                      'Geometric center (chordwise), m', 'Geometric center (thickness), m',
                                      'Shear center (chordwise), m', 'Shear center (thickness), m'])
        for i in range(len(beam_mass_per_length)):  # receive number of radial sections
            beam_prop_writer.writerow([str(radial_stations[i]), str(beam_mass_per_length[i,0]),
                                       str(beam_mass_center[i,0]), str(beam_mass_center[i,1]),
                                       str(beam_neutral_axes[i,0]), str(beam_neutral_axes[i,1]),
                                       str(beam_geometric_center[i,0]), str(beam_geometric_center[i,1]),
                                       str(beam_shear_center[i,0]), str(beam_shear_center[i,1])])

    csv_file.close()
    # -------------------------------------------------- #

    # Export stiffness matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_stiff]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_stiff)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_stiff[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()

    # -------------------------------------------------- #
    # Export mass matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_mass]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_inertia)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_inertia[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()
    # -------------------------------------------------- #
    return None



def anbax_export_beam_struct_properties(folder_str, job_str, radial_stations, coordsys, solver, beam_stiff, beam_inertia, beam_mass_per_length):

    if solver=='anbax':
        export_name_general = 'anbax_beam_properties_general.csv'
        export_name_stiff = 'anbax_beam_properties_stiff_matrices.csv'
        export_name_mass = 'anbax_beam_properties_mass_matrices.csv'
    else:
        print('Define correct solver name (vabs or anbax) when calling export_beam_struct_properties')

    # -------------------------------------------------- #
    # Export mass per unit length for the defined radial stations
    if os.path.isdir(folder_str + 'csv_export/') == False:
        os.mkdir(folder_str + 'csv_export/')
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_general]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if coordsys == 'BeamDyn':
            beam_prop_writer.writerow(['Coordinate system:', 'Beamdyn coordinates'])
        elif coordsys == 'VABS/SONATA':
            beam_prop_writer.writerow(['Coordinate system:', 'VABS/SONATA coordinates'])
        elif coordsys == 'ANBAX':
            beam_prop_writer.writerow(['Coordinate system:', 'VABS/SONATA coordinates'])
        else:
            beam_prop_writer.writerow(['Coordinate system:', 'to be verified'])

        beam_prop_writer.writerow(['section in r/R', 'Mass per unit length [kg/m]'])
        for i in range(len(beam_mass_per_length)):  # receive number of radial sections
            beam_prop_writer.writerow([str(radial_stations[i]), str(beam_mass_per_length[i,0])])

    csv_file.close()
    # -------------------------------------------------- #

    # Export stiffness matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_stiff]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_stiff)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_stiff[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()

    # -------------------------------------------------- #
    # Export mass matrices for the defined radial stations
    with open(''.join([folder_str + 'csv_export/' + job_str[0:-5] + '_' + export_name_mass]), mode='w') as csv_file:
        beam_prop_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for i in range(len(beam_inertia)):  # receive number of radial sections
            beam_prop_writer.writerow(' ')
            beam_prop_writer.writerow(['section in r/R', str(radial_stations[i])])

            for j in range(6):  # number of rows for each matrix
                beam_prop_writer.writerow(beam_inertia[i, j, :])
                # beam_prop_writer.writerow(job.beam_properties[i, 1].TS[j, :])  # can eventually be called as a standalone via the job.beam_properties object
    csv_file.close()
    # -------------------------------------------------- #
    return None


def strain_energy_eval(blade, MatID=None):
    """
    Calculate the strain energy for a blade with already calculated
    internal loads.

    Parameters
    ----------
    blade : SONATA.classBlade.Blade
        A blade object from SONATA. The blade response for applied internal
        forces and moments should already be calculated.
    MatID : int or None, optional
        Which material should be post processed.
        Cells with other materials will be ignored.
        Input is not fully tested yet.
        If None, then all materials will be included in results.
        Materials can be checked by looking at `blade.materials`
        The default is None.

    Returns
    -------
    total_energy : float
        Total energy in the material across the full length of the blade.
    directional_energy : (6,) numpy.ndarray
        Energy contributions from each direction (see Notes for ordering).
    strain_all : list of (6,Ni) numpy.ndarray
        Strain values for each direction and relevant element.
        Each list entry corresponds to a different section each of which has
        a numpy.ndarray. The first index of the numpy.narray is the direction
        (see Notes) and the second is the cell. Not all cells may be included
        if `MatID` is used to filter elements. In that case NaN is returned
        for cells that do not match `MatID`.
    energy_all : list of (6,Ni) numpy.ndarray
        Same format as `strain_all`, but records the total energy associated
        with the corresponding elements.

    Notes
    -----

    Stress and strain are processed in the material coordinates.
    For no rotation of the material, directions are
    1 - axial to the blade,
    2 - along the arclength of a given section,
    3 - through the thickness of layers.

    Stress and strain directions of shear webs have not been verified.

    Output stress/strain/directional components are in the order
    [11, 22, 33, 23, 13, 12]

    Trapezoid integration is used to calculate total energy using sections
    provided in `blade`.

    If the blade definition used kg,m,s units, then the output energy will be
    in units of J.

    Strain energy is integrated with respect to the spanwise-length of the
    blade.
    It may be more accurate to integrate with respect to the arclength of the
    reference line.

    """

    length = blade.yml['components']['blade']['outer_shape_bem'] \
                    ['reference_axis']['z']['values'][-1]

    # Energy calculated in material coordinates.
    energyM_length = np.zeros((len(blade.sections), 6))

    strain_all = len(blade.sections) * [None]
    energy_all = len(blade.sections) * [None]

    stations_array = np.array([sec[0] for sec in blade.sections])

    quad_length_weights = np.zeros(len(blade.sections))
    quad_length_weights[0]    = 0.5*(stations_array[1] - stations_array[0])
    quad_length_weights[-1]   = 0.5*(stations_array[-1] - stations_array[-2])
    quad_length_weights[1:-1] = 0.5*(stations_array[2:] - stations_array[0:-2])

    for sec_ind,section in enumerate(blade.sections):

        (x,cs) = section
        cells = cs.mesh

        stressM = np.nan*np.zeros((len(cells), 6))
        strainM = np.nan*np.zeros((len(cells), 6))

        area = np.nan*np.zeros(len(cells))

        for ind,c in enumerate(cells):

            if c.MatID == MatID or MatID == None:
                # In Material Coordinates
                # Stresses
                stressM[ind, 0] = c.stressM.sigma11
                stressM[ind, 1] = c.stressM.sigma22
                stressM[ind, 2] = c.stressM.sigma33

                stressM[ind, 3] = c.stressM.sigma23
                stressM[ind, 4] = c.stressM.sigma13
                stressM[ind, 5] = c.stressM.sigma12

                # Strains
                strainM[ind, 0] = c.strainM.epsilon11
                strainM[ind, 1] = c.strainM.epsilon22
                strainM[ind, 2] = c.strainM.epsilon33

                strainM[ind, 3] = c.strainM.epsilon23
                strainM[ind, 4] = c.strainM.epsilon13
                strainM[ind, 5] = c.strainM.epsilon12

                # Area
                area[ind] = c.calc_area()

        # mask out all of the other materials
        mask = np.logical_not(np.isnan(area))

        stressM = stressM[mask, :]
        strainM = strainM[mask, :]
        area = area[mask]

        # In material coordinates
        energyM_density = 0.5*stressM*strainM

        # The shear components have a factor of 2 in the contraction
        energyM_density[:, 3:] *= 2

        # Energy per length along the beam
        energyM_length[sec_ind, :] = energyM_density.T @ area

        # Collect statistics for distribution of energy v. strain.
        strain_all[sec_ind] = strainM
        energy_all[sec_ind] = energyM_density * area.reshape(-1,1) \
                                    * quad_length_weights[sec_ind] * length

    directional_energy = np.trapz(energyM_length.T, stations_array) * length

    total_energy = np.sum(directional_energy)

    return total_energy, directional_energy, strain_all, energy_all
