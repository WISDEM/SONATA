"""
Utility file for tesing
"""

import numpy as np


def load_bd_blade(fname):
    """
    Load matrices generated in beamdyn blade file to check in tests.

    Parameters
    ----------
    fname : str
        Filename and/or path of beamdyn blade file to load.

    Returns
    -------
    mass : list of (6,6) numpy.ndarrays
        Mass matrices for each station.
    stiff : list of (6,6) numpy.ndarrays
        Stiffness matrices for each station.

    """
    
    # Get number of stations
    with open(fname, 'r') as file:
        lines = file.readlines()
       
    n_stations = int(lines[3].split()[0])
    
    mass = n_stations * [None]
    stiff = n_stations * [None]
    
    for i in range(n_stations):

        stiff[i] = np.genfromtxt(lines[11+i*15:17+i*15])

        mass[i] = np.genfromtxt(lines[18+i*15:24+i*15])

    return mass, stiff

def compare_bd_blade(ref_path, test_path, tolerance=1e-9):
    """
    

    Parameters
    ----------
    ref_path : filepath
        Reference beamdyn blade file.
    test_path : filepath
        Test beamdyn blade file.
    tolerance : float
        This is mulitiplied by the largest element of the 6x6 matrix to set
        `atol` in `np.allclose`

    Returns
    -------
    None.

    """
    
    mass_ref, stiff_ref = load_bd_blade(ref_path)
    mass_test, stiff_test = load_bd_blade(test_path)
    
    for i in range(len(mass_ref)):
        
        assert np.allclose(stiff_ref[i], stiff_test[i],
                           atol=tolerance*stiff_ref[i].max()), \
            "Stiffness matrix does not match at station index {:}.".format(i)

        print("Stiffness error: {:}".format(
            np.abs(stiff_ref[i]-stiff_test[i]).max() / stiff_ref[i].max()))

        assert np.allclose(mass_ref[i], mass_test[i],
                           atol=tolerance*mass_ref[i].max()), \
            "Mass matrix does not match at station index {:}.".format(i)

        print("Mass error: {:}".format(
            np.abs(mass_ref[i]-mass_test[i]).max() / mass_ref[i].max()))
        
    return
    