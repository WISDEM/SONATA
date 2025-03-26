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
    