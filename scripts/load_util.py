# Utilities for loading data outputted in a PLUTO run
# Author: Tyler
# Date: July 10, 2019

import os
import sys
import pyPLUTO as pp
import astropy.constants as const

def load_user_params(w_dir):
    """
    Loads user-defined parameters.

    Parameters
    ----------
    w_dir : str
        The location where the file is to be loaded from
    """
    print("Loading user params")
    params = []
    with open(os.path.join(w_dir, "user_params.dat"), "r") as f:
        for line in f:
            params.append(float(line))
    return params

def load_unit_constants(w_dir):
    """
    Loads the scaling constants from a file.

    Parameters
    ----------
    w_dir : str
        The location where the file is to be loaded from
    """
    UNIT_DENSITY  = sys.maxint
    UNIT_VELOCITY = sys.maxint
    UNIT_LENGTH   = sys.maxint
    try:
        print("Attempting to load unit constants")
        with open(os.path.join(w_dir, "unit_constants.dat"), "r") as f:
            for line in f:
                l = line.split(" ")
                assert(len(l) == 2), "Error reading unit_constants.dat ({0})".format(l)
                if l[0].startswith("UNIT_DENSITY"):
                    UNIT_DENSITY = float(l[1])
                elif l[0].startswith("UNIT_VELOCITY"):
                    UNIT_VELOCITY = float(l[1])
                elif l[0].startswith("UNIT_LENGTH"):
                    UNIT_LENGTH = float(l[1])
    except IOError:
        print(e)
        print("WARNING: Using default unit constants (not loaded from run)")
        UNIT_DENSITY  = 1.0e-9
        UNIT_VELOCITY = 1.0e5
        UNIT_LENGTH   = 1.6*6.9911e9 / 1.2

    print("Using unit constants:\n\tUNIT_DENSITY = {0} [g/cm]\n"
          "\tUNIT_VELOCITY = {1} [cm/s]\n"
          "\tUNIT_LENGTH = {2} [cm]\n".format(
            UNIT_DENSITY, UNIT_VELOCITY, UNIT_LENGTH
        )
    )

    assert(UNIT_DENSITY != sys.maxint and
           UNIT_VELOCITY != sys.maxint and
           UNIT_LENGTH != sys.maxint), "Your unit constants don't make sense!"
    return UNIT_DENSITY, UNIT_VELOCITY, UNIT_LENGTH

def get_unit_constants():
    

def load_pluto_dataframes(w_dir):
    """
    Loads all the PLUTO dataframes in the specified directory into a list

    Parameters
    ----------
    w_dir : str
        The location where the files are to be loaded from
    """
    n_frames = len([file for file in os.listdir(w_dir) if "data" in file and ".dbl" in file])
    data = list() # List of data handles
    for i in range(n_frames):
        d = pp.pload(i, w_dir)
        data.append(d)
    return data