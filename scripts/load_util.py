# Utilities for loading data outputted in a PLUTO run
# Author: Tyler
# Date: July 10, 2019

import os
import sys
import pyPLUTO as pp
import numpy as np
from settings import *

def load_user_params(w_dir):
    """
    Loads user-defined parameters into the settings

    Parameters
    ----------
    w_dir : str
        The location where the file is to be loaded from
    """
    print("Loading user params")
    params = []
    try:
        with open(os.path.join(w_dir, "user_params.dat"), "r") as f:
            for line in f:
                params.append(float(line))
    except IOError:
        print("No user params found")
        return
    
    user_params = {}
    user_params['ALPHA']     = None
    user_params['VMAX']      = None
    user_params['EXP_DECAY'] = None
    user_params['TRELAX']    = None
    user_params['BSURFACE']  = None
    user_params['ETA']       = None
    
    for i in range(len(params)):
        if i == 0:
            user_params['ALPHA'] = params[i]
        elif i == 1:
            user_params['VMAX'] = params[i]
        elif i == 2:
            user_params['EXP_DECAY'] = params[i]
        elif i == 3:
            user_params['TRELAX'] = params[i]
        elif i == 4:
            user_params['BSURFACE'] = params[i]
        elif i == 5:
            user_params['ETA'] = params[i]
    
    print("User params:")
    for k,v in user_params.items():
        print("\t" + k + ": " + str(v))
    
    settings.user_params = user_params
    return user_params

def load_unit_constants(w_dir):
    """
    Loads the scaling constants from a file into the settings

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
    except IOError as e:
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
    
    settings.UNIT_DENSITY = UNIT_DENSITY
    settings.UNIT_VELOCITY = UNIT_VELOCITY
    settings.UNIT_LENGTH = UNIT_LENGTH
    return UNIT_DENSITY, UNIT_VELOCITY, UNIT_LENGTH

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

def load_vector_field_3d(ref_frame, w_dir, fname="eta_field.dat"):
    """
    Loads a vector field from a file into a numpy array. Order of iteration is
    k-dir, j-dir, i-dir, and each line contains the field components in order
    i-j-k.

    Parameters
    ----------
    ref_frame : pyPLUTO.pload
        Reference data frame, used for getting coordinates etc.
    w_dir : str
        The location where the file is to be loaded from
    fname : str
        The name of the file to load the resistivity field from
    """
    r_tot     = ref_frame.n1_tot
    theta_tot = ref_frame.n2_tot
    phi_tot   = ref_frame.n3_tot
    try:
        with open(os.path.join(w_dir, fname), "r") as f:
            data = f.readlines()
    except Exception as e:
        print(e)
        return
    expected_num_lines = r_tot * theta_tot * phi_tot
    assert(len(data) == expected_num_lines), (
        "Error reading {0} (got {1} lines, expected {2})".format(
            fname, len(data), expected_num_lines
        )
    )
    fx1 = np.zeros((r_tot, theta_tot, phi_tot))
    fx2 = np.zeros((r_tot, theta_tot, phi_tot))
    fx3 = np.zeros((r_tot, theta_tot, phi_tot))
    for k in range(phi_tot):
        for j in range(theta_tot):
            for i in range(r_tot):
                idx = k * theta_tot * r_tot + j * r_tot + i
                fijk = data[idx].split(" ")
                fijk = [float(e) for e in fijk]
                assert(len(fijk) == 3), \
                    "Wrong number of field components (i,j,k)=({0},{1},{2})".format(
                        i, j, k
                    )
                fx1[i,j,k] = fijk[0]
                fx2[i,j,k] = fijk[1]
                fx3[i,j,k] = fijk[2]
    return fx1, fx2, fx3