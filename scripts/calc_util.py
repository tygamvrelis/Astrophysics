# Calculation utilities
# Author: Tyler
# Date: July 10, 2019

import numpy as np
from field_util import get_field_qty
from settings import *

# TODO(tyler): vectorize computations. Terrible efficiency at the moment
def compute_ohmic_heating(frame, **kwargs):
    """
    Computes ohmic heating values for each time step directly from the data sets

    Parameters
    ----------
    frame : pyPLUTO.pload
        Data frame to get the coordinate axes and current values from. Must
        be in code units
    eta_field : np.ndarray
        Resistivity field array. MUST BE IN CODE UNITS!!
    eta_const : double
        Constant resistivity value. Must be in code units
    """
    assert('eta_field' in kwargs or 'eta_const' in kwargs), \
        "Must specify either eta_field or eta_const"
    assert(not ('eta_field' in kwargs and 'eta_const' in kwargs)), \
        "You must pass in only ONE of eta_field or eta_const"
    
    const_eta = 'eta_const' in kwargs
    if not const_eta:
        eta_f = kwargs['eta_field']
        eta_r = eta_f.ex1
        eta_theta = eta_f.ex2
        eta_phi = eta_f.ex3
    else:
        eta_c = kwargs['eta_const']

    r     = frame.x1
    theta = frame.x2
    phi   = frame.x3
    j_r     = get_field_qty('j_r', frame)
    j_theta = get_field_qty('j_theta', frame)
    j_phi   = get_field_qty('j_phi', frame)
    r_tot = frame.n1_tot
    theta_tot = frame.n2_tot
    phi_tot = frame.n3_tot
    I = 0
    for i in range(r_tot - 1):
        dVr = np.abs(r[i + 1]**3 - r[i]**3) / 3.0
        for j in range(theta_tot - 1):
            dmu = np.abs(np.cos(theta[j]) - np.cos(theta[j + 1]))
            for k in range(phi_tot - 1):
                # Volume element. Partially calculated before to reduce number
                # of repeated calculations
                d_phi = phi[k + 1] - phi[k]
                vol = dVr * dmu * d_phi # Volume element formula from PLUTO

                # Field quantity, choosing middle point in the spherical wedge
                j_r_s     = (j_r[i+1,j+1,k+1]     + j_r[i,j,k])     / 2
                j_theta_s = (j_theta[i+1,j+1,k+1] + j_theta[i,j,k]) / 2
                j_phi_s   = (j_phi[i+1,j+1,k+1]   + j_phi[i,j,k])   / 2

                # Update integral
                if const_eta:
                    j_sq = j_r_s**2 + j_theta_s**2 + j_phi_s**2
                    I += eta_c * j_sq * vol
                else:
                    eta_r_s     = (eta_r[i+1,j+1,k+1]     + eta_r[i,j,k])     / 2
                    eta_theta_s = (eta_theta[i+1,j+1,k+1] + eta_theta[i,j,k]) / 2
                    eta_phi_s   = (eta_phi[i+1,j+1,k+1]   + eta_phi[i,j,k])   / 2
                    eta_j_sq = eta_r_s*(j_r_s**2) + \
                               eta_theta_s*(j_theta_s**2) + \
                               eta_phi_s*(j_phi_s**2)
                    I += eta_j_sq * vol
    return I

def sph2cart(r, theta, phi):
    """
    Converts the inputted spherical coordinates to cartesian coordinates

    Parameters
    ----------
    r : double
        The radius
    theta : double
        The altitude
    phi : double
        The azimuth
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z
