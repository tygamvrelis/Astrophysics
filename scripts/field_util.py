# Scalar and vector field utilities
# Author: Tyler
# Date: July 10, 2019

import numpy as np
import astropy.constants as const
from load_util import load_vector_field_3d
from settings import *

def get_field_qty(field, df=None, bg_field=False):
    """
    Wrapper for retrieving field quantities from a PLUTO data frame

    Parameters
    ----------
    field : str
        Name of field
    df : pyPLUTO.pload
        PLUTO data frame
    bg_field : bool
        Enables the background magnetic field (if applicable) if True
    """
    plot_v = field == 'v_phi' or field == 'v_r' or field == 'v_theta'
    plot_b = field == 'b_phi' or field == 'b_r' or field == 'b_theta'
    plot_j = field == 'j_phi' or field == 'j_r' or field == 'j_theta'
    plot_eta = field == 'eta_phi' or field == 'eta_r' or field == 'eta_theta'
    assert(field == 'rho' or plot_v or plot_b or plot_j or plot_eta), \
        "field argument is invalid"
    assert(not (df == None and (plot_v or plot_b or plot_j))), \
        "Must specify data frame"
    
    # Find the field quantity
    if field == 'rho':
        qty = df.rho
    elif plot_v:
        if field == 'v_r':
            qty = df.vx1
        elif field == 'v_theta':
            qty = df.vx2
        elif field == 'v_phi':
            qty = df.vx3
    elif plot_b:
        # Accept both conventions since this has changed over time
        try:
            if field == 'b_r':
                qty = df.Bx1
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx1
            elif field == 'b_theta':
                qty = df.Bx2
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx2
            elif field == 'b_phi':
                qty = df.Bx3
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx3
        except:
            if field == 'b_r':
                qty = df.bx1
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx1
            elif field == 'b_theta':
                qty = df.bx2
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx2
            elif field == 'b_phi':
                qty = df.bx3
                if bg_field:
                    qty = np.copy(qty) + settings.B0.Bx3
    elif plot_j:
        # Accept both conventions since this has changed over time
        try:
            if field == 'j_r':
                qty = df.Jx1
            elif field == 'j_theta':
                qty = df.Jx2
            elif field == 'j_phi':
                qty = df.Jx3
        except:
            if field == 'j_r':
                qty = df.jx1
            elif field == 'j_theta':
                qty = df.jx2
            elif field == 'j_phi':
                qty = df.jx3
    elif plot_eta:
        if field == 'eta_r':
            qty = settings.eta.ex1
        elif field == 'eta_theta':
            qty = settings.eta.ex2
        elif field == 'eta_phi':
            qty = settings.eta.ex3
    return qty;

class EtaField:
    """
    Resistivity field (vector)
    """
    def __init__(self, fx1=0, fx2=0, fx3=0):
        self.ex1 = fx1
        self.ex2 = fx2
        self.ex3 = fx3

    def load_field(self, ref_frame, w_dir, fname="eta_field.dat"):
        """
        Loads the resistivity field from a file. Order of iteration is k-dir, j-dir,
        i-dir, and each line contains the field components in order i-j-k.
        
        Also converts to physical units (SI) from code units.

        Parameters
        ----------
        ref_frame : pyPLUTO.pload
            Reference data frame, used for getting coordinates etc.
        w_dir : str
            The location where the file is to be loaded from
        fname : str
            The name of the file to load the resistivity field from
        """
        ex1, ex2, ex3 = load_vector_field_3d(ref_frame, w_dir, fname)
        self.ex1 = eta_cgs_to_si(get_physical_eta_units(ex1))
        self.ex2 = eta_cgs_to_si(get_physical_eta_units(ex2))
        self.ex3 = eta_cgs_to_si(get_physical_eta_units(ex3))
    
    def assign_field(self, etax1, etax2, etax3):
        """
        Allows existing arrays to be assigned to the internal field.
        
        Also converts to physical units (SI) from code units.

        Parameters
        ----------
        etaxi : np.ndarray
            Field quantity about i = 1, 2, 3 coordinates
        """
        self.ex1 = eta_cgs_to_si(get_physical_eta_units(etax1))
        self.ex2 = eta_cgs_to_si(get_physical_eta_units(etax2))
        self.ex3 = eta_cgs_to_si(get_physical_eta_units(etax3))
        
    def __truediv__(self, other):
        """
        Division by scalar
        """
        return EtaField(self.ex1 / other, self.ex2 / other, self.ex3 / other)

    def __div__(self, other):
        """
        Division by scalar
        """
        return EtaField(self.ex1 / other, self.ex2 / other, self.ex3 / other)

    def __mul__(self, other):
        """
        Multiplication by scalar
        """
        return EtaField(self.ex1 * other, self.ex2 * other, self.ex3 * other)

class MagField:
    """
    Magnetic field (vector)
    """
    def __init__(self, fx1=0, fx2=0, fx3=0):
        self.Bx1 = fx1
        self.Bx2 = fx2
        self.Bx3 = fx3
        
    def init_dipole_field(self, ref_frame, B_surface):
        """
        Generates a dipole field with the specified surface value.

        Parameters
        ----------
        ref_frame : pyPLUTO.pload
            Reference data frame, used for getting coordinates etc.
        B_surface : double
            Value of B field at the edge of the atmosphere
        """
        self.Bx1 = np.zeros_like(ref_frame.Bx1)
        self.Bx2 = np.zeros_like(ref_frame.Bx2)
        self.Bx3 = np.zeros_like(ref_frame.Bx3)
        M = B_surface * (ref_frame.x1r[-1] ** 3)
        r_tot     = ref_frame.n1_tot
        theta_tot = ref_frame.n2_tot
        phi_tot   = ref_frame.n3_tot
        for i in range(r_tot):
            r_cubed = ref_frame.x1[i] ** 3
            for j in range(theta_tot):
                bg_Bx1 = 2.0 * M * np.cos(ref_frame.x2[j]) / r_cubed
                bg_Bx2 = M * np.sin(ref_frame.x2[j]) / r_cubed
                for k in range(phi_tot):
                    self.Bx1[i,j,k] = bg_Bx1
                    self.Bx2[i,j,k] = bg_Bx2
    
    def load_field(self, ref_frame, w_dir, fname="bg_field.dat"):
        """
        Loads the background field from a file. Order of iteration is k-dir, j-dir,
        i-dir, and each line contains the field components in order i-j-k.
        
        Also converts to physical units (cgs) from code units.

        Parameters
        ----------
        ref_frame : pyPLUTO.pload
            Reference data frame, used for getting coordinates etc.
        w_dir : str
            The location where the file is to be loaded from
        fname : str
            The name of the file to load the background field from
        """
        Bx1, Bx2, Bx3 = load_vector_field_3d(ref_frame, w_dir, fname)
        self.Bx1 = get_physical_b_units(Bx1)
        self.Bx2 = get_physical_b_units(Bx2)
        self.Bx3 = get_physical_b_units(Bx3)
    
    def assign_field(self, B0x1, B0x2, B0x3):
        """
        Allows existing arrays to be assigned to the internal field.
        
        Also converts to physical units (SI) from code units.

        Parameters
        ----------
        B0xi : np.ndarray
            Field quantity about i = 1, 2, 3 coordinates
        """
        self.Bx1 = get_physical_b_units(B0x1)
        self.Bx2 = get_physical_b_units(B0x2)
        self.Bx3 = get_physical_b_units(B0x3)
    
    def __truediv__(self, other):
        """
        Division by scalar
        """
        return MagField(self.Bx1 / other, self.Bx2 / other, self.Bx3 / other)

    def __div__(self, other):
        """
        Division by scalar
        """
        return MagField(self.Bx1 / other, self.Bx2 / other, self.Bx3 / other)

    def __mul__(self, other):
        """
        Multiplication by scalar
        """
        return MagField(self.Bx1 * other, self.Bx2 * other, self.Bx3 * other)

# TODO (tyler): use astropy.units for conversions as much as possible
def get_physical_b_units(b):
    '''
    Converts magnetic field from code units to Gauss as per section 5.1.1 in the
    PLUTO user manual

    Parameters
    ----------
    b : np.ndarray
        Magnetic field in code units
    '''
    
    return b * settings.UNIT_VELOCITY * np.sqrt(4 * np.pi * settings.UNIT_DENSITY)

def get_physical_eta_units(e):
    '''
    Converts resistivity from code units to cgs units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in code units
    '''
    return e * settings.UNIT_LENGTH * settings.UNIT_VELOCITY

def get_eta_code_units(e):
    '''
    Converts resistivity from cgs units to code units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in cgs units
    '''
    return e / get_physical_eta_units(1)

def eta_cgs_to_si(e):
    '''
    Converts resistivity from cgs units to SI units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in cgs units
    '''
    return e * 10**11 / const.c.to('cm/s').value**2

def eta_si_to_cgs(e):
    '''
    Converts resistivity from SI units to cgs units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in SI units
    '''
    return e / eta_cgs_to_si(1)