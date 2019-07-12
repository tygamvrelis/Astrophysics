# Scalar and vector field utilities
# Author: Tyler
# Date: July 10, 2019

import numpy as np
import astropy.constants as const
import astropy.units as u
from load_util import load_vector_field_3d
from settings import *

# TODO (tyler): use astropy.units for conversions as much as possible
def eta_cgs_to_si(e):
    '''
    Converts resistivity from cgs units to SI units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in cgs units
    '''
    return (e * 10**11 / const.c.to('cm/s')**2).value * u.S / u.m

def eta_si_to_cgs(e):
    '''
    Converts resistivity from SI units to cgs units

    Parameters
    ----------
    e : np.ndarray
        Resistivity field in SI units
    '''
    return (e / eta_cgs_to_si(1)).value * u.cm**2 / u.s

def get_field_qty(field, df=None, bg_field=True):
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

def get_field_units(field):
    """
    Wrapper for retrieving field quantity units in cgs.
    
    Returns a tuple whose first element indicates whether or not the field
    argument has a corresponding unit (True), and if True, returns the unit as
    the second element

    Parameters
    ----------
    field : str
        Name of field
    """
    get_v = field == 'v_phi' or field == 'v_r' or field == 'v_theta'
    get_b = field == 'b_phi' or field == 'b_r' or field == 'b_theta'
    get_j = field == 'j_phi' or field == 'j_r' or field == 'j_theta'
    get_eta = field == 'eta_phi' or field == 'eta_r' or field == 'eta_theta'
    if not (field == 'rho' or get_v or get_b or get_j or get_eta):
        return False, u.dimensionless_unscaled
    
    # Find the field quantity
    if field == 'rho':
        units = u.g / u.cm**3
    elif get_v:
        units = u.cm / u.s
    elif get_b:
        units = u.cm**(-1/2) * u.g * u.s**(-1)
    elif get_j:
        units = u.g**(1/2) * u.cm**(-3/2) * u.s**(-1)
    elif get_eta:
        units = u.cm**2 / u.s
    return True, units

def get_scaling_factor(field):
    """
    Gets a scaling factor that, when multiplied with the field quantity in
    dimensionless units, yields the field quantity in cgs units

    Returns a tuple whose first element indicates whether or not the field
    argument requires scaling to cgs (True), and if True, returns the scaling
    factor as the second element
    
    Parameters
    ----------
    field : str
        Name of field
    """
    get_v = field == 'v_phi' or field == 'v_r' or field == 'v_theta'
    get_b = field == 'b_phi' or field == 'b_r' or field == 'b_theta'
    get_j = field == 'j_phi' or field == 'j_r' or field == 'j_theta'
    get_eta = field == 'eta_phi' or field == 'eta_r' or field == 'eta_theta'
    get_time = field == 'SimTime' # TODO (tyler): check this...
    if not (field == 'rho' or get_v or get_b or get_j or get_eta or get_time):
        return False, 1.0
    
    # Find the field quantity
    if field == 'rho':
        factor = settings.UNIT_DENSITY
    elif get_v:
        factor = settings.UNIT_VELOCITY
    elif get_time:
        factor = settings.UNIT_TIME
    elif get_b:
        # Converts magnetic field from code units to cgs units (Gauss) as per
        # section 5.1.1 in the PLUTO user manual
        factor = settings.UNIT_VELOCITY * np.sqrt(4 * np.pi * settings.UNIT_DENSITY)
    elif get_j:
        # Converts current from code units to cgs units
        factor = settings.UNIT_DENSITY**(1/2) / settings.UNIT_TIME
    elif get_eta:
        # Converts resistivity from code units to cgs units
        factor = settings.UNIT_LENGTH * settings.UNIT_VELOCITY
    return True, factor

def convert_pluto_dataframes_to_cgs(data):
    """
    Scales the quantities in the dataframes to cgs and attaches units
    """
    for frame in data:
        for kv in frame.__dict__.items():
            # 1. Scale
            needs_scaling, scaling_factor = get_scaling_factor(kv[0])
            if needs_scaling:
                kv[1] *= scaling_factor
            
            # 2. Add units to all quantities
            has_units, units = get_field_units(kv[0])
            if has_units:
                kv[1] *= units

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
        Generates a dipole field with the specified surface value in Gauss

        Parameters
        ----------
        ref_frame : pyPLUTO.pload
            Reference data frame, used for getting coordinates etc.
        B_surface : double
            Value of B field at the edge of the atmosphere, in Gauss
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
        # Gauss
        self.Bx1 *= get_field_units('b_r')
        self.Bx2 *= get_field_units('b_theta')
        self.Bx3 *= get_field_units('b_phi')
    
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
        self.Bx1 = b_code_units_to_cgs(Bx1)
        self.Bx2 = b_code_units_to_cgs(Bx2)
        self.Bx3 = b_code_units_to_cgs(Bx3)

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
