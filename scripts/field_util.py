# Scalar and vector field utilities
# Author: Tyler
# Date: July 10, 2019

import numpy

def get_field_qty(field, df, bg_field=True):
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
    assert(field == 'rho' or plot_v or plot_b or plot_j or plot_eta), "field argument is invalid"
    
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
                qty = np.copy(qty) + B0.Bx1
        elif field == 'b_theta':
            qty = df.Bx2
            if bg_field:
                qty = np.copy(qty) + B0.Bx2
        elif field == 'b_phi':
            qty = df.Bx3
            if bg_field:
                qty = np.copy(qty) + B0.Bx3
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
            qty = eta.ex1
        elif field == 'eta_theta':
            qty = eta.ex2
        elif field == 'eta_phi':
            qty = eta.ex3
    return qty;

# TODO (tyler): make EtaField and MagField derive from a common vector field class
# TODO (tyler): create a load_vector_field function in load_util and call it here
class EtaField:
    def __init__(self, _ex1=0, _ex2=0, _ex3=0):
        self.ex1 = _ex1
        self.ex2 = _ex2
        self.ex3 = _ex3

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
        r_tot     = ref_frame.n1_tot
        theta_tot = ref_frame.n2_tot
        phi_tot   = ref_frame.n3_tot
        with open(os.path.join(w_dir, fname), "r") as f:
            data = f.readlines()
        expected_num_lines = r_tot * theta_tot * phi_tot
        assert(len(data) == expected_num_lines), "Error reading {0} (got {1} lines, expected {2})".format(fname, len(data), expected_num_lines)
        self.ex1 = np.zeros((r_tot, theta_tot, phi_tot))
        self.ex2 = np.zeros((r_tot, theta_tot, phi_tot))
        self.ex3 = np.zeros((r_tot, theta_tot, phi_tot))
        for k in range(phi_tot):
            for j in range(theta_tot):
                for i in range(r_tot):
                    idx = k * theta_tot * r_tot + j * r_tot + i
                    etaijk = data[idx].split(" ")
                    etaijk = [float(e) for e in etaijk]
                    assert(len(etaijk) == 3), "Wrong number of field components"
                    self.ex1[i,j,k] = etaijk[0]
                    self.ex2[i,j,k] = etaijk[1]
                    self.ex3[i,j,k] = etaijk[2]
        self.ex1 = eta_cgs_to_si(get_physical_eta_units(self.ex1))
        self.ex2 = eta_cgs_to_si(get_physical_eta_units(self.ex2))
        self.ex3 = eta_cgs_to_si(get_physical_eta_units(self.ex3))
        
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
    def __init__(self):
        self.Bx1 = 0
        self.Bx2 = 0
        self.Bx3 = 0
        
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
        
        Also converts to physical units (c.g.s.) from code units.

        Parameters
        ----------
        ref_frame : pyPLUTO.pload
            Reference data frame, used for getting coordinates etc.
        w_dir : str
            The location where the file is to be loaded from
        fname : str
            The name of the file to load the background field from
        """
        r_tot     = ref_frame.n1_tot
        theta_tot = ref_frame.n2_tot
        phi_tot   = ref_frame.n3_tot
        with open(os.path.join(w_dir, fname), "r") as f:
            data = f.readlines()
        expected_num_lines = r_tot * theta_tot * phi_tot
        assert(len(data) == expected_num_lines), "Error reading {0} (got {1} lines, expected {2})".format(fname, len(data), expected_num_lines)
        self.Bx1 = np.zeros_like(ref_frame.Bx1)
        self.Bx2 = np.zeros_like(ref_frame.Bx2)
        self.Bx3 = np.zeros_like(ref_frame.Bx3)
        for k in range(phi_tot):
            for j in range(theta_tot):
                for i in range(r_tot):
                    idx = k * theta_tot * r_tot + j * r_tot + i
                    B0ijk = data[idx].split(" ")
                    B0ijk = [float(e) for e in B0ijk]
                    assert(len(B0ijk) == 3), "Wrong number of field components"
                    self.Bx1[i,j,k] = B0ijk[0]
                    self.Bx2[i,j,k] = B0ijk[1]
                    self.Bx3[i,j,k] = B0ijk[2]
        self.Bx1 = get_physical_b_units(self.Bx1)
        self.Bx2 = get_physical_b_units(self.Bx2)
        self.Bx3 = get_physical_b_units(self.Bx3)
