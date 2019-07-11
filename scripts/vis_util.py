# Visualization utilities
# Author: Tyler
# Date: July 10, 2019

import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio
from field_util import *
from calc_util import *
from load_util import *

DEFAULT_FRAME_PERIOD = 0.25 # seconds
class image_set_manager:
    def __init__(self, w_dir):
        self.__im_set = []
        self.__wdir = w_dir
    
    def save_plt(self, fname, track=True):
        """
        Saves plots to the working directory
    
        Parameters
        ----------
        fname : str
            The name to be given to the image file
        track : bool
            True if this file should be tracked in the analysis set
        """
        name = fname + '.png'
        plt.savefig(os.path.join(self.__wdir, name))
        if track:
            self.__im_set.append(name)
    
    def end_set(self, frame_period=DEFAULT_FRAME_PERIOD, del_frames=True):
        """
        Makes a .gif animation of the current frames in the set
        """
        images = []
        for filename in self.__im_set:
            images.append(imageio.imread(os.path.join(self.__wdir, filename)))
        gif_name = os.path.join(self.__wdir, self.__im_set[-1][:-4] + '.gif')
        print(gif_name)
        imageio.mimsave(gif_name, images, duration=frame_period)
        # Delete individual frames to save space
        if del_frames:
            for filename in self.__im_set:
                os.remove(os.path.join(self.__wdir, filename))
        self.clear_set()
        
    def clear_set(self):
        """
        Clears the current image set
        """
        del self.__im_set[:]

def set_image_defaults():
    """
    Runs some one-time matplotlib configurations for image formatting
    """
    # Image size
    width = 8.5
    height = 11
    
    # One-time plot configuration
    plt.rcParams['figure.figsize'] = [width, height]
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5)

def plot_2d(
    data,
    im_mgr,
    const_elevation=False,
    theta=90,
    const_azimuth=False,
    phi=0,
    log=False,
    data_idx=-1,
    field='',
    arrows=False,
    bg_field=True,
    v_min=float('inf'),
    v_max=float('inf')
):
    """
    Plots scalar and vector field quantities.
    
    Parameters
    ----------
    data : list of pyPLUTO.pload
        List of data frames
    im_mgr : image_set_manager
        Manages image collection for animation
    const_elevation : bool
        If True, plots the density as a function of radius and elevation
        Exclusive to const_azimuth
    theta : double
        The elevation (latitude) to fix, if applicable. Given in degrees
    phi : double
        The longitude to fix, if applicable. Given in degrees
    const_azimuth : double
        If True, plots the density as a function of radius and azimuth
        Exclusive to const_elevation
    log : bool
        (Optional) Plots the logarithm for scalar fields (base 10)
    data_idx : int
        (Optional) Causes only the specified frame to be plotted
    field : str
        The field quantity to plot
            - "rho" for density
            - "v_phi" for azimuthal (zonal) velocity
            - "v_r" for radial velocity
            - "v_theta" for meridional velocity
            - And similarly for magnetic flux density (e.g. "b_phi", etc.)
            - And similarly for current (e.g. "j_theta", etc.)
            - And similar for resistivity (e.g. "eta_r", etc.)
    arrows : bool
        Plot field vectors if true
    bg_field : bool
        Sum the background magnetic field with the deviation field if true
    v_min : float
        Specifies what the colorbar's min color should be
    v_max : float
        Specifies what the colorbar's max color should be
    """
    if const_elevation == False and const_azimuth == False:
        const_elevation = True
    
    # Sanity check on inputs. Not trying to be exhaustive here, just
    # trying to reduce the probability that I make a mistake while using
    # this code
    plot_v = field == 'v_phi' or field == 'v_r' or field == 'v_theta'
    plot_b = field == 'b_phi' or field == 'b_r' or field == 'b_theta'
    plot_j = field == 'j_phi' or field == 'j_r' or field == 'j_theta'
    plot_eta = field == 'eta_phi' or field == 'eta_r' or field == 'eta_theta'
    is_scalar_field = field == 'rho'
    assert(field == '' or field == 'rho' or plot_v or plot_b or plot_j or plot_eta), \
        "field argument is invalid"
    assert(const_elevation != True or const_azimuth != True), \
        "EITHER the azimuth or elevation may be specified, not both"
    assert(not (field == '' and log == True)), "Cannot do log plot if no field quantity is specified"
    assert(data_idx <= len(data)), "data_idx must be <= len(data)"
    assert(not (is_scalar_field and arrows)), \
        "Cannot plot arrows for scalar field {0}".format(field)
    
    single_only = data_idx != -1
    nframes = len(data)
    nrows   = 1 if single_only else np.round(np.sqrt(nframes))
    ncols   = 1 if single_only else np.round(np.sqrt(nframes) + 1)
    theta_idx = (np.abs(data[0].x2 - theta * np.pi / 180)).argmin()
    phi_idx   = (np.abs(data[0].x3 - phi * np.pi / 180)).argmin()
    try:
        fig = plt.figure()
        for idx in range(nframes):
            if single_only and idx != data_idx:
                continue
            
            X = data[idx].x1
            # Set up quantities
            if not is_scalar_field:
                if plot_v:
                    U = get_field_qty("v_r", data[idx], bg_field=bg_field)
                elif plot_b:
                    U = get_field_qty("b_r", data[idx], bg_field=bg_field)
                elif plot_j:
                    U = get_field_qty("j_r", data[idx], bg_field=bg_field)
                elif plot_eta:
                    U = get_field_qty("eta_r", data[idx], bg_field=bg_field)
            if const_azimuth:
                Y = 90 - (180 / np.pi) * data[idx].x2
                if not is_scalar_field:
                    if plot_v:
                        V = get_field_qty("v_theta", data[idx], bg_field=bg_field)
                    if plot_b:
                        V = get_field_qty("b_theta", data[idx], bg_field=bg_field)
                    if plot_j:
                        V = get_field_qty("j_theta", data[idx], bg_field=bg_field)
                    if plot_eta:
                        V = get_field_qty("eta_theta", data[idx], bg_field=bg_field)
                    U = U[:,:,phi_idx].T
                    V = V[:,:,phi_idx].T
                    V = -1 * V # Needed for getting the theta vector direction right
                qty = get_field_qty(field, data[idx], bg_field=bg_field)[:,:,phi_idx].T
                annotations = [field, 'X1 (radius)', 'X2 (elevation)']
            elif const_elevation:
                Y = (180 / np.pi) * data[idx].x3
                if not is_scalar_field:
                    if plot_v:
                        V = get_field_qty("v_phi", data[idx], bg_field=bg_field)
                    if plot_b:
                        V = get_field_qty("b_phi", data[idx], bg_field=bg_field)
                    if plot_j:
                        V = get_field_qty("j_phi", data[idx], bg_field=bg_field)
                    if plot_eta:
                        V = get_field_qty("eta_phi", data[idx], bg_field=bg_field)
                    U = U[:,theta_idx,:].T
                    V = V[:,theta_idx,:].T
                qty = get_field_qty(field, data[idx], bg_field=bg_field)[:,theta_idx,:].T
                annotations = [field, 'X1 (radius)', 'X3 (azimuth)']
    
            # Plotting
            if single_only:
                plt.subplot(nrows, ncols, 1)
            else:
                plt.subplot(nrows, ncols, idx + 1)
            if field != '':
                if log:
                    annotations[0] = "Log10 " + annotations[0]
                    qty = np.log10(qty)
                if v_min == float('inf'):
                    v_min = np.floor(np.min(qty))
                if v_max == float('inf'):
                    v_max = np.ceil(np.max(qty))
                plt.pcolormesh(X, Y, qty, vmin=v_min, vmax=v_max)
                if (single_only or ((idx != 0 and (idx+1) % ncols == 0) or idx == nframes-1)):
                    plt.colorbar()
            if arrows:
                plt.quiver(X, Y, U, V)
            annotations[0] = annotations[0] + ' (time: ' + str(np.round(data[idx].SimTime, 2)) + ')'
            plt.title(annotations[0])
            plt.xlabel(annotations[1])
            if single_only or idx % ncols == 0:
                plt.ylabel(annotations[2])
            if single_only:
                break
    
        # Save image to working directory
        if const_azimuth:
            fname = "rad_and_theta_vs_" + field
        elif const_elevation:
            fname = "rad_and_phi_vs_" + field
        if arrows:
            if fname[-1] != '_':
                fname = fname + "_"
            fname = fname + "vector"
        fname = fname + "_max" + str(v_max)
        if bg_field and plot_b:
            fname = fname + "_bgfield"
        if not single_only:
            im_mgr.save_plt(fname, track=False)
        else:
            fname = fname + "_t" + str(np.round(data[idx].SimTime, 5))
            im_mgr.save_plt(fname)
        plt.close()
    except Exception as e:
        print(e)
        plt.close()

def subplot_radius_vs_eta_2d(data, theta, phi, log=False):
    """
    Plots eta profiles as a function of the radius for fixed theta and phi, for
    all time steps. Note that the values specified for theta and/or phi may not
    exist on the grid since it is discrete. The closest value will be used in
    such cases.

    Parameters
    ----------
    data : list of pyPLUTO.pload
        List of data frames
    theta : double
        The desired elevation angle (spherical coordinates)
    phi : double
        The desired azimuthal angle
    log : bool
        (Optional) Plots the logarithm of the density (base 10)
    """
    theta_idx = (np.abs(data[0].x2 - theta * np.pi / 180)).argmin()
    phi_idx = (np.abs(data[0].x3 - phi * np.pi / 180)).argmin()
    
    # Eta is only dumped out once at the beginning of the code, so no need to do
    # more than 1 plot
    try:
        fig = plt.figure()
        plt.subplot(1, 1, 1)
        qty = get_field_qty("eta_r")[:, theta_idx, phi_idx] # TODO (tyler): add way to specify which component...
        if log:
            qty = np.log10(qty)
        plt.scatter(data[0].x1, qty)
        plt.title('theta = ' + str(np.round(data[0].x2[theta_idx] * 180 / np.pi, decimals=1)) +
                 ', phi = ' + str(np.round(data[0].x3[phi_idx] * 180 / np.pi, decimals=1))
        )
        plt.xlabel('X1 (radius)')
        if log:
            plt.ylabel('Log 10 Resistivity [S/m]')
        else:
            plt.ylabel('Resistivity (SI units)')
    except Exception as e:
        print(e)
        plt.close()

def subplot_radius_vs_density_2d(data, theta, phi, log=False, first_only=False):
    """
    Plots density profiles as a function of the radius for fixed theta and phi,
    for all time steps. Note that the values specified for theta and/or phi may
    not exist on the grid since it is discrete. The closest value will be used
    in such cases.

    Parameters
    ----------
    data : list of pyPLUTO.pload
        List of data frames
    theta : double
        The desired elevation angle (spherical coordinates)
    phi : double
        The desired azimuthal angle
    log : bool
        (Optional) Plots the logarithm of the density (base 10)
    first_only : bool
        (Optional) Causes only the first frame to be plotted
    """
    theta_idx = (np.abs(data[0].x2 - theta * np.pi / 180)).argmin()
    phi_idx = (np.abs(data[0].x3 - phi * np.pi / 180)).argmin()
    
    nframes = 1 if first_only else len(data)
    nrows   = 1 if first_only else np.round(np.sqrt(nframes))
    ncols   = 1 if first_only else np.round(np.sqrt(nframes) + 1)
    for idx in range(nframes):
        plt.subplot(nrows, ncols, idx + 1)
        qty = data[idx].rho[:, theta_idx, phi_idx]
        if log:
            qty = np.log10(qty)
        plt.scatter(data[idx].x1, qty)
        plt.title('Time: ' + str(np.round(data[idx].SimTime, 2)) +
                  ' (theta = ' + str(np.round(data[idx].x2[theta_idx] * 180 / np.pi, decimals=1)) +
                  ', phi = ' + str(np.round(data[idx].x3[phi_idx] * 180 / np.pi, decimals=1)) + ')'
        )
        plt.xlabel('X1 (radius)')
        if log:
            plt.ylabel('Log 10 Density')
        else:
            plt.ylabel('Density')

def plot_extrema_over_time(data, im_mgr, field, bg_field=False):
    """
    Plots the min and max values for the specified field quantity, for each time
    step

    Parameters
    ----------
    data : pyPLUTO.pload
        PLUTO data frame
    im_mgr : image_set_manager
        Manages image collection for animation
    field : str
        Name of field
    bg_field : bool
        Enables the background magnetic field (if applicable) if True
    """
    try:
        fig, ax = plt.subplots()
        t = []
        field_max = []
        field_min = []
        for frame in data:
            t.append(frame.SimTime)
            qty = get_field_qty(field, frame, bg_field=bg_field)
            field_max.append(np.max(qty))
            field_min.append(np.min(qty))
        ax.scatter(t, field_max, c="blue", label="max")
        ax.scatter(t, field_min, c="red", label="min")
        ax.legend()
        plt.title(field)
        plt.xlabel('Time')
        plt.ylabel('{0}'.format(field)) # TODO (tyler): support adding units to the y-axis labels
        fname = field + "_extrema_vs_t"
        if bg_field:
            fname = fname + "_bgfield"
        im_mgr.save_plt(fname)
        print("Saved plot as " + fname)
        plt.close();
    except Exception as e:
        print(e)
        plt.close()

def plot_ohmic_heating_from_python(data, im_mgr, **kwargs):
    """
    Computes ohmic heating values for each time step directly from the data sets
    and plots it

    Parameters
    ----------
    data : np.ndarray
        Data frames to get the coordinate axes and current values from. Must be
        in code units
    im_mgr : image_set_manager
        Manages image collection for animation
    eta_field : np.ndarray
        Resistivity field array. MUST BE IN CODE UNITS!!
    eta_const : double
        Constant resistivity value. Must be in code units
    """
    assert('eta_field' in kwargs or 'eta_const' in kwargs), "Must specify either eta_field or eta_const"
    assert(not ('eta_field' in kwargs and 'eta_const' in kwargs)), "You must pass in only ONE of eta_field or eta_const"
    
    t = []
    I_set = []
    for frame in data[1:]:
        t.append(frame.SimTime)
        I = compute_ohmic_heating(frame, **kwargs)
        I_set.append(I)

    try:
        fig = plt.figure()
        plt.scatter(t, I_set)
        plt.title('Ohmic heating vs time')
        plt.xlabel('Time')
        plt.ylabel('Integral of $\eta|J|^2$ (code units)')
        fname = "ohmic_heating_vs_t_python"
        const_eta = 'eta_const' in kwargs
        if const_eta:
            fname += "_eta_const"
        else:
            fname += "_eta_field"
        im_mgr.save_plt(fname)
        print("Saved plot as " + fname)
        plt.close()
    except Exception as e:
        print(e)
        plt.close()

def plot_ohmic_heating_from_pluto(im_mgr, w_dir, fname="heating.dat"):
    """
    Loads ohmic heating data outputted by PLUTO and plots it

    Parameters
    ----------
    im_mgr : image_set_manager
        Manages image collection for animation
    w_dir : str
        The location where the file is to be loaded from
    fname : str
        Name of file to load ohmic heating data from
    """
    h = []
    t = []
    TIME_IDX = 1
    HEATING_IDX = 2
    try:
        print("Attempting to load {0}".format(fname))
        with open(os.path.join(w_dir, fname), "r") as f:
            f.next() # Skip first line since it's just labels for the data
            for line in f:
                step_time_heating = line.split(" ")
                t.append(float(step_time_heating[TIME_IDX]))
                h.append(float(step_time_heating[HEATING_IDX]))
    except Exception as e:
        print(e)
        print("An exception occurred while reading heating.dat. Are you sure"
              " the file exists?")
    try:
        fig = plt.figure()
        plt.scatter(t, h)
        plt.title('Ohmic heating vs time')
        plt.xlabel('Time')
        plt.ylabel('Integral of $\eta|J|^2$ (code units)')
        im_mgr.save_plt("ohmic_heating_vs_t_PLUTO")
        plt.close();
        print("Plotting complete")
    except Exception as e:
        print(e)
        plt.close()

# TODO(tyler): vectorize computations. Terrible efficiency at the moment
def pltSphData3D(d, qty, axes, r=-1, theta=-1, phi=-1, silent=True, f=-1,
    log=False, v_min=float('inf'), v_max=float('inf'), size=20
):
    '''
    Displays a 3D scatter plot of the specified field quantity with one of the
    coordinates held constant.

    Parameters
    ----------
    d : pyPLUTO.pload
        Interface to data object
    qty : str
        Name of field quantity (e.g. 'rho' or 'prs')
    axes : matplotlib.axes.SubplotBase
        Plotting interface
    r, theta, phi : int
        (Optional) Fixes a grid index for the radius, elevation, or azimuth. This
        quantity will be kept constant throughout the plot. Only one may be set at
        a time. Defaults to a plot of constant radius.
    silent : bool
        Prints info messages when True, otherwise no messages are printed
    log : bool
        (Optional) Plots the logarithm of the field quantity (base 10)
    v_min : float
        Specifies the min value for the color scale
    v_max : float
        Specifies the max value for the color scale
    size : float
        Specifies point size to pass into the plotter
    
    Example
    -------
    # This will use the pload object `dh` to plot the density for a constant radius.
    The radius setting of 20 fixes the radius at dh.x1[20].
    pltSphData3D(dh, 'rho', r=20)
    '''
    
    # We default to a constant radius (specifically, the innermost radius)
    if r == -1 and phi == -1 and theta == -1:
        r = 0
    
    # Sanity checks on inputs
    assert(
        (r >= 0  and phi == -1 and theta == -1) or
        (r == -1 and phi >= 0  and theta == -1) or
        (r == -1 and phi == -1 and theta >= 0 )
    ), "Only one of r, phi, and theta can be specified"
    
    plot_v = qty == 'v_phi' or qty == 'v_r' or qty =='v_theta'
    plot_b = qty == 'b_phi' or qty == 'b_r' or qty =='b_theta'
    assert(qty == 'rho' or plot_v or plot_b), "Invalid field quantity"
    
    # Coordinates
    #    x1 = radial
    #    x2 = latitudinal (theta)
    #    x3 = longitudinal (azimuthal, phi)
    x, y, z, c = [], [], [], []
    
    num_r = d.x1.shape[0]
    num_theta = d.x2.shape[0]
    num_phi = d.x3.shape[0]
    
    # Find the field quantity
    if qty == 'rho':
        field = d.rho
    elif plot_v:
        if qty == 'v_r':
            field = d.vx1
        elif qty == 'v_theta':
            field = d.vx2
        elif qty == 'v_phi':
            field = d.vx3
    elif plot_b:
        if qty == 'b_r':
            field = d.Bx1
        elif qty == 'b_theta':
            field = d.Bx2
        elif qty == 'b_phi':
            field = d.Bx3

    if log:
        field = np.log(field)
    
    # Plot the quantity
    if(r >= 0):
        assert(r <= num_r - 1), "r is too large (must be <= {0})".format(num_r - 1)
        r_l = d.x1[r]
        if not silent:
            print("Plotting sphere of radius {0}".format(r_l))
        for i in range(num_theta):
            theta_l = d.x2[i]
            for j in range(num_phi):
                phi_l = d.x3[j]
                x_p, y_p, z_p = sph2cart(r_l, theta_l, phi_l)
                x.append(x_p)
                y.append(y_p)
                z.append(z_p)
                c.append(field[r,i,j])
    elif(theta >= 0):
        assert(theta <= num_theta - 1), "theta is too large (must be <= {0})".format(num_theta - 1)
        theta_l = d.x2[theta]
        if not silent:
            print("Plotting cone with elevation angle {0} deg".format(theta_l * 180 / np.pi))
        for i in range(num_r):
            r_l = d.x1[i]
            for j in range(num_phi):
                phi_l = d.x3[j]
                x_p, y_p, z_p = sph2cart(r_l, theta_l, phi_l)
                x.append(x_p)
                y.append(y_p)
                z.append(z_p)
                c.append(field[i,theta,j])
    else:
        assert(phi <= num_phi - 1), "phi is too large (must be <= {0})".format(num_phi - 1)
        phi_l = d.x3[phi]
        if not silent:
            print("Plotting disk with azimuth {0} deg".format(phi_l * 180 / np.pi))
        for i in range(num_r):
            r_l = d.x1[i]
            for j in range(num_theta):
                theta_l = d.x2[j]
                x_p, y_p, z_p = sph2cart(r_l, theta_l, phi_l)
                x.append(x_p)
                y.append(y_p)
                z.append(z_p)
                c.append(field[i,j,phi])
    
    plt.xlabel('x')
    plt.ylabel('y')
    
    # Update the plot axes, but don't make them smaller
    axes.set_xlim3d(
        min(axes.get_xlim3d()[0], min(x)),
        max(axes.get_xlim3d()[1], max(x))
    )
    axes.set_ylim3d(
        min(axes.get_ylim3d()[0], min(y)),
        max(axes.get_ylim3d()[1], max(y))
    )
    axes.set_zlim3d(
        min(axes.get_zlim3d()[0], min(z)),
        max(axes.get_zlim3d()[1], max(z))
    )
    
    # Compute the max and min values used for the color
    if v_min == float('inf'):
        v_min = np.floor(np.min(field))
    if v_max == float('inf'):
        v_max = np.ceil(np.max(field))
    axes.scatter3D(x, y, z, c=c, vmin=v_min, vmax=v_max, s=size)