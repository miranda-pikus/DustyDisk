import matplotlib.pyplot as plt
import numpy as np
from DustyDisk.Functions import Grid
import DustyDisk.Constants as Constants
#plt.style.use('./PlotStyling.mplstyle')

def PlotGasDensity(whichGrid, which_ax, color='purple',label='gas density'):
    '''
    Args:
        whichGrid (Grid) : the grid class to plot
        which_ax (ax) : axes object to put plot on
        color (string): color of line
        label (string): label of line
    '''
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.sigma_gas, 
                  color=color, label=label)
    which_ax.set_ylabel(r'density (g cm$^{-3}$)')

def PlotGasPressure(whichGrid, which_ax, color='orange',label='gas pressure'):
    '''
    Args:
        whichGrid (Grid) : the grid class to plot
        which_ax (ax) : axes object to put plot on
        color (string): color of line
        label (string): label of line
    '''
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.Pressure, 
                  color=color, label=label)
    which_ax.set_ylabel(r'pressure (dyne)')

def PlotDriftVelocity(whichGrid, which_ax, color='forestgreen', label='drift velocity'):
    '''
    Args:
        whichGrid (Grid) : the grid class to plot
        which_ax (ax) : axes object to put plot on
        color (string): color of line
        label (string): label of line
    '''
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.vdrift(), 
                  color=color, label=label)
    which_ax.set_ylabel(r'v$_{drift}$ (cm/s)')

def PlotDustDensity(whichGrid, which_ax, color='black', label='dust density'):
    '''
    Args:
        whichGrid (Grid) : the grid class to plot
        which_ax (ax) : axes object to put plot on
        color (string): color of line
        label (string): label of line
    '''
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.dust_density(), 
                    color=color, label=label)
    which_ax.set_ylabel(r'norm. dust density (g cm$^{-3}$)')


def PlotQuantity(theGrid, which_Qs):
    '''
    plots a various quantity as a function of radius based on argument
    Arg: 
        theGrid (Grid) : class that contains information on system
        which_q (list of strings string) : what quantity/s to plot
            possible key words include: 
                'Gas_Density', 'Gas_Pressure', 
    '''
    fig, ax = plt.subplots(len(which_Qs),1, figsize=(8,2*len(which_Qs)))
    for qi, q in enumerate(which_Qs):
        if q == 'Gas_Density':
            PlotGasDensity(theGrid, ax[qi])
        elif q == 'Gas_Pressure':
            PlotGasPressure(theGrid, ax[qi])
        elif q == 'Drift_Velocity':
            PlotDriftVelocity(theGrid, ax[qi])
        elif q == 'Dust_Density':
            PlotDustDensity(theGrid, ax[qi])
        
    for axi in ax:
        axi.legend()
        axi.grid()
        axi.set_xlabel('radius [AU]')
        axi.set_yscale('log')
    plt.show()


def PlotSpherical2D_DustImage(theGrid, cutfrac=0, colormap='plasma'):
    '''
    Plots a 2D Image extrapolated from the 1D calculation for dust image assuming spherical symmetry.
    Args:
        theGrid (Grid) : the input Grid() model
        cutfrac (float): the fraction of the domain that you wish to cut out of the image 
                            (e.g. cutfrac=0.4 cuts 40% of the domain off both ends)
        colormap (string): the desired color map scheme to be used in image

    '''
    N = len(theGrid.radius)
    xvals = np.arange(-N, N)
    yvals = np.arange(-N, N)
    X, Y = np.meshgrid(xvals, yvals)
    R = np.sqrt(X**2 + Y**2)
    # convert the radial array into 1D coordinates
    R_index = R.astype(int)
    # if any indices are over, force them to be one less
    R_index[R_index >= N] = N-1
    DustDensityImage = theGrid.dust_density()[R_index]
    rmax_AU = np.max(theGrid.radius)/Constants.AU
    #plt.figure(figsize=(6,6))
    cutfrac = .5 # cut 50% of the domain in both directions to zoom in on the dust accumulation
    cut_i = int(cutfrac*N)
    #print(DustDensityImage[cut_i:-cut_i, cut_i:-cut_i].shape)
    plt.imshow(DustDensityImage[cut_i:-cut_i, cut_i:-cut_i], 
            extent=[-rmax_AU*(1-cutfrac), rmax_AU*(1-cutfrac), -rmax_AU*(1-cutfrac), rmax_AU*(1-cutfrac)],
            cmap='plasma')
    plt.colorbar()
    plt.xlabel('radius [AU]')
    plt.ylabel('radius [AU]')
    plt.annotate('normalized dust density for toyGrid', xy=(.05, 1.02), xycoords='axes fraction')

