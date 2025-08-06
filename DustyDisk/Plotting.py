import matplotlib.pyplot as plt
from Main import Grid
import Constants
plt.style.use('PlotStyling.mplstyle')

def PlotGasDensity(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.sigma_gas, 
                  color='blue', label='gas density')
    which_ax.set_ylabel(r'density (g cm$^{-3}$)')

def PlotGasPressure(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.Pressure, 
                  color='forestgreen', label='gas pressure')
    which_ax.set_ylabel(r'pressure (dyne)')

def PlotDriftVelocity(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.v_drift, 
                  color='purple', label='drift velocity')
    which_ax.set_ylabel(r'v$_{drift}$ (cm/s)')

def PlotDustDensity(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.rho_dust, 
                    color='black', label='dust density')
    which_ax.set_ylabel(r'density (g cm$^{-3}$)')


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


def Spherical2D_DustImage(theGrid):

    '''
    for when dust density is added
    rmax = 1
    N = 150
    rvals = np.linspace(0, rmax, N)
    height = 1
    width = .01*height
    pos = rmax*0.6
    rhodust = height * np.exp(-(rvals-pos)**2/(2*width**2))
    x = np.linspace(0, rmax, N)
    y = np.linspace(0, rmax, N)
    X,Y= np.meshgrid(x,y,indexing='ij')
    R = np.sqrt(X**2 + Y**2)
    RhoDust = height * np.exp(-(R-pos)**2/(2*width**2))
    plt.plot(rvals, rhodust)
    plt.figure()
    plt.imshow(RhoDust, extent=[0,rmax,0,rmax], norm='log',origin='lower',cmap='Blues')
    plt.colorbar()
    '''
