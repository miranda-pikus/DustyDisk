import matplotlib.pyplot as plt
from Main import Grid
import Constants
plt.style.use('PlotStyling.mplstyle')

def PlotGasDensity(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.sigma_gas, label='gas density')
    #which_ax.set_ylabel('')

def PlotGasPressure(whichGrid, which_ax):
    which_ax.plot(whichGrid.radius/Constants.AU, whichGrid.Pressure, label='gas pressure')

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
            PlotGasDensity(theGrid, ax.ravel()[qi])
        elif q == 'Gas_Pressure':
            PlotGasPressure(theGrid, ax.ravel()[qi])
        
    for axi in ax.ravel():
        axi.grid()
        axi.set_xlabel('radius [AU]')
        axi.set_yscale('log')
    plt.show()

