import numpy as np
import Constants 

class Grid():
    '''
    class that represents the grid that gets initialized by the user with
    an input radius, gas density, and gas temperature profile
    '''
    def __init__(self, radius, sigma_gas, Tgas, mu_gas, Mstar, unit_type):
        '''
        Args: 
            radius (array) : radial profile
            sigma_gas (array) : gas density
            Tgas (array) : gas temperature
            mu_gas (float) : mean molecular weight of gas
            Mstar (float) : mass of central star
            unit_type (string) : unit system (cgs everywhere for now)
        '''
        self.radius = radius 
        self.sigma_gas = sigma_gas 
        self.Tgas = Tgas 
        #self.unit_type = unit_type # unit system (string)

        # sound speed
        Cs = np.sqrt(Constants.k_B * Tgas * (radius/Constants.AU))**(-0.5) / (mu_gas*Constants.m_H) # cm/s
        # keplerian angular velocity
        Omega_K = np.sqrt(Constants.G * Mstar / radius**3)
        # keplerian velocity
        v_K = Omega_K * radius
        # pressure
        self.Pressure = sigma_gas * Cs * Omega_K / np.sqrt(2*np.pi)


def Initialize_System(radius, sigma_gas, Tgas, mu_gas, Mstar, unit_type='cgs'):
    '''
    Args: 
        radius (array of floats) -- radius distribution
        sigma_gas (array of floats) -- surface gas density 
        Tgas (array of floats) -- gas temperature 
    Returns: theGrid (Grid) with user-input radius, sigma_gas, Tgas
    '''
    theGrid = Grid(radius, sigma_gas, Tgas, mu_gas, Mstar, unit_type)
    return theGrid
