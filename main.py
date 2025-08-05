
import numpy as np
import Constants 

class Grid():
    '''
    class that represents the grid that gets initialized by the user with
    an input radius, gas density, and gas temperature profile
    '''
    def __init__(self, radius, sigma_gas, Tgas, mu_gas, Mstar, grain_size, unit_type):
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
        self.grain = grain_size
        self.rho_g = sigma_g /(np.sqrt(2*np.pi)*H)

        # sound speed
        Cs = np.sqrt(Constants.k_B * Tgas * (radius/Constants.AU))**(-0.5) / (mu_gas*Constants.m_H) # cm/s
        # keplerian angular velocity
        Omega_K = np.sqrt(Constants.G * Mstar / radius**3)
        # keplerian velocity
        self.v_K = Omega_K * radius
        H = Cs / Omega_K
        
        # pressure
        self.Pressure =rho_g * Cs **2
        self.St = (np.pi / 2) * rho_s * s / sigma_gas
        self.dpdr = np.gradient(self.Pressure, self.radius)
    def vdrift(self):
        eta = -0.5 *self.radius* self.dpdr / (self.rho_g * self.v_K**2) 
        v_drift= -2 * eta * v_K * St / (1 + St**2)
        return v_drift
   


def Initialize_System(radius, sigma_gas, Tgas, mu_gas, Mstar, grain_size,unit_type='cgs'):
    '''
    Args: 
        radius (array of floats) -- radius distribution
        sigma_gas (array of floats) -- surface gas density 
        Tgas (array of floats) -- gas temperature 
    Returns: theGrid (Grid) with user-input radius, sigma_gas, Tgas
    '''
    theGrid = Grid(radius, sigma_gas, Tgas, mu_gas, Mstar, grain_size,unit_type)
    return theGrid
