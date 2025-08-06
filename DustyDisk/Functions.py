
import numpy as np
import DustyDisk.Constants as Constants

class Grid():
    '''
    class that represents the grid that gets initialized by the user with
    an input radius, gas density, and gas temperature profile
    '''
    def __init__(self, radius, sigma_gas, Tgas, mu_gas, Mstar, grain_size, unit_type, dt=1e8, t_final =1e12):
        '''
        Initializes a Grid object

        Args: 
            radius (array) : radial profile
            sigma_gas (array) : gas density
            Tgas (array) : gas temperature
            mu_gas (float) : mean molecular weight of gas
            Mstar (float) : mass of central star
            unit_type (string) : unit system (cgs everywhere for now)

        Attributes: 
            radius 
            sigma_gas
            Tgas
            dt
            t_fin
            Nt
            grain_size
            Cs : sound speed
            v_K : keplerian velocity
            rho_g
            Pressure
            St
            dpdr : pressure gradient 
        '''
        self.radius = radius 
        self.sigma_gas = sigma_gas 
        self.Tgas = Tgas 
        self.dt = dt
        self.t_fin = t_final
        self.Nt =  int(t_final // dt)
        #self.unit_type = unit_type # unit system (string)
        self.grain_size = grain_size
        # sound speed
        self.Cs = np.sqrt(Constants.k_B * Tgas * (radius/Constants.AU)**(-0.5) / (mu_gas*Constants.m_H)) # cm/s
        # keplerian angular velocity
        Omega_K = np.sqrt(Constants.G * Mstar / radius**3)
        # keplerian velocity
        self.v_K = Omega_K * radius
        H = self.Cs / Omega_K

        self.rho_g = self.sigma_gas /(np.sqrt(2*np.pi)*H)

        
        # pressure
        self.Pressure =self.rho_g * self.Cs **2
        self.St = (np.pi / 2) * Constants.rho_s * self.grain_size / sigma_gas
        self.dpdr = np.gradient(self.Pressure, self.radius)

        
    def vdrift(self):
        eta = -0.5 *self.radius* self.dpdr / (self.rho_g * self.v_K**2) 
        v_drift= -2 * eta * self.v_K * self.St / (1 + self.St**2)

        return v_drift
    def dust_density(self):
        '''
        Args:
            dt (float): timestep (s)
            Nt (int): number of steps

        Returns:
            sigma_dust (array): normalized evolved dust surface density
        '''
        r = self.radius
        sigma_d = np.ones_like(self.sigma_gas)  
        
        v_drift = self.vdrift()
        for _ in range(self.Nt):
            F = sigma_d * v_drift
            dsigma_dt = -1 / r * np.gradient(r * F, r)
            sigma_d += dsigma_dt * self.dt
            sigma_d = np.maximum(sigma_d, 1e-20)

        sigma_d /= np.max(sigma_d)
        return sigma_d

    



def Initialize_System(radius: np.array, 
                      sigma_gas: np.array, 
                      Tgas, mu_gas, Mstar, grain_size,unit_type='cgs'):
    '''
    Initialize the environment to be handled by DustyDisk package.
    Depends on the user-defined input such as density, temperature, grain size, etc...

    Args: 
        radius (array): numpy vector. Radius values in the domain.
        sigma_gas (array) -- surface gas density 
        Tgas (array) -- gas temperature 

    Returns: 
        Grid: object that uses the density and temperature profiles for other calculations
    '''
    theGrid = Grid(radius, sigma_gas, Tgas, mu_gas, Mstar, grain_size,unit_type,dt=1e9, t_final=1e12)
    return theGrid
