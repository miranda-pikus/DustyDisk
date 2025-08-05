import numpy as np

class Grid():
    '''
    class that represents the grid that gets initialized by the user with
    an input radius, gas density, and gas temperature profile
    '''
    def __init__(self, radius, sigma_gas, Tgas, unit_type):
        self.radius = radius # array of radius
        self.sigma_gas = sigma_gas # array of gas density
        self.Tgas = Tgas # array of gas temperature 
        self.unit_type = unit_type # unit system (string)
        self.Pressure = np.zeros(len(radius))
        # integrate different units here eventually, for now assume cgs

    #def CalculatePressureGrid(self):




def Initialize_System(radius, sigma_gas, Tgas, unit_type='cgs'):
    '''
    Args: 
        radius (array of floats) -- radius distribution
        sigma_gas (array of floats) -- surface gas density 
        Tgas (array of floats) -- gas temperature 
    Returns: theGrid (Grid) with user-input radius, sigma_gas, Tgas
    '''
    if len(radius)*len(sigma_gas)*len(Tgas)/3 != len(radius):
        print('cannot initialize system, input arrays are different lengths! try again.')
        return 0.
    else: 
        theGrid = Grid(radius, sigma_gas, Tgas, unit_type)
        return theGrid
