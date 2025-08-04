import numpy as np
from scipy.signal import find_peaks

class PressureProfile(object):
    """
    1D array for pressure 
    either load in or analytical form of some kind? 
    """
    def __init__(self):
        self.pressure_arr = np.array([])
        self.radius_arr = np.array([])

def Find_Maxima(pressure_prof):
    """
    Args: pressure_prof (PressureProfile) 
    Returns: local maxima ([float]): list of floats
    """
    maxima, _ = find_peaks(pressure_prof.pressure_arr, height=0)
    rmaxima = pressure_prof.radius_arr[maxima]
    return rmaxima

    
