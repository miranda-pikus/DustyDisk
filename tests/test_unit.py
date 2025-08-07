import numpy as np
import os
from DustyDisk.Functions import *
from DustyDisk.Constants import AU as au
import pytest

def test_constant_pressure():
    '''
    checks the calculation of vdrift and density for dust given a constant pressure 
    '''
    r = np.linspace(0.1, 100, 1000) * au
    sigma_gas = 200 
    Tgas = 300 # Kelvin, assume constant temperature for now 
    mugas = 2.34  # mean molecular weight of gas
    Mstar = 3 * Constants.Msun
    constP_Grid = Initialize_System(r, sigma_gas, Tgas, mugas, Mstar, grain_size=1e-3,pressure_type='constant')
    # for constant pressure, we expect the drift velocity to be 0
    vdrift_vals = constP_Grid.vdrift()

    assert pytest.approx(vdrift_vals) == 0.0