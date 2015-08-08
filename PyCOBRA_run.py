import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, e, m_p

from PyCOBRA_base import Beam, TwissMap


n_macroparticles = 1000
intensity = 1e11
energy = 26e9
gamma = np.sqrt(1 + (energy/(m_p*c**2))**2)
beam = Beam(n_macroparticles, weight=intensity/n_macroparticles,
            gamma=gamma, epsn_x=2e-6, epsn_y=2e-6, epsn_z=0.3)

twiss_x = TwissMap()


