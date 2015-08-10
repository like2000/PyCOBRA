import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, e, m_p

import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
plt.ion()

from PyCOBRA_accelerator import TwissMap
from PyCOBRA_beam import Bunch, gaussian_generator


R = 6911./(2*np.pi)
energy = 26e9
gamma = np.sqrt(1 + (e*energy/(m_p*c**2))**2)
beta = np.sqrt(1-gamma**-2)
gamma_tr = 18
eta = gamma_tr**-2 - gamma**-2
Q_x = 20.13
Q_y = 20.18
Q_s = 0.017
beta_x = R/Q_x
beta_y = R/Q_y
beta_z = eta*R/Q_s
twiss_x = TwissMap('x', 0, beta_x, 0, beta_x, Q_x)
twiss_y = TwissMap('y', 0, beta_y, 0, beta_y, Q_y)
twiss_z = TwissMap('z', 0, beta_z, 0, beta_z, Q_s)


n_macroparticles = 2000
intensity = 1e11
bunch = Bunch(n_macroparticles,
              intensity/n_macroparticles, e, m_p, gamma,
              gaussian_generator(2e-6/(beta*gamma), ('x', 'xp'), 0, beta_x),
              gaussian_generator(2e-6/(beta*gamma), ('y', 'yp'), 0, beta_y),
              gaussian_generator(2.00, ('z', 'dp'), 0, beta_z))


fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12, 14), tight_layout=True)


one_turn_map = [twiss_x, twiss_y, twiss_z]
for k in range(60):
    for m in one_turn_map:
        m.kick(bunch)

    ax1.scatter(bunch.x, bunch.xp, c='g', lw=0.1)
    ax2.scatter(bunch.y, bunch.yp, c='g', lw=0.1)
    ax3.scatter(bunch.z, bunch.dp, c='g', lw=0.1)

    ax1.plot(np.mean(bunch.x), np.mean(bunch.xp), 'ro')
    ax2.plot(np.mean(bunch.y), np.mean(bunch.yp), 'ro')
    ax3.plot(np.mean(bunch.z), np.mean(bunch.dp), 'ro')

    [ax.set_xlim(-2e-2, 2e-2) for ax in [ax1, ax2]]
    [ax.set_ylim(-2e-4, 2e-4) for ax in [ax1, ax2]]
    ax3.set_xlim(-100, 100)
    ax3.set_ylim(-1, 1)

    plt.draw()
    [ax.cla() for ax in [ax1, ax2, ax3]]


plt.close()
