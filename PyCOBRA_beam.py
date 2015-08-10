from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.constants import e, m_p


def gaussian_generator(eps_geo, phase_space_tuple=('x', 'xp'), alpha=0, beta=1):

    sigma = np.sqrt(eps_geo)

    def generate(bunch):
        n_macroparticles = bunch.n_macroparticles
        x = np.random.normal(scale=sigma, size=n_macroparticles)
        xp = np.random.normal(scale=sigma, size=n_macroparticles)

        M = np.array([[np.sqrt(beta), 0],
                      [-alpha/np.sqrt(beta), 1./np.sqrt(beta)]])
        x, xp = M[0,0]*x + M[0,1]*xp, M[1,0]*x + M[1,1]*xp

        setattr(bunch, phase_space_tuple[0], x)
        setattr(bunch, phase_space_tuple[1], xp)

    return generate


class Bunch(object):

    def __init__(self, n_macroparticles,
                 weight=1, charge=e, mass=m_p, gamma=1,
                 *phase_space_generators):

        self.n_macroparticles = n_macroparticles

        self.weight = weight
        self.charge = charge
        self.mass   = mass
        self.gamma  = gamma

        [generate(self) for generate in phase_space_generators]

    def emittance_normalised(self, x, xp):

        return np.sqrt(self.gamma**2 - 1) * \
            np.sqrt( np.std(x**2)*np.std(xp**2) - np.std(x*xp)**2 )

    def epsn_x(self):
        return self.emittance_normalised(self.x, self.xp)

    def epsn_y(self):
        return self.emittance_normalised(self.y, self.yp)

    def epsn_z(self):
        return self.emittance_normalised(self.z, self.dp)


class Beam(object):

    def __init__(self, bunches_list):

        self.n_macroparticles = sum([b.n_macroparticles for b in bunches_list])

        self.weight = np.concatenate(b.weight for b in bunches_list)
        self.charge = np.concatenate(b.charge for b in bunches_list)
        self.mass   = np.concatenate(b.mass for b in bunches_list)
        self.gamma  = np.concatenate(b.gamma for b in bunches_list)

        self.x  = np.concatenate(b.x for b in bunches_list)
        self.xp = np.concatenate(b.xp for b in bunches_list)
        self.y  = np.concatenate(b.y for b in bunches_list)
        self.yp = np.concatenate(b.yp for b in bunches_list)
        self.z  = np.concatenate(b.z for b in bunches_list)
        self.dp = np.concatenate(b.dp for b in bunches_list)
