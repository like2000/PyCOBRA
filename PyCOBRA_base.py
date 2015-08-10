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
        return emittance_normalised(self.x, self.xp)

    def epsn_y(self):
        return emittance_normalised(self.y, self.yp)

    def epsn_z(self):
        return emittance_normalised(self.z, self.dp)


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


class MachineElement(object):

    __metaclass__ = ABCMeta

    @abstractmethod
    def kick(self, beam):
        pass


class TwissMap(MachineElement):

    def __init__(self, plane='x',
                 alpha_0=0, beta_0=100, alpha_1=0, beta_1=100, dmu=0,
                 *detuners):

        B     = np.array([[1./np.sqrt(beta_0), 0],
                          [alpha_0/np.sqrt(beta_0), np.sqrt(beta_0)]])
        R     = np.array([[np.cos(dmu),  np.sin(dmu)],
                          [-np.sin(dmu), np.cos(dmu)]])
        B_inv = np.array([[np.sqrt(beta_1), 0],
                          [-alpha_1/np.sqrt(beta_1), 1./np.sqrt(beta_1)]])

        I = np.array([[1, 0],
                      [0, 1]])
        S = np.array([[0, 1],
                      [-1, 0]])

        self.dmu = dmu
        self.M = np.dot(B_inv, np.dot(R, B))
        self.C = np.dot(B_inv, np.dot(I, B))
        self.S = np.dot(B_inv, np.dot(S, B))

    def kick(self, beam):

        if self.plane=='x':
            beam.x, beam.xp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * self.x \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) * self.xp, \
                              (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * self.x \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu)) * self.xp
        if self.plane=='y':
            beam.y, beam.yp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * self.y \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) *self.yp, \
                            + (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * self.y \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu))* self.yp
        if self.plane=='z':
            beam.z, beam.dp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * self.z \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) *self.dp, \
                            + (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * self.z \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu))* self.dp


class RFMap(MachineElement):

    def __init__(self, V, h, dphi):
        pass
