from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.constants import e, m_p


class Beam(object):

    def __init__(self, n_macroparticles,
                 weight=1, charge=e, mass=m_p, gamma=1,
                 epsn_x=1, epsn_y=1, epsn_z=1):

        self.n_macroparticles = n_macroparticles

        self.weight = weight
        self.charge = charge
        self.mass   = mass
        self.gamma  = gamma

        self.x  = np.zeros(n_macroparticles)
        self.xp = np.zeros(n_macroparticles)
        self.y  = np.zeros(n_macroparticles)
        self.yp = np.zeros(n_macroparticles)
        self.z  = np.zeros(n_macroparticles)
        self.dp = np.zeros(n_macroparticles)


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
            beam.x = np.dot(self.M, beam.x)
            beam.xp = np.dot(self.M, beam.xp)
        if self.plane=='y':
            beam.y = np.dot(self.M, beam.y)
            beam.yp = np.dot(self.M, beam.yp)
        if self.plane=='z':
            beam.z = np.dot(self.M, beam.z)
            beam.dp = np.dot(self.M, beam.dp)


class RFMap(MachineElement):

    def __init__(self, V, h, dphi):
        pass
