from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.constants import e, m_p


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

        self.plane = plane
        self.dmu = dmu
        self.M = np.dot(B_inv, np.dot(R, B))
        self.C = np.dot(B_inv, np.dot(I, B))
        self.S = np.dot(B_inv, np.dot(S, B))

    def kick(self, beam):

        if self.plane=='x':
            beam.x, beam.xp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * beam.x \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) * beam.xp, \
                              (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * beam.x \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu)) * beam.xp
        if self.plane=='y':
            beam.y, beam.yp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * beam.y \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) * beam.yp, \
                            + (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * beam.y \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu)) * beam.yp
        if self.plane=='z':
            beam.z, beam.dp = (self.C[0,0]*np.cos(self.dmu) + self.S[0,0]*np.sin(self.dmu)) * beam.z \
                            + (self.C[0,1]*np.cos(self.dmu) + self.S[0,1]*np.sin(self.dmu)) * beam.dp, \
                            + (self.C[1,0]*np.cos(self.dmu) + self.S[1,0]*np.sin(self.dmu)) * beam.z \
                            + (self.C[1,1]*np.cos(self.dmu) + self.S[1,1]*np.sin(self.dmu)) * beam.dp


class RFMap(MachineElement):

    def __init__(self, V, h, dphi):
        pass
