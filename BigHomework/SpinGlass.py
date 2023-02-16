import numpy as np
import matplotlib.pyplot as plt
import random
from Schrage_16807 import *

r = Schrage16807(seed_time())


#Simulated Tempering to simulate Spin Glass in Ising 2D model
class SpinGlass:
    def __init__(self, width=20, low_temperature=1.0, high_temperature=3.0):
        self.L = width
        self.T_l = low_temperature
        self.T_h = high_temperature
        N = width * width
        self.J = np.zeros((N,N))
        for i in range(N):
            for j in range(i,N):
                if i == j:
                    self.J[i,j] = 0
                else:
                    self.J[i,j] = r.normal(0,1)
                    self.J[j,i] = self.J[i,j]
        self.num_spins = N
        self.spin_config = np.empty(self.num_spins)
        for i in range(self.num_spins):
            self.spin_config[i] = 1 if r.rand() > 0.5 else -1
        self.nbr = {i : ((i // width) * width + (i + 1) % width, (i + width) % N,
                    (i // width) * width + (i - 1) % width, (i - width) % N) \
                                            for i in list(range(N))}
        self.average_energy = np.sum(self.get_energy())/self.num_spins
        self.Cv = []
        self.M = []

    def set_temperature(self,temperature):
        self.T = temperature;

    def get_energy(self):
        energy_ = np.zeros(np.shape(self.spin_config))
        idx = 0
        for spin in self.spin_config:
            energy_[idx] = -sum(spin*self.spin_config[n]*self.J[idx,n]
                                for n in self.nbr[idx])
            idx += 1
        return energy_/2
