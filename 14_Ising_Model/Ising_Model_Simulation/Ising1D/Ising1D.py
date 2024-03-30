import numpy as np
import matplotlib.pyplot as plt
import random
from Schrage_16807 import *


class Ising1D:
    def __init__(self, width=20, temperature=1.0, H=np.zeros(3), normal = (0,0,1)):
        self.T = temperature
        self.J = 1.0
        self.H = H  # 外磁场
        self.normal = np.asarray(normal)
        self.mu_B = 1.0
        self.r = Schrage16807(seed_time())
        N = width
        self.num_spins = N
        self.spin_config = np.empty(self.num_spins)
        for i in range(self.num_spins):
            self.spin_config[i] = 1 if self.r.rand() > 0.5 else -1
        self.nbr = {i: ((i + 1) % N, (i - 1) % N)
                    for i in list(range(N))}
        self.energy = self.get_energy()
        self.Cv = 0
        self.magnetization = self.get_magnetization()
        self.Cv_list = []
        self.T_list = []

    def init_state(self):
        self.spin_config = np.empty(self.num_spins)
        for i in range(self.num_spins):
            self.spin_config[i] = 1 if self.r.rand() > 0.5 else -1
        self.energy = self.get_energy()
        self.magnetization = self.get_magnetization()

    def set_temperature(self, temperature):
        self.T = temperature

    def set_H(self, H):
        self.H = H

    def get_local_energy(self, k,new_spin=None):
        if new_spin != None:
            spin = new_spin
        else:
            spin = self.spin_config[k]
        return -sum(spin*self.spin_config[n]*self.J + spin*np.dot(self.H,self.normal)*self.mu_B
                        for n in self.nbr[k])

    def get_energy(self):
        energy_ = 0
        idx = 0
        for idx in range(self.num_spins):  # calculate energy per spin
            energy_ += self.get_local_energy(idx)
        return energy_/2


    def get_magnetization(self):
        #单位体积内的磁矩为磁化强度
        return sum(self.spin_config)*self.mu_B/self.num_spins 

    def sweeps(self):
        beta = 1.0 / self.T
        spin_idx = list(range(self.num_spins))
        random.shuffle(spin_idx)
        for idx in spin_idx:
            #k = np.random.randint(0, N - 1)#randomly choose a spin
            if(self.r.rand() < 0.5):
                old_energy = self.get_local_energy(idx)
                new_spin = -self.spin_config[idx]
                new_energy = self.get_local_energy(idx, new_spin=new_spin)
                delta_E = new_energy - old_energy
                if self.r.rand() < np.exp(-beta * delta_E):
                    self.spin_config[idx] = new_spin

    def equilibrate(self, max_nsweeps=int(1e4), temperature=None, H=np.zeros(3), show=False):
        if temperature != None:
            self.T = temperature
        if H.any():
            self.H = H
        dic_thermal_t = {}
        dic_thermal_t['energy'] = []
        dic_thermal_t['magnetization'] = []
        beta = 1.0/self.T
        energy_temp = 0
        print("\nequilibrate: target temperature = %.3f" % self.T)

        for k in list(range(max_nsweeps)):
            if(k%5e3 == 0):
                print("loop : %d" % (k+1))
            self.sweeps()
            #list_M.append(np.abs(np.sum(S)/N))
            energy = self.get_energy()  # 总能量
            magnetization = self.get_magnetization()  # 磁化强度
            dic_thermal_t['energy'].append(energy)
            dic_thermal_t['magnetization'].append(magnetization)
            #print( abs(energy-energy_temp)/abs(energy))
            if show & (k % 1e3 == 0):
                print('#sweeps=%i' % (k+1))
                print('energy=%.2f' % energy)
                print('magnetization=%.2f' % np.linalg.norm(magnetization))
                self.show()
            if k == max_nsweeps-1:
                print('\nequilibrium state is reached at T=%.1f' % self.T)
                print('#sweep=%i' % k)
                print('energy=%.2f' % energy)
                print('magnetization=%.2f' % np.linalg.norm(magnetization))
                break
            energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        if(nstates < 4e4):
            nstates = int(nstates*2/3)
        elif(nstates <=1e5):
            nstates = int(-2e4)
        else:
            nstates = int(-5e4)
        # 取后一半的平均值,舍弃前一半的热化过程
        energy = np.average(dic_thermal_t['energy'][int(nstates):])
        self.energy = energy
        energy2 = np.average(
            np.power(dic_thermal_t['energy'][int(nstates):], 2))
        self.Cv_list.append((energy2-energy**2)*(beta**2))  # 计算热容
        print("Cv = %.3f" % self.Cv_list[-1])
        self.T_list.append(self.T)
        magnetization = np.average(
            dic_thermal_t['magnetization'][int(nstates):], axis=0)
        self.magnetization = np.asarray(magnetization)

    def simulate(self, max_nsweeps=int(1e4), temperature=None, H=np.zeros(3)):
        if temperature != None:
                self.T = temperature
        if H.any():
            self.H = H
        dic_thermal_t = {}
        dic_thermal_t['energy'] = []
        dic_thermal_t['magnetization'] = []
        beta = 1.0/self.T
        energy_temp = 0
        for k in list(range(max_nsweeps)):
            self.sweeps()
            #list_M.append(np.abs(np.sum(S)/N))
            energy = self.get_energy()  # 总能量
            magnetization = self.get_magnetization()  # 磁化强度
            dic_thermal_t['energy'].append(energy)
            dic_thermal_t['magnetization'].append(magnetization)
            if ((abs(energy-energy_temp)/abs(energy) < 1e-8) & (k > 500)) or k == max_nsweeps-1:
                print('#Ising model sweeps=%i' % (k+1))
                break
            energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        if(nstates < 4e4):
            nstates = int(nstates*2/3)
        elif(nstates <=1e5):
            nstates = int(-2e4)
        else:
            nstates = int(-5e4)
        # 取后一半的平均值,舍弃前一半的热化过程
        energy = np.average(dic_thermal_t['energy'][int(nstates/2):])
        self.energy = energy
        energy2 = np.average(
            np.power(dic_thermal_t['energy'][int(nstates/2):], 2))
        magnetization = np.average(
            dic_thermal_t['magnetization'][int(nstates/2):], axis=0)
        self.magnetization = np.asarray(magnetization)
        self.Cv = (energy2-energy**2)*(beta**2)  # 计算热容

    #模拟退火
    def annealing(self, T_init=2.5, T_final=0.1, n=20, show_equi=False):
        # initialize spins. Orientations are taken from 0 - 2pi randomly.
        #initialize spin configuration
        dic_thermal = {}
        dic_thermal['temperature'] = list(np.linspace(T_init, T_final, n))
        dic_thermal['energy'] = []
        dic_thermal['Cv'] = []
        for T in dic_thermal['temperature']:
            self.equilibrate(temperature=T)
            if show_equi:
                self.show()
            dic_thermal['energy'] += [self.energy]
            dic_thermal['Cv'] += [self.Cv]
        plt.plot(dic_thermal['temperature'], dic_thermal['Cv'], '*-')
        plt.ylabel(r'$C_v$')
        plt.xlabel('T')
        plt.show()
        plt.plot(dic_thermal['temperature'], dic_thermal['energy'], '*-')
        plt.ylabel(r'$\langle E \rangle$')
        plt.xlabel('T')
        plt.show()
        return dic_thermal

    ## visulize a configurtion
    #  input：S/ spin configuration in list form
    def show(self, colored=False):
        return

    def plot_Cv(self):
        plt.plot(self.T_list, self.Cv_list, '*-')
        plt.xlabel('Temperature')
        plt.ylabel('Specific heat')
        plt.title('Specific heat of 1D Ising model')
        plt.show()


if __name__ == '__main__':
    model = Ising1D(width=50, temperature=1, H=np.array([0, 0, 1]))
    #model.equilibrate()
    T_list = np.linspace(0.1, 1.5, 15)
    # T_list = np.append(T_list, np.linspace(0.5, 1, 20))
    # T_list = np.append(T_list, np.linspace(1, 3, 10))
    for t in T_list:
        model.equilibrate(max_nsweeps = int(1e5),temperature=t)
        model.init_state()
    model.plot_Cv()
    #model.annealing(T_init=2.5, T_final=0.1, n=20, show_equi=False)
