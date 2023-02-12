import numpy as np
import matplotlib.pyplot as plt
import random
from Schrage_16807 import *

r = Schrage16807(seed_time())

class XYModel:
    def __init__(self, width=20, temperature=1.0):
        self.L = width
        self.T = temperature
        self.J = 1.0
        N = width * width
        self.num_spins = N
        self.spin_config =self.spin_config = r.rand(self.num_spins)*2*np.pi
        self.nbr = {i : ((i // width) * width + (i + 1) % width, (i + width) % N,
                    (i // width) * width + (i - 1) % width, (i - width) % N) \
                                            for i in list(range(N))}
        self.energy = np.sum(self.get_energy())/self.num_spins
        self.Cv = []

    def set_temperature(self,temperature):
        self.T = temperature;

    def get_energy(self):
        energy_ = np.zeros(np.shape(self.spin_config))
        idx = 0
        for spin in self.spin_config:  # calculate energy per spin
            energy_[idx] = -sum(np.cos(spin-self.spin_config[n])
                                for n in self.nbr[idx])  # nearst neighbor of kth spin
            idx += 1
        return energy_/2

    def _mc_step(self):
        beta = 1.0 / self.T
        spin_idx = list(range(self.num_spins))
        random.shuffle(spin_idx)
        for idx in spin_idx:
            #k = np.random.randint(0, N - 1)#randomly choose a spin
            dtheta = r.uniform(-np.pi,np.pi)
            new_theta = self.spin_config[idx] + dtheta
            delta_E = -sum(np.cos(new_theta-self.spin_config[n])
                            - np.cos(self.spin_config[idx]-self.spin_config[n])
                            for n in self.nbr[idx]) 
            if r.uniform(0.0, 1.0) < np.exp(-beta * delta_E):
                self.spin_config[idx] += dtheta


    def equilibrate(self,max_nsweeps=int(1e4),temperature=None,H=None,show = False):
        if temperature != None:
            self.T = temperature
        dic_thermal_t = {}
        dic_thermal_t['energy']=[]
        beta = 1.0/self.T
        energy_temp = 0
        for k in list(range(max_nsweeps)):
            self._mc_step()
            #list_M.append(np.abs(np.sum(S)/N))
            energy = np.sum(self.get_energy())/self.num_spins #每个粒子平均能量
            dic_thermal_t['energy'] += [energy]
            #print( abs(energy-energy_temp)/abs(energy))
            if show  & (k%1e3 ==0):
                print('#sweeps=%i'% (k+1))
                print('energy=%.2f'%energy)
                self.show()
            if ((abs(energy-energy_temp)/abs(energy)<1e-4) & (k>500)) or k == max_nsweeps-1:
                print('\nequilibrium state is reached at T=%.1f'%self.T)
                print('#sweep=%i'%k)
                print('energy=%.2f'%energy)
                break
            energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        energy=np.average(dic_thermal_t['energy'][int(nstates/2):])#取后一半的平均值,舍弃前一半的热化过程
        self.energy = energy
        energy2=np.average(np.power(dic_thermal_t['energy'][int(nstates/2):],2))
        self.Cv=(energy2-energy**2)*beta**2#计算热容

    def simulate(self,max_nsweeps=int(1e4) ):
        for step in range(max_nsweeps):
            self._mc_step()
            #self.energies.append(self.energy())

    #模拟退火
    def annealing(self,T_init=2.5,T_final=0.1,n = 20,show_equi=False):
        # initialize spins. Orientations are taken from 0 - 2pi randomly.
        #initialize spin configuration  
        dic_thermal = {}
        dic_thermal['temperature']=list(np.linspace(T_init,T_final,n))
        dic_thermal['energy']=[]
        dic_thermal['Cv']=[]
        for T in dic_thermal['temperature']:
            self.equilibrate(temperature=T)
            if show_equi:
                self.show()
            dic_thermal['energy'] += [self.energy]
            dic_thermal['Cv'] += [self.Cv]
        plt.plot(dic_thermal['temperature'],dic_thermal['Cv'],'.')
        plt.ylabel(r'$C_v$')
        plt.xlabel('T')
        plt.show()
        plt.plot(dic_thermal['temperature'],dic_thermal['energy'],'.')
        plt.ylabel(r'$\langle E \rangle$')
        plt.xlabel('T')
        plt.show()
        return dic_thermal

    @staticmethod
    ## convert configuration inz list to matrix form
    def list2matrix(S):
        N = int(np.size(S))
        L = int(np.sqrt(N))
        S = np.reshape(S, (L, L))
        return S

    ## visulize a configurtion
    #  input：S/ spin configuration in list form
    def show(self, colored=False):
        config_matrix = self.list2matrix(self.spin_config)
        X, Y = np.meshgrid(np.arange(0, self.L), np.arange(0, self.L))
        U = np.cos(config_matrix)
        V = np.sin(config_matrix)
        plt.figure(figsize=(4, 4), dpi=100)
        #plt.title('Arrows scale with plot width, not view')
        Q = plt.quiver(X, Y, U, V, units='width')
        qk = plt.quiverkey(Q, 0.1, 0.1, 1, r'$spin$', labelpos='E',
                           coordinates='figure')
        plt.title('T=%.2f' % self.T+', #spins=' +
                  str(self.L)+'x'+str(self.L))
        plt.axis('off')
        plt.show()


    def plot_energy(self):
        plt.plot(self.energies)
        plt.xlabel('Monte Carlo steps')
        plt.ylabel('Energy')
        plt.show()


if __name__ == '__main__':
    model = XYModel(width=20, temperature=1)
    model.equilibrate()
    model.show()
