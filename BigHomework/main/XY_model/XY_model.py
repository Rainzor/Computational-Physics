import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from Schrage_16807 import *

class XYModel:
    def __init__(self, width=20, temperature=1.0, H=np.zeros(2)):
        self.L = width
        self.T = temperature
        self.J = 1.0
        self.mu_B = 1.0
        self.r = Schrage16807(seed_time())
        self.H = H#外磁场
        N = width * width
        self.num_spins = N
        self.spin_config = self.r.rand(self.num_spins)*2*np.pi
        self.nbr = {i : ((i // width) * width + (i + 1) % width, (i + width) % N,
                    (i // width) * width + (i - 1) % width, (i - width) % N) \
                                            for i in list(range(N))}
        self.energy = self.get_energy()
        self.magnetization = self.get_magnetization()
        self.Cv = 0
        self.Cv_list = []
        self.T_list = []
        self.energy_list = []

    def set_temperature(self,temperature):
        self.T = temperature;
    def set_H(self,H):
        self.H = H
    def set_normal(self,normal):
        self.normal = normal
    def get_local_energy(self, k, new_spin = None):
        if new_spin != None:
            spin = new_spin
        else:
            spin = self.spin_config[k]
        direc = np.array([np.cos(spin),np.sin(spin)])
        return -sum(np.cos(spin-self.spin_config[n])*self.J + np.dot(direc,self.H)*self.mu_B
                                for n in self.nbr[k]) 

    def get_energy(self):
        energy_ = 0
        idx = 0
        for idx in range(self.num_spins):  # calculate energy per spin
            energy_ += self.get_local_energy(idx)
        return energy_/2
    
    def init_state(self):
            self.spin_config = self.r.rand(self.num_spins)*2*np.pi
            self.energy = self.get_energy()
            self.magnetization  = self.get_magnetization()
            self.energy_list = []

    def get_magnetization(self):
        #单位体积内的磁矩为磁化强度
        return np.asarray([np.sum(np.cos(self.spin_config)),np.sum(np.sin(self.spin_config))])*self.mu_B/self.num_spins

    def sweeps(self):
        beta = 1.0 / self.T
        spin_idx = list(range(self.num_spins))
        random.shuffle(spin_idx)
        for idx in spin_idx:
            #k = np.random.randint(0, N - 1)#randomly choose a spin
            drec = np.array([np.cos(self.spin_config[idx]),np.sin(self.spin_config[idx])])
            old_energy = self.get_local_energy(idx)
            new_theta = self.spin_config[idx] + self.r.uniform(-np.pi/4,np.pi/4)
            new_energy = self.get_local_energy(idx,new_theta)
            delta_E = new_energy - old_energy
            if self.r.rand() < np.exp(-beta * delta_E):
                self.spin_config[idx] = new_theta%(2*np.pi)


    def equilibrate(self,max_nsweeps=int(1e4),temperature=None,H = np.zeros(2), show = False):
        if temperature != None:
            self.T = temperature
        if H.any():
            self.H = H
        dic_thermal_t = {}
        dic_thermal_t['energy']=[]
        dic_thermal_t['magnetization']=[]
        beta = 1.0/self.T
        energy_temp = 0
        for k in list(range(max_nsweeps)):
            self.sweeps()
            #list_M.append(np.abs(np.sum(S)/N))
            energy = self.get_energy()#总能量
            magnetization = self.get_magnetization()#磁化强度
            dic_thermal_t['energy'].append(energy)
            dic_thermal_t['magnetization'].append(magnetization)
            #print( abs(energy-energy_temp)/abs(energy))
            if show  & (k%1e3 ==0):
                print('#sweeps=%i'% (k+1))
                print('energy=%.2f'%energy)
                print('magnetization=%.2f'%np.linalg.norm(magnetization))
                self.show()
            if ((abs(energy-energy_temp)/abs(energy)<1e-5) & (k>1000)) or k == max_nsweeps-1:
                print('\nequilibrium state is reached at T=%.2f'%self.T)
                print('#sweep=%i'%k)
                print('energy=%.3f'%energy)
                print('magnetization=%.3f'%np.linalg.norm(magnetization))
                break
            energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        energy=np.average(dic_thermal_t['energy'][int(nstates/2):])#取后一半的平均值,舍弃前一半的热化过程
        self.energy = energy
        energy2=np.average(np.power(dic_thermal_t['energy'][int(nstates/2):],2))
        self.Cv_list.append((energy2-energy**2)*(beta**2))#计算热容
        self.Cv = self.Cv_list[-1]
        self.T_list.append(self.T)
        magnetization = np.average(dic_thermal_t['magnetization'][int(nstates/2):],axis=0)
        self.magnetization = np.asarray(magnetization)

    def simulate(self, max_nsweeps=int(1e4), temperature=None, H = np.zeros(2)):
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
            if(k%5e3 == 0):
                print("loop : %d" % (k+1))
            self.sweeps()
            #list_M.append(np.abs(np.sum(S)/N))
            energy = self.get_energy()  # 总能量
            magnetization = self.get_magnetization()  # 磁化强度
            dic_thermal_t['energy'].append(energy)
            self.energy_list.append(energy)
            dic_thermal_t['magnetization'].append(magnetization)
            if ((abs(energy-energy_temp)/abs(energy) < 1e-5) & (k > 1000)) or k == max_nsweeps-1:
                print('#XY model sweeps=%i' % (k+1))
                break
            energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        # 取后一半的平均值,舍弃前一半的热化过程
        energy = np.average(dic_thermal_t['energy'][int(nstates/2):])
        self.energy = energy
        energy2 = np.average(
            np.power(dic_thermal_t['energy'][int(nstates/2):], 2))
        magnetization = np.average(
            dic_thermal_t['magnetization'][int(nstates/2):], axis=0)
        self.magnetization = np.asarray(magnetization)
        self.Cv = (energy2-energy**2)*(beta**2)#计算热容

    #模拟退火
    def annealing(self,T_init=2.5,T_final=0.1,n = 20,show_equi=False,max_nsweeps=int(1e4)):
        # initialize spins. Orientations are taken from 0 - 2pi randomly.
        #initialize spin configuration  
        dic_thermal = {}
        dic_thermal['temperature']=list(np.linspace(T_init,T_final,n))
        dic_thermal['energy']=[]
        dic_thermal['Cv']=[]
        for T in dic_thermal['temperature']:
            self.equilibrate(max_nsweeps=max_nsweeps,temperature=T)
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
        plt.title('T=%.3f' % self.T+', #spins=' +
                  str(self.L)+'x'+str(self.L))
        plt.axis('off')
        plt.show()


    def plot_Cv(self):
        plt.plot(self.T_list, self.Cv_list, '*-')
        plt.xlabel('Temperature')
        plt.ylabel('Specific heat')
        plt.show()
    def plot_energy(self):
        plt.plot(self.energy_list)
        plt.xlabel('Steps')
        plt.ylabel('Energy')
        plt.show()


if __name__ == '__main__':
    model = XYModel(width=25, temperature=0.01)
    dic_thermal = {}
    dic_thermal["temperature"] = []
    dic_thermal["energy"] = []
    dic_thermal["Cv"] = []

    model.simulate(max_nsweeps=int(2e4))
    model.show()
    # model.plot_energy()
    #model.annealing(T_init=30,T_final=0.8,n = 10,show_equi=False,max_nsweeps=int(5e3))
    #model.show()
    # T_list = np.linspace(0.1, 1, 7,dtype=float)
    # T_list = np.append(T_list, np.linspace(1.1, 1.5, 10,dtype=float))
    # T_list = np.append(T_list, np.linspace(1.6, 2.5, 8, dtype=float))
    # for t in T_list:
    #     model.equilibrate(max_nsweeps = int(2e4),temperature=t)
    #     dic_thermal["temperature"].append(t)
    #     dic_thermal["energy"].append(model.energy)
    #     dic_thermal["Cv"].append(model.Cv)
    #     model.init_state()
    # model.plot_Cv()

    # dic_thermal = pd.DataFrame(dic_thermal)
    # dic_thermal.to_csv('XY_model.csv',index=False)

