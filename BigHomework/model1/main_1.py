import numpy as np
import matplotlib.pyplot as plt
import random
import time
import pandas as pd
from Schrage_16807 import *

r = Schrage16807(seed_time())
mu_B = 1.0

class Ising:
    def __init__(self, width=20, temperature=1.0, normal = (1,0,0)):
        self.L = width
        self.T = temperature
        self.J = 1.0
        N = width * width
        self.num_spins_xy = N
        self.num_spins_1D = width
        XY_list = r.rand(self.num_spins_xy)*2*np.pi
        Ising1D_list = np.empty(self.num_spins_1D)
        for i in range(width):
            Ising1D_list[i] = 1 if r.rand() > 0.5 else -1
        self.spin_config = {"XY": XY_list, "1D": Ising1D_list}

        nbr_XY = {i: ((i // width) * width + (i + 1) % width, (i + width) % N,
                    (i // width) * width + (i - 1) % width, (i - width) % N)
                    for i in list(range(N))}
        nbr_1D = {i: ((i + 1) % width, (i - 1) % width)
                    for i in list(range(width))}
        self.nbr = {"XY": nbr_XY, "1D": nbr_1D}

        nbr_1DtoXY = [(int)(N/2+i)%N for i in list(range(width))]
        nbr_XYto1D = {(int)(N/2+i)%N: i for i in list(range(width))}
        self.nbr_mix = {"XY": nbr_1DtoXY, "1D": nbr_XYto1D}

        self.normal = np.asarray(normal)
        self.normal2D = np.asarray([normal[0], normal[1]])
        self.energy = self.get_energy()
        self.magnetization  = self.get_magnetization()
        self.Cv = 0
        self.Cv_list = []
        self.T_list = []
        self.energy_list = []
        self.magnetization_list = []
    
    def init_state(self):
        XY_list = r.rand(self.num_spins_xy)*2*np.pi
        Ising1D_list = np.empty(self.num_spins_1D)
        for i in range(self.num_spins_1D):
            Ising1D_list[i] = 1 if r.rand() > 0.5 else -1
        self.spin_config = {"XY": XY_list, "1D": Ising1D_list}
        self.energy = self.get_energy()
        self.magnetization  = self.get_magnetization()
        self.energy_list = []
        self.magnetization_list = []
    def clear(self):
        self.energy_list = []
        self.magnetization_list = []
        self.Cv_list = []
        self.T_list = []
    def set_temperature(self, temperature):
        self.T = temperature

    def set_normal(self, normal):
        self.normal = normal

    def get_local_energy_XY(self, k, new_spin=None):
        if new_spin != None:
            spin = new_spin
        else:
            spin = self.spin_config["XY"][k]

        if k in self.nbr_mix["1D"]:
            index_1D = self.nbr_mix["1D"][k]
        else:
            index_1D = None

        direc = np.array([np.cos(spin), np.sin(spin)])
        sum_XY = -sum(np.cos(spin-self.spin_config["XY"][n])*self.J for n in self.nbr["XY"][k])
        if index_1D != None:
            return sum_XY - np.dot(direc, self.normal2D)*mu_B*self.spin_config["1D"][index_1D]
        else:
            return sum_XY
        
    def get_local_energy_1D(self, k, new_spin=None):
        if new_spin != None:
            spin = new_spin
        else:
            spin = self.spin_config['1D'][k]
        index_XY = self.nbr_mix["XY"][k]
        direc = np.array([np.cos(self.spin_config["XY"][index_XY]), np.sin(self.spin_config["XY"][index_XY])])
        sum_1D = -sum(spin*self.spin_config["1D"][n]*self.J for n in self.nbr["1D"][k])
        return sum_1D - np.dot(direc, self.normal2D)*mu_B*spin

    def get_energy(self):
        energy_ = 0
        idx = 0
        for idx in range(self.num_spins_xy):  # calculate energy per spin
            energy_ += self.get_local_energy_XY(idx)
        for idx in range(self.L):
            energy_ += self.get_local_energy_1D(idx)
        return energy_/2
    
    def get_magnetization(self):
        mag_1D = np.sum(self.spin_config["1D"])*self.normal
        return np.asarray([ np.sum(np.cos(self.spin_config["XY"]))+mag_1D[0],
                            np.sum(np.sin(self.spin_config["XY"]))+mag_1D[1],
                            mag_1D[2]])*mu_B


    def sweep_XY(self, idx):
        beta = 1.0 / self.T
        old_energy = self.get_local_energy_XY(idx)
        new_theta = self.spin_config["XY"][idx] + r.uniform(-np.pi/4, np.pi/4)
        new_energy = self.get_local_energy_XY(idx, new_theta)
        delta_E = new_energy - old_energy
        if r.rand() < np.exp(-beta * delta_E):
            self.spin_config["XY"][idx] = new_theta % (2*np.pi)

    def sweep_1D(self, idx):
        beta = 1.0 / self.T
        #k = np.random.randint(0, N - 1)#randomly choose a spin
        if(r.rand() < 0.5):
            old_energy = self.get_local_energy_1D(idx)
            new_spin = -self.spin_config["1D"][idx]
            new_energy = self.get_local_energy_1D(idx, new_spin=new_spin)
            delta_E = new_energy - old_energy
            if r.rand() < np.exp(-beta * delta_E):
                self.spin_config["1D"][idx] = new_spin

    def sweeps(self):
        spin_XY_idx = list(range(self.num_spins_xy))
        random.shuffle(spin_XY_idx)
        num_XY = 0
        spin_1D_idx = list(range(self.L))
        random.shuffle(spin_1D_idx)
        num_1D = 0

        if(r.rand() < 0.5) and (num_XY < self.num_spins_xy):
            self.sweep_XY(spin_XY_idx[num_XY])
            num_XY += 1
        elif num_1D < self.num_spins_1D:
            self.sweep_1D(spin_1D_idx[num_1D])
            num_1D += 1

    def equilibrate(self, max_nsweeps=int(4e3), temperature=None, show=False):
        if temperature != None:
            self.T = temperature
        try:
            if(temperature <= 0):
                raise ValueError
        except ValueError:
            print("temperature must be positive!")
            self.T = 0.01
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
            try:
                if(energy < -1e10):
                    raise ValueError
            except ValueError:
                self.init_state()
                print("\nenergy too large, reset state")
                k = 0
                continue
            else:
                magnetization = self.get_magnetization()  # 磁化强度
                dic_thermal_t['energy'].append(energy)
                self.energy_list.append(energy)
                dic_thermal_t['magnetization'].append(magnetization)
                self.magnetization_list.append(np.linalg.norm(magnetization))
                #print( abs(energy-energy_temp)/abs(energy))
                if show & (k % 1e3 == 0):
                    print('#sweeps=%i' % (k+1))
                    print('energy=%.2f' % energy)
                    print('magnetization=%.2f' % np.linalg.norm(magnetization))
                    self.show()
                error = abs(energy-energy_temp)
                #if ((error < 1e-6) & (k > 500)) or k == max_nsweeps-1:
                if k == max_nsweeps-1:
                    print('equilibrium state is reached at T=%.3f' % self.T)
                    print('#sweep=%i' % (k+1))
                    print('energy=%.8f' % energy)
                    print('magnetization=%.8f' % np.linalg.norm(magnetization))
                    print('error=%.10f' % error)
                    break
                energy_temp = energy
        nstates = len(dic_thermal_t['energy'])
        if(nstates < 4e4):
            nstates = int(nstates*2/3)
        elif(nstates <1e5):
            nstates = int(-2e4)
        else:
            nstates = int(-5e4)
        # 舍弃前面的热化过程
        energy = np.average(dic_thermal_t['energy'][int(nstates):])
        self.energy = energy
        energy2 = np.average(
            np.power(dic_thermal_t['energy'][int(nstates):], 2))
        self.Cv_list.append((energy2-energy**2)*(beta**2))  # 计算热容
        self.Cv = self.Cv_list[-1]
        print('Cv=%.8f' % self.Cv)
        self.T_list.append(self.T)
        magnetization = np.average(
            dic_thermal_t['magnetization'][int(nstates):], axis=0)
        self.magnetization = np.asarray(magnetization)



    #模拟退火
    def annealing(self, T_init=2.5, T_final=0.1, n=20, max_nsweeps=int(1e4),show_equi=False,**kwargs):
        # initialize spins. Orientations are taken from 0 - 2pi randomly.
        #initialize spin configuration
        dic_thermal = {}
        dic_thermal['temperature'] = list(np.linspace(T_init, T_final, n))
        dic_thermal['energy'] = []
        dic_thermal['Cv'] = []
        dic_thermal['magnetization'] = []
        
        for T in dic_thermal['temperature']:
            start_time =  time.time()
            if(T == T_final):
                self.equilibrate(max_nsweeps=int(max_nsweeps+5e4), temperature=T)
            else:
                self.equilibrate(max_nsweeps=max_nsweeps, temperature=T)
            #self.equilibrate(max_nsweeps=max_nsweeps, temperature=T)
            if show_equi:
                self.show()
            dic_thermal['energy'] += [self.energy]
            dic_thermal['Cv'] += [self.Cv]
            dic_thermal['magnetization'] += [self.magnetization]
            delta_time = time.time() - start_time
            print('equilibration time cost:%d m %.3f s' % (delta_time//60, delta_time%60))
        # plt.plot(dic_thermal['temperature'], dic_thermal['Cv'], '.-')
        # plt.ylabel(r'$C_v$')
        # plt.xlabel('T')
        # plt.title('Heat capacity in XY model with 1D Ising interaction')
        # plt.show()
        # plt.plot(dic_thermal['temperature'], dic_thermal['energy'], '.-')
        # plt.ylabel(r'$\langle E \rangle$')
        # plt.xlabel('T')
        # plt.title('Energy in XY model with 1D Ising interaction')
        # plt.show()
        # mag_inten = np.linalg.norm(dic_thermal['magnetization'], axis=1)
        # plt.plot(dic_thermal['temperature'],mag_inten, '.-')
        # plt.ylabel(r'$\langle M \rangle$')
        # plt.xlabel('T')
        # plt.title('Magnetization in XY model with 1D Ising interaction')
        # plt.show()

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
    def show(self, colored=False, save=False):
        config_matrix = self.list2matrix(self.spin_config["XY"])
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
        if save:
            plt.savefig('./img/T%.3f_model1.png' % self.T)
            plt.clf()
        else:
            plt.show()

    def plot_Cv(self):
        plt.plot(self.T_list, self.Cv_list, '*-')
        plt.xlabel('Temperature')
        plt.ylabel('Specific heat')
        plt.title('L=%i' % self.L)
        plt.show()

    def plot_energy(self,start=0,save=False):
        plt.plot(self.energy_list[start:])
        plt.xlabel('Steps')
        plt.ylabel('Energy')
        plt.title('T=%.3f' % self.T)
        if save:
            plt.savefig('./img/energy_T%.3f_model1.png' % self.T)       
            plt.clf()
        else:
            plt.show()
    
    def plot_magnetization(self,save=False):
        plt.plot(self.magnetization_list)
        plt.xlabel('Steps')
        plt.ylabel('Magnetization')
        plt.title('T=%.3f' % self.T)
        if save:
            plt.savefig('./img/magnetization_T%.3f_model1.png' % self.T)   
            plt.clf()
        else:
            plt.show()

if __name__ == '__main__':
    model = Ising(width=15, temperature=0.2)
    dic_thermal = {}
    dic_thermal["temperature"] = []
    dic_thermal["energy"] = []
    dic_thermal["Cv"] = []
    dic_thermal["magnetization"] = []

    T_list = np.asarray(
        [0.3,0.1,0.01])
    start_time = time.time()
    for T in T_list:
        model.equilibrate(max_nsweeps=int(3e5),temperature=T)
        #model.annealing(T_init=1.35, T_final=1.3, n=2, max_nsweeps=int(2e5),show_equi=False)
        delta_time = time.time()-start_time
        print('time cost: %d min %.2f sec ' % (delta_time//60, delta_time % 60))
        #model.show()
        model.plot_energy(save=True)
        #model.plot_energy(start=int(2.5e5))
        model.plot_magnetization(save=True)
        model.show(save=True)
        dic_thermal["temperature"].append(model.T)
        dic_thermal["energy"].append(model.energy)
        dic_thermal["Cv"].append(model.Cv)
        m = model.magnetization
        dic_thermal["magnetization"].append(np.linalg.norm(m))
        model.init_state()
    dic_thermal = pd.DataFrame(dic_thermal)
    dic_thermal.to_csv('data_backup.csv', index=False)


    # #target_T = np.arange(1.3, 0, -0.1)
    # target_T = np.asarray([0.1,0.4,0.8,1.4])
    # #target_T = np.arange(0.1, 1.3, 0.1)

    # start_time = time.time()
    # for t in target_T:
    #     model.annealing(T_init=t+0.05, T_final=t, n=2, max_nsweeps=int(8e4),show_equi=False)
    #     #model.equilibrate(max_nsweeps=int(10e4),temperature=t)
    #     dic_thermal["temperature"].append(model.T)
    #     dic_thermal["energy"].append(model.energy)
    #     dic_thermal["Cv"].append(model.Cv)
    #     dic_thermal["magnetization"].append(model.magnetization)
    #     model.plot_energy(save=True)
    #     model.plot_magnetization(save=True)
    #     model.init_state()
    
    # #dic_thermal = model.annealing(T_init=3.5, T_final=1, n=26, max_nsweeps=int(5e4),show_equi=False)
    # delta_time = time.time()-start_time
    # min = delta_time//60
    # sec = delta_time % 60
    # hour = min//60
    # min = min % 60
    # print('time cost:%d hours %d mins %.2f secs ' % (hour, min, sec))
    # dic_thermal["magnetization"] = np.linalg.norm(
    #     dic_thermal['magnetization'], axis=1)
    # dic_thermal = pd.DataFrame(dic_thermal)
    # dic_thermal.to_csv('data_backup.csv',index=False)
    # df = pd.read_csv('data.csv')
    # df = pd.concat([df,dic_thermal],axis=0)
    # df.sort_values(by='temperature',inplace=True)
    # df.index = range(len(df))
    # df = df.round({"temperature": 3})
    # df = df.groupby("temperature").mean()
    # df.reset_index(inplace=True)



    # t = df["temperature"].values
    # e = df["energy"].values
    # c = df["Cv"].values
    # m = df["magnetization"].values

    # fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    # ax[0, 0].scatter(t, e, s=0.5)
    # ax[0, 0].set_xlabel("temperature")
    # ax[0, 0].set_ylabel("energy")
    # ax[0, 0].set_title("energy")
    # ax[0, 1].scatter(t, c, s=2, color="red")
    # ax[0, 1].plot(t, c, linewidth=1, color="black")
    # ax[0, 1].set_ylabel("Cv")
    # ax[0, 1].set_title("Cv")
    # ax[1, 0].plot(t, m, linewidth=1)
    # ax[1, 0].set_xlabel("temperature")
    # ax[1, 0].set_ylabel("magnetization")
    # ax[1, 0].set_title("magnetization")
    # plt.show()

    # df.to_csv('data.csv',index=False)

