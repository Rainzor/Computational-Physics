from Schrage_16807 import*
import pandas as pd
from warnings import WarningMessage
import numpy as np
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod

#基类
class Growth_Model(ABC):
    """Growth models can be imagined as having a seed.
        Particles are then gradually added according to certain rules, 
        connecting them into clusters or clusters containing more particles.
        The growth model is a class of models that can be used to simulate the growth of clusters.
    
    Args:
        N (int, optional): The length of the side of the square. Defaults to 200.
        max_num (int, optional): The maximum number of particles. Defaults to 400.
        parameter (int, optional): The parameter of the model. Defaults to 1.

    Methods:
        _check_boundary: check if the point is in the boundary (TODO)
        run: run the model  (TODO)
        plot: plot the model (ACHIEVED)
        save: save the model data (ACHIEVED)
        data: return the model data (ACHIEVED)
        _check_point: check if the point is in the model (ACHIEVED)
        random_walk: random walk    (ACHIEVED)
    
    Attributes:
        N: the length of the image
        r: random number generator
        max_num: the max number of the particles
        _particles_set: the set of the particles
        parameter: the parameter of the model
        particle_number: the number of the particles
        model_name: the name of the model


    """
    def __init__(self, N=200, max_num=400,parameter=1):
        self.N = N  #图像边长
        self.r = Schrage16807(seed_time())
        if(max_num<N**2//4):            
            self.max_num = max_num #最大点数
        elif(max_num<N**2):
            raise WarningMessage("The max number of particals is too large, it will be set to N**2//4")
        self._particles_set = {(0, 0)}  # 集合
        self.parameter = parameter
        self.particle_number = 1#粒子数
        self.model_name = "Growth_Model"     
  
    @abstractmethod
    def _check_boundery(self, x, y):  
    # 检查点是否在边界内
        pass

    def _check_point(self, x, y):  
    # 检查点是否已在集合中
        return (x, y) in self._particles_set


    def random_walk(self,x,y):
    #随机漫步
        p = self.r.rand()
        if p < 1/4 :
            x -= 1
        elif p<1/2:
            x += 1
        elif p<3/4:
            y -= 1
        else:
            y += 1
        return x,y

    def add_particle(self,x,y):
    #生长一个粒子
        self._particles_set.add((x,y))

    @abstractmethod
    def update(self):
    #更新内部状态
        pass

    @abstractmethod
    def run(self):
    #运行程序
        pass

    @property
    def data(self):
    #返回数据
        x,y = zip(*self._particles_set)
        return x,y
    
    def plot(self, **kwargs):
    #画图
        size = kwargs['size']
        label = kwargs['label']
        x, y = self.data
        plt.style.use('dark_background')
        plt.scatter(x, y, s=size, marker=',', color='white',
                    edgecolors='none', label=label)
        plt.title(self.model_name)
        plt.legend()
        #隐藏坐标轴
        plt.axis('square')
        plt.xticks([])
        plt.yticks([])
        plt.show()

    def plot_data(self, filename, **kwargs):
    #读取数据
        df = pd.read_csv(filename)
        x = df['x']
        y = df['y']
        #画图
        size = kwargs['size']
        label = kwargs['label']
        plt.style.use('dark_background')
        plt.scatter(x, y, s=size, marker=',', color='white',
                    edgecolors='none', label=label)
        plt.title(filename)
        plt.legend()
        #隐藏坐标轴
        plt.axis('square')
        plt.xticks([])
        plt.yticks([])
        plt.show()
    def plot_GIF(self,**kwargs):
    #画动图
        """plot the GIF of the model

        """
        size = kwargs['size']
        label = kwargs['label']

        plt.style.use('dark_background')
        plt.figure(figsize=(10, 10))
        plt.ion()
        while self.particle_number < self.max_num:
            self.update()
            self.particle_number+=1

            plt.cla()
            plt.title(self.model_name+" created by Rainzor")
            x,y=zip(*self._particles_set)
            plt.scatter(x, y, s=size, marker=',', color='white',
                    edgecolors='none', label=label)
            plt.legend()
            #隐藏坐标轴
            plt.axis('square')
            plt.xticks([])
            plt.yticks([])
            plt.xlim(-self.N//2, self.N//2)
            plt.ylim(-self.N//2, self.N//2)

            plt.pause(0.01)
        plt.ioff()
        plt.show()


    def save(self):
    #保存数据
        x, y = self.data
        df = pd.DataFrame(data={'x': x, 'y': y})
        df.to_csv(self.model_name+'.csv', index=False)

    def read(self,filename):
    #读取数据
        df = pd.read_csv(filename)
        x = df['x']
        y = df['y']
        self._particles_set = set(zip(x,y))
        self.particle_number = len(self._particles_set)
        
        
#DLA模型
class DLA(Growth_Model):
    """DLA model
    The Diffusion Limited Aggregation (DLA) algorithm used a random walk to add cells.
    The random walk is biased by a potential field that attracts the walker to the cluster.

    For DLA, this has the stickiness parameter, which varies from 0 to 1

    We need to plot a picture to show the result，
    and create a 2D grid that is N by N cells in size

    Attributes:
        N: the length of the image

        parameter: the stickness parameter of the model

        r: random number generator

        model_name: DLA

        max_num: the max number of the particles

        radius: the max distance of the particles from the center

        _particles_set: the set of the particles

        particle_number: the number of the particles

    Methods:
        _check_boundery: check if the point is in the boundary 

        _check_point: check if the point is in the particles_set

        _check_around: check if the point is around the particles set

        _is_legal: check if the point is legal

        hausdorff_dim: calculate the hausdorff dimension of the model

        _sandbox: calculate the sandbox of the model

        _cov_density: calculate the coverage density of the model



    """
    def __init__(self,N,stickness=1,max_num=1000) -> None:
        self.N = N  #图像边长
        self.parameter = stickness #粘性参数
        if(max_num<N**2//4):            
            self.max_num = max_num #最大点数
        elif(max_num<N**2):
            raise WarningMessage("The max number of particals is too large, it will be set to N**2//4")

        self.model_name = "DLA"

        self.r = Schrage16807(seed_time())
        self.range_x = N//2+1
        self.range_y = N//2+1       
        self.radius = 20 #回转半径

        self.particle_number = 2      #记录点的个数
        self._particles_set = {(0,0),(1,0)} # store initial data of DLA

        self.dim = 2

    def _check_boundery(self,x,y):#检查点是否在边界内
        return -self.range_x< x <self.range_x and -self.range_y < y < self.range_y 

    # def _check_point(self,x,y):#检查点是否已在集合中
    #     return (x,y) in self._particles_set

    def _check_around(self,x,y):
        #检查周围8个点是否已被占据
        for i in range(-1,2):
            for j in range(-1,2):
                if (x+i, y+j) in self._particles_set:
                    return True
        return False
    def _check_around4(self,x,y):
        for i in [-1,1]:
            for j in [-1,1]:
                if (x, y+j) in self._particles_set:
                    return True
                if (x+i, y) in self._particles_set:
                    return True
        return False
    def _is_legal(self,x,y):
        return self._check_boundery(x,y) and  not self._check_point(x,y)

    @property
    def stickness(self):
        return self.parameter
#------------------------DLA数据生成--------------------------------
    def _rand_sample(self):
    #随机采样圆上的点或边界上的点
        radius = np.round(2*self.radius)
        if(radius<self.N//2):#如果直径小于N，就在圆上采样
            while 1:
                x = self.r.rand()*2-1
                y = self.r.rand()*2-1
                if x**2+y**2 <= 1:
                    break
            r = np.sqrt(x**2+y**2)
            x = np.round(x/r*radius)
            y = np.round(y/r*radius)
            return x,y
        else:#如果直径大于N，就在边界上采样

            N = self.N
            r_int = int(self.r.rand()*N-N//2) #向0取整
            p = self.r.rand()
            if p<1/4:
                x = N//2
                y = r_int
            elif p<1/2:
                x = -N//2 
                y = r_int
            elif p<3/4:
                x = r_int
                y = N//2
            else:
                x = r_int
                y = -N//2

            return x,y

    def generate_point(self):
    #创造一个随机点
        while 1:
            x, y = self._rand_sample()
            if self._is_legal(x,y):#不能在粒子集合内也不能超出边界
                break
        return  x, y

    def update(self):
        """
        添加一个粒子的步骤如下：

        1. 生成一个随机点,进入循环

        2. 如果随机游走的点超出边界，就重新生成一个随机点

        3. 随机漫步，直到碰到集合中的点

        4. 如果没有黏住，就继续随机游走下去

        5. 如果黏住，就将该点加入集合中,循环结束
        """
        x, y = self.generate_point()
        while 1:
            x, y = self.random_walk(x, y)
            if not self._check_boundery(x, y):  # 如果超出边界，就重新生成一个随机点
                x, y = self.generate_point()
            if self._check_around4(x, y):
                if self.r.rand() < self.stickness:  # 如果周围有点且黏住
                    self.add_particle(x, y)
                    if x**2+y**2 > self.radius**2:
                        self.radius = np.sqrt(x**2+y**2)#更新回转半径
                    return x,y

    def run(self):
        """Run the DLA model
        代码主体内容在update函数中
        """
        while self.particle_number < self.max_num:
            print("Particle number: ", self.particle_number)
            self.update()
            self.particle_number += 1

    def read(self,filename):
        #读取数据
        df = pd.read_csv(filename)
        x = df['x']
        y = df['y']
        r = x**2+y**2
        self.radius = np.sqrt(r.max())
        self._particles_set = set(zip(x,y))
        self.particle_number = len(self._particles_set)
#----------------------分形维数计算--------------------------------
    def hausdorff_dim(self, method):
        if(method=='sandbox'):
            return self._sandbox()
        elif(method == 'cov_density'):
            return self._cov_density()

    def _sandbox(self):
            #沙盒法计算维数
        radius = 4
        x,y = self.data
        x = np.array(x)
        y = np.array(y)
        num_list = []
        exp_list = []
        len_bound = self.radius*4/5
        while radius < len_bound:            
            #利用Boolean数组进行筛选
            num = sum((x<radius)&(x>-radius)&(y<radius)&(y>-radius))
            num_list.append(num)
            exp_list.append(np.log2(radius*2))
            radius = radius*np.sqrt(2)
        
        exp_list = np.array(exp_list)
        log_num_list = np.log2(num_list)
        coeff = np.polyfit(exp_list, log_num_list, 1)
        plt.plot(exp_list, log_num_list,'o',label='data points')
        plt.plot(exp_list, coeff[0]*exp_list+coeff[1],label='fitting line')
        plt.legend()
        plt.title('Sandbox method')
        plt.xlabel('log2(r)')
        plt.ylabel('log2(N)')
        plt.show()
        self.dim = coeff[0]
        return coeff[0]

    def _cov_density(self):
        #密度相关法计算维数
        set_data = self._particles_set
        num = len(set_data)
        set_l = set()
        set_r = set()
        set_u = set()
        set_d = set()

        radius = 16
        density_list = []
        exp_list = []
        len_bound = self.radius*4/5
        while radius < len_bound:
            #每次循环增长2倍半径
            set_l.clear()
            set_r.clear()
            set_u.clear()
            set_d.clear()
            for x,y in set_data:
                set_l.add((x-radius,y))
                set_r.add((x+radius,y))
                set_u.add((x,y+radius))
                set_d.add((x,y-radius))

            #利用集合的交集计算密度
            density = (len(set_data&set_l)+len(set_data&set_r)+len(set_data&set_u)+len(set_data&set_d))
            density_list.append(density)
            exp_list.append(np.log2(radius))

            # if(radius*np.sqrt(2) < len_bound):
            #     break

            #这里计算的是 sqrt2 倍的半径，得到更加精细的结果
            set_l.clear()
            set_r.clear()
            set_u.clear()
            set_d.clear()
            for x,y in set_data:
                set_l.add((x-radius,y+radius))
                set_r.add((x+radius,y-radius))
                set_u.add((x+radius,y+radius))
                set_d.add((x-radius,y-radius))

            #利用集合的交集计算密度
            density = (len(set_data&set_l)+len(set_data&set_r)+len(set_data&set_u)+len(set_data&set_d))
            density_list.append(density)
            exp_list.append(np.log2(radius)+0.5)
            radius = radius<<1

        exp_list = np.array(exp_list[:-1])
        log_num_list = np.log2(density_list[:-1])
        density_list = np.array(density_list)
        coeff = np.polyfit(exp_list, log_num_list, 1)
        plt.plot(exp_list, log_num_list, 'o', label='data points')
        plt.plot(exp_list, coeff[0]*exp_list+coeff[1], label='fitting line')
        plt.legend()
        plt.title('density related method')
        plt.xlabel(r'$\log_2(r)$')
        plt.ylabel(r'$\log_2(C(r))$')
        plt.show()

        self.dim = 2 + coeff[0]
        return self.dim


if __name__=="__main__":
    dla = DLA(N=500,stickness=1,max_num=10000)
    # 生成DLA数据 10000个点，时间大概需要5min
    # time_start = time.time()
    # dla.run()
    # time_end=time.time()
    # min = (time_end-time_start)//60
    # sec = (time_end-time_start)%60
    # print('DLA {} points time cost: {}min, {}sec'.format(dla.max_num,min,sec))
    # dla.plot(size=0.01,label='stickness={}'.format(dla.stickness))

    # 为了实验方便，直接读取数据，即可复现实验结果
    dla.read("DLA.csv")
    # 画图
    #dla.plot(size=0.02,label='stickness={}'.format(dla.stickness))
    
    #计算分形维数，两种方法，输出结果，并画图
    dim = dla.hausdorff_dim('sandbox')
    print("hausdorff dimension by sandbox: ", dim)

    dim = dla.hausdorff_dim('cov_density')
    print("hausdorff dimension by cov_density: ", dim)
