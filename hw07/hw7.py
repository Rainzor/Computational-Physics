from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#选择函数F(x)
def F_x(x):
    if x >= 2900 and x <= 2994:
        return 0.015
    elif x < 3006:
        return 0.096
    elif x <= 3013:
        return 0.005
    else:
        return 0

#the inverse function for cdf of F(x)
def F_ppf(p,x):
    if x<0 or x>1:
        print("Error!\n")
        return
    if p<0 or p>1:
        print("Error!\n")
        return
    elif p<=0.5429:
        return 94*x+2900
    elif p<0.9865:
        return 12*x+2994
    else:
        return 7*x+3006



class MyDistribution(object):
    """
    根据实验数据得到的离散分布

    Attributes:
        p_data:概率密度分布数据
        c_data:累积分布函数数据
        size:数据点个数

    Method:
        pdf(x):概率密度函数
        cdf(x):累积分布函数
        ppf(p):分位点
        rvs_direct(size):直接抽样法得到的随机变量
        rvs_select(size):舍选法得到的随机变量
    """
    def __init__(self,data) -> None:
        data = np.array(data)
        self.p_data = data[np.argsort(data[:, 0])]
        self.size = np.size(self.p_data,0)
        self.c_data = data.copy()
        for i in range(self.size):
            if i>0:
                self.c_data[i, 1] = self.c_data[i, 1]+self.c_data[i-1, 1]
            else:
                continue
    #概率密度函数
    def pdf(self,x):
        """Probability density function.

            Parameters:
                x(float): a random variable.

            Returns:
                p_x(float): the probility density of x
        """
        data_x = self.p_data[:,0]
        x_min = data_x[0]
        x_max = data_x[-1]
        p_x = self.p_data[:,1]
        if(x < x_min or x > x_max+1):
            return 0
        else:
            x = np.floor(x)
            return p_x[data_x == x]
    #累积分布函数
    def cdf(self, x):
        """Cumulative density function.

            Parameters:
			x(float): a random variable.
			
		Returns:
		    c_x(float): the Cumulative density probility at x
        """
        data_x = self.c_data[:, 0]
        x_min = data_x[0]
        x_max = data_x[-1]
        c_x = self.c_data[:, 1]
        if x < x_min:
            return 0
        elif x > x_max+1:
            return 1
        else:
            x = np.floor(x)
            return c_x[data_x == x]
    
    #分位点 (累积分布函数的反函数).
    def ppf(self, p):
        """Percent point function (inverse of cdf — percentiles).

            Parameters:
			p(float): the probility in cumulative density function

            Returns:
		    x(int): Percent point 
        """
        data_x = self.c_data[:, 0]
        x_min = data_x[0]
        x_max = data_x[-1]
        c_x = self.c_data[:, 1]
        if p>0 or p<1:
            index = self.__binary_search(p,c_x)
            return data_x[index]
        elif p == 0:
            return x_min 
        elif p == 1:
            return x_max
        else:
            print("error")
            return 


    #直接法得到size大小的随机数
    def rvs_direct(self,size=1):
        """Random variates using direct sampling method

            Parameters:
            size(int): the size of random variates as np.array format 

            Return:
            x_1(array(int)):random variates
        """
        seed = seed_time()
        r = Schrage16807(seed)
        xi = r.rand(size)
        x_1 = np.array([self.ppf(p) for p in xi])
        return x_1


    #舍选法得到size大小的随机数
    def rvs_select(self,size=1):
        """Random variates using selection sampling method

            Parameters:
            size(int): the size of random variates as np.array format 

            Return:
            result(array(float)):random variates

            print the efficiency of selection    
        """
        seed = seed_time()
        r = Schrage16807(seed)
        i = 0
        j = 0
        result = np.zeros(size)
        while(i<size):
            xi_x = F_ppf(p = r.rand(),x = r.rand())
            xi_y = r.rand()
            if self.__selection_method(xi_x,xi_y):
                result[i] = xi_x
                i = i+1
            j = j+1
        print("Efficiency of Sampling:%.5f\n"%(size/j))
        return result

    #舍选法判断函数
    def __selection_method(self,x,y):
        M = F_x(x)
        return M*y<=self.pdf(x)

    #二分查找法
    def __binary_search(self,x,data):
        low = 0
        high = self.size-1
        while low <= high:
            mid = (low+high)//2
            if data[mid] == x:
                return mid
            elif data[mid] > x:
                high = mid-1
            else:
                low = mid+1
        if  low < self.size:
            return low
        else:
            return 0



if __name__=="__main__":
    df = pd.read_table('data.TXT')
    dataset = np.array(df,'float')
    #归一化处理
    sum_n = np.sum(dataset[:,1])
    dataset[:,1] = dataset[:,1]/sum_n
    N = 100000

    #输入数据，得到相关的概率分布函数类
    my_dis = MyDistribution(dataset)

    #获得随机数据
    r_d = my_dis.rvs_direct(N)
    r_s = my_dis.rvs_select(N)

    #绘制直方图与概率密度分布图像
    num_bins = 113


    #绘制直接抽样法得到的图像
    fig, ax = plt.subplots()

    #the histogram of the data
    n, bins, patches = ax.hist(r_d, num_bins, density=True, label='Direct sampling data')

    # add a 'best fit' line
    ax.plot(dataset[:,0], dataset[:,1], '--', label="Best fit")
    ax.set_xlabel('Energy(eV)')
    ax.set_ylabel('Probability density')

    # Tweak spacing to prevent clipping of ylabel
    ax.legend()
    fig.tight_layout()
    plt.figure(1)


    #绘制舍选法得到的图像
    fig, ax = plt.subplots()

    #the histogram of the data
    n, bins, patches = ax.hist(r_s, num_bins, density=True, label='Selection sampling data')

    # add a 'best fit' line
    ax.plot(dataset[:,0], dataset[:,1], '--', label="Best fit")
    ax.set_xlabel('Energy(eV)')
    ax.set_ylabel('Probability density')

    # Tweak spacing to prevent clipping of ylabel
    ax.legend()
    fig.tight_layout()
    plt.figure(2)

    plt.show()

    df_rvs = pd.DataFrame({"Direct_Sampling":r_d,
                            "Selection_Sampling":r_s})
    df_rvs.to_csv("data_07.csv")






    
    