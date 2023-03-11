#Using UTF-8 Code Python
#Runze Wang
#PB20020480

import time
import statistics
from scipy.stats import chi2 #调用卡方分布
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Schrahe16807类
#随机输入初始值，反复调用rand(),得到随机数列
class Schrage16807(object):
    def __init__(self,seed = 42):
        self.__seed = seed
        self.__M = 2147483647
        self.__a = 16807

    def __iter_rand(self, i_n):#随机数迭代
        a = self.__a
        m = self.__M
        q = m//a
        r = m%a
        i_next = a * (i_n % q) - r * (i_n // q)
        if i_next < 0:
            return i_next + m
        else:
            return i_next

    def rand(self):
        self.__seed = self.__iter_rand(self.__seed)
        return self.__seed/self.__M

#取计算机时间为种子值
def seed_time():
    t = time.localtime(time.time())
    return t.tm_year-2000 + 70*(t.tm_mon+12*(t.tm_mday+31*(t.tm_hour+23*(t.tm_min+59*t.tm_sec))))


#计算随机数一致性，得到Chi2
def uniformity(data, k = 10):
    lp = np.linspace(0,1,k+1)
    nk = [0] * k
    size = len(data)
    for i in range(size):
        j =int((data[i]*k)//1)
        nk[j] = nk[j] + 1
        
    mk = size/k
    c = 0
    average = 0
    for i in range(k):
        c = c + ((nk[i]-mk)**2)/mk
        average = average 
    return c

#计算随机数相关性
def covariance2(data, l=2):
    x = np.array(data[:-l])
    y = np.array(data[l:])
    d = np.array(data)

    rcl = (x*y).mean()
    r1 = (d.mean())**2 
    r2 = (d*d).mean()

    return (rcl - r1)/(r2 - r1)

def makefile(data,filename):
    df = pd.DataFrame({'data':data,})
    df.to_csv(filename+'.csv',index=False)

def readfile(filename):
    df = pd.read_csv(filename)
    return np.array(df['data'])

if __name__ == "__main__":
    s = seed_time()
    test = Schrage16807(seed = s)

    #data = readfile(filename='random_data_1e7.csv') #读入原始数据 
    data = [test.rand() for i in range(10000001)] 
    data = np.array(data)

    # #绘图，散点图
    # x = data[0:5000]
    # y = data[2:5002]
    # plt.scatter(x,y, s=1, c='b')
    # plt.title("Plane distribution of random numbers(interval l=2 )")
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()

    #用k阶矩检验均匀性
    print("k阶矩检验均匀性:\n")
    error = []
    for k in range(6):
        for i in range(1,8):        
            xk = data[:int(10**i)]**k
            print("N:%.2e, k:%d, Average of xk:%.5f, Expectation:%.5f" %(10**i, k, xk.mean(), 1/(k+1)))
            if k==5:
                error.append(abs(xk.mean()-1/(k+1)))
        print("\n")
    for i in range(7):
         print("N:%.2e, k:5, error:%.5f" %(10**i,error[i]))
    print("\n")

    # #绘图
    # x = np.linspace(1,7,1001)
    # y = error[0]*np.sqrt(10)/(np.sqrt(10**x))
    # plt.plot(range(1,8),error,'o-',color = 'g',label="error")#o-:圆形
    # plt.plot(x,y,'--',color='b',label=r'$O(\frac{1}{\sqrt{N}})$')
    # plt.xlabel("lg(N)")
    # plt.ylabel("error")
    # plt.legend()
    # plt.show()



    #检验均匀性(卡方分布)
    print("卡方分布检验均匀性:\n")
    alpha = 0.05 #显著水平a(或置信度1-a)
    for i in range(4,8):
        for k in range(2,11):
            c2= uniformity(data[:int(10**i)],k)
            percent_point = chi2.ppf(1-alpha,k-1) 
            print("N:%.2e, k:%d, Statistics:%.5f, Percent_point:%.5f"%(10**i,k,c2,percent_point))
            
        print("\n")

    #检验关联性
    print("检验关联性:\n")
    for i in range(1,10):
        covar2 = covariance2(data=data,l=i)    
        print("l = %d, C(l) = %.10f"%(i,covar2))



