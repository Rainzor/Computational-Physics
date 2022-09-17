#Using UTF-8 Code Python
#Runze Wang
#PB20020480


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#Schrahe 16807类
#随机输入初始值，反复调用rand(),得到随机数列
class Schrage(object):
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



def uniformity(data, k = 10):#计算一致性
    lp = np.linspace(0,1,k+1)
    nk = [0] * k
    size = len(data)
    for i in range(size):
        j = (data[i]*k)%k
        nk[j] = nk[j] + 1
        # for j in range(k):
        #     if data[i] < lp[j+1] :
        #         continue
        #     else:
        #         nk[j] = nk[j]+1
        #         break
    mk = size/k
    chi = 0
    for i in range(k):
        chi = chi + ((nk[i]-mk)**2)/mk
    return chi

def covariance2(data, l=2):#计算相关性
    x = np.array(data[:-l])
    y = np.array(data[l:])
    d = np.array(data)

    rcl = (x*y).mean()
    r1 = (d.mean())**2 
    r2 = (d*d).mean()

    return (rcl - r1)/(r2 - r1)

def makefile(data,filename):
    df = pd.DataFrame({'data':data,})
    df.to_excel(filename+'xlsx',sheet_name='sheet1', index=False)


if __name__ == "__main__":
    test = Schrage(seed = 42)
    
    data = [test.rand() for i in range(100000)]
    print(data[:90]) 
    x = data[:-1]
    y = data[1:]

    #绘图，散点图
    plt.scatter(x,y, s=0.01, c='b')
    plt.show()

