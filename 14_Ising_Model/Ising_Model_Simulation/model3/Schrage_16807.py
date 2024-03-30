import time
import numpy as np
import pandas as pd


def seed_time():  # 取计算机时间为种子值
    t = time.localtime(time.time())
    return t.tm_year-2000 + 70*(t.tm_mon+12*(t.tm_mday+31*(t.tm_hour+23*(t.tm_min+59*t.tm_sec))))


class Schrage16807(object):
    def __init__(self, seed=42):
        self.__seed = seed
        self.__M = 2147483647
        self.__a = 16807

    def __iter_rand(self, i_n):  # Schrage 随机数迭代方法
        a = self.__a
        m = self.__M
        q = m//a
        r = m % a
        i_next = a * (i_n % q) - r * (i_n // q)
        if i_next < 0:
            return i_next + m
        else:
            return i_next

    # 返回一个随机数或者一个随机向量
    def __rand_list(self, dim=1) -> tuple[np.ndarray, float]:
        if dim == 1:  # 返回一个随机数
            self.__seed = self.__iter_rand(self.__seed)
            return self.__seed/self.__M
        elif dim > 1:  # 返回dim维的行向量随机数
            i_n = self.__seed
            l = np.zeros(dim)
            for i in range(dim):
                i_n = self.__iter_rand(i_n)
                l[i] = i_n
            self.__seed = i_n
            return l/self.__M

    def rand(self, d0=1, d1=None):
        if isinstance(d0, int):  # 保证是一个整数
            if d0 == 1:
                if d1 is None:
                    return self.__rand_list(1)  # 返回一个随机数
                elif isinstance(d1, int) and (d1 > 0):
                    temp = self.__rand_list(d1)
                    return temp  # 返回随机数或行向量
            elif d0 > 1:
                if d1 is None:
                    return self.__rand_list(d0)  # 返回行向量
                elif isinstance(d1, int) and (d1 > 0):
                    if d1 == 1:
                        temp = self.__rand_list(d0)
                        return temp.T  # 返回列向量
                    else:
                        temp = [self.__rand_list(d1) for item in range(d0)]
                        a = np.array(temp)
                        return a.reshape((d0, d1))  # 返回d0*d1的随机数矩阵
        return print("Error, please input a positive integer!\n")

    def __normal(self, size):
        """生成服从正态分布的随机数

        Parameters
        ----------
        loc : float
            正态分布的均值
        scale : float
            正态分布的标准差
        size : int or tuple of ints
            随机数的个数
        """
        i = 0
        tmp = np.zeros(size).flatten()
        l = len(tmp)

        while i < l:
            u = self.__rand_list(1)*2-1
            v = self.__rand_list(1)*2-1
            r = np.sqrt(u**2+v**2)
            if r >= 1:
                continue
            else:
                tmp[i] = (u/r)*np.sqrt(-4*np.log(r))
                if i+1 < l:
                    tmp[i+1] = (v/r)*np.sqrt(-4*np.log(r))
                i += 2
                # print((u/r)*np.sqrt(-4*np.log(r)), (v/r)*np.sqrt(-4*np.log(r)))


    def normal(self, loc=0, scale=1,size=None, ):
        if size is None:
            return self.__normal(1)*scale+loc
        else:
            return self.__normal(size)*scale+loc

    def uniform(self, low=0, high=1,size=None):
        if size is None:
            return self.__rand_list(1)*(high-low)+low
        else:
            return self.__rand_list(size)*(high-low)+low

    def randint(self, low=0, high=1,size=None):
        if size is None:
            return np.floor(self.__rand_list(1)*(high-low)+low)
        else:
            return np.floor(self.__rand_list(size)*(high-low)+low)


def makefile(data, filename):
    df = pd.DataFrame(data=data)
    df.to_csv(filename+'.csv', index=False)


def readfile(filename):
    df = pd.read_csv(filename)
    return np.array(df)
