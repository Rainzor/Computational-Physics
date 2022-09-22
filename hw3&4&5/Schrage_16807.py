import time
import numpy as np
def seed_time():
    t = time.localtime(time.time())
    return t.tm_year-2000 + 70*(t.tm_mon+12*(t.tm_mday+31*(t.tm_hour+23*(t.tm_min+59*t.tm_sec))))

class Schrage_16807(object):
    def __init__(self, seed=42):
        self.__seed = seed
        self.__M = 2147483647
        self.__a = 16807

    def __iter_rand(self, i_n):  # 随机数迭代
        a = self.__a
        m = self.__M
        q = m//a
        r = m % a
        i_next = a * (i_n % q) - r * (i_n // q)
        if i_next < 0:
            return i_next + m
        else:
            return i_next

    def rand(self,dim =1):
        if dim == 1:
            self.__seed = self.__iter_rand(self.__seed)
            return self.__seed/self.__M
        elif dim > 1:
            i_n =self.__seed
            l = np.zeros(dim)
            for i in range(dim):
                i_n=self.__iter_rand(i_n)
                l[i] = i_n;
            self.__seed = i_n
            return l/self.__M
         