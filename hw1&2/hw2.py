import imp
from operator import ge
from queue import Queue 
import numpy as np
from hw1 import Schrage

def decision(list):
    if list[0]>list[2] and list[2]>list[1]:
        return True
    else:
        return False

class Fibonacci(object):
    def __init(self, seed = 42):
        self.__seed = seed
        self.__M = int(2**32-5)
        self.__number_queue= self.__make_number44(seed)
    
    def __make_number44(self,seed):
        s = Schrage(seed);
        q1 = Queue(maxsize=21)
        q2 = Queue(maxsize=22)
        #使用16807随机生成前43个数
        for j in range(21):
            q1.put(s.rand())
        for j in range(22):
            q2.put(s.rand())
        q = (q1,q2)
        return q
    
    def rand(self):
        i_0 = self.__number_queue[0].get()
        i_21 = self.__number_queue[1].get()
        i_43 = i_21 - i_0 - 0;
        if temp<0:
            temp = temp + self.__M - 1
        self.__number_queue[0].put(i_21)
        self.__number_queue[1].put(temp)
        return temp

