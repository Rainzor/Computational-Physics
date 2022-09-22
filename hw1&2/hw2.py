#Using UTF-8 Code Python
#Runze Wang
#PB20020480

from operator import ge
from queue import Queue 
import numpy as np
from hw1 import Schrage_16807

def decision(list):
    if list[0]>list[2] and list[2]>list[1]:
        return 1
    else:
        return 0


class LaggedFibonacci(object):
    def __init__(self,seed=42):
        self.__seed = seed
        self.__M = int(2**32-5)
        self.__numberQueue= self.__make_number43()
    
    def __make_number43(self):
        s = Schrage_16807(self.__seed);
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
        i_0 = self.__numberQueue[0].get()
        i_21 = self.__numberQueue[1].get()
        i_43 = i_21 - i_0 - 0;
        if i_43<0:
            i_43 = i_43 + self.__M - 1
        self.__numberQueue[0].put(i_21)
        self.__numberQueue[1].put(i_43)
        return i_43/(2147483647) #i_43/(2**32-1)


if __name__ == "__main__":
    N = 50002
    schrage_16807 = Schrage_16807(seed=42)
    fibonacci = LaggedFibonacci(seed=42)

    data_16807 = [schrage_16807.rand() for i in range(N)]
    data_fib = [fibonacci.rand() for i in range(N)]
    N_16807 = 0
    N_fib = 0

    for i in range(1,N-1):
        N_16807 += decision(data_16807[i-1:i+2])
        N_fib += decision(data_fib[i-1:i+2])
    
    print("Seed = %d, Schrage 16807: %.6f, Fibonacci: %.6f" %( 42, N_16807/(N-2), N_fib/(N-2)))




