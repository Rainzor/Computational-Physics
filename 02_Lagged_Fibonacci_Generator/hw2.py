#Using UTF-8 Code Python
#Runze Wang
#PB20020480

from queue import Queue 
import numpy as np
from hw17.Schrage_16807 import*

#Fibonacci延迟产生器
class LaggedFibonacci(object):
    def __init__(self,seed=42):
        self.__seed = seed
        self.__numberQueue= self.__make_number43()

    def __iter_rand(self, i_n):#随机数迭代
        a = 16807
        m = 2147483647
        q = m//a
        r = m%a
        i_next = a * (i_n % q) - r * (i_n // q)
        if i_next < 0:
            return i_next + m
        else:
            return i_next

    def __make_number43(self):
        #创建队列数据结构
        q1 = Queue(maxsize=21)
        q2 = Queue(maxsize=22)
        t = self.__seed
        #使用16807随机生成前43个数
        for j in range(21):
            t = self.__iter_rand(t)
            q1.put(t)

        for j in range(22):
            t = self.__iter_rand(t)
            q2.put(t)

        self.__seed=t
        q = (q1,q2)
        return q
    
    def rand(self):
        i_0 = self.__numberQueue[0].get()
        i_21 = self.__numberQueue[1].get()
        i_43 = i_21 - i_0 - 0;
        #迭代循环
        if i_43<0:
            i_43 = i_43 + int(2**32-5) - 1
        self.__numberQueue[0].put(i_21)
        self.__numberQueue[1].put(i_43)
        return i_43/(4294967295) #i_43/(2**32-1)

def decision(list):
    if list[0]>list[2] and list[2]>list[1]:
        return 1
    else:
        return 0

        
if __name__ == "__main__":
    N = 10002
    # s = seed_time()
    # schrage_16807 = Schrage16807(seed=s)
    # data_16807 = [schrage_16807.rand() for i in range(N)]
    # fibonacci = LaggedFibonacci(seed=s)

    # data_fib = [fibonacci.rand() for i in range(N)]

    # df = pd.DataFrame()
    # df['16807']=data_16807
    # df['fibonacci']=data_fib
    # df.head()
    # df.to_csv("test_dataset.csv")

    iter_max = 10
    per_16807 = np.zeros(iter_max)
    per_fib = np.zeros(iter_max)
    for i in range(iter_max):
        s = seed_time()
        schrage_16807 = Schrage16807(seed=s)
        fibonacci = LaggedFibonacci(seed=s)
        data_16807 = [schrage_16807.rand() for i in range(N)]
        data_fib = [fibonacci.rand() for i in range(N)]
        N_16807 = 0
        N_fib = 0
        for j in range(1,N-1):
            N_16807 += decision(data_16807[j-1:j+2])
            N_fib += decision(data_fib[j-1:j+2])        
        per_16807[i] = N_16807/(N-2)
        per_fib[i] = N_fib/(N-2)
        print("Seed = %d, Schrage 16807: %.6f, Fibonacci: %.6f" %( s, per_16807[i], per_fib[i]))
        time.sleep(4)

    print("Averarge of Schrage 16807 : %.6f, Average of Fibonacci: %.6f" %( per_16807.mean(), per_fib.mean()))



