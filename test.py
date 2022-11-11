from audioop import add
from math import erf
import numpy as np
import matplotlib.pyplot as plt
from sympy import erfinv
import pandas as pd

from Schrage_16807 import Schrage16807
from Schrage_16807 import seed_time
#make random data of size, using 16807 method to generate random numbers and seed need to using time


class student():
    def __init__(self,id=0) -> None:
        self.id=id

def add_id(student):
    student.id+=1

if __name__=="__main__":

    score = [np.random.randint(0,10) for i in range(100)] # 此处随机生成一个数值列表
    score = pd.Series(score)
    se1 = pd.cut(score, [0,1,2,5,8,10],normalize=True) # 统计0-1,1-2依次类推各个区间的数值数量
    out= pd.value_counts(se1)
    print(out.values)
    print(out.index)
    
