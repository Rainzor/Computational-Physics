from audioop import add
from math import erf
import numpy as np
import matplotlib.pyplot as plt
from sympy import erfinv
import pandas as pd
from myStats import *
from Schrage_16807 import *
#make random data of size, using 16807 method to generate random numbers and seed need to using time


class student():
    def __init__(self,id=0) -> None:
        self.id=id

def add_id(student):
    student.id+=1

if __name__=="__main__":

    s = Schrage16807(seed=seed_time())
    g = s.normal(0,1,(2,100000))
    print(type(g))
    # n = norm()
    # h = n.rvs(100000)
    l = g.flatten()
    plt.hist(l, bins=100,range=(-5,5),density=True)
    plt.show()
    
