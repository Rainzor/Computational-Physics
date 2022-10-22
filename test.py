from audioop import add
from math import erf
import numpy as np
import matplotlib.pyplot as plt
from sympy import erfinv

from Schrage_16807 import Schrage16807
from Schrage_16807 import seed_time
#make random data of size, using 16807 method to generate random numbers and seed need to using time


class student():
    def __init__(self,id=0) -> None:
        self.id=id

def add_id(student):
    student.id+=1

if __name__=="__main__":
    s = student()
    add_id(s)
    print(s.id)
    add_id(s)
    print(s.id)