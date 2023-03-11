import multiprocessing as mp
from time import sleep
from Schrage_16807 import*
class cell:
    def __init__(self, x=0, y=0, potential=0):
        self.position = (x, y)
        self.potential = potential
        self.phi_eta = 0  # phi_eta = phi*eta
        self.partical_sum = 0

    def __hash__(self) -> int:
        return hash(self.position)

    def __eq__(self, __o: object) -> bool:
        return self.position == __o.position
    def __str__(self) -> str:
        return str(self.position)+":"+str(self.potential)

def get_rand(x,y):
    return Schrage16807(seed_time()).rand()*x-y

if __name__=="__main__":
    for i in range(20):
        print(seed_time())
        sleep(1)
