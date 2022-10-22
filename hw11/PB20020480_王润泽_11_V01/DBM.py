from sympy import re
from DLA import*


class cell:
    def __init__(self,x=0,y=0,potential=0):
        self.position = (x,y)
        self.potential = potential
        self.growing_speed = 0 #growing_speed = phi**eta
        self.partical_sum = 0

    def __hash__(self) -> int:
        return hash(self.position)

    def __eq__(self, __o: object) -> bool:
        return self.position == __o.position
    

class DBM(Growth_Model):
    """Dielectric breakdown(DBM) refers to the formation of electrically conducting 
    regions in an insulating material exposed to a strong electric field.

    At the heart of the DBM algorithm is selecting where to grow 
    the pattern based on the electric potential at the candidate cells

    Attributes:
        N (int): the size of the grid

        max_num (int): the max number of particles

        eta (float): the parameter of the DBM algorithm,it varies from zero to infinity. 
                    Valuse of eta near zero give thick patterns, and larger values give more skinny branches.

        model_name (str): the name of the model

        r (Schrage16807): the random number generator

        radius (int): the radius of the boundery circle

        _particles_set (set): the set of the particles

        particle_number (int): the number of the particles

        _candidate_set (set): the set of the candidates

        phi0 (int): the potential of the central point

    Methods:
        _check_boundery(x,y): check if the point is in the boundery

        _get_potentical(x,y): get a potential of the candidate by random walk once

        _set_candidate_potential(candidate): set the potential of the candidate

        choose_candidate(): choose a candidate from the candidate list

        add_candidates(x,y): add the candidates around central point (x,y) to candidate set

        delete_candidate(candidate): delete the candidate from the candidate set

        add_particle(x,y): add the coordinate of candidate to the particle set

        update_candidate_set_potential(): update the potential of the candidates in the candidate set

        run(): run the DBM algorithm

        plot(): plot the result of the DBM algorithm

        save(): save the result of the DBM algorithm

    """
    def __init__(self, N=100,eta=2,max_num=400):
        self.N = N
        if(max_num < N**2//4):
            self.max_num = max_num  # 最大点数
        elif(max_num < N**2):
            raise WarningMessage(
                "The max number of particals is too large, it will be set to N**2//4")
        self.parameter = eta
        self.phi0 = 1
        self.model_name = "DBM"
        self.r = Schrage16807(seed_time())
        self.radius = N//2

        self._particles_set =  {(0,0)}
        self.particle_number = 1
        self._candidate_set = set()
        self.add_candidates(0,0)
        self.update_candidate_set_potential()

    # def random_walk(self,x,y):
    
    # def _check_point(self,x,y):


    def _check_boundery(self,x,y):
        """check if the point is in the boundery

        Args:
            x (int): x coordinate
            y (int): y coordinate

        Returns:
            value: the potential of the boundery
                   or -1 means in the boundery
        """
        if self._check_point(x,y):
            return 1
        elif x**2+y**2 <= self.radius**2:
            return -1
        else: 
            return 0

    def _get_potentical(self, x, y):
        """get a potential of the candidate by random walk once

        Args:
            x (int): x coordinate
            y (int): y coordinate
        """
        while 1:
            x, y = self.random_walk(x, y)
            c = self._check_boundery(x, y)
            if c < 0:
                continue
            else:
                return c

    def _set_candidate_potential(self,candidate,iter_num=20):
        """update the potential of the cadidate

            random work iter_num times to get the potential of the candidate

        Args:
            candidate (candidate): candidate to be updated the potential
        
        Returns:
            value: the potential of the candidate
        """          
        potential = 0
        for i in range(iter_num):
            x,y = candidate.position
            potential += self._get_potentical(x,y)
        candidate.potential = potential/iter_num
        return candidate.potential


    @property
    def eta(self):
        return self.parameter

##########################################################################
    def choose_candidate(self):
        """choose a candidate from the candidate list

        Returns:
            the candidate 
        """
        if self.eta>0:
            sum = 0
            for candidate in self._candidate_set:
                candidate.growing_speed = abs(self.phi0-candidate.potential)**self.eta
                sum += candidate.growing_speed
                candidate.partical_sum = sum

            p = self.r.rand()*sum
            for candidate in self._candidate_set:
                if p < candidate.partical_sum:
                    return candidate
        else:
            sum = len(self._candidate_set)
            p = self.r.rand()*sum
            for candidate in self._candidate_set:
                sum -= 1
                if p > sum:
                    return candidate


    def add_candidates(self,x,y):
        """add the candidates around central point (x,y) to candidate set

        Args:
            the central point's coordinate (x,y)
        """
        for i in range(-1,2):
            for j in range(-1,2):
                if self._check_point(x+i,y+j):
                    continue
                else:
                    self._candidate_set.add(cell(x+i,y+j))

    def delete_candidate(self,candidate):
        """delete the candidate from the candidate set

        Args:
            candidate(cell): the candidate to be remove in the set
        Returns:
            value: the candidate
        """
        if candidate in self._candidate_set:
            self._candidate_set.remove(candidate)
        else:
            raise KeyError("the candidate is not in the candidate set")
        return candidate

        
    def add_particle(self,x,y):
        """add the coordinate of candidate to the particle set

        Args:
            x (int): x coordinate
            y (int): y coordinate
        """
        self._particles_set.add((x,y))

    def update_candidate_set_potential(self):
        """update the candidate set
        """
        for candidate in self._candidate_set:
            self._set_candidate_potential(candidate)
            
    
    def run(self):
        """run the DBM algorithm

        代码主体分成以下6步：

        1.依生长为概率，选择候选点candidate

        2.将候选点加入粒子集合particle_set

        3.将候选点周围的点加入候选集合candidate_set,不设置其势能potential

        4.删除第2步加入粒子集合的候选点candidate

        5.依照随机行走的方法,更新候选点集合candidate_set中所有的势能potential

        """
        while self.particle_number < self.max_num:
            print("Particle number: ",self.particle_number)
            candidate = self.choose_candidate()
            self.add_particle(*candidate.position)
            self.add_candidates(*candidate.position)
            self.delete_candidate(candidate)
            self.update_candidate_set_potential()
            self.particle_number += 1

    # @property
    # def data(self):

    # def plot(self,**kwargs):

    # def save(self):



if __name__ == "__main__":
    dbm = DBM(N=300,eta = 6,max_num=300)


    time_start = time.time()
    dbm.run()
    time_end=time.time()
    min = (time_end-time_start)//60
    sec = (time_end-time_start)%60
    print('DBM {} points time cost: {}min, {}sec'.format(dbm.max_num,min,sec))
    dbm.plot(size=25,label='eta={}'.format(dbm.eta))

    dbm.save()



