from DBM import *


class DBM_FAST(DBM):
    """
    This DBM using the faster method to calculate the potential

    Attributes:
        N (int): the size of the grid

        max_num (int): the max number of particles

        eta (int): the parameter of the DBM algorithm,it varies from zero to infinity. 
                   Valuse of eta near zero give thick patterns, and larger values give more skinny branches.

        model_name (str): the name of the model

        r (Schrage16807): the random number generator

        radius (int): the radius of the circle

        _particles_set (set): the set of the particles

        particle_number (int): the number of the particles

        _candidate_set (set): the set of the candidates

        phi0 (int): the potential of the central point

        max_potential (int): the max potential of the candidate set

        min_potential (int): the min potential of the candidate set

    Methods:
        _get_potential(x_e,y_e,x_p,y_p): get (x_p,y_p) potential(phi) from an electron at position (x_e,y_e) created

        _set_max_min_potential(p): set the max and min potential

        _set_candidate_growing_speed(candidate): set the growing_speed of the candidate

        _set_new_candidate_potential(new_candidate): set the potential of the candidate

        choose_candidate(): choose a candidate from the candidate list

        update_candidate_set_potential(x,y): update the candidate set value

        add_candidates(x,y): add the candidates of the position (x,y)

        add_particle(x,y): add a particle at the position (x,y)

        run(): run the DBM algorithm

        plot(): plot the result

        save(): save the result

        data(): return the result

    """

    def __init__(self, N=200, eta=2, max_num=400):
        self.N = N
        if(max_num < N**2//4):
            self.max_num = max_num  # 最大点数
        elif(max_num < N**2):
            raise WarningMessage(
                "The max number of particals is too large, it will be set to N**2//4")
        self.parameter = eta

        self.model_name = "DBM_FAST"
        self.r = Schrage16807(seed_time())
        self.radius = N//2
        self.R1 = 1/2
        self.R2 = N/2
        self._particles_set = {(0, 0)}
        self.particle_number = 1
        self._candidate_set = set()
        self.max_potential = 0
        self.min_potential = np.inf
        self.add_candidates(0, 0)

            

    def _get_potential(self,x_e,y_e,x_p,y_p):
        """get (x_p,y_p) potential(phi) from an electron at position (x_e,y_e) created

        Args:
            x_e (int): x coordinate of the electron
            y_e (int): y coordinate of the electron
            x_p (int): x coordinate of the position
            y_p (int): y coordinate of the position

        """
        r = np.sqrt((x_e-x_p)**2+(y_e-y_p)**2)
        #return 1-self.R1/r
        return self.R1/r

    def _set_max_min_potential(self,p):
        if p>self.max_potential:
            self.max_potential = p
        if p<self.min_potential:
            self.min_potential = p

    def _set_candidate_growing_speed(self,candidate):
        """set the growing_speed of the candidate

        """
        #p_e = (candidate.potential-self.min_potential)/(self.max_potential-self.min_potential)
        p_e = (self.max_potential-candidate.potential )/(self.max_potential-self.min_potential)
        candidate.growing_speed = p_e**self.eta
        return candidate.growing_speed

    # def _check_point(self, x, y):
    #     return (x, y) in self._particles_set


    def _set_new_candidate_potential(self,new_candidate):
        """set the potential of the candidate

        Args:
            new candidate(cell): the candidate to be set the potential
        """
        potential = 0
        x,y = new_candidate.position
        for p_x,p_y in self._particles_set:
            potential += self._get_potential(p_x,p_y,x,y)
        new_candidate.potential = potential
        return potential

############################################################################################################
    def choose_candidate(self):
        """choose a candidate from the candidate list

        Returns:
            a candidate 
        """
        if self.eta > 0:
            self.max_potential = 0
            self.min_potential = np.inf
            sum = 0
            for candidate in self._candidate_set:
                p = candidate.potential
                self._set_max_min_potential(p)

            for candidate in self._candidate_set:
                sum += self._set_candidate_growing_speed(candidate)
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



    # def add_particle(self, x, y):

    # def delete_candidate(self, candidate):

    def update_candidate_set_potential(self,x,y):
        """update the candidate set value
            
        """     
        for candidate in self._candidate_set:
            candidate.potential += self._get_potential(
                x,y, candidate.position[0], candidate.position[1])
        
    def add_candidates(self, x, y):
        """add the candidates around central point (x,y) to candidate set

        Args:
            the central point's coordinate (x,y)
        """
        for i in range(-1, 2):
            for j in range(-1, 2):
                if self._check_point(x+i, y+j) :#检查是否已经存在在粒子集合中
                    continue
                else:
                    new_candidate = cell(x+i, y+j)
                    self._set_new_candidate_potential(new_candidate)
                    self._candidate_set.add(new_candidate)
    


    def run(self):
        """run the DBM algorithm

        代码主体分成以下6步：

        1.依生长速率为概率, 选择候选点candidate

        2.将候选点加入粒子集合particle_set

        3.删除第2步加入粒子集合的候选点candidate

        4.更新候选点集合candidate_set中的势能potential

        5.将候选点周围的点加入候选集合candidate_set, 并且设置其位置与势能

        """
        while self.particle_number < self.max_num:
            print("Particle number: ", self.particle_number)
            candidate = self.choose_candidate()
            self.add_particle(*candidate.position)
            self.delete_candidate(candidate)
            self.update_candidate_set_potential(*candidate.position)
            self.add_candidates(*candidate.position)
            #self.update_remain_values()
            self.particle_number += 1

    #def data(self):

    #def plot(self):

    #def save(self):

if __name__ == "__main__":
    dbm_f = DBM_FAST(N=300,eta = 6,max_num=300)

    #计时
    time_start = time.time()
    dbm_f.run()
    time_end=time.time()
    min = (time_end-time_start)//60
    sec = (time_end-time_start)%60
    print('DBM FAST {} points time cost: {}min, {}sec'.format(dbm_f.max_num,min,sec))


    dbm_f.plot(size=10,label='eta={}'.format(dbm_f.eta))

    dbm_f.save()

    # df = pd.read_csv('hw11/DBM_FAST.csv')
    # x = df['x']
    # y = df['y']
    # #画图
    # size = 5
    # label = 'eta={}'.format(6)
    # plt.style.use('dark_background')
    # plt.scatter(x, y, s=size, marker=',', color='white',
    #             edgecolors='none', label=label)
    # plt.title("DBM_FAST")
    # plt.legend()
    # #隐藏坐标轴
    # plt.axis('square')
    # plt.xticks([])
    # plt.yticks([])
    # plt.show()
