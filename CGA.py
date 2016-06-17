# -*- coding: utf-8 -*-
import math, random
import matplotlib.pyplot as pl
ma = []
mi = []
ave = []

# chaotic genetic algorithm 混沌遗传算法
# 用logistics mapping 计算交叉算子， gray code

class Population:
    def __init__ (self, size, chrom_size, cp, mp, gen_max):
        # 种群信息
        self.individuals = []          # 个体集合
        self.fitness = []              # 个体适应度集合
        self.selector_probability = [] # 个体选择概率集合
        self.new_individuals = []      # 新一代个体集合
        
        self.elitist = {'chromosome':[0, 0], 'fitness':0, 'age':0} # 最佳个体信息
        
        self.size = size # 种群所包含的个体数
        self.chromosome_size = chrom_size # 个体的染色体长度
        self.crossover_probability = cp   # 个体之间的交叉概率
        self.mutation_probability = mp    # 个体之间的变异概率
        
        self.generation_max = gen_max # 种群进化的最大世代数
        self.age = 0                  # 种群当前所处世代
         
        # 随机产生初始个体集，并将新一代个体、适应度、选择概率等集合以 0 值进行初始化
        v = 2 ** self.chromosome_size - 1
        for i in range (self.size):
            self.individuals.append ([random.randint (0, v), random.randint (0, v),random.randint (0, v), random.randint (0, v)])
            self.new_individuals.append ([0, 0])
            self.fitness.append (0)
            self.selector_probability.append (0)
        
    def decode (self, interval, chromosome):
        d = interval[1] - interval[0]
        n = float (2 ** self.chromosome_size -1)
        return (interval[0] + chromosome * d / n)
    def encode (self, x):
        n = float (2 ** self.chromosome_size -1)
        return (int(n*x))
    
    def fitness_func (self, chrom1, chrom2):
        interval = [-10.0, 10.0]
        (x, y) = (self.decode (interval, chrom1), 
                  self.decode (interval, chrom2))
        n = lambda x, y: math.sin (math.sqrt (x*x + y*y)) ** 2 - 0.5
        d = lambda x, y: (1 + 0.001 * (x*x + y*y)) ** 2
        func = lambda x, y: 0.5 - n (x, y)/d (x, y)
        return func (x, y)
        
    def evaluate (self):
        sp = self.selector_probability
        for i in range (self.size):
            self.fitness[i] = self.fitness_func (self.individuals[i][0], 
                                                 self.individuals[i][1])
        ft_sum = sum (self.fitness)
        for i in range (self.size):
            sp[i] = self.fitness[i] / float (ft_sum)
        for i in range (1, self.size):
            sp[i] = sp[i] + sp[i-1]

    def select (self):
        (t, i) = (random.random (), 0)
        for p in self.selector_probability:
            if p > t:
                break
            i = i + 1
        return i
    
    def cross (self, chrom1, chrom2, mask10, mask20):
        p = random.random ()
        n = 2 ** self.chromosome_size -1
        (mask1, mask2) = (self.b2g(mask10), self.b2g(mask20))
        if chrom1 != chrom2 and p < self.crossover_probability:
            (h1, h2) = (chrom1 & mask1, chrom2 & mask2)
            (l1, l2) = (chrom1 & (~mask1), chrom2 & (~mask2))
            (chrom1, chrom2) = (h1 + l2, h2 + l1)
        return (chrom1, chrom2)
  
    def mutate (self, chrom):
        p = random.random ()  
        if p < self.mutation_probability:
            t = random.randint (1, self.chromosome_size)
            mask1 = 1 << (t - 1)
            mask2 = chrom & mask1
            if mask2 > 0:
                chrom = chrom & (~mask2)
            else:
                chrom = chrom ^ mask1
        return chrom
    
    def reproduct_elitist (self):
        # 与当前种群进行适应度比较，更新最佳个体
        j = 0
        for i in range (self.size):
            if self.elitist['fitness'] < self.fitness[i]:
                j = i
                self.elitist['fitness'] = self.fitness[i]
        if (j > 0):
            self.elitist['chromosome'][0] = self.individuals[j][0]
            self.elitist['chromosome'][1] = self.individuals[j][1]
            self.elitist['age'] = self.age
    def logistic(self, x):
        return 3.99*x*(1-x)
    def b2g(self, n):
        return n^(n >> 1)
    def g2b(self, n):
        mask = n >> 1
        while mask != 0:
            n = n ^ mask
            mask = mask >> 1
        return n

    def evolve (self):
        indvs = self.individuals
        new_indvs = self.new_individuals
        
        # 计算适应度及选择概率
        self.evaluate ()
        
        # 进化操作
        i = 0
        while True:
            # 选择两个个体，进行交叉与变异，产生新的种群
            idv1 = self.select ()
            idv2 = self.select ()
            
            # 交叉
            (idv1_x, idv1_y, mask01, mask02) = (indvs[idv1][0], indvs[idv1][1], indvs[idv1][2], indvs[idv1][3])
            (idv2_x, idv2_y, mask11, mask12) = (indvs[idv2][0], indvs[idv2][1], indvs[idv2][2], indvs[idv2][3])
            (idv1_x, idv2_x) = self.cross (idv1_x, idv2_x, mask01, mask02)
            (idv1_y, idv2_y) = self.cross (idv1_y, idv2_y, mask11, mask12)


            # 变异
            (idv1_x, idv1_y) = (self.mutate (idv1_x), self.mutate (idv1_y))
            (idv2_x, idv2_y) = (self.mutate (idv2_x), self.mutate (idv2_y))
            
            (new_indvs[i][0], new_indvs[i][1])     = (idv1_x, idv1_y)
            (new_indvs[i+1][0], new_indvs[i+1][1]) = (idv2_x, idv2_y)
            
            # 判断进化过程是否结束
            i = i + 2
            if i >= self.size:
                break
        
        # 最佳个体保留
        self.reproduct_elitist ()
        
        # 更新换代
        for i in range (self.size):
            self.individuals[i][0] = self.new_individuals[i][0]
            self.individuals[i][1] = self.new_individuals[i][1]
            a1 = self.individuals[i][2]
            a2 = self.individuals[i][3]
            inter = [0, 1]
            b1 = self.decode(inter, a1)
            b2 = self.decode(inter, a2)
            (c1, c2) = (self.logistic(b1), self.logistic(b2))
            (d1, d2) = (self.encode(c1), self.encode(c2))
            self.individuals[i][2] = self.b2g(d1)
            self.individuals[i][3] = self.b2g(d2)
  
    def run (self):
        for i in range (self.generation_max):
            self.age = i
            self.evolve ()
            ma.append(max (self.fitness))
            mi.append( min (self.fitness))
            ave.append(sum (self.fitness)/self.size)
            
if __name__ == '__main__':
    pop = Population (50, 25, 0.8, 0.1, 300)
    pop.run ()
    i = range(300)
    print max(ma)
    pl.plot(i, ma,'r', i, mi , 'b', i,ave,'y')
    pl.show()
