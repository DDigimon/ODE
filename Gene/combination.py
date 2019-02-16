import numpy as np
class combination():
    def __init__(self,gene_list):
        self.a=np.random.random()
        self.gene_list=gene_list
        self.n_list=[]
        self.k_list=[]
        for gene in gene_list:
            self.n_list.append(np.random.random())
            self.k_list.append(np.random.random())