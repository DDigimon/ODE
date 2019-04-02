import numpy as np
import random
class ParamManager():
    def __init__(self):
        self.act_k=np.random.random()
        self.inact_k=np.random.random()
        self.act_n=int(random.randint(1,2))
        self.inact_n=int(random.randint(1,2))
    def gene_normal_init(self,size,DNA_size,mu=0.5,sigma=1):
        while True:
            self.gene_value=sigma * np.random.standard_normal(size) + mu
            self.gene_value=np.around(self.gene_value,decimals=3)
            if np.min(self.gene_value)>0:
                break
        while True:
            self.DNA_gene_value = sigma * np.random.standard_normal(DNA_size) + mu
            self.DNA_gene_value = np.around(self.DNA_gene_value, decimals=3)
            if np.min(self.DNA_gene_value) > 0:
                break

# params=ParamManager()
# params.gene_normal_init(16,5)
# print(params.gene_value)