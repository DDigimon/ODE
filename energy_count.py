import pickle
import math
import scipy.integrate as intergrates
import numpy as np

with open('./data/final_data.pkl','rb') as f:
    gene_data=pickle.load(f)

class Energy():
    def __init__(self,gene_data,chosed_gene):
        self.gene_data=gene_data
        self.chosed_gene=chosed_gene
        self.gene_value={}


    # def count_based(self):
    #     for key in self.chosed_gene:
    #

    def _function_define(self,v,m,gama):
        return 1/(1+math.exp((m-v)/float(gama)))

    def energy_count(self):
        for gene in self.chosed_gene:
            self.gene_value[gene]={}
            self.gene_value[gene]['m']=np.mean(self.gene_data[gene])
            self.gene_value[gene]['v']=np.var(self.gene_data[gene])
            # TODO is the value is same in equation?
            self.gene_value[gene]['gama']=np.random.random()
        print(self.gene_value)


        # intergrates.fixed_quad(function(x), 0, 0)



energy_counter=Energy(gene_data,['Gata6','Sox1'])
energy_counter.energy_count()

