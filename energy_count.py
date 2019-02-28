import pickle
import math
with open('./data/final_data.pkl','rb') as f:
    gene_data=pickle.load(f)

class Energy():
    def __init__(self,gene_data,chosed_gene):
        self.gene_data=gene_data
        self.chosed_gene=chosed_gene

    # def count_based(self):
    #     for key in self.chosed_gene:
    #

    def _function_define(self,v,m,gama):
        return 1/(1+math.exp((m-v)/float(gama)))



energy_counter=Energy(gene_data,['Gata6','Sox1'])

