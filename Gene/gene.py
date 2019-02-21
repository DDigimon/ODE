import numpy as np
import random

class gene():
    def __init__(self,name,DNAgene=None):
        # TODO init method
        self.name=name
        self.value=round(np.random.random(),5)
        self.kbase=round(np.random.random(),5)
        self.act_value = round(np.random.random(),5)
        self.inact_value = round(np.random.random(),5)
        self.gama=round(np.random.random(),5)

        self.activate_link=[]
        self.inactivate_link=[]
        self.link_value={}
        self.DNAgene=DNAgene
        if self.DNAgene!=None:
            self.DNA_value=self.DNAgene.value
        else:
            self.DNA_value=1

        self.function_value = 0


    def init_gene_value(self):
        self.value=round(np.random.random(),5)

    def init_param(self):
        for i in self.activate_link:
            self.link_value[i.name]={}
            # TODO init method
            self.link_value[i.name]['n']=1
            self.link_value[i.name]['k']=round(np.random.random(),5)
        for i in self.inactivate_link:
            self.link_value[i.name]={}
            self.link_value[i.name]['n']=1
            self.link_value[i.name]['k']=round(np.random.random(),5)

    def one_tune_value(self):
        for i in self.activate_link:
            self.link_value[i.name]['v']=self.act(i)
        for i in self.inactivate_link:
            self.link_value[i.name]['v']=self.in_act(i)

    def act(self,gene):
        return round(gene.value ** self.link_value[gene.name]['n'] / \
               (self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] +
                gene.value ** self.link_value[gene.name]['n']),10)

    def in_act(self,gene):
        return round(self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] / \
               self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] + \
               gene.value ** self.link_value[gene.name]['n'],10)

    def or_op(self, gene_list):
        result = 0
        for gene in gene_list:
            result += self.link_value[gene.name]['v']
        return result

    def and_op(self, gene_list):
        result = 1
        for gene in gene_list:
            result *= self.link_value[gene.name]['v']
        return result



    def ODE_result(self):
        self.value=self.DNA_value*(self.kbase+self.function_value)-self.gama


