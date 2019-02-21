import numpy as np
import random

class gene():
    def __init__(self,name,DNAgene=None):
        # TODO init method
        self.name=name
        self.value=np.random.random()
        self.kbase=np.random.random()
        self.act_value = np.random.random()
        self.inact_value = np.random.random()
        self.gama=np.random.random()

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
        self.value=np.random.random()

    def init_param(self):
        for i in self.activate_link:
            self.link_value[i.name]={}
            # TODO init method
            self.link_value[i.name]['n']=random.random()
            self.link_value[i.name]['k']=random.random()
        for i in self.inactivate_link:
            self.link_value[i.name]={}
            self.link_value[i.name]['n']=random.random()
            self.link_value[i.name]['k']=random.random()

    def one_tune_value(self):
        for i in self.activate_link:
            self.link_value[i.name]['v']=self.act(i)
        for i in self.inactivate_link:
            self.link_value[i.name]['v']=self.in_act(i)

    def act(self,gene):
        return gene.value ** self.link_value[gene.name]['n'] / \
               (self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] +
                gene.value ** self.link_value[gene.name]['n'])

    def in_act(self,gene):
        return self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] / \
               self.link_value[gene.name]['k'] ** self.link_value[gene.name]['n'] + \
               gene.value ** self.link_value[gene.name]['n']

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


