import numpy as np
import random

class gene():
    def __init__(self,name,DNAgene=None):
        # TODO init method
        # TODO acc control
        self.type='float64'
        self.max_acc=3
        self.nor_mu=0
        self.nor_sigm=0.5
        self.name=name
        self.value=round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)
        self.kbase=round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)
        self.act_value = round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)
        self.inact_value = round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)
        self.gama=round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)

        self.activate_link=[]
        self.inactivate_link=[]
        self.link_value={}
        self.iter_link_value={}
        self.DNAgene=DNAgene
        if self.DNAgene!=None:
            self.DNA_value=self.DNAgene.value
        else:
            self.DNA_value=1

        self.function_value = 0

    def _n_times(self,a,n):
        for _ in range(n):
            a*=a
            a=round(a,self.max_acc)
        if a==0.:
            a=0.1**self.max_acc
        return a

    def init_gene_value(self):
        self.value=round(np.random.normal(self.nor_mu,self.nor_sigm),self.max_acc)

    def init_param(self,act_n=random.randint(1,4),inact_n=random.randint(1,4),
                   act_k=round(np.random.normal(0,0.5),3),inact_k=round(np.random.normal(0,0.5),3)):
        for i in self.activate_link:
            self.link_value[i.name]={}
            # TODO init method
            self.link_value[i.name]['n']=act_n
            self.link_value[i.name]['k']=act_k
        for i in self.inactivate_link:
            self.link_value[i.name]={}
            self.link_value[i.name]['n']=inact_n
            self.link_value[i.name]['k']=inact_k

    def one_tune_value(self):
        for i in self.activate_link:
            self.iter_link_value[i.name]=self.act(i)
        for i in self.inactivate_link:
            self.iter_link_value[i.name]=self.in_act(i)

    def act(self,gene):
        return round(self._n_times(gene.value , self.link_value[gene.name]['n']) / \
               (self._n_times(self.link_value[gene.name]['k'] , self.link_value[gene.name]['n']) +
                self._n_times(gene.value , self.link_value[gene.name]['n'])),self.max_acc)

    def in_act(self,gene):
        return round(self._n_times(self.link_value[gene.name]['k'] , self.link_value[gene.name]['n']) / \
               self._n_times(self.link_value[gene.name]['k'] , self.link_value[gene.name]['n']) + \
               self._n_times(gene.value , self.link_value[gene.name]['n']),self.max_acc)

    def or_op(self, gene_list):
        result = 0
        for gene in gene_list:
            result += self.iter_link_value[gene.name]
        return result

    def and_op(self, gene_list):
        result = 1
        for gene in gene_list:
            result *= self.iter_link_value[gene.name]
            result=round(result,5)
        return result



    def ODE_result(self):
        self.value=round(self.DNA_value*(self.kbase+self.function_value)-self.gama,self.max_acc)


