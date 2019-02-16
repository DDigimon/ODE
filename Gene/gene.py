import numpy as np

class gene():
    def __init__(self,DNAgene=None):
        # TODO init method
        self.value=np.random.random()
        self.activate_link=[]
        self.inactivate_link=[]
        self.function=0
        self.link_value={}
        self.DNAgene=DNAgene

    def init_param(self):
        for i in self.activate_link:
            self.link_value[i]={}
            # TODO init method
            self.link_value[i]['n']=np.random.random()
            self.link_value[i]['k']=np.random.random()
        for i in self.inactivate_link:
            self.link_value[i]={}
            self.link_value[i]['n']=np.random.random()
            self.link_value[i]['k']=np.random.random()

    def one_tune_value(self):
        for i in self.activate_link:
            self.link_value[i]['v']=self.act(i)
        for i in self.inactivate_link:
            self.link_value[i]['v']=self.in_act(i)

    def act(self,gene):
        return gene.value ** self.link_value[gene]['n'] / \
               (self.link_value[gene]['k'] ** self.link_value[gene]['n'] + gene.value ** self.link_value[gene]['n'])

    def in_act(self,gene):
        return self.link_value[gene]['k'] ** self.link_value[gene]['n'] / \
               self.link_value[gene]['k'] ** self.link_value[gene]['n'] + gene.value ** self.link_value[gene]['n']

    def or_op(self, gene_list):
        result = 0
        for gene in gene_list:
            result += self.link_value[gene]['v']
        return result

    def and_op(self, gene_list):
        result = 1
        for gene in gene_list:
            result *= self.link_value[gene]['v']
        return result


