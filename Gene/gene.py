import numpy as np

class gene():
    def __init__(self,DNAgene=None):
        # TODO init method
        self.init_value=np.random.random()
        self.activate_link=[]
        self.inactivate_link=[]
        self.function=0
        self.DNAgene=DNAgene
        self.link_value={}

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