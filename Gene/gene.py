import numpy as np

class gene():
    def __init__(self):
        # TODO init method
        self.init_value=np.random.random()
        self.activate_link=[]
        self.inactivate_link=[]