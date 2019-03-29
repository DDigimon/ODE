import numpy as np
import random
class ParamManager():
    def __init__(self):
        self.act_k=np.random.random()
        self.inact_k=np.random.random()
        self.act_n=int(random.randint(1,4))
        self.inact_n=int(random.randint(1,4))