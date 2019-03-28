import pickle
import math
import scipy.integrate as intergrates
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

with open('./data/final_data.pkl','rb') as f:
    gene_data=pickle.load(f)


class EnergyLandscape():
    '''
    1. u=-In(p)
    2. use Fokker–Planck (FP) equation and Wiener process get p=e^(-x^2/2t)/sqrt(2*PI*t)
    3. use normal distribution p=exp(-(x-u)^2/2*a^2)/sqrt(2*PI*a)
    4. each point of (x,y) have two solution of u according to 3
    5. u=ux*uy according to Gaussian distribution
    '''
    def __init__(self,gene_data,select_gene):
        self.gene_data=gene_data
        self.select_gene=select_gene

        self.gene_statistic_value={}
        for key in self.select_gene:
            self.gene_statistic_value[key]={}
            self.gene_statistic_value[key]['mean']=0
            self.gene_statistic_value[key]['var']=0

    def init_landscape(self,size,x_range,y_range):
        self.landscape=np.zeros((size))
        self.x_sample=np.linspace(x_range[0],x_range[1],size[0])
        self.y_sample=np.linspace(y_range[0],y_range[1],size[1])
        for key_gene in self.select_gene:
            # print(min(self.gene_data[key_gene]),max(self.gene_data[key_gene]))
            self.gene_statistic_value[key_gene]['mean']=np.mean(self.gene_data[key_gene])
            self.gene_statistic_value[key_gene]['var']=np.var(self.gene_data[key_gene])
        # print(self.landscape)
        # print(self.x_sample)
        # print(self.y_sample)

    def _quasi_potential(self,x,gene_name):
        p=math.exp(-(x-self.gene_statistic_value[gene_name]['mean'])**2/(2*self.gene_statistic_value[gene_name]['var']))\
          /math.sqrt(2*math.pi*self.gene_statistic_value[gene_name]['var'])
        u=-math.log(p)
        return u

    def quasi_landscape(self):
        self.ux=[]
        self.uy=[]
        for x in self.x_sample:
            self.ux.append(self._quasi_potential(x,self.select_gene[0]))
        for y in self.y_sample:
            self.uy.append(self._quasi_potential(y,self.select_gene[1]))

        for i in range(len(self.landscape)):
            for j in range(len(self.landscape[i])):
                self.landscape[i][j]=self.ux[i]*self.uy[j]


        self.X,self.Y=np.meshgrid(self.x_sample,self.y_sample)

        print(self.landscape)

    def draw_landscape(self):
        fig = plt.figure()
        axes3d = Axes3D(fig)
        axes3d.plot_surface(self.X, self.Y, self.landscape,cmap=cm.coolwarm)
        plt.show()


energy_land=EnergyLandscape(gene_data,['Sox1','Gata6'])
energy_land.init_landscape((20,20),(-0.5,1.5),(-0.5,1.5))
energy_land.quasi_landscape()
energy_land.draw_landscape()