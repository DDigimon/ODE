import pickle
import math
import scipy.integrate as intergrates
import numpy as np
from GPy.models.gplvm import GPLVM
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from final_data_generator import Generator

# with open('./data/final_data.pkl','rb') as f:
#     gene_data=pickle.load(f)


class EnergyLandscape():
    '''
    1. u=-In(p)
    2. use Fokker–Planck (FP) equation and Wiener process get p=e^(-x^2/2t)/sqrt(2*PI*t)
    3. use normal distribution p=exp(-(x-u)^2/2*a^2)/sqrt(2*PI*a)
    4. each point of (x,y) have two solution of u according to 3
    5. u=ux*uy according to Gaussian distribution
    '''
    def __init__(self,select_gene,gene_data_path,generator_path,mode=None,generator_final_data=False,
                 max_group=100,latent_dim=2,gene_sum_max=4,gene_sum_min=-4):
        self.select_gene=select_gene
        self.latent_dim=latent_dim
        self.mode=mode
        if len(select_gene)>2:
            self.mode='GPDM'

        generator_data=Generator(max_group=max_group,choose_gene=select_gene,gene_data_path=gene_data_path,
                                 result_path=generator_path,gene_sum_max=gene_sum_max,gene_sum_min=gene_sum_min)
        if generator_final_data==False:
            generator_data.load_gene_data()
        else:
            generator_data.generator()
        self.gene_data=generator_data.gene_dic

        if self.mode=='GPDM':
            self.gpdm_output=self.GPDM_solver()



    def init_landscape(self,size,x_range=None,y_range=None):
        self.landscape = np.ones((size))
        self.gene_statistic_value = {}
        for id,key_gene in enumerate(self.select_gene):
            if self.mode==None:
                self.gene_statistic_value[key_gene] = {}
                self.gene_statistic_value[key_gene]['mean'] = np.mean(self.gene_data[key_gene])
                self.gene_statistic_value[key_gene]['var'] = np.var(self.gene_data[key_gene])
            if self.mode=='GPDM':
                self.gene_statistic_value[id]={}
                self.gene_statistic_value[id]['mean']=np.mean(self.gene_data[key_gene])
                self.gene_statistic_value[id]['var']=np.var(self.gene_data[key_gene])


        if self.mode!='GPDM':
            self.x_sample = np.linspace(x_range[0], x_range[1], size[0])
            self.y_sample = np.linspace(y_range[0], y_range[1], size[1])


        else:
            self.x_value = []
            self.y_value = []
            self.landscape_point=np.zeros((size[0],size[1])).tolist()
            # print(self.gpdm_output)
            for i in range(len(self.select_gene)):
                self.x_value.append(self.gpdm_output.X[i][0])
                self.y_value.append(self.gpdm_output.X[i][1])
            self.x_value = np.array(self.x_value)
            self.y_value = np.array(self.y_value)
            x_start = self.x_value.min()+x_range[0]
            x_end = self.x_value.max()+x_range[1]
            y_start = self.y_value.min()+y_range[0]
            y_end = self.y_value.max()+y_range[1]
            self.x_sample=np.linspace(x_start,x_end,size[0])
            self.y_sample=np.linspace(y_start,y_end,size[1])
            for i in range(len(self.x_sample)):
                for j in range(len(self.y_sample)):
                    self.landscape_point[i][j]=self.gpdm_output.predict(np.array([[self.x_sample[i],self.y_sample[j]]]))[0][0]



    def _quasi_potential(self,x,mean,var):
        if var==0:
            var=0.001
        p=math.exp(-(x-mean)**2/(2*var))/math.sqrt(2*math.pi*var)
        if p==0:
            p=0.001
        u=-math.log(p)
        if u >=500:
            u=500
        return u

    def quasi_landscape(self):
        self.ux=[]
        self.uy=[]
        if self.mode=='GPDM':
            for i in range(len(self.x_sample)):
                for j in range(len(self.y_sample)):
                    u=[]
                    for id in range(len(self.landscape_point[i][j])):
                        u.append(self._quasi_potential(self.landscape_point[i][j][id],
                                                       self.gene_statistic_value[id]['mean'],
                                                       self.gene_statistic_value[id]['var']))

                    for u_value in u:
                        self.landscape[i][j]*=u_value
                        if self.landscape[i][j] >= 1000:
                            self.landscape[i][j] = 1000


        else:
            for x in self.x_sample:
                self.ux.append(self._quasi_potential(x,self.gene_statistic_value[self.select_gene[0]]['mean'],
                                                     self.gene_statistic_value[self.select_gene[0]]['var']))
            for y in self.y_sample:
                self.uy.append(self._quasi_potential(y,self.gene_statistic_value[self.select_gene[1]]['mean'],
                                                     self.gene_statistic_value[self.select_gene[1]]['var']))

            for i in range(len(self.landscape)):
                for j in range(len(self.landscape[i])):
                    self.landscape[i][j]=self.ux[i]*self.uy[j]
                    if self.landscape[i][j]>=1000:
                        self.landscape[i][j]=1000


        self.X,self.Y=np.meshgrid(self.x_sample,self.y_sample)

        # print(self.landscape)

    def draw_landscape(self):
        fig = plt.figure()
        axes3d = Axes3D(fig)
        axes3d.plot_surface(self.X, self.Y, self.landscape,cmap=cm.coolwarm)
        plt.show()

    def GPDM_solver(self):
        train_data=[]
        print(self.gene_data)
        for key_gene in self.select_gene:
            train_data.append(self.gene_data[key_gene])
        self.train_data=np.array(train_data).T
        output = GPLVM(self.train_data, self.latent_dim, init='PCA')
        output.optimize(messages=True,max_iters=20)
        return output


# gene_data_path='./data/result16.pkl'
# generator_path='./data/generator_result16.pkl'
#
# energy_land=EnergyLandscape(select_gene=['Oct4','Sox2','Nanog','Cdx2','Pax6','Sox1','Gata6','Myc','Klf4'],
#                             gene_data_path=gene_data_path,generator_path=generator_path,
#                             mode='GPDM',generator_final_data=False)
# # energy_land=EnergyLandscape(select_gene=['Sox1','Gata6'],
# #                             gene_data_path=gene_data_path,generator_path=generator_path,
# #                             generator_final_data=False)
# energy_land.init_landscape((30,30),(-0.5,1.5),(-0.5,1.5))
# energy_land.quasi_landscape()
# energy_land.draw_landscape()


import os
dir_path='./data/result6/'
file_lens=len(os.listdir(dir_path))
for file in os.listdir(dir_path):
    if 'init' not in file:
        print(file)
        gene_data_path = dir_path+file
        generator_path = './data/generator_result.pkl'

        # energy_land = EnergyLandscape(select_gene=['Oct4', 'Sox2', 'Nanog', 'Cdx2', 'Pax6', 'Sox1', 'Gata6', 'Myc', 'Klf4'],
        #                               gene_data_path=gene_data_path, generator_path=generator_path,
        #                               mode='GPDM', generator_final_data=True)
        energy_land=EnergyLandscape(select_gene=['Sox1','Gata6'],
                                    gene_data_path=gene_data_path,generator_path=generator_path,
                                    generator_final_data=True)
        energy_land.init_landscape((30, 30), (-0.5, 1.5), (-0.5, 1.5))
        energy_land.quasi_landscape()
        energy_land.draw_landscape()
        # energy_land=EnergyLandscape(select_gene=['Sox1','Oct4'],
        #                             gene_data_path=gene_data_path,generator_path=generator_path,
        #                             generator_final_data=False)
        # energy_land.init_landscape((30, 30), (-0.5, 1.5), (-0.5, 1.5))
        # energy_land.quasi_landscape()
        # energy_land.draw_landscape()
        #
        # energy_land = EnergyLandscape(select_gene=['Oct4', 'Sox2', 'Nanog', 'Cdx2', 'Pax6', 'Sox1', 'Gata6', 'Myc', 'Klf4'],
        #                               gene_data_path=gene_data_path, generator_path=generator_path,
        #                               mode='GPDM', generator_final_data=False)
        # energy_land.init_landscape((30, 30), (-0.5, 1.5), (-0.5, 1.5))
        # energy_land.quasi_landscape()
        # energy_land.draw_landscape()