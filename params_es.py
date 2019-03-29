from data import GeneData,DataSaver
from tqdm import tqdm
import numpy as np
import pickle
import math
max_iterations=100
max_temperature=100
# temperature_rate=0.95

class ES():
    def __init__(self,max_att_iter,temperature_rate,save_path):
        self.max_trajectories=max_att_iter
        self.temperature_rate=temperature_rate
        self.save_path=save_path
        self.attractor_value=None



    def attractor_count(self):
        # init
        iter_count=0
        gene_data=GeneData()
        # avoid data saver
        del gene_data
        gene_data=GeneData()
        self.data_saver=DataSaver(gene_data)
        self.attractor_value = gene_data.attractor_define()
        break_counter = 0

        while iter_count<self.max_trajectories:
            break_counter+=1
            data=gene_data.trajectories_op(100,128)
            gap=gene_data._sum_value_for_data(data,gene_data.steps)-\
                gene_data._sum_value_for_data(data,gene_data.steps-2)
            # for a stable state
            # print(break_counter)
            if gap>gene_data.gene_num*gene_data.max_acc:
                # attractor value
                tmp_attractor_value = gene_data.attractor_define()
                for i in range(len(self.attractor_value)):
                    self.attractor_value[i]+=tmp_attractor_value[i]

            if gap > gene_data.gene_num * gene_data.max_acc or break_counter >= 300:
                iter_count += 1
                break_counter = 0
                if iter_count%10==0 and iter_count!=0:
                    print(iter_count,'/',self.max_trajectories,end='\t')

            del gene_data
            gene_data=GeneData()
            gene_data=self.data_saver.release_data(gene_data)

            # if break_counter>=350:
            #     break


    def cost_count(self):
        self.cost=self.max_trajectories

        if self.attractor_value!=None:
            for i in self.attractor_value:
                self.cost+=abs(self.max_trajectories/len(self.attractor_value)-i)-i
        print(self.attractor_value)

    def SA(self,iterations,temperature):
        # init
        current_cost=self.max_trajectories

        for i in range(1,iterations):
            self.attractor_count()
            self.cost_count()
            new_cost=self.cost
            if new_cost<=current_cost or\
                math.exp((current_cost-new_cost)/float(temperature))>np.random.random():
                current_cost=new_cost
                self.data_saver.save_data(self.save_path)
                print(current_cost,'Update !')
            temperature*=self.temperature_rate
            temperature=round(temperature,3)
            print('Epoch:',i)









es=ES(50,temperature_rate=max_temperature,save_path='./data/result.pkl')
es.SA(max_iterations,max_temperature)
