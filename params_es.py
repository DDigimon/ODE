from data import GeneData,DataSaver
from tqdm import tqdm
import pickle
import math
max_iterations=100
max_temperature=100
temperature_rate=0.95

class ES():
    def __init__(self,max_att_iter):
        self.max_trajectories=max_att_iter

    def attractor_count(self):
        # init
        iter_count=0
        gene_data=GeneData()
        # avoid data saver
        del gene_data
        gene_data=GeneData()
        self.data_saver=DataSaver(gene_data)
        self.attractor_value = gene_data._attractor_define()

        while iter_count<self.max_trajectories:
            data=gene_data.trajectories_op(100,1,128)[0]
            gap=gene_data._sum_value_for_data(data,gene_data.steps)-\
                gene_data._sum_value_for_data(data,gene_data.steps-2)
            # for a stable state
            if gap>gene_data.gene_num*gene_data.max_acc:
                # attractor value
                tmp_attractor_value = gene_data._attractor_define()
                for i in range(len(self.attractor_value)):
                    self.attractor_value[i]+=tmp_attractor_value[i]
                iter_count += 1
            del gene_data
            gene_data=GeneData()
            gene_data=self.data_saver.release_data(gene_data)


    def cost_count(self):
        self.cost=self.max_trajectories
        for i in self.attractor_value:
            self.cost+=abs(self.max_trajectories/len(self.attractor_value)-i)-i

    def SA(self,iterations,temperature):
        # init
        self.attractor_count()
        self.cost_count()
        current_cost=self.cost
        current_data=self.data_saver
        for _ in tqdm(range(iterations)):
            self.attractor_count()
            self.cost_count()
            new_cost=self.cost
            new_data=self.data_saver
            if new_cost<=current_cost or\
                math.exp((current_cost-new_cost)/float(temperature)):
                current_cost=new_cost
                current_data=new_data
            temperature*=temperature_rate
            temperature=round(temperature,3)
        with open('./data/result.pkl','rb',encoding='utf-8') as f:
            pickle.dump(current_data,f)








es=ES(10)
es.attractor_count()
es.cost_count()
es.SA(max_iterations,max_temperature)
