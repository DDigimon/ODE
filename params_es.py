from data import GeneData
max_iterations=100
max_temperature=100

class ES():
    def __init__(self,max_att_iter):
        self.max_trajectories=max_att_iter

    def attractor_count(self):
        # init
        iter_count=0
        gene_data=GeneData()
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
                print(iter_count)
            del gene_data
            gene_data=GeneData()
        print(self.attractor_value)

    def cost_count(self):
        self.cost=self.max_trajectories
        for i in self.attractor_value:
            self.cost+=abs(self.max_trajectories/len(self.attractor_value)-i)-i
        print(self.cost)


es=ES(10)
es.attractor_count()
es.cost_count()
