import pickle
from data import GeneData,DataSaver
class Generator():
    def __init__(self,max_group,choose_gene,gene_data_path,result_path,gene_sum_min,gene_sum_max):
        self.step_recorder={}
        self.max_group=max_group
        self.gene_data_path=gene_data_path
        self.result_path=result_path
        self.gene_sum_min=gene_sum_min
        self.gene_sum_max=gene_sum_max
        self.choose_gene=choose_gene
        self.gene_dic={}
        for gene in choose_gene:
            self.gene_dic[gene] = []


    def generator(self):
        iter_count=0
        while iter_count < self.max_group:
            gene_data = GeneData()
            data_saver = DataSaver(gene_data)
            data_saver.load_data(self.gene_data_path)
            gene_data = data_saver.release_data(gene_data)
            data = gene_data.trajectories_op(100, 128)
            for gene in data.keys():
                self.step_recorder[gene]=[]
            # self.step_recorder.append(data)
            sum = 0
            for gene in self.choose_gene:
                sum += data[gene][len(data[gene]) - 1]

            del gene_data
            if sum > self.gene_sum_min and sum < self.gene_sum_max:
                iter_count += 1
                for gene in self.choose_gene:
                    self.gene_dic[gene].append(data[gene][len(data[gene]) - 1])
                for gene in data.keys():
                    self.step_recorder[gene].append(data[gene][len(data[gene])-1])
                if iter_count % 10 == 0:
                    print(iter_count, '/', self.max_group)
        print(self.gene_dic)
        self._save_gene_dic()

    def _save_gene_dic(self):
        with open(self.result_path, 'wb') as f:
            pickle.dump(self.step_recorder, f)
    def load_gene_data(self):
        with open(self.result_path,'rb') as f:
            self.gene_dic=pickle.load(f)
