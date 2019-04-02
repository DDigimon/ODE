from Gene.gene import gene
from Gene.merge_gene import merge_gene
from Gene.param_manager import ParamManager
import pickle

import numpy as np


class GeneData():
    def __init__(self):
        # TODO use special method to load data, some mode?
        self.gene_list=[]
        self.DNA_gene_list=[]
        self.gene_name_dic={}
        self.gene_id_dic={}
        self.trajectories_num=0
        self.gene_num=0
        self.max_acc=4

        self.Oct4_gene=gene('Oct4_gene')
        self.Sox2_gene=gene('Sox2_gene')
        self.Nanog_gene=gene('Nanog_gene')
        self.Gata6_gene=gene('Gata6_gene')
        self.Sox1_gene=gene('Sox1_gene')

        self.Oct4=gene('Oct4',DNAgene=self.Oct4_gene)
        self.Sox2=gene('Sox2',DNAgene=self.Sox2_gene)
        self.Nanog=gene('Nanog',DNAgene=self.Nanog_gene)
        self.Cdx2=gene('Cdx2')
        self.Sox1=gene('Sox1',DNAgene=self.Sox1_gene)
        self.Gcnf=gene('Gcnf')
        self.Pax6=gene('Pax6')
        self.Gata6=gene('Gata6',DNAgene=self.Gata6_gene)
        self.Myc=gene('Myc')
        self.Klf4=gene('Klf4')
        self.Prc2=gene('Prc2')
        self.Mbd3=gene('Mbd3')
        self.Nurd=gene('Nurd')
        self.EA1=gene('EA1')
        self.EA2=gene('EA2')
        self.EA3=gene('EA3')

        self.gene_list=[self.Oct4,self.Sox2,self.Nanog,self.Cdx2,self.Sox1,self.Gcnf,self.Pax6,self.Gata6,
                        self.Myc,self.Klf4,self.Prc2,self.Mbd3,self.Nurd,self.EA1,self.EA2,self.EA3]

        self.DNA_gene_list=[self.Oct4_gene,self.Nanog_gene,self.Gata6_gene,self.Sox1_gene,self.Sox2_gene]
        self.gene_num=len(self.gene_list)+len(self.DNA_gene_list)

        self.OC=merge_gene('OC',[self.Oct4,self.Cdx2])
        self.OS=merge_gene('OS',[self.Oct4,self.Sox2])
        self.OSN=merge_gene('OSN',[self.Oct4,self.Sox2,self.Nanog])
        self.OG=merge_gene('OG',[self.Oct4,self.Gata6])

        for each_gene in self.gene_list:
            self.gene_id_dic[each_gene] = each_gene.name
            self.gene_name_dic[each_gene.name] = each_gene
        for each_gene in self.DNA_gene_list:
            self.gene_id_dic[each_gene] = each_gene.name
            self.gene_name_dic[each_gene.name] = each_gene

        self.params=ParamManager()
        self._add_link()

    def _add_link(self):
        # TODO trajectories_num count
        self.Oct4.activate_link=[self.OS,self.Klf4,self.OSN,self.Myc]
        self.Oct4.inactivate_link=[self.Gata6,self.Sox1,self.OC,self.Gcnf]
        self.trajectories_num+=9

        self.Sox2.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Sox2.inactivate_link=[self.Sox1,self.Gata6,self.Oct4]
        self.trajectories_num+=8

        self.Nanog.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Nanog.inactivate_link=[self.OG,self.Sox1,self.Oct4]
        self.trajectories_num+=4

        self.Cdx2.activate_link=[self.Cdx2]
        self.Cdx2.inactivate_link=[self.Nanog,self.OC]
        self.trajectories_num+=1

        self.Gcnf.activate_link=[self.Cdx2,self.Gata6,self.Gcnf]
        self.trajectories_num+=1

        self.Pax6.activate_link=[self.Sox2,self.Pax6]
        self.Pax6.inactivate_link=[self.Oct4,self.Nanog]
        self.trajectories_num+=2

        self.Sox1.activate_link=[self.Sox2]
        self.Sox1.inactivate_link=[self.Gata6,self.Oct4,self.Nanog]
        self.trajectories_num+=4

        self.Gata6.activate_link=[self.Oct4]
        self.Gata6.inactivate_link=[self.Sox1,self.Sox2,self.Nanog]
        self.trajectories_num+=4

        self.Myc.activate_link=[self.Klf4]
        self.Myc.inactivate_link=[self.Sox1,self.Gata6,self.Cdx2,self.Gcnf]
        self.trajectories_num+=2

        self.Klf4.activate_link=[self.Klf4]
        self.Klf4.inactivate_link=[self.Sox1,self.Cdx2,self.Gata6,self.Gcnf]
        self.trajectories_num+=2

        self.Prc2.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.Prc2.inactivate_link=[self.Mbd3,self.Nurd]
        self.trajectories_num+=2

        self.Mbd3.activate_link=[self.Gata6]
        self.Mbd3.inactivate_link=[self.Prc2,self.Nurd]
        self.trajectories_num+=1

        self.Nurd.activate_link=[self.Sox1]
        self.Nurd.inactivate_link=[self.Prc2,self.Mbd3]
        self.trajectories_num+=1

        self.EA1.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.EA1.inactivate_link=[self.EA2,self.EA3]
        self.trajectories_num+=2

        self.EA2.activate_link=[self.Gata6]
        self.EA2.inactivate_link=[self.EA1,self.EA3]
        self.trajectories_num+=1

        self.EA3.activate_link=[self.Sox1]
        self.EA3.inactivate_link=[self.EA1,self.EA2]
        self.trajectories_num+=1

        for gene in self.gene_list:
            gene.init_param()

        self.Oct4_gene.activate_link=[self.EA1]
        self.Oct4_gene.inactivate_link=[self.Mbd3,self.Nurd]

        self.Sox2_gene.activate_link=[self.EA1]
        self.Sox2_gene.inactivate_link=[self.Mbd3,self.Nurd]

        self.Nanog_gene.activate_link=[self.EA1]
        self.Nanog_gene.inactivate_link=[self.Mbd3,self.Nurd]

        self.Gata6_gene.activate_link=[self.EA2]
        self.Gata6_gene.inactivate_link=[self.Prc2,self.Nurd]

        self.Sox1_gene.activate_link=[self.EA3]
        self.Sox1_gene.inactivate_link=[self.Prc2,self.Mbd3]

        for gene in self.DNA_gene_list:
            gene.init_param()
            # self.trajectories_num+=len(gene.activate_link)+len(gene.inactivate_link)

    def _sum_value_for_data(self,data,n):
        '''

        :param data:
        :param n: the n(th)
        :return:
        '''
        sum_value=0
        n-=1
        for gene in self.gene_list:
            sum_value+=data[gene.name][n]
        for gene in self.DNA_gene_list:
            sum_value+=data[gene.name][n]
        return sum_value

    def attractor_define(self):
        attractor=[]
        attractor.append(int(bool(self.Sox1.value>self.Oct4.value and
                                   self.Sox1.value>self.Sox2.value and
                                   self.Sox1.value>self.Gata6.value)))
        attractor.append(int(bool(self.Gata6.value>self.Oct4.value and
                                   self.Gata6.value>self.Sox2.value and
                                   self.Gata6.value>self.Sox1.value)))
        attractor.append(int(bool(self.Oct4.value>self.Sox1.value and
                                   self.Oct4.value>self.Gata6.value and
                                   self.Sox2.value>self.Gata6.value and
                                   self.Sox2.value>self.Sox1.value)))

        return attractor


    def one_tune(self):
        # for gene in self.gene_list:
        #     gene.one_tune_value()
        #
        # for gene in self.DNA_gene_list:
        #     gene.one_tune_value()
        #
        def count_func(n, k, gene_list, mode='a'):
            value = 1
            for gene in gene_list:
                value *= gene.value
            if (k ** n + value ** n)==0:
                return 0.001
            if mode == 'a':
                return round(value ** n / (k ** n + value ** n), 3)
            if mode == 'i':
                return round(k ** n / (k ** n + value ** n), 3)

        # a1,a2,a3,a4,a5,a6,a7,a8,a9,a10=0.4802,0.7065,0.6033,0.2405,0.8269,0.0386,0.7158,0.3503,0.6858,0.7
        # n1,n2 ,n3,n4=1,5,7,2
        # k1,k2,k3,k4=0.1,0.4,0.0252,0.4086
        # a=0.1
        # b=0.01
        a=0.1
        a1,a2,a3,a4,a5,a6,a7,a8,a9,a10=0.4802,0.7065,0.6033,0.2405,0.8269,0.0386,0.7158,0.3503,0.6858,0.7
        # n1,n2 ,n3,n4=1,5,7,2
        # k1,k2,k3,k4=0.1,0.4,0.0252,0.4086
        ka,ki=0.1,0.4
        na,ni=1,5
        a=0.1
        aa=0.7
        b=0.01
        kb=0.11
        self.Oct4.function_value=a*(count_func(na,ka,[self.Oct4,self.Sox2])+count_func(na,k1,[self.Klf4])+
                                    count_func(na,k1,[self.Oct4,self.Sox2,self.Nanog])+count_func(na,k1,[self.Myc]))*\
                                 count_func(ni,k2,[self.Gata6],mode='i')*count_func(n2,k2,[self.Sox1],mode='i')+\
                                 (count_func(n2,k2,[self.Oct4,self.Cdx2],mode='i')+count_func(n2,k2,[self.Gcnf],mode='i')*b)

        self.Sox2.function_value=a*(count_func(na,k1,[self.Oct4,self.Sox2])+
                                    count_func(na,k1,[self.Oct4,self.Sox2,self.Nanog])+
                                    count_func(na,k1,[self.Myc])+count_func(na,k1,[self.Klf4]))*\
                                 count_func(n2,k2,[self.Sox1],mode='i')*count_func(ni,k2,[self.Gata6],mode='i')+\
                                 count_func(n3,k2,[self.Oct4],mode='i')*b

        self.Nanog.function_value=a*(count_func(na,k1,[self.Oct4,self.Sox2])+
                                     count_func(na,k1,[self.Oct4,self.Sox2,self.Nanog])+count_func(na,k1,[self.Myc])+
                                     count_func(na,k1,[self.Klf4]))*count_func(n2,k2,[self.Oct4,self.Gata6],mode='i')*\
                                  count_func(n2,k2,[self.Sox1],mode='i')+count_func(n2,k2,[self.Oct4],mode='i')*b

        self.Cdx2.function_value=a10*count_func(na,k1,[self.Cdx2])*count_func(n3,k2,[self.Nanog],mode='i')*\
                                 count_func(n3,k2,[self.Oct4,self.Cdx2],mode='i')*b

        self.Gcnf.function_value=a*a10*(count_func(na,k1,[self.Cdx2])+count_func(na,k2,[self.Nanog])+count_func(na,k2,[self.Gcnf]))

        self.Pax6.function_value=a*a10*(count_func(na,k1,[self.Sox2])+count_func(na,k1,[self.Gata6]))*\
                                 count_func(n3,k2,[self.Oct4],mode='i')*count_func(n3,k2,[self.Nanog],mode='i')*b

        self.Sox1.function_value=a*a10*count_func(na,k1,[self.Sox2])*count_func(ni,k2,[self.Gata6],mode='i')*\
                                 count_func(n3,k2,[self.Oct4],mode='i')+count_func(n3,k2,[self.Nanog],mode='i')*b

        self.Gata6.function_value=a*a10*count_func(na,k1,[self.Oct4])*count_func(n3,k2,[self.Sox1],mode='i')*\
                                  count_func(n3,k2,[self.Sox2],mode='i')+count_func(n3,k3,[self.Nanog])*b

        self.Myc.function_value=a*count_func(na,k1,[self.Klf4])*count_func(n3,k2,[self.Sox1],mode='i')*\
                                count_func(n3,k2,[self.Gata6],mode='i')*count_func(n3,k2,[self.Cdx2],mode='i')+\
                                count_func(n3,k2,[self.Gcnf],mode='i')*b

        self.Klf4.function_value=a*count_func(na,k1,[self.Klf4])*count_func(n3,k2,[self.Sox1],mode='i')*\
                                 count_func(n3,k2,[self.Cdx2],mode='i')*count_func(n3,k2,[self.Gata6],mode='i')+\
                                 count_func(n3,k2,[self.Gcnf],mode='i')*b

        self.Prc2.function_value=a1*(count_func(n4,k3,[self.Oct4])+count_func(n4,k3,[self.Sox2])+
                                     count_func(n4,k4,[self.Myc])+count_func(n4,k4,[self.Klf4]))+\
                                 a2*(count_func(n4,k4,[self.Mbd3])+count_func(n4,k4,[self.Nurd]))

        self.Mbd3.function_value=a3*count_func(n4,k3,[self.Gata6])+\
                                 a4*count_func(n4,k4,[self.Prc2],mode='i')+count_func(n4,k4,[self.Mbd3],mode='i')

        self.Nurd.function_value=a3*count_func(n4,k3,[self.Sox1])+\
                                 a4*count_func(n4,k4,[self.Prc2],mode='i')+count_func(n4,k4,[self.Nurd],mode='i')

        self.EA1.function_value=a5*(count_func(n4,k3,[self.Oct4])+count_func(n4,k3,[self.Sox2])+
                                    count_func(n4,k4,[self.Myc])+count_func(n4,k4,[self.Klf4]))+\
                                a6*(count_func(n4,k4,[self.EA2],mode='i')+count_func(n4,k4,[self.EA3],mode='i'))

        self.EA2.function_value=a7*(count_func(n4,k3,[self.Gata6])+count_func(n4,k4,[self.EA1],mode='i')+
                                    count_func(n4,k4,[self.EA3],mode='i'))

        self.EA3.function_value=a7*(count_func(n4,k3,[self.Sox1])+count_func(n4,k4,[self.EA1],mode='i')+
                                    count_func(n4,k4,[self.EA2],mode='i'))




        self.Oct4_gene.value=0

        self.Sox2_gene.value=0

        self.Nanog_gene.value=0

        self.Gata6_gene.value=0

        self.Sox1_gene.value=0
        # self.Oct4.function_value=a*(count_func(n1,k1,[self.Oct4,self.Sox2])+count_func(n1,k1,[self.Klf4])+
        #                             count_func(n1,k1,[self.Oct4,self.Sox2,self.Nanog])+count_func(n1,k1,[self.Myc]))*\
        #                          count_func(n2,k2,[self.Gata6],mode='i')*count_func(n2,k2,[self.Sox1],mode='i')+\
        #                          (count_func(n2,k2,[self.Oct4,self.Cdx2],mode='i')+count_func(n2,k2,[self.Gcnf])*b)
        #
        # self.Sox2.function_value=a*(count_func(n1,k1,[self.Oct4,self.Sox2])+
        #                             count_func(n1,k1,[self.Oct4,self.Sox2,self.Nanog])+
        #                             count_func(n1,k1,[self.Myc])+count_func(n1,k1,[self.Klf4]))*\
        #                          count_func(n2,k2,[self.Sox1],mode='i')*count_func(n2,k2,[self.Gata6],mode='i')+\
        #                          count_func(n3,k2,[self.Oct4],mode='i')*b
        #
        # self.Nanog.function_value=a*(count_func(n1,k1,[self.Oct4,self.Sox2])+
        #                              count_func(n1,k1,[self.Oct4,self.Sox2,self.Nanog])+count_func(n1,k1,[self.Myc])+
        #                              count_func(n1,k1,[self.Klf4]))*count_func(n2,k2,[self.Oct4,self.Gata6],mode='i')*\
        #                           count_func(n2,k2,[self.Sox1],mode='i')+count_func(n2,k2,[self.Oct4],mode='i')*b
        #
        # self.Cdx2.function_value=a10*count_func(n1,k1,[self.Cdx2])*count_func(n3,k2,[self.Nanog],mode='i')*\
        #                          count_func(n3,k2,[self.Oct4,self.Cdx2],mode='i')*b
        #
        # self.Gcnf.function_value=a*a10*(count_func(n1,k1,[self.Cdx2])+count_func(n3,k2,[self.Nanog])+count_func(n3,k2,[self.Gcnf]))
        #
        # self.Pax6.function_value=a*a10*(count_func(n1,k1,[self.Sox2])+count_func(n1,k1,[self.Gata6]))*\
        #                          count_func(n3,k2,[self.Oct4],mode='i')*count_func(n3,k2,[self.Nanog],mode='i')*b
        #
        # self.Sox1.function_value=a*a10*count_func(n1,k1,[self.Sox2])*count_func(n3,k2,[self.Gata6],mode='i')*\
        #                          count_func(n3,k2,[self.Oct4],mode='i')+count_func(n3,k2,[self.Nanog],mode='i')*b
        #
        # self.Gata6.function_value=a*a10*count_func(n1,k1,[self.Oct4])*count_func(n3,k2,[self.Sox1],mode='i')*\
        #                           count_func(n3,k2,[self.Sox2],mode='i')+count_func(n3,k3,[self.Nanog])*b
        #
        # self.Myc.function_value=a*count_func(n1,k1,[self.Klf4])*count_func(n3,k2,[self.Sox1],mode='i')*\
        #                         count_func(n3,k2,[self.Gata6],mode='i')*count_func(n3,k2,[self.Cdx2],mode='i')+\
        #                         count_func(n3,k2,[self.Gcnf],mode='i')*b
        #
        # self.Klf4.function_value=a*count_func(n1,k1,[self.Klf4])*count_func(n3,k2,[self.Sox1],mode='i')*\
        #                          count_func(n3,k2,[self.Cdx2],mode='i')*count_func(n3,k2,[self.Gata6],mode='i')+\
        #                          count_func(n3,k2,[self.Gcnf],mode='i')*b
        #
        # self.Prc2.function_value=a1*(count_func(n4,k3,[self.Oct4])+count_func(n4,k3,[self.Sox2])+
        #                              count_func(n4,k4,[self.Myc])+count_func(n4,k4,[self.Klf4]))+\
        #                          a2*(count_func(n4,k4,[self.Mbd3])+count_func(n4,k4,[self.Nurd]))
        #
        # self.Mbd3.function_value=a3*count_func(n4,k3,[self.Gata6])+\
        #                          a4*count_func(n4,k4,[self.Prc2],mode='i')+count_func(n4,k4,[self.Mbd3],mode='i')
        #
        # self.Nurd.function_value=a3*count_func(n4,k3,[self.Sox1])+\
        #                          a4*count_func(n4,k4,[self.Prc2],mode='i')+count_func(n4,k4,[self.Nurd],mode='i')
        #
        # self.EA1.function_value=a5*(count_func(n4,k3,[self.Oct4])+count_func(n4,k3,[self.Sox2])+
        #                             count_func(n4,k4,[self.Myc])+count_func(n4,k4,[self.Klf4]))+\
        #                         a6*(count_func(n4,k4,[self.EA2],mode='i')+count_func(n4,k4,[self.EA3],mode='i'))
        #
        # self.EA2.function_value=a7*(count_func(n4,k3,[self.Gata6])+count_func(n4,k4,[self.EA1],mode='i')+
        #                             count_func(n4,k4,[self.EA3],mode='i'))
        #
        # self.EA3.function_value=a7*(count_func(n4,k3,[self.Sox1])+count_func(n4,k4,[self.EA1],mode='i')+
        #                             count_func(n4,k4,[self.EA2],mode='i'))
        #
        #
        #
        #
        # self.Oct4_gene.value=a8*count_func(n4,k3,[self.Mbd3],mode='i')+count_func(n4,k3,[self.Nurd],mode='i')+\
        #                      count_func(n4,k4,[self.EA1])
        #
        # self.Sox2_gene.value=a8*count_func(n4,k3,[self.Mbd3],mode='i')+count_func(n4,k3,[self.Nurd],mode='i')+\
        #                      count_func(n4,k4,[self.EA1])
        #
        # self.Nanog_gene.value=a8*count_func(n4,k3,[self.Mbd3],mode='i')+count_func(n4,k3,[self.Nurd],mode='i')+\
        #                       count_func(n4,k4,[self.EA1],mode='i')
        #
        # self.Gata6_gene.value=a9*count_func(n4,k3,[self.Prc2],mode='i')+count_func(n4,k3,[self.Nurd],mode='i')+\
        #                       count_func(n4,k4,[self.EA2],mode='i')
        #
        # self.Sox1_gene.value=a9*count_func(n4,k3,[self.Prc2],mode='i')+count_func(n4,k3,[self.Mbd3],mode='i')+\
        #                       count_func(n4,k4,[self.EA3],mode='i')
        #
        for gene in self.gene_list:
            gene.function_value = round(gene.function_value, gene.max_acc)
            gene.ODE_result()
        for gene in self.DNA_gene_list:
            gene.value=round(gene.value,gene.max_acc)




    def reinit_gene_value(self):
        for gene in self.gene_list:
            gene.init_gene_value()

        self.Oct4_gene.init_gene_value()
        self.Nanog_gene.init_gene_value()
        self.Sox2_gene.init_gene_value()
        self.Gata6_gene.init_gene_value()
        self.Sox1_gene.init_gene_value()
        # print(len(self.gene_list))

    def trajectories_op(self, max_temperature, R):
        '''
        # TODO how to define n and R
        cost count in the same function
        cost function need to be improved
        :param max_temperature:
        :param R:
        :return:
        '''
        n=max_temperature*(2^8)
        self.steps=int(n/R)
        data=[]
        gene_dic = {}
        gene_name_dic = {}
        # TODO whether need a function
        # init
        for each_gene in self.gene_list:
            gene_dic[each_gene.name] = []
            gene_name_dic[each_gene.name] = each_gene
        for each_gene in self.DNA_gene_list:
            gene_dic[each_gene.name] = []
            gene_name_dic[each_gene.name] = each_gene


        for times in range(self.steps):
            for key in gene_dic:
                # dataset save
                gene_dic[key].append(gene_name_dic[key].value)
            self.one_tune()

        data.append(gene_dic)
        return data[0]


class DataSaver():
    def __init__(self,gene_data):
        self.gene_data=gene_data

    def release_data(self,gene_data):
        for gene in gene_data.gene_name_dic:
            gene_data.gene_name_dic[gene].link_value=self.gene_data.gene_name_dic[gene].link_value
            gene_data.gene_name_dic[gene].act_value=self.gene_data.gene_name_dic[gene].act_value
            gene_data.gene_name_dic[gene].inact_value=self.gene_data.gene_name_dic[gene].inact_value
        return gene_data

    def load_data(self,load_path):
        with open(load_path,'rb') as f:
            self.gene_data=pickle.load(f)

    def save_data(self,save_path):
        with open(save_path,'wb') as f:
            pickle.dump(self.gene_data,f)






#
# gene_data=GeneData()
# data_saver=DataSaver(gene_data)
#
# gene_data.reinit_gene_value()
# print(gene_data.trajectories_op(100, 1, 128))
# del gene_data
# gene_data=GeneData()
# data_saver.release_data(gene_data)
# # gene_data.attractor_count(10)
# for gene in gene_data.gene_list:
#     print(gene.name,gene.value)
# gene_data.one_tune()
# print(gene_data.Oct4.value)
# print(gene_data.Oct4.value, gene_data.Nanog.value)

