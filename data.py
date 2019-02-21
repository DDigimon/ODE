from Gene.gene import gene
from Gene.merge_gene import merge_gene


class GeneData():
    def __init__(self):
        # TODO use special method to load data, some mode?
        self.gene_list=[]
        self.DNA_gene_list=[]
        self.gene_dic={}
        self.gene_name_dic={}

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


        self.OC=merge_gene('OC',[self.Oct4,self.Cdx2])
        self.OS=merge_gene('OS',[self.Oct4,self.Sox2])
        self.OSN=merge_gene('OSN',[self.Oct4,self.Sox2,self.Nanog])
        self.OG=merge_gene('OG',[self.Oct4,self.Gata6])

        # TODO whether need a function

        for each_gene in self.gene_list:
            self.gene_dic[each_gene.name]=[]
            self.gene_name_dic[each_gene.name]=each_gene
        for each_gene in self.DNA_gene_list:
            self.gene_dic[each_gene.name]=[]
            self.gene_name_dic[each_gene.name]=each_gene

    def add_link(self):
        self.Oct4.activate_link=[self.OS,self.Klf4,self.OSN,self.Myc]
        self.Oct4.inactivate_link=[self.Gata6,self.Sox1,self.OC,self.Gcnf]

        self.Sox2.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Sox2.inactivate_link=[self.Sox1,self.Gata6,self.Oct4]

        self.Nanog.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Nanog.inactivate_link=[self.OG,self.Sox1,self.Oct4]

        self.Cdx2.activate_link=[self.Cdx2]
        self.Cdx2.inactivate_link=[self.Nanog,self.OC]

        self.Gcnf.activate_link=[self.Cdx2,self.Gata6,self.Gcnf]

        self.Pax6.activate_link=[self.Sox2,self.Pax6]
        self.Pax6.inactivate_link=[self.Oct4,self.Nanog]

        self.Sox1.activate_link=[self.Sox2]
        self.Sox1.inactivate_link=[self.Gata6,self.Oct4,self.Nanog]

        self.Gata6.activate_link=[self.Oct4]
        self.Gata6.inactivate_link=[self.Sox1,self.Sox2,self.Nanog]
        self.Gata6.init_param()

        self.Myc.activate_link=[self.Klf4]
        self.Myc.inactivate_link=[self.Sox1,self.Gata6,self.Cdx2,self.Gcnf]

        self.Klf4.activate_link=[self.Klf4]
        self.Klf4.inactivate_link=[self.Sox1,self.Cdx2,self.Gata6,self.Gcnf]

        self.Prc2.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.Prc2.inactivate_link=[self.Mbd3,self.Nurd]

        self.Mbd3.activate_link=[self.Gata6]
        self.Mbd3.inactivate_link=[self.Prc2,self.Nurd]

        self.Nurd.activate_link=[self.Sox1]
        self.Nurd.inactivate_link=[self.Prc2,self.Mbd3]

        self.EA1.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.EA1.inactivate_link=[self.EA2,self.EA3]

        self.EA2.activate_link=[self.Gata6]
        self.EA2.inactivate_link=[self.EA1,self.EA3]

        self.EA3.activate_link=[self.Sox1]
        self.EA3.inactivate_link=[self.EA1,self.EA2]

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

    def one_tune(self):
        for gene in self.gene_list:
            gene.one_tune_value()

        for gene in self.DNA_gene_list:
            gene.one_tune_value()

        self.Oct4.function_value= self.Oct4.act_value*\
                                  self.Oct4.or_op([self.OS, self.Klf4, self.OSN, self.Myc]) * \
                                  self.Oct4.and_op([self.Gata6,self.Sox1]) + \
                                  self.Oct4.inact_value*self.Oct4.or_op([self.OC,self.Gcnf])

        self.Sox2.function_value=self.Sox2.act_value*\
                                 (self.Sox2.or_op([self.OS,self.OSN,self.Myc,self.Klf4]))*\
                                 self.Sox2.and_op([self.Sox1,self.Gata6])+\
                                 self.Sox2.inact_value*self.Sox2.and_op([self.Oct4])

        self.Nanog.function_value=self.Nanog.act_value*\
                                  self.Nanog.or_op([self.OS,self.OSN,self.Myc,self.Klf4])*\
                                  self.Nanog.and_op([self.OG,self.Sox1])+\
                                  self.Nanog.inact_value*self.Nanog.and_op([self.Oct4])

        self.Cdx2.function_value=self.Cdx2.act_value*\
                                 self.Cdx2.and_op([self.Cdx2,self.Nanog,self.OC])*self.Cdx2.inact_value

        self.Gcnf.function_value=self.Gcnf.act_value*self.Gcnf.or_op([self.Cdx2,self.Gata6,self.Gcnf])

        self.Pax6.function_value=self.Pax6.act_value*self.Pax6.or_op([self.Sox2,self.Pax6])*\
                                 self.Pax6.and_op([self.Oct4,self.Nanog])*self.Pax6.inact_value

        self.Sox1.function_value=self.Sox1.act_value*self.Sox1.and_op([self.Sox2,self.Gata6,self.Oct4])+\
                                 self.Sox1.and_op([self.Nanog])*self.Sox1.inact_value

        self.Gata6.function_value=self.Gata6.act_value*self.Gata6.and_op([self.Oct4,self.Sox1,self.Sox2])+\
                                  self.Gata6.inact_value*self.Gata6.and_op([self.Nanog])

        self.Myc.function_value=self.Myc.act_value*self.Myc.and_op([self.Klf4,self.Sox1,self.Gata6,self.Cdx2])+\
                                self.Myc.inact_value*self.Myc.and_op([self.Gcnf])

        self.Klf4.function_value=self.Klf4.act_value*self.Klf4.and_op([self.Klf4,self.Sox1,self.Cdx2,self.Gata6])+\
                                 self.Klf4.inact_value*self.Klf4.and_op([self.Gcnf])

        self.Prc2.function_value=self.Prc2.act_value*self.Prc2.or_op([self.Oct4,self.Sox2,self.Myc,self.Klf4])+\
                                 self.Prc2.inact_value*self.Prc2.and_op([self.Mbd3,self.Nurd])

        self.Mbd3.function_value=self.Mbd3.act_value*self.Mbd3.or_op([self.Gata6])+\
                                 self.Mbd3.inact_value*self.Mbd3.or_op([self.Prc2,self.Nurd])

        self.Nurd.function_value=self.Nurd.act_value*self.Nurd.or_op([self.Sox1])+\
                                 self.Nurd.inact_value*self.Nurd.or_op([self.Prc2,self.Mbd3])

        self.EA1.function_value=self.EA1.act_value*self.EA1.or_op([self.Oct4,self.Sox2,self.Myc,self.Klf4])+\
                                self.EA1.inact_value*self.EA1.or_op([self.EA2,self.EA3])

        self.EA2.function_value=self.EA2.act_value*self.EA2.or_op([self.Gata6])+\
                                self.EA2.inact_value*self.EA2.or_op([self.EA1,self.EA3])

        self.EA3.function_value=self.EA3.act_value*self.EA3.or_op([self.Sox1])+\
                                self.EA3.inact_value*self.EA3.or_op([self.EA1,self.EA2])

        for gene in self.gene_list:
            gene.ODE_result()

        self.Oct4_gene.value=self.Oct4_gene.act_value*self.Oct4_gene.or_op([self.EA1])+\
                             self.Oct4_gene.inact_value*self.Oct4_gene.or_op([self.Mbd3,self.Nurd])

        self.Sox2_gene.value=self.Sox2_gene.act_value*self.Sox2_gene.or_op([self.EA1])+\
                             self.Sox2_gene.inact_value*self.Sox2_gene.or_op([self.Mbd3,self.Nurd])

        self.Nanog_gene.value=self.Nanog_gene.act_value*self.Nanog_gene.or_op([self.EA1])+\
                              self.Nanog_gene.inact_value*self.Nanog_gene.or_op([self.Mbd3,self.Nurd])

        self.Gata6_gene.value=self.Gata6_gene.act_value*self.Gata6_gene.or_op([self.EA2])+\
                              self.Gata6_gene.inact_value*self.Gata6_gene.or_op([self.Prc2,self.Nurd])

        self.Sox1_gene.value=self.Sox1_gene.act_value*self.Sox1_gene.or_op([self.EA3])+\
                             self.Sox1_gene.inact_value*self.Sox1_gene.or_op([self.Prc2,self.Mbd3])



    def reinit_gene_value(self):
        for gene in self.gene_list:
            gene.init_gene_value()

        self.Oct4_gene.init_gene_value()
        self.Nanog_gene.init_gene_value()
        self.Sox2_gene.init_gene_value()
        self.Gata6_gene.init_gene_value()
        self.Sox1_gene.init_gene_value()

    def trajectories_dataset(self,max_temperature,max_iter,R):
        '''
        # TODO how to define n and R
        :param max_temperature:
        :param R:
        :return:
        '''
        n=max_temperature*(2^8)
        steps=int(n/R)
        print(steps)
        data=[]
        for id in range(max_iter):
            for _ in range(steps):
                for key in self.gene_dic:
                    self.gene_dic[key].append(self.gene_name_dic[key].value)
                self.one_tune()
            data.append(self.gene_dic)
        return data


#
# gene_data=GeneData()
# gene_data.add_link()
# print(gene_data.trajectories_dataset(100,1,128))
# gene_data.one_tune()
# print(gene_data.Oct4.value)
# print(gene_data.Oct4.value, gene_data.Nanog.value)

