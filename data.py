from Gene.gene import gene
from Gene.DNAgene import DNAgene
from Gene.merge_gene import merge_gene


class GeneData():
    def __init__(self):
        self.Oct4_gene=DNAgene()
        self.Sox2_gene=DNAgene()
        self.Nanog_gene=DNAgene()
        self.Gata6_gene=DNAgene()
        self.Sox1_gene=DNAgene()

        self.Oct4=gene('Oct4',DNAgene=self.Oct4_gene)
        self.Sox2=gene('Sox2')
        self.Nanog=gene('Nanog')
        self.Cdx2=gene('Cdx2')
        self.Sox1=gene('Sox1')
        self.Gcnf=gene('Gcnf')
        self.Pax6=gene('Pax6')
        self.Gata6=gene('Gata6')
        self.Myc=gene('Myc')
        self.Klf4=gene('Klf4')
        self.Prc2=gene('Prc2')
        self.Mbd3=gene('Mbd3')
        self.Nurd=gene('Nurd')
        self.EA1=gene('EA1')
        self.EA2=gene('EA2')
        self.EA3=gene('EA3')


        self.OC=merge_gene('OC',[self.Oct4,self.Cdx2])
        self.OS=merge_gene('OS',[self.Oct4,self.Sox2])
        self.OSN=merge_gene('OSN',[self.Oct4,self.Sox2,self.Nanog])
        self.OG=merge_gene('OG',[self.Oct4,self.Gata6])

    def add_link(self):
        self.Oct4.activate_link=[self.OS,self.Klf4,self.OSN,self.Myc]
        self.Oct4.inactivate_link=[self.Gata6,self.Sox1,self.OC,self.Gcnf]
        self.Oct4.init_param()

        self.Sox2.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Sox2.inactivate_link=[self.Sox1,self.Gata6,self.Oct4]
        self.Sox2.init_param()

        self.Nanog.activate_link=[self.OS,self.OSN,self.Myc,self.Klf4]
        self.Nanog.inactivate_link=[self.OG,self.Sox1,self.Oct4]
        self.Nanog.init_param()

        self.Cdx2.activate_link=[self.Cdx2]
        self.Cdx2.inactivate_link=[self.Nanog,self.OC]
        self.Cdx2.init_param()

        self.Gcnf.activate_link=[self.Cdx2,self.Gata6,self.Gcnf]
        self.Gcnf.init_param()

        self.Pax6.activate_link=[self.Sox2,self.Pax6]
        self.Pax6.inactivate_link=[self.Oct4,self.Nanog]
        self.Pax6.init_param()

        self.Sox1.activate_link=[self.Sox2]
        self.Sox1.inactivate_link=[self.Gata6,self.Oct4,self.Nanog]
        self.Sox1.init_param()

        self.Gata6.activate_link=[self.Oct4]
        self.Gata6.inactivate_link=[self.Sox1,self.Sox2,self.Nanog]
        self.Gata6.init_param()

        self.Myc.activate_link=[self.Klf4]
        self.Myc.inactivate_link=[self.Sox1,self.Gata6,self.Cdx2,self.Gcnf]
        self.Myc.init_param()

        self.Klf4.activate_link=[self.Klf4]
        self.Klf4.inactivate_link=[self.Sox1,self.Cdx2,self.Gata6,self.Gcnf]
        self.Klf4.init_param()

        self.Prc2.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.Prc2.inactivate_link=[self.Mbd3,self.Nurd]
        self.Prc2.init_param()

        self.Mbd3.activate_link=[self.Gata6]
        self.Mbd3.inactivate_link=[self.Prc2,self.Nurd]
        self.Mbd3.init_param()

        self.Nurd.activate_link=[self.Sox1]
        self.Nurd.inactivate_link=[self.Prc2,self.Mbd3]
        self.Nurd.init_param()

        self.EA1.activate_link=[self.Oct4,self.Sox2,self.Myc,self.Klf4]
        self.EA1.inactivate_link=[self.EA2,self.EA3]
        self.EA1.init_param()

        self.EA2.activate_link=[self.Gata6]
        self.EA2.inactivate_link=[self.EA1,self.EA3]
        self.EA2.init_param()

        self.EA3.activate_link=[self.Sox1]
        self.EA3.inactivate_link=[self.EA1,self.EA2]
        self.EA3.init_param()

    def one_tune(self):
        self.Oct4.one_tune_value()
        self.Oct4.function_value= self.Oct4.act_value*\
                                  self.Oct4.or_op([self.OS, self.Klf4, self.OSN, self.Myc]) * \
                                  self.Oct4.and_op([self.Gata6,self.Sox1]) + \
                                  self.Oct4.inact_value*self.Oct4.or_op([self.OC,self.Gcnf])
        self.Oct4.ODE_result()

        self.Sox2.one_tune_value()
        self.Sox2.ODE_result()









gene_data=GeneData()
gene_data.add_link()
print(gene_data.Oct4.value)
gene_data.one_tune()
print(gene_data.Oct4.value)
print(gene_data.Oct4.value, gene_data.Nanog.value)
