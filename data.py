from Gene.gene import gene
from Gene.DNAgene import DNAgene
from Gene.merge_gene import merge_gene
from Gene.combination import combination as cp

class GeneData():
    def __init__(self):
        self.Oct4=gene()
        self.Sox2=gene()
        self.Nanog=gene()
        self.Cdx2=gene()
        self.Sox1=gene()
        self.Gcnf=gene()
        self.Pax6=gene()
        self.Gata6=gene()
        self.Myc=gene()
        self.Klf4=gene()
        self.Prc2=gene()
        self.Mbd3=gene()
        self.Nurd=gene()
        self.EA1=gene()
        self.EA2=gene()
        self.EA3=gene()

        self.Oct4_gene=DNAgene()
        self.Sox2_gene=DNAgene()
        self.Nanog_gene=DNAgene()
        self.Gata6_gene=DNAgene()
        self.Sox1_gene=DNAgene()

        self.OC=merge_gene()
        self.OS=merge_gene()
        self.OSN=merge_gene()
        self.OG=merge_gene()

    def add_link(self):
        # TODO multiple between link combination, while add in the combination
        # self.Oct4.activate_link.append(cp([self.OS,self.Klf4,self.OSN,self.Myc]))
        # self.Oct4.inactivate_link.append(cp([self.Gata6]))
        # self.Oct4.inactivate_link.append(cp([self.Sox1]))
        #
        # self.Sox2.activate_link.append(cp([self.OS,self.OSN,self.Myc,self.Klf4]))

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





gene_data=GeneData()
gene_data.add_link()
print(gene_data.Oct4.init_value,gene_data.Nanog.init_value)
