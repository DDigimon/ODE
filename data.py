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
        self.Oct4.activate_link.append(cp([self.OS,self.Klf4,self.OSN,self.Myc]))
        self.Oct4.inactivate_link.append(cp([self.Gata6]))
        self.Oct4.inactivate_link.append(cp([self.Sox1]))

        self.Sox2.activate_link.append(cp([self.OS,self.OSN,self.Myc,self.Klf4]))



gene_data=GeneData()
gene_data.add_link()
print(gene_data.Oct4.init_value,gene_data.Nanog.init_value)
