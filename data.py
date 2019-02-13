from Gene.gene import gene
from Gene.DNAgene import DNAgene
from Gene.merge_gene import merge_gene

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

gene_data=GeneData()
print(gene_data.Oct4.init_value,gene_data.Nanog.init_value)