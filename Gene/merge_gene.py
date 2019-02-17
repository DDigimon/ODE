class merge_gene():
    def __init__(self,name,gene_list):
        self.value=1
        self.name=name
        for i in gene_list:
            self.value*=i.value