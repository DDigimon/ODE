class merge_gene():
    def __init__(self,gene_list):
        self.value=1
        for i in gene_list:
            self.value*=i.value