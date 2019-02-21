from data import GeneData
max_iterations=100
max_temperature=100

gene_data=GeneData()
gene_data.add_link()
for _ in range(max_iterations):
    gene_data.one_tune()