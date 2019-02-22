from data import GeneData
max_iterations=100
max_temperature=100

result_gene_data=GeneData()
for _ in range(max_iterations):
    gene_data = GeneData()
    data=gene_data.trajectories_op(max_temperature, 1, 128)
