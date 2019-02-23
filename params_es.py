from data import GeneData
max_iterations=100
max_temperature=100

# result_gene_data=GeneData()
# for _ in range(max_iterations):
#     gene_data = GeneData()
#     data=gene_data.trajectories_op(max_temperature, 1, 128)

def attractor_count(max_att_iter):
    # init
    iter_count=0
    gene_data=GeneData()
    attractor_value = gene_data._attractor_define()

    while iter_count<max_att_iter:
        data=gene_data.trajectories_op(100,1,128)[0]
        gap=gene_data._sum_value_for_data(data,gene_data.steps)-\
            gene_data._sum_value_for_data(data,gene_data.steps-2)
        # for a stable state
        if gap>gene_data.gene_num*gene_data.max_acc or gap==float('nan'):
            # attractor value
            tmp_attractor_value = gene_data._attractor_define()
            for i in range(len(attractor_value)):
                attractor_value[i]+=tmp_attractor_value[i]
            iter_count += 1
            print(iter_count)
        del gene_data
        gene_data=GeneData()
    print(attractor_value)

attractor_count(10)
