def act(gene_value,n,k):
    return gene_value**n/(k**n+gene_value**n)

def in_act(gene_value,n,k):
    return k**n/k**n+gene_value**n