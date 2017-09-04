import algorithms

def min_e(G, m=0, return_index=False):

    x = G[:,0]
    y = G[:,1]

    ix_sel = algorithms.min_e_approximation(x,y,m)

    if(return_index):
        return ix_sel
    else:
        return G[ix_sel,:]

def min_num(G, epsilon=0.0, return_index=False):

    x = G[:,0]
    y = G[:,1]

    ix_sel = algorithms.min_num_approximation(x,y,epsilon)

    if(return_index):
        return ix_sel
    else:
        return G[ix_sel,:]
