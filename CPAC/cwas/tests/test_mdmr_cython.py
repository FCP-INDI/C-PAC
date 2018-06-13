import os

def test_mdmr():

    from CPAC.cwas.mdmr import distance, mdmr
    import numpy as np

    X = np.genfromtxt(os.path.join(os.path.dirname(__file__), 'X.csv'), delimiter=',')
    Y = np.genfromtxt(os.path.join(os.path.dirname(__file__), 'Y.csv'), delimiter=',')

    X = np.random.uniform(size=(5, 10))
    D = mdmr(X, 20)
    
    print
    for i in range(D.n):
        print
        for j in range(D.n):
            print D[i, j],