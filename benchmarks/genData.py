import numpy as np
import h5py

nPoints = [100, 100, 100]
nDimensions = [200, 500, 1000]
nClusters = [10, 10, 10]

for (ps, ds, k) in zip(nPoints, nDimensions, nClusters):
    filename = f"tests/test_{ps*k}p_{ds}d_{k}k.h5"



    means = 200 * np.random.rand(k, ds) - 100
    x=[]
    for k in range(k):
        x.append( means[k] +  np.random.randn(ps, ds) )
    x = np.vstack(x)
    np.random.shuffle(x)

    with h5py.File(filename, "w") as data_file:
        data_file.create_dataset("DBSCAN", data=x, dtype='float32')

    print(x.shape)


