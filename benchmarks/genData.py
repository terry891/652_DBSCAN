import numpy as np
import h5py

numPts = 20
numDims = 3
numK = 4
filename = "test1.h5"

means = 200 * np.random.rand(numK, numDims) - 100
x=[]
for k in range(numK):
	x.append( means[k] +  np.random.randn(numPts, numDims) )
x = np.vstack(x)
np.random.shuffle(x)

with h5py.File(filename, "w") as data_file:
    data_file.create_dataset("DBSCAN", data=x, dtype='float32')

print(x.shape)


