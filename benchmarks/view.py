import h5py
import matplotlib.pyplot as plt
import numpy as np





filename = 'tests/iris.h5'



with h5py.File(filename, 'r') as f:
    print("Keys: %s" % list(f.keys())[0])
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    print(type(f[a_group_key])) 

    # If a_group_key is a group name, 
    # this gets the object names in the group and returns as a list
    data = list(f[a_group_key])

    # preferred methods to get dataset values:
    ds_obj = f[a_group_key]      # returns as a h5py dataset object
    ds_arr = f[a_group_key][()]  # returns as a numpy array

    # get the shape of the dataset
    print(ds_obj.shape)
    # print(data)



    # Plot 
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []
    for i in data:
        xs.append(i[0])
        ys.append(i[1])
        zs.append(i[2])
    
    ax.scatter(xs, ys, zs, marker='o')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()




