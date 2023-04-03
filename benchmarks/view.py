import h5py


filename = 'iris.h5'

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
    print(data)




