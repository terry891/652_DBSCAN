import numpy as np



def one_data_point(vars, coefs):
    return np.dot(vars, coefs)  

# Generate n random pointss
n = 10
x1 = np.random.rand(n)
x2 = np.random.rand(n)
x3 = np.random.rand(n)

# Compute y values using the function
y = f(x1, x2, x3)

# Print the points
for i in range(n):
    print("Point {}: ({}, {}, {}) -> y = {}".format(i+1, x1[i], x2[i], x3[i], y[i]))
