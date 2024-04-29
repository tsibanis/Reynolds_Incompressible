import numpy as np

# Load data from file
data = np.loadtxt("filename.txt")

# Extract first column
col1 = data[:,0]

# Print as array
print(col1)