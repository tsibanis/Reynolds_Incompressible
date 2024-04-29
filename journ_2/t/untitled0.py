
import numpy as np

data = np.loadtxt('t_1000_1000.dat')
column1 = data[:, 1]
column1_str = ','.join(str(x) for x in column1)

print(column1_str)