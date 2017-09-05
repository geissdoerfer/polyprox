#!/usr/bin/env python
"""Small demo for visualization of algorithms

Both approximations (by max error/ by number of sampling points) are applied
to an exemplary trajectory implemented as random walk.
"""

import polyprox
import numpy as np
import rdp
import time
import matplotlib.pyplot as plt

__author__ = "Kai Geissdoerfer"
__copyright__ = "Copyright (c) 2017, CSIRO"
__credits__ = ["Kun Zhao"]
__license__ = "CSIRO Open Source Software License Agreement"
__version__ = "0.1"
__maintainer__ = "Kai Geissdoerfer"
__email__ = "kai.geissdoerfer@mailbox.tu-berlin.de"
__status__ = "Production"

## Random walk
N = 2000

phi = np.random.uniform(0.0,np.pi,N)
s = np.random.uniform(1.0,5.0,N)

x = np.cumsum(s*np.cos(phi))
y = np.cumsum(s*np.sin(phi))

G = np.array((x,y)).transpose()

## Min-num problem
# Set maximum allowable error
epsilon = 20.0

# Time execution of min_num method provided in this package
t_start = time.clock()
G_pp =  polyprox.min_num(G,epsilon)
t_exec_pp = time.clock()-t_start

# Time execution of Hirschmann's implementation
t_start = time.clock()
G_rdp =  rdp.rdp(G,epsilon)
t_exec_rdp = time.clock()-t_start

print("Same result as rdp: {}".format(np.array_equal(G_rdp,G_pp)))
print("Speedup: {:.2f}%".format((t_exec_rdp-t_exec_pp)/t_exec_rdp*100.0))

plt.plot(G[:,0], G[:,1], label= 'Groundtruth')
plt.plot(G_pp[:,0], G_pp[:,1], 'g--o', label = 'Approximation')
plt.legend()
plt.show()

## Min-e problem
# Number of points by which to approximate
m = 10
G_pp =  polyprox.min_e(G,m)
t_exec_pp = time.clock()-t_start

print("Requested: {} Got: {}".format(m, len(G_pp)-2))

plt.plot(G[:,0], G[:,1], label= 'Groundtruth')
plt.plot(G_pp[:,0], G_pp[:,1], 'r--o', label = 'Approximation')
plt.legend()
plt.show()
