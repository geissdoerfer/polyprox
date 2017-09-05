#!/usr/bin/env python
"""Provides python interface to underlying C algorithms"""

import algorithms

__author__ = "Kai Geissdoerfer"
__copyright__ = "Copyright (c) 2017, CSIRO"
__credits__ = ["Kun Zhao"]
__license__ = "CSIRO Open Source Software License Agreement"
__version__ = "0.1"
__maintainer__ = "Kai Geissdoerfer"
__email__ = "kai.geissdoerfer@mailbox.tu-berlin.de"
__status__ = "Production"

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
