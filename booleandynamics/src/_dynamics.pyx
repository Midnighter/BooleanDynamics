# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False


#"""
#===============
#Boolean Updates
#===============
#
#:Author:
#    Moritz Emanuel Beber
#:Date:
#    2014-09-27
#:Copyright:
#    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
#:File:
#    _dynamics.pyx
#
#.. |c| unicode:: U+A9
#"""


cimport cdynamics as dyn


ctypedef unsigned char UChar


def time_series(UChar[:] states, int[:] inc_adj, int[:] inc_reg, int[:] inc_ptr,
        int num_nodes, UChar[:] series, int steps):
    # randomize initial state
    dyn.time_series(&states[0], &inc_adj[0], &inc_reg[0], &inc_ptr[0], num_nodes,
            &series[0], steps)
    return series

