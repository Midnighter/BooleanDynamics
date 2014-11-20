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
#    cdynamics.pxd
#
#.. |c| unicode:: U+A9
#"""


cdef extern from "dynamics.h":
    void attractor_finder(unsigned char *states, int *inc_adj, int *inc_reg, int *inc_ptr,
            const int num_nodes, size_t *attractor_id, const size_t num_states)
    void time_series(unsigned char *states, int *inc_adj, int *inc_reg, int *inc_ptr,
            const int num_nodes, unsigned char *series, const size_t time)

