# -*- coding: utf-8 -*-


"""
========================
Boolean Dynamics Classes
========================

:Author:
    Moritz Emanuel Beber
:Date:
    2014-09-27
:Copyright:
    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
:File:
    classes.py

.. |c| unicode:: U+A9
"""


from __future__ import (absolute_import, unicode_literals)


__all__ = ["BooleanDynamicsError", "BooleanDynamics"]


import logging

import numpy as np

from . import _dynamics as dyn


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


class BooleanDynamicsError(Exception):
    pass


class BooleanDynamics(object):
    """
    """

    def __init__(self, net, node2id=None, function="regulatory"):
        """
        """
        self.num_nodes = None
        self.nodes = None
        self.node2id = None
        self.incidence_ptr = None
        self.incidence_adj = None
        self.incidence_func = None
        self.from_trn(net, node2id, function)

    def from_trn(self, trn, node2id=None, function="regulatory"):
        """
        Prepare data structures for boolean dynamics using a transcriptional
        regulatory network.

        Parameters
        ----------
        trn: nx.(Multi-)DiGraph
            The transcriptional regulatory network (TRN).
        node2id: dict (optional)
            A mapping from nodes in the network to indices running from 0 to (N - 1).
        function: hashable (optional)
            In case of a DiGraph, the keyword for edge data that describes the
            regulatory function. For a MultiDiGraph the key is used.
        """
        if len(trn) < 2 or trn.size() < 2:
            raise BooleanDynamicsError("aborting due to small network size")
        self.nodes = sorted(trn.nodes())
        self.num_nodes = len(self.nodes)
        if node2id is None:
            self.node2id = {n: i for (i, n) in enumerate(self.nodes)}
        else:
            self.node2id = node2id
        self.incidence_ptr = np.zeros(self.num_nodes + 1, dtype=np.int32)
        self.incidence_adj = np.zeros(trn.size(), dtype=np.int32)
        self.incidence_func = np.zeros(trn.size(), dtype=np.int32)
        if trn.is_multigraph():
            for n in self.nodes:
                i = self.node2id[n]
                i_ptr = self.incidence_ptr[i]
                in_edges = trn.in_edges(n, keys=True)
                self.incidence_ptr[i + 1] = i_ptr + len(in_edges)
                for (j_ptr, (u, v, k)) in enumerate(in_edges, start=i_ptr):
                    j = self.node2id[u]
                    self.incidence_adj[j_ptr] = j
                    self.incidence_func[j_ptr] = k
        else:
            for n in self.nodes:
                i = self.node2id[n]
                i_ptr = self.incidence_ptr[i]
                in_edges = trn.in_edges(n, data=True)
                self.incidence_ptr[i + 1] = i_ptr + len(in_edges)
                for (j_ptr, (u, v, data)) in enumerate(in_edges, start=i_ptr):
                    j = self.node2id[u]
                    self.incidence_adj[j_ptr] = j
                    self.incidence_func[j_ptr] = data[function]

    def time_series(self, steps, states=None, asynchronous=False, seed=None):
        """
        Follow the ON/OFF states of all nodes for a number of time steps.

        Parameters
        ----------
        steps: int
            Number of steps in the time series.
        states: bool (optional)
            Provide initial ON/OFF states for each node.
        asynchronous: bool (optional)
            Switch between a synchronous (default) and an asynchronous state update.
        seed: hashable (optional)
            Set the state of the seed for the random number generator.

        Returns
        -------
        A two dimensional array, where the first dimension corresponds to the
        nodes and the second to the time points. The number of time points is
        `steps` + 1, i.e., initial state and all time steps.
        """
        steps = steps + 1
        if seed is not None:
            np.random.seed(seed)
        if states is None:
            states = np.random.random_integers(0, 1,
                    self.num_nodes).astype(np.ubyte)
        else:
            states = np.asarray(states, dtype=np.ubyte)
        series = np.zeros(self.num_nodes * steps, dtype=np.ubyte)
        dyn.time_series(states, self.incidence_adj, self.incidence_func,
                self.incidence_ptr, self.num_nodes, series, steps)
        return series.reshape((steps, self.num_nodes)).T

