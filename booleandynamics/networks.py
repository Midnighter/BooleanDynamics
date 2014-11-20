# -*- coding: utf-8 -*-


"""
=========================
Boolean Dynamics Networks
=========================

:Author:
    Moritz Emanuel Beber
:Date:
    2014-11-19
:Copyright:
    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
:File:
    networks.py

.. |c| unicode:: U+A9
"""


from __future__ import (absolute_import, unicode_literals)


__all__ = ["random_regulatory", "test_feedback"]


import logging
import random

import networkx as nx


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def add_self_inhibition(net, function="function"):
    if net.is_multigraph():
        for node in net:
            if not any(triple[2] == -1 for triple in\
                    net.in_edges_iter(node, keys=True)):
                net.add_edge(node, node, key=-1)
    else:
        for node in net:
            if not any(data[function] == -1 for (_, _, data) in\
                    net.in_edges_iter(node, data=True)):
                net.add_edge(node, node, **{function: -1})

def random_regulatory(num_nodes, num_activating, num_inhibiting,
        function="function", seed=None):
    net = nx.gnm_random_graph(num_nodes, num_activating + num_inhibiting,
            seed=seed, directed=True)
    edges = net.edges()
    activating = random.sample(edges, num_activating)
    for (u, v) in activating:
        net[u][v][function] = 1
    for (u, v) in list(set(edges) - set(activating)):
        net[u][v][function] = -1
    add_self_inhibition(net, function)
    return net

def test_feedback():
    net = nx.DiGraph()
    net.add_edge("A", "B", function=1)
    net.add_edge("B", "C", function=1)
    net.add_edge("C", "A", function=1)
    add_self_inhibition(net)
    return net

# TODO: adapt to this type of network
def draw(self, filename, output_format="pdf", layout_program="fdp",
            layout_args="", distinct=False):
    import pygraphviz as pgv
    OPTIONS = misc.OptionsManager.get_instance()
    net = pgv.AGraph(directed=True, name=filename, strict=True)
    node_attr= dict()
    link_attr= dict()
    # add compound nodes
    indeces = dict(itertools.izip(self.compounds, itertools.count()))
    for (cmpd, i) in indeces.iteritems():
        net.add_node(i, label=str(cmpd), shape="ellipse", **node_attr)
    # add reactions
    indeces.update(itertools.izip(self.reactions,
            itertools.count(len(self.compounds))))
    i = len(self.compounds) + len(self.reactions)
    for rxn in self.reactions:
        net.add_node(indeces[rxn], label=str(rxn), shape="box", **node_attr)
        # add forward reaction links
        for cmpd in self.predecessors(rxn):
            net.add_edge(indeces[cmpd], indeces[rxn], **link_attr)
        for cmpd in self.successors(rxn):
            net.add_edge(indeces[rxn], indeces[cmpd], **link_attr)
        if rxn.reversible:
            if distinct:
                rev = pymet.BasicReaction(str(rxn) + OPTIONS.reversible_suffix)
                indeces[rev] = i
                net.add_node(i, label=str(rev), shape="box", **node_attr)
                # add backward reaction links
                for cmpd in self.predecessors(rxn):
                    net.add_edge(indeces[rev], indeces[cmpd], **link_attr)
                for cmpd in self.successors(rxn):
                    net.add_edge(indeces[cmpd], indeces[rev], **link_attr)
                i += 1
            else:
                # add backward reaction links
                for cmpd in self.predecessors(rxn):
                    net.add_edge(indeces[rxn], indeces[cmpd],
                            style="dotted", **link_attr)
                for cmpd in self.successors(rxn):
                    net.add_edge(indeces[cmpd], indeces[rxn],
                            style="dotted", **link_attr)
    filename = "%s.%s" % (filename, output_format)
    net.draw(filename, prog=layout_program, args=layout_args)

