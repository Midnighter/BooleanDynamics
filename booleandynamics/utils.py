# -*- coding: utf-8 -*-


"""
========================
Boolean Dynamics Utility
========================

:Author:
    Moritz Emanuel Beber
:Date:
    2014-10-30
:Copyright:
    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
:File:
    utils.py

.. |c| unicode:: U+A9
"""


from __future__ import (absolute_import, unicode_literals)


__all__ = ["to_expression", "stitch_series"]


import logging

import numpy as np
from numpy.lib.stride_tricks import as_strided


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def contiguous_blocks(matrix, size, overlap=0):
    if overlap < 0:
        raise ValueError("overlap must be zero or positive")
    if overlap >= size:
        raise ValueError("overlap between windows cannot be equal to or"\
                " greater than the windows themselves")
    step = size - overlap
    num_cols = matrix.shape[1] // step
    num_rows = matrix.shape[0]
    shape = (num_cols, num_rows, size)
    strides = ((step - overlap) * matrix.strides[-1],) + matrix.strides
    return as_strided(matrix, shape, strides)

def to_expression(series, size, norm="window"):
    """
    Compute the activity per gene within a window normalized by the total
    activity in that window.

    Parameters
    ----------
    series: numpy.ndarray
        Boolean array of dimension (number of nodes x number of time steps).
    size: int
        Size of the windows over which to compute the activity.
    norm: str
        Normalize by activity per window or by total activity ('window',
        'total').

    Returns
    -------
    """
    num_cols = series.shape[1] // size
    num_rows = series.shape[0]
    expression = np.zeros((num_rows, num_cols))
    if norm == "total":
        for (i, block) in enumerate(contiguous_blocks(series, size)):
            expression[:, i] = block.sum(axis=1).astype(np.double)
        expression /= series.sum().astype(np.double)
    elif norm == "window":
        for (i, block) in enumerate(contiguous_blocks(series, size)):
            expression[:, i] = block.sum(axis=1).astype(np.double) / block.sum(
                    ).astype(np.double)
    else:
        raise ValueError("unknown normalization method '{}'".format(norm))
    return expression

def stitch_series(rbn, repeat, steps, seed=None):
    """
    Stitch together multiple simulated node activity series for different
    starting conditions.

    Parameters
    ----------
    rbn: BooleanDynamics
        An instance of BooleanDynamics with prepared topology.
    repeat: int
        How many different starting conditions to put together.
    steps: int
        The number of steps of activity to follow per run.

    Returns
    -------
    2D array where the first dimension corresponds to the nodes and the second
    to ``repeat * steps``.
    """
    if seed is not None:
        np.random.seed(seed)
    long_series = np.zeros((rbn.num_nodes, repeat * steps), dtype=np.ubyte)
    known_states = set()
    for i in range(repeat):
        states = tuple(np.random.random_integers(0, 1, rbn.num_nodes))
        while states in known_states:
            states = tuple(np.random.random_integers(0, 1, rbn.num_nodes))
        t = i * steps
        long_series[:, t: (t + steps)] = rbn.time_series(steps - 1, states=states)
    return long_series

