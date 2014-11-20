#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
============================
Boolean Dynamics on Networks
============================

:Author:
    Moritz Emanuel Beber
:Date:
    2014-09-27
:Copyright:
    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
:File:
    __init__.py

.. |c| unicode:: U+A9
"""


from __future__ import (absolute_import, unicode_literals)

import sys
import logging
import multiprocessing
import argparse
from logging.config import dictConfig

import numpy as np
from sqlalchemy import create_engine
from progressbar import (ProgressBar, Timer, SimpleProgress, Bar, Percentage, ETA)

import booleandynamics as bd
import pyorganism.regulation as pyreg


LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def worker(job):
    generator = getattr(bd, job.generator)
    net = generator(job.num_nodes, job.num_activating, job.num_inhibiting,
            function=job.function, seed=job.seed)
    rbn = bd.BooleanDynamics(net, function=job.function)
    series = bd.stitch_series(rbn, job.num_repeats, job.steps, seed=job.seed)
    expr = bd.to_expression(series, job.win_size)
    control = pyreg.ContinuousControl()
    control.from_trn(net, function=job.function)
    trans_expr = np.ascontiguousarray(expr.T) # necessary for C code
    results = list()
    job_id = job.id
    (z_ca, _, _) = control.series_ctc(trans_expr, "absolute",
            random_num=job.random_num)
    # need to add delay: None, otherwise later values are ignored in insert
    results.extend({"job_id": job_id, "ctc": val, "delay": None,
            "measure": "absolute"} for val in z_ca)
    (z_cd, _, _) = control.series_ctc(trans_expr, "difference",
            random_num=job.random_num)
    results.extend({"job_id": job_id, "ctc": val, "delay": None,
            "measure": "difference"} for val in z_cd)
    (z_cf, _, _) = control.series_ctc(trans_expr, "functional",
            random_num=job.random_num)
    results.extend({"job_id": job_id, "ctc": val, "delay": None,
            "measure": "functional"} for val in z_cf)
    (z_cs, _, _) = control.series_ctc(trans_expr, "functional-comparison",
            random_num=job.random_num)
    results.extend({"job_id": job_id, "ctc": val, "delay": None,
            "measure": "functional-comparison"} for val in z_cs)
    (z_ct1, _, _) = control.series_ctc(trans_expr, "delayed-functional",
            random_num=job.random_num, delay=1)
    results.extend({"job_id": job_id, "ctc": val, "delay": 1,
            "measure": "delayed-functional"} for val in z_ct1)
    (z_ct2, _, _) = control.series_ctc(trans_expr, "delayed-functional",
            random_num=job.random_num, delay=2)
    results.extend({"job_id": job_id, "ctc": val, "delay": 2,
            "measure": "delayed-functional"} for val in z_ct2)
    (z_ct3, _, _) = control.series_ctc(trans_expr, "delayed-functional",
            random_num=job.random_num, delay=3)
    results.extend({"job_id": job_id, "ctc": val, "delay": 3,
            "measure": "delayed-functional"} for val in z_ct3)
    return (job_id, results)

def main(args):
    engine = create_engine(args.engine)
    bd.Base.metadata.bind = engine
    bd.Session.configure(bind=engine)
    session = bd.Session()
    tasks = session.query(bd.Job).filter(~bd.Job.complete).all()
    if len(tasks) == 0:
        LOGGER.warn("Nothing to do")
        return
    pool = multiprocessing.Pool(args.nproc)
    result_it = pool.imap_unordered(worker, tasks)
    bar = ProgressBar(maxval=len(tasks), widgets=[Timer(), " ",
            SimpleProgress(), " ", Percentage(), " ", Bar(), " ",
            ETA()]).start()
    for (job_id, result) in result_it:
        try:
            job = session.query(bd.Job).filter_by(id=job_id).one()
            # more low-level method for speed
            session.execute(bd.ControlResult.__table__.insert(), result)
            job.complete = True
            session.commit()
        except Exception:
            session.rollback()
        bar += 1
    bar.finish()
    session.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-v", "--version", action="version", version="0.1")
    parser.add_argument("--log-level", dest="log_level", default="INFO",
            help="Log level, i.e., DEBUG, INFO, WARN, ERROR, CRITICAL (default: %(default)s)")
    parser.add_argument("--encoding", dest="encoding", default="utf-8",
            help="File encoding to assume (default: %(default)s)")
    parser.add_argument("-n", "--nproc", dest="nproc",
            default=multiprocessing.cpu_count(), type=int,
            help="Number of processors to use (default: %(default)s)")
    parser.add_argument("engine",
            help="Database connection string, e.g., 'sqlite+pysqlite:///file.db'")
    args = parser.parse_args()
    logger = multiprocessing.log_to_stderr()
    dictConfig({"version": 1, "incremental": True,
            "root": {"level": args.log_level}})
    if args.log_level != "DEBUG":
        logging.getLogger(logger.name).setLevel(logging.WARN)
        logging.getLogger("pyorganism").setLevel(logging.WARN)
    try:
        sys.exit(main(args))
    except: # we want to catch everything
        (err, msg, trace) = sys.exc_info()
        # do something
        raise err, msg, trace
    finally:
        logging.shutdown()

