# -*- coding: utf-8 -*-


"""
================================
Boolean Dynamics Database Models
================================

:Author:
    Moritz Emanuel Beber
:Date:
    2014-11-19
:Copyright:
    Copyright |c| 2014, Jacobs University Bremen gGmbH, all rights reserved.
:File:
    models.py

.. |c| unicode:: U+A9
"""


from __future__ import (absolute_import, unicode_literals)


__all__ = ["Base", "Session", "Job", "ControlResult"]


import logging

from sqlalchemy import (Column, ForeignKey, Integer, String, Sequence,
        Float, Boolean, select)
from sqlalchemy.orm import (sessionmaker, relationship)
from sqlalchemy.ext.declarative import declarative_base
from pandas import DataFrame

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


Base = declarative_base()
Session = sessionmaker()


class Job(Base):
    __tablename__ = "job"
    id = Column(Integer, Sequence("job_id_seq"), primary_key=True)
    generator = Column(String()) # e.g., simple
    num_nodes = Column(Integer)
    num_activating = Column(Integer)
    num_inhibiting = Column(Integer)
    num_repeats = Column(Integer)
    steps = Column(Integer)
    win_size = Column(Integer)
    random_num = Column(Integer) # often 1E04
    function = Column(String()) # e.g., regulatory
    seed = Column(Integer)
    complete = Column(Boolean, default=False, nullable=False)
    results = relationship("ControlResult", backref="job")


class ControlResult(Base):
    __tablename__ = "controlresult"
    id = Column(Integer, Sequence("controlresult_id_seq"), primary_key=True)
    job_id = Column(Integer, ForeignKey("job.id"))
    ctc = Column(Float)
    measure = Column(String(30)) # describe the point in a series
    delay = Column(Integer) # usually between 1 and 6 with delay measures

    @classmethod
    def load_frame(cls, session):
        """
        Load part of the table into a well-formatted pandas.DataFrame.

        session can be any object with the execute method.
        """
        result = cls.__table__
        stmt = select([result.c.id, result.c.ctc,
                result.c.measure, result.c.delay])
        result = session.execute(stmt)
        df = DataFrame(iter(result), columns=result.keys())
        df.set_index("id", inplace=True)
        return df

