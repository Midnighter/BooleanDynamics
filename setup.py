#!/usr/bin/env python
# -*- coding: utf-8 -*-


from os.path import join


if __name__ == '__main__':
    from setuptools import (setup, Extension, find_packages)
    from Cython.Distutils import build_ext
    from numpy import get_include

    # dynamics
    sources = ["_dynamics.pyx", "dynamics.c"]
    c_path = join("booleandynamics", "src")
    dynamics = Extension("booleandynamics._dynamics",
        sources=[join(c_path, src) for src in sources],
        include_dirs=[c_path, get_include()]
    )

    setup(
        name="booleandynamics",
        version="0.0.1",
        author="Moritz Emanuel Beber",
        license="BSD 3-Clause License",
        packages=find_packages(),
        ext_modules=[dynamics],
        cmdclass={"build_ext": build_ext}
    )

