#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import re

from setuptools import setup, find_packages


classifiers = """\
    Development Status :: 4 - Beta
    Operating System :: OS Independent
    Programming Language :: Python
    Programming Language :: Python :: 2
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.4
    Programming Language :: Python :: 3.5
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: 3.7
"""


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_version():
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read('cooler', '_version.py'),
        re.MULTILINE).group(1)
    return version


def get_long_description():
    return _read('README.md')


install_requires = [
    'six',
    'numpy>=1.9',
    'scipy>=0.16',
    'pandas>=0.19',
    'h5py>=2.5',
    'click>=6.6',
    'pysam>0.8',
    'cytoolz',
    'multiprocess',
    'biopython',
    'pyfaidx',
    'pypairix',
    'asciitree',
    'pyyaml'
]


tests_require = [
    'pytest',
    'mock'
]


extras_require = {
    'docs': [
        'Sphinx>=1.1',
        'numpydoc>=0.5'
    ]
}


setup(
    name='cooler',
    author='Nezar Abdennur',
    author_email='nezar@mit.edu',
    version=get_version(),
    license='BSD',
    description='Sparse binary format for genomic interaction matrices',
    long_description=get_long_description(),
    long_description_content_type='text/markdown',
    keywords=['genomics', 'bioinformatics', 'Hi-C', 'contact', 'matrix', 'format', 'hdf5'],
    url='https://github.com/mirnylab/cooler',
    packages=find_packages(),
    zip_safe=False,
    classifiers=[s.strip() for s in classifiers.split('\n') if s],
    install_requires=install_requires,
    tests_require=tests_require,
    extras_require=extras_require,
    entry_points={
        'console_scripts': [
            'cooler = cooler.cli:cli',
        ]
    }
)
