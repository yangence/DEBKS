#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from debks.version import __version__

setup(name='DEBKS',
      version=__version__,
      description='Analysis differential back-spliced in of ribo-zero RNA-seq',
      author='Zelin Liu',
      author_email='zlliu@bjmu.edu.cn',
      url='https://github.com/yangence/DEBKS',
      license='GPL3',
      keywords='circular RNAs',
      python_requires=">=3",
      packages=find_packages(),
      install_requires=[
          'scipy>=1.2.1',
          'pysam>=0.15.2',
          'pandas>=0.24.2',
          'numpy>=1.16.6',
          'docopt>=0.6.2'
      ],
      entry_points={
          'console_scripts': [
              'DEBKS=debks.DEBKS_main:main'
          ],
      },
      )
