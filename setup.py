from setuptools import setup
import os

version_loc = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           'madic/version.py'))
exec(open(version_loc).read())  # supplies __version__

long_description = """

Overview
========
madic is a Python package for detecting and correcting interference in targeted
mass spectrometry data using high-performance Pandas data structures.

Purpose
=======
Multiple reaction monitoring (MRM, also known as selected reaction monitoring /
SRM) data acquired during targeted mass spectrometry experiments can suffer
from complex interference patterns that hinder accurate quantification. This
software identifies and corrects interference, as well increases the confidence
of peptide calls according to three qualities of an ideal chromatogram:
relative transition ratios similar to synthetic standards, co-elution of all
individual target transitions with respective stable isotope labeled (SIL)
transitions, and a distinguished peak relative to the background. Additionally,
these criteria should be satisfied for all replicate injections of a given
sample.

More information
================
Follow the github URL below for more information.
"""

setup(name='madic',
      packages=['madic'],
      version=__version__,
      description='MAtrix-Dependent Interference Correction of targeted mass '
      'spectrometry data',
      long_description=long_description,
      keywords='mass spectrometry SRM MRM Skyline LC-MS/MS',
      author='Derek Croote',
      author_email='dcroote@stanford.edu',
      url='https://github.com/dcroote/madic',
      install_requires=[
          'numpy',
          'pandas>=0.20',
          'scipy>=0.14',
          'matplotlib'
      ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      entry_points={
          'console_scripts': ['madic=madic.command_line:main']},
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: MacOS',
          'Operating System :: POSIX',
          'Operating System :: Unix']
      )
