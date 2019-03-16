madic: MAtrix-Dependent Interference Correction
===============================================

madic is a Python package for detecting and correcting interference in targeted mass spectrometry data using Pandas data structures.

|travis-build-status| |appveyor-build-status| |coverage-status|

Purpose
-------

Multiple reaction monitoring (MRM, also known as selected reaction monitoring / SRM) data acquired during targeted mass spectrometry experiments can suffer from complex interference patterns that hinder accurate quantification. This software identifies and corrects interference, as well increases the confidence of peptide calls according to three qualities of an ideal chromatogram: relative transition ratios similar to synthetic standards, co-elution of all individual target transitions with respective stable isotope labeled (SIL) transitions, and a distinguished peak relative to the background. Additionally, these criteria should be satisfied for all replicate injections of a given sample.

See our `biorxiv preprint`_ for an illustration of the challenges of matrix interference when detecting allergenic protein in commercial food products.

See the `Jupyter Notebook tutorial`_ for a visual walk-through of interference correction in action.

.. image:: static/overview.png

Documentation and Tutorial
--------------------------

See `Read the Docs`_ for documentation and the `Jupyter Notebook tutorial`_ for a visual walk-through of interference correction.


Installation
------------

madic depends on a functioning Python 2.7 or 3.5+ installation (I recommend installing Anaconda_, which works on Windows, Linux, and macOS) and the following python packages: scipy_, numpy_, pandas_, and matplotlib_. For example:

.. code-block:: bash
    
    conda create -n madic_env python=3.5 numpy scipy pandas matplotlib

Installing the current release using the `pip` python package manager:

.. code-block:: bash

    pip install madic

Installing the development version (requires git):

.. code-block:: bash

    pip install git+https://github.com/dcroote/madic.git#egg=madic


Usage
-----

This tool is capable of processing SRM / MRM data from any instrument vendor so long as the raw data is capable of being loaded into Skyline_ version >= 4.1 (you can check the program version by going to Help --> About).

It is strongly preferred that injections contain stable isotope (also known as "heavy") labeled versions of all measured peptides. In the absence of these heavy peptides madic will fail to evaluate transition retention time variability and overall levels of confidence in the results may be somewhat compromised.

.. _data-prep-label:

Data Preparation
################

1. You will need to download the following Skyline report file to export data from Skyline with the correct formatting:
    + `Skyline transition report`_

2. Load the report file into Skyline using the following steps:
    + File --> Export --> Report --> Edit list --> Import

3. Two reports need to be exported:
    1. Sample data - this report should originate from a Skyline document that contains samples to be analyzed
    2. Reference data - this report should originate from a Skyline document that contains only standards or QC samples that provide reference transition ratios.

Program execution
#################

MADIC can be used as a standalone python executable or be imported into an interactive jupyter notebook for visualization (see *Tutorial* section below).

Command line usage:
::

    madic.py [options] csv ref_csv delimiter delimiter_pos

:code:`delimiter` and :code:`delimiter_pos` are necessary to identify samples from within replicate injection names. For example, :code:`_` (an underscore) for the :code:`delimiter` and :code:`1` for the :code:`delimiter_pos` would correctly identify the two samples (sampleA and sampleB) from injections with the following names:

    + Study1_sampleA_injection1
    + Study1_sampleA_injection2
    + Study1_sampleB_injection1
    + Study1_sampleB_injection2

(Remember Python uses zero-based indexing, so after splitting each replicate name into an array of three elements using the underscore, the second element is referred to as ``delimiter_pos`` one.)


Testing
-------

Uses :code:`pytest`. To test madic, run :code:`python setup.py test` in the source directory.

Development
-----------

Please submit issues to the `issue tracker`_.

Pull requests welcome!

License
-------

`BSD 3-Clause License`_


.. _Skyline transition report: https://raw.githubusercontent.com/dcroote/madic/master/static/madic_skyline_transition_report_and_chromatogram.skyr
.. _Skyline: https://skyline.ms/
.. _biorxiv preprint: https://www.biorxiv.org/content/early/2017/12/12/231266
.. _Anaconda: https://conda.io/docs/user-guide/install/index.html
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/
.. _pandas: https://pandas.pydata.org/
.. _matplotlib: https://matplotlib.org/
.. _Jupyter Notebook tutorial: https://github.com/dcroote/madic/blob/master/examples/tutorial.ipynb
.. _issue tracker: https://github.com/dcroote/madic/issues
.. _BSD 3-Clause License: https://github.com/dcroote/madic/blob/master/LICENSE
.. _Read the Docs: http://madic.rtfd.io/

.. |travis-build-status| image:: https://travis-ci.org/dcroote/madic.svg?branch=master
    :alt: travis build status
    :target: https://travis-ci.org/dcroote/madic

.. |appveyor-build-status| image:: https://ci.appveyor.com/api/projects/status/oiir9453oqluvmpm/branch/master?svg=true
    :alt: appveyor build stats
    :target: https://ci.appveyor.com/project/dcroote/madic/branch/master

.. |coverage-status| image:: https://coveralls.io/repos/github/dcroote/madic/badge.svg?branch=master
    :alt: coverage
    :target: https://coveralls.io/github/dcroote/madic?branch=master
