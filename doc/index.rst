===============================
eartH2Observe Downscaling tools
===============================

.. note::
      This documentation is for version |version| of e2o_dstools, release |release|
      This documentation was generated |today|

Introduction
============

e2o_downscaling-tools consists of a number of python programs and procedures that facilitate local
application of the earth2observe global water resources reanalysis. The tools
can connect directly to the project's data server and save (resampled) data to a local computer
for further analysis or direct application. The current first versions
of the tool focusses on downscaling the global forcing dataset used in
the project :cite:`weedonwfdei2014`.



Quick start using the examples
==============================

The steps below should allow you to run the examples on a windows computer.
Linux instructions are missing at the moment but in principle just running "python setup.py install" should get you started.

1. Download the release: https://github.com/earth2observe/downscaling-tools/releases/download/20015.1/e2o_dstools-64-bit-2015.1.zip
2. Download the examples: https://github.com/earth2observe/downscaling-tools/releases/download/20015.1/examples.zip
3. Unzip both files into an empty directory
4. Goto the examples\\getvar directory and run the examplerun1.bat file
5. the result are stored as PCRaster maps in the output directory. These can be opened in QGIS or aguila




The e2o_getvar script
=====================
.. toctree::
   :maxdepth: 2

   e2o_getvar

The radiation module
====================
.. toctree::
   :maxdepth: 2

   radiation

The evaporation module
======================
.. toctree::
   :maxdepth: 2

   evaporation




Examples and tests
==================
.. toctree::
   :maxdepth: 2
   

FAQ
===
.. toctree::
   :maxdepth: 2

   faq

   
Release notes
=============
.. toctree::
   :maxdepth: 2

   release-notes


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



TODO
====

.. todolist::


.. only:: html

   .. rubric:: References

.. bibliography:: e2o.bib
   :style: alpha
