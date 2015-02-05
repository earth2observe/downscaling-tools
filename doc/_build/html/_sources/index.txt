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


The figure below shows the steps used to generated down-scale reference evaporation.

.. digraph:: steps

    "Clear sky radiation maps" [shape=box];    
    "Reference evaporation" [shape=box];    
    "e2o_radiation.py" -> "Clear sky radiation maps" [label =" Correct for aspect and slope with DEM"];
     "e2o_calculateEvaporation.py" -> "Reference evaporation" [label =" Downscale using DEM and clear-sky maps"]
    "Clear sky radiation maps" -> "e2o_calculateEvaporation.py"

    dpi=69;


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
