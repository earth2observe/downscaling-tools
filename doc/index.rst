=============================
Introduction and installation
=============================

Introduction
============

.. note::
      This documentation is for version |version| of e2o_dstools, release |release|
      This documentation was generated |today|


e2o_downscaling-tools consists of a number of python programs and procedures that facilitate local
application of the earth2observe global water resources reanalysis (http://www.earth2observe.eu). The tools
can connect directly to the project's data server (http://wci.earth2observe.eu) and save (resampled)
data to a local computer for further analysis or direct application. The current first versions
of the tool focusses on downscaling the global forcing dataset used in
the project.


Installation
============

The latest release can be downloaded from github: https://github.com/earth2observe/downscaling-tools/releases

Installing the window exe distribution
--------------------------------------

Download the file e2o_downscale-RELEASE-normal-win32-64.zip (or similar) from the releases tab to you computer. Please note that
this file will only work on windows computers with a 64 bit operating system.

Make a directory in which to store the program (e.g c:\\dstools) and unzip the contents of the zip file
into this directory, you should now have something like this:

::

    c:\dstools\e2o_downscale-2016-3-normal-win32-64\
    c:\dstools\examples\



Now we can check if the program is installed successfully by runnin gone of the programs from the command prompt. First
open a windows command prompt (press the windows key and type cmd). Next type in the name (full path) of the program to
run (in this case e2o_getvar.exe) as shown below:

::

    c:\>c:\dstools\e2o_downscale-2016-3-normal-win32-64\e2o_getvar.exe

    Get a variable from the forcing data from the e2o server for a specific region and time range

    usage:

        e2o_getvar.py -I inifile [-l loglevel][-h]

        -I inifile - ini file with settings which data to get
        -l loglevel (must be one of DEBUG, WARNING, ERROR)


if you see the above you have installed the program succesfully!

To save you from having to type the full path to the program every time you want to run it you can add the
directory c:\\dstools\\e2o_downscale-2016-3-normal-win32-64\\ to you computer's search path (see e.g.
http://www.howtogeek.com/118594/how-to-edit-your-system-path-for-easy-command-line-access/ on how to do this)

If you want to downscale precipitation data you must first download the required data
(https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec)
and extract the data in the e2o_downscale-2016-3-normal-win32-64\\data directory.


Installing the python distribution
----------------------------------

As a first step make sure all the required packages are installed:

- pcraster python extensions - www.pcraster.eu
- netCDF4
- numpy
- gdal
- matplotlib
- scipy

The easiest option we have found is to use the Anaconda python distribution and use the conda program
to install the requirement packages. pcraster needs to be installed manually.

Download the zip file with the source code from github (see above) and unzip the file into an empty directory and unzip
it's contents. If you want to downscale precipitation data you must first download the required data
(https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec) and extract the data in the
e2o_downscale-2016-3-normal-win32-64/data directory.

In the example below we have unzipped the file
onto the root of the C-drive. We can see the following files:

::

    c:\src\downscaling-tools-2016.3\data\
    c:\src\downscaling-tools-2016.3\doc\
    c:\src\downscaling-tools-2016.3\e2o_dstools\
    c:\src\downscaling-tools-2016.3\examples\
    c:\src\downscaling-tools-2016.3\setup.py
    c:\src\downscaling-tools-2016.3\make_exe.py
    c:\src\downscaling-tools-2016.3\README.rst

Open a windows command prompt (press the windows key and type cmd) and navigate to the directory in which the file
setup.py is located and run the setup script:

::

    c:>cd  src\downscaling-tools-2016.3
    c:>python setup.py install

This well give you a lot of output, ending with something similar to what is shown below:

::

    removing 'build\bdist.win-amd64\egg' (and everything under it)
    Processing e2o_dstoools-0.1-py2.7.egg
    Removing c:\anaconda\lib\site-packages\e2o_dstoools-0.1-py2.7.egg
    Copying e2o_dstoools-0.1-py2.7.egg to c:\anaconda\lib\site-packages
    e2o-dstoools 0.1 is already the active version in easy-install.pth
    Installing e2o_getvar.py script to C:\Anaconda\Scripts
    Installing e2o_calculateEvaporation.py script to C:\Anaconda\Scripts

    Installed c:\anaconda\lib\site-packages\e2o_dstoools-0.1-py2.7.egg
    Processing dependencies for e2o-dstoools==0.1
    Finished processing dependencies for e2o-dstoools==0.1





===========
User manual
===========

.. toctree::
   :maxdepth: 3

   manual.rst


================
Reference manual
================

The e2o_getvar script
=====================
.. toctree::
   :maxdepth: 3

   e2o_getvar

The e2o_radiation script
========================
.. toctree::
   :maxdepth: 2

   radiation

The e2o_calculateEvaporation script
===================================
.. toctree::
   :maxdepth: 2

   evaporation

===
FAQ
===

.. toctree::
   :maxdepth: 2

   faq

=============
Release notes
=============
.. toctree::
   :maxdepth: 2

   release-notes

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


====
TODO
====

.. todolist::


.. only:: html

   .. rubric:: References

.. bibliography:: e2o.bib
   :style: alpha
