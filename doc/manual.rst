e2o-downscaling tools - manual
==============================


Introduction
------------

Within this manual we introduce a Python based meteorological downscaling tool that allows one to:

+ retrieve all meteorological variables that are part of the eartH2Observe WRR1 and WRR2 datasets
 for a region of interest on a user defined grid extent and resolution;
+ downscale the meteorological variables temperature and air pressure using a DEM based elevation correction;
+ calculate potential evaporation from the WRR1 and WRR2 datasets using the Penman-Moneith, Priestley-Taylor
    or Hargreaves equation, optionally considering elevation corrections for temperature, air pressure and
    radiation and shading corrections for radiation.

This document is merely a technical manual. For background reading on the scientific concepts used we
will in the text refer to the online documentation.


Region specific user-defined settings
-------------------------------------

In the previous section we have introduced the default way of running the eartH2Observe downscaling tool. In this
section we discuss all options that are user adjustable and can be specified for the region, resolution and period of
 interest. Options will be discussed following the order in the .ini file.

User specified DEM of area of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The user can specify his region of interest, and the desired regular grid resolution by providing a high resolution
DEM covering the region of interest. The file should be provided in GEOTiff format with the name DEM.tif in the
following directory:

::

    ../e2o_downscaling/e2o_dstools/highresdem

The DEM will be used for the definition of the extend and resolution of the generated meteorological output files. If
 the option ‘downscaling’ is turned on (see section 3.2.4) the altitudes in the DEM will be used to spatially
 downscale temperature, air pressure and radiation.


User settings in the e2o_calculateEvaporation.ini file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below a description of  the content of the .ini file is given. Options can be turned on or off by:

+ for multiple options: 	adding / deleting a hashtag in front fo the option
+ for simple on or off: 	indication True (on) or False (off)


*Server location*

The server location refers to the online location where the eartH2Observe WRR datasets are hosted. At this stage
there are two datasets available:

1.	the WRR1 dataset which is directly derived from the global WATCH-Forcing-Data-ERA-Interim dataset (Weedon et al.,
2014) that was developed as processor of the WATCH forcing dataset developed in the EU FP6 project EU-WATCH (Weedon
et al., 2012). This dataset has a spatial resolution of 0.5 degrees and a 3-hourly temporal resolution.

2.	the WRR2 dataset in which several improvements have been made compared to the WRR1 dataset (see REF Dutra et al ).
 This dataset has a spatial resolution of 0.25 degrees and a 3-hourly temporal resolution.

The serverroot for both datasets is:

serverroot = http://wci.earth2observe.eu/thredds/dodsC/

The wrrsetroot depends on the dataset a user prefers to use:

+ For the WRR1 : 	wrrsetroot = ecmwf/met_forcing_v0/

+ For the WRR2 : 	wrrsetroot = ecmwf/met_forcing_v1/

*Evaporation options*

To turn the generation of evaporation grids on or off the calculateEvap option is given.

calculateEvap = True	-- evaporation grids will be generated
calculateEvap = False	-- no evaporation grids will be generated (reduces run times when no evaporation is needed).

