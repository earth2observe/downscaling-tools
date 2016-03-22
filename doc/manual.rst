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

3.2.2	Evaporation options
To turn the generation of evaporation grids on or off the calculateEvap option is given.

calculateEvap = True	 	 evaporation grids will be generated
calculateEvap = False	 no evaporation grids will be generated (reduces run times when no evaporation is needed).


Three methods to calculate evaporation have been implemented in the down-scaling tool:

+ Penman-Monteith : a physically based equation considering most relevant atmospheric processes

+ Priestley-Taylor : a substitute of the Penman-Monteith equation where the aerodynamic term has been replaced by an
empirical multiplier

+ Hargreaves : A simplified form of the Penman-Monteith equation using temperature and an annual radiation cycle as input

For more information on the implemented equations and a brief comparison we refer to Sperna Weiland et al. (2015).

The equation to be used can simply be selected by removing the hashdeck in front of the specific method making sure
the hashdeck is present before all other methods. In the example below Penman-Monteith will be used.

::

    # Choose one of the three methods below
    evapMethod = PenmanMonteith
    #evapMethod = Hargreaves
    #evapMethod = PriestleyTaylor

*Resampling and downscaling*

With the option resampling the user can select whether the data needs to be resampled to the resolution of the by the
 user provided DEM (see section 3.1).

If the option resampling is set to True a second optimization can be chosen with the downscaling option. If the
downscaling option is set to True temperature and air pressure will be corrected based on the difference in
altitude in the high-resolution user specified DEM and the low-resolution DEM that belongs to the WRR1 or WRR2
datasets.

These DEMs are located in the folder : ../e2o_downscaling/e2o_dstools/lowresdem and are called demWRR1.tif and
demWRR2.tif. The downscaling tool automatically selects the correct DEM based on the selected meteorological forcing
(met_forcing_v0 or met_forcing_v1) defined in the .ini file at wrrsetroot.   When both downscaling and resampling are
 set to false the maximum spatial extend required for the data to be read from the netCDFs file can be set by
 defining the corners of the area of interest: latmin, latmax, lonmin and lonmax.

If one is for example only interested in data for Australia the process can be accelerated by avoiding the reading of
 the full world maps from the netCDFs file by setting an extend slightly larger than the Australian continent.

::

    # Specify bounding box to download from server. Should be a bit bigger that the DEM
    latmin = -45
    latmax = -4
    lonmin = 110
    lonmax = 155

*Variable lapse rate*

For the downscaling of temperature, air pressure and radiation from the WRR1 dataset only a constant lapse rate of -0
.006 degrees/m can be used. The WRR2 datasets contains monthly fields of spatially and temporal varying lapse rates
 – derived from atmospheric conditions. To use these varying lapse rate fields for the downscaling the following
 option should be set to True:

::

    # useVarLapseRate = True -> use spatial and temporal varying lapse rate provided as part of the WRR2 forcing dataset
    # in stead of the default value of -0.006
    useVarLapseRate = True


*Time period*

The WRR1 and WRR2 datasets are available for the period 01-01-1979 to 31-12-2012. The user can specify the period of
interest, see the example below for 1979:

::

    # Start and end-year, month and day of the evaporation calculations
    startyear = 1979
    endyear= 1979
    startmonth = 1
    endmonth = 12
    startday = 1
    endday = 31


3.2.7	Radiation correction

The WRR1 and WRR2 provide potential solar radiation which is the radiation of an unobstructed or cloudless sky. The
magnitude of this potential solar radiation that reaches the earth surface depends on the position of the sun the
solar altitude or solar angleduring the day, the inclination of the solar rays with the earth’s surface, the amount
of radiation at the outer layer of the earth’s atmosphere, the transmissivity of the sky and the altitude of the
earth’s surface.

With the high resolution DEM the potential solar radiation can be corrected for aspect and shading. The correction
for cloudiness and other back scatter is derived from the transmissivity of the air and the path length radiation
needs to travel before reaching the earth’s surface.

The coefficient for radiation correction are calculated in the radiation sub-routine which will be described in
section 4. The directory where the correction files are located should be defined in the ini file:

::

    [downscaling]
    # Where to find the output of the e2o_radiation script
    radcordir=output_rad

Below you will find the filenames and there content:
COR00000.??? - Total clear sky radiation on DEM
SUN00000.??? - Nr of time intervals a pixel was in the sun
FLAT0000.??? - Total clear sky radiation on a flat surface
CORDIR00.??? - Direct clear sky radiation on DEM
FLATDIR0.??? - Direct clear sky radiation on a flat surface

For full details see: REF to documentation Jaap

*Output*

The user can specify the format of the output files – any of the gdal formats can be selected. These can, together
with their shortnames, be found at:
http://www.gdal.org/formats_list.html

::

    [output]
    # Gdal output format string
    # See: http://www.gdal.org/formats_list.html
    # examples: AAIGrid, PCRaster, GTiff etc
    format=	GTiff

The user can indicate the output location where all files should be stored:

::

    directory=output/


The first letters of the evaporation output files are set with the prefix:

::

    prefix=PET

If all other meteorological variables need to be saved the “saveall” option should be set to true.

::

    # If saveall is true all variables used are saved instead of only the PET
    saveall=1



Example e2o_calculateevaporation ini file:

::
 
    [url]
    # Server location and location of the WRR forcing
    serverroot = http://wci.earth2observe.eu/thredds/dodsC/
    wrrsetroot = ecmwf/met_forcing_v1/

    [selection]
    # What to do
    calculateEvap = True
    # Choose one of the three methods below
    evapMethod = PenmanMonteith
    #evapMethod = Hargreaves
    #evapMethod = PriestleyTaylor

    # Specifye box to download from server. Should be a bit bigger that the DEM
    latmin = -90
    latmax = +90
    lonmin = -180
    lonmax = 180

    # Start and end-year, month and day of the evaporation calculations
    startyear = 1979
    endyear= 1979
    startmonth = 1
    endmonth = 12
    startday = 1
    endday = 31

    [downscaling]
    # location of original DEM (WFDEI) and the local high resolution DEM
    highResDEM=highresdem\DEM.tif
    # Resampling = True -> resample to resolution of dEM specified in downscaling section
    # Downscaling = True -> also apply DEM based correctiosn of T, Radiation, Pressure
    resampling  = True
    downscaling = True
    # useVarLapseRate = True -> use spatial and temporal varying lapse rate provided as part of the WRR2 forcing dataset iso the default value of -0.006
    useVarLapseRate = True
    # Wher to fine the output of the e2o_radiation script
    radcordir=output_rad

    [output]
    # Gdal output format string
    # See: http://www.gdal.org/formats_list.html
    # examples: AAIGrid, PCRaster, GTiff etc
    format=	GTiff
    directory=output/
    prefix=PET
    # Is saveall is true all variables used are saved instead of only the PET
    saveall=1