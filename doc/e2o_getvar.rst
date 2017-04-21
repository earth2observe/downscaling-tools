e2o_getvar - retrieve variables
===============================

e2o_getval gets variables from the server for a specific region (specified in the ini file)
and optionally resamples those to a specified Digital Elevation Model (DEM).
e2o_getvar outputs daily values in the original units. The script does
physical downscaling only for Temperature and Precipitation. The other variables are interpolated only.
In order to downscale the precipitation data you must first download the data from:
https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec.
We use the WorldClim data, see Hijmans et. al. 2005. If you have other
high resolution climatology data you can use that instead of the worldclim data.

.. note::

    Use the e2o_calculateEvaporation script for more elaborate downscaling options for evaporation
    and the variables needed to calculate it.



Usage:
::

    e2o_getvar.py -I the_ini_file


ini file configuration
----------------------


The .ini file below shows the available options

.. literalinclude:: _static/examplerun1.ini



Table: Variables and names to be used in the ini file

==================================== ====================== ==================== =========================================
Variable                             Filenames WRR2         Filenames WRR1       standard_name
==================================== ====================== ==================== =========================================
Temperature                          Tair_daily_EI_025_     Tair_daily_E2OBS_    air_temperature
DownwellingLongWaveRadiation         LWdown_daily_EI_025_   LWdown_daily_E2OBS_  surface_downwelling_longwave_flux_in_air
SurfaceAtmosphericPressure           PSurf_daily_EI_025_    PSurf_daily_E2OBS_   surface_air_pressure
NearSurfaceSpecificHumidity          Qair_daily_EI_025_     Qair_daily_E2OBS_    specific_humidity
Rainfall                             Rainf_daily_MSWEP_025_ Rainf_daily_E2OBS_   rainfall_flux
SurfaceIncidentShortwaveRadiation    SWdown_daily_EI_025_   SWdown_daily_E2OBS_  surface_downwelling_shortwave_flux_in_air
SnowfallRate                         Snowf_daily_MSWEP_025_ Snowf_daily_E2OBS_   snowfall_flux
TotalPrecipitation                   -                      -                    precipitation_flux
NearSurfaceWindSpeed                 Wind_daily_EI_025_     Wind_daily_E2OBS_    wind_speed
LapseRate                            lapseM_EI_025_                 -            air_temperature_lapse_rate
==================================== ====================== ==================== =========================================

TotalPrecipitation is a virtual name and will instruct the script to sum Rainfall and SnowfallRate.


Using the examples
------------------

examplerun1.ini
~~~~~~~~~~~~~~~
This example downloads and resamples MSWEP rainfall data (WRR2). No downscaling is applied. All output
is stored in sequential geotiff files in output/Rainfall/P0000000.???, one for each timestep.

examplerun2.ini
~~~~~~~~~~~~~~~
This example downloads and resamples MSWEP rainfall data (WRR2). Downscaling is applied using nearest
interpolation. Therefore you need to have downloaded and installed the Prec.zip file
at: https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec
All output is stored in sequential geotiff files in output/Rainfall/PDS00000.???, one for each timestep.

examplerun3.ini
~~~~~~~~~~~~~~~
This example downloads and resamples MSWEP rainfall data (WRR2). Downscaling is applied using linear
interpolation. Therefore you need to have downloaded and installed the Prec.zip file
at: https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec
All output is stored in sequential geotiff files in output/Rainfall/PDSL0000.???, one for each timestep.


examplerun4.ini
~~~~~~~~~~~~~~~
As example 3 but the output is saved in a netcdf file (output/Rainfall.nc)

examplerun5.ini
~~~~~~~~~~~~~~~
This example downloads and resamples MSWEP snowfall data (WRR2). Downscaling is applied using linear
interpolation. Therefore you need to have downloaded and installed the Prec.zip file
at: https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec
All output is stored in sequential geotiff files in output/SnowfallRate/SNO00000.???, one for each timestep.


examplerun6.ini
~~~~~~~~~~~~~~~
This example downloads and resamples WFDEI TotalPrecipitation (Rain and Snow) data (WRR1). No downscaling is applied but
only resampling  using linear interpolation.
All output is stored in sequential geotiff files in output/TotalPrecipitation/PWWW10000.???, one for each timestep.



examplerun7.ini
~~~~~~~~~~~~~~~
This example downloads and resamples WFDEI Air Temperature  data (WRR1). Downscaling is applied using a fixed laps rate and
 linear interpolation.
All output is stored in sequential geotiff files in output/Temperature/TEMP0000.???, one for each timestep.






Implementation
--------------



.. automodule:: e2o_dstools.e2o_getvar
    :members:


References
----------

Hijmans, R.J., S.E. Cameron, J.L. Parra, P.G. Jones and A. Jarvis, 2005.
Very high resolution interpolated climate surfaces for global land areas.
International Journal of Climatology 25: 1965-1978.