e2o_getvar - retrieve variables
===============================

e2o_getval gets variables from the server for a specific region (specified in the ini file)
and optionally resamples those to a specified Digital Elevation Model (DEM).
e2o_getvar outputs daily values in the original units. The script does
physical downscaling only for Temperature and Precipitation. The other variables are interpolated only.
In order to downscale the precipitation data you must first download the data from:
https://github.com/earth2observe/downscaling-tools/releases/tag/2016-data-prec.
Next, check if the resolution you want is available in the data/Prec
directory ( we use the WorldClim data, see Hijmans et. al. 2005) and make files for the new resolution.
You can use the newres.bat script as an example. If you have other
high resolution climatology data you can use that instead of the worldclim data.

Use the e2o_calculateEvaporation script for more elaborate downscaling options for radiation.

Usage:
::

    e2o_getvar.py -I the_ini_file


ini file configuration
----------------------


The .ini file below shows the available options

.. literalinclude:: _static/e2o_getvar.ini


An example ini file be found :download:`here. <_static/e2o_getvar.ini>`




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
NearSurfaceWindSpeed                 Wind_daily_EI_025_     Wind_daily_E2OBS_    wind_speed
LapseRate                            lapseM_EI_025_                 -            air_temperature_lapse_rate
==================================== ====================== ==================== =========================================


Implementation
--------------



.. automodule:: e2o_dstools.e2o_getvar
    :members:


References
----------

Hijmans, R.J., S.E. Cameron, J.L. Parra, P.G. Jones and A. Jarvis, 2005.
Very high resolution interpolated climate surfaces for global land areas.
International Journal of Climatology 25: 1965-1978.