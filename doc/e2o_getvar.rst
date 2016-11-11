e2o_getvar - retrieve variables
===============================

e2o_getval gets variables from the server for a specific region (specified in the ini file)
and optionally resamples those to a specified Digital Elevation Model (DEM).
e2o_getvar outputs daily values in the original units. The script does not do
any physical downscaling other than simple interpolation. Use the e2o_calculateEvaporation script
for more elaborate downscaling options.

Usage:
::

    e2o_getvar.py -I the_ini_file


ini file configuration
----------------------

| Variable                          | Filenames WRR2       | Filenames WRR1      | standard_name                             |
|-----------------------------------|----------------------|---------------------|-------------------------------------------|
| Temperature                       | Tair_daily_EI_025_   | Tair_daily_E2OBS_   | air_temperature                           |
| DownwellingLongWaveRadiation      | LWdown_daily_EI_025_ | LWdown_daily_E2OBS_ | surface_downwelling_longwave_flux_in_air  |
| SurfaceAtmosphericPressure        | PSurf_daily_EI_025_  | PSurf_daily_E2OBS_  | surface_air_pressure                      |
| NearSurfaceSpecificHumidity       | Qair_daily_EI_025_   | Qair_daily_E2OBS_   | specific_humidity                         |
| Rainfall                          | Rainf_daily_EI_025_  | Rainf_daily_E2OBS_  | rainfal_flux                              |
| SurfaceIncidentShortwaveRadiation | SWdown_daily_EI_025_ | SWdown_daily_E2OBS_ | surface_downwelling_shortwave_flux_in_air |
| SnowfallRate                      | Snowf_daily_EI_025_  | Snowf_daily_E2OBS_  | snowfall_flux                             |
| NearSurfaceWindSpeed              | Wind_daily_EI_025_   | Wind_daily_E2OBS_   | wind_speed                                |
| LapseRate                         | lapseM_EI_025_       |                     | air_temperature_lapse_rate                |



The .ini file below shows the available options

.. literalinclude:: _download/example1.ini


An example ini file be found :download:`here. <_download/example1.ini>`

An example config that downloads 1 yr of Precipitation and Temperature data can be found in the
examples/getvar directory.

Implementation
--------------



.. automodule:: e2o_dstools.e2o_getvar
    :members:

