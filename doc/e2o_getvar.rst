e2o_getvar - retrieve variables
===============================

e2o_getval gets variables from the server for a specific region (specified in the ini file)
and optionally resamples those to a specified Digital Elevation Model (DEM).
e2o_getvar outputs daily values in the original units.

Usage:
::

    e2o_getvar.py -I the_ini_file


ini file configuration
----------------------
The .ini file below shows the available options

.. literalinclude:: _download/e2o_getvar.ini


An example ini file be found :download:`here. <_download/e2o_getvar.ini>`

An example config that downloads 1 yr of Precipitation and Temperature data can be found in the
examples/getvar directory.

Implementation
--------------



.. automodule:: e2o_dstools.e2o_getvar
    :members:

