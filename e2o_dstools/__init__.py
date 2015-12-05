__all__ = ["e2o_utils","e2o_getvar","e2o_calculateEvaporation"]

__author__ = 'schelle'
__version__ = '0.1'
__release__ = "2015"


import osgeo.gdal as gdal
import netCDF4
import netcdftime
import scipy
import scipy.interpolate
import scipy.special