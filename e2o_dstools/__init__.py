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
import sys
import os


if hasattr(sys, "frozen"):
    _ROOT = os.path.abspath(os.path.dirname(__file__)).split("library.zip")[0]
    #os.environ['GDAL_DATA'] = os.path.join(_ROOT,'gdal-data')
else:
    _ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)