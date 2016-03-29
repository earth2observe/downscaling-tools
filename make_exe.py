from distutils.core import setup
from bbfreeze import Freezer
import ctypes
import os
import glob, shutil
import netCDF4_utils

nrbits = str(ctypes.sizeof(ctypes.c_voidp) * 8)
data_files=[]

thename = "e2o_dstools-"+ str(nrbits)+"-bit"

f = Freezer("e2o_dstools-" + str(nrbits) +"-bit",\
            includes=('numpy','scipy','scipy.special._ufuncs_cxx','scipy.sparse.csgraph._validation'))
f.addScript("e2o_dstools/e2o_radiation.py")
f.addScript("e2o_dstools/e2o_getvar.py")
f.addScript("e2o_dstools/e2o_calculateEvaporation.py")
f()    # starts the freezing process



ddir = "c:/pcraster4-64/lib/"
data_files.append((".", glob.glob(ddir + "/*.dll")))

ddir = 'e2o_dstools/data/'
data_files.append(("./data", glob.glob(ddir + "/*.*")))


gdaldata = os.getenv("GDAL_DATA")
data_files.append(("./gdal-data", glob.glob(gdaldata + "/*.*")))

print data_files
print "Copying extra data files..."
for dirr in data_files:
    timake = os.path.join(thename ,dirr[0])
    print timake
    if not os.path.exists(timake):
        os.makedirs(timake)
    for tocp in dirr[1]:
        print tocp
        shutil.copy(tocp,timake)
