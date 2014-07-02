"""
Get a variable from the forcing data from the e2o server for a specific region and time range


usage:

    e2o_getvar.py -I inifile

    -I inifile - ini file with settings which data to get
"""

#TODO: Make initial version

import getopt, sys, os, netCDF4
import osgeo.gdal as gdal
from osgeo.gdalconst import *




serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
wrrsetroot = "ecmwf/met_forcing_v0/"
variable = "Tair_daily_E2OBS_"
startyear = 1979
endyear= 1980
startmonth = 1
endmonth = 12
latmin = 51
latmax = 52
lonmin = 3
lonmax = 4

ncurl = "http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/1980/Tair_daily_E2OBS_198001.nc"

def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2

    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data,0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'



class ncdatset():
    """
    Wrapper around the nc object to simplify things
    
    Opens the dataset and determines the number of dimensions
    3 = T, Lat Lon
    4 = T heigth, Lat Lon
    """
    def __init__(self,ncurl):

        self.nc = netCDF4.Dataset(ncurl)
        self.lat = self.getlat(self.nc)
        self.lon = self.getlon(self.nc)
        self.heigth = self.getheigth(self.nc)
        self.dimensions = 3 if self.heigth == None else 4
        self.time =   self.gettime(self.nc)
        self.timesteps = self.time.shape[0]


    def __del__(self):
        self.nc.close()
        
        
    def getlat(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'latitude':
                return ncdataset.variables[a]

        return None

    def getvarbyname(self,name):
        """
        """

        for a in self.nc.variables:
            print self.nc.variables[a].standard_name
            if  self.nc.variables[a].standard_name == name:
                return self.nc.variables[a]

        return None

    def gettime(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'time':
                return ncdataset.variables[a]

        return None
        

    def getlon(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'longitude':
                return ncdataset.variables[a]

        return None


    def getheigth(self,ncdataset):
        """
        """

        for a in ncdataset.variables:
            if  ncdataset.variables[a].standard_name == 'height':
                return ncdataset.variables[a]

        return None





a = "http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/1979/Tair_daily_E2OBS_197901.nc"
#def main(argv=None):
"""
Perform command line execution of the model.
"""

#if argv is None:
#    argv = sys.argv[1:]
#    if len(argv) == 0:
#        usage()
#        return
#
argv =["lll"]
try:
    opts, args = getopt.getopt(argv, 'I:')
except getopt.error, msg:
    usage(msg)

for o, a in opts:
    if o == '-I': inifile = a


# Generate the data range
years = arange(startyear,endyear+1,1)
months = arange(startmonth, endmonth + 1,1)

ncflist = []
for year in years:
    for month in months:
        ncflist.append(serverroot + wrrsetroot + "%d" % (year) + "/" + variable + "%d%02d.nc" % (year,month))

# first assumption: all files have the same geo dimentsion
BB = dict(
       lon=[ lonmin, lonmax],
       lat= [ latmin, latmax]
       )
    
print ncflist[0]
print a
nc = ncdatset(ncflist[0])
lat = nc.lat[:]
lon = nc.lat[:]
(latidx,) = logical_and(lat >= BB['lat'][0], lat < BB['lat'][1]).nonzero()
(lonidx,) = logical_and(lon >= BB['lon'][0], lon < BB['lon'][1]).nonzero()
print lat
print lon
print latidx
print lonidx
 
myvar = nc.getvarbyname("air_temperature")   
    
if nc.dimensions == 3:
    window = myvar[:,latidx.min():latidx.max(),lonidx.min():lonidx.max()] 
else:
    window = myvar[:,0,latidx.min():latidx.max(),lonidx.min():lonidx.max()]     
  

ts_mean = window.mean(axis=2).mean(axis=1)   
    
#main(argv="ggg")
