"""
Get a variable from the forcing data from the e2o server for a specific region and time range


usage:

    e2o_getvar.py -I inifile

    -I inifile - ini file with settings which data to get
"""

#TODO: Add local cache
#TODO: Fix Hargreaves problem

import getopt, sys, os, netCDF4
import osgeo.gdal as gdal
from osgeo.gdalconst import *
import datetime
from numpy import *
import numpy as np
from e2o_utils import *
import pdb
#import pcraster as pcr


#pcr.setglobaloption("radians")



#ncurl = "http://wci.earth2observe.eu/thredds/dodsC/ecmwf/met_forcing_v0/1980/Tair_daily_E2OBS_198001.nc"

def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)





class ncdatset():
    """
    Wrapper around the nc object to simplify things
    
    Opens the dataset and determines the number of dimensions
    3 = T, Lat Lon
    4 = T heigth, Lat Lon
    """
    def __init__(self,ncurl,logger):

        self.nc = netCDF4.Dataset(ncurl)
        self.logger = logger
        self.lat = self.getlat(self.nc)
        if self.lat == None:
            self.logger.error("No lat information found!")
        self.lon = self.getlon(self.nc)
        if self.lon == None:
            self.logger.error("No lon information found!")

        self.heigth = self.getheigth(self.nc)
        self.dimensions = 3 if self.heigth == None else 4
        self.time =   self.gettime(self.nc)
        self.timesteps = self.time.shape[0]
        self.logger.debug(self.nc)


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





def get_times_daily(startdate,enddate, serverroot, wrrsetroot, filename,logger):
    """
    generate a dictionary with date/times and the NC files in which the data resides
    """

    numdays = enddate - startdate
    dateList = []
    filelist = {}
    for x in range (0, numdays.days + 1):
        dateList.append(startdate + datetime.timedelta(days = x))
    
    for thedate in dateList:
        ncfile = serverroot + wrrsetroot + "%d" % (thedate.year) + "/" + filename + "%d%02d.nc" % (thedate.year,thedate.month)
        filelist[str(thedate)] = ncfile
    
    return filelist, dateList


def get_times(startdate,enddate, serverroot, wrrsetroot, filename, timestepSeconds, logger):
    """
    generate a dictionary with date/times and the NC files in which the data resides for flexible timestep
    """

    numdays = enddate - startdate
    dateList = []
    filelist = {}
    #days
    for x in range (0, numdays.days + 1):
        #selected time-step in seconds
        delta = 0
        while delta < 86400:
            dateList.append(startdate + datetime.timedelta(seconds = delta))
            delta += timestepSeconds
    
    for thedate in dateList:
        ncfile = serverroot + wrrsetroot + "%d" % (thedate.year) + "/" + filename + "%d%02d.nc" % (thedate.year,thedate.month)
        filelist[str(thedate)] = ncfile
    
    return filelist, dateList


class getstepdaily():
    """
    class to get data from a set of NC files
    Initialise with a list of netcdf files and a variable name (standard_name)

    """
    
    def __init__(self,nclist,BB,varname,logger):
        """
        """
        self.o_nc_files = []
        self.list = nclist
        self.varname = varname
        self.BB = BB
        self.latidx = []
        self.lonidx =[]
        self.lat=[]
        self.lon=[]
        self.data = []
        self.logger = logger

    def getdate(self,thedate):

        datestr = str(thedate)
        lat = None
        lon = None
        window = None
        
        if datestr in self.list.keys():
            self.dset = ncdatset(self.list[datestr],self.logger)
            data = self.dset.getvarbyname(self.varname)
            lat = flipud(self.dset.lat[:])
            
            lon = self.dset.lon[:]
            (latidx,) = logical_and(lat >= self.BB['lat'][0], lat < self.BB['lat'][1]).nonzero()
            (lonidx,) = logical_and(lon >= self.BB['lon'][0], lon < self.BB['lon'][1]).nonzero()

            time = self.dset.time
            timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)
            
            dpos = thedate.day -1
            
            if self.dset.dimensions ==3:
                window = data[dpos,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[dpos,0,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            
            self.lat = lat[latidx]
            self.lon = lon[lonidx]
    
        else:
            self.logger.error( "cannot find: " + datestr)
            
        return self.lat, self.lon, window

    def getdates(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files for speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.debug("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret
        
    def getdates_seconds(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files fro speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.debug("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret
        
class getstep():
    """
    class to get data from a set of NC files for user defined timestep in seconds
    Initialise with a list of netcdf files and a variable name (standard_name)

    """
    
    def __init__(self,nclist,BB,varname,timestepSeconds,logger):
        """
        """
        self.o_nc_files = []
        self.list = nclist
        self.varname = varname
        self.BB = BB
        self.latidx = []
        self.lonidx =[]
        self.lat=[]
        self.lon=[]
        self.data = []
        self.logger = logger

    def getdate(self,thedate):

        datestr = str(thedate)
        lat = None
        lon = None
        window = None
        
        if datestr in self.list.keys():
            self.dset = ncdatset(self.list[datestr],self.logger)
            data = self.dset.getvarbyname(self.varname)
            lat = flipud(self.dset.lat[:])
            
            lon = self.dset.lon[:]
            (latidx,) = logical_and(lat >= self.BB['lat'][0], lat < self.BB['lat'][1]).nonzero()
            (lonidx,) = logical_and(lon >= self.BB['lon'][0], lon < self.BB['lon'][1]).nonzero()

            time = self.dset.time
            timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)
            
            if timestepSecond < 3600:
                dpos = thedate.second -1
            elif timestepSecond < 86400:
                dpos = thedate.hour -1
            else:
                dpos = thedate.day -1
            
            if self.dset.dimensions ==3:
                window = data[dpos,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[dpos,0,latidx.min():latidx.max()+1,lonidx.min():lonidx.max()+1]
            
            self.lat = lat[latidx]
            self.lon = lon[lonidx]
    
        else:
            self.logger.error( "cannot find: " + datestr)
            
        return self.lat, self.lon, window

    def getdates(self,alldates):
        """
        Does not work yet
        """
        lat = None
        lon = None
        ret = []
        
        lastnc = None
        
        # here we loop over nc files fro speed reasons
        
        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)
            
            time = self.dset.time
            tar = time[:]
            timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            #print timeObj
            spos = nonzero(timeObj == alldates[0])[0]
            if len(spos) != 1:
                spos = 0
            else:
                spos = int(spos)
 
            epos = nonzero(timeObj == alldates[-1])[0]
            if len(epos) != 1:
                epos = len(tar)
            else:
                epos = int(epos + 1)
            
            
            self.logger.debug("Processing url: " + theone )
            
            data = self.dset.getvarbyname(self.varname)
            
            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )
                
            lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()


            if self.dset.dimensions ==3:
                window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]
            if self.dset.dimensions ==4:
                window = data[spos:epos,0,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

            
            self.lat = lat[self.latidx]
            self.lon = lon[self.lonidx]
            
            if len(ret) == 0:
                ret = window.copy()
            else:
                ret = vstack((ret,window))
        
        return ret
                         
def save_as_mapsstack(lat,lon,data,times,directory,prefix="E2O",oformat="PCRaster"):        
    
    cnt = 0
    if not os.path.exists(directory):
        os.mkdir(directory)
    for a in times:
            below_thousand = cnt % 1000
            above_thousand = cnt / 1000
            mapname  = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)
            #print "saving map: " + os.path.join(directory,mapname)
            writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[cnt,:,:]),-999.0)
            cnt = cnt + 1    

def save_as_mapsstack_per_day(lat,lon,data,ncnt,directory,prefix="E2O",oformat="PCRaster"):        
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    below_thousand = ncnt % 1000
    above_thousand = ncnt / 1000
    mapname  = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)
    #print "saving map: " + os.path.join(directory,mapname)
    writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[:,:]),-999.0)
    
def PenmanMonteith(relevantDataFields, Tmax, Tmin):
    
    """
    relevantDataFields : ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']
    """
    Tmean   =  relevantDataFields[0]
    Rlin    =  relevantDataFields[1]
    Pres    =  relevantDataFields[2]
    Q       =  relevantDataFields[3]
    Rsin    =  relevantDataFields[4]
    Wsp     =  relevantDataFields[5]
              
    """
    Computes Penman-Monteith reference evaporation
    Inputs:
        Rsin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) incoming shortwave radiation [W m-2]
        Rlin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) incoming longwave radiation [W m-2]
        Tmean:       netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily mean temp [K]
        Tmax:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily max. temp [K]
        Tmin:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) daily min. temp [K]
        Wsp:         netCDF obj or NumPy array   -- 3D array (time, lat, lon) wind speed [m s-2]
        Q:           netCDF obj or NumPy array   -- 3D array (time, lat, lon) specific humididy [kg kg-1]
        Pres:        netCDF obj or NumPy array   -- 3D array (time, lat, lon) Surface air pressure [Pa]
        Pet:         netCDF obj or NumPy array   -- 3D array (time, lat, lon) for target data
    Outputs:
        trg_var:    netCDF obj or NumPy array   -- 3D array (time, lat, lon) for target data, updated with computed values
    """
    cp           = 1013         # specific heat of air 1013 [J kg-1 K-1]
    TimeStepSecs = 86400        # timestep in seconds
    karman       = 0.41         # von Karman constant [-]
    vegh         = 0.12         # vegetation height [m] 
    refh         = 2            # reference height where wind speed is measured
    alpha        = 0.23         # albedo, 0.23 [-]
    rs           = 70           # surface resistance, 70 [s m-1]
    R            = 287.058      # Universal gas constant [J kg-1 K-1]
    convmm       = 1000*TimeStepSecs # conversion from meters to millimeters
    sigma        = 5.6704e-8    # stephan boltzmann [W m-2 K-4]
    eps          = 0.622        # ratio of water vapour/dry air molecular weights [-]
    # outgoing long-wave [W m-2]: assume that surface temperature is equal to temperature close to the ground. Assume emissivity = 1.
    # net radiation [W m-2]
    Rlout = sigma*Tmean**4. # emissivity assumed 1
    Rlnet = Rlin - Rlout
    Rnet  = (1-alpha)*Rsin + Rlnet
    
    # saturation vapour pressure [Pa]
    es = lambda T:610.8*np.exp((17.27*(Tmean-273.15))/((Tmean-273.15)+237.3))
    es_min  = es(Tmin)
    es_max  = es(Tmax)
    es_mean = (es_min+es_max)/2.

    # actual vapour pressure
    
    """
        eps*ea
    q = ------
        P-((1-eps)*ea)
        
            q*P    
    ea =- ---------------
           (eps-1)*q-eps
    """
    ea = lambda Pres, Q, eps: -(Q*Pres)/((eps-1)*Q-eps)
    ea_mean = ea(Pres, Q, eps)
    
#    ea = lambda rel_hum, es: rel_hum/100.*es
#    ea_min  = ea(Hr, es_min)
#    ea_max  = ea(Hr, es_max)
#    ea_mean = (ea_min + ea_max)/2
    
    # vapour pressure deficit
    vpd = np.maximum(es_mean - ea_mean, 0.)
    
    # density of air [kg m-3]
    rho = Pres/(Tmean*R)
    
    # Latent heat [J kg-1]
    Lheat = (2.501-(0.002361*(Tmean-273.15)))*1e6

    # slope of vapour pressure [Pa K-1]
    deltop  = 4098. *(610.8*np.exp((17.27*(Tmean-273.15))/((Tmean-273.15)+237.3)))
    delbase = ((Tmean-273.15)+237.3)**2
    delta   = deltop/delbase

    # psychrometric constant
    gamma   = cp*Pres/(eps*Lheat)
    
    # aerodynamic resistance
    z = 10 # height of wind speed variable (10 meters above surface)
    Wsp_2 = Wsp*4.87/(np.log(67.8*z-5.42))
    ra = 208./Wsp_2
    
    PETtop  = delta*Rnet + rho*cp*vpd/ra
    PETbase = delta + gamma*(1+rs/ra)
    PET     = np.maximum(PETtop/PETbase, 0)
    PETmm   = PET/Lheat*TimeStepSecs
    return PETmm #, Lheat, Rlnet, ea_mean, es_mean, Wsp_2



def hargreaves(lat, currentdate, relevantDataFields, Tmax, Tmin):
    
    
    #get day of year
    tt  = currentdate.timetuple()
    JULDAY = tt.tm_yday
#    #Latitude radians
    LatRad= lat*np.pi/180.0
    test = np.tan(LatRad)
#    ### water euivalent extraterrestial radiation ###    
#    # declination (rad)
    declin = 0.4093*(np.sin(((2.0*pi*JULDAY)/365.0)-1.405))
#    # sunset hour angle
    arccosInput = (-(np.tan(LatRad))*(np.tan(declin)))
#    
    arccosInput = np.minimum(1,arccosInput)
    arccosInput = np.maximum(-1,arccosInput)
    sunangle = np.arccos(arccosInput)
#    # distance of earth to sun
    distsun = 1+0.033*(np.cos((2*pi*JULDAY)/365.0))
#    # SO = water equivalent extra terrestiral radiation in mm/day
    Ra = 15.392*distsun*(sunangle*(np.sin(LatRad))*(np.sin(declin))+(np.cos(LatRad))*(np.cos(declin))*(np.sin(sunangle)))
    strDay       = str(JULDAY)
#    while len(strDay) < 3:
#        strDay  = str(0)+ strDay
#    fileName    = 'so000000.' + strDay
#    RaMap       = pcr.readmap(os.path.join('maps',fileName))
#    Ra          = np.flipud(pcr.pcr2numpy(RaMap,0))
    
    #print 'test2'
        
    airT = relevantDataFields[0]
    PETmm = 0.0023*Ra*((np.maximum(0,(airT-273.0))) + 17.8)*sqrt(np.maximum(0,(Tmax-Tmin)))
 
    return PETmm, Ra, LatRad, sunangle, declin

#### MAIN ####

def main(argv=None):

    serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
    wrrsetroot = "ecmwf/met_forcing_v0/"
    
    #available variables with corresponding file names and standard_names as in NC files
    variables = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','Rainfall','SurfaceIncidentShortwaveRadiation','SnowfallRate','NearSurfaceWindSpeed']
    filenames = ["Tair_daily_E2OBS_","LWdown_daily_E2OBS_","PSurf_daily_E2OBS_","Qair_daily_E2OBS_",\
                    "Rainf_daily_E2OBS_","SWdown_daily_E2OBS_","Snowf_daily_E2OBS_","Wind_daily_E2OBS_"]
    standard_names = ['air_temperature','surface_downwelling_longwave_flux_in_air','surface_air_pressure','specific_humidity',\
                        'rainfal_flux','surface_downwelling_shortwave_flux_in_air','snowfall_flux','wind_speed']
    
    #defaults in absence of ini file
    filename = "Tair_daily_E2OBS_"
    standard_name ='air_temperature'
    startyear = 1979
    endyear= 1980
    startmonth = 1
    endmonth = 12
    latmin = 51.25
    latmax = 51.75
    lonmin = 5.25
    lonmax = 5.75
    startday = 1
    endday = 1
    getDataForVar = True
    calculateEvap = False
    evapMethod = None
    
    
    #argv = ["-I","e2o_getvar.ini"]
    
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            exit()   
    try:
        opts, args = getopt.getopt(argv, 'I:')
    except getopt.error, msg:
        usage(msg)
    
    for o, a in opts:
        if o == '-I': inifile = a
            
    logger, ch = setlogger("e2o_getvar.log","e2o_getvar",level=logging.INFO)
    logger.debug("Reading settings from ini: " + inifile)
    theconf = iniFileSetUp(a)
    
    # Read period from file
    startyear = int(configget(logger,theconf,"selection","startyear",str(startyear)))
    endyear = int(configget(logger,theconf,"selection","endyear",str(endyear)))
    endmonth = int(configget(logger,theconf,"selection","endmonth",str(endmonth)))
    startmonth = int(configget(logger,theconf,"selection","startmonth",str(startmonth)))
    endday = int(configget(logger,theconf,"selection","endday",str(endday)))
    startday = int(configget(logger,theconf,"selection","startday",str(startday)))
    start = datetime.datetime(startyear,startmonth,startday)
    end = datetime.datetime(endyear,endmonth,endday)
    
    #read remaining settings from in file        
    lonmax = float(configget(logger,theconf,"selection","lonmax",str(lonmax)))
    lonmin = float(configget(logger,theconf,"selection","lonmin",str(lonmin)))
    latmax = float(configget(logger,theconf,"selection","latmax",str(latmax)))
    latmin = float(configget(logger,theconf,"selection","latmin",str(latmin)))
    BB = dict(
           lon=[ lonmin, lonmax],
           lat= [ latmin, latmax]
           )
    #serverroot = configget(logger,theconf,"url","serverroot",serverroot)
    #wrrsetroot = configget(logger,theconf,"url","wrrsetroot",wrrsetroot)
    oformat = configget(logger,theconf,"output","format","PCRaster")
    #odir = configget(logger,theconf,"output","directory","output/")
    oprefix = configget(logger,theconf,"output","prefix","E2O")
       
    # Check whether evaporation should be calculated
    calculateEvap   = configget(logger,theconf,"selection","calculateEvap",calculateEvap)
    
    if calculateEvap == 'True':
        evapMethod      = configget(logger,theconf,"selection","evapMethod",evapMethod)
        
    if evapMethod == 'PenmanMonteith':
        relevantVars = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','SurfaceIncidentShortwaveRadiation','NearSurfaceWindSpeed']        
    elif evapMethod == 'Hargreaves':
        relevantVars = ['Temperature']        
    
    currentdate = start
    ncnt = 0
    while currentdate <= end:
        # Get all daily datafields needed and aad to list
        relevantDataFields = []
        for i in range (0,len(variables)):
            odir = configget(logger,theconf,"output","directory","output/")
            if variables[i] in relevantVars:
                logger.info("Getting data field: " + filename)
                filename = filenames[i]
                standard_name = standard_names[i]
    
                tlist, timelist = get_times_daily(currentdate,currentdate,serverroot, wrrsetroot, filename,logger )             
                ncstepobj = getstepdaily(tlist,BB,standard_name,logger)
                mstack = ncstepobj.getdates(timelist)
    
                relevantDataFields.append(mstack)
                
        if evapMethod == 'PenmanMonteith':
            # retrieve 3 hourly Temperature and calculate max and min Temperature            
            filename = 'Tair_E2OBS_'
            standard_name = 'air_temperature'
            timestepSeconds = 10800
            
            tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )       
            ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)       
            mstack = ncstepobj.getdates(timelist)                       
            tmin = mstack.min(axis=0)            
            tmax = mstack.max(axis=0)
            
            PETmm = PenmanMonteith(relevantDataFields, tmax, tmin)
            
        if evapMethod == 'Hargreaves':
            # retrieve 3 hourly Temperature and calculate max and min Temperature       
            filename = 'Tair_E2OBS_'
            standard_name = 'air_temperature'
            timestepSeconds = 10800


            tlist, timelist = get_times(currentdate,currentdate,serverroot, wrrsetroot, filename,timestepSeconds,logger )
            ncstepobj = getstep(tlist,BB,standard_name,timestepSeconds,logger)
            #ncstepobj = getstepdaily(tlist,BB,standard_name,logger)
    
            mstack = ncstepobj.getdates(timelist)   
            latitude = ncstepobj.lat[:]
            #assuming a resolution of 0.5 degrees
            LATITUDE = np.ones(((2*(latmax-latmin)),(2*(lonmax-lonmin))))
            for i in range (0,int((2*(lonmax-lonmin)))):
                LATITUDE[:,i]=LATITUDE[:,i]*latitude
            tmin = mstack.min(axis=0)            
            tmax = mstack.max(axis=0)

            PETmm, Ra, dst, angle, dec = hargreaves(LATITUDE,currentdate,relevantDataFields, tmax, tmin)
            dst = dst * 180.0/pi
            save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,Ra,int(ncnt),odir,prefix="RA",oformat=oformat)
            save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,dst,int(ncnt),odir,prefix="DST",oformat=oformat)
            save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,angle,int(ncnt),odir,prefix="ANG",oformat=oformat)
            #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,dec,int(ncnt),odir,prefix="DEC",oformat=oformat)

        logger.info("Saving PET data for: " +str(currentdate))
        #save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,PETmm[0],int(ncnt),odir,prefix=oprefix,oformat=oformat)  
        save_as_mapsstack_per_day(ncstepobj.lat,ncstepobj.lon,PETmm[0],int(ncnt),odir,prefix=oprefix,oformat=oformat)  
        

        currentdate += datetime.timedelta(days=1)
        ncnt +=1        
    
    logger.info("Done.")

if __name__ == "__main__":
    main()



