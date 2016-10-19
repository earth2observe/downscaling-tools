# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 07:36:06 2014

@author: schelle
"""

import osgeo
import osgeo.ogr
import pyproj

from osgeo import gdal, gdalconst
import logging, netCDF4
import logging.handlers
import ConfigParser
import os
from numpy import *
import sys
import datetime
import datetime as dt
import time
from scipy import interpolate

try:
    import netCDF4.utils
except:
    import netCDF4_utils

import netcdftime



globmetadata = {}
globmetadata['title'] = 'e2o_dstools output'
globmetadata['institution'] = 'Deltares'
globmetadata['source'] = 'eartH2Observe'
globmetadata['history'] = time.ctime()
globmetadata['references'] = 'https://github.com/earth2observe'
globmetadata['Conventions'] = 'CF-1.4'



class ncdatset():
    """
    Wrapper around the nc object to simplify things

    Opens the dataset and determines the number of dimensions
    3 = T, Lat Lon
    4 = T heigth, Lat Lon
    """
    def __init__(self,ncurl,logger):
        self.logger = logger
        try:
            self.nc = netCDF4.Dataset(ncurl)
        except:
            self.logger.error("Failed to open remote file: " + ncurl)
            sys.exit(12)

        self.lat = self.getlat(self.nc)
        if self.lat == None:
            self.logger.error("No lat information found!")
        self.lon = self.getlon(self.nc)
        if self.lon == None:
            self.logger.error("No lon information found!")

        self.heigth = self.getheigth(self.nc)
        self.dimensions = 3 if self.heigth == None else 4
        self.time =   self.gettime(self.nc)
        #self.timesteps = self.time.shape[0]
        self.logger.debug(self.nc)


    def __del__(self):
        self.nc.close()


    def getlat(self,ncdataset):
        """
        Check for standard name latitude. If that is not present assume the variable name
        """

        for a in self.nc.variables:
            if hasattr(self.nc.variables[a],'standard_name'):
                if  self.nc.variables[a].standard_name == 'latitude':
                    return self.nc.variables[a]
            else:
                if a == 'latitude':
                    return self.nc.variables[a]
                if a == 'lat':
                    return self.nc.variables[a]

        return None

    def getvarbyname(self,name):
        """
        Check for standard name. If that is not present assume the variable name
        """

        for a in self.nc.variables:
            if hasattr(self.nc.variables[a],'standard_name'):
                if  self.nc.variables[a].standard_name == name:
                    return self.nc.variables[a]
            else:
                if  a  == name:
                    return self.nc.variables[a]


        return None

    def gettime(self,ncdataset):
        """
        Check for standard name. If that is not present assume the variable name
        """

        for a in ncdataset.variables:
            if hasattr(self.nc.variables[a],'standard_name'):
                if 'time' in ncdataset.variables[a].standard_name:
                    return ncdataset.variables[a]
            elif hasattr(self.nc.variables[a],'units'):
                if 'since' in ncdataset.variables[a].units:
                    return ncdataset.variables[a]

        return None


    def getlon(self,ncdataset):
        """
        Check for standard name. If that is not present assume the variable name
        """

        for a in self.nc.variables:
            if hasattr(self.nc.variables[a],'standard_name'):
                if  self.nc.variables[a].standard_name == 'longitude':
                    return self.nc.variables[a]
            else:
                if a == 'longitude':
                    return self.nc.variables[a]
                if a == 'lon':
                    return self.nc.variables[a]

        return None


    def getheigth(self,ncdataset):
        """
        """

        for a in self.nc.variables:
            if hasattr(self.nc.variables[a],'standard_name'):
                if  self.nc.variables[a].standard_name == 'height':
                    return self.nc.variables[a]
            else:
                if a == 'height':
                    return self.nc.variables[a]

        return None



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
        flipped = 0
        ret = []

        lastnc = None

        # here we loop over nc files for speed reasons

        for theone in  unique(self.list.values()):
            self.dset = ncdatset(theone,self.logger)

            time = self.dset.time
            tar = time[:]
            try:
                timeObj = netCDF4.num2date(tar, units=time.units, calendar=time.calendar)
            except:
                timeObj = netCDF4.num2date(tar, units=time.units, calendar='standard')
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


            self.logger.info("Processing url: " + theone )

            data = self.dset.getvarbyname(self.varname)

            if data == None:
                self.logger.error("dataset with standard_name " + self.varname + " not found" )


            if self.dset.lat[0] > self.dset.lat[-1]:
                lat = self.dset.lat[::-1]
                flipped = 1
            else:
                lat = self.dset.lat[:]
            lon = self.dset.lon[:]

            (self.latidx,) = logical_and(lat >= self.BB['lat'][0], lat <= self.BB['lat'][1]).nonzero()
            (self.lonidx,) = logical_and(lon >= self.BB['lon'][0], lon <= self.BB['lon'][1]).nonzero()

            if self.dset.dimensions ==3:
                if 'lat' in data.dimensions[0]:

                    window = zeros((data.shape[2],self.latidx.max() - self.latidx.min() + 1,self.lonidx.max() - self.lonidx.min() + 1))
                    for i in arange(spos,epos):
                        tt = flipud(data[:,:,i])
                        window[i,:,:] = tt[self.latidx.min():self.latidx.max() + 1,
                                 self.lonidx.min():self.lonidx.max() + 1]
                else:
                    window = data[spos:epos,self.latidx.min():self.latidx.max()+1,self.lonidx.min():self.lonidx.max()+1]

                if flipped:
                    window = flipud(window).copy()

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


            self.logger.info("Processing url: " + theone )

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
        self.timestepSeconds = timestepSeconds

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

            if self.timestepSeconds < 3600:
                dpos = thedate.second -1
            elif self.timestepSeconds < 86400:
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


            self.logger.info("Processing url: " + theone )

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


def prepare_nc(trgFile, timeList, x, y, metadata, logger, EPSG="EPSG:4326", units=None,
               calendar='gregorian', Format="NETCDF4", complevel=9, zlib=True, least_significant_digit=None):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 1900-01-01 00:00:00'
    """
    import datetime as dt

    logger.info('Setting up netcdf output: ' + trgFile)

    if units == None: # Use start of the run
        epoch = timeList[0]
        units = 'seconds since %04d-%02d-%02d %02d:%02d:%02d.0 00:00' % (
        epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second)

    startDayNr = netCDF4.date2num(timeList[0].replace(tzinfo=None), units=units, calendar=calendar)
    endDayNr = netCDF4.date2num(timeList[-1].replace(tzinfo=None), units=units, calendar=calendar)

    timeAR = linspace(startDayNr, endDayNr, num=len(timeList))

    nc_trg = netCDF4.Dataset(trgFile, 'w', format=Format, zlib=zlib, complevel=complevel)

    logger.info(
        'Setting up dimensions and attributes. Steps: ' + str(len(timeList)) + ' lat: ' + str(len(y)) + " lon: " + str(
            len(x)))
    if len(timeAR) == 1:
        nc_trg.createDimension('time', 1)
    else:
        nc_trg.createDimension('time', 0)  # NrOfDays*8

    DateHour = nc_trg.createVariable('time', 'f8', ('time',), fill_value=-9999., zlib=zlib, complevel=complevel)
    DateHour.units = units
    DateHour.calendar = calendar
    DateHour.standard_name = 'time'
    DateHour.long_name = 'time'
    DateHour.axis = 'T'
    DateHour[:] = timeAR

    # make a proj4 string
    srs = osgeo.osr.SpatialReference()
    res = srs.ImportFromEPSG(int(EPSG[5:]))
    if res != 0:
        logger.error("EPGS not converted correctly: " + EPSG + ". Is the GDAL_DATA environment variable set correctly?")
        exit(1)

    projStr = srs.ExportToProj4()
    proj_src = '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'


    if srs.IsProjected() == 0: # ONly lat lon needed
        nc_trg.createDimension('lat', len(y))
        nc_trg.createDimension('lon', len(x))
        y_var = nc_trg.createVariable('lat', 'f4', ('lat',), fill_value=-9999., zlib=zlib, complevel=complevel)
        y_var.standard_name = 'latitude'
        y_var.long_name = 'latitude'
        y_var.units = 'degrees_north'
        y_var.axis = 'Y'
        x_var = nc_trg.createVariable('lon', 'f4', ('lon',), fill_value=-9999., zlib=zlib, complevel=complevel)
        x_var.standard_name = 'longitude'
        x_var.long_name = 'longitude'
        x_var.units = 'degrees_east'
        x_var.axis = 'X'
        y_var[:] = y
        x_var[:] = x
        crs = nc_trg.createVariable('crs', 'c')
        crs.long_name = 'wgs84'
        crs.proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        crs.grid_mapping_name = 'latitude_longitude'
    else:  # Assume regular grid in m
        nc_trg.createDimension('y', len(y))
        nc_trg.createDimension('x', len(x))
        y_var = nc_trg.createVariable('y', 'f4', ('y',), fill_value=-9999., zlib=zlib, complevel=complevel)
        y_var.standard_name = 'projection_y_coordinate'
        y_var.long_name = 'y-coordinate in Cartesian system'
        y_var.units = 'm'
        y_var.axis = 'Y'
        x_var = nc_trg.createVariable('x', 'f4', ('x',), fill_value=-9999., zlib=zlib, complevel=complevel)
        x_var.standard_name = 'projection_x_coordinate'
        x_var.long_name = 'x-coordinate in Cartesian system'
        x_var.units = 'm'
        x_var.axis = 'X'
        y_var[:] = y
        x_var[:] = x
        crs = nc_trg.createVariable('crs', 'c')
        crs.long_name = EPSG
        crs.grid_mapping_name = 'universal_transverse_mercator'
        crs.utm_zone_number = srs.GetUTMZone()
        crs.semi_major_axis = srs.GetSemiMajor()
        crs.inverse_flattening = srs.GetInvFlattening()
        crs._CoordinateTransformType = "Projection"
        crs._CoordinateAxisTypes = "y x"
        crs.proj4_params = projStr
        # Also write lat lon fields
        XI,YI = meshgrid(x,y)
        lon_vals,lat_vals = convertCoord(projStr, proj_src, XI, YI)
        # Need to create lat-lon fields
        lat = nc_trg.createVariable('lat','f4',('y','x',))
        lat.standard_name = 'latitude'
        lat.long_name = 'latitude coordinate'
        lat.units = 'degrees_north'
        lat.coordinates = 'lat lon'
        lat.grid_mapping = 'wgs84'
        #lat._CoordinateAxisType = "Lat"
        lat[:,:] = lat_vals
        lon = nc_trg.createVariable('lon','f4',('y','x',))
        lon.standard_name = 'longitude'
        lon.long_name = 'longitude coordinate'
        lon.units = 'degrees_east'
        lon.coordinates = 'lat lon'
        lon.grid_mapping = 'wgs84'
        #lon._CoordinateAxisType = "Lon"
        lon[:,:] = lon_vals

    crs.EPSG_code = EPSG

    # now add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()



class netcdfoutput():
    def __init__(self, netcdffile, x, y, logger, starttime, timesteps, EPSG="EPSG:4326", timestepsecs=86400,
                 metadata={}, zlib=True, Format="NETCDF4",
                 maxbuf=25, least_significant_digit=None):
        """
        Under construction
        """

        self.EPSG = EPSG
        self.zlib = zlib
        self.Format = Format
        self.least_significant_digit = least_significant_digit

        def date_range(start, end, timestepsecs):
                r = int((end + dt.timedelta(seconds=timestepsecs) - start).total_seconds()/timestepsecs)
                return [start + dt.timedelta(seconds=(timestepsecs * i)) for i in range(r)]

        self.logger = logger
        # Do not allow a max buffer larger than the number of timesteps
        self.maxbuf = maxbuf if timesteps >= maxbuf else timesteps
        self.ncfile = netcdffile
        self.timesteps = timesteps

        # Shift one timestep as we output at the end
        #starttime = starttime + dt.timedelta(seconds=timestepsecs)
        end = starttime + dt.timedelta(seconds=timestepsecs * (self.timesteps))

        timeList = date_range(starttime, end, timestepsecs)
        self.timestepbuffer = zeros((self.maxbuf, len(y), len(x)))
        self.bufflst = {}

        globmetadata.update(metadata)

        prepare_nc(self.ncfile, timeList, x, y, globmetadata, logger, Format=self.Format, EPSG=EPSG,zlib=self.zlib,
                   least_significant_digit=self.least_significant_digit)

    def savetimestep(self, timestep, data, unit="mm", var='P', name="Precipitation"):
        """
        save a single timestep for a variable

        input:
            - timestep - current timestep
            - data - pcraster map to save
            - unit - unit string
            - var - variable string
            - name - name of the variable
        """
        # Open target netCDF file
        var = os.path.basename(var)
        self.nc_trg = netCDF4.Dataset(self.ncfile, 'a', format=self.Format, zlib=self.zlib, complevel=9)
        self.nc_trg.set_fill_off()
        # read time axis and convert to time objects
        # TODO: use this to append time
        # time = self.nc_trg.variables['time']
        # timeObj = netCDF4.num2date(time[:], units=time.units, calendar=time.calendar)

        idx = timestep - 1

        buffreset = (idx + 1) % self.maxbuf
        bufpos = (idx) % self.maxbuf

        try:
            nc_var = self.nc_trg.variables[var]
        except:
            self.logger.debug("Creating variable " + var + " in netcdf file. Format: " + self.Format)
            if self.EPSG.lower() == "epsg:4326":
                nc_var = self.nc_trg.createVariable(var, 'f4', ('time', 'lat', 'lon',), fill_value=-9999.0, zlib=self.zlib,
                                                    complevel=9, least_significant_digit=self.least_significant_digit)
                nc_var.coordinates = "lat lon"
            else:
                nc_var = self.nc_trg.createVariable(var, 'f4', ('time', 'y', 'x',), fill_value=-9999.0, zlib=self.zlib,
                                                    complevel=9, least_significant_digit=self.least_significant_digit)
                nc_var.coordinates = "lat lon"
                nc_var.grid_mapping = "crs"

            nc_var.units = unit
            nc_var.standard_name = name
            self.nc_trg.sync()

        miss = float(nc_var._FillValue)

        if self.bufflst.has_key(var):
            self.bufflst[var][bufpos, :, :] = data
        else:
            self.bufflst[var] = self.timestepbuffer.copy()
            self.bufflst[var][bufpos, :, :] = data

        # Write out timestep buffer.....

        if buffreset == 0 or idx == self.maxbuf - 1 or self.timesteps <= timestep:
            spos = idx - bufpos
            self.logger.debug(
                "Writing buffer for " + var + " to file at: " + str(spos) + " " + str(int(bufpos) + 1) + " timesteps")
            nc_var[spos:idx + 1, :, :] = self.bufflst[var][0:bufpos + 1, :, :]
            self.nc_trg.sync()

    def finish(self):
        """
        Flushes and closes the netcdf file

        :return: Nothing
        """
        if hasattr(self, "nc_trg"):
            self.nc_trg.sync()
            self.nc_trg.close()



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
        #ncfile = serverroot + wrrsetroot + filename + "%d" % (thedate.year) + ".nc"
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



def readMap(fileName, fileFormat,logger):
    """
    ead geographical file into memory

    :param fileName: Name of the file to read
    :param fileFormat: Gdal format string
    :param logger: logger object
    :return: resX, resY, cols, rows, x, y, data, FillVal
    """


    #Open file for binary-reading

    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        logger.error('Could not open ' + fileName + '. Something went wrong!! Shutting down')
        sys.exit(1)
    # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    del ds
    return resX, resY, cols, rows, x, y, data, FillVal



def correctTemp(Temp,LapseRate,elevationCorrection):

    """
    Elevation based correction of temperature

    inputs:
    Temperature             = daily mean, min or max temperature (degrees Celcius)
    Elevation correction    = difference between high resolution and low resolution (0.5 degrees) DEM  [m]


    """

    #apply elevation correction

    Temp_cor   = Temp + LapseRate * elevationCorrection

    return Temp_cor


def getmapnamemonth(yearday,prefix):
    """
    generate a pcraster type mapname (month) based on a yearday as input
    :var number: number of the mape
    :var prefix: prefix for the map

    :return: Name
    """
    month = (datetime.date(2000,1,1) + datetime.timedelta(days=yearday)).timetuple().tm_mon

    below_thousand = month % 1000
    above_thousand = month / 1000

    realprefix =  os.path.basename(prefix)



    mapname = os.path.join(os.path.dirname(prefix),
                           str(realprefix + '%0' + str(8-len(realprefix)) + '.f.%03.f') % (above_thousand, below_thousand))

    return mapname

def getmapname(number,prefix):
    """
    generate a pcraster type mapname based on timestep and prefix
    :var number: number of the mape
    :var prefix: prefix for the map

    :return: Name
    """

    below_thousand = number % 1000
    above_thousand = number / 1000
    mapname = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)

    return mapname

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """
    Write geographical data into file. Also replace NaN by FillVall

    :param fileName:
    :param fileFormat:
    :param x:
    :param y:
    :param data:
    :param FillVal:
    :return:
    """


    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    data[isnan(data)] = FillVal
    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
        print "Output format: " + fileFormat
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2

    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data.astype(float32),0,0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(FillVal)

    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    if fileFormat == 'GTiff':
        outDataset = driver2.CreateCopy(fileName, TempDataset, 0 ,options = ['COMPRESS=LZW'])
    else:
        outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'



def configget(log,config,section,var,default):
    """   
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value 
    with the default in the config-file
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default string
        
    Returns:
        - string - either the value from the config file or the default value
    """
    

    try:
        ret = config.get(section, var)
    except:
        ret = default
        log.info( "returning default (" + str(default) + ") for " + section + ":" + var)
        configset(config,section,var,str(default), overwrite=False)
    

    return ret       

def configset(config,section,var,value, overwrite=False):
    """   
    Sets a string in the in memory representation of the config object
    Deos NOT overwrite existing values if overwrite is set to False (default)
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to set
        - value - the value to set
        - overwrite (optional, default is False)
   
    Returns:
        - nothing
        
    """
    
    if not config.has_section(section):
        config.add_section(section)
        config.set(section,var,value)
    else:     
        if not config.has_option(section,var):
            config.set(section,var,value)
        else:
            if overwrite:
                config.set(section,var,value)

def iniFileSetUp(configfile):
    """
    Reads .ini file and sets default values if not present
    """
    # TODO: clean up wflwo specific stuff
    #setTheEnv(runId='runId,caseName='caseName)
    # Try and read config file and set default options
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config

def setlogger(logfilename,loggername, level=logging.INFO):
    """
    Set-up the logging system and return a logger object. Exit if this fails
    """

    try:
        #create logger
        logger = logging.getLogger(loggername)
        if not isinstance(level, int):
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(level)
        ch = logging.FileHandler(logfilename,mode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(2)

def closeLogger(logger, ch):
    """
    Closes the logger
    """
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch


def getmapname(number,prefix):
    """
    generate a pcraster type mapname based on timestep and prefix
    :var number: number of the mape
    :var prefix: prefix for the map

    :return: Name
    """
    below_thousand = number % 1000
    above_thousand = number / 1000
    mapname = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)

    return mapname


def resample_grid(gridZ_in,Xin,Yin,Xout,Yout,method='nearest',FillVal=1E31):
    """
    Resample a regular grid be supplying original and new x, y row,col coordinates
    Missing values is set to 1E31, the PCRaster standard. If the output grid
    has a lower resolution the interpolation is done in two steps by first applying
    a moving window.

    :param gridZ_in: datablock of original data (e.g. from readMap)
    :param Xin: X-coordinates of all columns (e.g. from readMap)
    :param Yin: Y-coordinates of all columns (e.g. from readMap)
    :param Xout: X-coordinates of all columns in the new map (e.g. from readMap)
    :param Yout: Y-coordinates of all columns in the new map (e.g. from readMap)
    :param method: linear, nearest, [cubic, quintic <- switched off for now]
    :param FillVal: Value for the missing values

    :return: datablock of new grid.
    """
    
    # we need to sort the y data (must be ascending)
    # and thus flip the image
    Ysrt = sort(Yin)

    # First average if the output has a lower resolution than the input grid
    insize = min(diff(Xin))
    outsize = min(diff(Xout))
    #if insize < outsize:
    #    xsize = outsize/insize
    #    gridZ_in = ndimage.filters.percentile_filter(gridZ_in,50,size=xsize)

    if method in 'nearest linear':
        interobj = interpolate.RegularGridInterpolator((Ysrt,Xin), flipud(gridZ_in), method=method ,bounds_error=False,fill_value=float32(FillVal))
        _x, _y = meshgrid(Xout, Yout)
        yx_outpoints = transpose([_y.flatten(), _x.flatten()])

        # interpolate
        res = interobj(yx_outpoints)
        result = reshape(res, (len(Yout), len(Xout)))
        result[isnan(result)] = FillVal
    # elif method in 'cubic quintic':
    #     interobj = interpolate.interp2d(Xin, Ysrt, gridZ_in,  kind=method, bounds_error=False)
    #     res = interobj(Xout, Yout)
    #     result = reshape(res, (len(Yout), len(Xout)))
    #     result[isnan(result)] = FillVal
    #     res = flipud(res)
    else:
        raise ValueError("Interpolation method " + method + " not known.")

    #reshape to the new grid

    
    return result
