"""
Get a variable from the forcing data from the e2o server for a specific region and time range


usage:

    e2o_getvar.py -I inifile [-l loglevel][-h]

    -I inifile - ini file with settings which data to get
    -l loglevel (must be one of DEBUG, WARNING, ERROR)

    See: http://e2o-downscaling-tools.readthedocs.io/en/latest/e2o_getvar.html
"""


import getopt, sys, os
import datetime
from numpy import *
from e2o_dstools.e2o_utils import *
import e2o_dstools
import math

def usage(*args):
    """
    Print usage information

    -  *args: command line arguments given
    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)







def save_as_mapsstack(lat,lon,data,times,directory,prefix="E2O",oformat="PCRaster"):
    """
    Save a data matrix (multiple times) as a stack of (pcraster) maps.

    :param lat:
    :param lon:
    :param data:
    :param times:
    :param directory:
    :param prefix:
    :param oformat:
    :return:
    """
    
    cnt = 0
    if not os.path.exists(directory):
        os.mkdir(directory)
    for a in times:
            below_thousand = cnt % 1000
            above_thousand = cnt / 1000
            mapname  = str(prefix + '%0' + str(8-len(prefix)) + '.f.%03.f') % (above_thousand, below_thousand)
            print "saving map: " + os.path.join(directory,mapname)
            writeMap(os.path.join(directory,mapname),oformat,lon,lat[::-1],flipud(data[cnt,:,:]),-999.0)
            cnt = cnt + 1    



def dict_split(d, chunk_size=1):
    return [
            dict(item for item in sorted(d.items())[i:i+chunk_size])
            for i in range(0, len(d.items()), chunk_size)
           ]


def downscale(variable, data, datalow, wrrversion,serverroot,wrrsetroot,BB,currentdate,XH,YH,XL,YL,
              interpolmethod, ncoutfillval, hiresdem, resLowResDEM, interpolmode, logger):
    """

    :param variable:
    :param data:
    :param wrrversion:
    :return:
    """

    if variable == 'Temperature':
        if wrrversion == 2:
            standard_name = 'air_temperature_lapse_rate'
            tlist, timelist = get_times_daily(currentdate, currentdate, serverroot, wrrsetroot,
                                              "lapseM_EI_025_", logger)

            ncstepobj = getstepdaily(tlist, BB, standard_name, logger)
            mmstack = ncstepobj.getdates(timelist)
            lapse_rate = flipud(mmstack.mean(axis=0))
            lapse_rate = resample_grid(lapse_rate, ncstepobj.lon, ncstepobj.lat, XH, YH,
                                       method=interpolmethod, FillVal=ncoutfillval)
        else:
            lapse_rate = -0.006
        retdata = data + lapse_rate * (hiresdem - resLowResDEM)
    elif variable == "Rainfall":
        if wrrversion == 2:
            # First read P data in WRR2 resolution
            yday = (currentdate.date() - datetime.date(currentdate.year, 1, 1)).days + 1
            PclimMapName = getmapnamemonth(yday, "prec")
            PclimWRR2 = e2o_dstools.get_data(os.path.join('Prec/0.2500/', PclimMapName))
            resX, resY, cols, rows, CLIMLON, CLIMLAT, PlowClim, FillVal = \
                readMap(PclimWRR2, 'GTiff', logger)
            PlowClim = PlowClim.astype(float32)
            PlowClim[PlowClim == FillVal] = NaN

            # Now read the P climate data for the current resolution
            Resstr = "%1.4f" % diff(XH).mean() # find current resolution
            PclimCUR = e2o_dstools.get_data(os.path.join('Prec',Resstr,PclimMapName))
            resX, resY, cols, rows, HILON, HILAT, PHiClim, FillVal = \
                readMap(PclimCUR, 'GTiff', logger)
            PHiClim = PHiClim.astype(float32)
            PHiClim[PHiClim==FillVal] = NaN

            # Now determine diff between current data at original resolution and the climatology at the same resolutions
            # resample climatology first as it may not cover the whole earth
            _PlowClim= resample_grid(PlowClim,  CLIMLON, CLIMLAT,XL, YL,
                                       method='nearest', FillVal=NaN)

            datalow[datalow<=0] = NaN
            multlow = _PlowClim/datalow
            multhi = resample_grid(multlow, XL, YL, XH, YH,
                                       method=interpolmode, FillVal=NaN)

            _PHiClim = resample_grid(PHiClim,  HILON, HILAT,XH, YH,
                                       method='nearest', FillVal=NaN)
            multhi[multhi <= 0] = NaN
            retdata = _PHiClim/multhi
            # fill with interpolated original data
            retdata[~isfinite(retdata)] = data[~isfinite(retdata)]
        else:
            retdata = data
    else:
        retdata = data


    return retdata





def main(argv=None):

    serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
    wrrsetroot = "ecmwf/met_forcing_v0/"
    variable = "Tair_daily_E2OBS_"
    # defaults, overwritten by info from ini file
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
    getDataForVar = False
    loglevel=logging.INFO
    downscaling = "False"
    resampling = "True"
    metadata = {}
    

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    try:
        opts, args = getopt.getopt(argv, 'I:l:h')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-I': inifile = a
        if o == '-h': usage()
        if o == '-l': exec "loglevel = logging." + a


    logger = setlogger("e2o_getvar.log","e2o_getvar",level=loglevel)
    logger.debug("Reading settings from in: " + inifile)
    theconf = iniFileSetUp(inifile)
    interpolmethod ='linear'
    ncoutfillval=-9999.0


    lonmax = float(configget(logger,theconf,"selection","lonmax",str(lonmax)))
    lonmin = float(configget(logger,theconf,"selection","lonmin",str(lonmin)))
    latmax = float(configget(logger,theconf,"selection","latmax",str(latmax)))
    latmin = float(configget(logger,theconf,"selection","latmin",str(latmin)))

    BB = dict( lon=[ lonmin, lonmax], lat= [ latmin, latmax])

    startyear = int(configget(logger,theconf,"selection","startyear",str(startyear)))
    endyear = int(configget(logger,theconf,"selection","endyear",str(endyear)))
    endmonth = int(configget(logger,theconf,"selection","endmonth",str(endmonth)))
    startmonth = int(configget(logger,theconf,"selection","startmonth",str(startmonth)))
    endday = int(configget(logger,theconf,"selection","endday",str(endday)))
    startday = int(configget(logger,theconf,"selection","startday",str(startday)))
    serverroot = configget(logger,theconf,"url","serverroot",serverroot)
    wrrsetroot = configget(logger,theconf,"url","wrrsetroot",wrrsetroot)

    start = datetime.datetime(startyear,startmonth,startday)
    end = datetime.datetime(endyear,endmonth,endday)

    oformat = configget(logger,theconf,"output","format","PCRaster")
    oodir = configget(logger,theconf,"output","directory","output/")
    oprefix = configget(logger,theconf,"output","prefix","E2O")
    resampling  = configget(logger,theconf,"selection","resampling",resampling)
    FNhighResDEM = configget(logger,theconf,"downscaling","highResDEM","downscaledem.map")
    netcdfout = configget(logger, theconf, "output", "netcdfout", "None")
    rwbuffer = int(configget(logger, theconf, "output", "rwbuffer", "10"))
    ncoutfillval = float(configget(logger, theconf, "output", "ncoutfillval", "-9999.9"))
    least_significant_digit = int(configget(logger, theconf, "output", "least_significant_digit", "6"))
    interpolmethod = configget(logger, theconf, "downscaling", "interpolmethod", interpolmethod)
    downscaling = configget(logger, theconf, "downscaling", "downscaling", downscaling)
    wrrversion = int(configget(logger, theconf, "selection", "wrrversion", '2'))
    timestepsecs = float(configget(logger, theconf, "selection", "timestepsecs", '86400'))
    timeflatten =configget(logger, theconf, "conversion", "intime", 'mean')

    logger.debug("Done reading settings.")

    if  wrrversion ==2:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR2.tif')
    else:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR1.tif')



    variable = configget(logger, theconf, "selection", "variable", 'None')
    filename = configget(logger, theconf, "selection", "filename", 'None')
    standard_name = configget(logger, theconf, "selection", "standard_name", 'None')



    FNlowResDEM = configget(logger, theconf, "downscaling", "lowResDEM", FNlowResDEM)

    if downscaling =="True" or resampling == "True":
        resX, resY, cols, rows, xhires, yhires, hiresdem, FillVal = readMap(FNhighResDEM,'PCRaster',logger)
        resX, resY, cols, rows, xlres, ylres, lowresdem, FillVal = readMap(FNlowResDEM, 'PCRaster', logger)
        # Resample orid dem to new resolution using nearest
        x = xhires
        y = yhires
        demmask = hiresdem != FillVal
        mismask = hiresdem == FillVal
        Ldemmask = lowresdem != FillVal
        Lmismask = lowresdem == FillVal
        # Fille gaps in high res DEM with Zeros for ineterpolation purposes
        lowresdem[Lmismask] = 0.0
        resLowResDEM = resample_grid(lowresdem, xlres, ylres, xhires, yhires, method=interpolmethod,
                                     FillVal=FillVal)

        lowresdem[Lmismask] = FillVal
        BB = dict(lon=[min(x), max(x)], lat=[min(y), max(y)])
    else:
        resX, resY, cols, rows, xlres, ylres, lowresdem, FillVal = readMap(FNlowResDEM, 'PCRaster', logger)
        x = xlres
        y = ylres

    EndStep = (end - start).days + 1
    StartStep = 1



    # Setup netcdf output plus meta information
    if netcdfout != 'None':
        try:
            globalmetadata = getmetadatafromini(inifile,'netcdf_attributes')
            metadata.update(globalmetadata)
        except:
            logger.warn("No netcdf metadata found in ini file. Expected section: netcdf_attributes")
        ncout = netcdfoutput(netcdfout,x,y,logger,start,EndStep - StartStep,metadata=metadata,maxbuf=rwbuffer,
                             least_significant_digit=least_significant_digit,FillVal=ncoutfillval)


    #Add options for multiple variables
    #for i in range (0,len(variables)):
    getDataForVar = False
    try:
        varmetadata = getmetadatafromini(inifile, 'netcdf_attributes_' + variable)
    except:
        varmetadata = {}
        logger.warn("Could not find section with variable metadata: netcdf_attributes_" + variable)


    valid_max = float(configget(logger, theconf, 'selection', "valid_max", '1E31'))
    valid_min = float(configget(logger, theconf, 'selection', "valid_min", '-1E31'))

    start = datetime.datetime(startyear, startmonth, startday)
    end = datetime.datetime(endyear, endmonth, endday)

    cnt = 0
    odir = os.path.join(oodir,variable)
    if not os.path.exists(odir):
        os.makedirs(odir)

    allsteps=(end-start).days
    #for tlist,timelist in zip(chunks,lchunks):
    for thisstep in arange(0,allsteps +1):
        currentdate=start+datetime.timedelta(days=thisstep)
        logger.info("Processing date: " + str(currentdate))
        if variable == "Rainfall":
            tlist, timelist = get_times_P(currentdate, currentdate, serverroot, wrrsetroot, filename,
                                        timestepsecs, logger)
        else:
            tlist, timelist = get_times(currentdate, currentdate, serverroot, wrrsetroot, filename,
                                        timestepsecs, logger)
        ncstepobj = getstep(tlist,BB,standard_name,timestepsecs,logger)
        # get the steps for this time


        mstack = ncstepobj.getdates(timelist)
        exec "thevar = mstack." + timeflatten + "(axis=0)"

        thevar[thevar<valid_min]=NaN
        thevar[thevar > valid_max] = NaN
        arcnt = 0

        mapname = getmapname(cnt+1,oprefix)
        convstr = configget(logger, theconf, "conversion", variable, 'none')
        if convstr != 'none':
            convstr = convstr.replace(variable, 'thevar')
            try:
                exec "thevar =  " + convstr
            except:
                logger.error("Conversion string not valid: " + convstr)

        if resampling == "True":
            newdata = resample_grid(flipud(thevar),ncstepobj.lon,ncstepobj.lat, xhires,
                                    yhires,method=interpolmethod,FillVal=NaN)



            # Process temperature, downscale and use laps rate if possible
            # needs to be fixed
            if downscaling == 'True':
                newdata = downscale(variable, newdata,flipud(thevar),wrrversion,serverroot,wrrsetroot,BB,currentdate, xhires, yhires,ncstepobj.lon,ncstepobj.lat,
                              interpolmethod, ncoutfillval, hiresdem, resLowResDEM,interpolmethod,logger)




            newdata[isnan(newdata)] = ncoutfillval
            newdata[~isfinite(newdata)] = ncoutfillval
            if netcdfout != 'None':
                logger.info("Saving step to netcdf: " + str(thisstep))
                ncout.savetimestep(cnt + 1, newdata, name=standard_name, var=variable, metadata=varmetadata)
            else:
                logger.info("Writing map: " + os.path.join(odir, mapname))
                writeMap(os.path.join(odir,mapname),oformat,xhires,yhires,newdata,ncoutfillval)
        else:
            newdata = flipud(thevar).copy()
            if convstr != 'none':
                convstr = convstr.replace(variable,'newdata')
                try:
                    exec "newdata =  " + convstr
                except:
                    logger.error("Conversion string not valid: " + convstr)

            if netcdfout != 'None':
                logger.info("Saving step to netcdf: " + str(thisstep))
                ncout.savetimestep(cnt + 1, newdata, name=standard_name, var=variable, metadata=varmetadata)
            else:
                logger.info("Writing map: " + os.path.join(odir, mapname))
                writeMap(os.path.join(odir,mapname),oformat,ncstepobj.lon,ncstepobj.lat[::-1],newdata,ncoutfillval)


        cnt = cnt + 1
        arcnt = arcnt +1
        del ncstepobj


    logger.info("Done.")



if __name__ == "__main__":
    main()



