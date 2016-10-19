"""
Get a variable from the forcing data from the e2o server for a specific region and time range


usage:

    e2o_getvar.py -I inifile [-l loglevel][-h]

    -I inifile - ini file with settings which data to get
    -l loglevel (must be one of DEBUG, WARNING, ERROR)
"""


import getopt, sys, os
import datetime
from numpy import *
from e2o_dstools.e2o_utils import *
import e2o_dstools


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




def main(argv=None):


    serverroot = "http://wci.earth2observe.eu/thredds/dodsC/"
    wrrsetroot = "ecmwf/met_forcing_v0/"
    variable = "Tair_daily_E2OBS_"
    
    #available variables with corresponding file names and standard_names as in NC files
    variables = ['Temperature','DownwellingLongWaveRadiation','SurfaceAtmosphericPressure',\
                    'NearSurfaceSpecificHumidity','Rainfall','SurfaceIncidentShortwaveRadiation','SnowfallRate','NearSurfaceWindSpeed']
    filenameswrr1 = ["Tair_daily_E2OBS_","LWdown_daily_E2OBS_","PSurf_daily_E2OBS_","Qair_daily_E2OBS_",\
                    "Rainf_daily_E2OBS_","SWdown_daily_E2OBS_","Snowf_daily_E2OBS_","Wind_daily_E2OBS_"]
    filenameswrr2 = ["Tair_daily_EI_025_", "LWdown_daily_EI_025_", "PSurf_daily_EI_025_", "Qair_daily_EI_025_", \
                 "Rainf_daily_EI_025_", "SWdown_daily_EI_025_", "Snowf_daily_EI_025_", "Wind_daily_EI_025_"]

    standard_names = ['air_temperature','surface_downwelling_longwave_flux_in_air','surface_air_pressure','specific_humidity',\
                        'rainfal_flux','surface_downwelling_shortwave_flux_in_air','snowfall_flux','wind_speed']
    
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
    interpolmethod = configget(logger, theconf, "downscaling", "interpolmethod", interpolmethod)
    downscaling = configget(logger, theconf, "downscaling", "downscaling", downscaling)

    logger.debug("Done reading settings.")

    if 'met_forcing_v1' in wrrsetroot:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR2.tif')
        filenames = filenameswrr2
    else:
        FNlowResDEM     = e2o_dstools.get_data('DEM-WRR1.tif')
        filenames = filenameswrr1

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
                                     FillVal=0.0)

        lowresdem[Lmismask] = FillVal
        BB = dict(lon=[min(x), max(x)], lat=[min(y), max(y)])
    else:
        resX, resY, cols, rows, xlres, ylres, lowresdem, FillVal = readMap(FNlowResDEM, 'PCRaster', logger)
        x = xlres
        y = ylres

    EndStep = (end - start).days + 1
    StartStep = 1



    if netcdfout != 'None':
            ncout = netcdfoutput(netcdfout,x,y,logging,start,EndStep - StartStep)

    #Add options for multiple variables
    for i in range (0,len(variables)):
        getDataForVar = False
        # Check whether variable exists in ini file
        getDataForVar = configget(logger,theconf,"selection",variables[i],"False")
        # If variable is True read timeseries from file
        if getDataForVar == 'True':
            filename = filenames[i]
            standard_name = standard_names[i]

            start = datetime.datetime(startyear,startmonth,startday)
            end = datetime.datetime(endyear,endmonth,endday)
            odir = os.path.join(oodir,variables[i])
            if not os.path.exists(odir):
                os.makedirs(odir)
            tlist, timelist = get_times_daily(start,end,serverroot, wrrsetroot, filename,logger )       

            ncstepobj = getstepdaily(tlist,BB,standard_name,logger)

            #print unique(tlist.values())
            mstack = ncstepobj.getdates(timelist)

            logger.info("Saving " + ncstepobj.varname + " to mapstack " + odir + oprefix)


            cnt = 0
            for a in timelist:
                mapname = getmapname(cnt+1,oprefix)
                convstr = configget(logger, theconf, "conversion", variables[i], 'none')
                if resampling == "True":
                    newdata = resample_grid(flipud(mstack[cnt,:,:]),ncstepobj.lon,ncstepobj.lat, xhires,
                                            yhires,method=interpolmethod)
                    if convstr != 'none':
                        convstr = convstr.replace(variables[i],'newdata')
                        try:
                            exec "newdata =  " + convstr
                        except:
                            logger.error("Conversion string not valid: " + convstr)

                    # Process temperature, downscale and use laps rate if possible
                    if variables[i] == 'Temperature':
                        if 'met_forcing_v1' in wrrsetroot:
                            if downscaling == 'True':
                                standard_name = 'air_temperature_lapse_rate'
                                tlist, timelist = get_times_daily(a, a, serverroot, wrrsetroot,
                                                                  "lapseM_EI_025_", logger)

                                ncstepobj = getstepdaily(tlist, BB, standard_name, logger)
                                mmstack = ncstepobj.getdates(timelist)
                                lapse_rate = flipud(mmstack.mean(axis=0))
                                lapse_rate = resample_grid(lapse_rate, ncstepobj.lon, ncstepobj.lat, xhires, yhires,
                                                           method=interpolmethod, FillVal=FillVal)
                            else:
                                lapse_rate = -0.006
                        else:
                            if downscaling == 'True':
                                lapse_rate = -0.006

                        newdata =  newdata + lapse_rate * (hiresdem - resLowResDEM)
                    writeMap(os.path.join(odir,mapname),oformat,xhires,yhires,newdata,-999.0)
                else:
                    newdata = flipud(mstack[cnt,:,:]).copy()
                    if convstr != 'none':
                        convstr = convstr.replace(variables[i],'newdata')
                        try:
                            exec "newdata =  " + convstr
                        except:
                            logger.error("Conversion string not valid: " + convstr)
                    writeMap(os.path.join(odir,mapname),oformat,ncstepobj.lon,ncstepobj.lat[::-1],newdata,-999.0)

                if netcdfout != 'None':
                    ncout.savetimestep(cnt+1,newdata,name=standard_name,var=variables[i])
                cnt = cnt + 1


    logger.info("Done.")



if __name__ == "__main__":
    main()



