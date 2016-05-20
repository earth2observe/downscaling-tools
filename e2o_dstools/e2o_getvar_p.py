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
    variables = ['Rainfall']
    filenameswrr1 = [""]

    standard_names = ['precipitation']
    
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

    oformat = configget(logger,theconf,"output","format","PCRaster")
    oodir = configget(logger,theconf,"output","directory","output/")
    oprefix = configget(logger,theconf,"output","prefix","E2O")
    resampling  = configget(logger,theconf,"selection","resampling",resampling)
    FNhighResDEM = configget(logger,theconf,"downscaling","highResDEM","downscaledem.map")

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
        interpolmethod=configget(logger,theconf,"downscaling","interpolmethod",interpolmethod)
        # Resample orid dem to new resolution using nearest


    #Add options for multiple variables
    for i in range (0,len(variables)):
        getDataForVar = False
        # Check whether variable exists in ini file
        getDataForVar = configget(logger,theconf,"selection",variables[i],"False")
        # If variable is True read timeseries from file
        if getDataForVar == 'True':
            filename = filenames[i]
            print filename
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
                    newdata = resample_grid(flipud(mstack[cnt,:,:].astype(float32)),ncstepobj.lon,ncstepobj.lat, \
                                            xhires,yhires,method=interpolmethod,FillVal=nan)
                    if convstr != 'none':
                        convstr = convstr.replace(variables[i],'newdata')
                        try:
                            exec "newdata =  " + convstr
                        except:
                            logger.error("Conversion string not valid: " + convstr)
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

                cnt = cnt + 1
                #save_as_mapsstack(ncstepobj.lat,ncstepobj.lon,mstack,timelist,odir,prefix=oprefix,oformat=oformat)

    logger.info("Done.")



if __name__ == "__main__":
    main()



