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
              interpolmethod, ncoutfillval, hiresdem, resLowResDEM, interpolmode, logger,localdatadir=''):
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
    elif variable == "Rainfall" or 'SnowfallRate' or "TotalPrecipitation":
        if wrrversion == 2:
            # First read P data in WRR2 resolution
            yday = (currentdate.date() - datetime.date(currentdate.year, 1, 1)).days + 1
            PclimMapName = getmapnamemonth(yday, "prec")
            PclimWRR2 = e2o_dstools.get_data(os.path.join('Prec/0.2500/', PclimMapName))
            if not os.path.exists(PclimWRR2):
                logger.error('Cannot find WorlClim data for downscaling: ' + PclimWRR2)
                logger.error('Please download and process the worldclim data at resolution: ' + PclimWRR2)
                exit(1)

            resX, resY, cols, rows, CLIMLON, CLIMLAT, PlowClim, FillVal = \
                readMap(PclimWRR2, 'GTiff', logger)
            PlowClim = PlowClim.astype(float32)
            PlowClim[PlowClim == FillVal] = NaN

            # Now read the P climate data for the current resolution
            Resstr = "%1.4f" % diff(XH).mean() # find current resolution
            PclimCUR = e2o_dstools.get_data(os.path.join('Prec',Resstr,PclimMapName))
            if not os.path.exists(PclimCUR):
                userdir = os.path.join(localdatadir,'Prec_Clim_For_Downscale',Resstr)
                userfile = os.path.join(userdir,PclimMapName)
                if not os.path.exists(userfile):
                    logger.info('Cannot find WorlClim data for downscaling: ' + PclimCUR)
                    logger.info('Resampling 0.0083 to : ' + Resstr)
                    PclimMapName = getmapnamemonth(yday, "prec")
                    PclimBaseRes = e2o_dstools.get_data(os.path.join('Prec/0.0083/', PclimMapName))
                    resXX, resYX, cols_, rows, baseCLIMLON, baseCLIMLAT, baseClim, FillVal = \
                        readMap(PclimBaseRes, 'GTiff', logger)
                    PHiClim = resample_grid(baseClim,baseCLIMLON,baseCLIMLAT,XH,YH,method='nearest', FillVal=NaN)
                    HILON = XH
                    HILAT = YH
                    #exit(1)
                    if not os.path.exists(userdir):
                        logger.info('Making directory:  ' + userdir)
                        os.makedirs(userdir)
                    writeMap(userfile,"GTiff",XH,YH,PHiClim.astype(float32),1E31)
                else:
                    resX, resY, cols, rows, HILON, HILAT, PHiClim, FillVal = \
                        readMap(userfile, 'GTiff', logger)
                    PHiClim = PHiClim.astype(float32)
                    PHiClim[PHiClim == FillVal] = NaN
            else:
                resX, resY, cols, rows, HILON, HILAT, PHiClim, FillVal = \
                    readMap(PclimCUR, 'GTiff', logger)
                PHiClim = PHiClim.astype(float32)
                PHiClim[PHiClim==FillVal] = NaN

            # Now determine diff between current data at original resolution and the climatology at the same resolutions
            # resample climatology first as it may not cover the whole earth

            if YL[0] < YL[-1]:
                yas = YL[::-1]
            else:
                yas =  YL
            _PlowClim= resample_grid(PlowClim,  CLIMLON, CLIMLAT,XL, yas,
                                       method='nearest', FillVal=NaN)

            datalow[datalow<=0] = NaN
            multlow = _PlowClim/datalow
            multhi = resample_grid(multlow, XL, yas, XH, YH,
                                       method=interpolmode, FillVal=NaN)

            _PHiClim = resample_grid(PHiClim,  HILON, HILAT,XH, YH,
                                       method='nearest', FillVal=NaN)
            multhi[multhi <= 0] = NaN
            retdata = _PHiClim/multhi
            # fill with interpolated original data
            retdata[~isfinite(retdata)] = data[~isfinite(retdata)]
        else:
            # First read P data in WRR1 resolution
            yday = (currentdate.date() - datetime.date(currentdate.year, 1, 1)).days + 1
            PclimMapName = getmapnamemonth(yday, "prec")
            PclimWRR1 = e2o_dstools.get_data(os.path.join('Prec/0.5000/', PclimMapName))
            if not os.path.exists(PclimWRR1):
                logger.error('Cannot find WorlClim data for downscaling: ' + PclimWRR1)
                logger.error('Please download and process the worldclim data at resolution: ' + PclimWRR1)
                exit(1)

            resX, resY, cols, rows, CLIMLON, CLIMLAT, PlowClim, FillVal = \
                readMap(PclimWRR1, 'GTiff', logger)
            PlowClim = PlowClim.astype(float32)
            PlowClim[PlowClim == FillVal] = NaN

            # Now read the P climate data for the current resolution
            Resstr = "%1.4f" % diff(XH).mean() # find current resolution
            PclimCUR = e2o_dstools.get_data(os.path.join('Prec',Resstr,PclimMapName))
            if not os.path.exists(PclimCUR):
                userdir = os.path.join(localdatadir,'Prec_Clim_For_Downscale',Resstr)
                userfile = os.path.join(userdir,PclimMapName)
                if not os.path.exists(userfile):
                    logger.info('Cannot find WorlClim data for downscaling: ' + PclimCUR)
                    logger.info('Resampling 0.0083 to : ' + Resstr)
                    PclimMapName = getmapnamemonth(yday, "prec")
                    PclimBaseRes = e2o_dstools.get_data(os.path.join('Prec/0.0083/', PclimMapName))
                    resXX, resYX, cols_, rows, baseCLIMLON, baseCLIMLAT, baseClim, FillVal = \
                        readMap(PclimBaseRes, 'GTiff', logger)
                    PHiClim = resample_grid(baseClim,baseCLIMLON,baseCLIMLAT,XH,YH,method='nearest', FillVal=NaN)
                    HILON = XH
                    HILAT = YH
                    #exit(1)
                    if not os.path.exists(userdir):
                        logger.info('Making directory:  ' + userdir)
                        os.makedirs(userdir)
                    writeMap(userfile,"GTiff",XH,YH,PHiClim.astype(float32),1E31)
                else:
                    resX, resY, cols, rows, HILON, HILAT, PHiClim, FillVal = \
                        readMap(userfile, 'GTiff', logger)
                    PHiClim = PHiClim.astype(float32)
                    PHiClim[PHiClim == FillVal] = NaN
            else:
                resX, resY, cols, rows, HILON, HILAT, PHiClim, FillVal = \
                    readMap(PclimCUR, 'GTiff', logger)
                PHiClim = PHiClim.astype(float32)
                PHiClim[PHiClim==FillVal] = NaN

            # Now determine diff between current data at original resolution and the climatology at the same resolutions
            # resample climatology first as it may not cover the whole earth

            if YL[0] < YL[-1]:
                yas = YL[::-1]
            else:
                yas =  YL
            _PlowClim= resample_grid(PlowClim,  CLIMLON, CLIMLAT,XL, yas,
                                       method='nearest', FillVal=NaN)

            datalow[datalow<=0] = NaN
            multlow = _PlowClim/datalow
            multhi = resample_grid(multlow, XL, yas, XH, YH,
                                       method=interpolmode, FillVal=NaN)

            _PHiClim = resample_grid(PHiClim,  HILON, HILAT,XH, YH,
                                       method='nearest', FillVal=NaN)
            multhi[multhi <= 0] = NaN
            retdata = _PHiClim/multhi
            # fill with interpolated original data
            retdata[~isfinite(retdata)] = data[~isfinite(retdata)]
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

    netcdfout = configget(logger, theconf, "output", "netcdfout", "None")
    rwbuffer = int(configget(logger, theconf, "output", "rwbuffer", "10"))
    ncoutfillval = float(configget(logger, theconf, "output", "ncoutfillval", "-9999.9"))
    least_significant_digit = int(configget(logger, theconf, "output", "least_significant_digit", "6"))
    NetCDFFormat = configget(logger, theconf, "output", "netCDFFormat", "NETCDF4")
    interpolmethod = configget(logger, theconf, "downscaling", "interpolmethod", interpolmethod)
    downscaling = configget(logger, theconf, "downscaling", "downscaling", downscaling)
    wrrversion = int(configget(logger, theconf, "selection", "wrrversion", '2'))
    timestepsecs = float(configget(logger, theconf, "selection", "timestepsecs", '86400'))
    timeflatten =configget(logger, theconf, "conversion", "intime", 'mean')
    FNhighResDEM = configget(logger, theconf, "downscaling", "highResDEM", "downscaledem.map")

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
        # Get the extent from the high res DEM
        lonmin, latmin, lonmax, latmax = get_extent(FNhighResDEM,)
        BB = dict(lon=[lonmin, lonmax], lat=[latmin, latmax])

        resX, resY, cols, rows, xhires, yhires, hiresdem, FillVal = readMap(FNhighResDEM,'PCRaster',logger)
        resX, resY, cols, rows, xlres, ylres, lowresdem, FillVal = readMap(FNlowResDEM, 'PCRaster', logger)
        lonmin, latmin, lonmax, latmax = get_extent(FNhighResDEM,resX)
        BB = dict(lon=[lonmin, lonmax], lat=[latmin, latmax])
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
        #BB = dict(lon=[min(x), max(x)], lat=[min(y), max(y)])
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
                             least_significant_digit=least_significant_digit,FillVal=ncoutfillval,Format=NetCDFFormat)


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

    # Check for pseudo variable and expand
    if variable == "TotalPrecipitation":
        if wrrversion == 1:
            _variable = ['Rainfall','SnowfallRate']
            _standard_name = ['rainfall_flux','snowfall_flux']
            _filename = ['Rainf_daily_EOBS_', 'Snowf_daily_EOBS_']
        else:
            _variable = ['Rainfall', 'SnowfallRate']
            _standard_name = ['rainfal_flux', 'snowfall_flux']
            _filename = ['Rainf_daily_MSWEP_025_','Snowf_daily_MSWEP_025_']



    #for tlist,timelist in zip(chunks,lchunks):
    for thisstep in arange(0,allsteps +1):
        _tlist = []
        _timelist = []

        mapname = getmapname(cnt + 1, oprefix)

        if os.path.exists(os.path.join(odir,mapname)):
            logger.info('Skipping map: ' + mapname)
        else:
            currentdate=start+datetime.timedelta(days=thisstep)
            logger.info("Processing date: " + str(currentdate))

            if wrrversion == 1:
                if standard_name == "rainfall_flux": # Hack to support wrong names in WRR1 for now
                    standard_name = 'rainfal_flux'
            if wrrversion == 2 and (variable == "Rainfall" or variable == "SnowfallRate"):
                if variable == "TotalPrecipitation":
                    _tlist_, _timelist_ = get_times_P(currentdate, currentdate, serverroot, wrrsetroot, _filename[0],
                                                  timestepsecs, logger)
                    _tlist.append(_tlist_)
                    _timelist.append(_timelist_)
                    _tlist_, _timelist_ = get_times_P(currentdate, currentdate, serverroot, wrrsetroot,
                                                          _filename[1],timestepsecs, logger)
                    _tlist.append(_tlist_)
                    _timelist.append(_timelist_)
                else:
                    tlist, timelist = get_times_P(currentdate, currentdate, serverroot, wrrsetroot, filename,
                                            timestepsecs, logger)
            else:
                if variable == "TotalPrecipitation":
                    _tlist_, _timelist_ = get_times_P(currentdate, currentdate, serverroot, wrrsetroot, _filename[0],
                                                  timestepsecs, logger)
                    _tlist.append(_tlist_)
                    _timelist.append(_timelist_)
                    _tlist_, _timelist_ = get_times_P(currentdate, currentdate, serverroot, wrrsetroot,
                                                          _filename[1],timestepsecs, logger)
                    _tlist.append(_tlist_)
                    _timelist.append(_timelist_)
                else:
                    tlist, timelist = get_times(currentdate, currentdate, serverroot, wrrsetroot, filename,
                                                timestepsecs, logger)

            if len(_tlist) > 1: # Add Snow and rain
                ncstepobj = getstep(_tlist[0],BB,_standard_name[0],timestepsecs,logger)
                # get the steps for this time
                mstack = ncstepobj.getdates(_timelist[0])
                ncstepobj1 = getstep(_tlist[1], BB, _standard_name[1], timestepsecs, logger)
                mstack1 = ncstepobj1.getdates(_timelist[1])
                exec "thevar = mstack." + timeflatten + "(axis=0)"
                exec "thevar1 = mstack1." + timeflatten + "(axis=0)"
                thevar[thevar<valid_min]=NaN
                thevar[thevar > valid_max] = NaN
                thevar1[thevar1 < valid_min] = NaN
                thevar1[thevar1 > valid_max] = NaN
                thevar = thevar + thevar1
            else:
                ncstepobj = getstep(tlist,BB,standard_name,timestepsecs,logger)
                # get the steps for this time
                mstack = ncstepobj.getdates(timelist)
                exec "thevar = mstack." + timeflatten + "(axis=0)"
                thevar[thevar<valid_min]=NaN
                thevar[thevar > valid_max] = NaN


            #mapname = getmapname(cnt+1,oprefix)
            convstr = configget(logger, theconf, "conversion", variable, 'none')
            if convstr != 'none':
                convstr = convstr.replace(variable, 'thevar')
                varmetadata['conversionlog'] = convstr
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
                                  interpolmethod, ncoutfillval, hiresdem, resLowResDEM,interpolmethod,logger,localdatadir=odir)


                newdata[isnan(newdata)] = ncoutfillval
                newdata[~isfinite(newdata)] = ncoutfillval
                if netcdfout != 'None':
                    logger.info("Saving step to netcdf: " + str(thisstep))
                    if variable == 'TotalPrecipitation':
                        ostdname = 'precipitation_flux'
                    else:
                        ostdname = standard_name
                    ncout.savetimestep(cnt + 1, newdata, name=ostdname, var=variable, metadata=varmetadata)
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

                del ncstepobj

        cnt = cnt + 1




    logger.info("Done.")



if __name__ == "__main__":
    main()



