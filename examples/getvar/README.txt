The directory contain three example config files. If the tools are installed 
correctly you can run them using the following commands:

1: Resample MSWEP precipitation (Rainfall only!) for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun1.ini

2: Resample and downscale MSWEP precipitation using nearest interpolation (Rainfall only!) 
for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun2.ini

3:Resample and downscale MSWEP precipitation using linear interpolation (Rainfall only!) 
for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun3.ini

4:Resample and downscale MSWEP precipitation using linear interpolation (Rainfall only!) 
for the Rhine basin to 0.0366666666 degree and save into a netcdf file
e2o_getvar.exe -I examplerun4.ini

5:Resample and downscale MSWEP precipitation using linear interpolation (Snowfall only!) 
for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun5.ini

6:Resample and downscale WFDEI precipitation using linear interpolation (Rainfall only!) 
for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun6.ini

7:Resample and downscale MSWEP precipitation+snowfall (TotalPrecipitation) using linear interpolation 
for the Rhine basin to 0.0366666666 degree
e2o_getvar.exe -I examplerun7.ini