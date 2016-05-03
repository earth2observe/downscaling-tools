
This example determines reference potential ET for the rhine basin for 5  days

1) use the runrad.bat file (on windows) to run the radiation
correction example. The maps are saved in the output_rad dir. These maps are needed to run
the second step.

2) run the runevap.bat file to download the WDFEI forcing data from the earth2observe 
server and calculate reference ET (using the radiation maps to correct for the DEM) and
apply elevation based downscaling.