@echo off
rem normally run a whole year, for this example only 10 days!
..\..\e2o_downscale-2017-1-normal-win32-64\e2o_radiation.exe -S 1 -E 10 -D highResDEM\wflow_dem.map -d WRR2 -O output_rad -s 5 -e 22 -l DEBUG
