@echo off
rem normally run a whole year, for this example only 10 days!
..\..\e2o_dstools-64-bit\e2o_radiation.exe -S 1 -E 10 -D highResDEM\wflow_dem.map -d WRR1 -O output_rad -s 5 -e 22 -l DEBUG
