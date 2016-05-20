@echo off
echo Get precipitation and Temperature for the first 15 days of the first month of 2010 and 
echo resamples the data to a local dem using linear interpolation.
pause Press enter to tun...
..\..\e2o_dstools-64-bit\e2o_getvar.exe -I examplerun1.ini
echo results can be found in the output directory

