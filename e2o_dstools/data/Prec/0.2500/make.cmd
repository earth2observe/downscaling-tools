@echo off
for %%f in (..\0.0833\*.bil) do (
	echo %%f
	gdalwarp -tr 0.25 0.25 -r average -of GTiff -co COMPRESS=DEFLATE %%f %%~nf.tif
)
