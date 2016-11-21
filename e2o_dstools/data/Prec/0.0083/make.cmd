@echo off
for %%f in (*.bil) do (
	echo %%f
	gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9 %%f %%~nf.tif
)
