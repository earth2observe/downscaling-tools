@echo off
mkdir WorldClim

move prec_5m_bil.zip WorldClim\
move prec_30s_bil.zip WorldClim\
move prec_2-5m_bil.zip WorldClim\
move prec_10m_bil.zip WorldClim\

cd WorldClim
unzip prec_10m_bil.zip
mkdir prec_10m_bil
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec1.bil prec_10m_bil\prec0000.001
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec10.bil prec_10m_bil\prec0000.010
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec11.bil prec_10m_bil\prec0000.011
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec12.bil prec_10m_bil\prec0000.012
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec2.bil prec_10m_bil\prec0000.002
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec3.bil prec_10m_bil\prec0000.003
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec4.bil prec_10m_bil\prec0000.004
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec5.bil prec_10m_bil\prec0000.005
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec6.bil prec_10m_bil\prec0000.006
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec7.bil prec_10m_bil\prec0000.007
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec8.bil prec_10m_bil\prec0000.008
gdal_translate -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  prec9.bil prec_10m_bil\prec0000.009
del *.bil
del *.hdr
@echo off

set resolution=0.2500

mkdir ..\%resolution%

for %%f in (prec_10m_bil\prec????.???) do (
	echo %%f
	gdalwarp -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  -tr %resolution% %resolution% %%f   ..\%resolution%\%%~nf
)


cd ..