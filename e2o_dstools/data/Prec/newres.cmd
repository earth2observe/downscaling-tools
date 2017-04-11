
@echo off
set resolution=0.5000

mkdir %resolution%

for %%f in (0.0083/prec*.*) do (
	gdalwarp -r average -of GTiff -co COMPRESS=DEFLATE -co ZLEVEL=9  -tr %resolution% %resolution% 0.0083\%%f   %resolution%\%%f
)
