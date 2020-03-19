# prepare_data_for_WRF-VPRM
this code help prepare data for WRF-VPRM simulation using MODIS datasets.

## input dataset needed
tile based EVI and LSWI from MODIS

tile based landcover map from MODIS

This code merge and reproject the MODIS tile data to 8km Lambert conformal conic projection. The merged and reprojected data is then combined and generate a netcdf file.

GDAL library is needed.



