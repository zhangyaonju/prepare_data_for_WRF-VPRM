
######first I need to run get the select tiles and merged for this year
#for EVI first
year<-2018
dir.create(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/EVI"),recursive = T)
dir.create(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/LSWI"),recursive = T)
dir.create(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/reproj_EVI"),recursive = T)
dir.create(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/reproj_LSWI"),recursive = T)

for (i in (1:46)*8-7){
  system(paste("gdal_merge.py -o /data/ifs/VPM/request/CA_VPRM",year,"/EVI/EVI.",year,".",
               formatC(i,flag="0",width = 3),".tif /data/ifs/VPM/driving_data/MOD09A1_006_BISE_SG/",year,"/",
               "{h21v03,h22v03,h23v03,h24v03,h25v03,h26v03,h22v04,h23v04,h24v04,h25v04,h26v04,h27v04,",
               "h23v05,h24v05,h25v05,h26v05,h27v05,h28v05,h23v06,h24v06,h25v06,h26v06,h27v06,h28v06,",
               "h29v06,h24v07,h25v07,h26v07,h27v07,h28v07,h29v07}/MOD09A1.A",year,"",
               formatC(i,flag="0",width = 3),"*.tif",sep=""))
}

system(paste("for i in /data/ifs/VPM/request/CA_VPRM",year,"/EVI/*.tif; do gdalwarp -t_srs ",
             "'+proj=lcc +lat_1=36.138999938964844 +lat_2=36.138999938964844 +lat_0=36.13899230957031 ",
             "+lon_0=100.1500015258789 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs' ",
             "-tr 20000 20000 -te -2850000 -2340000 2850000 2340000 -r average -srcnodata -32768 $i ",
             "/data/ifs/VPM/request/CA_VPRM",year,"/reproj_EVI/${i##*/}; done",sep=""))


#for LSWI then
for (i in (1:46)*8-7){
  system(paste("gdal_merge.py -o /data/ifs/VPM/request/CA_VPRM",year,"/LSWI/LSWI.",year,".",
               formatC(i,flag="0",width = 3),".tif /data/ifs/modis/products_006/mod09a1/geotiff/lswi/",year,"/",
               "{h21v03,h22v03,h23v03,h24v03,h25v03,h26v03,h22v04,h23v04,h24v04,h25v04,h26v04,h27v04,",
               "h23v05,h24v05,h25v05,h26v05,h27v05,h28v05,h23v06,h24v06,h25v06,h26v06,h27v06,h28v06,",
               "h29v06,h24v07,h25v07,h26v07,h27v07,h28v07,h29v07}/MOD09A1.A",year,"",
               formatC(i,flag="0",width = 3),"*.tif",sep=""))
}

system(paste("for i in /data/ifs/VPM/request/CA_VPRM",year,"/LSWI/*.tif; do gdalwarp -t_srs ",
             "'+proj=lcc +lat_1=36.138999938964844 +lat_2=36.138999938964844 +lat_0=36.13899230957031 ",
             "+lon_0=100.1500015258789 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs' ",
             "-tr 20000 20000 -te -2850000 -2340000 2850000 2340000 -r average -srcnodata -32768 $i ",
             "/data/ifs/VPM/request/CA_VPRM",year,"/reproj_LSWI/${i##*/}; done",sep=""))


system(paste("gdal_merge.py -o /data/ifs/VPM/request/CA_VPRM",year,"/LSWImax.",year,".tif",
             " /data/ifs/VPM/driving_data/LSWImax_MA_06/",year-2,"/",
             "{h21v03,h22v03,h23v03,h24v03,h25v03,h26v03,h22v04,h23v04,h24v04,h25v04,h26v04,h27v04,",
             "h23v05,h24v05,h25v05,h26v05,h27v05,h28v05,h23v06,h24v06,h25v06,h26v06,h27v06,h28v06,",
             "h29v06,h24v07,h25v07,h26v07,h27v07,h28v07,h29v07}.",year-2,".maxLSWI_MA.tif",
             sep=""))


system(paste("gdalwarp -t_srs ",
             "'+proj=lcc +lat_1=36.138999938964844 +lat_2=36.138999938964844 +lat_0=36.13899230957031 ",
             "+lon_0=100.1500015258789 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs' ",
             "-tr 20000 20000 -te -2850000 -2340000 2850000 2340000 -r average -srcnodata -32768 ",
             "/data/ifs/VPM/request/CA_VPRM",year,"/LSWImax.",year,".tif ",
             "/data/ifs/VPM/request/CA_VPRM",year,"/LSWImax_reproj.",year,".tif -overwrite",sep=""))

library(raster)
library(rgdal)
library(ncdf4)
year=2018
setwd(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/"))
#setwd("V:/VPM/request/CA_VPRM/")
EVI_file<-list.files("./reproj_EVI/",pattern="*.tif",full.names = T)
EVIstack<-stack(EVI_file)
LSWI_file<-list.files("./reproj_LSWI/",pattern="*.tif",full.names = T)
LSWIstack<-stack(LSWI_file)

minLSWI<-calc(LSWIstack,fun = function(x) {quantile(x,probs = 0.05,na.rm=TRUE)})
writeRaster(minLSWI,paste0("./LSWImin_reproj.",year,".tif"),overwrite=T)

minEVI<-calc(EVIstack,min,na.rm=T)
writeRaster(minEVI,paste0("./EVImin_reproj.",year,".tif"))

maxEVI<-calc(EVIstack,max,na.rm=T)
writeRaster(maxEVI,paste0("./EVImax_reproj.",year,".tif"))


lc_files<-list.files("/data/ifs/VPM/request/CA_VPRM/LC_new/",full.names=T)
if (length(lc_files)<8){
  
  system(paste("gdalwarp -t_srs ",
               "'+proj=lcc +lat_1=36.138999938964844 +lat_2=36.138999938964844 +lat_0=36.13899230957031 ",
               "+lon_0=100.1500015258789 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs' ",
               "-tr 500 500 -te -2850000 -2340000 2850000 2340000 -r near -srcnodata -32768 -overwrite ",
               "/data/ifs/VPM/request/CA_VPRM/MCD12Q1.500.2016.tif ",
               "/data/ifs/VPM/request/CA_VPRM/MCD12Q1.reproj.500.2016.tif",sep=""))
  
  lc_data<-raster('/data/ifs/VPM/request/CA_VPRM/MCD12Q1.reproj.500.2016.tif')
  EVI_data<-raster("/data/ifs/VPM/request/CA_VPRM/EVImax_reproj.2016.tif")
  
  lc<-c('ENF','EBF','DNF','DBF','MF','CSH','OSH','WSA','SAV','GRA','WET','CRO','URB','NVM','SNO','BAR','WAT')
  si<-40
  
  percent<-array(NA,dim=c(234,285,16))
  
  for (i in 1:234){
    for (j in 1:285){
      data<-lc_data[(1+floor((i-1)*si)):floor(i*si),(1+floor((j-1)*si)):floor(j*si)]
      n<-length(data)
      for (k in 1:16){
        percent[i,j,k]<-sum(data==k)/n
      }
      #percent[i,j,17]<-sum(data==0)/n
    }
  }
  for (i in 1:16){
    out<-setValues(EVI_data,as.vector(t(percent[,,i])))
    writeRaster(out,paste('/data/ifs/VPM/request/CA_VPRM/LC/',formatC(i,format='d',width=2,flag='0'),
                          lc[i],'landcover_percent.2016.tif'),overwrite=T)
  }
  
  d<-c(1,1,2,2,3,4,4,5,5,7,7,6,8,7,8,7)
  lc8<-c("EF","DF","MF","SH","SA","CR","GR","OT")
  LC_f<-list.files("/data/ifs/VPM/request/CA_VPRM/LC/",pattern="*.tif",full.names = T)
  for (i in 1:8){
    biome_f<-LC_f[d==i]
    biome_stack<-stack(biome_f)
    biome_pct<-sum(biome_stack)
    writeRaster(biome_pct,paste("/data/ifs/VPM/request/CA_VPRM/LC_new/",i,"_",lc8[i],".tif",sep=''),overwrite=T)
  }
}


#create and write the netCDF file
library(raster)
#library(rgdal)
library(ncdf4)

setwd(paste0("/data/ifs/VPM/request/CA_VPRM",year,"/"))

x3<-seq(-2840,2840,20)
y3<-seq(-2330,2330,20)
t3<-seq(0,360,8)
xdim<-ncdim_def("x","km_east",as.double(x3))
ydim<-ncdim_def("y","km_north",as.double(y3))
timedim<-ncdim_def("time","days since 2016/01/01",as.double(t3))
lcdim<-ncdim_def("8 land cover type","percent of land cover for EF DF MF SH SA CR GR OT",as.double(1:8))
###
#define variables
fillvalue<-1e-9
dlname<-"EVI minimum"
EVImin_def<-ncvar_def("EVImin",'NA',list(xdim,ydim),fillvalue,dlname,prec = "single")
dlname<-"EVI maximum"
EVImax_def<-ncvar_def("EVImax",'NA',list(xdim,ydim),fillvalue,dlname,prec = "single")
dlname<-"EVI"
EVI_def<-ncvar_def("EVI",'NA',list(xdim,ydim,timedim),fillvalue,dlname,prec = "single")
dlname<-"LSWI minimum"
LSWImin_def<-ncvar_def("LSWImin",'NA',list(xdim,ydim),fillvalue,dlname,prec = "single")
dlname<-"LSWI maximum"
LSWImax_def<-ncvar_def("LSWImax",'NA',list(xdim,ydim),fillvalue,dlname,prec = "single")
dlname<-"LSWI"
LSWI_def<-ncvar_def("LSWI",'NA',list(xdim,ydim,timedim),fillvalue,dlname,prec = "single")
################
dlname<-"landcover"
Landcover_def<-ncvar_def("landcover",'NA',list(xdim,ydim,lcdim),fillvalue,dlname,prec = "single")

# create netCDF file and put arrays
ncfname <- paste0("./all_data_",year,".nc")
ncout <- nc_create(ncfname,list(EVImin_def,EVImax_def,EVI_def,LSWImin_def,LSWImax_def,LSWI_def,Landcover_def))

#put variables
EVI_F<-list.files("./reproj_EVI/",pattern='*.tif',full.names = T)
EVI<-stack(EVI_F)
evi_array<-getValues(EVI)/10000
evi_array[evi_array==-32768]=-2000
dim(evi_array)<-c(285,234,46)
ncvar_put(ncout,EVI_def,evi_array)

LSWI_F<-list.files("./reproj_LSWI/",pattern='*.tif',full.names = T)
LSWI<-stack(LSWI_F)
lswi_array<-getValues(LSWI)/10000
dim(lswi_array)<-c(285,234,46)
ncvar_put(ncout,LSWI_def,lswi_array)

LC_F<-list.files("/data/ifs/VPM/request/CA_VPRM/LC_new/",pattern='*.tif',full.names = T)
LC<-stack(LC_F)
lc_array<-getValues(LC)
dim(lc_array)<-c(285,234,8)
ncvar_put(ncout,Landcover_def,lc_array)


EVImin<-raster(paste0("./EVImin_reproj.",year,".tif"))
evimin_array<-getValues(EVImin)/10000
dim(evimin_array)<-c(285,234)
ncvar_put(ncout,EVImin_def,evimin_array)

EVImax<-raster(paste0("./EVImax_reproj.",year,".tif"))
evimax_array<-getValues(EVImax)/10000
dim(evimax_array)<-c(285,234)
ncvar_put(ncout,EVImax_def,evimax_array)

LSWImin<-raster(paste0("./LSWImin_reproj.",year,".tif"))
LSWImin_array<-getValues(LSWImin)/10000
dim(LSWImin_array)<-c(285,234)
ncvar_put(ncout,LSWImin_def,LSWImin_array)

LSWImax<-raster(paste0("./LSWImax_reproj.",year,".tif"))
LSWImax_array<-getValues(LSWImax)/10000
dim(LSWImax_array)<-c(285,234)
ncvar_put(ncout,LSWImax_def,LSWImax_array)

#put additional atrributes into dimension and data variables
ncatt_put(ncout,"x", 'axis', "X")
ncatt_put(ncout,"y", 'axis', "Y")
ncatt_put(ncout,"time", 'axis', "T")

ncatt_put(ncout,0,"title","data for VRPM simulation")
nc_close(ncout)




