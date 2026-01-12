# Load library
library(gstat)
library (sf)
library (sp)
library (ggplot2)
library(automap)
library(sfheaders)
library(terra)
library(dplyr)

#Define Filepath
dir=choose.dir()
csv.path=paste0(dir,"\\atl08p_2_2025-06-03_13_56_19_889.parquet.csv")
egm.path=paste0(dir,"\\us_nga_egm2008_1.tif")
bound.path=paste0(dir,"\\polygon_138.geojson")
output.dir=paste0(dir,"\\")
transect.dir=paste0(dir,"\\Transect.geojson")
demnas.dir=paste0(dir,"\\DEMNAS_Clipped.tif")
als.dir=paste0(dir,"\\Polygon_138_utm_49S_dtm.idw.tif")


#Load dataset
df=read.csv(csv.path)
df$x=1:nrow(df)
egm=rast(egm.path)
bound_sf<-geojsonsf::geojson_sf(bound.path)
bound_sf<-st_transform(bound_sf,4326)
egm.p=terra::crop(egm, bound_sf, mask=TRUE)
egm.p=terra::project(egm.p,"EPSG:32749")
bound_sf<-st_transform(bound_sf,32749)
bound_sf<-st_buffer(bound_sf,200)
tr.sf=st_read(transect.dir)
demnas=rast(demnas.dir)
als=rast(als.dir)

#Convert tabular to spatial
sf=sf_multipoint(
  df,
  x = "longitude",
  y = "latitude",
  multipoint_id = "x",
  keep = TRUE
)

st_crs(sf)<-4326
sf<-st_transform(sf,32749)

#Filter Data
Q1 <- quantile(df$h_te_median,.25)
Q3 <- quantile(df$h_te_median, .75)
IQR <- IQR(df$h_te_median)
filter_sf <- subset(sf, sf$h_te_median> (Q1 - 2*IQR) & sf$h_te_median< (Q3 + 2*IQR))

#Add Orthometric Height
egm_sf=extract(egm.p,vect(filter_sf))
colnames(egm_sf)[1]="x"
join=filter_sf %>% left_join(egm_sf,by="x") 
join$ortho<-join$h_te_median-join$geoid_undulation
join<-join%>% na.omit("ortho")

#Kriging Interpolation
sf_grid <- bound_sf %>% 
  st_make_grid(cellsize = 200, what = "centers") %>% # grid of points
  st_intersection(bound_sf)

print(paste("Process Started: ",Sys.time()))
time_start=Sys.time()

kriging_result = autoKrige(ortho~1, join, sf_grid)

time_finish=Sys.time()
time_take=time_finish-time_start

print(paste("Process Finished: ",Sys.time()))
print(paste("Time Taken: ",time_take))

res<-st_as_sf(kriging_result$krige_output)

#Convert Grid points to raster
krg_raster <- rast(data.frame(
  x = st_coordinates(res)[,1],
  y = st_coordinates(res)[,2],
  pred_value = res$var1.pred
), type = "xyz", crs = "EPSG: 32749")

template_resamp <- rast(extent=ext(krg_raster), crs = crs(krg_raster), resolution = 1)
krg_raster.1  <- resample(krg_raster, template_resamp, method = "bilinear")

#Visualize Result
plot(krg_raster.1)
title("DTM ICESat-2 - Polygon 138")

#Extract raster based in transect check points
smp=st_line_sample(tr.sf,density = 1/100)
vect=vect(smp)
sample.ics <- terra::extract(krg_raster.1, vect)
sample.als <- terra::extract(als, vect)
sample.dns <- terra::extract(demnas, vect)
sample.df=cbind(sample.ics,sample.als,sample.dns)
colnames(sample.df)[2]="ics"
colnames(sample.df)[4]="als"
colnames(sample.df)[6]="dns"
colnames(sample.df)[3]="ID_2"
colnames(sample.df)[5]="ID_3"

#Calculate Evaluation Metrics
rmse.ics=sqrt(mean((sample.df$als - sample.df$ics )^2))
rmse.dns=sqrt(mean((sample.df$als  - sample.df$dns )^2))

#Export Result
st_write(res,paste0(output.dir,"krig_res.shp"),append=FALSE)
writeRaster(krg_raster.1,paste0(output.dir,"kriging_prediction_1m.tif"),overwrite=TRUE)
st_write(sample.df,paste0(output.dir,"sample_df.shp"),append=FALSE)
st_write(smp,paste0(output.dir,"sample_pts.shp"),append=FALSE)

#Plot Graph
plot.ics <- ggplot(sample.df,aes(x=als, y=ics)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=5, y=8, label=paste0("RMSE: ",round(rmse.ics,2)," m"),
           color="black")+
  labs(title ="ICESat DTM",x = "Actual (m)", y = "Predicted (m)", hjust=0.5) +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,10) +
  ylim(0,10)

plot.dns <- ggplot(sample.df,aes(x=als, y=dns)) +
  geom_point(shape=21,colour = "blue", fill = "white",size=1,stroke = 1.5) + 
  geom_abline(a=0,b=0,linetype=1,colour="red")+
  annotate(geom="text", x=5, y=8, label=paste0("RMSE: ",round(rmse.dns,2)," m"),
           color="black")+
  labs(title ="DEMNAS",x = "Actual (m)", y = "Predicted (m)") +
  coord_fixed()+  # Set aspect ratio t  o be equal
  xlim(0,10) +
  ylim(0,10)

grid.arrange(plot.ics+ theme(plot.title = element_text(hjust = 0.5)),
             plot.dns+ theme(plot.title = element_text(hjust = 0.5)),
             ncol=2)
