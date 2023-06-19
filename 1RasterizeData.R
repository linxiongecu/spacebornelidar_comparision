### Script to raserize globally 
is1 <- glas_veg[, c('lon', 'lat', 'ht')]

#coordinates(is1) <- ~ lon + lat
#### map region 
# Set a new extent for the raster
#extent <- extent(xmin, xmax, ymin, ymax) 

#r <- raster(e, ncol=2880, nrow=1440, crs = '+proj=longlat +datum=WGS84') ### 0.125 degree 1/8
e <- extent(-180, 180, -90, 90)
r <- raster(e, ncol=11520, nrow= 5760, crs = '+proj=longlat +datum=WGS84') ### 1/32 degree ---3km

r <- raster(e, ncol=23040, nrow= 11520, crs = '+proj=longlat +datum=WGS84') ### 1/64 degree ---1.7km
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(is1[, 1:2], r, is1[,3], fun=mean) # min, max, or mean

##### save to rasterized IS1
# Save the raster to a file
#writeRaster(x, filename = "IS1_raster_14km_98percent.tif", format = "GTiff")
writeRaster(x, filename = "IS1_raster_1dot7km.tif", format = "GTiff")