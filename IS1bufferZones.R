########################Create a buffer =500 m and save in gpkg format 
library(sf)
library(sp)
load('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/icesat_glas_umd_v1.RData')
is1 <- glas_veg[, c('lon', 'lat', 'ht')]

########### how many 0.125 degree grids ????
# grid IS1 first, export to gridded xyz.
#############



for (i in 1:2880-1) {
   lon = -180 + i*0.125
   #print(lon)
   for (j in 1:832-1) { 
   lat = -52 + j*0.125
   
   
   #lon = -54.54
   #lat = -6.16
   
   ### for each grid 
   # Step 1: filter is1 points 
   ### 
   # Define the center point and resolution
   # center <- c(13.3125,49.1875)      # Center coordinates (longitude, latitude)
   resolution <- 0.125     # Resolution in degrees
   
   # Calculate the half-length of the square in each dimension
   half_length <- resolution / 2
   
   # Calculate the coordinates of the square vertices
   x_coord <- lon + c(-1, 1, 1, -1, -1) * half_length
   y_coord <- lat + c(-1, -1, 1, 1, -1) * half_length
   # Filter points by grid polygon 
   in_flag <- point.in.polygon(is1$lon, is1$lat, x_coord, y_coord)
   filtered_points <- data.frame(cbind(is1$lon, is1$lat, 
                                       is1$ht, in_flag))
   filtered_points <- filtered_points[filtered_points$in_flag > 0, ]
   if (nrow(filtered_points) == 0) {next}
   print(paste(lon, lat))
   ### Step 2 buffer points 
   dat_sf <- st_as_sf(filtered_points, coords = c("V1", "V2"), crs = 4326)
   dat_circles <- st_buffer(dat_sf, dist = 500) # too long 
   dat_circles_union <- st_union(dat_circles)
   #plot(dat_circles_union)
   bufferName = paste0("Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Out/is1_buffer500m_", i,'_', j, '.gpkg')
   print(bufferName)
   st_write(dat_circles_union, bufferName, layer = "icesat", delete_layer = TRUE) # overwrites)
   }
}


# Save sf object to GeoPackage file
#st_write(buffered_points, "Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/is1_buffer500m.gpkg", layer = "icesat")

