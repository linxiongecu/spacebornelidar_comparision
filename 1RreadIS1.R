#### read is1 data
#### 06/06/2023
load('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/icesat_glas_umd_v1.RData')

glas_veg[1,]


#st_growing ed_growing
# if st < ed,   growing days [st, ed]
# if st > ed, growing days [0, ed] & [st, 365]

## filtering growing season data
## jday03 , glas launched in January 12, 2003

# Function to select column based on value
select_row <- function(row) {
  
  st <- row['st_growing']
  ed <- row['ed_growing']
  d  <- (row['jday03'] + 12) %% 365

  if (st < ed) {
    if (d > st && d < ed) {
      return (TRUE)
    }
    
  } else {
    if (d > st || d < ed) {
      return (TRUE)
    }
    
  }
  return (FALSE)
}
# Randomly select 1000 samples
# Randomly select 1000 rows
random_rows <- sample(nrow(glas_veg), size = 1000, replace = FALSE)
# Select the randomly sampled rows from the data frame
selected_df <- glas_veg[random_rows, ]
#apply(X, MARGIN, FUN)
# MARGIN: a matrix 1 indicates rows, 2 indicates columns
test <- apply(selected_df , 1, select_row)
selected_df['growing_flag'] <- test
# Clear row numbers
rownames(selected_df) <- NULL

### filter by flag value 
# Filter using logical indexing
filtered_df_grow<- selected_df[selected_df$growing_flag == TRUE, ]




# Apply the function to each row get a growing flag ; 
glas_veg['growing_flag']  <- apply(glas_veg, 1, select_row)
filtered_glas_veg<- glas_veg[glas_veg$growing_flag == TRUE, ]
is1_grow <- filtered_glas_veg[, c('lon', 'lat', 'ht')]






#### rasterize  data 
library(raster)
is1 <- glas_veg[, c('lon', 'lat', 'ht')]


write.csv(is1, file = "C:/Users/lxiong/Desktop/Data/IS1.csv", row.names = FALSE)


################################
is1 <- is1_grow
################################

#coordinates(is1) <- ~ lon + lat
#### map region 
# Set a new extent for the raster
#extent <- extent(xmin, xmax, ymin, ymax) 
e <- extent(-180, 180, -90, 90)
r <- raster(e, ncol=2880, nrow=1440, crs = '+proj=longlat +datum=WGS84') ### 0.125 degree
### 0.008 degree
r <- raster(e, ncol=45000, nrow= 22500, crs = '+proj=longlat +datum=WGS84')

# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(is1[, 1:2], r, is1[,3], fun=mean) # min, max, or mean


##### 98% height 
#x <- rasterize(is1[, 1:2], r, is1[,3], 
#               fun = function(i,...) quantile(i, probs=0.98, na.rm=T))

# Rasterize the points by count
is1_density <- rasterize(is1[, 1:2], r, is1[,3],fun = 'count')


# Merge the rasters
merged_raster <- stack(x, is1_density)
names(merged_raster) <- c('height', 'count')

##### save to rasterized IS1
# Save the raster to a file
#writeRaster(x, filename = "IS1_raster_14km_98percent.tif", format = "GTiff")
writeRaster(x, filename = "IS1_raster_14km_mean_growingSeason.tif", format = "GTiff")


### export to csv 
# Convert raster to XYZ format
xyz_data <- rasterToPoints(x)
# Write the data frame to a CSV file
write.csv(xyz_data, file = "is1_global.csv", row.names = FALSE)

library(rasterVis)
# Interpolate the raster
slope <- terrain(x, opt='slope')
aspect <- terrain(x, opt='aspect')
hill <- hillShade(slope, aspect, 40, 270)
plot(hill,  xlim = c(-180,180), ylim = c(-90, 90),
     col=grey(0:100/100), legend=FALSE, main='IS1 canopy')
plot(x, col=rainbow(25, alpha=0.35),add=TRUE)
#################
########









