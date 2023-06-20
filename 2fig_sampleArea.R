### Fig for sampling area
library(raster)
library(sf)

############each 0.125 degree grid, what is the sampling area ????
############ common grid have samples;

# ###### Step 1 get command grids
# gedi_gridded <- read.csv("gedi_calval_shot_filtered_gridded_optimal_filterDisturbance.csv")
# # Extract raster values at the specified coordinates
# x <- gedi_gridded$x
# y <- gedi_gridded$y
# is1 <- raster('IS1_raster_14km.tif')
# is1_values <- extract(is1 , cbind(x, y))
# gedi_rh98 <-   gedi_gridded$layer
# result <- cbind(x, y,  gedi_rh98, is1_values)
# #result <- cbind(x, y,  gedi_rh98, is1_values_98percent)
# result <- data.frame(result)
# # Remove rows with NA values
# result <- result[complete.cases(result), ]
##### 722 command grids 
cmgrids <- read.csv("Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Out/result_gridded_06192023.csv")

####for each grid , get sample area 
### Step 1 filter data by the cell
gedi_shots_filtered_final <- read.csv('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Out/gedi_calval_shot_filtered_06192023.csv')
#gedi_shots_filtered_final<- as.data.frame(lapply(gedi_shots_filtered_final, as.numeric))
load('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/icesat_glas_umd_v1.RData')
is1 <- glas_veg[, c('lon', 'lat', 'ht')]



library(sf)

############Step 2 in each grid, only keep shots of GEDI/IS1 are within 500 m.
# Create an empty data frame
samples_df <- data.frame()
#### LOOP THROUGH rows
#### 
# Get the number of rows in the data frame
num_rows <- nrow(cmgrids)
#num_rows <- 10
# Loop through each row
for (i in 1:num_rows) {
  # Access the row using 'i'
  current_row <- cmgrids[i, ]
  print(paste('Processing Grid Number: ', i))
  
  # Perform operations on the current row
  # ...
  
  
  # Define the center point and resolution
  # center <- c(13.3125,49.1875)      # Center coordinates (longitude, latitude)
  resolution <- 0.125     # Resolution in degrees
  
  # Calculate the half-length of the square in each dimension
  half_length <- resolution / 2
  
  # Calculate the coordinates of the square vertices
  x_coord <- current_row$x + c(-1, 1, 1, -1, -1) * half_length
  y_coord <- current_row$y + c(-1, -1, 1, 1, -1) * half_length
  #xym <- data.frame(cbind(x_coord, y_coord))
  # Create an sf object
  # library(sp)
  # p = Polygon(xym)
  # ps = Polygons(list(p),1)
  # sps = SpatialPolygons(list(ps))
  # plot(sps)
  
  # Filter points by grid polygon 
  in_flag <- point.in.polygon(gedi_shots_filtered_final$lon, gedi_shots_filtered_final$lat, x_coord, y_coord)
  filtered_points <- data.frame(cbind(gedi_shots_filtered_final$lon, gedi_shots_filtered_final$lat, 
                                      gedi_shots_filtered_final$rh98, in_flag))
  filtered_points <- filtered_points[filtered_points$in_flag > 0, ]
  
  if (nrow(filtered_points) == 0) {next}
  
  # For IS1 data
  in_flag_is1 <- point.in.polygon(is1$lon, is1$lat, x_coord, y_coord)
  filtered_points_is1 <- data.frame(cbind(is1$lon, is1$lat, is1$ht, in_flag_is1))
  filtered_points_is1 <- filtered_points_is1[filtered_points_is1$in_flag_is1 > 0, ]
  if (nrow(filtered_points_is1) == 0) {next}

  

  
  #### another method  
  #st_is_within_distance
  gedi_sf <- st_as_sf(filtered_points, coords = c("V1", "V2"), crs = 4326)
  is1_sf <- st_as_sf(filtered_points_is1, coords = c("V1", "V2"), crs = 4326)
  
  # Check if points are within the distance threshold from the reference point
  is_within_distance <- st_is_within_distance(gedi_sf, is1_sf, dist = 500) #500 m
  
  ####  Here we only filter GEDI, can I filter IS1?
  ####  I want both shots <= 1km?
  
  filtered_points$flag <- apply(is_within_distance , 1, any)
  filtered_points <- filtered_points[filtered_points$flag == TRUE, ]
  
  ##update GEDI sf
  gedi_sf <- st_as_sf(filtered_points, coords = c("V1", "V2"), crs = 4326)
  is_within_distance <- st_is_within_distance(is1_sf, gedi_sf, dist = 500) #500 m
  filtered_points_is1$flag <- apply(is_within_distance , 1, any)
  filtered_points_is1 <- filtered_points_is1[filtered_points_is1$flag == TRUE, ]

  
  if (nrow(filtered_points) == 0) {next}
  if (nrow(filtered_points_is1) == 0) {next}
  
  plot(filtered_points$V1,filtered_points$V2 )
   ###Add the second point cloud to the plot
  points(filtered_points_is1$V1,filtered_points_is1$V2, col = "red", pch = 16)

  
  diff = mean(filtered_points$V3) - mean(filtered_points_is1$V3)
  
  samples_df <- rbind(samples_df, data.frame(Grid = i,
                                             N_GEDI = nrow(filtered_points), 
                                             N_IS1 = nrow(filtered_points_is1),
                                             GEI_h = mean(filtered_points$V3),
                                             IS1_h = mean(filtered_points_is1$V3),
                                             diff = diff
                                             ))
  
  
  
  
  ### Step 2 each point get a circle 
  # 
  # # Buffer circles by 100m
  # dat_sf <- st_as_sf(filtered_points, coords = c("V1", "V2"), crs = 4326)
  # dat_circles <- st_buffer(dat_sf, dist = 12.5) # too long 
  # 
  # # Find overlaps and union features
  # parts <- st_cast(st_union(dat_circles),"POLYGON")
  # plot(parts)
  # 
  # 
  # #plot(dat_sf)
  # #plot(dat_circles)
  # # Calculate the area of the polygon
  # #area <- st_area(dat_circles)
  # 
  # 
  # ### Step 3 polygons union 
  # # Plot the sf object using ggplot2
  # ggplot() +
  # geom_sf(data = data_com, fill = "red")+ 
  # geom_sf(data = dat_sf, color = "blue")+
  # xlim(13.31, 13.32)+ 
  # ylim(49.06, 49.07)
  #   
  # 
  # #plot(st_union(dat_circles))
  # area <- st_area(totalsurface)
  # 
  # 
  # cell_sf <- st_polygon(
  #   list(
  #     cbind(
  #       x_coord, 
  #       y_coord)
  #   )
  # )
  # 
  # cell_sf = st_sfc(cell_sf, crs=4326)
  # plot(cell_sf )
  # area_total <- st_area(cell_sf)
  # 
  # ### from GEE
  # ## Rectangle area: m^2
  # ### 126268001.87933438
  # 
  # ### Step 4 sample area percentage 
  # Spercent <- area  / area_total * 100
  # # Print or store the results
  # # ...
  # print(Spercent)
}


#################fig 

samples_df_order <- samples_df[order(samples_df$N_GEDI), ]

# Example data frames
df1 <- data.frame(
  x = seq(1,nrow(samples_df_order)),
  y = samples_df_order$N_GEDI,
  group = "GEDI"
)

df2 <- data.frame(
  x = seq(1,nrow(samples_df_order)),
  y = samples_df_order$N_IS1,
  group = "ICESat"
)

# Combine the data frames
combined_df <- rbind(df1, df2)



# Plot the values using ggplot
ggplot(combined_df , aes(x , y =y , fill = group)) +
  geom_bar(stat = "identity", position = "identity") +
  xlab("Grid index") +
  ylab("Number of GEDI/IS1 shots per grid") + 
  coord_cartesian(clip = "off", ylim= c(0,1000))
  


 
samples_df$diff_number  <-  samples_df$N_GEDI - samples_df$N_IS1

ggplot(samples_df , aes(Grid , diff)) +
  geom_bar(stat = "identity") +
  xlab("Grids number") +
  ylab("GEDI - IS1 diff") +
  ylim(-10, 10)

## his
ggplot(samples_df , aes(diff )) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.7) +
  xlab("Difference") +
  ylab("GEDI/IS1 number")


#### median plot 

#### 
result <- cbind(samples_df$IS1_h,samples_df$GEI_h)
result <- data.frame(result)  
colnames(result) <- c('is1_values', 'gedi_rh98')
# Define height class intervals
height_bins <- cut(result$is1_values, breaks = seq(0.5,60.5))

# Group data by height class and calculate mean weight
gedi_median <- aggregate(result$gedi_rh98, by = list(height_class = height_bins), FUN = median)
gedi_mean <- aggregate(result$gedi_rh98, by = list(height_class = height_bins), FUN = mean)
gedi_median_x <- as.character(gedi_median$height_class)


# Load the stringr package
library(stringr)
gedi_median_x <-  str_extract(gedi_median_x, "(?<=\\,).*?(?=\\.)")
gedi_median_x <- as.numeric(gedi_median_x)
# Calculate the number of samples in each class
height_bins <- cut(result$is1_values, breaks = seq(0.5,60.5,2))

result$height_class <- height_bins 

sample_counts <- result %>% count(height_class)
x_sample <-  str_extract(sample_counts$height_class, "(?<=\\,).*?(?=\\.)")

sample_counts$x <- as.numeric(x_sample)



df_1 <- data.frame(gedi_median_x, gedi_median$x ) # median
colnames(df_1) <- c('x', 'y')

df_1<- merge(df_1, sample_counts, by = "x", all = TRUE)
df_1['type'] = 'median'
df_2 <- data.frame(gedi_median_x, gedi_mean$x) # mean
colnames(df_2) <- c('x', 'y')


df_2<- merge(df_2, sample_counts, by = "x", all = TRUE)
df_2['type'] = 'mean'
# Combine the data frames


combined_df <- rbind(df_1, df_2)
# Remove NA rows
combined_df <- na.omit(combined_df)

df_1 <- na.omit(df_1 )

p2 <- ggplot(combined_df,aes(x , y , color = type, fill = type)) +
  geom_point(aes(shape=type),size =3) +
  scale_color_manual(values = c("mean" = "red", "median" = "blue") )+
  #scale_fill_manual(values = c("mean" = "orange", "median" = "lightblue") )+
  scale_shape_manual(values = c("mean" = 16, "median" = 17))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        #legend.position = c(0.1, 0.8),
  )+
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed",linewidth = 1 )+
  xlim(c(0,30)) +
  ylim(c(0,30))+
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
  # geom_text( aes(
  #                                 label = format(n, big.mark = ",", 
  #                                 scientific = FALSE)), vjust = -1, size = 5)+
  labs(x = "Icesat-1 [m]", y = "GEDI_cal/val_shots [m]", title = "Median of GEDI at 1/8° from ICESat")

p2


# Save the plot as a TIFF file
tiff("plot_median_within500m.tiff", width = 6, height = 5, units = "in", res = 300)
p2
dev.off()





