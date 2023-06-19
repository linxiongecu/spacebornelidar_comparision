#### Compare GEDI calval vs IS1

### read gridded 0.125 degree IS1 data

# Load required packages
library(raster)
library(ggplot2)
library(dplyr)
# Load required packages
library(arrow)
# 
# gedi_fsbd <- read.csv('GEDIcalvalMetrics.csv')
# 
# ############# grid gedi data to 0.125 degree grid first 
# 
# e <- extent(-180, 180, -90, 90)
# r <- raster(e, ncol=2880, nrow=1440, crs = '+proj=longlat +datum=WGS84') ### 0.125 degree
# # you need to provide a function 'fun' for when there are multiple points per cell
# x <- rasterize(gedi_fsbd[, 7:8], r, gedi_fsbd[,4], fun=mean) # min, max, or mean
# 
# ### export to csv 
# # Convert raster to XYZ format
# xyz_data <- rasterToPoints(x)
# # Write the data frame to a CSV file
# write.csv(xyz_data, file = "gedi_fsbd_gridded.csv", row.names = FALSE)
# 
# ##########
# gedi_gridded <- read.csv("gedi_fsbd_gridded.csv")
# 
# 
# #### given x and y ----give me is1 raster value 
# # Extract raster values at the specified coordinates
# x <- gedi_gridded$x
# y <- gedi_gridded$y
# is1_values <- extract(is1 , cbind(x, y))
# 
# 
# #### result 
# 
# gedi_rh98 <- gedi_gridded$layer
# result <- cbind(x, y,  gedi_rh98, is1_values)
# result <- data.frame(result)
# 
# # Remove rows with NA values
# result <- result[complete.cases(result), ]
# 
# 
# 
# #### statistics 
# #### 
# N = nrow(result)
# Bias = mean(result$gedi_rh98-result$is1_values)
# MAE = mean(abs(result$gedi_rh98-result$is1_values))
# r <- cor(result$is1_values, result$gedi_rh98)
# R2 <- r^2
# # Calculate RMSE
# rmse <- sqrt(mean((result$is1_values - result$gedi_rh98)^2))
# rmse_percent <- rmse/mean(result$gedi_rh98)*100
# 
# ##### add text in fig 
# # Create example data
# lx <- c(40, 40, 40, 50, 40)
# ly <- c(20, 15, 10, 10, 5)
# labels <- c(paste('N =',round(N)), paste('Bias =',round(Bias,1), 'MAE =', round(MAE, 1)), 
#             paste('r =',round(r, 2)), paste('R^2 =', round(R2,2)),
#             paste('rmse =', round(rmse, 1), '%rmse =', round(rmse_percent, 1))  )
# 
# 
# # Create a scatterplot with density
# #geom_text() adds only text to the plot. 
# #geom_label() draws a rectangle behind the text, making it easier to read.
# 
# 
# p <- ggplot(result, aes(x = result$is1_values, y = result$gedi_rh98)) +
#   geom_bin2d(binwidth = c(1.5, 1.5)) +
#   scale_fill_continuous(type = "viridis", limits=c(0,10),oob = scales::squish) +
#   #geom_point() +
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank())+
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed",linewidth = 2 )+
#   xlim(c(0,60)) +
#   ylim(c(0,60))+
#   #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
#   labs(x = "Icesat-1 [m]", y = "GEDI_calval_simulated [m]", title = "Canopy height")+
#   annotate("text", x = lx, y = ly, label = labels, hjust = 0) 
# p
# 
# 
# # Save the plot as a TIFF file
# tiff("plot.tiff", width = 6, height = 5, units = "in", res = 300)
# p
# dev.off()


is1_raster <- raster('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/IS1_raster_1dot7km.tif')

################################################
##################
#################Compare with GEDI real shots 
#####################################################################


# Read Parquet data
#df_gedi_Parquet<- read_parquet("Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/calval_20230417.parquet") ### 10G data  ## takes time to read
## file in pc 
#df_gedi_Parquet<- read_parquet("C:/Users/lxiong/Desktop/Data/calval_20230417.parquet")

### data L2A
#  geolocation/rh_a1_098
#  lat_lowestmode
#  lon_lowestmode
colnames(df_gedi_Parquet)

#/gedi/gedi_rh_a1/gedi_rh098
#/geolocation/latitude
# /geolocation/longitude
# /gedi/selected_algorithm

gedi_shots<- cbind(df_gedi_Parquet$`/geolocation/latitude`,            #X1
                   df_gedi_Parquet$`/geolocation/longitude`,           #X2
                   df_gedi_Parquet$`/gedi/gedi_rh_a1/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a2/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a3/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a4/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a5/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a6/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a7/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a8/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a9/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/gedi_rh_a10/gedi_rh098`,
                   df_gedi_Parquet$`/gedi/selected_algorithm`, # the one use to a shot    X13
                   df_gedi_Parquet$`/gedi/optimal_algorithm`, # the one select right ground  X14
                   df_gedi_Parquet$`/quality/algorithmrun_flag`,
                   df_gedi_Parquet$`/quality/degrade_flag`, # Flag 1 indicating degraded stat, 1 true  X16
                   df_gedi_Parquet$`/quality/quality_flag`, #  (1=valid, 0=invalid)   X17
                   df_gedi_Parquet$`/quality/surface_flag`,  #   Flag 1 is surface?  
                   df_gedi_Parquet$`/gedi/sensitivity`, # should > 0.95  X19
                   df_gedi_Parquet$`/geolocation/shot_number`,   #x20
                   df_gedi_Parquet$`/simulation/als_rh098`,   #x21
                   df_gedi_Parquet$`/land_cover_data/leaf_on_doy_1km`,   #x22
                   df_gedi_Parquet$`/simulation/als_project`, # X23   char
                   #df_gedi_Parquet$`/gedi/beam`,  #strong???  
                   #df_gedi_Parquet$`beam_number`,  #strong???  
                   df_gedi_Parquet$`geolocation/delta_time`, # X24 # Time delta since Jan 1 00:00 2018.
                   df_gedi_Parquet$`/land_cover_data/leaf_off_doy_1km` #X25 leaf off start days
                   #df_gedi_Parquet$`time_of_day`,   # 
                   #df_gedi_Parquet$`date_time`     # 
                   )
###               /simulation/als_project 


gedi_shots<- data.frame(gedi_shots)
# Convert all columns to numeric     ## be careful about shot number conversion 
#gedi_shots[, 1:22] <- as.data.frame(lapply(gedi_shots[, 1:22], as.numeric)) # takes too much time
#summary(gedi_shots)

### do we have duplicate shots here ?
gedi_shots %>% filter( X20 == 19660500200096920)

## outlier
gedi_shots %>% filter( X20 == 78470600300588016)
gedi_shots %>% filter( X20 == 111830200200267280)
outlier <- gedi_shots %>% filter( X20 == 210660800200257248)





# save the data
#write.csv(gedi_shots, file = "gedi_shots.csv", row.names = FALSE)


 gedi_shots$X1 <- as.numeric(gedi_shots$X1)
 gedi_shots$X2 <- as.numeric(gedi_shots$X2)
 gedi_shots$X3 <- as.numeric(gedi_shots$X3)
######
#### omit disturbed shots, degrage shots, surface shots, etc. 

#gedi_shots <- read.csv("gedi_shots.csv")


# Filter dataframe based on a specific value
filtered_df <- gedi_shots[gedi_shots$X19 >  0.95 
                          & gedi_shots$X16 <1  
                          & gedi_shots$X17 > 0
                         # & gedi_shots$X13 == gedi_shots$X14
                          ,  ]
#### fitler by project name 
#Read ALS projet name file 
als_name <- readLines('als_name.txt')

# Remove multiple occurrences of "\t"
als_name <- gsub("\t+", "", als_name)

# Remove spaces
als_name  <- gsub("\\s", "", als_name)

als_name <- as.list(als_name)

# Search string
# Check if the list contains the search name
contains_name <- 'neon_teak2021' %in% als_name

# Print the result
print(contains_name)

# Filter rows based on names in the list
filtered_df <- filtered_df[filtered_df$X23 %in% als_name, ]


#### leaf on and off 
##### Similar to IS1 

### convert delta time to DOY 
### X24  delta ime 	s

### Time delta since Jan 1 00:00 2018

# Convert seconds to POSIXlt datetime object
datetime <- as.POSIXlt(as.numeric(filtered_df$X24), origin = "2018-01-01")

# Extract the day of the year
filtered_df['day_of_year'] <- format(datetime, "%j")


##########data is leaf on ? or off ? 

##### leaf-on vs. leaf-off filter
####  X22   vs X25
# Function to select column based on value
select_row <- function(row) {
  
  st <- row['X22']
  ed <- row['X25']
  d  <- row['day_of_year']
  
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

filtered_df['growing_flag'] <- apply(filtered_df , 1, select_row)
# Filter using logical indexing
filtered_df_leafOn<- filtered_df[filtered_df$growing_flag == TRUE, ]
filtered_df_leafOff<- filtered_df[filtered_df$growing_flag == FALSE, ]





####################### filter by selected_algorithm X13


# Function to select column based on value
select_column <- function(row) {
  
  column_name <- paste0('X',as.character(as.numeric(row["X14"])+2))   # X13
  
  
  return(row[column_name])
}

# Apply the function to each row
selected_columns <- apply(filtered_df, 1, select_column)


filtered_gedi_shots <- data.frame(cbind(filtered_df$X2, filtered_df$X1, selected_columns, 
                                        filtered_df$X21 ,filtered_df$X20,  filtered_df$growing_flag ))

colnames(filtered_gedi_shots) <- c('lon', 'lat', 'rh98', 'rh98_als', 'X20', 'leaf_flag')

head(filtered_gedi_shots)

######### filtered by forest disturbance ; get a disturbance image 

### go to GEE, sampeRegions with forest cover. 
### come back not working, 
### use h3 command 


######### read  disturbance data ####
loss <- read_parquet("Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Calval/forestLoss/loss.parquet") ### extract using gh3
loss_dt <- cbind(loss$shot_number, loss$loss)  ## 0 false 1: true
loss_dt<- data.frame(loss_dt)
#loss_dt <- as.data.frame(lapply(loss_dt, as.numeric))
colnames(loss_dt) <- c('X20', 'X23_loss')
loss_dt$X20 <- as.character(loss_dt$X20)

########## read updated disturbance data

hansen_forest_loss <- read_parquet("Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Calval/calval_hansen_forest_loss_20230425_string.parquet"
                                   ) ### shared by Tiago

hansen_forest_loss <- data.frame(hansen_forest_loss )

head(hansen_forest_loss)
# Rename the 'old_column_name' to 'new_column_name'
hansen_forest_loss <- hansen_forest_loss %>%
  rename(X20 = shot_number) ## to match previous data



############### joint 
library(dplyr)
# Merge the two Parquet files by common column

options(scipen = 999) # shotnumber is too big 
# x , y
# Remove duplicate values in the "value" column
# same shotnumber, different rows = different data;

# Merge data frames based on ID column
#merged_df <- merge(df1, df2, by = "ID", all.x = TRUE)


## here shot number is character
#filtered_gedi_shots_rmDuplicates <- distinct(filtered_gedi_shots, X20, .keep_all = TRUE) 


#loss_dt_rmDuplicates <- distinct(loss_dt, X20, .keep_all = TRUE)
##############
## A left_join() keeps all observations in x.
## inner_join() only keeps observations from x that have a matching key in y.

merged_file <- left_join(filtered_gedi_shots, hansen_forest_loss, by = "X20")  # by shot number join two data frames. 
head(merged_file)

#loss_dt %>% filter( X20 == 19660500200096920)

##########################################
##############after all filtering ---####
#########################################



dt <- merged_file[merged_file$loss < 1, ]
## remove NA rows
# Remove rows with NA values
dt  <- dt [complete.cases(dt ), ]
write.csv(dt, file = "gedi_calval_shot_filtered_06142023.csv", row.names = FALSE)



#################################################GEDI grid at 0.125 degree
merged_file <- read.csv('gedi_calval_shot_filtered_06142023.csv')
head(merged_file)




merged_file$lon <- as.numeric(merged_file$lon)
merged_file$lat <- as.numeric(merged_file$lat)
merged_file$rh98 <- as.numeric(merged_file$rh98)

############# grid gedi data to 0.125 degree grid first 

e <- extent(-180, 180, -90, 90)
#r <- raster(e, ncol=2880, nrow=1440, crs = '+proj=longlat +datum=WGS84') ### 0.125 degree

r <- raster(e, ncol=23040, nrow= 11520, crs = '+proj=longlat +datum=WGS84') ### 1/64 degree ---1.7km


# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(merged_file[, 1:2], r, merged_file[,3], fun=mean) # min, max, or mean
##### 98% height 
#x <- rasterize(filtered_gedi_shots[, 1:2], r, filtered_gedi_shots[,3], 
#               fun = function(i,...) quantile(i, probs=0.98, na.rm=T))



### export to csv 
# Convert raster to XYZ format
xyz_data <- rasterToPoints(x)
# Write the data frame to a CSV file
write.csv(xyz_data, file = "gedi_calval_shot_filtered_gridded_optimal_filterDisturbance_1dot7km.csv", row.names = FALSE)





##########
gedi_gridded <- read.csv("gedi_calval_shot_filtered_gridded_optimal_filterDisturbance_1dot7km.csv")



GEDI_raster <- raster('Z:/vclgp/xiongl/GEDIglobal/global/global_img.tif')

library(rasterVis)
# Plot the DEM using levelplot function from rasterVis
levelplot(GEDI_raster, margin = FALSE, col.regions = rainbow(100),
           at = seq(0, 4000, length.out = 100))

# Convert raster to XYZ format
xyz_data <- rasterToPoints(GEDI_raster)

x<- xyz_data[,1]
y<- xyz_data[,2]



#### given x and y ----give me is1 raster value 
# Extract raster values at the specified coordinates
#x <- gedi_gridded$x
#y <- gedi_gridded$y




#### GEDI calval, not gridded 
#gedi_not_gridded <- read.csv("gedi_calval_shot_filtered_06142023.csv")
#x <- gedi_not_gridded$lon
#y <- gedi_not_gridded$lat


#is1 <- raster('IS1_raster_14km.tif')
is1_values <- extract(is1_raster , cbind(x, y))


#is1_98percent <- raster('IS1_raster_14km_98percent.tif')
#is1_values_98percent <- extract(is1_98percent, cbind(x, y))



#### result 
#gedi_rh98 <-   gedi_not_gridded$rh98              #gedi_gridded$layer

#gedi_rh98 <-   gedi_gridded$layer


gedi_rh98 <-   xyz_data[,3]/100  # cm to m

#### swith to see

is1_values<-  extract(is1_raster , cbind(x, y))


result <- cbind(x, y,  gedi_rh98, is1_values)


#result <- cbind(x, y,  gedi_rh98, is1_values_98percent)



result <- data.frame(result)

# Remove rows with NA values
result <- result[complete.cases(result), ]




#### statistics 
#### 
N = nrow(result)
N = format(N, big.mark = ",", scientific = FALSE)     # Apply format function
# "10,000,000"

Bias = mean(result$gedi_rh98-result$is1_values)
MAE = mean(abs(result$gedi_rh98-result$is1_values))
r <- cor(result$is1_values, result$gedi_rh98)
R2 <- r^2
# Calculate RMSE
rmse <- sqrt(mean((result$is1_values - result$gedi_rh98)^2))
rmse_percent <- rmse/mean(result$gedi_rh98)*100

##### add text in fig 
# Create example data
lx <- c(0, 0, 0,  0)
ly <- c(60, 55, 50, 45)
labels <- c(paste('N =',N), 
            paste('Bias =',round(Bias,1), ', MAE =', round(MAE, 1)), 
            paste('r =',round(r, 2),', R^2 =', round(R2,2)),
            paste('rmse =', round(rmse, 1), ', %rmse =', round(rmse_percent, 1))  )


# Create a scatterplot with density
#geom_text() adds only text to the plot. 
#geom_label() draws a rectangle behind the text, making it easier to read.


p <- ggplot(result, aes(x = result$is1_values, y = result$gedi_rh98)) +
  geom_bin2d(binwidth = c(1.5, 1.5)) +
  scale_fill_continuous(type = "viridis", limits=c(0, 500),oob = scales::squish) +
  #geom_point() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed",linewidth = 2 )+
  xlim(c(0,62)) +
  ylim(c(0,62))+
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
  labs(x = "Icesat-1 [m]", y = "GEDI_calval_shots [m]", title = "Canopy height comparison (mean_gridded)")+
  annotate("rect", xmin = 0, xmax = 27, ymin = 43, ymax = 62,
             alpha = 0.9, fill = 'gray')+
  annotate("text", x = lx, y = ly, label = labels, hjust = 0) 
 
p


# Save the plot as a TIFF file
tiff("plot_filtered_gridded_98.tiff", width = 6, height = 5, units = "in", res = 300)
p
dev.off()




###############median height comparison 

## group result by IS-1  height class
##switch
#tmp1 <- result$is1_values
#tmp2 <- result$gedi_rh98

#result$gedi_rh98 <- tmp1

#result$is1_values <-tmp2

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

p2 <- ggplot(df_1,aes(x , y , color = type, fill = type)) +
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
      xlim(c(0,42)) +
      ylim(c(0,42))+
      #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
      # geom_text( aes(
      #                                 label = format(n, big.mark = ",", 
      #                                 scientific = FALSE)), vjust = -1, size = 5)+
      labs(x = "Icesat-1 [m]", y = "GEDI_calval_shots [m]", title = "Median of GEDI at 1/8° from ICESat")

p2
# Save the plot as a TIFF file
tiff("plot_median_mean_filtered_within1km.tiff", width = 6, height = 5, units = "in", res = 300)
p2
dev.off()


### Fig. box lots of diff
# Fig. box plots; Define height class intervals

result$diff <- result$gedi_rh98 - result$is1_values


height_bins <- cut(result$is1_values, breaks = seq(0.5,60.5,2))

result$height_class <- height_bins 


# Calculate the number of samples in each class
sample_counts <- result %>% count(height_class)


pbox <- ggplot(result, aes(x = height_class, y = diff)) +
  # coef = 0
  geom_boxplot(outlier.shape = NA,  fill = "lightblue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs( x = "ICESat height class [m]", y = "GEDI_rh98 - IS1_gridded [m]")+
  ylim(-15, 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
  geom_text(data = sample_counts, aes(x = height_class, y = 9, 
            label = format(n, big.mark = ",", scientific = FALSE)),angle = 45, vjust = -1, size =2)

# Save the plot as a TIFF file
tiff("plot_GEDI_IS1_boxplot_06142023.tiff", width = 6, height = 3, units = "in", res = 300)
pbox
dev.off()











##############################################################################
###############GEDI ALS comparison############################################
##############################################################################

#head(filtered_gedi_shots)
head(merged_file)

dt <- merged_file[merged_file$loss < 1, ]
## remove NA rows
# Remove rows with NA values
result <- result[complete.cases(result), ]


##########leaf on 
### leaf on 
### 
dt_leafOn<- merged_file[merged_file$leaf_flag == TRUE, ]
dt_leafOff<- merged_file[merged_file$leaf_flag == FALSE, ]



#dt <- data.frame(cbind(filtered_gedi_shots$V4, filtered_gedi_shots$selected_columns))
#colnames(dt) <- c('ALSsim', 'GEDI')

###############gridded and compare ????



e <- extent(-180, 180, -90, 90)
r <- raster(e, ncol=2880, nrow=1440, crs = '+proj=longlat +datum=WGS84') ### 0.125 degree
# you need to provide a function 'fun' for when there are multiple points per cell
x_gedi <- rasterize(dt[, 1:2], r, dt[,3], fun=mean) # min, max, or mean
x_als <- rasterize(dt[, 1:2], r, as.numeric(dt[,4]), fun=mean)
##### 98% height 
#x <- rasterize(filtered_gedi_shots[, 1:2], r, filtered_gedi_shots[,3], 
#               fun = function(i,...) quantile(i, probs=0.98, na.rm=T))
xyz_data_gedi <- rasterToPoints(x_gedi)
xyz_data_als <- rasterToPoints(x_als)

dt2 <- data.frame(cbind(xyz_data_als[,3], xyz_data_gedi[,3]))
colnames(dt2) <- c('rh98_als', 'rh98')

############################################################################################

#### statistics 
#### 
result <- dt2

result <- dt[, c('rh98_als', 'rh98')]

#result <- dt_leafOn[, c('rh98_als', 'rh98')]
#result <- dt_leafOff[, c('rh98_als', 'rh98')]


# Remove rows with NA values
result <- result[complete.cases(result), ]
result <- as.data.frame(lapply(result, as.numeric))    
colnames(result) <- c('is1_values', 'gedi_rh98')

result <- result[result$is1_values < 60, ]

N = nrow(result)
N = format(N, big.mark = ",", scientific = FALSE)     # Apply format function
# "10,000,000"

Bias = mean(result$gedi_rh98-result$is1_values)
MAE = mean(abs(result$gedi_rh98-result$is1_values))
r <- cor(result$is1_values, result$gedi_rh98)
R2 <- r^2
# Calculate RMSE
rmse <- sqrt(mean((result$is1_values - result$gedi_rh98)^2))
rmse_percent <- rmse/mean(result$gedi_rh98)*100

##### add text in fig 
# Create example data
lx <- c(0, 0, 0,  0)
ly <- c(60, 55, 50, 45)
labels <- c(paste('N =',N), 
            paste('Bias =',round(Bias,1), ', MAE =', round(MAE, 1)), 
            paste('r =',round(r, 2),', R^2 =', round(R2,2)),
            paste('rmse =', round(rmse, 1), ', %rmse =', round(rmse_percent, 1))  )
###################################################################


# c('is1_values', 'gedi_rh98')

p3_all <- ggplot(result, aes(x = is1_values, y = gedi_rh98)) +
  geom_bin2d(binwidth = c(1.5, 1.5)) +
  scale_fill_continuous(type = "viridis", limits=c(0,1000),oob = scales::squish) +
  #scale_fill_gradient(low = "white", high = "darkred", limits=c(0,800), oob = scales::squish )+
  #geom_point() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"),
        legend.position = "bottom",
        legend.key.width = unit(dev.size()[1] / 10, "inches"),
        legend.key.height = unit(dev.size()[1] / 45, "inches")
        )+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed",linewidth = 2 )+
  xlim(c(0,62)) +
  ylim(c(0,62))+
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
  labs(x = "ALSsim [m]", y = "GEDI_calval_shots [m]")+ #, title = "GEDI vs ALS simulation")+
  annotate("rect", xmin = 0, xmax = 40, ymin = 43, ymax = 62,
           alpha = 0.9, fill = 'gray')+
  annotate("text", x = lx, y = ly, label = labels, hjust = 0, size = 2.5)

p3_off

############ all on off fig.


library(ggpubr)
# Create a new plot device
tiff("plot_GEDI_ALS.tiff", width = 8, height = 3.5, units = "in", res = 300)


ggarrange(p3_all, p3_on, p3_off, 
             nrow = 1, #ncol=3,
             labels = c("a) All", "b) Leaf on", " c) Leaf off"),
             common.legend = TRUE,
             vjust = 1,
            # heights = c(2,1,1),
             legend = "bottom")
# Save the plot
dev.off()





# Save the plot as a TIFF file
tiff("plot_GEDI_ALSsim.tiff", width = 6, height = 5, units = "in", res = 300)
p3
dev.off()


## group result by IS-1  height class




# Define height class intervals
height_bins <- cut(result$is1_values, breaks = seq(1.5,60.5))

# Group data by height class and calculate mean weight
gedi_median <- aggregate(result$gedi_rh98, by = list(height_class = height_bins), FUN = median)
gedi_mean <- aggregate(result$gedi_rh98, by = list(height_class = height_bins), FUN = mean)
gedi_median_x <- as.character(gedi_median$height_class)


# Load the stringr package
library(stringr)
gedi_median_x <-  str_extract(gedi_median_x, "(?<=\\,).*?(?=\\])")
gedi_median_x <- as.numeric(gedi_median_x)
df_1 <- data.frame(gedi_median_x, gedi_median$x) # median
colnames(df_1) <- c('x', 'y')
df_1['type'] = 'median'
df_2 <- data.frame(gedi_median_x, gedi_mean$x) # mean
colnames(df_2) <- c('x', 'y')
df_2['type'] = 'mean'
# Combine the data frames
combined_df <- rbind(df_1, df_2)

combined_df <- df_1



p2 <- ggplot(combined_df,aes(x , y, color = type, shape = type)) +
  scale_color_manual(values = c("mean" = "red", "median" = "blue") )+
  #scale_fill_manual(values = c("mean" = "orange", "median" = "lightblue") )+
  scale_shape_manual(values = c("mean" = 16, "median" = 17))+
  geom_point(aes(shape=type),size =3) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        #legend.position = c(0.1, 0.8),
  )+
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed",linewidth = 1 )+
  xlim(c(0,42)) +
  ylim(c(0,42))+
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
  labs(x = "ALSsim [m]", y = "GEDI_calval_shots [m]", title = "Median of GEDI from ALS")

p2
# Save the plot as a TIFF file
tiff("plot_GEDI_ALS_median_filtered_06142023.tiff", width = 6, height = 5, units = "in", res = 300)
p2
dev.off()



#######GEDI-ALS vs ALS
combined_df <- rbind(df_1, df_2)
combined_df['gedi_als'] <- combined_df$y - combined_df$x

p_GEDI_ALS <- ggplot(combined_df,aes(x , gedi_als, color = type, fill = type)) +
  scale_color_manual(values = c("mean" = "red", "median" = "blue") )+
  #scale_fill_manual(values = c("mean" = "orange", "median" = "lightblue") )+
  scale_shape_manual(values = c("mean" = 16, "median" = 17))+
  geom_point(aes(shape=type),size =3) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        #legend.position = c(0.1, 0.8),
  )+
  geom_abline(slope = 0, intercept = 0, color = "black", linetype = "dashed",linewidth = 1 )+
  xlim(c(0,42)) +
  ylim(c(-10,10))+
  #geom_smooth(method = "lm", linewidth = 2, fullrange = TRUE, se = TRUE) +
  labs(x = "ALSsim [m]", y = "GEDI - ALS [m]")

p_GEDI_ALS 
# Save the plot as a TIFF file
tiff("plot_GEDI_ALS_ALS_06142023.tiff", width = 6, height = 3, units = "in", res = 300)
p_GEDI_ALS 
dev.off()








##### what happened to height class (0, 1]

outlier <- dt[dt$rh98_als <= 1 & dt$rh98_als > 0,]


##### GEDI - ALS vs ALS plot 
GEDI_ALS <- result$gedi_rh98 - result$is1_values

dt_boxplot <- data.frame(cbind(result$is1_values,GEDI_ALS ))


### scatter plot 

p2d <- ggplot(dt_boxplot , aes(V1, GEDI_ALS)) +
  geom_bin2d(binwidth = c(1,1)) +
  #geom_point() +
  #geom_density_2d() +
  scale_fill_continuous(type = "viridis", limits=c(0,1000),oob = scales::squish) +
  #scale_fill_gradient(low = "white", high = "darkred", limits=c(0,800), oob = scales::squish )+
  #geom_point() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = margin(0.5, 0.1, 0.1, 0.1, "cm"))+
  xlim(0, 60)+
  ylim(-15, 15)+
  labs(x = "ALSsim [m]", y = "GEDI - ALSsim [m]", title = "2D density plot (bin size (m)= 1 x 1)")
# Save the plot as a TIFF file
tiff("plot_GEDI_ALS_ALS_density.tiff", width = 6, height = 3, units = "in", res = 300)
p2d 
dev.off()



######   20 points  ################













# Fig. box plots; Define height class intervals
height_bins <- cut(dt_boxplot$V1, breaks = seq(1,60,4))

dt_boxplot$height_class <- height_bins 

ggplot(dt_boxplot, aes(x = height_class, y = GEDI_ALS)) +
  geom_boxplot() +
  labs(title = "Boxplot by Height Class", x = "Height", y = "Value")+
  ylim(-10, 10)




