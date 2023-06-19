#### Read calval database
#### Export x, y and rh100



# Load required packages
library(sp)
library(rgdal)


##### Read FSBD data base  data 
##### Step 1 get list of data files 
#### database path:  
directory <- "Z:/vclgp/data/gedi/l4_calval/gedi_fsbd"
#### 
# List all the files in the directory
file_list <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)


# Create an empty dataframe
df <- data.frame()

###############loop through file list 
for (i in 1:length(file_list)){
  f = file_list[i]
  #print(f)
  data <- readRDS(f)
  # Check if the dataframe contains a specific column
  if ("fpdata" %in% names(data)) {
    print(paste("The dataframe contains the column:", "fpdata"))
    d <- data$fpdata
    ### utm east, north, rh98, epsg, 
    d1 <- cbind(d$g.x, d$g.y,d$l.epsg,d$rhReal98, d$date, d$obstime)
    ### convert utm to wgs84
    d1<- data.frame(d1)
    ###
    x <- as.numeric(d1$X1)
    y <- as.numeric(d1$X2)
    ### one field site, the epsg should be the same;
    epsg <- as.numeric(d1$X3[1])
    # Create a SpatialPoints object with UTM coordinates and set the UTM CRS
    utm_points <- SpatialPoints(
      coords = data.frame(x,y),
      proj4string = CRS(paste0("+init=epsg:", epsg))
    )
    # Transform UTM coordinates to WGS84 using the WGS84 EPSG code
    wgs84_points <- spTransform(utm_points, CRS(paste0("+init=epsg:", 4326)))
    
    # Extract the WGS84 coordinates
    d1$lon <- coordinates(wgs84_points)[, 1]
    d1$lat <- coordinates(wgs84_points)[, 2]
    
    df <- rbind(df, d1)
  
  }
  
}

colnames(df) <- c('east', 'north', 'epsg', 'rhReal98', 'date', 'time', 'lon', 'lat')



# Write the dataframe to a CSV file
write.csv(df, "GEDIcalvalMetrics.csv", row.names = FALSE)


dt <- read.csv("GEDIcalvalMetrics.csv")






############################
##########   TEST ##########
############################
############################

# Load the RDS file

data <- readRDS("gedicalval_usa_neonsoap_20220623_r03.rds")
test <- data$fpdata



#In the fpdata dataframe, 
#g.x and g.y give the footprint coordinates. 
#rhReal100 is the preferred version of rh100



# Explore the data
print(data)



### test convertion utm to wgs84 here

###
x <- as.numeric(d1$X1)
y <- as.numeric(d1$X2)
### one field site, the epsg should be the same;
epsg <- as.numeric(d1$X3[1])
# Create a SpatialPoints object with UTM coordinates and set the UTM CRS
utm_points <- SpatialPoints(
  coords = data.frame(x,y),
  proj4string = CRS(paste0("+init=epsg:", epsg))
)
# Transform UTM coordinates to WGS84 using the WGS84 EPSG code
wgs84_points <- spTransform(utm_points, CRS(paste0("+init=epsg:", 4326)))

# Extract the WGS84 coordinates
d1$lon <- coordinates(wgs84_points)[, 1]
d1$lat <- coordinates(wgs84_points)[, 2]




