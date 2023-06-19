### get ISS track data
# https://spotthestation.nasa.gov/trajectory_data.cfm
#install.packages("asteRisk")  # For reading OEM data

#install.packages("ggplot2")  # For plotting
library(asteRisk)
library(ggplot2)

testOEM_ISS <- readOEM('Z:/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/ISS.OEM_J2K_EPH.txt')
testOEM_ISS$header
testOEM_ISS$dataBlocks[[1]]$objectName
testOEM_ISS$dataBlocks[[1]]$referenceFrame
head(testOEM_ISS$dataBlocks[[1]]$ephemerides)

# GCRFtoLATLON Convert coordinates from GCRF to geodetic latitude, longitude and altitude
# GCRFtoLATLON(position_GCRF, dateTime, degreesOutput=TRUE)

position_GCRF <- data.frame(
  cbind(
    testOEM_ISS$dataBlocks[[1]]$ephemerides$position_X * 1000,
    testOEM_ISS$dataBlocks[[1]]$ephemerides$position_Y * 1000,
    testOEM_ISS$dataBlocks[[1]]$ephemerides$position_Z * 1000
    
  )
)
colnames(position_GCRF) <- c('X', 'Y', 'Z')
dateTime <- data.frame(cbind(testOEM_ISS$dataBlocks[[1]]$ephemerides$epoch))

# Load the lubridate package
library(lubridate)
# Parse the date and time string
datetime <- ymd_hms(dateTime$cbind.testOEM_ISS.dataBlocks..1...ephemerides.epoch.)

#install.packages("asteRiskData", repos="https://rafael-ayala.github.io/drat/")
library(asteRiskData)

getLatestSpaceData(targets="all")


result <- data.frame()


for (i in 1: nrow(position_GCRF)) {
  row <- position_GCRF[i, ]# Access current row using index 'i'
  print(row)
  # Convert dataframe to array
  array <- as.matrix(row)
  t <- datetime[i]
  if (i  == 4128) {t <- "2023-06-28 00:00:00 UTC"}
  if (i  == 4488) {t <- "2023-06-29 00:00:00 UTC"}
  if (i  == 4848) {t <- "2023-06-30 00:00:00 UTC"}
  if (i  == 5208) {t <- "2023-07-01 00:00:00 UTC"}
  #print(t)
  out <- GCRFtoLATLON(array , t, degreesOutput=TRUE)
  #print(out)
  result <- rbind(result, out)
  
}


# Create a scatter plot
ggplot(result, aes(x =X.119.716310990453 , y = X9.38681001909205, color = X414317.231569037)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "Longitude", y = "Latitude", color = "Height") +
  theme_minimal()


# Create a histogram of elevation
ggplot(result, aes(x = X414317.231569037)) +
  geom_histogram(binwidth = 100, fill = "blue", color = "black") +
  labs(x = "Elevation", y = "Count") +
  theme_minimal()

#### time and elevation 

dt <- data.frame(cbind(dateTime, result$X414317.231569037))
colnames(dt) <- c('t', 'h')  
# Subset the first 1000 columns
dt <- dt[ 4000:4100, ]
ggplot(dt, aes(x =t , y = h, color = h)) +
  geom_point(size=5) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "t", y = "h", color = "Height") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



