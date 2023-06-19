


##########
gedi_gridded <- read.csv("gedi_calval_shot_filtered_gridded_98.csv")


#### given x and y ----give me is1 raster value 
# Extract raster values at the specified coordinates
x <- gedi_gridded$x
y <- gedi_gridded$y

is1_98percent <- raster('IS1_raster_14km_98percent.tif')
is1_values_98percent <- extract(is1_98percent, cbind(x, y))



#### result 

gedi_rh98 <- gedi_gridded$layer
result <- cbind(x, y,  gedi_rh98, is1_values_98percent)
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
lx <- c(0, 0, 0, 10, 0)
ly <- c(60, 55, 50, 50, 45)
labels <- c(paste('N =',N), paste('Bias =',round(Bias,1), 'MAE =', round(MAE, 1)), 
            paste('r =',round(r, 2)), paste('R^2 =', round(R2,2)),
            paste('rmse =', round(rmse, 1), '%rmse =', round(rmse_percent, 1))  )


# Create a scatterplot with density
#geom_text() adds only text to the plot. 
#geom_label() draws a rectangle behind the text, making it easier to read.


p <- ggplot(result, aes(x = result$is1_values, y = result$gedi_rh98)) +
  geom_bin2d(binwidth = c(1.5, 1.5)) +
  scale_fill_continuous(type = "viridis", limits=c(0,20),oob = scales::squish) +
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




# Define height class intervals
height_bins <- cut(result$is1_values, breaks = seq(0,60))

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


p2 <- ggplot(combined_df,aes(x , y , color = type, fill = type)) +
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
  labs(x = "Icesat-1 [m]", y = "GEDI_calval_shots [m]", title = "Median and mean of GEDI at 0.125° from ICESat")

p2
# Save the plot as a TIFF file
tiff("plot_median_mean_filtered_98.tiff", width = 6, height = 5, units = "in", res = 300)
p2
dev.off()