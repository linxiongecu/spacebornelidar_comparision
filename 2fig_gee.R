### read csv from GEE result
dt <- read.csv('c:/Users/lxiong/Downloads/ee-chart (11).csv')
colnames(dt) <- c('is1_values','gedi_rh98')


result <- data.frame(dt)

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
  scale_fill_continuous(type = "viridis", limits=c(0, 100),oob = scales::squish) +
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



###############median height comparison 

## group result by IS-1  height class
##switch
#tmp1 <- result$is1_values
#tmp2 <- result$gedi_rh98

#result$gedi_rh98 <- tmp1

#result$is1_values <-tmp2

#### 
#result <- cbind(samples_df$IS1_h,samples_df$GEI_h)
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
  labs(x = "Icesat-1 [m]", y = "GEDI_shots [m]", title = "Median of GEDI at 1/8° from ICESat")

p2

