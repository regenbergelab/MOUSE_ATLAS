library(chngpt)
library(ggplot2)
library(sp)
library(lmtest)
library(plotrix)
library(ellipse)
library(plot3D)
library(raster)
library(rgl)
library(plotly)
library(openxlsx)

install.packages("segmented")
library(segmented)
setwd("C:/Users/vzl597/Desktop/Xue Analysis")
cortex3 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 4)
cortex12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 5)
cortex22 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 6)

cortex3$age <- "Young"
cortex12$age <- "Adult"
cortex22$age <- "Old"

hip3 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 7)
hip12 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 8)
hip22 <- read_excel("Mb that were found in the 3 month and 12 month liver, hippocampus and cortex.xlsx", 9)

hip3$age <- "Young"
hip12$age <- "Adult"
hip22$age <- "Old"

data <- rbind(cortex3, cortex12, cortex22)
dataH <- rbind(hip3, hip12, hip22)

data
data.frame(dataH)

ggplot(hip22, aes(x = log_norm_count_per_MB, y = log_ecc_number_per_MB, color= log_norm_count_per_MB)) +
  geom_point(color="#BFBEBE") +
  xlab("Log Normalized Count per MB") +
  ylab("Log ECC Number per MB") +
  ggtitle("LOESS Regression") +
  ylim(0, 10)+
  stat_density_2d(alpha=0.08,geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "white",bins = 20)+
  geom_smooth(mapping = aes(x = log_norm_count_per_MB, y = log_ecc_number_per_MB), method="loess")


random_data <- rnorm(5000, mean = 2.5, sd = 2)
random_data1 <- rnorm(5000, mean = 13, sd = 2)
plot(random_data,random_data1)
random <- data.frame(log_norm_count_per_MB = random_data1, log_ecc_number_per_MB = random_data)

ggplot(random, aes(x = log_norm_count_per_MB, y = log_ecc_number_per_MB, color= log_norm_count_per_MB)) +
  geom_point(color="#BFBEBE") +
  xlab("Log Normalized Count per MB") +
  ylab("Log ECC Number per MB") +
  ggtitle("LOESS Regression") +
  ylim(0, 10)+
  stat_density_2d(alpha=0.08,geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "white",bins = 20)+
  geom_smooth(mapping = aes(x = log_norm_count_per_MB, y = log_ecc_number_per_MB), method="loess")



loess <- loess(genesah$log_ecc_number_per_MB  ~ genesah$log_norm_count_per_MB)
summary(loess)


# Score the model on a grid
grid <- seq(from = 1, to = 20, length.out = 1000)
scored_values <- predict(loess, newdata = grid)
# Display the scored values in a table
score_table <- data.frame(grid, scored_values)
plot(score_table)


# Ratio (eccDNA/RNA level) through aging
genesy$Ratio <- genesy$log_ecc_number_per_MB /genesy$log_norm_count_per_MB 
genesa$Ratio <- genesa$log_ecc_number_per_MB /genesa$log_norm_count_per_MB 
geneso$Ratio <- geneso$log_ecc_number_per_MB /geneso$log_norm_count_per_MB

liver <- rbind(genesy, genesa, geneso)
liver$eccDNA <- 10^(liver$log_ecc_number_per_MB/1000000/liver$slen)
liver$Ratio <- liver$eccDNA /liver$log_norm_count_per_MB
liver$logRatio <- log(liver$Ratio,2)

ggplot(liver, aes(x = age, y= Ratio)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#999999", "#3e7b7d", "#a7dcdd")) +
  scale_y_continuous() +
  ylim(c(0.5,1.25))+
  labs(x = "Age and Condition", y = "eccDNA count / RNA level", title = "Intron density") 
