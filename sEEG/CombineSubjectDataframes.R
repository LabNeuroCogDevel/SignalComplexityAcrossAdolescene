# Libraries ----

library(LNCDR)
library(data.table)
library(dplyr)
library(factoextra)
library(ggplot2)
library(e1071)
library(caret)
attach(mtcars)
library(grid)
library(gridExtra)
library(plotrix)
library(mgcv)
library(readxl)
library(lme4)
library(lubridate)
library(checkmate)
library(lmerTest)
library(tidyr)
library(jtools)

# Initialize directory and empty dataframes ----

# Set your working directory to the folder containing your CSV files
setwd('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/individualSubFiles/')

# List all CSV files in the directory
csv_files <- list.files(pattern = "MultiScaleEntropy20_referenced.csv")

# Initialize an empty data frame to store the combined data
combined_data <- data.frame()

# Loop through each CSV file, read it, and combine it ----
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Change header argument if needed
  combined_data <- rbind(combined_data, data)
}


write.csv(combined_data, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/5Mins_MSE20_allSubjects.csv', row.names = F)
