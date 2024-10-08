

library(LNCDR)
library(data.table)
library(dplyr)
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
library(tvem)
library(interactions)
library(akima)
library(mice)


entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/5Mins_MSE20_allSubjects.csv')

nonEpContacts <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/nonEPcontacts.csv')

# List of subjects
subjects <-intersect(unique(nonEpContacts$Subject), unique(entropy$Subject))


# Initialize an empty data frame to store the results
result <- data.frame()

# Loop through subjects
for (subject in subjects) {
  # Select contacts for the current subject
  selectContacts <- nonEpContacts %>%
    filter(Subject == subject) %>%
    select(nonEpContacts)
  
  # Create a formatted string
  formattedString <- paste(selectContacts$nonEpContacts, collapse = '|')
  
  # Filter rows in entropy for the current subject
  filteredData <- entropy[grep(formattedString, entropy$label), ] %>%
    subset(., !grepl('ET', label)) %>%
    subset(., !grepl('EEG', label)) %>%
    filter(Subject == subject)
  
  # Append the result to the overall result data frame
  subjectResult <- bind_rows(result, filteredData)
  
  result <- rbind(result, subjectResult)
}



