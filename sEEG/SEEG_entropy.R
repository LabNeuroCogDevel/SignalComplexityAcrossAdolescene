

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
library(eegUtils)
library(tvem)
library(interactions)
library(akima)
library(mice)
library(gratia)


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 1.5))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

## MultiScale entropy ----

MSentropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/5Mins_MSE20_allSubjects.csv')
MSentropy <- separate(MSentropy, label, into = c("Letter", "Number"), sep = "_", convert = TRUE) 
MSentropy$Number <- sub("^0", "", MSentropy$Number)

## All electrodes, not averaged 
MSentropy_long <- MSentropy %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)")

MSentropy_long <- unite(MSentropy_long, "label", Letter, Number, sep = "_") # need the labels to match the TEI since i removed the 0 in front of their signal digits 

# Find each subjects max entropy value on the MSE curve 
subs = unique(MSentropy_long$Subject)

maxValues <- data.frame()

for (subject in subs) {
  
  subDataFrame <- MSentropy_long %>% filter(Subject == subject)
  chans = unique(subDataFrame$label)
  
  for (chan in chans) {
    # Check if the subject has the visit number in the dataset
    
    subChanData <- subDataFrame %>% filter(label == chan)
    maxSubchanEntropy <- max(subChanData$MSx)
    
    subChanInfo <- data.frame(Subject = subject, label = chan, maxEntropy = maxSubchanEntropy)
    maxValues <- rbind(maxValues, subChanInfo)
  }
}

maxValuesAge <- merge(maxValues, MSentropy_long, by = c("Subject", "label")) %>% merge(., agefile, by = "Subject")

lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(Subject=~1))
summary(gam.model$gam)


lunaize(ggplot(data = maxValuesAge, aes(x = age, y = Var1)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = maxValuesAge, random=list(Subject=~1))
summary(gam.model$gam)



### Parse out non EP electrodes ----
nonEpContacts <- read_excel('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/nonEPcontacts.xlsx')

# List of subjects
subjects <-intersect(unique(nonEpContacts$Subject), unique(MSentropy$Subject))

result <- data.frame()

# Loop through subjects
for (subject in subjects) {
  
  # Initialize an empty data frame to store the results
  letters_list <- c()
  numbers_list <- c()
  
  
  # Select contacts for the current subject
  selectContacts <- nonEpContacts %>%
    filter(Subject == subject) %>%
    select(nonEpContacts)
  
  # Create a formatted string
  formattedString <- paste(selectContacts$nonEpContacts, collapse = '|')
  
  # Split the string by pipe (|)
  segments <- strsplit(formattedString, "\\|")[[1]]
  
  # Function to expand ranges like "1-5" to individual numbers 1, 2, 3, 4, 5
  expand_range <- function(range) {
    range_parts <- as.numeric(strsplit(range, "-")[[1]])
    if (length(range_parts) == 2 && all(!is.na(range_parts))) {
      return(range_parts[1]:range_parts[2])
    } else {
      warning("Invalid range:", range)
      return(NULL)
    }
  }
  
  # Process each segment
  for (segment in segments) {
    # Split each segment into letter and numbers
    parts <- strsplit(segment, " ")[[1]]
    if (length(parts) > 1) {
      letter <- parts[1]
      ranges <- parts[-1]
      
      # Expand ranges and append to the lists
      expanded_ranges <- unlist(sapply(ranges, expand_range))
      
      if (!any(is.null(expanded_ranges))) {
        letters_list <- c(letters_list, rep(letter, length(expanded_ranges)))
        numbers_list <- c(numbers_list, expanded_ranges)
      }
    }
  }
  
  # Create a data frame with letters and numbers
  letterNumberList <- data.frame(
    Letter = letters_list,
    Number = numbers_list
  )
  
  
  filteredData <- merge(MSentropy, letterNumberList, by = c("Letter", "Number")) %>%
    subset(., !grepl('EEG', Letter)) %>%
    filter(Subject == subject)
  
  
  result <- rbind(result, filteredData)
}


nonEpAge <- merge(result, agefile, by = "Subject")

subject_averages <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, MSx11, MSx12, MSx13, MSx14, MSx15, MSx16, MSx17, MSx18, MSx19, MSx20, Var1, age)~ Subject, data = nonEpAge, FUN = mean)

ggplot(data = nonEpAge, aes(x = age, y = MSx7)) + geom_point() + 
  geom_smooth(aes(group = 1, color = "black"), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2) 

lunaize(ggplot(data = subject_averages, aes(x = age, y = Var1)) + geom_point() + 
  geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + 
    ylab("AUC")


nonEpAge_long <- nonEpAge %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>% mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults')))) 

lunaize(ggplot(data = nonEpAge_long, aes(x = as.numeric(timeScale), y = MSx, color = ageGroup))  + 
  geom_smooth(aes(group = ageGroup), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2) )

subAvg_long <- subject_averages %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>% mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults')))) 


lunaize(ggplot(data = subAvg_long, aes(x = as.numeric(timeScale), y = MSx, color = ageGroup))  + 
          geom_smooth(aes(group = ageGroup), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


### All non EP Electrodes, not averaged ----
nonEpAge_long <- nonEpAge %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>% mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))

nonEpAge_long <- unite(nonEpAge_long, "label", Letter, Number, sep = "_") # need the labels to match the TEI since i removed the 0 in front of their signal digits 

# Find each subjects max entropy value on the MSE curve 
subs = unique(nonEpAge_long$Subject)

maxValues <- data.frame()

for (subject in subs) {
  
  subDataFrame <- nonEpAge_long %>% filter(Subject == subject)
  chans = unique(subDataFrame$label)
  
  for (chan in chans) {
    # Check if the subject has the visit number in the dataset
    
    subChanData <- subDataFrame %>% filter(label == chan)
    maxSubchanEntropy <- max(subChanData$MSx)
    
    subChanInfo <- data.frame(Subject = subject, label = chan, maxEntropy = maxSubchanEntropy)
    maxValues <- rbind(maxValues, subChanInfo)
  }
}

maxValuesAge <- merge(maxValues, nonEpAge_long, by = c("Subject", "label"))

maxValuesAge_outlier <- maxValuesAge %>%
  mutate(across(c("maxEntropy", "Var1"), naoutlier))

lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(Subject=~1))
summary(gam.model$gam)


lunaize(ggplot(data = maxValuesAge, aes(x = age, y = Var1)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = maxValuesAge, random=list(Subject=~1))
summary(gam.model$gam)



### Average across all non EP electrodes ----
subject_averages <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, MSx11, MSx12, MSx13, MSx14, MSx15, MSx16, MSx17, MSx18, MSx19, MSx20, Var1, age)~ Subject, data = nonEpAge, FUN = mean)

subject_averages_long <- subject_averages %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))

# Find each subjects max entropy value on the MSE curve 
subs = unique(subject_averages_long$Subject)

maxValues <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset

      subData <- subject_averages_long %>% filter(Subject == subject)
      maxSubEntropy <- max(subData$MSx)
      
      subInfo <- data.frame(Subject = subject, maxEntropy = maxSubEntropy)
      maxValues <- rbind(maxValues, subInfo)
}

maxValuesAge <- merge(maxValues, agefile, by = c("Subject"))

write.csv(maxValuesAge, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/sEEG/Results/allSubjectsMaxMSE_referenced.csv', row.names = F)

lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge)
summary(gam.model$gam)


lunaize(ggplot(data = subject_averages, aes(x = age, y = Var1)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = subject_averages)
summary(gam.model$gam)


## sEEG FOOOF ----

fooof <- read.csv('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/fooof/allSubjectsFOOOFmeasures2.csv')
fooof <- separate(fooof, Channel, into = c("Letter", "Number"), sep = "_", convert = TRUE) 
fooof$Number <- sub("^0", "", fooof$Number)


### Parse out non EP electrodes ----
rois <- read.csv('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE_Lab/entropy/allSubjects_rois.csv')
rois$rois5 <- rois$rois6
rois <- separate(rois, rois6, into = c("Letter", "Number"), sep = "_", convert = TRUE) 
rois$Number <- sub("^0", "", rois$Number)
rois <- rename(rois, roi = rois1, x = rois2, y = rois3, z = rois4, Subject = rois7)

# List of subjects
subjects <-intersect(unique(nonEpContacts$Subject), unique(fooof$Subject))

result <- data.frame()


# Loop through subjects
for (subject in subjects) {
  
  # Initialize lists to store letters and numbers
  letters_list <- c()
  numbers_list <- c()
  
  # Select contacts for the current subject
  selectContacts <- nonEpContacts %>%
    filter(Subject == subject)
  
  # Process each contact segment
  for (segment in selectContacts$nonEpContacts) {
    # Split each segment into letter and numbers part
    parts <- strsplit(segment, " ")[[1]]
    if (length(parts) > 1) {
      letter <- parts[1]
      ranges_string <- paste(parts[-1], collapse = " ")  # Join remaining parts
      
      # Split ranges by comma to handle multiple ranges
      ranges <- strsplit(ranges_string, ", ")[[1]]
      
      # Function to expand ranges like "1-5" to individual numbers 1, 2, 3, 4, 5
      expand_range <- function(range) {
        range_parts <- as.numeric(strsplit(range, "-")[[1]])
        if (length(range_parts) == 2 && all(!is.na(range_parts))) {
          return(range_parts[1]:range_parts[2])
        } else {
          warning("Invalid range:", range)
          return(NULL)
        }
      }
      
      # Process each range in the segment
      for (range in ranges) {
        # Expand range into individual numbers
        expanded_numbers <- expand_range(range)
        if (!is.null(expanded_numbers)) {
          # Append letter and expanded numbers to lists
          letters_list <- c(letters_list, rep(letter, length(expanded_numbers)))
          numbers_list <- c(numbers_list, expanded_numbers)
        }
      }
    }
  }
  
  # Create a data frame with letters and numbers
  letterNumberList <- data.frame(
    Letter = letters_list,
    Number = numbers_list
  )
  
  # Merge with 'fooof' data and filter
  filteredData <- merge(fooof, letterNumberList, by = c("Letter", "Number")) %>%
    subset(., !grepl('EEG', Letter)) %>%
    filter(Subject == subject)
  
  # Bind filtered data to result
  result <- rbind(result, filteredData)
}


result <- merge(result, rois, by = c("Subject", "Letter", "Number")) %>% subset(., !grepl('Background', roi)) 

fooof_nonEpAge <- merge(result, agefile, by = "Subject")

write.csv(fooof_nonEpAge, file = '/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/fooof/allSubjects_fooof_nonEPelectrodes.csv', row.names = F)


#### Looking at regions ----
fooof_nonEpAge_group <- fooof_nonEpAge %>%
  mutate(
    region = ifelse(grepl("\\(GapMap\\)", roi),  # Check if "(GapMap)" is present
                   gsub("\\s*\\(GapMap\\).*", "", roi),  # Extract text before "(GapMap)" and remove any leading whitespace
                   gsub(".*\\((.*?)\\).*", "\\1", roi)),  # Extract text within parentheses
    side = substr(roi, nchar(roi), nchar(roi))  # Extract last character (L or R)
  ) %>% mutate(ageGroup = as.factor(cut(age, c(0,17,Inf), labels = c('Adol','Adults'))))


df_by_region <- fooof_nonEpAge_group %>%
  group_by(region) %>%
  summarise(df = n_distinct(age) - 1)  # Calculate df based on unique age values

# Set a threshold for minimum degrees of freedom
min_df_threshold <- 3  # Adjust this threshold as needed

# Filter regions based on degrees of freedom
regions_with_sufficient_df <- df_by_region %>%
  filter(df >= min_df_threshold) %>%
  pull(region)

# Filter broadbandPower data to include only regions with sufficient df
filtered_data <- fooof_nonEpAge_group %>%
  filter(region %in% regions_with_sufficient_df)

lunaize(ggplot(data = filtered_data, 
               aes(x = age, y = Exponent)) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.05, linewidth=2)) +
  ggtitle("Exponent All Regions with DF > 3") + facet_wrap(~region)


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3) + region, data = filtered_data, random=list(Subject=~1))
summary(gam.model$gam)


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3), data = filtered_data %>% filter(region == 'Temporal-to-Parietal'), random=list(Subject=~1))
summary(gam.model$gam)



lunaize(ggplot(data = filtered_data, 
               aes(x = age, y = Offset)) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.05, linewidth=2)) +
  ggtitle("Offset All Regions with DF > 3") + facet_wrap(~region)

gam.model <-  mgcv::gamm(Offset ~ s(age, k = 3), data = filtered_data %>% filter(region == 'STS'), random=list(Subject=~1))
summary(gam.model$gam)



(ggplot(data = fooof_nonEpAge_group %>% group_by(region, ageGroup) %>% summarize(Exponent = mean(Exponent, na.rm=T)), aes(x = reorder(region, Exponent), y = Exponent, fill = ageGroup)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(x = "Regions", y = "Exponent", fill = "Age Group") +
  scale_fill_manual(values = c("Adol" = "skyblue", "Adults" = "orange")) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels by 45 degrees and align to the right
)


lunaize(ggplot(data = fooof_nonEpAge_group %>% group_by(region) %>% summarize(Exponent = mean(Exponent, na.rm=T)), aes(x = reorder(region, Exponent), y = Exponent)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, fill = "skyblue") +
    labs(x = "Regions", y = "Exponent"))+
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels by 45 degrees and align to the right



lunaize(ggplot(data = fooof_nonEpAge_group %>% group_by(region) %>% summarize(Offset = mean(Offset, na.rm=T)), aes(x = reorder(region, Offset), y = Offset)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.8, fill = "skyblue") +
          labs(x = "Regions", y = "Offset"))+
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Set y-axis limits to start from zero or the minimum value
  # Rotate labels by 45 degrees and align to the right


fooof_nonEpAge_outlier <- fooof_nonEpAge %>% group_by(Subject) %>%
  mutate(across(c("Exponent", "Offset"), naoutlier))



lunaize(ggplot(data = fooof_nonEpAge_group, aes(x = age, y = Exponent)) + geom_point() +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


lunaize(ggplot(data = fooof_nonEpAge_group, aes(x = age, y = Exponent)) + geom_point() +
          geom_smooth(aes(group = 1), method="lm", alpha=0.4, linewidth=2))


lunaize(ggplot(data = fooof_nonEpAge, aes(x = age, y = Error)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

lunaize(ggplot(data = fooof_nonEpAge, aes(x = age, y = Offset)) + geom_point() +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) 


## Merge MSE and FOOOF ----
fooof_nonEpAge <- unite(fooof_nonEpAge, "label", Letter, Number, sep = "_") # need the labels to match the TEI since i removed the 0 in front of their signal digits 
fooofMSE_allElectrodes <- merge(fooof_nonEpAge, maxValuesAge_outlier, by = c("Subject", "label", "age"))


lunaize(ggplot(data = fooofMSE_allElectrodes, aes(x = maxEntropy, y = Exponent)) + geom_point() +
          geom_smooth(aes(group = 1), method="loess", alpha=0.4, linewidth=2))


lunaize(ggplot(data = fooofMSE_allElectrodes, aes(x = Var1, y = Exponent)) + geom_point() +
          geom_smooth(aes(group = 1), method="loess", alpha=0.4, linewidth=2))


lunaize(ggplot(data = fooofMSE_allElectrodes, aes(x = maxEntropy, y = Offset)) + geom_point() +
          geom_smooth(aes(group = 1), method="loess", alpha=0.4, linewidth=2))

lunaize(ggplot(data = fooofMSE_allElectrodes, aes(x = Var1, y = Offset)) + geom_point() +
          geom_smooth(aes(group = 1), method="loess", alpha=0.4, linewidth=2))


### Average non EP electrodes ----
fooofMSE_averages <- aggregate(cbind(Exponent, Offset, maxEntropy, Var1, age)~ Subject, data = fooofMSE_allElectrodes, FUN = mean)


lunaize(ggplot(data = fooofMSE_averages, aes(x = maxEntropy, y = Exponent)) + geom_point() +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

lunaize(ggplot(data = fooofMSE_averages, aes(x = maxEntropy, y = Offset)) + geom_point() +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

#### By timescale ----

MSEfooof <- merge(fooof_nonEpAge, nonEpAge_long, by = c("Subject", "label", "age"))%>% mutate(timeScale = (as.numeric(timeScale)))


lunaize(ggplot(data = MSEfooof %>% filter(timeScale >=14 & timeScale<=20), aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))

lunaize(ggplot(data = MSEfooof  %>% filter(timeScale <=10), aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))


model <- lmer(MSx ~ Exponent + age + timeScale + (1|Subject), data = MSEfooof )
summary(model)

model <- lmer(MSx ~ Exponent *age + timeScale + (1|Subject), data = MSEfooof)
summary(model)

interactionmodel <- lmer(MSx ~ ageGroup + Exponent*timeScale + (1|Subject), data = MSEfooof )
summary(model)

prediction<- ggpredict(interactionmodel, terms = c("Exponent", "timeScale[2,4,6,8,10,14,16,18,20]", "ageGroup"))

plot(prediction)

ggeffects::hypothesis_test(interactionmodel,  terms = c("Exponent", "timeScale[14,15,16, 17, 18, 19,20]", "ageGroup"), test = NULL)


#### aggregate electrodes ----

fooofMSE_averages <- aggregate(cbind(Exponent, Offset, Var1,MSx, age)~ Subject + timeScale, data = fooofMSE_allElectrodes, FUN = mean) %>% mutate(ageGroup = cut(age, c(0,14,17,20,23,Inf), labels = c('10-13','14-16','17-19','20-22','23-30')))%>% mutate(timeScale = (as.numeric(timeScale)))


lunaize(ggplot(data = fooofMSE_averages %>% filter(timeScale >=14 & timeScale<=20), aes(x = Exponent, y = MSx))+ geom_point()+
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))

lunaize(ggplot(data = fooofMSE_averages  %>% filter(timeScale <=10), aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))



model <- lmer(MSx ~ Exponent + age + timeScale + (1|Subject), data = fooofMSE_averages )
summary(model)

model <- lmer(MSx ~ Exponent *age + timeScale + (1|Subject), data = fooofMSE_averages)
summary(model)

interactionmodel <- lmer(MSx ~ ageGroup + Exponent*timeScale + (1|Subject), data = fooofMSE_averages )
summary(model)

prediction<- ggpredict(interactionmodel, terms = c("Exponent", "timeScale[2,4,6,8,10,14,16,18,20]", "ageGroup"))

plot(prediction)

ggeffects::hypothesis_test(interactionmodel,  terms = c("Exponent", "timeScale[14,15,16, 17, 18, 19,20]", "ageGroup"), test = NULL)

