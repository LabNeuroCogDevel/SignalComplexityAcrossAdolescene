# Load Libraries 

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


outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)

# Load dataframes ----

acw <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/acw/allSubjectsACW_newACWscript.csv') %>% mutate(Condition = 'eyesClosed')
fooof <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsAllChannelsFooofMeasures_20230911.csv')


# Merge ACW and FOOOF ----

fooof_acw <- merge(acw, fooof, by = c("Subject", "Channel","Condition")) 

# Outlier detection ----

fooof_acw_outlier <- fooof_acw %>% group_by(Subject) %>%
  mutate(across(c("ACW_old","ACW_0","ACW_50", "Exponent", "Offset"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


# ACW old vs new 

lunaize(ggplot(data = fooof_acw_outlier, aes(x = ACW_50, y = ACW_old, color = Channel)) + 
          geom_smooth(aes(group = Channel), method="lm", alpha=0.02, linewidth=2))

# ACW vs age ----

lunaize(ggplot(data = fooof_acw_outlier, aes(x = age, y = ACW_50)) + 
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.02, linewidth=2))


gam.model <-  mgcv::gamm(ACW ~ s(age, k = 3) + Channel, data = fooof_acw_outlier, random=list(lunaID=~1))
summary(gam.model$gam)


# ACW vs FOOOF ----

lunaize(ggplot(data = fooof_acw_outlier, aes(x = Exponent, y = ACW_50)) + 
          geom_smooth(aes(group = 1), method="lm", alpha=0.02, linewidth=2))


gam.model <-  mgcv::gamm(ACW ~ s(age, k = 3) + Channel, data = fooof_acw_outlier, random=list(lunaID=~1))
summary(gam.model$gam)


## Assign broadmanns areas ----



