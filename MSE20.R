
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

outliers <- function(x) {
  (abs(x - mean(x, na.rm= T)) > (sd(x, na.rm= T) * 2))
}
naoutlier <- function(x) ifelse(outliers(x), NA, x)


# MSE Entropy 20 ----

entropy <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_allChans_MSE20.csv') %>% select(c(-type))

# Entropy Outlier Detection ----
entropy_outlier <- entropy %>% group_by(Subject) %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                  "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% ungroup()


merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
ageValues <- merge7t[c("lunaid","eeg.date","visitno","eeg.age")]
colnames(ageValues) <- c("lunaID", "visitDate","visitno", "age")
ageValues$Subject <- paste(ageValues$lunaID, ageValues$visitDate, sep = "_")

entropyAge <- merge(entropy_outlier, ageValues, by = "Subject")

# Entropy averaged across all electrodes ----
entropyAgeAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, MSx11, MSx12, MSx13, MSx14, MSx15, MSx16, MSx17, MSx18, MSx19, MSx20, Var1, age) ~ Subject + visitno, data = entropyAge, FUN = mean)
write.csv(entropyAgeAvg, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_MSE20_allChansAvg.csv')



entropyAgeAvg_outlier <- entropyAgeAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                  "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = MSx6)) + geom_point() + 
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),
               names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,13,16,19,Inf), labels = c('10-12','13-15','16-18','Adults'))) %>%
  mutate(timeScale = as.numeric(timeScale))


lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = as.numeric(timeScale), y = MSx)) + geom_point(aes(color = ageGroup)) +
          geom_smooth(aes(group = ageGroup, color = ageGroup), method=mgcv::"gam", 
                      formula = y ~ s(x, k = 5, fx = T), alpha=0.4, linewidth=1)) + xlab("Time Scale")



entropy_with_se <- entropyAgeAvg_outlier_long %>%
  group_by(timeScale, ageGroup) %>%
  summarise(
    MSx_mean = mean(MSx, na.rm = TRUE),
    MSx_se = sd(MSx, na.rm = TRUE) / sqrt(n())
  )

lunaize(ggplot(data = entropy_with_se, aes(x = as.factor(timeScale), y = MSx_mean)) + geom_point(aes(color=ageGroup), position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = MSx_mean - MSx_se, ymax = MSx_mean + MSx_se, color = ageGroup), position = position_dodge(width = 0.5), width = 0.2)) + 
  xlab("Time Scale") + 
  ylab("Mean MSx")

lunaize(ggplot(data = entropyAgeAvg_outlier, aes(x = age, y = MSx13)) + geom_point()+
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(aes(group = 1), method="loess", alpha=0.2, linewidth=1))


lunaize(ggplot(data = entropyAgeAvg_outlier_long %>% filter(timeScale >14 & timeScale <=20), aes(x = age, y = MSx, color = timeScale)) +
          geom_smooth(aes(group = timeScale, color = timeScale), method="loess", alpha=0.2, linewidth=1))

lunaize(ggplot(data = entropyAgeAvg_outlier_long %>% filter(timeScale >=1 & timeScale <=5), aes(x = age, y = MSx, color = timeScale)) +
          geom_smooth(aes(group = timeScale, color = timeScale), method="loess", alpha=0.2, linewidth=1))


lunaize(ggplot(data = entropyAgeAvg_outlier_long %>% arrange(timeScale), 
               aes(x = as.numeric(timeScale), y = MSx,
                   group = interaction(lunaID, visitno), color = age)) + 
          geom_line(alpha = 0.4)) + xlab("Time Scale")

lunaize(ggplot(data = entropyAgeAvg_outlier_long, aes(x = age, y = Var1)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AUC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3), data = entropyAgeAvg_outlier_long, random=list(lunaID=~1))
summary(gam.model$gam)

## predict timescale by MSE for each ageGroup ----
entropyAgeAvg_outlier_long <- entropyAgeAvg_outlier_long %>% mutate(timeScale_factor = as.factor(timeScale))

interactionmodel <- lmer(MSx ~ ageGroup*(timeScale_factor) + (1|lunaID), data = entropyAgeAvg_outlier_long)
summary(interactionmodel)

prediction<- ggpredict(interactionmodel, terms = c("timeScale_factor", "ageGroup"))
plot(prediction)

lunaize(ggplot(data = prediction %>% filter(x %in% seq(2, 20, by = 2)), 
       aes(x = x, y = predicted, group = group, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size =2) +  # Dodge the points
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.2, 
                position = position_dodge(width = 0.5)) ) +xlab("Time Scale") +ylab("predicted MSE")


ggeffects::hypothesis_test(interactionmodel,  terms = c("timeScale_factor[2,18]", "ageGroup"))



# Find each subjects max entropy value on the MSE curve ----
subs = unique(entropyAgeAvg_outlier_long$lunaID)
visits = unique(entropyAgeAvg_outlier_long$visitno)


maxValues <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  for (visitNum in visits) {
    if (any(entropyAgeAvg_outlier_long$lunaID == subject & entropyAgeAvg_outlier_long$visitno == visitNum)) {
      subData <- entropyAgeAvg_outlier_long %>% filter(lunaID == subject) %>% filter(visitno == visitNum)
      if (all(is.na(subData$MSx))) {
        next
      }
      
      maxSubEntropy <- max(subData$MSx, na.rm = T)
      timescale <- subData$timeScale[which(subData$MSx == maxSubEntropy)]
      
      subInfo <- data.frame(lunaID = subject, visitno = visitNum, maxEntropy = maxSubEntropy, maxTimeScale = timescale)
      maxValues <- rbind(maxValues, subInfo)
    }
  }
}

maxValuesAge <- merge(maxValues, ageValues, by = c('lunaID', "visitno"))
write.csv(maxValues, file = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Entropy/Results/allSubjects_maxEntropy_allChansAvg.csv')

lunaize(ggplot(data = maxValuesAge, aes(x = age, y = maxEntropy)) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(lunaID=~1))
summary(gam.model$gam)


lunaize(ggplot(data = maxValuesAge, aes(x = age, y = as.numeric(maxTimeScale))) + geom_point() +
          geom_line(aes(group= interaction(lunaID)), alpha = 0.2) +
          geom_smooth(group = 1, method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy TimeScale")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(lunaID=~1))
summary(gam.model$gam)

# MSE Entropy vs FOOOF ----

## Load and avg fooof across all electrodes ----
fooofAllchans <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Aperiodic_MRS_Development/Results/allSubjectsAllChannelsFooofMeasures_20230911.csv')

fooofAllchans <- fooofAllchans %>% 
  rename(labels = Channel)

fooofAvg <- aggregate(cbind(Offset, Exponent) ~ Subject + Condition, data = fooofAllchans, FUN = mean)

fooofAvg_outlier <- fooofAvg %>% group_by(Condition) %>%
  mutate(across(c("Offset", "Exponent"), naoutlier)) %>% ungroup() %>% separate(Subject, c('lunaID','vdate'))

## Merge MSE and fooof ----

MSEfooof <- merge(fooofAvg_outlier, entropyAgeAvg_outlier_long, by = c('lunaID','vdate')) %>%
  merge(.,  maxValuesAge, by = c("lunaID", "visitno", "age")) %>% mutate(timeScale = (as.numeric(timeScale)))


lunaize(ggplot(data = MSEfooof %>% filter(Condition == 'eyesOpen') %>% filter(timeScale >=14 & timeScale<=20), aes(x = Exponent, y = MSx))+ geom_point()+
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))

lunaize(ggplot(data = MSEfooof %>% filter(Condition == 'eyesOpen') %>% filter(timeScale <=10), aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2))


model <- lmer(MSx ~ Exponent + age + timeScale + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(model)

model <- lmer(MSx ~ Exponent *age + timeScale + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(model)

interactionmodel <- lmer(MSx ~ ageGroup + Exponent*timeScale + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(model)

prediction<- ggpredict(interactionmodel, terms = c("Exponent", "timeScale[2,4,6,8,10,14,16,18,20]", "ageGroup"))

plot(prediction)

ggeffects::hypothesis_test(interactionmodel,  terms = c("Exponent", "timeScale[2,6,10,15,20]", "ageGroup"), test = NULL)


### AUC ----
lunaize(ggplot(data = MSEfooof, aes(y = Exponent, x = Var1)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2) +
          geom_smooth(aes(group = Condition, color = Condition),
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))+xlab("AUC")


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3)+ Condition + Var1, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

lunaize(ggplot(data = MSEfooof, aes(y = Offset, x = Var1))+ geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2)+
          geom_smooth(aes(group = Condition, color = Condition), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))+xlab("AUC")

gam.model <-  mgcv::gamm(Offset ~ s(age, k = 3)+ Condition + Var1, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

### Max Entropy ----

lunaize(ggplot(data = MSEfooof, aes(y = Exponent, x = maxEntropy)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2) +
          geom_smooth(aes(group = Condition, color = Condition),
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3)+ Condition + maxEntropy, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

lunaize(ggplot(data = MSEfooof, aes(y = Offset, x = maxEntropy))+ geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2)+
          geom_smooth(aes(group = Condition, color = Condition), 
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

gam.model <-  mgcv::gamm(Offset ~ s(age, k = 3)+ Condition + maxEntropy, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

# Regions ----
regionalEntropy <- entropyAge %>%
  mutate(region = case_when(
    startsWith(labels, "F") ~ "Frontal",
    startsWith(labels, "C") ~ "Central",
    startsWith(labels, "P") ~ "Parietal",
    startsWith(labels, "O") ~ "Occipital"
  ))


regionalEntropyAvg <- aggregate(cbind(MSx1, MSx2, MSx3, MSx4, MSx5, MSx6, MSx7, MSx8, MSx9, MSx10, MSx11, MSx12, MSx13, MSx14, MSx15, MSx16, MSx17, MSx18, MSx19, MSx20, Var1, age) ~ Subject + visitno + region, data = regionalEntropy, FUN = mean)

regionalEntropyAvg_outlier <- regionalEntropyAvg %>%
  mutate(across(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                  "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                  "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20", "Var1"), naoutlier)) %>% separate(Subject, c('lunaID','vdate'))


regionalEntropyAvg_outlier_long <- regionalEntropyAvg_outlier %>%
  pivot_longer(cols = starts_with(c("MSx1", "MSx2", "MSx3", "MSx4", "MSx5", "MSx6", 
                                    "MSx7", "MSx8", "MSx9", "MSx10", "MSx11", "MSx12", "MSx13", "MSx14", 
                                    "MSx15", "MSx16", "MSx17", "MSx18", "MSx19", "MSx20")),
               names_to = c(".value", "timeScale"),names_pattern = "(\\D+)(\\d+)") %>%  mutate(ageGroup = cut(age, c(0,13,16,19,Inf), labels = c('10-12','13-15','16-18','Adults'))) %>%
  mutate(timeScale = as.numeric(timeScale))

lunaize(ggplot(data = regionalEntropyAvg_outlier_long, aes(x = age, y = MSx, color = as.numeric(timeScale))) + 
          geom_smooth(aes(group = as.numeric(timeScale), color = as.numeric(timeScale)), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) +
  facet_wrap(~region)

regionalEntropyAvg_outlier_long_se <- regionalEntropyAvg_outlier_long %>%
  group_by(timeScale, ageGroup,region) %>%
  summarise(
    MSx_mean = mean(MSx, na.rm = TRUE),
    MSx_se = sd(MSx, na.rm = TRUE) / sqrt(n())
  )

lunaize(ggplot(data = regionalEntropyAvg_outlier_long_se%>% filter(timeScale %in% seq(2, 20, by = 2)), aes(x = as.factor(timeScale), y = MSx_mean))) + geom_point(aes(color=ageGroup, shape = region), position = position_dodge(width = 0.7))+
          geom_errorbar(aes(ymin = MSx_mean - MSx_se, ymax = MSx_mean + MSx_se, color = ageGroup), position = position_dodge(width = 0.7), width = 0.2) + 
  xlab("Time Scale") + 
  ylab("Mean MSx")


## AUC ----
lunaize(ggplot(data = regionalEntropyAvg_outlier_long, aes(x = age, y = Var1, color = region)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, region)), alpha = 0.2) +
          geom_smooth(aes(group = region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("AOC")  


gam.model <-  mgcv::gamm(Var1 ~ s(age, k = 3) + region, data = regionalEntropyAvg_outlier_long, random=list(lunaID=~1))
summary(gam.model$gam)

## Max entropy ----
# Find each subjects max entropy value on the MSE curve 
subs = unique(regionalEntropyAvg_outlier_long$lunaID)
visits = unique(regionalEntropyAvg_outlier_long$visitno)
regions = unique(regionalEntropyAvg_outlier_long$region)



maxValuesRegions <- data.frame()

for (subject in subs) {
  # Check if the subject has the visit number in the dataset
  for (visitNum in visits) {
    for (regionNow in regions) {
    if (any(regionalEntropyAvg_outlier_long$lunaID == subject & regionalEntropyAvg_outlier_long$visitno == visitNum & regionalEntropyAvg_outlier_long$region == regionNow)) {
      subData <- regionalEntropyAvg_outlier_long %>% filter(lunaID == subject) %>% filter(visitno == visitNum)%>% filter(region == regionNow)
      if (all(is.na(subData$MSx))) {
        next
      }
      
      maxSubEntropy <- max(subData$MSx, na.rm = T)
      timescale <- subData$timeScale[which(subData$MSx == maxSubEntropy)]
      
      subInfo <- data.frame(lunaID = subject, visitno = visitNum, maxEntropy = maxSubEntropy, maxTimeScale = timescale, region = regionNow)
      maxValuesRegions <- rbind(maxValuesRegions, subInfo)
    }
    }
  }
}

maxValuesRegions <- merge(maxValuesRegions, ageValues, by = c('lunaID', "visitno"))


lunaize(ggplot(data = maxValuesRegions, aes(x = age, y = maxEntropy, color = region)) + geom_point(aes(color = region)) +
          geom_line(aes(group= interaction(lunaID, region)), alpha = 0.2) +
          geom_smooth(aes(group = region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3) + region, data = maxValuesRegions, random=list(lunaID=~1))
summary(gam.model$gam)


lunaize(ggplot(data = maxValuesRegions, aes(x = age, y = as.numeric(maxTimeScale), color = region)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, region)), alpha = 0.2) +
          geom_smooth(aes(group = region), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2)) + ylab("Max Entropy TimeScale")  


gam.model <-  mgcv::gamm(maxEntropy ~ s(age, k = 3), data = maxValuesAge, random=list(lunaID=~1))
summary(gam.model$gam)


# Regions vs FOOOF ----


MSEfooof <- merge(fooofAvg_outlier, regionalEntropyAvg_outlier_long, by = c('lunaID','vdate')) %>%
  merge(.,  maxValuesAge, by = c("lunaID", "visitno", "age")) %>% mutate(timeScale = (as.numeric(timeScale)))


lunaize(ggplot(data = MSEfooof %>% filter(Condition == 'eyesOpen') %>% filter(timeScale >=14 & timeScale<=20), 
               aes(x = Exponent, y = MSx))+ 
          geom_smooth(aes(group = timeScale, color = timeScale), 
                      method= "loess", alpha=0.2, linewidth=2)) + facet_wrap(~region)


model <- lmer(MSx ~ Exponent + age + timeScale + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(model)

model <- lmer(MSx ~ Exponent *age + timeScale + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(model)

interactionmodel <- lmer(MSx ~ ageGroup + Exponent*timeScale +region + (1|lunaID), data = MSEfooof %>% filter(Condition == 'eyesOpen'))
summary(interactionmodel)

prediction<- ggpredict(interactionmodel, terms = c("Exponent", "timeScale[2,4,6,8,10,14,16,18,20]", "ageGroup", "region"))

plot(prediction)

ggeffects::hypothesis_test(interactionmodel,  terms = c("Exponent", "timeScale[2,6,10,15,20]", "ageGroup"), test = NULL)


### AUC ----
lunaize(ggplot(data = MSEfooof, aes(y = Exponent, x = Var1)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2) +
          geom_smooth(aes(group = Condition, color = Condition),
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))+xlab("AUC")


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3)+ Condition + Var1, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

lunaize(ggplot(data = MSEfooof, aes(y = Offset, x = Var1))+ geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2)+
          geom_smooth(aes(group = Condition, color = Condition), method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))+xlab("AUC")

gam.model <-  mgcv::gamm(Offset ~ s(age, k = 3)+ Condition + Var1, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

### Max Entropy ----

lunaize(ggplot(data = MSEfooof, aes(y = Exponent, x = maxEntropy)) + geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2) +
          geom_smooth(aes(group = Condition, color = Condition),
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))


gam.model <-  mgcv::gamm(Exponent ~ s(age, k = 3)+ Condition + maxEntropy, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

lunaize(ggplot(data = MSEfooof, aes(y = Offset, x = maxEntropy))+ geom_point() +
          geom_line(aes(group= interaction(lunaID, Condition)), alpha = 0.2)+
          geom_smooth(aes(group = Condition, color = Condition), 
                      method=mgcv::"gam", formula = y ~ s(x, k = 3, fx = T), alpha=0.4, linewidth=2))

gam.model <-  mgcv::gamm(Offset ~ s(age, k = 3)+ Condition + maxEntropy, data = MSEfooof, random=list(lunaID=~1))
summary(gam.model$gam)

