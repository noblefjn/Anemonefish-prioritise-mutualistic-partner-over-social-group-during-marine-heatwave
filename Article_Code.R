####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
##########Anemonefish prioritise mutualistic partner over social group during marine heatwave#################################
##########Code for data analysis #################################################################

#Data Prep ####
data<-read.csv("Completed_database.2.csv", header =  TRUE)
#load required packages
library(ggsignif)
library(rmarkdown)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(DHARMa)
library(sjPlot)
library(emmeans)
library(performance)
library(ggeffects)
library(lme4)
library(ggsignif)
library(scales)
library(ggeffects)
library(tidyr)
library(gghalves)
library(patchwork)
library(lubridate)
library(ggdist)
library(grid)
library(lmerTest)

#Ensure that all factors are in correct format 
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
data$Bleached_Overall<-as.factor(data$Bleached_Overall)
data$Behaviour<-as.factor(data$Behaviour)
data$Subject<-as.factor(data$Subject)
data$Heatwave_Status<- as.factor(data$Heatwave_Status)

#Create column for Proportion of time visible
data <- data %>%
  mutate(Time_vis_proportion = Time_vis / Total_duration_.s.)
data$Anemone<-as.factor(data$Anemone)

#Filter out anemones where video data is incomplete from main database
data <- data %>% 
  filter(!Anemone %in% c("B62", "B46","B64", "B83", "LY353","LY399","SY517","B142","B28","B49", "B68"))
#remove empty subjects
data <- data %>%
  droplevels(data$Subject)
data %>%
  distinct(Anemone, Bleached_Overall) %>%  # get unique anemone-bleaching status count
  count(Bleached_Overall)

#merges all aggressive categories into one
data <- data %>%
  mutate(Behavioural_Category = case_when(
    Behavioural_Category %in% c("Aggressive - Display", "Aggressive - Overt", "Aggressive Display") ~ "Aggressive",
    TRUE ~ Behavioural_Category
  ))
#do same for neutral
data <- data %>%
  mutate(Behavioural_Category = case_when(
    Behavioural_Category %in% c("Aggressive - Display", "Aggressive - Overt", "Aggressive Display") ~ "Aggressive",
    Behavioural_Category %in% c("Neutral", "Neutral - Bidirectional", "Neutral - Unidirectional") ~ "Neutral",
    TRUE ~ Behavioural_Category
  ))
#scale continuous variables for use in models
data$Day_Temp_SC <- scale(data$Day_Temp, center = TRUE, scale = TRUE)
data$N_fish_SC <- scale(data$N_fish, center = TRUE, scale = TRUE)
data$Anemone_Area_SC <- scale(data$Anemone_Area, center = TRUE, scale = TRUE)
data$ThermalExposure_SC <-scale(data$ThermalExposure, center = TRUE, scale = TRUE)
#Create separate dataframes for each individual behavioural category to allow for easier analysis
data$FishID <- paste(data$Anemone, data$Subject, sep = "_")

data <- data %>%
  mutate(Date = as.Date(Date)) %>%
  group_by(Anemone) %>%
  mutate(
    Days_Before_During = if_else(
      Heatwave_Status == "Before",
      0,
      as.numeric(
        Date[Heatwave_Status == "During"][1] -
          Date[Heatwave_Status == "Before"][1]
      )
    )
  ) %>%
  ungroup()
data$DateDifference<- data$Days_Before_During

data <- data %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(
    Days_since_earliest = as.numeric(Date - min(Date, na.rm = TRUE))
  ) %>%
  ungroup()

data <- data %>%
  group_by(Anemone) %>%
  mutate(
    Day_Temp_Difference = if_else(
      Heatwave_Status == "Before",
      0,
      Day_Temp[Heatwave_Status == "During"][1] -
        Day_Temp[Heatwave_Status == "Before"][1]
    )
  ) %>%
  ungroup()



datamain<- data %>%
  filter(Behavioural_Category %in% c("Maintenance"))
dataagg<- data %>%
  filter(Behavioural_Category %in% c("Aggressive"))
datadef<- data %>%
  filter(Behavioural_Category %in% c("Defence"))
dataneu<- data %>%
  filter(Behavioural_Category %in% c("Neutral"))
datasub<- data%>%
  filter(Behavioural_Category %in% c("Submissive"))
dataall <- data
dataall$FishID <- paste(dataall$Anemone, data$Subject, sep = "_")
data$FishID <- paste(data$Anemone, data$Subject, sep = "_")
head(data)

#Activity model####
##Ready DB####
datavis1 <- data %>%
  group_by(Observation_id, Subject) %>%
  summarise(
    Time_vis_proportion = first(Time_vis_proportion),
    Heatwave_Status = first(Heatwave_Status),
    Bleached_Overall = first(Bleached_Overall),
    Behaviours_per_second = first(Behaviours_per_second),
    N_fish = first(N_fish),
    Anemone_Area_SC = first(Anemone_Area_SC),
    Anemone = first(Anemone),
    FishID = first(FishID),
    DateDifference = first(DateDifference),
    Days_since_earliest = first(Days_since_earliest),
    ThermalExposure = first(ThermalExposure),
    ThermalExposure_SC = first(ThermalExposure_SC), 
    Day_Temp = first(Day_Temp),
    Day_Temp_SC = first (Day_Temp_SC),
    Day_Temp_Difference = first(Day_Temp_Difference),
    Time_vis = first(Time_vis),
    .groups = "drop"
  )

datavis1_rate <- datavis1 %>%
  filter(Time_vis > 0)

# Check
sum(datavis1_rate$Time_vis <= 0, na.rm = TRUE)
nrow(datavis1_rate)

n_distinct(datavis1_rate$FishID)
n_distinct(datavis1_rate$Anemone)
##Model####
datavisTP<-lmer(Time_vis_proportion~Heatwave_Status*Bleached_Overall+Behaviours_per_second+
                  Subject+N_fish+Anemone_Area_SC+Day_Temp_SC+
                  (1|Anemone) +(1|FishID),
                data=datavis1_rate)

summary(datavisTP)
tab_model(datavisTP, digits = 3)
r2(datavisTP)
VarCorr(datavisTP)
performance::check_singularity(datavisTP) 
check_collinearity(datavisTP)
#check to see if VIF caused by interaction
datavisTPnull<-lmer(Time_vis_proportion~Heatwave_Status+Bleached_Overall+Behaviours_per_second+
                      Subject+N_fish+Anemone_Area_SC+Day_Temp_SC+
                      (1|Anemone) +(1|FishID),
                    data=datavis1_rate)
check_collinearity(datavisTPnull) #VIF caused by interaction, ok to go with original 
anova(datavisTP,datavisTPnull)

sim_ratedatavisTP<- simulateResiduals(datavisTP)
testDispersion(sim_ratedatavisTP) #ok
testUniformity(sim_ratedatavisTP) #ok
testZeroInflation(sim_ratedatavisTP) #ok
testOutliers(sim_ratedatavisTP, type = "bootstrap") #ok
plot(sim_ratedatavisTP)#ok
##Emmeans Contrast####
emm_vis <- emmeans(
  datavisTP,
  ~ Heatwave_Status * Bleached_Overall,
  type = "response"
)

emm_vis
emm_time <- emmeans(
  datavisTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

pairs(emm_time, reverse = TRUE, ,
      infer = c(TRUE, TRUE))

emm_int <- emmeans(
  datavisTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_int,
  interaction = "revpairwise",
  type = "response",
  adjust = "none",
  infer = c(TRUE, TRUE)
)
emmeans(
  datavisTP,
  ~ Heatwave_Status * Bleached_Overall,
  type = "response"
)
emm_vis_df <- as.data.frame(emm_vis)
percentage_changevis <- emm_vis_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = emmean[Heatwave_Status == "Before"],
    During = emmean[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )

percentage_changevis
#No sig dif between groups or tp
###Percentage Changes####
emm_bleachvis <- emmeans(
  datavisTP,
  ~ Bleached_Overall | Heatwave_Status
)

pairs(emm_bleachvis, reverse = TRUE,
      infer = c(TRUE, TRUE))
emm_visbleach_df <- as.data.frame(emm_bleachvis)

percentage_changevisbleach <- emm_visbleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = emmean[Bleached_Overall == FALSE],
    Bleached = emmean[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )

percentage_changevisbleach

###Fish Rank Pairwise####
emm_subjectvis <- emmeans(
  datavisTP,
  ~ Subject
)

emm_subjectvis
emm_subjectvis_df <- as.data.frame(emm_subjectvis)

emm_subject <- emmeans(
  datavisTP,
  ~ Subject
)

pairs(emm_subject, reverse = TRUE,
      infer = c(TRUE, TRUE))
percentage_change_subject <- emm_subjectvis_df |>
  dplyr::summarise(
    Rank1 = emmean[Subject == "Rank 1"],
    Rank2 = emmean[Subject == "Rank 2"],
    Rank3 = emmean[Subject == "Rank 3"],
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    Rank2_vs_Rank3 = round(
      ((Rank2 - Rank3) / Rank1) * 100, 1
    )
  )


percentage_change_subject

  #Behaviour Rates####
##Aggression ####
###Data Prep/Check####
dataagg <- dataagg %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()
# Remove observations with zero observation time
dataagg_rate <- dataagg %>%
  filter(Time_vis > 0)

# Check
sum(dataagg_rate$Time_vis <= 0, na.rm = TRUE)
nrow(dataagg_rate)


dataagg_rate_total <- dataagg_rate %>%
  group_by(
    Observation_id,
    Subject,
    Anemone,
    Heatwave_Status,
    Bleached_Overall,
    N_fish,
    Anemone_Area_SC,
    Day_Temp_SC,
    Time_vis,
    FishID
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

###Model ####
model_agg_rateTP <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = dataagg_rate_total,
  family = nbinom1
)
summary(model_agg_rateTP)
r2(model_agg_rateTP)
tab_model(model_agg_rateTP, digits = 3)

sim_agg <- simulateResiduals(model_agg_rateTP)

testZeroInflation(sim_agg)
testDispersion(sim_agg)
testUniformity(sim_agg)
performance::check_singularity(model_agg_rateTP) 
plot(sim_agg)
testQuantiles(sim_agg) 
check_collinearity(model_agg_rateTP) #colinearity low in null 


model_agg_rateTPnull <- glmmTMB(
  Count ~ Heatwave_Status + Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = dataagg_rate_total,
  family = nbinom1
)
anova(model_agg_rateTP,model_agg_rateTPnull) #sig
check_collinearity(model_agg_rateTPnull) #colinearity low in null 

###Emmeans Contrasts#####
emm_agg <- emmeans(
  model_agg_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

emm_agg
contrast(
  emm_agg,
  method = "pairwise",
  adjust = "none"
)
contrast(
  emm_agg,
  method = "revpairwise",
  type = "response",
  adjust = "none",  
  infer = c(TRUE, TRUE)
)

emm_int <- emmeans(
  model_agg_rateTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_int,
  interaction = c("revpairwise", "revpairwise"),
  type = "response",   infer = c(TRUE, TRUE)

)

emm_agg_df <- as.data.frame(emm_agg)
emm_agg_df
percentage_changeagg <- emm_agg_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = response[Heatwave_Status == "Before"],
    During = response[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )
percentage_changeagg

emm_bleachagg <- emmeans(
  model_agg_rateTP,
  ~ Bleached_Overall | Heatwave_Status,
  type = "response"
)

pairs(emm_bleachagg, reverse = TRUE,   infer = c(TRUE, TRUE)
)

emm_aggbleach_df <- as.data.frame(emm_bleachagg)
emm_aggbleach_df

percentage_changeaggbleach <- emm_aggbleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = response[Bleached_Overall == FALSE],
    Bleached = response[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )

percentage_changeaggbleach
###Rank Pairwise####
emm_subjectagg <- emmeans(
  model_agg_rateTP,
  ~ Subject,
  type = "response"
)
pairs(
  emm_subjectagg,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",  
  infer = c(TRUE, TRUE)

)
emm_subjectagg_df <- as.data.frame(emm_subjectagg)
emm_subjectagg_df

percentage_change_subjectagg <- emm_subjectagg_df |>
  dplyr::summarise(
    Rank1 = response[Subject == "Rank 1"],
    Rank2 = response[Subject == "Rank 2"],
    Rank3 = response[Subject == "Rank 3"],
    
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank2 = round(
      ((Rank3 - Rank2) / Rank2) * 100, 1
    )
  )

percentage_change_subjectagg
# Pairwise comparisons between ranks
pairs(
  emm_subjectagg,
  reverse = TRUE,
  type = "response",
  adjust = "tukey"
)
##Neutral####
###Data Prep/Check####
dataneu <- dataneu %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone, Days_since_earliest,Day_Temp_SC) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()
# Remove observations with zero observation time
dataneu_rate <- dataneu %>%
  filter(Time_vis > 0)

# Check
sum(dataneu_rate$Time_vis <= 0, na.rm = TRUE)
nrow(dataneu_rate)

dataneu_rate_total <- dataneu_rate %>%
  group_by(
    Observation_id,
    Subject,
    Anemone,
    Heatwave_Status,
    Bleached_Overall,
    N_fish,
    Anemone_Area_SC,
    Day_Temp_SC,
    Time_vis,
    FishID
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

###Model####

model_neu_rateTP <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = dataneu_rate_total,
  family = nbinom1
)
summary(model_neu_rateTP)
r2(model_neu_rateTP)
performance::r2(
  model_neu_rateTP,
  tolerance = 1e-10
)
performance::check_singularity(model_neu_rateTP) #True = FishID
VarCorr(model_neu_rateTP)
####Sensitivity Check for removing FishID ####
model_neu_rateTPsens <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) +
    offset(log(Time_vis)),
  data = dataneu_rate_total,
  family = nbinom1
)
summary(model_neu_rateTPsens)
anova(model_neu_rateTPsens, model_neu_rateTP)
#exactly the same, ok to leave in Fish ID as random and keep structure identical 


check_collinearity(model_neu_rateTP) #VIF High

model_neu_rateTPnull <- glmmTMB(
  Count ~ Heatwave_Status + Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = dataneu_rate_total,
  family = nbinom1
)


check_collinearity(model_neu_rateTPnull) #moderate VIF Heatwave and Day_Temp 5.45


sim_neu <- simulateResiduals(model_neu_rateTP)

testZeroInflation(sim_neu)#ok
testDispersion(sim_neu)#ok
testUniformity(sim_neu)#ok
testOutliers(sim_neu)
plot(sim_neu)
####Check residuals as top quantile deviates####
plotResiduals(sim_neu)
plotResiduals(sim_neu, form = fitted(model_neu_rateTP))
plotResiduals(sim_neu, form = dataneu_rate_total$Heatwave_Status)
plotResiduals(sim_neu, form = dataneu_rate_total$Bleached_Overall)
plotResiduals(sim_neu, form = dataneu_rate_total$N_fish)
plotResiduals(sim_neu, form = dataneu_rate_total$Day_Temp_SC)

resneu <- residuals(sim_neu)

cor.test(
  resneu,
  dataneu_rate_total$Day_Temp_SC,
  method = "spearman"
)#Residuals not correlated with Day Temp, ok to have some deviation in upper#quantile as all other test statistics are ok

###Emmeans Contracts####


emm_neu <- emmeans(
  model_neu_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

emm_neu
contrast(
  emm_neu,
  method = "pairwise",
  adjust = "none"
)
contrast(
  emm_neu,
  method = "revpairwise",
  type = "response",
  adjust = "none",  
  infer = c(TRUE, TRUE)
)

emm_intneu <- emmeans(
  model_neu_rateTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_intneu,
  interaction = c("revpairwise", "revpairwise"),
  type = "response"
)

emm_neu_df <- as.data.frame(emm_neu)
emm_neu_df
percentage_changeneu <- emm_neu_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = response[Heatwave_Status == "Before"],
    During = response[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )
percentage_changeneu

emm_bleachneu <- emmeans(
  model_neu_rateTP,
  ~ Bleached_Overall | Heatwave_Status,
  type = "response"
)

pairs(emm_bleachneu, reverse = FALSE, ,  
      infer = c(TRUE, TRUE))

emm_neubleach_df <- as.data.frame(emm_bleachneu)
emm_neubleach_df

percentage_changeneubleach <- emm_neubleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = response[Bleached_Overall == FALSE],
    Bleached = response[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )
percentage_changeneubleach

###Rank Pairwise####
emm_subjectneu <- emmeans(
  model_neu_rateTP,
  ~ Subject,
  type = "response"
)
pairs(
  emm_subjectneu,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",
  infer = c(TRUE, TRUE)
)
emm_subjectneu_df <- as.data.frame(emm_subjectneu)
emm_subjectneu_df

percentage_change_subjectneu <- emm_subjectneu_df |>
  dplyr::summarise(
    Rank1 = response[Subject == "Rank 1"],
    Rank2 = response[Subject == "Rank 2"],
    Rank3 = response[Subject == "Rank 3"],
    
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank2 = round(
      ((Rank3 - Rank2) / Rank2) * 100, 1
    )
  )

percentage_change_subjectneu
# Pairwise comparisons between ranks
pairs(
  emm_subjectagg,
  reverse = TRUE,
  type = "response",
  adjust = "tukey"
)
##Submissive####
###Date Prep/Check
datasub <- datasub %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone, Days_since_earliest,Day_Temp_SC) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()
# Remove observations with zero observation time
datasub_rate <- datasub %>%
  filter(Time_vis > 0)

# Check
sum(datasub_rate$Time_vis <= 0, na.rm = TRUE)
nrow(datasub_rate)

datasub_rate_total <- datasub_rate %>%
  group_by(
    Observation_id,
    Subject,
    Anemone,
    Heatwave_Status,
    Bleached_Overall,
    N_fish,
    Anemone_Area_SC,
    Day_Temp_SC,
    Time_vis,
    FishID
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

###Model####
model_sub_rateTP <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datasub_rate_total,
  family = nbinom1
)
summary(model_sub_rateTP)
r2(model_sub_rateTP)
performance::r2(
  model_sub_rateTP,
  tolerance = 1e-10
)
performance::check_singularity(model_sub_rateTP) #True = FishID
check_collinearity(model_sub_rateTP) 

model_sub_rateTPnull <- glmmTMB(
  Count ~ Heatwave_Status + Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datasub_rate_total,
  family = nbinom1
)
check_collinearity(model_sub_rateTPnull) #low vif without interaction

sim_sub <- simulateResiduals(model_sub_rateTP)

testZeroInflation(sim_sub)#ok
testDispersion(sim_sub)#ok
testUniformity(sim_sub)#ok
plot(sim_sub)#ok
testQuantiles(sim_sub)
###Emmeans Contrasts####
emm_sub <- emmeans(
  model_sub_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

emm_sub
contrast(
  emm_sub,
  method = "pairwise",
  adjust = "none"  ,  
  infer = c(TRUE, TRUE)
)
contrast(
  emm_sub,
  method = "revpairwise",
  type = "response",
  adjust = "none",  
  infer = c(TRUE, TRUE),
)

sub_contrasts <- contrast(
  emm_sub,
  method = "pairwise",
  adjust = "none",
  infer = c(TRUE, TRUE)
)

sub_contrasts_df <- as.data.frame(sub_contrasts)

sub_contrasts_df$lower.CL <- round(sub_contrasts_df$lower.CL, 3)
sub_contrasts_df$upper.CL <- round(sub_contrasts_df$upper.CL, 3)

sub_contrasts_df


emm_intsub <- emmeans(
  model_sub_rateTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_intsub,
  interaction = c("revpairwise", "revpairwise"),
  type = "response"
)

emm_sub_df <- as.data.frame(emm_sub)
emm_sub_df
percentage_changesub <- emm_sub_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = response[Heatwave_Status == "Before"],
    During = response[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )
percentage_changesub

emm_bleachsub <- emmeans(
  model_sub_rateTP,
  ~ Bleached_Overall | Heatwave_Status,
  type = "response"
)

pairs(emm_bleachsub, reverse = FALSE,,  
      infer = c(TRUE, TRUE))

emm_subbleach_df <- as.data.frame(emm_bleachsub)
emm_subbleach_df

percentage_changesubbleach <- emm_subbleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = response[Bleached_Overall == FALSE],
    Bleached = response[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )
percentage_changesubbleach

###Rank Pairwise####
emm_subjectsub <- emmeans(
  model_sub_rateTP,
  ~ Subject,
  type = "response"
)
pairs(
  emm_subjectsub,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",  
  infer = c(TRUE, TRUE)
)
emm_subjectsub_df <- as.data.frame(emm_subjectsub)
emm_subjectsub_df

percentage_change_subjectsub <- emm_subjectsub_df |>
  dplyr::summarise(
    Rank1 = response[Subject == "Rank 1"],
    Rank2 = response[Subject == "Rank 2"],
    Rank3 = response[Subject == "Rank 3"],
    
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank2 = round(
      ((Rank3 - Rank2) / Rank2) * 100, 1
    )
  )

percentage_change_subjectsub
# Pairwise comparisons between ranks
pairs(
  emm_subjectsub,
  reverse = TRUE,
  type = "response",
  adjust = "tukey"
)


##Maintenance####
###Data Prep/Check
datamain <- datamain %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone, Days_since_earliest,Day_Temp_SC) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()
# Remove observations with zero observation time
datamain_rate <- datamain %>%
  filter(Time_vis > 0)

# Check
sum(datamain_rate$Time_vis <= 0, na.rm = TRUE)
nrow(datamain_rate)

datamain_rate_total <- datamain_rate %>%
  group_by(
    Observation_id,
    Subject,
    Anemone,
    Heatwave_Status,
    Bleached_Overall,
    N_fish,
    Anemone_Area_SC,
    Day_Temp_SC,
    Time_vis,
    FishID
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

###Model####
model_main_rateTP <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datamain_rate_total,
  family = nbinom1
)
summary(model_main_rateTP)
r2(model_main_rateTP)
performance::r2(
  model_main_rateTP,
  tolerance = 1e-10
)
performance::check_singularity(model_main_rateTP) 

check_collinearity(model_main_rateTP) 
model_main_rateTPnull <- glmmTMB(
  Count ~ Heatwave_Status + Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datamain_rate_total,
  family = nbinom1
)
check_collinearity(model_main_rateTPnull) #low vif with no interaction


sim_main <- simulateResiduals(model_main_rateTP)

testZeroInflation(sim_main)#ok
testDispersion(sim_main)#ok
testUniformity(sim_main)#ok
plot(sim_main)
testQuantiles(sim_main)

####Check residuals as top quantile deviates####
plotResiduals(sim_main)
plotResiduals(sim_main, form = fitted(model_main_rateTP))
plotResiduals(sim_main, form = datamain_rate_total$Heatwave_Status)
plotResiduals(sim_main, form = datamain_rate_total$Bleached_Overall)
plotResiduals(sim_main, form = datamain_rate_total$N_fish)
plotResiduals(sim_main, form = datamain_rate_total$Day_Temp_SC)

resmain <- residuals(sim_main)

cor.test(
  resmain,
  datamain_rate_total$Day_Temp_SC,
  method = "spearman"
)
#Residuals not correlated with Day Temp, ok to have some deviation in upper
#quantile as all other test statistics are ok

###Emmeans Contrasts####
emm_main <- emmeans(
  model_main_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

emm_main
contrast(
  emm_main,
  method = "pairwise",
  adjust = "none",  
  infer = c(TRUE, TRUE)
)
contrast(
  emm_main,
  method = "revpairwise",
  type = "response",
  adjust = "none"
)

emm_intmain <- emmeans(
  model_main_rateTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_intmain,
  interaction = c("revpairwise", "revpairwise"),
  type = "response"
)

emm_main_df <- as.data.frame(emm_main)
emm_main_df
percentage_changemain <- emm_main_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = response[Heatwave_Status == "Before"],
    During = response[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )
percentage_changemain

emm_bleachmain <- emmeans(
  model_main_rateTP,
  ~ Bleached_Overall | Heatwave_Status,
  type = "response"
)

pairs(emm_bleachmain, reverse = FALSE,  
      infer = c(TRUE, TRUE))

emm_mainbleach_df <- as.data.frame(emm_bleachmain)
emm_mainbleach_df

percentage_changemainbleach <- emm_mainbleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = response[Bleached_Overall == FALSE],
    Bleached = response[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )
percentage_changemainbleach

###Rank Pairwise####
emm_subjectmain <- emmeans(
  model_main_rateTP,
  ~ Subject,
  type = "response"
)
pairs(
  emm_subjectmain,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",  
  infer = c(TRUE, TRUE)
)
emm_subjectmain_df <- as.data.frame(emm_subjectmain)
emm_subjectmain_df

percentage_change_subjectmain <- emm_subjectmain_df |>
  dplyr::summarise(
    Rank1 = response[Subject == "Rank 1"],
    Rank2 = response[Subject == "Rank 2"],
    Rank3 = response[Subject == "Rank 3"],
    
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank2 = round(
      ((Rank3 - Rank2) / Rank2) * 100, 1
    )
  )

percentage_change_subjectmain
# Pairwise comparisons between ranks
pairs(
  emm_subjectmain,
  reverse = TRUE,
  type = "response",
  adjust = "tukey"
)


##Defence####
###Data Prep/Check####
datadef <- datadef %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone, Days_since_earliest,Day_Temp_SC) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()
# Remove observations with zero observation time
datadef_rate <- datadef %>%
  filter(Time_vis > 0)

# Check
sum(datadef_rate$Time_vis <= 0, na.rm = TRUE)
nrow(datadef_rate)

datadef_rate_total <- datadef_rate %>%
  group_by(
    Observation_id,
    Subject,
    Anemone,
    Heatwave_Status,
    Bleached_Overall,
    N_fish,
    Anemone_Area_SC,
    Day_Temp_SC,
    Time_vis,
    FishID
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )

###Model####

model_def_rateTP <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datadef_rate_total,
  family = nbinom1
)
summary(model_def_rateTP)
r2(model_def_rateTP)
performance::r2(
  model_def_rateTP,
  tolerance = 1e-10
)
performance::check_singularity(model_def_rateTP) 
check_collinearity(model_def_rateTP) 

model_def_rateTPnull <- glmmTMB(
  Count ~ Heatwave_Status + Bleached_Overall +
    N_fish +
    Anemone_Area_SC +
    Subject +
    Day_Temp_SC +
    (1 | Anemone) + (1|FishID)+
    offset(log(Time_vis)),
  data = datadef_rate_total,
  family = nbinom1
)
check_collinearity(model_def_rateTPnull) #moderate vif = 5.73

sim_def <- simulateResiduals(model_def_rateTP)

testZeroInflation(sim_def)
testDispersion(sim_def)
testUniformity(sim_def)
plot(sim_def)
testOutliers(sim_def)
testQuantiles(sim_def)
####Check residuals as top quantile deviates####


plotResiduals(sim_def)
plotResiduals(sim_def, form = fitted(model_def_rateTP))
plotResiduals(sim_def, form = datadef_rate_total$Heatwave_Status)
plotResiduals(sim_def, form = datadef_rate_total$Bleached_Overall)
plotResiduals(sim_def, form = datadef_rate_total$N_fish)
plotResiduals(sim_def, form = datadef_rate_total$Day_Temp_SC)

resdef <- residuals(sim_def)

cor.test(
  resdef,
  datadef_rate_total$Day_Temp_SC,
  method = "spearman"
)
#Residuals not correlated with Day Temp, ok to have some deviation in upper
#quantile as all other test statistics are ok

###Emmeans Contrasts####
emm_def <- emmeans(
  model_def_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  type = "response"
)

emm_def
contrast(
  emm_def,
  method = "pairwise",
  adjust = "none",  
  infer = c(TRUE, TRUE)
)
contrast(
  emm_def,
  method = "revpairwise",
  type = "response",
  adjust = "none"
)

emm_intdef <- emmeans(
  model_def_rateTP,
  ~ Heatwave_Status * Bleached_Overall
)

contrast(
  emm_intdef,
  interaction = c("revpairwise", "revpairwise"),
  type = "response",  
  infer = c(TRUE, TRUE)
)
emm_def_df <- as.data.frame(emm_def)
percentage_changedef <- emm_def_df |>
  dplyr::group_by(Bleached_Overall) |>
  dplyr::summarise(
    Before = response[Heatwave_Status == "Before"],
    During = response[Heatwave_Status == "During"],
    Percentage_change = round(
      ((During - Before) / Before) * 100, 1
    )
  )
percentage_changedef

emm_bleachdef <- emmeans(
  model_def_rateTP,
  ~ Bleached_Overall | Heatwave_Status,
  type = "response"
)

pairs(emm_bleachdef, reverse = FALSE ,  
      infer = c(TRUE, TRUE))

emm_defbleach_df <- as.data.frame(emm_bleachdef)
emm_defbleach_df

percentage_changedefbleach <- emm_defbleach_df |>
  dplyr::group_by(Heatwave_Status) |>
  dplyr::summarise(
    Unbleached = response[Bleached_Overall == FALSE],
    Bleached = response[Bleached_Overall == TRUE],
    Percentage_change = round(
      ((Unbleached - Bleached) / Bleached) * 100, 1
    )
  )
percentage_changedefbleach

###Rank Pairwise####
emm_subjectdef <- emmeans(
  model_def_rateTP,
  ~ Subject,
  type = "response"
)
pairs(
  emm_subjectdef,
  reverse = TRUE,
  type = "response",
  adjust = "tukey",  
  infer = c(TRUE, TRUE)
)
emm_subjectdef_df <- as.data.frame(emm_subjectdef)
emm_subjectdef_df

percentage_change_subjectdef <- emm_subjectdef_df |>
  dplyr::summarise(
    Rank1 = response[Subject == "Rank 1"],
    Rank2 = response[Subject == "Rank 2"],
    Rank3 = response[Subject == "Rank 3"],
    
    Rank2_vs_Rank1 = round(
      ((Rank2 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank1 = round(
      ((Rank3 - Rank1) / Rank1) * 100, 1
    ),
    
    Rank3_vs_Rank2 = round(
      ((Rank3 - Rank2) / Rank2) * 100, 1
    )
  )

percentage_change_subjectdef
# Pairwise comparisons between ranks
pairs(
  emm_subjectdef,
  reverse = TRUE,
  type = "response",
  adjust = "tukey"
)



#Plots####
##Activity Plots####
# Extract predicted values from the model

emm_time <- emmeans(
  datavisTP,
  ~ Heatwave_Status | Bleached_Overall
)

contrast_time <- pairs(
  emm_time,
  reverse = TRUE
)

contrast_time
datavis1_rate$predicted <- predict(datavisTP)
ggplot(
  datavis1_rate,
  aes(
    x = Heatwave_Status,
    y = predicted,
    fill = Bleached_Overall
  )
) +
  
  geom_half_violin(
    side = "l",
    alpha = 0.6,
    trim = FALSE,
    position = position_dodge(width = 0.8)
  ) +
  
  geom_jitter(
    aes(color = Bleached_Overall),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.8
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  
  labs(
    x = "Time Point",
    y = "Predicted Proportion of Time Visible",
    fill = "Bleaching Status",
    color = "Bleaching Status"
  ) +
  
  scale_fill_manual(
    values = c(
      "FALSE" = "#009E73",
      "TRUE" = "#D55E00"
    ),
    labels = c(
      "FALSE" = "Unbleached",
      "TRUE" = "Bleached"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "#009E73",
      "TRUE" = "#D55E00"
    ),
    labels = c(
      "FALSE" = "Unbleached",
      "TRUE" = "Bleached"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  ) +
  
  # Unbleached: Before vs During
  annotate(
    "segment",
    x = 0.8, xend = 1.8,
    y = 0.85, yend = 0.85,
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 0.8, xend = 0.8,
    y = 0.85, yend = 0.83,
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 1.8, xend = 1.8,
    y = 0.85, yend = 0.83,
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 1,              # moved left
    y = 0.90,
    label = "n.s.",
    size = 5
  ) +
  
  # Bleached: Before vs During
  annotate(
    "segment",
    x = 1.2, xend = 2.2,
    y = 0.92, yend = 0.92,
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 1.2, xend = 1.2,
    y = 0.92, yend = 0.90,
    linewidth = 1.2
  ) +
  annotate(
    "segment",
    x = 2.2, xend = 2.2,
    y = 0.92, yend = 0.90,
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = 1.55,              # moved left
    y = 0.97,
    label = "n.s.",
    size = 5
  )


##Rate Plots ####

green_palette <- c("#009E73", "#33A373", "#66B373", "#99C373", "#CCE373")
vermillion_palette <- c("#FFB380", "#E57A33", "#D55E00", "#A64500", "#7A2A00")
blue_palette <- c("#80B3FF", "#3390E5", "#0072B2", "#005490", "#003D70")
skyblue_palette <- c("#56B4E9", "#7ABDEB", "#9EC5ED", "#C2CDF0", "#E6D5F2")
orange_palette <- c("#E69F00", "#FFB84D", "#FF9933", "#FF8000", "#CC6600")


###Maintenance####
# Obtain estimated rates at 1 second of observation time
emm_main_rate <- emmeans(
  model_main_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  at = list(Time_vis = 1),
  type = "response"
)

# Convert from behaviours/second to behaviours/minute
emm_main_rate_df <- as.data.frame(emm_main_rate) %>%
  mutate(
    response = response * 60,
    asymp.LCL = asymp.LCL * 60,
    asymp.UCL = asymp.UCL * 60
  )
####Brackets####
bracketsmain <- data.frame(
  Bleached_Overall = c(FALSE, TRUE),
  x = c(1, 1),
  xend = c(2, 2),
  y = c(2.6, 2.3),
  yend = c(2.6, 2.3),
  y_lower = c(2.4, 2.1),
  label = c("n.s.", "n.s."),
  label_y = c(2.4, 2.1)
)
Rateplotmain <- ggplot(
  emm_main_rate_df,
  aes(
    x = Heatwave_Status,
    y = response,
    group = Bleached_Overall
  )
) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    alpha = 0.2,
    color = NA,
    fill = green_palette[1]
  ) +
  
  # Model-estimated mean
  geom_line(
    linewidth = 1,
    color = green_palette[1]
  ) +
  
  # Model-estimated mean points
  geom_point(
    size = 2.5,
    color = green_palette[1]
  ) +
  
  # Separate panels for bleaching status
  facet_wrap(
    ~ Bleached_Overall,
    labeller = labeller(
      Bleached_Overall = c(
        "FALSE" = "Unbleached",
        "TRUE" = "Bleached"
      )
    )
  ) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    name = "Estimated behaviours per minute",
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated behaviours per minute"
  )
####Add brackets to plot####
Rateplotmain <- Rateplotmain +
  
  # Horizontal part
  geom_segment(
    data = bracketsmain,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Left vertical
  geom_segment(
    data = bracketsmain,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Right vertical
  geom_segment(
    data = bracketsmain,
    aes(
      x = xend,
      xend = xend,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Significance labels
  geom_text(
    data = bracketsmain,
    aes(
      x = (x + xend) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

Rateplotmain
###Defence####
# Obtain estimated rates at 1 second of observation time
emm_def_rate <- emmeans(
  model_def_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  at = list(Time_vis = 1),
  type = "response"
)

# Convert from behaviours/second to behaviours/minute
emm_def_rate_df <- as.data.frame(emm_def_rate) %>%
  mutate(
    response = response * 60,
    asymp.LCL = asymp.LCL * 60,
    asymp.UCL = asymp.UCL * 60
  )
####Brackets####
bracketsdef <- data.frame(
  Bleached_Overall = c(FALSE, TRUE),
  x = c(1, 1),
  xend = c(2, 2),
  y = c(2.5, 2),
  yend = c(2.5, 2),
  y_lower = c(2.2, 1.7),
  label = c("n.s.", "*"),
  label_y = c(2.2, 1.7)
)
Rateplotdef <- ggplot(
  emm_def_rate_df,
  aes(
    x = Heatwave_Status,
    y = response,
    group = Bleached_Overall
  )
) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    alpha = 0.2,
    color = NA,
    fill = vermillion_palette[1]
  ) +
  
  # Model-estimated mean
  geom_line(
    linewidth = 1,
    color = vermillion_palette[1]
  ) +
  
  # Model-estimated mean points
  geom_point(
    size = 2.5,
    color = vermillion_palette[1]
  ) +
  
  # Separate panels for bleaching status
  facet_wrap(
    ~ Bleached_Overall,
    labeller = labeller(
      Bleached_Overall = c(
        "FALSE" = "Unbleached",
        "TRUE" = "Bleached"
      )
    )
  ) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    name = "Estimated behaviours per minute",
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated behaviours per minute"
  )
####Add brackets to plot####
Rateplotdef <- Rateplotdef +
  
  # Horizontal part
  geom_segment(
    data = bracketsdef,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Left vertical
  geom_segment(
    data = bracketsdef,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Right vertical
  geom_segment(
    data = bracketsdef,
    aes(
      x = xend,
      xend = xend,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Significance labels
  geom_text(
    data = bracketsdef,
    aes(
      x = (x + xend) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

Rateplotdef

###Aggressive####
# Obtain estimated rates at 1 second of observation time
emm_agg_rate <- emmeans(
  model_agg_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  at = list(Time_vis = 1),
  type = "response"
)

# Convert from behaviours/second to behaviours/minute
emm_agg_rate_df <- as.data.frame(emm_agg_rate) %>%
  mutate(
    response = response * 60,
    asymp.LCL = asymp.LCL * 60,
    asymp.UCL = asymp.UCL * 60
  )
####Brackets####
bracketsagg <- data.frame(
  Bleached_Overall = c(FALSE, TRUE),
  x = c(1, 1),
  xend = c(2, 2),
  y = c(0.05, 0.085),
  yend = c(0.05, 0.085),
  y_lower = c(0.045, 0.08),
  label = c("n.s.", "**"),
  label_y = c(0.045, 0.09)
)


Rateplotagg <- ggplot(
  emm_agg_rate_df,
  aes(
    x = Heatwave_Status,
    y = response,
    group = Bleached_Overall
  )
) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    alpha = 0.2,
    color = NA,
    fill = blue_palette[1]
  ) +
  
  # Model-estimated mean
  geom_line(
    linewidth = 1,
    color = blue_palette[1]
  ) +
  
  # Model-estimated mean points
  geom_point(
    size = 2.5,
    color = blue_palette[1]
  ) +
  
  # Separate panels for bleaching status
  facet_wrap(
    ~ Bleached_Overall,
    labeller = labeller(
      Bleached_Overall = c(
        "FALSE" = "Unbleached",
        "TRUE" = "Bleached"
      )
    )
  ) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    name = "Estimated behaviours per minute",
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated behaviours per minute"
  ) 
####Add brackets to plot####
Rateplotagg <- Rateplotagg +
  
  # Horizontal part
  geom_segment(
    data = bracketsagg,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Left vertical
  geom_segment(
    data = bracketsagg,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Right vertical
  geom_segment(
    data = bracketsagg,
    aes(
      x = xend,
      xend = xend,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Significance labels
  geom_text(
    data = bracketsagg,
    aes(
      x = (x + xend) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

Rateplotagg


###Neutral####
# Obtain estimated rates at 1 second of observation time
emm_neu_rate <- emmeans(
  model_neu_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  at = list(Time_vis = 1),
  type = "response"
)

# Convert from behaviours/second to behaviours/minute
emm_neu_rate_df <- as.data.frame(emm_neu_rate) %>%
  mutate(
    response = response * 60,
    asymp.LCL = asymp.LCL * 60,
    asymp.UCL = asymp.UCL * 60
  )
####Brackets####
bracketsneu <- data.frame(
  Bleached_Overall = c(FALSE, TRUE),
  x = c(1, 1),
  xend = c(2, 2),
  y = c(0.55, 0.5),
  yend = c(0.55, 0.5),
  y_lower = c(0.5, 0.45),
  label = c("n.s.", "n.s."),
  label_y = c(0.5, 0.45)
)
Rateplotneu <- ggplot(
  emm_neu_rate_df,
  aes(
    x = Heatwave_Status,
    y = response,
    group = Bleached_Overall
  )
) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    alpha = 0.2,
    color = NA,
    fill = skyblue_palette[1]
  ) +
  
  # Model-estimated mean
  geom_line(
    linewidth = 1,
    color = skyblue_palette[1]
  ) +
  
  # Model-estimated mean points
  geom_point(
    size = 2.5,
    color = skyblue_palette[1]
  ) +
  
  # Separate panels for bleaching status
  facet_wrap(
    ~ Bleached_Overall,
    labeller = labeller(
      Bleached_Overall = c(
        "FALSE" = "Unbleached",
        "TRUE" = "Bleached"
      )
    )
  ) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    name = "Estimated behaviours per minute",
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated behaviours per minute"
  )
####Add brackets to plot####
Rateplotneu <- Rateplotneu +
  
  # Horizontal part
  geom_segment(
    data = bracketsneu,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Left vertical
  geom_segment(
    data = bracketsneu,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Right vertical
  geom_segment(
    data = bracketsneu,
    aes(
      x = xend,
      xend = xend,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Significance labels
  geom_text(
    data = bracketsneu,
    aes(
      x = (x + xend) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

Rateplotneu

###Submissive####
# Obtain estimated rates at 1 second of observation time
emm_sub_rate <- emmeans(
  model_neu_rateTP,
  ~ Heatwave_Status | Bleached_Overall,
  at = list(Time_vis = 1),
  type = "response"
)

# Convert from behaviours/second to behaviours/minute
emm_sub_rate_df <- as.data.frame(emm_sub_rate) %>%
  mutate(
    response = response * 60,
    asymp.LCL = asymp.LCL * 60,
    asymp.UCL = asymp.UCL * 60
  )
####Brackets####
bracketssub<- data.frame(
  Bleached_Overall = c(FALSE, TRUE),
  x = c(1, 1),
  xend = c(2, 2),
  y = c(0.55, 0.4),
  yend = c(0.55, 0.4),
  y_lower = c(0.5, 0.35),
  label = c("n.s.", "n.s."),
  label_y = c(0.5, 0.35)
)
Rateplotsub <- ggplot(
  emm_sub_rate_df,
  aes(
    x = Heatwave_Status,
    y = response,
    group = Bleached_Overall
  )
) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    alpha = 0.2,
    color = NA,
    fill = orange_palette[1]
  ) +
  
  # Model-estimated mean
  geom_line(
    linewidth = 1,
    color = orange_palette[1]
  ) +
  
  # Model-estimated mean points
  geom_point(
    size = 2.5,
    color = orange_palette[1]
  ) +
  
  # Separate panels for bleaching status
  facet_wrap(
    ~ Bleached_Overall,
    labeller = labeller(
      Bleached_Overall = c(
        "FALSE" = "Unbleached",
        "TRUE" = "Bleached"
      )
    )
  ) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    name = "Estimated behaviours per minute",
    labels = scales::label_number(accuracy = 0.01)
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated behaviours per minute"
  )
####Add brackets to plot####
Rateplotsub <- Rateplotsub+
  
  # Horizontal part
  geom_segment(
    data = bracketssub,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Left vertical
  geom_segment(
    data = bracketssub,
    aes(
      x = x,
      xend = x,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Right vertical
  geom_segment(
    data = bracketssub,
    aes(
      x = xend,
      xend = xend,
      y = y,
      yend = y_lower
    ),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  # Significance labels
  geom_text(
    data = bracketssub,
    aes(
      x = (x + xend) / 2,
      y = label_y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5
  )

Rateplotsub

###Combining Rate Plots####
Rateplotagg <- Rateplotagg +
  labs(tag = "A. Aggressive")

Rateplotneu <- Rateplotneu +
  labs(tag = "B. Neutral")

Rateplotsub <- Rateplotsub +
  labs(tag = "C. Submissive")

Rateplotmain <- Rateplotmain +
  labs(tag = "D. Maintenance")

Rateplotdef <- Rateplotdef +
  labs(tag = "E. Defence")




CombinedRatePlot <- 
  Rateplotagg /
  Rateplotneu /
  Rateplotsub /
  Rateplotmain /
  Rateplotdef &
  theme(
    plot.tag = element_text(
      size = 16,
      face = "bold"
    ),
    plot.tag.position = c(0.1, 0.98)
  )

CombinedRatePlot


##Temperature Plots ####
library(sf)
library(ncdf4)
# Read your Kimbe Bay polygon
kimbe_polygon <- st_read("kimbe.kml")

# Make sure it is WGS84
kimbe_polygon <- st_transform(kimbe_polygon, 4326)

# Check it
plot(st_geometry(kimbe_polygon))
st_bbox(kimbe_polygon)

nc <- nc_open("subset1.nc")
print(nc)
nc$dim$time
ncvar_get(nc, "time")
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals

sst <- ncvar_get(nc, "sst")

dates <- as.Date(
  time,
  origin = "1800-01-01"
)

nc_close(nc)

range(dates)
range(lon)
range(lat)
dim(sst)

grid <- expand.grid(
  lon = lon,
  lat = lat
)

grid_sf <- st_as_sf(
  grid,
  coords = c("lon", "lat"),
  crs = 4326
)


inside <- st_within(
  grid_sf,
  kimbe_polygon,
  sparse = FALSE
)

grid_sf$inside <- apply(inside, 1, any)

kimbe_grid <- grid_sf %>%
  filter(inside)

nrow(kimbe_grid)


plot(
  st_geometry(kimbe_polygon),
  col = "lightgrey",
  border = "red"
)

plot(
  st_geometry(kimbe_grid),
  add = TRUE,
  pch = 16,
  col = "blue"
)


kimbe_coords <- st_coordinates(kimbe_grid)

kimbe_grid$lon <- kimbe_coords[, 1]
kimbe_grid$lat <- kimbe_coords[, 2]

sst_kimbe <- lapply(1:nrow(kimbe_grid), function(i) {
  
  lon_i <- which.min(abs(lon - kimbe_grid$lon[i]))
  lat_i <- which.min(abs(lat - kimbe_grid$lat[i]))
  
  data.frame(
    lon = lon[lon_i],
    lat = lat[lat_i],
    date = dates,
    SST = sst[lon_i, lat_i, ]
  )
  
}) %>%
  bind_rows()

sst_weekly <- sst_kimbe %>%
  group_by(date) %>%
  summarise(
    SST = mean(SST, na.rm = TRUE),
    .groups = "drop"
  )
head(sst_weekly)
summary(sst_weekly$SST)


sst_monthly <- sst_weekly %>%
  mutate(
    Year = year(date),
    Month = month(date),
    Month_name = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(Year, Month, Month_name) %>%
  summarise(
    SST = mean(SST, na.rm = TRUE),
    .groups = "drop"
  )


gallDB<-read.csv("RF_data_Nov2022.csv")
gallDB



nearshore_mean_temp <- mean(
  gallDB$MEAN_ANNUAL_TEMP[gallDB$REEFTYPE == "Nearshore"],
  na.rm = TRUE
)

annualplot <- ggplot(
  sst_monthly,
  aes(
    x = Month,
    y = SST,
    group = Year,
    colour = factor(Year)
  )
) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 2.5) +
  
  # Mean annual temperature of nearshore reefs
  geom_hline(
    yintercept = nearshore_mean_temp,
    linetype = "dashed",
    linewidth = 1,
    colour = "black"
  ) +
  
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  
  scale_y_continuous(
    name = "Mean SST (°C)"
  ) +
  
  labs(
    x = "Month",
    colour = "Year"
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

annualplot
annualplot <- annualplot +
  annotate(
    "text",
    x = 7,
    y = nearshore_mean_temp,
    label = "Galbraith et al 2022 annual inshore mean temperature",
    hjust = 1,
    vjust = -0.5,
    size = 5
  )
annualplot
names(nc$var)
unique(gallDB$REEFTYPE)

#SST anomaly plot

nc3 <- nc_open("subset3.nc")

print(nc3)

names(nc3$var)


names(nc3$var)
nc3$dim$time$units
range(nc3$dim$time$vals)

lon <- nc3$dim$lon$vals
lat <- nc3$dim$lat$vals
time <- nc3$dim$time$vals

anom <- ncvar_get(nc3, "anom")

dates <- as.Date(
  time,
  origin = "1800-01-01"
)

nc_close(nc3)

range(dates)

grid <- expand.grid(
  lon = lon,
  lat = lat
)

grid_sf <- st_as_sf(
  grid,
  coords = c("lon", "lat"),
  crs = 4326
)

inside <- st_within(
  grid_sf,
  kimbe_polygon,
  sparse = FALSE
)

grid_sf$inside <- apply(inside, 1, any)

kimbe_grid <- grid_sf %>%
  filter(inside)
nrow(kimbe_grid)

kimbe_coords <- st_coordinates(kimbe_grid)

kimbe_grid$lon <- kimbe_coords[, 1]
kimbe_grid$lat <- kimbe_coords[, 2]

anom_kimbe <- lapply(1:nrow(kimbe_grid), function(i) {
  
  lon_i <- which.min(abs(lon - kimbe_grid$lon[i]))
  lat_i <- which.min(abs(lat - kimbe_grid$lat[i]))
  
  data.frame(
    lon = lon[lon_i],
    lat = lat[lat_i],
    date = dates,
    SST_anomaly = anom[lon_i, lat_i, ]
  )
  
}) %>%
  bind_rows()

anom_daily <- anom_kimbe %>%
  group_by(date) %>%
  summarise(
    SST_anomaly = mean(SST_anomaly, na.rm = TRUE),
    .groups = "drop"
  )

head(anom_daily)
summary(anom_daily$SST_anomaly)

anom_monthly <- anom_daily %>%
  mutate(
    Year = year(date),
    Month = month(date),
    Month_name = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(Year, Month, Month_name) %>%
  summarise(
    SST_anomaly = mean(SST_anomaly, na.rm = TRUE),
    .groups = "drop"
  )

anomplot<-ggplot(
  anom_monthly,
  aes(
    x = Month,
    y = SST_anomaly,
    group = Year,
    colour = factor(Year)
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 2.5) +
  
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  
  scale_y_continuous(
    name = "SST anomaly (°C)"
  ) +
  
  labs(
    x = "Month",
    colour = "Year"
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )


CombTempBayPlot <- 
  (annualplot | anomplot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

CombTempBayPlot

###Time Point Temperature Comparison ####
ggplot(data, aes(x = Heatwave_Status, y = Day_Temp, fill = Heatwave_Status)) +
  
  # Half violin (cloud)
  geom_half_violin(
    side = "l",
    alpha = 0.6,
    trim = FALSE,
    position = position_dodge(width = 0.8)
  ) +
  
  # Jittered points (rain)
  geom_jitter(
    aes(color = Heatwave_Status),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.8
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  
  # Boxplot (handle)
  geom_boxplot(
    width = 0.15,
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  
  # Labels
  labs(
    x = "Time Point",
    y = "Temperature (°C)"
  ) +
  
  # Match Plot 1 x-axis relabeling
  scale_x_discrete(
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  # Match Plot 1 colour scheme (Okabe–Ito)
  scale_fill_manual(
    values = c(
      "Before" = "#E69F00",   # cool blue
      "During" = "#56B4E9"    # warm orange
    ),
    name = NULL
  ) +
  
  scale_color_manual(
    values = c(
      "Before" = "#E69F00",
      "During" = "#56B4E9"
    ),
    name = NULL
  )+
  
  # Theme (same as Plot 1)
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    # Add axis lines
    axis.line = element_line(color = "black"),
    
    # Make axis text larger
    axis.text = element_text(size = 14),
    
    # Make axis titles larger
    axis.title = element_text(size = 16)
  )


####Bleached Numbers per time point ####


data_unique <- data %>%
  distinct(Anemone, Heatwave_Status, Bleached_Status)
unique(data_unique$Bleached_Status)
unique(data_unique$Heatwave_Status)

data_unique$Bleached_Status <- factor(
  data_unique$Bleached_Status,
  levels = c("Non-bleached", "Bleached"),  # keep order
  labels = c("Unbleached", "Bleached")     # new labels for legend
)

# Stacked bar chart
ggplot(data_unique, aes(x = Heatwave_Status, fill = Bleached_Status)) +
  
  geom_bar(
    position = "stack",
    width = 0.7
  ) +
  
  labs(
    x = "Time Point",
    y = "Number of Anemones",
    fill = "Bleaching Status"
  ) +
  
  scale_x_discrete(
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Unbleached" = "#009E73",
      "Bleached"   = "#D55E00"
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


data_unique %>%
  count(Heatwave_Status, Bleached_Status)




