##########Anemonefish prioritise mutualistic partner over social group during marine heatwave – Code for data analysis ##########




##Data Prep ####
data<-read.csv("Article_Database.csv", header =  TRUE)
#load required packages
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

#Ensure that all factors are in correct format 
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
data$Bleached_Overall<-as.factor(data$Bleached_Overall)
data$Behaviour<-as.factor(data$Behaviour)
data$Subject<-as.factor(data$Subject)

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


##Activity plots ####
#Create model to analyse activity levels 
datavis<-lmer(Time_vis_proportion~Heatwave_Status*Bleached_Overall+ThermalExposure_SC+Behaviours_per_second+
                 Subject+N_fish+Anemone_Area_SC+
                 (1|Anemone),
               data=data)
#Model Outputs 
anova(datavis)

summary(datavis)

r2(datavis)
tab_model(datavis) #Easier to read in tab, used throughout code
check_model(datavis) #check model performance

#Posthoc Testing to allow analysis of interactiions 
emmeans_interactionvis <- emmeans(datavis,pairwise ~ Heatwave_Status | Bleached_Overall)
emmeans_interactionvis

emmeans_interactionvis1 <- emmeans(datavis,pairwise ~ Bleached_Overall | Heatwave_Status)
emmeans_interactionvis1
emmip(datavis, Bleached_Overall ~ Heatwave_Status) #shows effect size 

tapply(data$Time_vis_proportion, data$Heatwave_Status, quantile)
tapply(data$Time_vis_proportion, data$Bleached_Overall, quantile)
quantile(data$Time_vis_proportion)


emm_df <- as.data.frame(emmeans_interactionvis)



###Raincloud plots ####
data$predicted <- predict(datavis) #extract predicted values from the model

ggplot(data, aes(x = Heatwave_Status, y = predicted, fill = Bleached_Overall)) +
  
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
    fill = "Bleaching Status",   # legend title for fill
    color = "Bleaching Status"   # legend title for points
  ) +
  
  # Okabe–Ito colours
  scale_fill_manual(
    values = c("FALSE" = "#009E73", "TRUE" = "#D55E00"),
    labels = c("FALSE" = "Unbleached", "TRUE" = "Bleached")
  ) +
  
  scale_color_manual(
    values = c("FALSE" = "#009E73", "TRUE" = "#D55E00"),
    labels = c("FALSE" = "Unbleached", "TRUE" = "Bleached")
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
    # Add axis lines
    axis.line = element_line(color = "black"),
    # Make axis text larger
    axis.text = element_text(size = 14),
    # Make axis titles larger
    axis.title = element_text(size = 16)
  )

##Proportionality Data ####
###Model ####
dataall <- dataall %>%
  group_by(Subject, Heatwave_Status, Bleached_Overall, Anemone) %>% 
  mutate(Total_Count_All = sum(Count)) %>%
  ungroup()

colnames(dataall)

model <- glmmTMB(
  Count ~ Heatwave_Status * Bleached_Overall * Behavioural_Category 
  + N_fish + Anemone_Area_SC+Subject + ThermalExposure_SC + 
   # (1|Anemone)+
   offset(log(Total_Count_All)) ,  
  data = dataall,
  family = poisson
)
#Model Outputs
summary(model)
r2(model)
performance::check_singularity(model)


# Calculate estimated marginal means for the 3-way interaction
emm_model <- emmeans(
  model,
  ~ Heatwave_Status * Bleached_Overall * Behavioural_Category,
  type = "response"   # because Poisson link; gives predicted counts
)

# Create an interaction plot using emmip
emmip(
  model,
  Bleached_Overall ~ Heatwave_Status | Behavioural_Category,
  type = "response",
  CIs = TRUE,
  at = list(N_fish = mean(dataall$N_fish),
            Anemone_Area_SC = mean(dataall$Anemone_Area_SC),
            ThermalExposure_SC = mean(dataall$ThermalExposure_SC))
)

#thermal exposure crashes model, works when scaled
emm <- emmeans(model, ~ Bleached_Overall | Heatwave_Status * Behavioural_Category, type = "response")
pairs(emm)  
check_model(model)
emm1 <- emmeans(model, ~ Heatwave_Status | Bleached_Overall * Behavioural_Category, type = "response")
pairs(emm1) 


emm_df <- as.data.frame(pairs(emm1))

emm_df <- emm_df %>%
  mutate(
    percent_change = (1 / ratio - 1) * 100
  )

emm_df <- emm_df %>%
  mutate(
    direction = case_when(
      percent_change > 0 ~ "Increase",
      percent_change < 0 ~ "Decrease",
      TRUE ~ "No change"
    )
  )

emm_df_clean <- emm_df %>%
  select(Bleached_Overall, Behavioural_Category, contrast,
         ratio, percent_change, direction, p.value)
emm_df_clean
emm_bleach_behav <- emmeans(
  model,
  ~ Bleached_Overall | Behavioural_Category,
  type = "response"
)

pairs(emm_bleach_behav)
tab_model(model)

###Plots ####
#Nowdo summariseandproportions onfiltereddata
behaviour_prop <-data%>%
  group_by(Heatwave_Status,Bleached_Overall,Behavioural_Category) %>%
  summarise(Total_Count= sum(Count), .groups= "drop") %>%
  group_by(Heatwave_Status,Bleached_Overall) %>%
  
  mutate(Proportion= Total_Count / sum(Total_Count)) %>%
  mutate(Group= paste0(Heatwave_Status, "|", ifelse(as.logical(Bleached_Overall), "Bleached","Unbleached")))

#Calculatelabelpositionsandformattedlabels
behaviour_prop <-behaviour_prop %>%
  group_by(Group) %>%
  arrange(Group,Behavioural_Category) %>%
  mutate(
    label_pos= cumsum(Proportion)-(Proportion / 2),
    Label= scales::percent(Proportion, accuracy= 1)
  )
#Plotusingfilteredsummary
#Shows different proportions of behaviours, not used in paper
ggplot(behaviour_prop, aes(x = Group, y = Proportion, fill = Behavioural_Category)) +
  
  geom_bar(stat = "identity") +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0,1)
  ) +
  
  scale_fill_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  labs(
    x = "Time Point | Bleaching Status",
    y = "Proportion of Behaviour Observed",
    fill = "Behaviour"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # removes major gridlines
      panel.grid.minor = element_blank()   # removes minor gridlines
  ) +
  
  geom_text(
    aes(label = scales::percent(Proportion, accuracy = 1)),
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 3
  )

behaviour_prop %>%
  arrange(Heatwave_Status, Bleached_Overall, desc(Proportion)) %>%
  mutate(Proportion = round(Proportion, 3)) %>%
  print(n = Inf)   # show all rows

#Code to create figure showing amount and direction of change
ggplot(behaviour_prop, aes(x = Heatwave_Status, 
                           y = Proportion, 
                           color = Behavioural_Category, 
                           group = Behavioural_Category)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_color_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
  ) +
  
  labs(
    y = "Percentage of behaviours displayed",
    color = "Behaviour category",
  )

## Specific Behaviour Models ####
###Mutualism ####
####Maintenance####
#creating and displaying model
hurdmain  <- glmmTMB(Count ~ Heatwave_Status * Bleached_Overall*Behaviour+ N_fish + Anemone_Area_SC +
                       ThermalExposure_SC+ (1|Anemone)
                     + offset(log(Time_vis)),
                     ziformula = ~Heatwave_Status+ Bleached_Overall+Behaviour,
                     family = ziGamma(link="log"),
                     data = datamain)
summary(hurdmain)
check_model(hurdmain)

plot_model(hurdmain, type = "int", terms = c("Heatwave_Status", "Bleached_Overall" , "Behaviour"))
#posthoc analysis 
emmeans_interactionmain <- emmeans(hurdmain,pairwise ~ Heatwave_Status | Bleached_Overall | Behaviour)
emmeans_interactionmain
plot(emmeans_interactionmain)
emmeans_flipmain <- emmeans(hurdmain,pairwise ~ Bleached_Overall | Behaviour|Heatwave_Status)
emmeans_flipmain

emmip(hurdmain,Bleached_Overall ~ Heatwave_Status *Behaviour)


emm_df1<- as.data.frame(emmeans_interactionmain$emmeans)

emm_df1 <- emm_df1 %>%
  mutate(response = exp(emmean),
         lower = exp(asymp.LCL),
         upper = exp(asymp.UCL))

ggplot(emm_df1, aes(x = Behaviour, y = response, fill = Bleached_Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ Heatwave_Status) +
  labs(y = "Estimated Count (on response scale)", x = "Behaviour", fill = "Bleached") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#palettes for indiviudal plots
green_palette <- c("#009E73", "#33A373", "#66B373", "#99C373", "#CCE373")
vermillion_palette <- c("#FFB380", "#E57A33", "#D55E00", "#A64500", "#7A2A00")
blue_palette <- c("#80B3FF", "#3390E5", "#0072B2", "#005490", "#003D70")
skyblue_palette <- c("#56B4E9", "#7ABDEB", "#9EC5ED", "#C2CDF0", "#E6D5F2")
orange_palette <- c("#E69F00", "#FFB84D", "#FF9933", "#FF8000", "#CC6600")
colnames(emm_df1)
# Line plot
Mainplot<-ggplot(emm_df1, aes(
  x = Heatwave_Status,
  y = response,           # or emmean
  color = Behaviour,
  group = Behaviour,
  fill = Behaviour        # for ribbon
)) +
  # Ribbon for ±SE
  geom_ribbon(aes(ymin = response - SE, ymax = response + SE), alpha = 0.2, color = NA) +
  
  # Line for mean
  geom_line(size = 1) +
  
  # Points for mean
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c("Before" = "Time Point 1",
               "During" = "Time Point 2")
  ) +
  
  scale_color_manual(values = blue_palette) +
  scale_fill_manual(values = blue_palette) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),   # ← draws both x and y axes
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Time Point",
    y = "Estimated count",
    color = "Behaviour",
    fill = "Behaviour"
  )
Mainplot
####Defence####
hurddef <- glmmTMB(Count ~ Heatwave_Status* Behaviour* Bleached_Overall+ N_fish+ Anemone_Area_SC+
                     ThermalExposure_SC+ (1|Anemone)
                   + offset(log(Time_vis)),
                   ziformula = ~Heatwave_Status+ Bleached_Overall+Behaviour,
                   family = ziGamma(link="log"),
                   data = datadef)
summary(hurddef)
check_model(hurddef)

emmeans_interactiondef <- emmeans(hurddef,pairwise ~ Heatwave_Status | Bleached_Overall | Behaviour)
emmeans_interactiondef
plot(emmeans_interactiondef)
emmeans_flipdef <- emmeans(hurddef,pairwise ~ Bleached_Overall | Behaviour|Heatwave_Status)
emmeans_flipdef

emmip(hurddef,Bleached_Overall ~ Heatwave_Status *Behaviour)

emm_df3<- as.data.frame(emmeans_interactiondef$emmeans)

emm_df3 <- emm_df3 %>%
  mutate(response = exp(emmean),
         lower = exp(asymp.LCL),
         upper = exp(asymp.UCL))

ggplot(emm_df3, aes(x = Behaviour, y = response, fill = Bleached_Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ Heatwave_Status) +
  labs(y = "Estimated Count (on response scale)", x = "Behaviour", fill = "Bleached") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#plotting defence behaviours
Defplot<-ggplot(emm_df3, aes(
  x = Heatwave_Status,
  y = response,
  color = Behaviour,
  group = Behaviour,
  fill = Behaviour   # needed for ribbon
)) +
  
  # Ribbon for ±SE
  geom_ribbon(aes(ymin = response - SE, ymax = response + SE), alpha = 0.2, color = NA) +
  
  # Line for mean
  geom_line(size = 1) +
  
  # Points for mean
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c("Before" = "Time Point 1",
               "During" = "Time Point 2")
  ) +
  
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  
  scale_color_manual(values = vermillion_palette) +
  scale_fill_manual(values = vermillion_palette) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),   # ← draws both x and y axes
    strip.text = element_text(face = "bold")
  ) +

  labs(
    x = "Time Point",
    y = "Estimated count",
    color = "Behaviour",
    fill = "Behaviour"
  )
Defplot
###Social####
####Aggressive####
hurdagg <- glmmTMB(Count ~ Heatwave_Status*Bleached_Overall+Behaviour + N_fish+ Anemone_Area_SC+
                      ThermalExposure_SC+(1|Anemone)
                   + offset(log(Time_vis)),
                   ziformula = ~Heatwave_Status+ Bleached_Overall+Behaviour,
                   family = ziGamma(link="log"),
                   data = dataagg)
summary(hurdagg)
check_model(hurdagg)

emmeans_interactionagg <- emmeans(hurdagg,pairwise ~ Heatwave_Status | Bleached_Overall | Behaviour)
emmeans_interactionagg
plot(emmeans_interactionagg)
emmeans_flipagg <- emmeans(hurdagg,pairwise ~ Bleached_Overall | Behaviour|Heatwave_Status)
emmeans_flipagg

emmip(hurdmain,Bleached_Overall ~ Heatwave_Status *Behaviour)

emm_df2<- as.data.frame(emmeans_interactionagg$emmeans)

emm_df2 <- emm_df2 %>%
  mutate(response = exp(emmean),
         lower = exp(asymp.LCL),
         upper = exp(asymp.UCL))

ggplot(emm_df2, aes(x = Behaviour, y = response, fill = Bleached_Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ Heatwave_Status) +
  labs(y = "Estimated Count (on response scale)", x = "Behaviour", fill = "Bleached") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Plotting aggressive behaviour
Aggplot<-ggplot(emm_df2, aes(
  x = Heatwave_Status,
  y = response,
  color = Behaviour,
  group = Behaviour,
  fill = Behaviour   # needed for ribbon
)) +
  
  # Ribbon for ±SE
  geom_ribbon(aes(ymin = response - SE, ymax = response + SE), alpha = 0.2, color = NA) +
  
  # Line for mean
  geom_line(size = 1) +
  
  # Points for mean
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c("Before" = "Time Point 1",
               "During" = "Time Point 2")
  ) +
  
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  
  scale_color_manual(values = green_palette) +
  scale_fill_manual(values = green_palette) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),   # ← draws both x and y axes
    strip.text = element_text(face = "bold")
  ) +

  labs(
    x = "Time Point",
    y = "Estimated count",
    color = "Behaviour",
    fill = "Behaviour"
  )
Aggplot
####Neutral####
hurdneu  <- glmmTMB(Count ~ Heatwave_Status * Bleached_Overall*Behaviour+ N_fish +  Anemone_Area_SC+
                      ThermalExposure_SC+ (1|Anemone)
                    + offset(log(Time_vis)),
                    ziformula = ~Heatwave_Status+ Bleached_Overall+Behaviour,
                    family = ziGamma(link="log"),
                    data = dataneu)
summary(hurdneu)

check_model(hurdneu)

emmeans_interactionneu <- emmeans(hurdneu,pairwise ~ Heatwave_Status | Bleached_Overall | Behaviour)
emmeans_interactionneu
plot(emmeans_interactionneu)
emmeans_flipneu <- emmeans(hurdneu,pairwise ~ Bleached_Overall | Behaviour|Heatwave_Status)
emmeans_flipneu


emmip(hurdneu,Bleached_Overall ~ Heatwave_Status *Behaviour)

emm_df4<- as.data.frame(emmeans_interactionneu$emmeans)

emm_df4 <- emm_df4 %>%
  mutate(response = exp(emmean),
         lower = exp(asymp.LCL),
         upper = exp(asymp.UCL))

ggplot(emm_df4, aes(x = Behaviour, y = response, fill = Bleached_Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ Heatwave_Status) +
  labs(y = "Estimated Count (on response scale)", x = "Behaviour", fill = "Bleached") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



Neuplot<-ggplot(emm_df4, aes(
  x = Heatwave_Status,
  y = response,
  color = Behaviour,
  group = Behaviour,
  fill = Behaviour   # needed for ribbon
)) +
  
  # Ribbon for ±SE
  geom_ribbon(aes(ymin = response - SE, ymax = response + SE), alpha = 0.2, color = NA) +
  
  # Line for mean
  geom_line(size = 1) +
  
  # Points for mean
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    )
  ) +
  
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  
  scale_color_manual(values = skyblue_palette) +
  scale_fill_manual(values = skyblue_palette) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),   # ← draws both x and y axes
    strip.text = element_text(face = "bold")
  ) +

  labs(
    x = "Time Point",
    y = "Estimated count",
    color = "Behaviour",
    fill = "Behaviour"
  )
Neuplot
####Submissive####
hurdsub  <- glmmTMB(Count ~ Heatwave_Status * Bleached_Overall+Behaviour + N_fish + Anemone_Area_SC+
                       ThermalExposure_SC+(1|Anemone)
                    + offset(log(Time_vis)),
                    ziformula = ~Heatwave_Status+ Bleached_Overall+Behaviour,
                    family = ziGamma(link="log"),
                    data = datasub)

summary(hurdsub)


check_model(hurdsub)

emmeans_interactionsub <- emmeans(hurdsub,pairwise ~ Heatwave_Status | Bleached_Overall | Behaviour)
emmeans_interactionsub
plot(emmeans_interactionsub)
emmeans_flipsub <- emmeans(hurdsub,pairwise ~ Bleached_Overall | Behaviour|Heatwave_Status)
emmeans_flipsub


emmip(hurdsub,Bleached_Overall ~ Heatwave_Status *Behaviour)

emm_df5<- as.data.frame(emmeans_interactionsub$emmeans)

emm_df5 <- emm_df5 %>%
  mutate(response = exp(emmean),
         lower = exp(asymp.LCL),
         upper = exp(asymp.UCL))

ggplot(emm_df5, aes(x = Behaviour, y = response, fill = Bleached_Overall)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~ Heatwave_Status) +
  labs(y = "Estimated Count (on response scale)", x = "Behaviour", fill = "Bleached") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot
Subplot<-ggplot(emm_df5, aes(
  x = Heatwave_Status,
  y = response,           # or emmean
  color = Behaviour,
  group = Behaviour,
  fill = Behaviour        # for ribbon
)) +
  
  # Ribbon for ±SE
  geom_ribbon(aes(ymin = response - SE, ymax = response + SE), alpha = 0.2, color = NA) +
  
  # Line for mean
  geom_line(size = 1) +
  
  # Points for mean
  geom_point(size = 2) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c("Before" = "Time Point 1",
               "During" = "Time Point 2")
  ) +
  
  scale_color_manual(values = orange_palette) +
  scale_fill_manual(values = orange_palette) +
  
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),   # ← draws both x and y axes
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = "Time Point",
    y = "Estimated count",
    color = "Behaviour",
    fill = "Behaviour"
  )
Subplot

#Raincloud temperature plots####


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



#to add variance to prop line graph ####

behaviour_prop <- data %>%
  group_by(Heatwave_Status, Bleached_Overall, Behavioural_Category, Anemone) %>%
  summarise(Total_Count = sum(Count), .groups = "drop") %>%
  group_by(Heatwave_Status, Bleached_Overall, Anemone) %>%
  complete(Behavioural_Category, fill = list(Total_Count = 0)) %>%
  mutate(Proportion = Total_Count / sum(Total_Count)) %>%
  mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) %>%  # <- fix NaN
  mutate(Group = paste0(Heatwave_Status, "|",
                        ifelse(as.logical(Bleached_Overall), "Bleached","Unbleached"),
                        "|", Anemone))



ggplot(behaviour_prop, aes(
  x = Heatwave_Status, 
  y = Proportion, 
  color = Behavioural_Category, 
  group = Behavioural_Category,
  fill = Behavioural_Category  # needed for ribbons
)) +
  
  # Ribbon for mean ± SEM
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA) +
  
  # Line connecting mean points
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  
  # Mean points
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    ),
    expand = expansion(add = 0.5)  # increase spacing between labels
  ) +
  
  scale_y_continuous(
    name = "Percentage of behaviours displayed",
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  scale_color_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),   # Draw axes
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  
  labs(
    color = "Behaviour category",
    fill = "Behaviour category"
  )



behaviour_prop <- data %>%
  group_by(Heatwave_Status, Bleached_Overall, Behavioural_Category, Anemone) %>%
  summarise(Total_Count = sum(Count), .groups = "drop") %>%
  group_by(Heatwave_Status, Bleached_Overall, Anemone) %>%
  complete(Behavioural_Category, fill = list(Total_Count = 0)) %>%
  mutate(Proportion = Total_Count / sum(Total_Count)) %>%
  mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) %>%  # fix NaN
  mutate(Group = paste0(Heatwave_Status, "|",
                        ifelse(as.logical(Bleached_Overall), "Bleached","Unbleached"),
                        "|", Anemone))


#Proportion y axis
ggplot(behaviour_prop, aes(
  x = Heatwave_Status, 
  y = Proportion, 
  color = Behavioural_Category, 
  group = Behavioural_Category,
  fill = Behavioural_Category  # needed for ribbons
)) +
  
  # Ribbon for mean ± SEM
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA) +
  
  # Line connecting mean points
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  
  # Mean points
  stat_summary(fun = mean, geom = "point", size = 4) +  # slightly bigger points
  
  facet_wrap(~Bleached_Overall,
             labeller = labeller(
               Bleached_Overall = c(
                 "FALSE" = "Unbleached",
                 "TRUE" = "Bleached"
               )
             )) +
  
  scale_x_discrete(
    name = "Time Point",
    labels = c(
      "Before" = "Time Point 1",
      "During" = "Time Point 2"
    ),
    expand = expansion(add = 0.5)
  ) +
  
  scale_y_continuous(
    name = "Proportion of behaviours displayed"
  ) +
  
  scale_color_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "#009E73",
      "#D55E00",
      "#0072B2",
      "#56B4E9",
      "#E69F00"
    )
  ) +
  
  theme_minimal(base_size = 18) +  # base font size larger
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  
  labs(
    color = "Behaviour category",
    fill = "Behaviour category"
  )

#Combining all single behaiviour plots into one figure ####


(Mainplot | Defplot) /
  (Aggplot | Neuplot) /
  Subplot +
  plot_annotation(
    tag_levels = 'A','B','C','D','E'   # Labels plots A, B, C...
  )

class(Mainplot)
class(Defplot)
class(Aggplot)
class(Neuplot)
class(Subplot)



is.ggplot(Mainplot)
is.ggplot(Defplot)
is.ggplot(Aggplot)
is.ggplot(Neuplot)
is.ggplot(Subplot)




theme_pub <- theme_minimal(base_size = 16) +
  theme(
    # Remove gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Text
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    
    # Facet labels
    strip.text = element_text(size = 14, face = "bold"),
    
    # Legend
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    
    # Spacing
    plot.margin = margin(10, 10, 10, 10)
  )

Mainplot <- Mainplot + theme_pub
Defplot  <- Defplot  + theme_pub
Aggplot  <- Aggplot  + theme_pub
Neuplot  <- Neuplot  + theme_pub
Subplot  <- Subplot  + theme_pub


labs(y = "Estimated count")

x_theme <- list(
  scale_x_discrete(
    labels = c(
      "Before" = "Time\nPoint 1",
      "During" = "Time\nPoint 2"
    )
  ),
  labs(x = NULL),
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 12
    )
  )
)

Mainplot <- Mainplot + x_theme
Defplot  <- Defplot  + x_theme
Aggplot  <- Aggplot  + x_theme
Neuplot  <- Neuplot  + x_theme
Subplot  <- Subplot  + x_theme



combined_plot <- (Mainplot | Defplot|plot_spacer()) / 
  (Aggplot | Neuplot | Subplot) +
  plot_annotation(
    tag_levels = list(c(
      "A: Maintenance",
      "B: Defensive",
      "C: Aggressive",
      "D: Neutral",
      "E: Submissive"
    )))

combined_plot
