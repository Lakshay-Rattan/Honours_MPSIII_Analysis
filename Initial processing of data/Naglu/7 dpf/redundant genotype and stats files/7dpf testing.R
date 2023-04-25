# data manipulation
library(tidyverse)
library(magrittr)
library(readxl)

# statistical analusis
library(lme4)
library(broom)
library(car)
library(performance)
library(emmeans)
library(MASS)

# data visualstion
library(scales)
library(ggpubr)
library(ggeasy)
library(ggfortify)
library(ggbeeswarm)
library(kableExtra)

theme_set(theme_bw())

#deleted the trial 3 data as the trial was made by mistake, and renamed the subsequent trials 3, 4 and 5 instead of being called 4, 5 and 6.
data <-
  read_excel("/Users/Lakshay/Desktop/Honours/R data/7 dpf/Genotype 7dpf only 49-120.xlsx") %>%
  dplyr::select(Trial = 2,
                Position =`...3`,
                Bin = `...4`,
                Distance_travelled = 'Distance moved center-point Total mm',
                mean_velocity = "Velocity center-point Mean mm/s",
                zone_intercepts = "In zone Centre / center-point Frequency",
                time_at_middle = "In zone Centre / center-point Cumulative Duration s") %>%
  #time bins
  mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")),
         Trial = str_remove(Trial, pattern = "Trial     ") %>%
           as.factor()
  )




#---------------------------TURNING THINGS INTO NUMERICS
data$Distance_travelled <- as.numeric(data$Distance_travelled)
data$mean_velocity <- as.numeric(data$mean_velocity)
data$zone_intercepts <- as.numeric(data$zone_intercepts)
data$time_at_middle <- as.numeric(data$time_at_middle)
#------------------------------------------------------------


#------------------ meta
meta <-
  read_excel("/Users/Lakshay/Desktop/Honours/R data/7 dpf/Genotype 7dpf.xlsx") %>%
  mutate(Genotype = factor(Genotype, levels = c("wt", "het", "hom")),
       Trial = as.factor(Trial),
       Position = as.factor(Position)
)

#------- Joining meta to data
data <- meta %>%
  left_join(data, by = c("Trial", "Position")) %>%
  na.omit

#----- Proportions of Genotype
meta %>%
  na.omit %>%
  ggplot(aes(x = Genotype, fill = Genotype)) +
  geom_bar(colour = "black") +
  scale_fill_viridis_d() + # pretty viridis colour scheme
  ggtitle("Number of larvae by genotype")+
  ylab("Number of larvae")

#---------------- TESTING

#-- analysis
# TOTAL DISTANCE TRAVELLED
data %>%
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, totaldist) %>%
  unique() %>%
  ggplot(aes(x = Genotype, y = totaldist)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA
  ) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (mm)") +
  ggtitle("Total distance moved by larvae in one hour(mm)")



data %>%
  group_by(fish_id) %>% # tell R I want to group the 3 time bins per fish
  mutate(totaldist = sum(Distance_travelled)) %>% # then calculate the total distnace per fish over the hour
  dplyr::select(fish_id, Trial, Genotype, totaldist) %>%
  unique() %>%  # as there were 3 replicates because of the bins
  ggplot(aes(x = Genotype, y = totaldist)) +
  stat_summary(fun=mean,
               geom="bar",
               alpha = 0.5,
               aes(group = Genotype,
                   fill = Genotype)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (mm)") +
  ggtitle("Mean Total distance moved by larvae in one hour(mm)")



# make a new object which only contains the sum of the total distance travelled.
distancedata <- data %>%
  group_by(fish_id) %>%
  mutate(Distance_travelled = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
  ungroup() %>%
  unique() %>%
  droplevels()

# Looking at a histogram to visually see the frequency of different distances in bins.
# want something resembling a normal dist so that we can pply a lin model.
#  Looks fairly normal to me, a bit weird at the middle
with(distancedata, hist(Distance_travelled))

# Fit a linear mixed effect model
model_2 <- data %>%
  group_by(fish_id) %>%
  mutate(Distance_travelled = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
  ungroup() %>%
  unique() %>%
  droplevels() %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
  lmer(formula = (Distance_travelled)^(1/3) ~ Genotype + Trial + (1|gTrial))

# check model assumptions
plot(residuals(model_2) ~ fitted(model_2))
qqPlot(residuals(model_2))

# Test the effect of genotype using a type of chi sq test.
Anova(model_2) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()

# Vis the data
data %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(1:4, time_at_middle) %>% # select the first 4 cols which are the metadata
  unique() %>% #as there were 3 replicates because of the bins
  ggplot(aes(x = Genotype, y = time_at_middle)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA # do not show outliers, they are already plotted by geom_point
  ) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Time spent in middle (s)") +
  ggtitle("Total time in centre")

# time in the middle
timedata <- data %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(fish_id, Trial, Genotype, time_at_middle) %>%
  ungroup() %>%
  unique() %>%
  droplevels()

with(timedata, hist(time_at_middle))


model_time <-
  timedata %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
  lmer(formula = (time_at_middle)^(1/3) ~ Genotype + Trial + (1|gTrial))



# check model assumptions
plot(residuals(model_time) ~ fitted(model_time))
qqPlot(residuals(model_time))


# Test the effect of genotype using a type of chi sq test.
Anova(model_time) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()


#--------------------------- lets try trial per time point

timebins <- data %>%
  group_by(fish_id) %>%
  dplyr::select(fish_id, Trial, Genotype, Bin, Distance_travelled)

ggplot(data = timebins, aes(x = Genotype, y = Distance_travelled)) +
  geom_jitter(aes(colour = Genotype,
                shape = Trial),
                size = 2) +

  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA
  ) +
  facet_wrap(~Bin) +

  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (mm)") +
  ggtitle("Total distance moved by larvae in one hour(mm)")


##
#angel

lm(Distance_travelled ~ Genotype + Trial + Genotype:Trial, data = data[data$Bin == "Start-0:20:00",]) %>% anova()
# this says that genotype and trial are significant, but the interaction is not so:
model <- lm(Distance_travelled ~ Genotype + Trial, data = data[data$Bin == "Start-0:20:00",]) %>% aov()

TukeyHSD(model)

data$geno[data$Genotype == "wt"] <- 1
data$geno[data$Genotype == "het"] <- 2
data$geno[data$Genotype == "hom"] <- 3

lm(geno ~ Distance_travelled + mean_velocity + zone_intercepts + time_at_middle, data = data[data$Bin == "0:20:0:40:00",]) %>%
  summary()

# second time bin looks good for distance travelled?

lm(Distance_travelled ~ Genotype + Trial + Genotype:Trial, data = data[data$Bin == "0:20:00-0:40:00",]) %>% anova()
# this says that genotype and trial are significant, but the interaction is not so:
model <- lm(Distance_travelled ~ Genotype + Trial, data = data[data$Bin == "0:20:00-0:40:00",]) %>% aov()

TukeyHSD(model)
#0.051 and 0.65 P val het-wt and hom-wt respectively. 0.99 hom het so wt is clear outlier, sus.
