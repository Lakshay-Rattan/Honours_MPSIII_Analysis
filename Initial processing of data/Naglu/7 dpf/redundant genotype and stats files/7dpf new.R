# file.rename(from ='Untitled2', to = 'sussy'). This script will aim to separate activity by time bins

# data manipulation
library(tidyverse)
library(magrittr)
library(readxl)
library(dplyr)

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



data.naglu.7 <-
  read_excel('final data after fixing sensitivity issue.xlsx')
  dplyr::rename(
    "Trial" = 2,
    "Position" = 3,
    Bin = 4,
    Distance_travelled = 'Distance moved center-point Total mm',
    mean_velocity = 'Velocity center-point Mean mm/s',
    zone_intercepts = "In zone centre / center-point Frequency",
    time_at_middle = "In zone centre / center-point Cumulative Duration s") %>%


   #time bins
  mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")),
         Trial = str_remove(Trial, pattern = "Trial     ") %>%
          as.factor())



#mutate turns values into a variable (adding a bunch of data into a single value), factor turns the variable into a categorical variable (like wt, hom, het),
#level turn categories into ordered categories (Jan, Feb, March), c(combine) makes all of the three factors get sussed.
#without c, it would go up to wt, the program will go back to the start of the function, resulting in an error.
#have to include all factors for leveling to work, self explanatory.


meta.naglu.7 <- read_excel("/Users/Lakshay/Desktop/Honours/R data/7 dpf/Genotype 7dpf only 49-120.xlsx") %>%
  mutate(Genotype = factor(Genotype, levels = c("wt","het","hom")),
#apparently as.factor is better than factor() for these, maybe because of the levels in factor()?
         Trial = as.factor(Trial),
         Position = as.factor(Position)
  )
#Replacing names of things in the table
#meta <- meta %>%
  #mutate(Position = str_replace(Position, "Arena ", "A"))



#at this point, the name of the trials are different between the data and meta datasets. Lachie compared using data$Trial and meta$Trial.
#Then, he removed the word "Trial     " so now you can left join.
fulldata <- data_new %>%
  left_join(meta)

fulldata <- fulldata[-5]


fulldata <- fulldata %>% na.omit

#The daniovision calibration length was set at 1000 mm instead of 120. correcting it by multiplying by 0.012
fulldata$distance_travelled <- fulldata$distance_travelled * 0.012


# genotype count
ggplot(data = fulldata) +
  geom_bar(mapping = aes(x = Genotype))

# Distance travelled vs genotype
#have to group by fish_id as there are multiple time boxes per sample.


# grouping the time intervals by fish ID, then the mutate function adds a column. as you've grouped by fish id,
#the sum will be the distance travelled per fish id. IMPORTANT: NEED TO DO <- to save your progress otherwise it will
# only print it at the bottom like a bitch rather than saving that shit. Alternatively just do fulldata %<>% to show off.
#use unique to get rid of the duplicate  fish id (1,1,1 to 1).

distancedata <- fulldata %>%
  group_by(fish_id) %>%
  mutate(totaldist = sum(distance_travelled)) %>%
  dplyr::select(Genotype, fish_id, Trial, totaldist) %>%

#ask Karissa why I need to ungroup, apparently it can fuck shit up later on in the code
  ungroup() %>%
  unique() %>%
  droplevels()

# Look at the histogram.
# want something resembling a normal dist so that we can pply a lin model.
#  Looks fairly normal to me, a bit weird at the middle
with(distancedata, hist(distance_travelled))

#total distance by genotype in an hour
distancedata %>%
ggplot(aes(x = Genotype, y = totaldist)) +
  ylim(c(0,1300)) +
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


#mean distance by genotype in an hour, ask karissa about stats summary
  fulldata %>%
    group_by(fish_id) %>%
    mutate(totaldist = sum(distance_travelled)) %>%
    dplyr::select(fish_id, Trial, Genotype, totaldist) %>%
    unique() %>%
    ggplot(aes(x = Genotype, y = totaldist)) +
    ylim(c(0,1300)) +
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
  ggtitle("Mean total distance moved by larvae in one hour (mm)")


  # checking if data is hella symetrical, outliers giving the thing a tail.
  with(distancedata, hist(totaldist))

  model_2 <- fulldata %>%
    group_by(fish_id) %>%
    mutate(Distance_travelled = sum(distance_travelled)) %>%
    dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
    ungroup() %>%
    unique() %>%
    droplevels() %>%
    mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
    lmer(formula = (Distance_travelled)^(1/3) ~ Genotype + Trial + (1|gTrial))

  # check model assumptions
  plot(residuals(model_2) ~ fitted(model_2))
  qqPlot(residuals(model_2)) #data is too peaked in the middle

  # Test the effect of genotype using a type of chi sq test.
  Anova(model_2) %>%
    as.data.frame() %>%
    dplyr::rename(pval = `Pr(>Chisq)`) %>%
    kable(caption = "Type II Wald chisquare test of LME") %>%
    kable_styling()
#####no effect of genotype or trial#######, probably as all trials were done in afternoon

  #lets suss out the time in middle
middle <-
  fulldata %>%
    group_by(fish_id) %>%
    mutate(time_at_middle = sum(time_at_middle)) %>%
    dplyr::select(fish_id, Trial, Bin, Genotype, time_at_middle) %>%
    unique()

fiddle <- middle %>% unique()

middle %>%
  unique() %>%
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

#data analysis with middle/fiddle/diddle data

with(middle, hist(time_at_middle))

model_time <- middle %>%
  group_by(fish_id) %>%
  dplyr::select(fish_id, Trial, Genotype, time_at_middle) %>%
  ungroup() %>%
  unique() %>%
  droplevels() %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
  lmer(formula = (time_at_middle)^(1/3) ~ Genotype + Trial + (1|gTrial))

# check model assumptions
plot(residuals(model_time) ~ fitted(model_time))
qqPlot(residuals(model_time)) #data is too peaked in the middle

# Test the effect of genotype using a type of chi sq test.
Anova(model_time) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()




