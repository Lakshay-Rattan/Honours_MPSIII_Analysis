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

#deleted the trial 3 data on excel file as the trial was made by mistake, and renamed the subsequent trials 3, 4 and 5 instead of being called 4, 5 and 6.
data.naglu.7 <-
  read_excel("/Users/Lakshay/Desktop/Honours/R data/7 dpf/Statistics.xlsx") %>%
  dplyr::select(Trial = 2,
                Position =`...3`,
                Bin = `...4`,
                Distance_travelled = 'Distance moved center-point Total mm',
                mean_velocity = "Velocity center-point Mean mm/s",
                zone_intercepts = "In zone centre / center-point Frequency",
                time_at_middle = "In zone centre / center-point Cumulative Duration s",
                -`Independent Variable Independent Variable Treatment`) %>%
  #time bins
  mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")),
         Trial = str_remove(Trial, pattern = "Trial     ") %>%
           as.factor()
  )

data.20.naglu.7 %>%
mutate(totaldist = sum(Distance_travelled)) %>%
ggplot() +
  geom_density(aes(x = totaldist, fill = Genotype),
               alpha = 0.5)


  data.naglu.7 <-
    data.naglu.7 %>%
    mutate(Distance_travelled = as.numeric(Distance_travelled))

# read in metadata
meta.naglu.7 <-
  read_excel("/Users/Lakshay/Desktop/Honours/R data/7 dpf/Old geno and stats stuff/Genotype 7dpf different name of arenas.xlsx") %>%
  mutate(Genotype = factor(Genotype, levels = c("wt", "het", "hom")),
         Trial = as.factor(Trial),
         Position = as.factor(Position)
  )

# join the spreadsheets together
data.naglu.7 <-
  meta.naglu.7 %>%
  left_join(data.naglu.7, by = c("Trial", "Position")) %>%
  na.omit

#----- Proportions of Genotype
meta.naglu.7 %>%
  na.omit %>%
  ggplot(aes(x = Genotype, fill = Genotype)) +
  geom_bar(colour = "black") +
  scale_fill_viridis_d() + # pretty viridis colour scheme
  ggtitle("Number of larvae by genotype") +
  ylab("Number of larvae")

 ggsave(filename = "proportions 7dpf.png", width = 20, height = 15, units = "cm", dpi = 600, scale = 1)


#-- analysis
# TOTAL DISTANCE TRAVELLED
data.naglu.7 %>%
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, totaldist) %>%
  unique() %>%
  ggplot(aes(x = Genotype, y = totaldist)) +
  geom_jitter(aes(colour = genotype,
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

#--- distance travelled BETWEEN TESTS
data.naglu.7 %>%
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
  ggtitle("Total distance moved by larvae in one hour(mm)") +
  facet_wrap(~Trial, nrow = 1)

data.naglu.7 %>%
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, totaldist) %>%
  subset(Trial == c(1,2)) %>%
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
  ggtitle("Total distance moved by larvae in one hour(mm)") +
  facet_wrap(~Trial)



