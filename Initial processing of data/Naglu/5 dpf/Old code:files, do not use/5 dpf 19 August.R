#libs for data manipulation

library(tidyverse)
library(magrittr)
library(readxl)

# libs for stat analysis

library(lme4)
library(broom)
library(car)
library(performance)
library(emmeans)
library(MASS)


# Libraries for data visualisation

library(scales)
library(ggpubr)
library(ggeasy)
library(ggfortify)
library(ggbeeswarm)
library(ggbeeswarm)
library(kableExtra)

# to make defualt theme of ggplots black and white instead of grey background
theme_set(theme_bw())

# naming zebrafish open field test results as 'data'
#data <- read_excel("Users/Lakshay/Downloads/MPSIII_larvae_behaviourExample/naglu het x het low light 5dpf 3Aug21/Statistics-2021_July_30_naglu_het_x_het_5dpf_light.xlsx") %>%
data <- 
  read_excel("/Users/Lakshay/Desktop/Honours/R data/final 5dpf.xlsx") %>% 

  #renaming
  dplyr::select(Trial = `...2`,
                Position = `...3`,
                Bin = `...4`,
                Distance_travelled = "Distance moved center-point Total mm",
                mean_velocity = "Velocity center-point Mean mm/s",
                time_moving = "Movement Moving / center-point Cumulative Duration s",
                time_at_middle =  "In zone Centre / center-point Cumulative Duration s", 
                zone_transitions = "In zone Centre / center-point Frequency",
                ) %>%
mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")), 
       Trial = str_remove(Trial, pattern = "Trial     ") %>% 
         as.factor()
)
meta <-
  read_excel("/Users/Lakshay/Desktop/Honours/R data/5 dpf/5 dpf 72 samples genotype.xlsx") %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "HET", "HOM")),
  Trial = as.factor(Trial),                
  Position = as.factor(Position)
  )

data <- meta %>%
  left_join(data, by= c("Trial", "Position")) %>%
  na.omit


meta %>%
  na.omit %>%
ggplot(aes(x = Genotype, fill = Genotype)) +
  geom_bar(colour = "black") +
  scale_fill_viridis_d() +
  ggtitle("Number of larval samples with differing genotypes for naglu") +
  ylab("Number of larvae")
#---------------------
 data %>% 
  group_by(fish_id) %>%
mutate(totaldist = sum(Distance_travelled))%>%
dplyr :: select(fish_id, Trial, Genotype, totaldist) %>% 
unique() %>%  
ggplot(aes(x=Genotype, y = totaldist)) +
geom_jitter(aes(colour = Genotype,
                shape = Trial),
            size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (cm)")
  ggsave(filename = "geno.png", width = 13, height = 15, units = "cm", dpi = 400, scale = 0.9)
#ggtitle("Total distance moved by 5 dpf larvae in an hour (cm)")
          
# trial results 
data %>% 
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled))%>%
  dplyr :: select(fish_id, Trial, Genotype, totaldist) %>% 
  unique() %>%  
  ggplot(aes(x= Trial, y = totaldist)) +
  geom_jitter(aes(shape = Genotype), 
              size = 3
              ) +
  geom_boxplot(aes(fill = Trial),
               alpha = 0.25,
               outlier.shape = NA) +
  scale_fill_viridis_d(end = 0.75, option = "inferno") +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (cm)") 
  ggsave(filename = "trial.png", width = 13, height = 15, units = "cm", dpi = 400, scale = 0.8)
  
            
               
              
# plotting as mean data
distancedata <-data %>%
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
  ungroup() %>%
  unique() %>%
  droplevels()

#making histo 
#skewed to the right 
with(distancedata,hist(Distance_travelled))

#linear fixed effect model
model_2 <- data %>%
  group_by(fish_id) %>%
  mutate(Distance_travelled = sum(Distance_travelled)) %>%
dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
  ungroup() %>%
  unique() %>%
  droplevels() %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>%
  lmer(formula = (Distance_travelled)^(1/3) ~ Genotype + Trial + 1|gTrial)  

#plotting the model - doesnt look all that linear, still doing parrametric test though
plot(residuals(model_2) ~ fitted(model_2))
qqPlot(residuals(model_2))

#effect of genotype using chi sq test
Anova(model_2) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type 2 Wald chisquare test of LME") %>%
  kable_styling() 

#time against the wall

#plotting
data %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(1:4, time_at_middle) %>% #  the first 4 cols  are metadata
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
  ggtitle("Total time larvae spent in the middle zone of the well in a one hour trial")







#------------------ poster plot
ggarrange(
  data %>% 
    group_by(fish_id) %>%
    mutate(totaldist = sum(Distance_travelled))%>%
    dplyr :: select(fish_id, Trial, Genotype, totaldist) %>% 
    unique() %>%  
    ggplot(aes(x=Genotype, y = totaldist)) +
    geom_jitter(aes(colour = Genotype,
                    shape = Trial),
                size = 2) +
    geom_boxplot(aes(fill = Genotype),
                 alpha = 0.25,
                 outlier.shape = NA) +
    scale_fill_viridis_d(end = 0.75) +
    scale_colour_viridis_d(end = 0.75) +
    scale_y_continuous(labels = comma, name = "Distance (mm)"),
  #ggtitle("Total distance moved by 5 dpf larvae in an hour (cm)"),

  
data %>% 
  group_by(fish_id) %>%
  mutate(totaldist = sum(Distance_travelled))%>%
  dplyr :: select(fish_id, Trial, Genotype, totaldist) %>% 
  unique() %>%  
  ggplot(aes(x= Trial, y = totaldist)) +
  geom_jitter(aes(shape = Genotype), 
              size = 2
  ) +
  geom_boxplot(aes(fill = Trial),
               alpha = 0.25,
               outlier.shape = NA) +
  scale_fill_viridis_d(end = 0.75, option = 'inferno') +
  scale_colour_viridis_d(end = 0.75, option = 'inferno') +
  scale_y_continuous(labels = comma, name = "Distance (cm)"),


data %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(1:4, time_at_middle) %>% 
  unique() %>% 
  ggplot(aes(x = Genotype, y = time_at_middle)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 2) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA 
  ) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Time spent in middle (s)"),
ncol = 1,
nrow = 3
) +
ggsave(filename = "yesallposter.png", width = 14, height = 30, units = "cm", dpi = 600, scale = 0.7)



