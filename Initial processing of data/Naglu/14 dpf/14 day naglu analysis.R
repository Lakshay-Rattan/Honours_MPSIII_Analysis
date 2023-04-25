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
library(dlookr)

# data visualstion
library(scales)
library(ggpubr)
library(ggeasy)
library(ggfortify)
library(ggbeeswarm)
library(kableExtra)

theme_set(theme_bw())


Total_dist %>%
  group_by(Genotype) %>%
#  dplyr::mutate(Total_dist = as.numeric(Total_dist)) %>%
  summarise(
    n = n(),
    mean = mean(Total_distance),
    median = median(Total_distance))



data <-
  read_excel("Statistics-Retesting 14dpf naglu 20 min bins.xlsx")
data <-
  subset(data, select = -c(`...1`,`Independent Variable Independent Variable Treatment`))

data <-
  data %>%
  rename(

    "Trial" = `...2`,
          "Position" = `...3`,
          "Bin" = `...4`,
          "min dist" = "Distance moved center-point Minimum mm",
          "Distance_travelled" = "Distance moved center-point Total mm",
          "Mean_Velocity" = "Velocity center-point Mean mm/s",
          "Freq_in_centre" = "In zone Centre / center-point Frequency" ,
          "Time_in_centre" = "In zone Centre / center-point Cumulative Duration s"
    ) %>%
  mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")),
       Trial = str_remove(Trial, pattern = "Trial     ") %>%
         as.factor())


meta <- read_excel("14dpf naglu genotype.xlsx") %>%
  subset(select = -`...5`) %>%
  mutate(Genotype = factor(Genotype, levels = c("wt","het","hom")),
         Trial = as.factor(Trial),
         Position = as.factor(Position)
  )

fulldata <-
  meta %>%
  left_join(data, by = c("Trial", "Position")) %>%
  na.omit

#subsetting T4 and 5 as they have compounding variable (arised from different pairmates)
  fulldata <-
fulldata %>%
  subset(!(Trial==4) &
         !(Trial==5))



# porportions
Proportion <-
ggplot(data = fulldata) +
        geom_bar(mapping = aes(x = Genotype))

#Karissa's function
plotGenoProps <- function(input, xaxis) {
  ggplot(data = input) +
    geom_bar(aes(x = Genotype))
}

plotGenoProps(input = fulldata, xaxis = Genotype)



#distance travelled----------------------------------------------------------------------------------------------------------------------

#this is the distance travelled in the hour
Total_dist <-
  fulldata %>%
  group_by(fish_id) %>%
  mutate(Total_distance = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Position, Genotype, Total_distance, Freq_in_centre) %>%
  dplyr::distinct(Total_distance, .keep_all = TRUE) %>%
  ungroup() %>%
  droplevels()



Total_dist %>% saveRDS("B.14.totaldist.rds")

fulldata %>% saveRDS("B.14.fulldata.rds")


#all things after this point are redundant as I repeat them in the final MPSI_III_Analysis project.




# Fit a linear mixed effect model
model_distance <- Total_dist %>%
  group_by(fish_id) %>%
  ungroup() %>%
  unique() %>%
  droplevels() %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>%
  lmer(formula = (log(Total_distance) ~ Genotype + (1|gTrial)))



# check model assumptionq
plot(residuals(model_distance) ~ fitted(model_distance)) #visual version of levene test.
qqPlot(residuals(model_distance)) #visually used to find normality of data

# Test the effect of genotype using a type of chi square test, major effect of trial.
Anova(model_distance) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()



#emeans total distance travelled----
print(emmeans(model_distance, ~ Genotype), type ="emmean") %>%
  as_tibble() %>%
  ggplot(aes(x = Genotype, y = emmean, colour = Genotype)) +
  geom_col(aes(fill = Genotype),
           alpha = 0.5,
           width = 0.75,
           position = position_dodge()) +
  geom_errorbar(aes(ymin = lower.CL, ymax =upper.CL),
                width = 0.125,
                size = 1,
                position = position_dodge(width = 0.75)) +
  scale_color_viridis_d(end = 0.75, option = "plasma") +
  scale_fill_viridis_d(end = 0.75, option = "plasma") +
  labs(y = "model-predicted distance travelled",
       x = "Genotype") +
  theme(text = element_text(size = 20))


#raw total distance travelled----

ggplot(data = Total_dist, aes(x = Genotype, y = Total_distance)) +
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
  ggtitle("Total distance moved by larvae in one hour (mm)") +
  theme(text = element_text(size = 20))


  ggsave(filename = "14dpf total distance travelled.png", width = 20, height = 15, units = "cm", dpi = 300, scale = 1.2)





  lm_bin <- fulldata %>%
    mutate(gTrial = interaction(Genotype, Trial, drop = TRUE),
           fish_id = as.character(fish_id)) %>% # these are to define interaction effects in the model
    lmer(formula = sqrt(Distance_travelled) ~ (Genotype*Bin) + Trial + (1|gTrial) + (1|fish_id),
         data = .)

  plot(residuals(lm_bin) ~ fitted(lm_bin))
  qqPlot(residuals(lm_bin))



  Anova(lm_bin) %>%
    as.data.frame() %>%
    dplyr::rename(pval = `Pr(>Chisq)`) %>%
    kable(caption = "Type II Wald chisquare test of LME") %>%
    kable_styling()


  #distance travelled by bin raw and emmean----

  print(emmeans(lm_bin, ~ Genotype * Bin), type = "response") %>%
    as_tibble() %>%
    mutate(binforvis = case_when(
      Bin == 1 ~ "Start-0:20:00",
      Bin == 2 ~ "0:20:00-0:40:00",
      Bin == 3 ~ "0:40:00-1:00:00",
    )) %>%
    ggplot(aes(x = binforvis, y = response, colour = Genotype)) +
    geom_col(aes(fill = Genotype),
             alpha = 0.5,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  size = 1,
                  position = position_dodge()) +
    theme(axis.text.x = element_text(hjust = 1,
                                     vjust = 1,
                                     angle = 45),
          legend.position = "bottom") +
    scale_color_viridis_d(end = 0.75) +
    scale_fill_viridis_d(end = 0.75) +
    facet_wrap(~Bin)



  #BINNED DISTANCE DATA

  fulldata %>%
  ggplot(aes(x = Genotype, y = Distance_travelled)) +
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
    ggtitle("Total distance moved by larvae in one hour (mm)") +
    theme(text = element_text(size = 15)) +
    facet_wrap(~Bin)

  ggsave(filename = "14dpf  distance travelled over an hour per 20 minute time bin.png", width = 30, height = 10, units = "cm", dpi = 600, scale = 1.2)






  # DISTANCE DATA BY TRIAL

  distancedata %>%
    subset(Distance_travelled) %>%
    ggplot(aes(x = Trial, y = Distance_travelled)) +
    geom_jitter(aes(
                    shape = Genotype),
                size = 4) +
    geom_boxplot(aes(fill = Trial),
                 alpha = 0.25,
                 outlier.shape = NA
    ) +
    scale_fill_viridis_d(option = "B", end = 0.75) +
    scale_colour_viridis_d(option = "B", end = 0.75) +
    scale_y_continuous(labels = comma, name = "Distance (mm)") +
    facet_wrap(~Trial, scales = "free_x", nrow = 1) +
    ggtitle("Total distance moved by larvae in one hour in each trial (mm)") +
    theme(text = element_text(size = 15))



  ggsave(filename = "14dpf total distance by trial.png", width = 20, height = 10, units = "cm", dpi = 300, scale = 1.5)



#Time in centre ------------------------------------------------------------------------------------------------------------------
#time in centre graph

centre_time_freq_graph <-
  fulldata %>%
  group_by(fish_id)%>%
  dplyr::select(Trial, Position, Genotype, fish_id, Bin, Time_in_centre,Freq_in_centre) %>%
  ungroup() %>%
  droplevels()


  centre_time_freq_graph %>%
  ggplot(aes(x = Genotype, y = Time_in_centre)) +
    geom_jitter(aes(colour = Genotype,
                    shape = Trial),
                size = 4) +
    geom_boxplot(aes(fill = Genotype),
                 alpha = 0.25,
                 outlier.shape = NA
    ) +
    scale_fill_viridis_d(end = 0.75) +
    scale_colour_viridis_d(end = 0.75) +
    scale_y_continuous(labels = comma, name = "Time in centre (s)") +
    ggtitle("Total time (in seconds) spent in the central zone for 14dpf larvae in one hour") +
    facet_wrap (~Bin)

  ggsave(filename = "14dpf total time spent in the centre per 20 minute time bin.png", width = 20, height = 10, units = "cm", dpi = 600, scale = 1.2)



  #this is the distance travelled in the hour
  Total_time <-
    fulldata %>%
    group_by(fish_id) %>%
    mutate(Total_time = sum(Time_in_centre)) %>%
    dplyr::select(Trial, Position, Genotype, Total_time, fish_id) %>%
    # another way but a bit more dodgy
    dplyr::distinct(Total_time, .keep_all = TRUE) %>%
    ungroup() %>%
    droplevels()




  # Fit a linear mixed effect model for SQRT DATA
  model_time <- Total_time %>%
    group_by(fish_id) %>%
    ungroup() %>%
    unique() %>%
    droplevels() %>%
    mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>%
    lmer(formula = (sqrt(Total_time) ~ Genotype + Trial + (1|gTrial)))



  # check model assumptionq
  plot(residuals(model_time) ~ fitted(model_time)) #visual version of levene test.
  qqPlot(residuals(model_time)) #visually used to find normality of data

  # Test the effect of genotype using a type of chi square test, major effect of trial.
  Anova(model_time) %>%
    as.data.frame() %>%
    dplyr::rename(pval = `Pr(>Chisq)`) %>%
    kable(caption = "Type II Wald chisquare test of LME") %>%
    kable_styling()



  lm_time <- hgsnat.5 %>%
    mutate(gTrial = interaction(Genotype, Trial, drop = TRUE),
           fish_id = as.character(fish_id)) %>%
    lmer(formula = sqrt(Time_in_centre) ~ (Genotype*Bin) + Trial + (1|gTrial) + (1|fish_id),
         data = .)


  plot(residuals(lm_time) ~ fitted(lm_time))
  qqPlot(residuals(lm_time))

  #----  overall time spent in centre
  print(emmeans(lm_time, ~ Genotype), type ="response") %>%
    as_tibble() %>%
    ggplot(aes(x = Genotype, y = response, colour = Genotype)) +
    geom_col(aes(fill = Genotype),
             alpha = 0.5,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = lower.CL, ymax =upper.CL),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_viridis_d(end = 0.75, option = "plasma") +
    scale_fill_viridis_d(end = 0.75, option = "plasma") +
    labs(y = "model-predicted time in centre",
         x = "Genotype") +
    theme(text = element_text(size = 20))


  #---- time spent in centre by bin.

  print(emmeans(lm_time, ~ Genotype * Bin), type = "response") %>%
    as_tibble() %>%
    mutate(binforvis = case_when(
      Bin == 1 ~ "Start-0:20:00",
      Bin == 2 ~ "0:20:00-0:40:00",
      Bin == 3 ~ "0:40:00-1:00:00",
    )) %>%
    ggplot(aes(x = binforvis, y = response, colour = Genotype)) +
    geom_col(aes(fill = Genotype),
             alpha = 0.5,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  size = 1,
                  position = position_dodge()) +
    theme(axis.text.x = element_text(hjust = 1,
                                     vjust = 1,
                                     angle = 45),
          legend.position = "bottom") +
    scale_color_viridis_d(end = 0.75) +
    scale_fill_viridis_d(end = 0.75) +
    facet_wrap(~Bin)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





# Freq in centre ----

  centre_time_freq_graph %>%
    ggplot(aes(x = Genotype, y = Freq_in_centre)) +
    geom_jitter(aes(colour = Genotype,
                    shape = Trial),
                size = 4) +
    geom_boxplot(aes(fill = Genotype),
                 alpha = 0.25,
                 outlier.shape = NA
    ) +
    scale_fill_viridis_d(end = 0.75) +
    scale_colour_viridis_d(end = 0.75) +
    scale_y_continuous(labels = comma, name = "Number of crosses") +
    ggtitle("Number of crosses into the central zone for 14dpf larvae in one hour") +
    facet_wrap (~Bin)

  ggsave(filename = "Number of crosses into the central zone for 14dpf larvae in one hour.png", width = 20, height = 10, units = "cm", dpi = 600, scale = 1.5)




#mixed effect model for the frequency in centre

  #---NO GOOD
  glm_crosses <-
  fulldata %>%
    glmer(Freq_in_centre ~ Genotype + Trial + (1|fish_id),
          data = .,
          family = poisson(link = "log"))

  check_overdispersion(glm_crosses) #hella overdispersed



  #negative binomial does not have assumption of equal dispersion so I will use that one instead
  glm_crosses <-
    fulldata %>%
    glmer.nb(Freq_in_centre ~ Genotype + Trial + Bin + (1|fish_id),
             data = .)

  Anova(glm_crosses) %>%
    as.data.frame() %>%
    dplyr::rename(pval = `Pr(>Chisq)`) %>%
    kable(caption = "Type II Wald chisquare test of LME") %>%
    kable_styling()


  #overall crosses into centre
  print(emmeans(glm_crosses, ~ Genotype), type ="response") %>%
    as_tibble() %>%
    ggplot(aes(x = Genotype, y = response, colour = Genotype)) +
    geom_col(aes(fill = Genotype),
             alpha = 0.5,
             width = 0.75,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax =asymp.UCL),
                  width = 0.125,
                  size = 1,
                  position = position_dodge(width = 0.75)) +
    scale_color_viridis_d(end = 0.75) +
    scale_fill_viridis_d(end = 0.75) +
    labs(y = "model-predicted crosses",
         x = "Genotype") +
    theme(text = element_text(size = 20))

  #---- Setting it up by bin



  print(emmeans(glm_crosses, ~ Genotype * Bin), type = "response") %>%
    as_tibble() %>%
    mutate(binforvis = case_when(
      Bin == 1 ~ "Start-0:20:00",
      Bin == 2 ~ "0:20:00-0:40:00",
      Bin == 3 ~ "0:40:00-1:00:00",

    )) %>%
    ggplot(aes(x = binforvis, y = response, colour = Genotype)) +
    geom_col(aes(fill = Genotype),
             alpha = 0.5,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  size = 1,
                  position = position_dodge()) +
    theme(axis.text.x = element_text(hjust = 1,
                                     vjust = 1,
                                     angle = 45),
          legend.position = "bottom") +
    scale_color_viridis_d(end = 0.75) +
    scale_fill_viridis_d(end = 0.75) +
    facet_wrap(~Bin)






#frequency in centre raw data
ggplot(data = centre_time_freq_graph, aes(x = Genotype, y = Freq_in_centre)) +
  geom_jitter(aes(colour = Genotype,
                     shape = Trial),
                         size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA
  ) +
  facet_wrap(~Bin)



#ggplate----

library(ggplate)
fulldata %>%
plate_plot(
  value = Distance_travelled,
  label = Distance_travelled,
  plate_size = 24,
  position = Position,
  plate_type = "round",
  scale = 1.5,
  label_size = 1.5
) +
  facet_wrap(~Trial ~Bin)














