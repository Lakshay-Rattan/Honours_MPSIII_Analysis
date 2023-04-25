
# Packages

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
library(dlookr)


# data visualstion
library(scales)
library(ggpubr)
library(ggeasy)
library(ggfortify)
library(ggbeeswarm)
library(kableExtra)
library(see)
library(ggplate)

theme_set(theme_bw())




distancedata.naglu.five %>%
  group_by(Genotype) %>%
  summarise(
    n = n(),
    mean = mean(Total_distance),
    median = median(Total_distance))


# import data ------------------------------------------------------------

# The data outputted from the Daniovision,
# as well as the metadata about each larvae are found in the naglu het x het... folder
# here i will import the data into R

data.naglu.five <-
  read_excel("5dpf EthoVision statistics after changing threshold.xlsx") %>%
  # renaming the column names to more informatative
  dplyr::select(Trial = `...2`,
                Position = `...3`,
                Bin = `...4`,
                Distance_travelled = "Distance moved center-point Total mm",
                cum_movement = "Movement Moving / center-point Cumulative Duration s",
                Freq_in_centre = "In zone Centre / center-point Frequency",
                Time_in_centre =  "In zone Centre / center-point Cumulative Duration s"



  ) %>%
  # here, I'm converting these cols factors, and tidying up the names
  mutate(Bin = factor(Bin, levels = c("Start-0:20:00", "0:20:00-0:40:00", "0:40:00-1:00:00")),
         Trial = str_remove(Trial, pattern = "Trial     ") %>%
           as.factor()
  )

# read in metadata
meta.naglu.five <-
  read_excel("5 dpf 72 samples genotype.xlsx") %>%
  mutate(Genotype = factor(Genotype, levels = c("wt", "het", "hom")),
         Trial = as.factor(Trial),
         Position = as.factor(Position)
  )

# join the spreadsheets together
data.naglu.five <- meta.naglu.five %>%
  left_join(data.naglu.five, by = c("Trial", "Position")) %>%
  na.omit




# check genotype propotions --------

#  The first check I do isto make sure that the fish in the family follow
# typical Mendelian ratios. No deviations from this are obsvered

meta.naglu.five %>%
  na.omit %>%
  ggplot(aes(x = Genotype, fill = Genotype)) +
  geom_bar(colour = "black") +
  scale_fill_viridis_d() + # pretty viridis colour scheme
  #ggtitle("Number of embryos in the analysis family by genotype") +
  ylab("Number of larvae")

ggsave(filename = "5dpf genotype proportion.png", width = 10, height = 10, units = "cm", dpi = 300, scale = 1)

# Analysis ------------
## Total distance travelled  ------------
#  In the data object, there are 3 timepoints per fish.
# I want to sum this up before I plot the graph

# Plot the data as a boxplot and raw data points



distancedata.naglu.five %>%
  (aes(x = Genotype, y = Total_distance)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA # do not show outliers, they are already plotted by geom_point
  ) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (mm)") +
  theme(text = element_text(size = 15))
# ggtitle("Total distance moved by larvae in one hour(mm)")

ggsave(filename = "Total distance in an hour.png", width = 15, height = 10, units = "cm", dpi = 300, scale = 1.5)




# make a new object which only contains the sum of the total distance travelled.
distancedata.naglu.five <- data.naglu.five %>%
  group_by(fish_id) %>%
  mutate(Total_distance = sum(Distance_travelled)) %>%
  dplyr::distinct(fish_id, .keep_all = TRUE) %>%
  dplyr::select(fish_id, Trial, Position, Genotype, Total_distance, Freq_in_centre) %>%
  ungroup() %>%
  unique() %>%
  droplevels()

#exporting the data as an rds file----


distancedata.naglu.five %>% saveRDS("/Users/Lakshay/Desktop/Honours/R data/One file to rule them all/b.5.totaldist.rds")

data.naglu.five %>% saveRDS("/Users/Lakshay/Desktop/Honours/R data/One file to rule them all/B.5.fulldata.rds")



# Fit a linear mixed effect model ---------------------
model_distance.naglu.five <- distancedata.naglu.five %>%
  group_by(fish_id) %>%
  #mutate(Distance_travelled = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Trial, Genotype, Distance_travelled) %>%
  ungroup() %>%
  #unique() %>%
  droplevels() %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
  lmer(formula = sqrt(Distance_travelled) ~ Genotype + Trial + (1|gTrial))


check_model(model_distance.naglu.five)
#checking homogeneity of variances using robust (levene) and strictly normal data tests.


# check model assumptions - not partocualrly perfect
plot(residuals(model_distance.naglu.five) ~ fitted(model_distance.naglu.five))
qqPlot(residuals(model_distance.naglu.five))



# Test the effect of genotype using a type of chi sq test.
Anova(model_distance.naglu.five) %>%
  as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()





#distance data per time bin----


data.naglu.five %>%
  group_by(fish_id) %>%
  dplyr::select(fish_id, Trial, Genotype, Distance_travelled, Bin) %>%
  ggplot(aes(x = Genotype, y = Distance_travelled)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.25,
               outlier.shape = NA) +
  scale_fill_viridis_d(end = 0.75) +
  scale_colour_viridis_d(end = 0.75) +
  scale_y_continuous(labels = comma, name = "Distance (mm)") +
  #ggtitle("Total distance travelled for larvae over three 20 minute time bins") +
  facet_wrap (~Bin) +
  theme(text = element_text(size = 15))

ggsave(filename = "Distance per bin.png", width = 20, height = 14, units = "cm", dpi = 300, scale = 1.5)



#next big bit is me trying to get normally distributed data. non of my methods worked,
#changing approach by performing analysis per individual bins - stopped doing that too.

with(data.naglu.five, hist(Distance_travelled))
#data.naglu.five <-
#trans <-
#data.naglu.five %>%
# mutate(sqrtdist = sqrt(Distance_travelled),
#      logdist = log(Distance_travelled),
#      b.cox = (Distance_travelled^(-2/3) - 1) / (-2/3), #refer to lambda calculation
#    invdist = (1/Distance_travelled))

#     trans %>%
#    with(hist(b.cox))

# normality(trans, Distance_travelled)
#normality(trans, sqrtdist)
#normality(trans, logdist)
#normality(trans, b.cox)
#normality(trans, invdist)

#working out for lamba for box cox transformation.
#b <- boxcox(lm(trans$Distance_travelled ~ 1))

#lambda <- b$x[which.max(b$y)]

#to recap, my data is not normally distributed, therefore, I will use Kruskal-Wallis test.

#so p value is 0.3206 BUT this does not account for the random effect of trial. Cant make linear model tho :/
kruskal.test(Distance_travelled ~ Genotype, data = data.naglu.five)








#Checking assumptions for time bin 1-------

time1dist <-
  data.naglu.five %>%
  filter(Bin == "Start-0:20:00")

#not normally distributed
time1dist %>%
  normality(Distance_travelled)

time1dist <-
  time1dist%>%
  mutate(sqrtdist = sqrt(Distance_travelled),
         logdist = log(Distance_travelled),
         b.cox = (Distance_travelled^(-2/3) - 1) / (-2/3), #refer to lambda calculation
         invdist = (1/Distance_travelled)) #%>%
#with(hist(b.cox))

#normality(b.cox)

b1 <- boxcox(lm(time1dist$Distance_travelled ~ 1))

lambda1 <- b1$x[which.max(b$y)]

#BOX COX IS NORMAL IM GOOD TO GO

#checking homogeneity of variances using robust (levene) and barlett test looking good so far.
leveneTest(b.cox ~ Genotype, time1dist)
bartlett.test(b.cox ~ Genotype, time1dist)



# Fitting a normal linear model. no need for mixed effect as ramdomness of trials not a factor
model_time1dist <- time1dist %>%
  group_by(fish_id) %>%
  #mutate(Distance_travelled = sum(Distance_travelled)) %>%
  dplyr::select(fish_id, Genotype, b.cox) %>%
  ungroup() %>%
  #unique() %>%
  droplevels() %>%


  #ok this only shows results for two levels (hom and het), need to use multiple linear model with separate columns for each genotype
  geno_time1dist <-
  time1dist %>%
  mutate(WT = Genotype == 'WT')
# OK I GIVE UP, ASK KARISSA WHATS GOING ON
# probs going wrong direction,



Chi_model_time1dist <-
  Anova(b.cox ~ Genotype, model_time1dist)
as.data.frame()
dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of anova") %>%
  kable_styling()





#ended up stopping this because s size may be too small.
# Time against the wall ----------

# Vis the data
rata %>%
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
  scale_fill_viridis_d('0.75') +
  scale_colour_viridis_d('0.75') +
  scale_y_continuous(labels = comma, name = "Time spent in middle (s)")

# Plot the data as mean bar and raw data points
rata %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(1:4, time_at_middle) %>%
  unique() %>% #as there were 3 replicates because of the bins
  ggplot(aes(x = Genotype, y = time_at_middle)) +
  stat_summary(fun=mean ,
               geom="bar",
               alpha = 0.5,
               aes(group = Genotype,
                   fill = Genotype)) +
  geom_jitter(aes(colour = Genotype,
                  shape = Trial),
              size = 4) +
  scale_fill_viridis_d("magma") +
  scale_colour_viridis_d('magma') +
  scale_y_continuous(labels = comma, name = "Time spent in middle (s)")
ggsave(filename = "av_middle.png", width = 13, height = 15, units = "cm", dpi = 400, scale = 0.9)
#ggtitle("Mean time larvae spent in the middle zone of the well in a one hour trial")

# make a new object which only contains the sum of the total distance travelled.
timedata <- rata %>%
  group_by(fish_id) %>%
  mutate(time_at_middle = sum(time_at_middle)) %>%
  dplyr::select(fish_id, Trial, Genotype, time_at_middle) %>%
  ungroup() %>%
  unique() %>%
  droplevels()

# Look at the histogram.
# want something resembling a normal dist so that we can pply a lin model.
#  Looks fairly normal to me, a bit L skew
with(timedata, hist(time_at_middle))

# DISCLAIMER: I'm not overly sure this is the correct statistical test.
# Need to consult a statistician.

# Fit a linear mixed effect model
model_time <-
  timedata %>%
  mutate(gTrial = interaction(Genotype, Trial, drop = TRUE)) %>% # these are to define interaction effects in the model
  lmer(formula = (time_at_middle)^(1/3) ~ Genotype + Trial + (1|gTrial))

# check model assumptions
plot(residuals(model_time) ~ fitted(model_time))
qqPlot(residuals(model_time))

# not too bad. A linear model should be appropriate

# Test the effect of genotype using a type of chi sq test.
Anova(model_time)
as.data.frame() %>%
  dplyr::rename(pval = `Pr(>Chisq)`) %>%
  kable(caption = "Type II Wald chisquare test of LME") %>%
  kable_styling()
