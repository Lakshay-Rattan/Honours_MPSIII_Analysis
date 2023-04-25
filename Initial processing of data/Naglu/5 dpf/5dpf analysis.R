
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
  read_excel("5dpf EthoVision statistics.xlsx") %>%
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


distancedata.naglu.five %>% saveRDS("b.5.totaldist.rds")

data.naglu.five %>% saveRDS("B.5.fulldata.rds")



