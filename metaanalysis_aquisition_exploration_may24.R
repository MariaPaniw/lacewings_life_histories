#####################################################
##                                                 ##
##        Neuroptera temperature life-history      ##
##                                                 ##
##              Cleaning + Exploration             ##
##                                                 ##
##                JJ - 13/05/2024                  ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)

#_______________________________________________________________________________
#### 1. load data ####

neuroptera_raw <- read_csv("LH_Neuroptera_red.csv")

#_______________________________________________________________________________
#### 2. looking at data spread per variable ####

nsum <- neuroptera_raw %>% 
  dplyr::select(-(tmin_E:Food), -Dorm.) %>% 
  mutate(Pr_sex_rat = as.numeric(gsub(",", ".", Pr_sex_rat)),
         Intr_rate_nat_incr = as.numeric(gsub(",", ".", Intr_rate_nat_incr))) %>% 
  pivot_longer(cols = Preovi_per:Vert_distr_M) %>% 
  filter(is.na(temp) == F & is.na(value) == F) %>% 
  group_by(name) %>% 
  summarise(n_obs = n(), n_species = n_distinct(sp.)) %>% 
  arrange(-n_obs) 

nsum$name <- factor(nsum$name, levels = nsum$name)

a <-
  ggplot(nsum, aes(x = name, y = n_obs)) +
  geom_col() +
  coord_flip() +
  labs(y = "Number of Observations", x = "Column name") +
  theme_bw()

b <- 
  ggplot(nsum, aes(x = name, y = n_species)) +
  geom_col() +
  coord_flip() +  
  labs(y = "Number of Species", x = "Column name") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

# ggsave(a + b + plot_layout(widths = c(4,3)), filename = "output/summary_plot.jpeg",
#        width = 17, height = 15, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 3. Cleaning + exploration ####

neuroptera <- neuroptera_raw %>% 
  mutate(ID = 1:n()) %>% 
  dplyr::select(ID:temp, Dev_1st_inst:L_surv, Repr_rate) %>% 
  rename(observation = ID, study = St_ID, species = sp., continent = `?ontinent`, 
         place = Place, habitat = Hab., temperature = temp, 
         development_time_1inst = Dev_1st_inst,
         development_time_2inst = Dev_2nd_inst, 
         development_time_3inst = Dev_3rd_inst,
         development_time_pupae = `Dev_ P`, 
         larval_survival = L_surv, 
         reproductive_rate = Repr_rate) %>% 
  filter(is.na(temperature) == F) %>% 
  mutate(species_lab = gsub("_", " ", species))

ht <- neuroptera %>% 
  dplyr::select(observation, species, development_time_1inst:reproductive_rate) %>% 
  pivot_longer(-c(observation, species)) %>% 
  na.omit() %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(show.legend = F, colour = "black", bins = 20, linewidth = 0.25) +
  facet_wrap(~ name, scales = "free") +
  labs(x = "Value", y = "Frequency") +
  theme_bw()

# ggsave(ht, filename = "output/histograms.jpeg", width = 19, height = 13, units = "cm", dpi = 600)

## plots
d1 <- ggplot(neuroptera, aes(x = temperature, y = development_time_1inst, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "1st Instar Development time", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

d2 <- ggplot(neuroptera, aes(x = temperature, y = development_time_2inst, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "2nd Instar Development time", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

d3 <- ggplot(neuroptera, aes(x = temperature, y = development_time_3inst, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "3rd Instar Development time", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

dp <- ggplot(neuroptera, aes(x = temperature, y = development_time_pupae, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "Pupae Development time", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

ls <- ggplot(neuroptera, aes(x = temperature, y = larval_survival, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "Larval survival (%)", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

rr <- ggplot(neuroptera, aes(x = temperature, y = reproductive_rate, colour = study)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Temperature (\u00B0C)", y = "Female reproductive rate", colour = "Study ID") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

#_______________________________________________________________________________
#### 4. Summary information per study ####

study_sum <- neuroptera %>% 
  mutate(study_f = 
           factor(study, 
                  levels = unique(.$study[order(as.numeric(gsub(pattern = "ID_", "", .$study)))]))) %>% 
  group_by(study_f) %>% 
  summarise(n_species = n_distinct(species),
            n_temperatures = n_distinct(temperature),
            n_obs = n())

sa <- ggplot(study_sum, aes(x = study_f, y = n_obs)) +
  geom_col() +
  coord_flip() +
  labs(y = "Number of Observations", x = "Study ID") +
  theme_bw()

sb <- ggplot(study_sum, aes(x = study_f, y = n_species)) +
  geom_col() +
  coord_flip() +
  labs(y = "Number of Species", x = "Study ID") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

sc <- ggplot(study_sum, aes(x = study_f, y = n_temperatures)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0,16,by = 2)) +
  coord_flip() +
  labs(y = "Number of Temperatures", x = "Study ID") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

# ggsave(sa + sb + sc + plot_layout(widths = c(9,8,8)),
#        filename = "output/study_summary.jpeg", width = 19, 
#        height = 16, units = "cm", dpi = 600)

#_______________________________________________________________________________
#### 5. Save ####

save(neuroptera, file = "neuroptera.RData")






