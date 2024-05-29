#####################################################
##                                                 ##
##        Neuroptera temperature life-history      ##
##                                                 ##
##             Coefficients per study              ##
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

load("neuroptera.RData", verbose = T)

#_______________________________________________________________________________
#### 2. Identifying studies to retain ####

study_sum <- neuroptera %>% 
  group_by(study) %>% 
  summarise(n_species = n_distinct(species),
            n_temperatures = n_distinct(temperature),
            n_obs = n(),
            n_1st = length(which(is.na(development_time_1inst) == F)),
            n_2nd = length(which(is.na(development_time_2inst) == F)),
            n_3rd = length(which(is.na(development_time_3inst) == F)),
            n_pup = length(which(is.na(development_time_pupae) == F)),
            n_surv = length(which(is.na(larval_survival) == F)),
            n_repr = length(which(is.na(reproductive_rate) == F)))

## 
study_sum2 <- neuroptera %>% 
  group_by(study) %>% 
  summarise(n_species = n_distinct(species),
            n_temperatures = n_distinct(temperature),
            n_obs = n(),
            n_1st = length(which(is.na(development_time_1inst) == F)),
            n_2nd = length(which(is.na(development_time_2inst) == F)),
            n_3rd = length(which(is.na(development_time_3inst) == F)),
            n_pup = length(which(is.na(development_time_pupae) == F)),
            n_surv = length(which(is.na(larval_survival) == F)),
            n_repr = length(which(is.na(reproductive_rate) == F)))

## Currently keeping only studies with:
## > 1 temperature treatment
## > 3 observations

retained_studies <- study_sum %>% 
  mutate(development_time_1inst = if_else(n_temperatures > 1 & n_1st > 3, 1, 0),
         development_time_2inst = if_else(n_temperatures > 1 & n_2nd > 3, 1, 0),
         development_time_3inst = if_else(n_temperatures > 1 & n_3rd > 3, 1, 0),
         development_time_pupae = if_else(n_temperatures > 1 & n_pup > 3, 1, 0),
         larval_survival = if_else(n_temperatures > 1 & n_surv > 3, 1, 0),
         reproductive_rate = if_else(n_temperatures > 1 & n_repr > 3, 1, 0)) %>% 
  dplyr::select(study, development_time_1inst:reproductive_rate) %>% 
  pivot_longer(-study) %>% 
  mutate(id = paste0(study, "_", name)) %>% 
  filter(value == 1)

#_______________________________________________________________________________
#### 3. Histograms for retained studies ####

neuroptera %>% 
  dplyr::select(study, development_time_1inst:reproductive_rate) %>% 
  pivot_longer(-study) %>% 
  mutate(id = paste0(study, "_", name)) %>% 
  filter(id %in% retained_studies$id) %>% 
  na.omit() %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(bins = 20, show.legend = F) +
  facet_wrap(~ name, scales = "free")
  
#_______________________________________________________________________________
#### 4. Iterating through studies and LH variables ####

# development time + reproduction - log transformed linear regression
# survival - logistic regression

neuroptera_coefficients <- retained_studies %>% 
  dplyr::select(study, name) %>% 
  mutate(coefficient = 0, n = 0, n_insect = 0, se = 0, se_coef = 0) %>% 
  # Studies with all NAs except for one temperature
  filter(study %in% c("ID_28", "ID_39") == F)
  
resid_list <- vector(mode = "list")

for(i in 1:nrow(neuroptera_coefficients)){
  
  # print(i)

  c_study = neuroptera_coefficients$study[i]
  c_resp = neuroptera_coefficients$name[i]
  
  # get right data
  c_dat_raw = filter(neuroptera, study == c_study) 
  
  c_dat = filter(neuroptera, study == c_study) %>% 
    dplyr::select(response = paste(c_resp), temperature) %>% 
    na.omit()
  
  # binomial if a survival model
  if(c_resp == "larval_survival"){mod = glm(response/100 ~ temperature, data = c_dat)}
  else{mod = lm(log(response + 1) ~ temperature, data = c_dat)}
  
  neuroptera_coefficients[i,]$coefficient <- summary(mod)$coefficients[2,1]
  neuroptera_coefficients[i,]$se_coef <- summary(mod)$coefficients[2,2]
  
  neuroptera_coefficients[i,]$n <- nrow(c_dat)
  neuroptera_coefficients[i,]$n_insect <- floor(sum(c_dat_raw$Instances_Count))
  neuroptera_coefficients[i,]$se <- sd(c_dat$response)/sqrt(nrow(c_dat))
  
  resid_list[[i]] <- tibble(fitted = fitted(mod), residual = resid(mod))
  
}

#_______________________________________________________________________________
#### 5. Summary plot - basis for metaregression ####

coef_study <- neuroptera_coefficients %>% 
  mutate(name = case_when(name == "development_time_1inst" ~ "D 1st instar",
                          name == "development_time_2inst" ~ "D 2nd instar",
                          name == "development_time_3inst" ~ "D 3rd instar",
                          name == "development_time_pupae" ~ "D Pupae",
                          name == "larval_survival" ~ "S pupae-adult",
                          name == "reproductive_rate" ~ "#eggs/female"),
         name = factor(name, levels = c("D 1st instar", "D 2nd instar", "D 3rd instar",
                                        "D Pupae", "S pupae-adult", "#eggs/female"))) %>% 
  na.omit() # taking away single study without n insects recorded

coef_summary <- coef_study %>% 
  ggplot(aes(x = coefficient, y = study)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmax = coefficient + 1.96*se_coef, 
                    xmin = coefficient - 1.96*se_coef), width = 0) +
  geom_point(aes(size = n_insect)) +
  scale_size_continuous(breaks = c(50, 150, 200, 400, 1500), 
                        range = c(2,6)) +
  facet_wrap(~ name, scales = "free_x", ncol = 6) +
  labs(x = "Temperature coefficient (linear predictor)", 
       y = "Study", size = "Number of\nInsects") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank())

# Plotting - commented out for ease
# ggsave(coef_summary, filename = "output/linear_coefficients.jpeg", 
#        width = 20, height = 18, units = "cm", dpi = 700)

# jpeg(filename = "output/linear_coefficient_residuals.jpeg", width = 40, height = 28, units = "cm", res = 700)
# par(mfrow = c(6,10))
# for(i in 1:length(resid_list)){
#   plot(resid_list[[i]]$residual ~ resid_list[[i]]$fitted, 
#        xlab = "Fitted", ylab = "Residual", 
#        main = paste0(neuroptera_coefficients[i,]$study, "_", neuroptera_coefficients[i,]$name))
#   abline(b = 0, a = 0)
# }
# dev.off()

save(coef_study, file = "coef_study.RData")

