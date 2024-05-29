#####################################################
##                                                 ##
##        Neuroptera temperature life-history      ##
##                                                 ##
##             Meta-analysis per trait             ##
##                                                 ##
##                JJ - 19/03/2024                  ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(meta)
library(metaviz)

load("coef_study.RData")

## Place to store mean data
res_meta <- data.frame(name = c("D 1st instar", "D 2nd instar", "D 3rd instar",
                            "D Pupae", "S pupae-adult", "#eggs/female"),
                   lwr = 0, mn = 0, upr = 0) %>% 
  mutate(name = factor(name, levels = c("D 1st instar", "D 2nd instar", "D 3rd instar",
                                        "D Pupae", "S pupae-adult", "#eggs/female")))

#_______________________________________________________________________________
#### 1. Meta-analyses for D 1st instar ####

d1inst <- coef_study %>% filter(name == "D 1st instar")

d1inst_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                        studlab = study, data = d1inst,
                        fixed = TRUE, random = FALSE, method.tau = "REML", 
                        title = "D 1st instar")
d1_sum <- summary(d1inst_meta)

jpeg(width = 20, height = 10, units = "cm", res = 600, file = "output/d1inst_meta.jpeg")
forest(d1inst_meta, 
       leftlabs = c("Study", "Total Insects"), 
       rightlabs = c("Coef", "95%-CI", "Weight"),
       col.square = "lightblue",
       smlab = "D 1st instar")
dev.off()

res_meta[1,2:4] <- c(d1_sum$common$lower, d1_sum$common$TE, d1_sum$common$upper)

#_______________________________________________________________________________
#### 2. Meta-analyses for D 2nd instar ####

d2inst <- coef_study %>% filter(name == "D 2nd instar")

d2inst_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                        studlab = study, data = d2inst,
                        fixed = TRUE, random = FALSE, method.tau = "REML", 
                        method.random.ci = "HK", title = "D 2nd instar")
d2_sum <- summary(d2inst_meta)

jpeg(width = 20, height = 10, units = "cm", res = 600, file = "output/d2inst_meta.jpeg")
forest(d2inst_meta, 
       leftlabs = c("Study", "Total Insects"), 
       rightlabs = c("Coef", "95%-CI", "Weight"),
       col.square = "lightblue",
       smlab = "D 2nd instar")
dev.off()

res_meta[2,2:4] <- c(d2_sum$common$lower, d2_sum$common$TE, d2_sum$common$upper)

#_______________________________________________________________________________
#### 3. Meta-analyses for D 3rd instar ####

d3inst <- coef_study %>% filter(name == "D 3rd instar")

d3inst_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                       studlab = study, data = d3inst,
                       fixed = TRUE, random = FALSE, method.tau = "REML", 
                       method.random.ci = "HK", title = "D 3rd instar")
d3_sum <- summary(d3inst_meta)

# jpeg(width = 20, height = 10, units = "cm", res = 600, file = "output/d3inst_meta.jpeg")
# forest(d3inst_meta, 
#        leftlabs = c("Study", "Total Insects"), 
#        rightlabs = c("Coef", "95%-CI", "Weight"),
#        col.square = "lightblue",
#        smlab = "D 3rd instar")
# dev.off()

res_meta[3,2:4] <- c(d3_sum$common$lower, d3_sum$common$TE, d3_sum$common$upper)

#_______________________________________________________________________________
#### 4. Meta-analyses for D Pupae ####

dpup <- coef_study %>% filter(name == "D Pupae")

dpup_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                     studlab = study, data = dpup,
                     fixed = TRUE, random = FALSE, method.tau = "REML", 
                     method.random.ci = "HK", title = "D Pupae")
dp_sum <- summary(dpup_meta)

# jpeg(width = 20, height = 10, units = "cm", res = 600, file = "output/dpup_meta.jpeg")
# forest(dpup_meta, 
#        leftlabs = c("Study", "Total Insects"), 
#        rightlabs = c("Coef", "95%-CI", "Weight"),
#        col.square = "lightblue",
#        smlab = "D Pupae")
# dev.off()

res_meta[4,2:4] <- c(dp_sum$common$lower, dp_sum$common$TE, dp_sum$common$upper)

#_______________________________________________________________________________
#### 5. Meta-analyses for S pupae-adult ####

sadult <- coef_study %>% filter(name == "S pupae-adult")

sadult_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                       studlab = study, data = sadult,
                       fixed = TRUE, random = FALSE, method.tau = "REML", 
                       method.random.ci = "HK", title = "S pupae-adult")
sa_sum <- summary(sadult_meta)

# jpeg(width = 20, height = 10, units = "cm", res = 600, file = "output/sadult_meta.jpeg")
# forest(sadult_meta, 
#        leftlabs = c("Study", "Total Insects"), 
#        rightlabs = c("Coef", "95%-CI", "Weight"),
#        col.square = "lightblue",
#        smlab = "S pupae-adult")
# dev.off()

res_meta[5,2:4] <- c(sa_sum$common$lower, sa_sum$common$TE, sa_sum$common$upper)

#_______________________________________________________________________________
#### 6. Meta-analyses for reproduction ####

repro <- coef_study %>% filter(name == "#eggs/female")

repro_meta <- metacor(n = n_insect, cor = coefficient, sm = "ZCOR", 
                       studlab = study, data = repro,
                       fixed = TRUE, random = FALSE, method.tau = "REML", 
                       method.random.ci = "HK", title = "#eggs/female")
rep_sum <- summary(repro_meta)

# jpeg(width = 20, height = 5, units = "cm", res = 600, file = "output/repro_meta.jpeg")
# forest(repro_meta, 
#        leftlabs = c("Study", "Total Insects"), 
#        rightlabs = c("Coef", "95%-CI", "Weight"),
#        col.square = "lightblue",
#        smlab = "#eggs/female")
# dev.off()

res_meta[6,2:4] <- c(rep_sum$common$lower, rep_sum$common$TE, rep_sum$common$upper)

#_______________________________________________________________________________
#### 7. Figure ####

cols_keep <- c("name", "studlab", "n", "TE", "lower", "upper")
meta_results_bind <- 
  bind_rows(data.frame(name = "D 1st instar", d1_sum)[,cols_keep],
          data.frame(name = "D 2nd instar", d2_sum)[,cols_keep],
          data.frame(name = "D 3rd instar", d3_sum)[,cols_keep],
          data.frame(name = "D Pupae", dp_sum)[,cols_keep],
          data.frame(name = "S pupae-adult", sa_sum)[,cols_keep],
          data.frame(name = "#eggs/female", rep_sum)[,cols_keep]) %>% 
  mutate(name = factor(name, levels = c("D 1st instar", "D 2nd instar", "D 3rd instar",
                                        "D Pupae", "S pupae-adult", "#eggs/female")))

fig_plot <- ggplot(meta_results_bind, aes(x = TE, y = studlab)) +
  geom_vline(xintercept = 0) +
  geom_vline(data = res_meta, aes(xintercept = mn), colour = "steelblue") +
  geom_rect(data = res_meta, aes(x = NULL, xmin = lwr, xmax = upr, 
                                 y = NULL, ymin = 0, ymax = Inf),
            fill = "lightsteelblue", alpha = 0.5) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
  geom_point(aes(size = n)) +
  scale_size_continuous(breaks = c(50, 150, 200, 400, 1500), 
                        range = c(2,6)) +
  #scale_x_continuous(breaks = seq(-0.4, 0.4, by = 0.2), limits = c(-0.4,0.4)) +
  facet_wrap(~ name, scales = "free_x", ncol = 6) +
  labs(x = "Temperature sensitivity", 
       y = "Study", size = "Number of\nInsects") +
  theme_bw(base_size = 15) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) 

# ggsave(fig_plot, filename = "output/meta_analysis_results.jpg",
#        width = 32, height = 12, units = "cm", dpi = 800)


