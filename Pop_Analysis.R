rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
               robumeta, ggpmisc, ggridges, ggbeeswarm, gridExtra)

source("./func.R")

# Importing Data Set
data <- read.csv("./Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

# Phylogenetic covariance matrix
tree <- ape::read.tree("./tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

##### Overall Model #####
priors <-  prior(student_t(3, 0, 20), class = "sd")

system.time(
  overall <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                       ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|obs),
                       data = data,
                       family = gaussian(),
                       data2 = list(A = A), 
                       chains = 4, 
                       cores = 4,
                       iter = 12000,
                       warmup = 2000,
                       thin = 5,
                       prior = priors,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       file = "./overall_model",
                       file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_overall <- as_draws_df(overall, variable = "b_Intercept")
b_overall <- data.frame(b_overall$b_Intercept)

sd <- as_draws_df(overall, variable = c("sd_Measurement__Intercept", "sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd <- data.frame("sd_Measurement__Intercept" = sd$sd_Measurement__Intercept, 
                 "sd_obs__Intercept" = sd$sd_obs__Intercept, 
                 "sd_phylo__Intercept" = sd$sd_phylo__Intercept, 
                 "sd_Study_ID__Intercept" = sd$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b <-  sapply(b_overall, mean)
ci_b <- as.vector(HPDinterval(as.mcmc(b_overall)))
pMCMC_b <- 2*(1 - max(table(b_overall<0) / nrow(b_overall)))

# Absolute magnitude
b_abs <- folded_norm(b_overall[,1], rowSums(sd))
mean_abs_b <- mean(b_abs)
ci.abs <- HPDinterval(as.mcmc(b_abs))

# Heterogeneity
overall_i2 <- i2(sd, data$Variance_Adjusted) 

##### Overall Model - Plasticity Mechanism Meta-regression #####
Plasticity_Exploration <- data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Plasticity_Exploration) <- Plasticity_Exploration$Plasticity_Category

Plasticity_Species_Count <- data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                            table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                            select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Plasticity_Species_Count) <- Plasticity_Species_Count$Plasticity_Category

Plasticity_Study_Count <- data %>% select("Study_ID", "Plasticity_Category") %>% 
                          table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                          select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Plasticity_Study_Count) <- Plasticity_Study_Count$Plasticity_Category

system.time(
  overall_plastic <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                               ~ Plasticity_Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Plasticity_Category, cor = FALSE)),
                               data = data,
                               family = gaussian(),
                               data2 = list(A = A), 
                               chains = 4, 
                               cores = 4,
                               iter = 12000,
                               warmup = 2000,
                               thin = 5,
                               prior = priors,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               file = "./overall_plastic_model",
                               file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_plastic <- as_draws_df(overall_plastic, variable = c("b_Intercept", 
                                                       "b_Plasticity_CategoryDevelopmental", 
                                                       "b_Plasticity_CategoryTransgenerational"))
b_plastic <- data.frame("b_Acclimation" = b_plastic$b_Intercept, 
                        "b_Developmental" = b_plastic$b_Plasticity_CategoryDevelopmental + b_plastic$b_Intercept, 
                        "b_Transgenerational" = b_plastic$b_Plasticity_CategoryTransgenerational + b_plastic$b_Intercept)


sd_plastic <- as_draws_df(overall_plastic, variable = c("sd_obs__Intercept:Plasticity_CategoryAcclimation",
                                                        "sd_obs__Intercept:Plasticity_CategoryDevelopmental",
                                                        "sd_obs__Intercept:Plasticity_CategoryTransgenerational",
                                                        "sd_phylo__Intercept", "sd_Study_ID__Intercept", 
                                                        "sd_Measurement__Intercept"))
sd_plastic <- data.frame("sd_Acclimation" = sd_plastic$`sd_obs__Intercept:Plasticity_CategoryAcclimation`,
                         "sd_Developmental" = sd_plastic$`sd_obs__Intercept:Plasticity_CategoryDevelopmental`, 
                         "sd_Transgenerational" = sd_plastic$`sd_obs__Intercept:Plasticity_CategoryTransgenerational`,
                         "sd_phylo__Intercept" = sd_plastic$`sd_phylo__Intercept`, 
                         "sd_Study_ID__Intercept" = sd_plastic$`sd_Study_ID__Intercept`, 
                         "sd_Measurement__Intercept" = sd_plastic$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
overall_plastic_means <- apply(b_plastic, 2, mean)
overall_plastic_cis <- apply(b_plastic, 2, function(x) HPDinterval(as.mcmc(x)))
overall_plastic_pMCMC <- apply(b_plastic, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_plastic_acc <- folded_norm(b_plastic$b_Acclimation, sqrt(rowSums(sd_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_plastic[, "sd_Acclimation"]^2)))
b_abs_plastic_dev <- folded_norm(b_plastic$b_Developmental, sqrt(rowSums(sd_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_plastic[, "sd_Developmental"]^2)))
b_abs_plastic_trans <- folded_norm(b_plastic$b_Transgenerational, sqrt(rowSums(sd_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_plastic[, "sd_Transgenerational"]^2)))
mean_abs_b_plastic_acc <- mean(b_abs_plastic_acc)
mean_abs_b_plastic_dev <- mean(b_abs_plastic_dev)
mean_abs_b_plastic_trans <- mean(b_abs_plastic_trans)
ci.abs_acc <- HPDinterval(as.mcmc(b_abs_plastic_acc))
ci.abs_dev <- HPDinterval(as.mcmc(b_abs_plastic_dev))
ci.abs_trans <- HPDinterval(as.mcmc(b_abs_plastic_trans))

# Heterogeneity
overall_i2_plastic <- i2(sd_plastic, data$Variance_Adjusted) 

# Overall_Plasticity_Summary
overall_plastic_means_list <- c(mean_abs_b_plastic_acc, mean_abs_b_plastic_dev, mean_abs_b_plastic_trans)
overall_plastic_low_ci <- c(ci.abs_acc[1], ci.abs_dev[1], ci.abs_trans[1])
overall_plastic_high_ci <- c(ci.abs_acc[2], ci.abs_dev[2], ci.abs_trans[2])
overall_plastic_categories <- c("Acclimation", "Developmental Plasticity", "Transgenerational Effects")

overall_plastic_summary <- matrix(c(overall_plastic_means_list, overall_plastic_low_ci, overall_plastic_high_ci), 
                                  nrow = 3, ncol = 3, byrow = FALSE, 
                                  dimnames = list(c(overall_plastic_categories), 
                                                  c("Mean", "Low_CI", "High_CI")))
overall_plastic_summary <- data.frame(overall_plastic_summary)

# Preparing Graph - Combined

Plasticity_rnames <- c("Acclimation", "Developmental Plasticity", "Transgenerational Effects")

Plasticity_k <- data.frame("k" = c(Plasticity_Exploration["Acclimation", "Freq"], 
                                   Plasticity_Exploration["Developmental", "Freq"], 
                                   Plasticity_Exploration["Transgenerational", "Freq"]), 
                           row.names = Plasticity_rnames)

Plasticity_group_no <- data.frame("Spp No." = c(Plasticity_Species_Count["Acclimation", "Freq"], 
                                                Plasticity_Species_Count["Developmental", "Freq"], 
                                                Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                  row.names = Plasticity_rnames)

Plasticity_study <- data.frame("Study" = c(Plasticity_Study_Count["Acclimation", "Freq"], 
                                           Plasticity_Study_Count["Developmental", "Freq"], 
                                           Plasticity_Study_Count["Transgenerational", "Freq"]), 
                               row.names = Plasticity_rnames)

Plasticity_table <- data.frame(estimate = overall_plastic_summary[,"Mean"], 
                               lowerCL = overall_plastic_summary[,"Low_CI"], 
                               upperCL = overall_plastic_summary[,"High_CI"], 
                               K = Plasticity_k[,1], 
                               group_no = Plasticity_group_no[,1], 
                               row.names = Plasticity_rnames)
Plasticity_table$name <- row.names(Plasticity_table)

Plasticity_raw_mean <- c(b_abs_plastic_acc, b_abs_plastic_dev, b_abs_plastic_trans)

Plasticity_raw_name <- c(replicate(8000, "Acclimation"), 
                         replicate(8000, "Developmental Plasticity"), 
                         replicate(8000, "Transgenerational Effects"))

Plasticity_raw_df <- data.frame("Model" = Plasticity_raw_name, 
                                "Effect" = Plasticity_raw_mean)

# Graph code - Combined

Plasticity_Order <- c("Transgenerational Effects", "Developmental Plasticity", "Acclimation")

density_plasticity <- Plasticity_table %>% mutate(name = fct_relevel(name, Plasticity_Order)) %>%
                      ggplot() +
                      geom_density_ridges(data = Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Plasticity_Order)), 
                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                      geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table)[1], 1)), xmin = max(Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table)[1], 1)), xmin = min(Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                     size = 1) +
                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                      size = 1, fatten = 2) +
                      theme_bw() +
                      guides(fill = "none", colour = "none") +
                      labs(x = TeX("Effect Size (PRRD)"), y = "") +
                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                       vjust = c(-0.8, -0.8, -2.7))) +
                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                      theme(axis.ticks = element_blank()) +
                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                      coord_cartesian(xlim = c(-0.01, 1.25)) +
                      annotate('text',  x = 1.25, y = (seq(1, dim(Plasticity_table)[1], 1)+0.4),
                      label= paste("italic(k)==", c(Plasticity_table["Transgenerational Effects", "K"], 
                                                    Plasticity_table["Developmental Plasticity", "K"], 
                                                    Plasticity_table["Acclimation", "K"]), "~","(", 
                                                  c(Plasticity_table["Transgenerational Effects", "group_no"], 
                                                    Plasticity_table["Developmental Plasticity", "group_no"], 
                                                    Plasticity_table["Acclimation", "group_no"]), 
                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                      geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_plastic_trans)-1)*100, 2), nsmall = 2), "%"), 
                                             paste(format(round(mean(exp(b_abs_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                             paste(format(round(mean(exp(b_abs_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                 x = rev(Plasticity_table$estimate+0.2), y = (seq(1, dim(Plasticity_table)[1], 1)+0.4)), size = 3.5)

density_plasticity #(400x320)

# Preparing Graph - Part 1

Plasticity_rnames_1 <- c("Acclimation", "Developmental Plasticity")

Plasticity_k_1 <- data.frame("k" = c(Plasticity_Exploration["Acclimation", "Freq"], 
                                     Plasticity_Exploration["Developmental", "Freq"]), 
                             row.names = Plasticity_rnames_1)

Plasticity_group_no_1 <- data.frame("Spp No." = c(Plasticity_Species_Count["Acclimation", "Freq"], 
                                                  Plasticity_Species_Count["Developmental", "Freq"]), 
                                    row.names = Plasticity_rnames_1)

Plasticity_study_1 <- data.frame("Study" = c(Plasticity_Study_Count["Acclimation", "Freq"], 
                                             Plasticity_Study_Count["Developmental", "Freq"]), 
                                 row.names = Plasticity_rnames_1)

overall_plastic_summary_Reorder_1 <- overall_plastic_summary[c("Acclimation", "Developmental Plasticity"), ]

Plasticity_table_1 <- data.frame(estimate = overall_plastic_summary_Reorder_1[,"Mean"], 
                                 lowerCL = overall_plastic_summary_Reorder_1[,"Low_CI"], 
                                 upperCL = overall_plastic_summary_Reorder_1[,"High_CI"], 
                                 K = Plasticity_k_1[,1], 
                                 group_no = Plasticity_group_no_1[,1], 
                                 row.names = Plasticity_rnames_1)
Plasticity_table_1$name <- row.names(Plasticity_table_1)

Plasticity_raw_mean_1 <- c(b_abs_plastic_acc, b_abs_plastic_dev)

Plasticity_raw_name_1 <- c(replicate(8000, "Acclimation"), 
                           replicate(8000, "Developmental Plasticity"))

Plasticity_raw_df_1 <- data.frame("Model" = Plasticity_raw_name_1, 
                                  "Effect" = Plasticity_raw_mean_1)

# Graph code - Part 1

Plasticity_Order_1 <- c("Developmental Plasticity", "Acclimation")

density_plasticity_1 <- Plasticity_table_1 %>% mutate(name = fct_relevel(name, Plasticity_Order_1)) %>%
                        ggplot() +
                        geom_density_ridges(data = Plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Plasticity_Order_1)), 
                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                       size = 1) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_1)[1], 1)), xmin = max(Plasticity_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                       size = 1) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_1)[1], 1)), xmin = min(Plasticity_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                       size = 1) +
                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                        size = 1, fatten = 2) +
                        theme_bw() +
                        guides(fill = "none", colour = "none") +
                        labs(x = TeX("Effect Size (PRRD)"), y = "") +
                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                         vjust = c(-0.8, -2.7))) +
                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                        theme(axis.ticks = element_blank()) +
                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                        scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                        scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                        coord_cartesian(xlim = c(-0.01, 1.25)) +
                        annotate('text',  x = 1.25, y = (seq(1, dim(Plasticity_table_1)[1], 1)+0.4),
                        label= paste("italic(k)==", c(Plasticity_table_1["Developmental Plasticity", "K"], 
                                                      Plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                    c(Plasticity_table_1["Developmental Plasticity", "group_no"], 
                                                      Plasticity_table_1["Acclimation", "group_no"]), 
                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                        geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                               paste(format(round(mean(exp(b_abs_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                       x = rev(Plasticity_table_1$estimate+0.2), y = (seq(1, dim(Plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_plasticity_1 #(400x240)

# Preparing Graph - Part 2

Plasticity_rnames_2 <- c("Transgenerational Effects")

Plasticity_k_2 <- data.frame("k" = c(Plasticity_Exploration["Transgenerational", "Freq"]), 
                             row.names = Plasticity_rnames_2)

Plasticity_group_no_2 <- data.frame("Spp No." = c(Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                    row.names = Plasticity_rnames_2)

Plasticity_study_2 <- data.frame("Study" = c(Plasticity_Study_Count["Transgenerational", "Freq"]), 
                                 row.names = Plasticity_rnames_2)

overall_plastic_summary_Reorder_2 <- overall_plastic_summary[c("Transgenerational Effects"), ]

Plasticity_table_2 <- data.frame(estimate = overall_plastic_summary_Reorder_2[,"Mean"], 
                                 lowerCL = overall_plastic_summary_Reorder_2[,"Low_CI"], 
                                 upperCL = overall_plastic_summary_Reorder_2[,"High_CI"], 
                                 K = Plasticity_k_2[,1], 
                                 group_no = Plasticity_group_no_2[,1], 
                                 row.names = Plasticity_rnames_2)
Plasticity_table_2$name <- row.names(Plasticity_table_2)

Plasticity_raw_mean_2 <- c(b_abs_plastic_trans)

Plasticity_raw_name_2 <- c(replicate(8000, "Transgenerational Effects"))

Plasticity_raw_df_2 <- data.frame("Model" = Plasticity_raw_name_2, 
                                  "Effect" = Plasticity_raw_mean_2)

# Graph code - Part 2

Plasticity_Order_2 <- c("Transgenerational Effects")

density_plasticity_2 <- Plasticity_table_2 %>% mutate(name = fct_relevel(name, Plasticity_Order_2)) %>%
                        ggplot() +
                        geom_density_ridges(data = Plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Plasticity_Order_2)), 
                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                            scale = 0.02, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                       size = 1) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_2)[1], 1)), xmin = max(Plasticity_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                       size = 1) +
                        geom_linerange(aes(y = rev(seq(1, dim(Plasticity_table_2)[1], 1)), xmin = min(Plasticity_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                       size = 1) +
                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                        size = 1, fatten = 2) +
                        theme_bw() +
                        guides(fill = "none", colour = "none") +
                        labs(x = TeX("Effect Size (PRRD)"), y = "") +
                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                         vjust = c(-0.8))) +
                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                        theme(axis.ticks = element_blank()) +
                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                        scale_colour_manual(values = c("#5D7AA1")) +
                        scale_fill_manual(values = c("#5D7AA1")) +
                        coord_cartesian(xlim = c(-0.01, 1.25)) +
                        annotate('text',  x = 1.25, y = (seq(1, dim(Plasticity_table_2)[1], 1)+0.4),
                        label= paste("italic(k)==", c(Plasticity_table_2["Transgenerational Effects", "K"]), "~","(", 
                                                    c(Plasticity_table_2["Transgenerational Effects", "group_no"]), 
                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                        geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_plastic_trans)-1)*100, 2), nsmall = 2), "%")), 
                                       x = rev(Plasticity_table_2$estimate+0.2), y = (seq(1, dim(Plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_plasticity_2 #(400x160)

##### Overall Model - Trait Category Meta-regression #####
Trait_Exploration <- data %>% select("Category") %>% table() %>% data.frame()
rownames(Trait_Exploration) <- Trait_Exploration$Category

Trait_Species_Count <- data %>% select("Scientific_Name", "Category") %>% 
                       table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                       select("Category") %>% table() %>% data.frame()
rownames(Trait_Species_Count) <- Trait_Species_Count$Category

Trait_Study_Count <- data %>% select("Study_ID", "Category") %>% 
                     table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                     select("Category") %>% table() %>% data.frame()
rownames(Trait_Study_Count) <- Trait_Study_Count$Category

system.time(
  overall_trait <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                             ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Category, cor = FALSE)),
                             data = data,
                             family = gaussian(),
                             data2 = list(A = A), 
                             chains = 4, 
                             cores = 4,
                             iter = 12000,
                             warmup = 2000,
                             thin = 5,
                             prior = priors,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             file = "./overall_trait_model",
                             file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_trait <- as_draws_df(overall_trait, variable = c("b_Intercept", "b_CategoryBiochemicalAssay", 
                                                   "b_CategoryGeneExpression", "b_CategoryLifeMHistoryTraits", 
                                                   "b_CategoryMorphology", "b_CategoryPhysiological", 
                                                   "b_CategoryTolerance"))
b_trait <- data.frame("b_Behavioural" = b_trait$b_Intercept, 
                      "b_BiochemicalAssay" = b_trait$b_CategoryBiochemicalAssay + b_trait$b_Intercept, 
                      "b_GeneExpression" = b_trait$b_CategoryGeneExpression + b_trait$b_Intercept, 
                      "b_LifeHistory" = b_trait$b_CategoryLifeMHistoryTraits + b_trait$b_Intercept, 
                      "b_Morphology" = b_trait$b_CategoryMorphology + b_trait$b_Intercept, 
                      "b_Physiological" = b_trait$b_CategoryPhysiological + b_trait$b_Intercept, 
                      "b_Tolerance" = b_trait$b_CategoryTolerance + b_trait$b_Intercept)


sd_trait <- as_draws_df(overall_trait, variable = c("sd_obs__Intercept:CategoryBehavioural", "sd_obs__Intercept:CategoryBiochemicalAssay", 
                                                    "sd_obs__Intercept:CategoryGeneExpression", "sd_obs__Intercept:CategoryLife-HistoryTraits", 
                                                    "sd_obs__Intercept:CategoryMorphology", "sd_obs__Intercept:CategoryPhysiological", 
                                                    "sd_obs__Intercept:CategoryTolerance", "sd_phylo__Intercept", 
                                                    "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_trait <- data.frame("sd_Behavioural" = sd_trait$`sd_obs__Intercept:CategoryBehavioural`,
                       "sd_BiochemicalAssay" = sd_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`, 
                       "sd_GeneExpression" = sd_trait$`sd_obs__Intercept:CategoryGeneExpression`, 
                       "sd_LifeHistory" = sd_trait$`sd_obs__Intercept:CategoryLife-HistoryTraits`, 
                       "sd_Morphology" = sd_trait$`sd_obs__Intercept:CategoryMorphology`, 
                       "sd_Physiological" = sd_trait$`sd_obs__Intercept:CategoryPhysiological`, 
                       "sd_Tolerance" = sd_trait$`sd_obs__Intercept:CategoryTolerance`,
                       "sd_phylo__Intercept" = sd_trait$`sd_phylo__Intercept`, 
                       "sd_Study_ID__Intercept" = sd_trait$`sd_Study_ID__Intercept`, 
                       "sd_Measurement__Intercept" = sd_trait$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
overall_trait_means <- apply(b_trait, 2, mean)
overall_trait_cis <- apply(b_trait, 2, function(x) HPDinterval(as.mcmc(x)))
overall_trait_pMCMC <- apply(b_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_trait_behavioural <- folded_norm(b_trait$b_Behavioural, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_Behavioural"]^2)))
b_abs_trait_biochem <- folded_norm(b_trait$b_BiochemicalAssay, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_trait_gene <- folded_norm(b_trait$b_GeneExpression, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_GeneExpression"]^2)))
b_abs_trait_life <- folded_norm(b_trait$b_LifeHistory, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_LifeHistory"]^2)))
b_abs_trait_morphology <- folded_norm(b_trait$b_Morphology, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_Morphology"]^2)))
b_abs_trait_physiological <- folded_norm(b_trait$b_Physiological, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_Physiological"]^2)))
b_abs_trait_tolerance <- folded_norm(b_trait$b_Tolerance, sqrt(rowSums(sd_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_trait[, "sd_Tolerance"]^2)))
mean_abs_b_trait_behavioural <- mean(b_abs_trait_behavioural)
mean_abs_b_trait_biochem <- mean(b_abs_trait_biochem)
mean_abs_b_trait_gene <- mean(b_abs_trait_gene)
mean_abs_b_trait_life <- mean(b_abs_trait_life)
mean_abs_b_trait_morphology <- mean(b_abs_trait_morphology)
mean_abs_b_trait_physiological <- mean(b_abs_trait_physiological)
mean_abs_b_trait_tolerance <- mean(b_abs_trait_tolerance)
ci.abs_behavioural <- HPDinterval(as.mcmc(b_abs_trait_behavioural))
ci.abs_biochem <- HPDinterval(as.mcmc(b_abs_trait_biochem))
ci.abs_gene <- HPDinterval(as.mcmc(b_abs_trait_gene))
ci.abs_life <- HPDinterval(as.mcmc(b_abs_trait_life))
ci.abs_morphology <- HPDinterval(as.mcmc(b_abs_trait_morphology))
ci.abs_physiological <- HPDinterval(as.mcmc(b_abs_trait_physiological))
ci.abs_tolerance <- HPDinterval(as.mcmc(b_abs_trait_tolerance))

# Heterogeneity
overall_i2_trait <- i2(sd_trait, data$Variance_Adjusted) 

# Overall_trait_summary

overall_trait_means_list <- c(mean_abs_b_trait_behavioural, mean_abs_b_trait_biochem, mean_abs_b_trait_gene, mean_abs_b_trait_life,
                              mean_abs_b_trait_morphology, mean_abs_b_trait_physiological, mean_abs_b_trait_tolerance)
overall_trait_low_ci <- c(ci.abs_behavioural[1], ci.abs_biochem[1], ci.abs_gene[1], ci.abs_life[1], ci.abs_morphology[1], ci.abs_physiological[1], 
                          ci.abs_tolerance[1])
overall_trait_high_ci <- c(ci.abs_behavioural[2], ci.abs_biochem[2], ci.abs_gene[2], ci.abs_life[2], ci.abs_morphology[2], ci.abs_physiological[2],
                           ci.abs_tolerance[2])
overall_trait_categories <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-History Traits",
                              "Morphology", "Physiological", "Tolerance")

overall_trait_summary <- matrix(c(overall_trait_means_list, overall_trait_low_ci, overall_trait_high_ci), 
                                nrow = 7, ncol = 3, byrow = FALSE, 
                                dimnames = list(c(overall_trait_categories), 
                                                c("Mean", "Low_CI", "High_CI")))
overall_trait_summary <- data.frame(overall_trait_summary)

# Preparing Graph - Combined

Trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits", 
                  "Morphological", "Physiological", "Tolerance")

Trait_k <- data.frame("k" = c(Trait_Exploration["Behavioural", "Freq"], 
                              Trait_Exploration["Biochemical Assay", "Freq"], 
                              Trait_Exploration["Gene Expression", "Freq"], 
                              Trait_Exploration["Life-History Traits", "Freq"], 
                              Trait_Exploration["Morphology", "Freq"], 
                              Trait_Exploration["Physiological", "Freq"], 
                              Trait_Exploration["Tolerance", "Freq"]), 
                      row.names = Trait_rnames)

Trait_group_no <- data.frame("Spp No." = c(Trait_Species_Count["Behavioural", "Freq"], 
                                           Trait_Species_Count["Biochemical Assay", "Freq"], 
                                           Trait_Species_Count["Gene Expression", "Freq"], 
                                           Trait_Species_Count["Life-History Traits", "Freq"], 
                                           Trait_Species_Count["Morphology", "Freq"], 
                                           Trait_Species_Count["Physiological", "Freq"], 
                                           Trait_Species_Count["Tolerance", "Freq"]), 
                             row.names = Trait_rnames)

Trait_study <- data.frame("Study" = c(Trait_Study_Count["Behavioural", "Freq"], 
                                      Trait_Study_Count["Biochemical Assay", "Freq"], 
                                      Trait_Study_Count["Gene Expression", "Freq"], 
                                      Trait_Study_Count["Life-History Traits", "Freq"], 
                                      Trait_Study_Count["Morphology", "Freq"], 
                                      Trait_Study_Count["Physiological", "Freq"], 
                                      Trait_Study_Count["Tolerance", "Freq"]), 
                          row.names = Trait_rnames)

overall_trait_summary_Reorder <- overall_trait_summary[c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-History Traits", 
                                                         "Morphology", "Physiological", "Tolerance"), ]

Trait_table <- data.frame(estimate = overall_trait_summary_Reorder[,"Mean"], 
                          lowerCL = overall_trait_summary_Reorder[,"Low_CI"], 
                          upperCL = overall_trait_summary_Reorder[,"High_CI"], 
                          K = Trait_k[,1], 
                          group_no = Trait_group_no[,1], 
                          row.names = Trait_rnames)
Trait_table$name <- row.names(Trait_table)

Trait_raw_mean <- c(b_abs_trait_behavioural, b_abs_trait_biochem, b_abs_trait_gene, 
                    b_abs_trait_life, b_abs_trait_morphology, b_abs_trait_physiological, 
                    b_abs_trait_tolerance)

Trait_raw_name <- c(replicate(8000, "Behavioural"), 
                    replicate(8000, "Biochemical Assay"), 
                    replicate(8000, "Gene Expression"), 
                    replicate(8000, "Life-history Traits"), 
                    replicate(8000, "Morphological"), 
                    replicate(8000, "Physiological"), 
                    replicate(8000, "Tolerance"))

Trait_raw_df <- data.frame("Model" = Trait_raw_name, 
                           "Effect" = Trait_raw_mean)

# Graph code - Combined

Trait_Order <- c("Tolerance", "Physiological", "Morphological", "Life-history Traits", 
                 "Gene Expression", "Biochemical Assay", "Behavioural")

density_trait <- Trait_table %>% mutate(name = fct_relevel(name, Trait_Order)) %>%
                 ggplot() +
                 geom_density_ridges(data = Trait_raw_df %>% mutate(Model = fct_relevel(Model, Trait_Order)), 
                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                 geom_linerange(aes(y = rev(seq(1, dim(Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(Trait_table)[1], 1)), xmin = max(Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(Trait_table)[1], 1)), xmin = min(Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                size = 1) +
                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                 size = 1, fatten = 2) +
                 theme_bw() +
                 guides(fill = "none", colour = "none") +
                 labs(x = TeX("Effect Size (PRRD)"), y = "") +
                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                  vjust = c(-2.7, -2.7, -2.7, -0.8, -0.8, -0.8, -2.7))) +
                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                 theme(axis.ticks = element_blank()) +
                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                 scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                 scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                 coord_cartesian(xlim = c(-0.01, 1.25)) +
                 annotate('text',  x = 1.25, y = (seq(1, dim(Trait_table)[1], 1)+0.4),
                 label= paste("italic(k)==", c(Trait_table["Tolerance", "K"], 
                                               Trait_table["Physiological", "K"], 
                                               Trait_table["Morphological", "K"], 
                                               Trait_table["Life-history Traits", "K"], 
                                               Trait_table["Gene Expression", "K"], 
                                               Trait_table["Biochemical Assay", "K"], 
                                               Trait_table["Behavioural", "K"]), "~","(", 
                                             c(Trait_table["Tolerance", "group_no"], 
                                               Trait_table["Physiological", "group_no"], 
                                               Trait_table["Morphological", "group_no"], 
                                               Trait_table["Life-history Traits", "group_no"], 
                                               Trait_table["Gene Expression", "group_no"], 
                                               Trait_table["Biochemical Assay", "group_no"], 
                                               Trait_table["Behavioural", "group_no"]), 
                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                 geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(b_abs_trait_physiological)-1)*100, 2), nsmall = 2), "%"),
                                        paste(format(round(mean(exp(b_abs_trait_morphology)-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(b_abs_trait_life)-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(b_abs_trait_gene)-1)*100, 2), nsmall = 2), "%"),
                                        paste(format(round(mean(exp(b_abs_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(b_abs_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                            x = rev(Trait_table$estimate+c(0.2, 0.2, 0.2, 0.2, 0.2, -0.2, 0.2)), y = (seq(1, dim(Trait_table)[1], 1)+0.4)), size = 3.5)

density_trait #(400x560)

# Preparing Graph - Part 1

Trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits")

Trait_k_1 <- data.frame("k" = c(Trait_Exploration["Behavioural", "Freq"], 
                                Trait_Exploration["Biochemical Assay", "Freq"], 
                                Trait_Exploration["Gene Expression", "Freq"], 
                                Trait_Exploration["Life-History Traits", "Freq"]), 
                        row.names = Trait_rnames_1)

Trait_group_no_1 <- data.frame("Spp No." = c(Trait_Species_Count["Behavioural", "Freq"], 
                                             Trait_Species_Count["Biochemical Assay", "Freq"], 
                                             Trait_Species_Count["Gene Expression", "Freq"], 
                                             Trait_Species_Count["Life-History Traits", "Freq"]), 
                               row.names = Trait_rnames_1)

Trait_study_1 <- data.frame("Study" = c(Trait_Study_Count["Behavioural", "Freq"], 
                                        Trait_Study_Count["Biochemical Assay", "Freq"], 
                                        Trait_Study_Count["Gene Expression", "Freq"], 
                                        Trait_Study_Count["Life-History Traits", "Freq"]), 
                            row.names = Trait_rnames_1)

overall_trait_summary_Reorder_1 <- overall_trait_summary[c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-History Traits"), ]

Trait_table_1 <- data.frame(estimate = overall_trait_summary_Reorder_1[,"Mean"], 
                            lowerCL = overall_trait_summary_Reorder_1[,"Low_CI"], 
                            upperCL = overall_trait_summary_Reorder_1[,"High_CI"], 
                            K = Trait_k_1[,1], 
                            group_no = Trait_group_no_1[,1], 
                            row.names = Trait_rnames_1)
Trait_table_1$name <- row.names(Trait_table_1)

Trait_raw_mean_1 <- c(b_abs_trait_behavioural, b_abs_trait_biochem, b_abs_trait_gene, b_abs_trait_life)

Trait_raw_name_1 <- c(replicate(8000, "Behavioural"), 
                      replicate(8000, "Biochemical Assay"), 
                      replicate(8000, "Gene Expression"), 
                      replicate(8000, "Life-history Traits"))

Trait_raw_df_1 <- data.frame("Model" = Trait_raw_name_1, 
                             "Effect" = Trait_raw_mean_1)

# Graph code - Part 1

Trait_Order_1 <- c("Life-history Traits", "Gene Expression", "Biochemical Assay", "Behavioural")

density_trait_1 <- Trait_table_1 %>% mutate(name = fct_relevel(name, Trait_Order_1)) %>%
                   ggplot() +
                   geom_density_ridges(data = Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Trait_Order_1)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_1)[1], 1)), xmin = max(Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_1)[1], 1)), xmin = min(Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                   size = 1, fatten = 2) +
                   theme_bw() +
                   guides(fill = "none", colour = "none") +
                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                    vjust = c(-0.8, -0.8, -0.8, -2.7))) +
                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                   theme(axis.ticks = element_blank()) +
                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                   coord_cartesian(xlim = c(-0.01, 1.25)) +
                   annotate('text',  x = 1.25, y = (seq(1, dim(Trait_table_1)[1], 1)+0.4),
                   label= paste("italic(k)==", c(Trait_table_1["Life-history Traits", "K"], 
                                                 Trait_table_1["Gene Expression", "K"], 
                                                 Trait_table_1["Biochemical Assay", "K"], 
                                                 Trait_table_1["Behavioural", "K"]), "~","(", 
                                               c(Trait_table_1["Life-history Traits", "group_no"], 
                                                 Trait_table_1["Gene Expression", "group_no"], 
                                                 Trait_table_1["Biochemical Assay", "group_no"], 
                                                 Trait_table_1["Behavioural", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_trait_life)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_trait_gene)-1)*100, 2), nsmall = 2), "%"),
                                           paste(format(round(mean(exp(b_abs_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                 x = rev(Trait_table_1$estimate+c(0.2, 0.2, 0.2, 0.2)), y = (seq(1, dim(Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_trait_1 #(400x400)

# Preparing Graph - Part 2

Trait_rnames_2 <- c("Morphological", "Physiological", "Tolerance")

Trait_k_2 <- data.frame("k" = c(Trait_Exploration["Morphology", "Freq"], 
                                Trait_Exploration["Physiological", "Freq"], 
                                Trait_Exploration["Tolerance", "Freq"]), 
                        row.names = Trait_rnames_2)

Trait_group_no_2 <- data.frame("Spp No." = c(Trait_Species_Count["Morphology", "Freq"], 
                                             Trait_Species_Count["Physiological", "Freq"], 
                                             Trait_Species_Count["Tolerance", "Freq"]), 
                               row.names = Trait_rnames_2)

Trait_study_2 <- data.frame("Study" = c(Trait_Study_Count["Morphology", "Freq"], 
                                        Trait_Study_Count["Physiological", "Freq"], 
                                        Trait_Study_Count["Tolerance", "Freq"]), 
                            row.names = Trait_rnames_2)

overall_trait_summary_Reorder_2 <- overall_trait_summary[c("Morphology", "Physiological", "Tolerance"), ]

Trait_table_2 <- data.frame(estimate = overall_trait_summary_Reorder_2[,"Mean"], 
                            lowerCL = overall_trait_summary_Reorder_2[,"Low_CI"], 
                            upperCL = overall_trait_summary_Reorder_2[,"High_CI"], 
                            K = Trait_k_2[,1], 
                            group_no = Trait_group_no_2[,1], 
                            row.names = Trait_rnames_2)
Trait_table_2$name <- row.names(Trait_table_2)

Trait_raw_mean_2 <- c(b_abs_trait_morphology, b_abs_trait_physiological, b_abs_trait_tolerance)

Trait_raw_name_2 <- c(replicate(8000, "Morphological"), 
                      replicate(8000, "Physiological"), 
                      replicate(8000, "Tolerance"))

Trait_raw_df_2 <- data.frame("Model" = Trait_raw_name_2, 
                             "Effect" = Trait_raw_mean_2)

# Graph code - Part 2

Trait_Order_2 <- c("Tolerance", "Physiological", "Morphological")

density_trait_2 <- Trait_table_2 %>% mutate(name = fct_relevel(name, Trait_Order_2)) %>%
                   ggplot() +
                   geom_density_ridges(data = Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Trait_Order_2)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_2)[1], 1)), xmin = max(Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(Trait_table_2)[1], 1)), xmin = min(Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                   size = 1, fatten = 2) +
                   theme_bw() +
                   guides(fill = "none", colour = "none") +
                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                    vjust = c(-2.7, -2.7, -2.7))) +
                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                   theme(axis.ticks = element_blank()) +
                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                   scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                   scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                   coord_cartesian(xlim = c(-0.01, 1.25)) +
                   annotate('text',  x = 1.25, y = (seq(1, dim(Trait_table_2)[1], 1)+0.4),
                   label= paste("italic(k)==", c(Trait_table_2["Tolerance", "K"], 
                                                 Trait_table_2["Physiological", "K"], 
                                                 Trait_table_2["Morphological", "K"]), "~","(", 
                                               c(Trait_table_2["Tolerance", "group_no"], 
                                                 Trait_table_2["Physiological", "group_no"], 
                                                 Trait_table_2["Morphological", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_trait_physiological)-1)*100, 2), nsmall = 2), "%"),
                                           paste(format(round(mean(exp(b_abs_trait_morphology)-1)*100, 2), nsmall = 2), "%")), 
                                 x = rev(Trait_table_2$estimate+c(0.2, -0.2, 0.2)), y = (seq(1, dim(Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_trait_2 #(400x320)

##### Overall Model - Taxonomy Category Meta-regression #####
Taxonomy_Exploration <- data %>% select("Class") %>% table() %>% data.frame()
Taxonomy_Exploration <- Taxonomy_Exploration %>% filter(Freq > 10)
rownames(Taxonomy_Exploration) <- Taxonomy_Exploration$Class

Taxonomy_Data <- data %>% filter(Class == "Actinopteri"| 
                                 Class == "Amphibia"| 
                                 Class == "Branchiopoda"|
                                 Class == "Gastropoda"| 
                                 Class == "Insecta")

Taxonomy_Species_Count <- Taxonomy_Data %>% select("Scientific_Name", "Class") %>% 
                          table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                          select("Class") %>% table() %>% data.frame()
rownames(Taxonomy_Species_Count) <- Taxonomy_Species_Count$Class

Taxonomy_Study_Count <- Taxonomy_Data %>% select("Study_ID", "Class") %>% 
                        table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                        select("Class") %>% table() %>% data.frame()
rownames(Taxonomy_Study_Count) <- Taxonomy_Study_Count$Class

Taxonomy_Species <- Taxonomy_Data %>% select("phylo") %>% unique()

Taxonomy_A <- as.data.frame(A)
Taxonomy_A <- Taxonomy_A[c(Taxonomy_Species$phylo), c(Taxonomy_Species$phylo)]
Taxonomy_A <- as.matrix(Taxonomy_A)

system.time(
  overall_taxonomy <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                                ~ Class + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Class, cor = FALSE)),
                                data = Taxonomy_Data,
                                family = gaussian(),
                                data2 = list(A = Taxonomy_A), 
                                chains = 4, 
                                cores = 4,
                                iter = 12000,
                                warmup = 2000,
                                thin = 5,
                                prior = priors,
                                control = list(adapt_delta = 0.99, max_treedepth = 15),
                                file = "./overall_taxonomy_model",
                                file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_taxonomy <- as_draws_df(overall_taxonomy, variable = c("b_Intercept", "b_ClassAmphibia", 
                                                         "b_ClassBranchiopoda", "b_ClassGastropoda", 
                                                         "b_ClassInsecta"))
b_taxonomy <- data.frame("b_Actinopteri" = b_taxonomy$b_Intercept, 
                         "b_Amphibia" = b_taxonomy$b_ClassAmphibia + b_taxonomy$b_Intercept, 
                         "b_Branchiopoda" = b_taxonomy$b_ClassBranchiopoda + b_taxonomy$b_Intercept, 
                         "b_Gastropoda" = b_taxonomy$b_ClassGastropoda + b_taxonomy$b_Intercept, 
                         "b_Insecta" = b_taxonomy$b_ClassInsecta + b_taxonomy$b_Intercept)


sd_taxonomy <- as_draws_df(overall_taxonomy, variable = c("sd_obs__Intercept:ClassActinopteri", "sd_obs__Intercept:ClassAmphibia", 
                                                          "sd_obs__Intercept:ClassBranchiopoda", "sd_obs__Intercept:ClassGastropoda", 
                                                          "sd_obs__Intercept:ClassInsecta", "sd_phylo__Intercept", 
                                                          "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_taxonomy <- data.frame("sd_Actinopteri" = sd_taxonomy$`sd_obs__Intercept:ClassActinopteri`,
                          "sd_Amphibia" = sd_taxonomy$`sd_obs__Intercept:ClassAmphibia`, 
                          "sd_Branchiopoda" = sd_taxonomy$`sd_obs__Intercept:ClassBranchiopoda`, 
                          "sd_Gastropoda" = sd_taxonomy$`sd_obs__Intercept:ClassGastropoda`, 
                          "sd_Insecta" = sd_taxonomy$`sd_obs__Intercept:ClassInsecta`, 
                          "sd_phylo__Intercept" = sd_taxonomy$`sd_phylo__Intercept`, 
                          "sd_Study_ID__Intercept" = sd_taxonomy$`sd_Study_ID__Intercept`, 
                          "sd_Measurement__Intercept" = sd_taxonomy$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
overall_taxonomy_means <- apply(b_taxonomy, 2, mean)
overall_taxonomy_cis <- apply(b_taxonomy, 2, function(x) HPDinterval(as.mcmc(x)))
overall_taxonomy_pMCMC <- apply(b_taxonomy, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_taxonomy_actinopteri <- folded_norm(b_taxonomy$b_Actinopteri, sqrt(rowSums(sd_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_taxonomy[, "sd_Actinopteri"]^2)))
b_abs_taxonomy_amphibia <- folded_norm(b_taxonomy$b_Amphibia, sqrt(rowSums(sd_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_taxonomy[, "sd_Amphibia"]^2)))
b_abs_taxonomy_branchiopoda <- folded_norm(b_taxonomy$b_Branchiopoda, sqrt(rowSums(sd_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_taxonomy[, "sd_Branchiopoda"]^2)))
b_abs_taxonomy_gastropoda <- folded_norm(b_taxonomy$b_Gastropoda, sqrt(rowSums(sd_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_taxonomy[, "sd_Gastropoda"]^2)))
b_abs_taxonomy_insecta <- folded_norm(b_taxonomy$b_Insecta, sqrt(rowSums(sd_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_taxonomy[, "sd_Insecta"]^2)))
mean_abs_b_taxonomy_actinopteri <- mean(b_abs_taxonomy_actinopteri)
mean_abs_b_taxonomy_amphibia <- mean(b_abs_taxonomy_amphibia)
mean_abs_b_taxonomy_branchiopoda <- mean(b_abs_taxonomy_branchiopoda)
mean_abs_b_taxonomy_gastropoda <- mean(b_abs_taxonomy_gastropoda)
mean_abs_b_taxonomy_insecta <- mean(b_abs_taxonomy_insecta)
ci.abs_taxonomy_actinopteri <- HPDinterval(as.mcmc(b_abs_taxonomy_actinopteri))
ci.abs_taxonomy_amphibia <- HPDinterval(as.mcmc(b_abs_taxonomy_amphibia))
ci.abs_taxonomy_branchiopoda <- HPDinterval(as.mcmc(b_abs_taxonomy_branchiopoda))
ci.abs_taxonomy_gastropoda <- HPDinterval(as.mcmc(b_abs_taxonomy_gastropoda))
ci.abs_taxonomy_insecta <- HPDinterval(as.mcmc(b_abs_taxonomy_insecta))

# Heterogeneity
overall_i2_taxonomy <- i2(sd_taxonomy, Taxonomy_Data$Variance_Adjusted) 

# Overall_trait_summary
overall_taxonomy_means_list <- c(mean_abs_b_taxonomy_actinopteri, mean_abs_b_taxonomy_amphibia, mean_abs_b_taxonomy_branchiopoda, 
                                 mean_abs_b_taxonomy_gastropoda, mean_abs_b_taxonomy_insecta)
overall_taxonomy_low_ci <- c(ci.abs_taxonomy_actinopteri[1], ci.abs_taxonomy_amphibia[1], ci.abs_taxonomy_branchiopoda[1], 
                             ci.abs_taxonomy_gastropoda[1], ci.abs_taxonomy_insecta[1])
overall_taxonomy_high_ci <- c(ci.abs_taxonomy_actinopteri[2], ci.abs_taxonomy_amphibia[2], ci.abs_taxonomy_branchiopoda[2], 
                              ci.abs_taxonomy_gastropoda[2], ci.abs_taxonomy_insecta[2])
overall_taxonomy_categories <- c("Actinopteri", "Amphibia", "Branchiopoda", 
                                 "Gastropoda", "Insecta")

overall_taxonomy_summary <- matrix(c(overall_taxonomy_means_list, overall_taxonomy_low_ci, overall_taxonomy_high_ci), 
                                     nrow = 5, ncol = 3, byrow = FALSE, 
                                     dimnames = list(c(overall_taxonomy_categories), 
                                                     c("Mean", "Low_CI", "High_CI")))
overall_taxonomy_summary <- data.frame(overall_taxonomy_summary)

# Preparing Graph - Combined

Taxonomy_rnames <- c("Actinopteri", "Amphibia", "Branchiopoda", 
                     "Gastropoda", "Insecta")

Taxonomy_k <- data.frame("k" = c(Taxonomy_Exploration["Actinopteri", "Freq"], 
                                 Taxonomy_Exploration["Amphibia", "Freq"], 
                                 Taxonomy_Exploration["Branchiopoda", "Freq"], 
                                 Taxonomy_Exploration["Gastropoda", "Freq"], 
                                 Taxonomy_Exploration["Insecta", "Freq"]), 
                         row.names = Taxonomy_rnames)

Taxonomy_group_no <- data.frame("Spp No." = c(Taxonomy_Species_Count["Actinopteri", "Freq"], 
                                              Taxonomy_Species_Count["Amphibia", "Freq"], 
                                              Taxonomy_Species_Count["Branchiopoda", "Freq"], 
                                              Taxonomy_Species_Count["Gastropoda", "Freq"], 
                                              Taxonomy_Species_Count["Insecta", "Freq"]), 
                                row.names = Taxonomy_rnames)

Taxonomy_study <- data.frame("Study" = c(Taxonomy_Study_Count["Actinopteri", "Freq"], 
                                         Taxonomy_Study_Count["Amphibia", "Freq"], 
                                         Taxonomy_Study_Count["Branchiopoda", "Freq"], 
                                         Taxonomy_Study_Count["Gastropoda", "Freq"], 
                                         Taxonomy_Study_Count["Insecta", "Freq"]), 
                             row.names = Taxonomy_rnames)

Taxonomy_table <- data.frame(estimate = overall_taxonomy_summary[,"Mean"], 
                             lowerCL = overall_taxonomy_summary[,"Low_CI"], 
                             upperCL = overall_taxonomy_summary[,"High_CI"], 
                             K = Taxonomy_k[,1], 
                             group_no = Taxonomy_group_no[,1], 
                             row.names = Taxonomy_rnames)
Taxonomy_table$name <- row.names(Taxonomy_table)

Taxonomy_raw_mean <- c(b_abs_taxonomy_actinopteri, b_abs_taxonomy_amphibia, b_abs_taxonomy_branchiopoda, 
                       b_abs_taxonomy_gastropoda, b_abs_taxonomy_insecta)

Taxonomy_raw_name <- c(replicate(8000, "Actinopteri"), 
                       replicate(8000, "Amphibia"), 
                       replicate(8000, "Branchiopoda"), 
                       replicate(8000, "Gastropoda"), 
                       replicate(8000, "Insecta"))

Taxonomy_raw_df <- data.frame("Model" = Taxonomy_raw_name, 
                              "Effect" = Taxonomy_raw_mean)

# Graph code - Combined

Taxonomy_Order <- c("Insecta", "Gastropoda", "Branchiopoda", 
                    "Amphibia", "Actinopteri")

density_taxonomy <- Taxonomy_table %>% mutate(name = fct_relevel(name, Taxonomy_Order)) %>%
                    ggplot() +
                    geom_density_ridges(data = Taxonomy_raw_df %>% mutate(Model = fct_relevel(Model, Taxonomy_Order)), 
                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                        scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                    geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                   size = 1) +
                    geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table)[1], 1)), xmin = max(Taxonomy_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                   size = 1) +
                    geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table)[1], 1)), xmin = min(Taxonomy_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                   size = 1) +
                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Taxonomy_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                    size = 1, fatten = 2) +
                    theme_bw() +
                    guides(fill = "none", colour = "none") +
                    labs(x = TeX("Effect Size (PRRD)"), y = "") +
                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                     vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                    theme(axis.ticks = element_blank()) +
                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                    scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                    scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                    coord_cartesian(xlim = c(-0.01, 1.25)) +
                    annotate('text',  x = 1.25, y = (seq(1, dim(Taxonomy_table)[1], 1)+0.4),
                    label= paste("italic(k)==", c(Taxonomy_table["Insecta", "K"], 
                                                  Taxonomy_table["Gastropoda", "K"], 
                                                  Taxonomy_table["Branchiopoda", "K"], 
                                                  Taxonomy_table["Amphibia", "K"], 
                                                  Taxonomy_table["Actinopteri", "K"]), "~","(", 
                                                c(Taxonomy_table["Insecta", "group_no"], 
                                                  Taxonomy_table["Gastropoda", "group_no"], 
                                                  Taxonomy_table["Branchiopoda", "group_no"], 
                                                  Taxonomy_table["Amphibia", "group_no"], 
                                                  Taxonomy_table["Actinopteri", "group_no"]), 
                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_taxonomy_insecta)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_taxonomy_gastropoda)-1)*100, 2), nsmall = 2), "%"),
                                           paste(format(round(mean(exp(b_abs_taxonomy_branchiopoda)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_taxonomy_amphibia)-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(b_abs_taxonomy_actinopteri)-1)*100, 2), nsmall = 2), "%")), 
                                   x = rev(Taxonomy_table$estimate+0.2), y = (seq(1, dim(Taxonomy_table)[1], 1)+0.4)), size = 3.5)

density_taxonomy #(400x480)

# Preparing Graph - Part 1

Taxonomy_rnames_1 <- c("Actinopteri", "Amphibia", "Branchiopoda")

Taxonomy_k_1 <- data.frame("k" = c(Taxonomy_Exploration["Actinopteri", "Freq"], 
                                   Taxonomy_Exploration["Amphibia", "Freq"], 
                                   Taxonomy_Exploration["Branchiopoda", "Freq"]), 
                           row.names = Taxonomy_rnames_1)

Taxonomy_group_no_1 <- data.frame("Spp No." = c(Taxonomy_Species_Count["Actinopteri", "Freq"], 
                                                Taxonomy_Species_Count["Amphibia", "Freq"], 
                                                Taxonomy_Species_Count["Branchiopoda", "Freq"]), 
                                  row.names = Taxonomy_rnames_1)

Taxonomy_study_1 <- data.frame("Study" = c(Taxonomy_Study_Count["Actinopteri", "Freq"], 
                                           Taxonomy_Study_Count["Amphibia", "Freq"], 
                                           Taxonomy_Study_Count["Branchiopoda", "Freq"]), 
                               row.names = Taxonomy_rnames_1)

overall_taxonomy_summary_Reorder_1 <- overall_taxonomy_summary[c("Actinopteri", "Amphibia", "Branchiopoda"), ]

Taxonomy_table_1 <- data.frame(estimate = overall_taxonomy_summary_Reorder_1[,"Mean"], 
                               lowerCL = overall_taxonomy_summary_Reorder_1[,"Low_CI"], 
                               upperCL = overall_taxonomy_summary_Reorder_1[,"High_CI"], 
                               K = Taxonomy_k_1[,1], 
                               group_no = Taxonomy_group_no_1[,1], 
                               row.names = Taxonomy_rnames_1)
Taxonomy_table_1$name <- row.names(Taxonomy_table_1)

Taxonomy_raw_mean_1 <- c(b_abs_taxonomy_actinopteri, b_abs_taxonomy_amphibia, b_abs_taxonomy_branchiopoda)

Taxonomy_raw_name_1 <- c(replicate(8000, "Actinopteri"), 
                         replicate(8000, "Amphibia"), 
                         replicate(8000, "Branchiopoda"))

Taxonomy_raw_df_1 <- data.frame("Model" = Taxonomy_raw_name_1, 
                                "Effect" = Taxonomy_raw_mean_1)

# Graph code - Part 1

Taxonomy_Order_1 <- c("Branchiopoda", "Amphibia", "Actinopteri")

density_taxonomy_1 <- Taxonomy_table_1 %>% mutate(name = fct_relevel(name, Taxonomy_Order_1)) %>%
                      ggplot() +
                      geom_density_ridges(data = Taxonomy_raw_df_1 %>% mutate(Model = fct_relevel(Model, Taxonomy_Order_1)), 
                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_1)[1], 1)), xmin = max(Taxonomy_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_1)[1], 1)), xmin = min(Taxonomy_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                     size = 1) +
                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Taxonomy_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                      size = 1, fatten = 2) +
                      theme_bw() +
                      guides(fill = "none", colour = "none") +
                      labs(x = TeX("Effect Size (PRRD)"), y = "") +
                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                       vjust = c(-2.7, -2.7, -2.7))) +
                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                      theme(axis.ticks = element_blank()) +
                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                      scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                      scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                      coord_cartesian(xlim = c(-0.01, 1.25)) +
                      annotate('text',  x = 1.25, y = (seq(1, dim(Taxonomy_table_1)[1], 1)+0.4),
                      label= paste("italic(k)==", c(Taxonomy_table_1["Branchiopoda", "K"], 
                                                    Taxonomy_table_1["Amphibia", "K"], 
                                                    Taxonomy_table_1["Actinopteri", "K"]), "~","(", 
                                                  c(Taxonomy_table_1["Branchiopoda", "group_no"], 
                                                    Taxonomy_table_1["Amphibia", "group_no"], 
                                                    Taxonomy_table_1["Actinopteri", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                      geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_taxonomy_branchiopoda)-1)*100, 2), nsmall = 2), "%"), 
                                             paste(format(round(mean(exp(b_abs_taxonomy_amphibia)-1)*100, 2), nsmall = 2), "%"), 
                                             paste(format(round(mean(exp(b_abs_taxonomy_actinopteri)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Taxonomy_table_1$estimate+0.2), y = (seq(1, dim(Taxonomy_table_1)[1], 1)+0.4)), size = 3.5)

density_taxonomy_1 #(400x320)

# Preparing Graph - Part 2

Taxonomy_rnames_2 <- c("Gastropoda", "Insecta")

Taxonomy_k_2 <- data.frame("k" = c(Taxonomy_Exploration["Gastropoda", "Freq"], 
                                   Taxonomy_Exploration["Insecta", "Freq"]), 
                           row.names = Taxonomy_rnames_2)

Taxonomy_group_no_2 <- data.frame("Spp No." = c(Taxonomy_Species_Count["Gastropoda", "Freq"], 
                                                Taxonomy_Species_Count["Insecta", "Freq"]), 
                                  row.names = Taxonomy_rnames_2)

Taxonomy_study_2 <- data.frame("Study" = c(Taxonomy_Study_Count["Gastropoda", "Freq"], 
                                           Taxonomy_Study_Count["Insecta", "Freq"]), 
                               row.names = Taxonomy_rnames_2)

overall_taxonomy_summary_Reorder_2 <- overall_taxonomy_summary[c("Gastropoda", "Insecta"), ]

Taxonomy_table_2 <- data.frame(estimate = overall_taxonomy_summary_Reorder_2[,"Mean"], 
                               lowerCL = overall_taxonomy_summary_Reorder_2[,"Low_CI"], 
                               upperCL = overall_taxonomy_summary_Reorder_2[,"High_CI"], 
                               K = Taxonomy_k_2[,1], 
                               group_no = Taxonomy_group_no_2[,1], 
                               row.names = Taxonomy_rnames_2)
Taxonomy_table_2$name <- row.names(Taxonomy_table_2)

Taxonomy_raw_mean_2 <- c(b_abs_taxonomy_gastropoda, b_abs_taxonomy_insecta)

Taxonomy_raw_name_2 <- c(replicate(8000, "Gastropoda"), 
                         replicate(8000, "Insecta"))

Taxonomy_raw_df_2 <- data.frame("Model" = Taxonomy_raw_name_2, 
                                "Effect" = Taxonomy_raw_mean_2)

# Graph code - Part 2

Taxonomy_Order_2 <- c("Insecta", "Gastropoda")

density_taxonomy_2 <- Taxonomy_table_2 %>% mutate(name = fct_relevel(name, Taxonomy_Order_2)) %>%
                      ggplot() +
                      geom_density_ridges(data = Taxonomy_raw_df_2 %>% mutate(Model = fct_relevel(Model, Taxonomy_Order_2)), 
                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_2)[1], 1)), xmin = max(Taxonomy_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                     size = 1) +
                      geom_linerange(aes(y = rev(seq(1, dim(Taxonomy_table_2)[1], 1)), xmin = min(Taxonomy_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                     size = 1) +
                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Taxonomy_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                      size = 1, fatten = 2) +
                      theme_bw() +
                      guides(fill = "none", colour = "none") +
                      labs(x = TeX("Effect Size (PRRD)"), y = "") +
                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                       vjust = c(-2.7, -2.7))) +
                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                      theme(axis.ticks = element_blank()) +
                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                      coord_cartesian(xlim = c(-0.01, 1.25)) +
                      annotate('text',  x = 1.25, y = (seq(1, dim(Taxonomy_table_2)[1], 1)+0.4),
                      label= paste("italic(k)==", c(Taxonomy_table_2["Insecta", "K"], 
                                                    Taxonomy_table_2["Gastropoda", "K"]), "~","(", 
                                                  c(Taxonomy_table_2["Insecta", "group_no"], 
                                                    Taxonomy_table_2["Gastropoda", "group_no"]), 
                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                      geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_taxonomy_insecta)-1)*100, 2), nsmall = 2), "%"), 
                                             paste(format(round(mean(exp(b_abs_taxonomy_gastropoda)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Taxonomy_table_2$estimate+0.2), y = (seq(1, dim(Taxonomy_table_2)[1], 1)+0.4)), size = 3.5)

density_taxonomy_2 #(400x240)

##### Overall Model - Treatment Meta-regression #####
Treatment_Exploration <- data %>% select("Type") %>% table() %>% data.frame()
Treatment_Exploration <- Treatment_Exploration %>% filter(Freq > 10)
rownames(Treatment_Exploration) <- Treatment_Exploration$Type

Treatment_Data <- data %>% filter(Type == "Diet"| 
                                  Type == "Humidity"| 
                                  Type == "Predator"|
                                  Type == "Salinity"| 
                                  Type == "Temperature"|
                                  Type == "Water Level")

Treatment_Species_Count <- Treatment_Data %>% select("Scientific_Name", "Type") %>% 
                           table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                           select("Type") %>% table() %>% data.frame()
rownames(Treatment_Species_Count) <- Treatment_Species_Count$Type

Treatment_Study_Count <- Treatment_Data %>% select("Study_ID", "Type") %>% 
                         table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                         select("Type") %>% table() %>% data.frame()
rownames(Treatment_Study_Count) <- Treatment_Study_Count$Type

Treatment_Species <- Treatment_Data %>% select("phylo") %>% unique()

Treatment_A <- as.data.frame(A)
Treatment_A <- Treatment_A[c(Treatment_Species$phylo), c(Treatment_Species$phylo)]
Treatment_A <- as.matrix(Treatment_A)

system.time(
  overall_treatment <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                             ~ Type + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Type, cor = FALSE)),
                             data = Treatment_Data,
                             family = gaussian(),
                             data2 = list(A = Treatment_A), 
                             chains = 4, 
                             cores = 4,
                             iter = 12000,
                             warmup = 2000,
                             thin = 5,
                             prior = priors,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             file = "./overall_treatment_model",
                             file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_treatment <- as_draws_df(overall_treatment, variable = c("b_Intercept", "b_TypeHumidity", 
                                                           "b_TypePredator", "b_TypeSalinity", 
                                                           "b_TypeTemperature", "b_TypeWaterLevel"))
b_treatment <- data.frame("b_Diet" = b_treatment$b_Intercept, 
                          "b_Humidity" = b_treatment$b_TypeHumidity + b_treatment$b_Intercept, 
                          "b_Predator" = b_treatment$b_TypePredator + b_treatment$b_Intercept, 
                          "b_Salinity" = b_treatment$b_TypeSalinity + b_treatment$b_Intercept, 
                          "b_Temperature" = b_treatment$b_TypeTemperature + b_treatment$b_Intercept, 
                          "b_WaterLevel" = b_treatment$b_TypeWaterLevel + b_treatment$b_Intercept)


sd_treatment <- as_draws_df(overall_treatment, variable = c("sd_obs__Intercept:TypeDiet", "sd_obs__Intercept:TypeHumidity", 
                                                            "sd_obs__Intercept:TypePredator", "sd_obs__Intercept:TypeSalinity", 
                                                            "sd_obs__Intercept:TypeTemperature", "sd_obs__Intercept:TypeWaterLevel", 
                                                            "sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_treatment <- data.frame("sd_Diet" = sd_treatment$`sd_obs__Intercept:TypeDiet`,
                           "sd_Humidity" = sd_treatment$`sd_obs__Intercept:TypeHumidity`, 
                           "sd_Predator" = sd_treatment$`sd_obs__Intercept:TypePredator`, 
                           "sd_Salinity" = sd_treatment$`sd_obs__Intercept:TypeSalinity`, 
                           "sd_Temperature" = sd_treatment$`sd_obs__Intercept:TypeTemperature`, 
                           "sd_WaterLevel" = sd_treatment$`sd_obs__Intercept:TypeWaterLevel`,
                           "sd_phylo__Intercept" = sd_treatment$`sd_phylo__Intercept`, 
                           "sd_Study_ID__Intercept" = sd_treatment$`sd_Study_ID__Intercept`, 
                           "sd_Measurement__Intercept" = sd_treatment$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
overall_treatment_means <- apply(b_treatment, 2, mean)
overall_treatment_cis <- apply(b_treatment, 2, function(x) HPDinterval(as.mcmc(x)))
overall_treatment_pMCMC <- apply(b_treatment, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_treatment_diet <- folded_norm(b_treatment$b_Diet, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_Diet"]^2)))
b_abs_treatment_humidity <- folded_norm(b_treatment$b_Humidity, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_Humidity"]^2)))
b_abs_treatment_predator <- folded_norm(b_treatment$b_Predator, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_Predator"]^2)))
b_abs_treatment_salinity <- folded_norm(b_treatment$b_Salinity, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_Salinity"]^2)))
b_abs_treatment_temperature <- folded_norm(b_treatment$b_Temperature, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_Temperature"]^2)))
b_abs_treatment_waterlevel <- folded_norm(b_treatment$b_WaterLevel, sqrt(rowSums(sd_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_treatment[, "sd_WaterLevel"]^2)))
mean_abs_b_treatment_diet <- mean(b_abs_treatment_diet)
mean_abs_b_treatment_humidity <- mean(b_abs_treatment_humidity)
mean_abs_b_treatment_predator <- mean(b_abs_treatment_predator)
mean_abs_b_treatment_salinity <- mean(b_abs_treatment_salinity)
mean_abs_b_treatment_temperature <- mean(b_abs_treatment_temperature)
mean_abs_b_treatment_waterlevel <- mean(b_abs_treatment_waterlevel)
ci.abs_treatment_diet <- HPDinterval(as.mcmc(b_abs_treatment_diet))
ci.abs_treatment_humidity <- HPDinterval(as.mcmc(b_abs_treatment_humidity))
ci.abs_treatment_predator <- HPDinterval(as.mcmc(b_abs_treatment_predator))
ci.abs_treatment_salinity <- HPDinterval(as.mcmc(b_abs_treatment_salinity))
ci.abs_treatment_temperature <- HPDinterval(as.mcmc(b_abs_treatment_temperature))
ci.abs_treatment_waterlevel <- HPDinterval(as.mcmc(b_abs_treatment_waterlevel))

# Heterogeneity
overall_i2_treatment <- i2(sd_treatment, Treatment_Data$Variance_Adjusted) 

# Overall_trait_summary

overall_treatment_means_list <- c(mean_abs_b_treatment_diet, mean_abs_b_treatment_humidity, mean_abs_b_treatment_predator, 
                                  mean_abs_b_treatment_salinity, mean_abs_b_treatment_temperature, mean_abs_b_treatment_waterlevel)
overall_treatment_low_ci <- c(ci.abs_treatment_diet[1], ci.abs_treatment_humidity[1], ci.abs_treatment_predator[1], 
                              ci.abs_treatment_salinity[1], ci.abs_treatment_temperature[1], ci.abs_treatment_waterlevel[1])
overall_treatment_high_ci <- c(ci.abs_treatment_diet[2], ci.abs_treatment_humidity[2], ci.abs_treatment_predator[2], 
                               ci.abs_treatment_salinity[2], ci.abs_treatment_temperature[2], ci.abs_treatment_waterlevel[2])
overall_treatment_categories <- c("Diet", "Humidity", "Predator", 
                                  "Salinity", "Temperature", "Water Level")

overall_treatment_summary <- matrix(c(overall_treatment_means_list, overall_treatment_low_ci, overall_treatment_high_ci), 
                                    nrow = 6, ncol = 3, byrow = FALSE, 
                                    dimnames = list(c(overall_treatment_categories), 
                                                    c("Mean", "Low_CI", "High_CI")))
overall_treatment_summary <- data.frame(overall_treatment_summary)

# Preparing Graph - Combined

Treatment_rnames <- c("Diet", "Humidity", "Predator", 
                      "Salinity", "Temperature", "Water Level")

Treatment_k <- data.frame("k" = c(Treatment_Exploration["Diet", "Freq"], 
                                  Treatment_Exploration["Humidity", "Freq"], 
                                  Treatment_Exploration["Predator", "Freq"], 
                                  Treatment_Exploration["Salinity", "Freq"], 
                                  Treatment_Exploration["Temperature", "Freq"], 
                                  Treatment_Exploration["Water Level", "Freq"]), 
                          row.names = Treatment_rnames)

Treatment_group_no <- data.frame("Spp No." = c(Treatment_Species_Count["Diet", "Freq"], 
                                               Treatment_Species_Count["Humidity", "Freq"], 
                                               Treatment_Species_Count["Predator", "Freq"], 
                                               Treatment_Species_Count["Salinity", "Freq"], 
                                               Treatment_Species_Count["Temperature", "Freq"], 
                                               Treatment_Species_Count["Water Level", "Freq"]), 
                                 row.names = Treatment_rnames)

Treatment_study <- data.frame("Study" = c(Treatment_Study_Count["Diet", "Freq"], 
                                          Treatment_Study_Count["Humidity", "Freq"], 
                                          Treatment_Study_Count["Predator", "Freq"], 
                                          Treatment_Study_Count["Salinity", "Freq"], 
                                          Treatment_Study_Count["Temperature", "Freq"], 
                                          Treatment_Study_Count["Water Level", "Freq"]), 
                              row.names = Treatment_rnames)

Treatment_table <- data.frame(estimate = overall_treatment_summary[,"Mean"], 
                              lowerCL = overall_treatment_summary[,"Low_CI"], 
                              upperCL = overall_treatment_summary[,"High_CI"], 
                              K = Treatment_k[,1], 
                              group_no = Treatment_group_no[,1], 
                              row.names = Treatment_rnames)
Treatment_table$name <- row.names(Treatment_table)

Treatment_raw_mean <- c(b_abs_treatment_diet, b_abs_treatment_humidity, b_abs_treatment_predator, 
                        b_abs_treatment_salinity, b_abs_treatment_temperature, b_abs_treatment_waterlevel)

Treatment_raw_name <- c(replicate(8000, "Diet"), 
                        replicate(8000, "Humidity"), 
                        replicate(8000, "Predator"), 
                        replicate(8000, "Salinity"), 
                        replicate(8000, "Temperature"), 
                        replicate(8000, "Water Level"))

Treatment_raw_df <- data.frame("Model" = Treatment_raw_name, 
                               "Effect" = Treatment_raw_mean)

# Graph code - Combined

Treatment_Order <- c("Water Level", "Temperature", "Salinity", 
                     "Predator", "Humidity", "Diet")

density_treatment <- Treatment_table %>% mutate(name = fct_relevel(name, Treatment_Order)) %>%
                     ggplot() +
                     geom_density_ridges(data = Treatment_raw_df %>% mutate(Model = fct_relevel(Model, Treatment_Order)), 
                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                     geom_linerange(aes(y = rev(seq(1, dim(Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(Treatment_table)[1], 1)), xmin = max(Treatment_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(Treatment_table)[1], 1)), xmin = min(Treatment_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                    size = 1) +
                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                     theme_bw() +
                     guides(fill = "none", colour = "none") +
                     labs(x = TeX("Effect Size (PRRD)"), y = "") +
                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                      vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7, -2.7))) +
                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                     theme(axis.ticks = element_blank()) +
                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                     scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                     scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                     coord_cartesian(xlim = c(-0.01, 1.25)) +
                     annotate('text',  x = 1.25, y = (seq(1, dim(Treatment_table)[1], 1)+0.4),
                     label= paste("italic(k)==", c(Treatment_table["Water Level", "K"], 
                                                   Treatment_table["Temperature", "K"], 
                                                   Treatment_table["Salinity", "K"], 
                                                   Treatment_table["Predator", "K"], 
                                                   Treatment_table["Humidity", "K"], 
                                                   Treatment_table["Diet", "K"]), "~","(", 
                                                 c(Treatment_table["Water Level", "group_no"], 
                                                   Treatment_table["Temperature", "group_no"], 
                                                   Treatment_table["Salinity", "group_no"], 
                                                   Treatment_table["Predator", "group_no"], 
                                                   Treatment_table["Humidity", "group_no"], 
                                                   Treatment_table["Diet", "group_no"]), 
                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                     geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_treatment_waterlevel)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_treatment_temperature)-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(b_abs_treatment_salinity)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_treatment_predator)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_treatment_humidity)-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(b_abs_treatment_diet)-1)*100, 2), nsmall = 2), "%")), 
                                    x = rev(Treatment_table$estimate+0.2), y = (seq(1, dim(Treatment_table)[1], 1)+0.4)), size = 3.5)

density_treatment #(400x560)

# Preparing Graph - Part 1

Treatment_rnames_1 <- c("Diet", "Humidity", "Predator")

Treatment_k_1 <- data.frame("k" = c(Treatment_Exploration["Diet", "Freq"], 
                                    Treatment_Exploration["Humidity", "Freq"], 
                                    Treatment_Exploration["Predator", "Freq"]), 
                            row.names = Treatment_rnames_1)

Treatment_group_no_1 <- data.frame("Spp No." = c(Treatment_Species_Count["Diet", "Freq"], 
                                                 Treatment_Species_Count["Humidity", "Freq"], 
                                                 Treatment_Species_Count["Predator", "Freq"]), 
                                   row.names = Treatment_rnames_1)

Treatment_study_1 <- data.frame("Study" = c(Treatment_Study_Count["Diet", "Freq"], 
                                            Treatment_Study_Count["Humidity", "Freq"], 
                                            Treatment_Study_Count["Predator", "Freq"]), 
                                row.names = Treatment_rnames_1)

overall_treatment_summary_Reorder_1 <- overall_treatment_summary[c("Diet", "Humidity", "Predator"), ]

Treatment_table_1 <- data.frame(estimate = overall_treatment_summary_Reorder_1[,"Mean"], 
                                lowerCL = overall_treatment_summary_Reorder_1[,"Low_CI"], 
                                upperCL = overall_treatment_summary_Reorder_1[,"High_CI"], 
                                K = Treatment_k_1[,1], 
                                group_no = Treatment_group_no_1[,1], 
                                row.names = Treatment_rnames_1)
Treatment_table_1$name <- row.names(Treatment_table_1)

Treatment_raw_mean_1 <- c(b_abs_treatment_diet, b_abs_treatment_humidity, b_abs_treatment_predator)

Treatment_raw_name_1 <- c(replicate(8000, "Diet"), 
                          replicate(8000, "Humidity"), 
                          replicate(8000, "Predator"))

Treatment_raw_df_1 <- data.frame("Model" = Treatment_raw_name_1, 
                                 "Effect" = Treatment_raw_mean_1)

# Graph code - Part 1

Treatment_Order_1 <- c("Predator", "Humidity", "Diet")

density_treatment_1 <- Treatment_table_1 %>% mutate(name = fct_relevel(name, Treatment_Order_1)) %>%
                       ggplot() +
                       geom_density_ridges(data = Treatment_raw_df_1 %>% mutate(Model = fct_relevel(Model, Treatment_Order_1)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_1)[1], 1)), xmin = max(Treatment_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_1)[1], 1)), xmin = min(Treatment_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                       scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                       annotate('text',  x = 1.25, y = (seq(1, dim(Treatment_table_1)[1], 1)+0.4),
                       label= paste("italic(k)==", c(Treatment_table_1["Predator", "K"], 
                                                     Treatment_table_1["Humidity", "K"], 
                                                     Treatment_table_1["Diet", "K"]), "~","(", 
                                                   c(Treatment_table_1["Predator", "group_no"], 
                                                     Treatment_table_1["Humidity", "group_no"], 
                                                     Treatment_table_1["Diet", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_treatment_predator)-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(b_abs_treatment_humidity)-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(b_abs_treatment_diet)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Treatment_table_1$estimate+0.2), y = (seq(1, dim(Treatment_table_1)[1], 1)+0.4)), size = 3.5)

density_treatment_1 #(400x320)

# Preparing Graph - Part 2

Treatment_rnames_2 <- c("Salinity", "Temperature", "Water Level")

Treatment_k_2 <- data.frame("k" = c(Treatment_Exploration["Salinity", "Freq"], 
                                    Treatment_Exploration["Temperature", "Freq"], 
                                    Treatment_Exploration["Water Level", "Freq"]), 
                            row.names = Treatment_rnames_2)

Treatment_group_no_2 <- data.frame("Spp No." = c(Treatment_Species_Count["Salinity", "Freq"], 
                                                 Treatment_Species_Count["Temperature", "Freq"], 
                                                 Treatment_Species_Count["Water Level", "Freq"]), 
                                   row.names = Treatment_rnames_2)

Treatment_study_2 <- data.frame("Study" = c(Treatment_Study_Count["Salinity", "Freq"], 
                                            Treatment_Study_Count["Temperature", "Freq"], 
                                            Treatment_Study_Count["Water Level", "Freq"]), 
                                row.names = Treatment_rnames_2)

overall_treatment_summary_Reorder_2 <- overall_treatment_summary[c("Salinity", "Temperature", "Water Level"), ]

Treatment_table_2 <- data.frame(estimate = overall_treatment_summary_Reorder_2[,"Mean"], 
                                lowerCL = overall_treatment_summary_Reorder_2[,"Low_CI"], 
                                upperCL = overall_treatment_summary_Reorder_2[,"High_CI"], 
                                K = Treatment_k_2[,1], 
                                group_no = Treatment_group_no_2[,1], 
                                row.names = Treatment_rnames_2)
Treatment_table_2$name <- row.names(Treatment_table_2)

Treatment_raw_mean_2 <- c(b_abs_treatment_salinity, b_abs_treatment_temperature, b_abs_treatment_waterlevel)

Treatment_raw_name_2 <- c(replicate(8000, "Salinity"), 
                          replicate(8000, "Temperature"), 
                          replicate(8000, "Water Level"))

Treatment_raw_df_2 <- data.frame("Model" = Treatment_raw_name_2, 
                                 "Effect" = Treatment_raw_mean_2)

# Graph code - Part 2

Treatment_Order_2 <- c("Water Level", "Temperature", "Salinity")

density_treatment_2 <- Treatment_table_2 %>% mutate(name = fct_relevel(name, Treatment_Order_2)) %>%
                       ggplot() +
                       geom_density_ridges(data = Treatment_raw_df_2 %>% mutate(Model = fct_relevel(Model, Treatment_Order_2)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_2)[1], 1)), xmin = max(Treatment_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Treatment_table_2)[1], 1)), xmin = min(Treatment_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                       annotate('text',  x = 1.25, y = (seq(1, dim(Treatment_table_2)[1], 1)+0.4),
                       label= paste("italic(k)==", c(Treatment_table_2["Water Level", "K"], 
                                                     Treatment_table_2["Temperature", "K"], 
                                                     Treatment_table_2["Salinity", "K"]), "~","(", 
                                                   c(Treatment_table_2["Water Level", "group_no"], 
                                                     Treatment_table_2["Temperature", "group_no"], 
                                                     Treatment_table_2["Salinity", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                        geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_treatment_waterlevel)-1)*100, 2), nsmall = 2), "%"), 
                                               paste(format(round(mean(exp(b_abs_treatment_temperature)-1)*100, 2), nsmall = 2), "%"),
                                               paste(format(round(mean(exp(b_abs_treatment_salinity)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Treatment_table_2$estimate+0.2), y = (seq(1, dim(Treatment_table_2)[1], 1)+0.4)), size = 3.5)

density_treatment_2 #(400x320)

##### Terrestrial Model #####
Terrestrial_Data <- data %>% filter(Ecosystem == "Terrestrial")
Terrestrial_Species <- Terrestrial_Data %>% select("phylo") %>% unique()

Terrestrial_A <- as.data.frame(A)
Terrestrial_A <- Terrestrial_A[c(Terrestrial_Species$phylo), c(Terrestrial_Species$phylo)]
Terrestrial_A <- as.matrix(Terrestrial_A)

system.time(
  terrestrial <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                           ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|obs),
                           data = Terrestrial_Data,
                           family = gaussian(),
                           data2 = list(A = Terrestrial_A), 
                           chains = 4, 
                           cores = 4,
                           iter = 12000,
                           warmup = 2000,
                           thin = 5,
                           prior = priors,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           file = "./terrestrial_model",
                           file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_terrestrial <- as_draws_df(terrestrial, variable = "b_Intercept")
b_terrestrial <- data.frame(b_terrestrial$b_Intercept)

sd_terrestrial <- as_draws_df(terrestrial, variable = c("sd_Measurement__Intercept", "sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_terrestrial <- data.frame("sd_Measurement__Intercept" = sd_terrestrial$sd_Measurement__Intercept, 
                             "sd_obs__Intercept" = sd_terrestrial$sd_obs__Intercept, 
                             "sd_phylo__Intercept" = sd_terrestrial$sd_phylo__Intercept, 
                             "sd_Study_ID__Intercept" = sd_terrestrial$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_terrestrial <-  sapply(b_terrestrial, mean)
ci_b_terrestrial <- as.vector(HPDinterval(as.mcmc(b_terrestrial)))
pMCMC_b_terrestrial <- 2*(1 - max(table(b_terrestrial<0) / nrow(b_terrestrial)))

# Absolute magnitude
b_abs_terrestrial <- folded_norm(b_terrestrial[,1], rowSums(sd_terrestrial))
mean_abs_b_terrestrial <- mean(b_abs_terrestrial)
ci.abs_terrestrial <- HPDinterval(as.mcmc(b_abs_terrestrial))

# Heterogeneity
terrestrial_i2 <- i2(sd_terrestrial, Terrestrial_Data$Variance_Adjusted) 

##### Terrestrial Model - Plasticity Mechanism Meta-regression #####
Terrestrial_Plasticity_Exploration <- Terrestrial_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Exploration) <- Terrestrial_Plasticity_Exploration$Plasticity_Category

Terrestrial_Plasticity_Data <- Terrestrial_Data %>% filter(Plasticity_Category == "Acclimation"| 
                                                           Plasticity_Category == "Developmental")

Terrestrial_Plasticity_Species_Count <- Terrestrial_Plasticity_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                        table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                        select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Species_Count) <- Terrestrial_Plasticity_Species_Count$Plasticity_Category

Terrestrial_Plasticity_Study_Count <- Terrestrial_Plasticity_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                      table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                      select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Study_Count) <- Terrestrial_Plasticity_Study_Count$Plasticity_Category

Terrestrial_Plasticity_Species <- Terrestrial_Plasticity_Data %>% select("phylo") %>% unique()

Terrestrial_Plasticity_A <- as.data.frame(A)
Terrestrial_Plasticity_A <- Terrestrial_Plasticity_A[c(Terrestrial_Plasticity_Species$phylo), c(Terrestrial_Plasticity_Species$phylo)]
Terrestrial_Plasticity_A <- as.matrix(Terrestrial_Plasticity_A)


system.time(
  terrestrial_plastic <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                                   ~ Plasticity_Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Plasticity_Category, cor = FALSE)),
                                   data = Terrestrial_Plasticity_Data,
                                   family = gaussian(),
                                   data2 = list(A = Terrestrial_Plasticity_A), 
                                   chains = 4, 
                                   cores = 4,
                                   iter = 12000,
                                   warmup = 2000,
                                   thin = 5,
                                   prior = priors,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                                   file = "./terrestrial_plastic_model",
                                   file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_terrestrial_plastic <- as_draws_df(terrestrial_plastic, variable = c("b_Intercept", 
                                                                       "b_Plasticity_CategoryDevelopmental"))
b_terrestrial_plastic <- data.frame("b_Acclimation" = b_terrestrial_plastic$b_Intercept, 
                                    "b_Developmental" = b_terrestrial_plastic$b_Plasticity_CategoryDevelopmental + b_terrestrial_plastic$b_Intercept)


sd_terrestrial_plastic <- as_draws_df(terrestrial_plastic, variable = c("sd_obs__Intercept:Plasticity_CategoryAcclimation",
                                                                        "sd_obs__Intercept:Plasticity_CategoryDevelopmental",
                                                                        "sd_phylo__Intercept", "sd_Study_ID__Intercept", 
                                                                        "sd_Measurement__Intercept"))
sd_terrestrial_plastic <- data.frame("sd_Acclimation" = sd_terrestrial_plastic$`sd_obs__Intercept:Plasticity_CategoryAcclimation`,
                                     "sd_Developmental" = sd_terrestrial_plastic$`sd_obs__Intercept:Plasticity_CategoryDevelopmental`,
                                     "sd_phylo__Intercept" = sd_terrestrial_plastic$`sd_phylo__Intercept`, 
                                     "sd_Study_ID__Intercept" = sd_terrestrial_plastic$`sd_Study_ID__Intercept`, 
                                     "sd_Measurement__Intercept" = sd_terrestrial_plastic$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
terrestrial_plastic_means <- apply(b_terrestrial_plastic, 2, mean)
terrestrial_plastic_cis <- apply(b_terrestrial_plastic, 2, function(x) HPDinterval(as.mcmc(x)))
terrestrial_plastic_pMCMC <- apply(b_terrestrial_plastic, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_terrestrial_plastic_acc <- folded_norm(b_terrestrial_plastic$b_Acclimation, sqrt(rowSums(sd_terrestrial_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_plastic[, "sd_Acclimation"]^2)))
b_abs_terrestrial_plastic_dev <- folded_norm(b_terrestrial_plastic$b_Developmental, sqrt(rowSums(sd_terrestrial_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_plastic[, "sd_Developmental"]^2)))
mean_abs_b_terrestrial_plastic_acc <- mean(b_abs_terrestrial_plastic_acc)
mean_abs_b_terrestrial_plastic_dev <- mean(b_abs_terrestrial_plastic_dev)
ci.abs_terrestrial_acc <- HPDinterval(as.mcmc(b_abs_terrestrial_plastic_acc))
ci.abs_terrestrial_dev <- HPDinterval(as.mcmc(b_abs_terrestrial_plastic_dev))

# Heterogeneity
terrestrial_i2_plastic <- i2(sd_terrestrial_plastic, Terrestrial_Plasticity_Data$Variance_Adjusted) 

# Overall_Plasticity_Summary
terrestrial_plastic_means_list <- c(mean_abs_b_terrestrial_plastic_acc, mean_abs_b_terrestrial_plastic_dev)
terrestrial_plastic_low_ci <- c(ci.abs_terrestrial_acc[1], ci.abs_terrestrial_dev[1])
terrestrial_plastic_high_ci <- c(ci.abs_terrestrial_acc[2], ci.abs_terrestrial_dev[2])
terrestrial_plastic_categories <- c("Acclimation", "Developmental Plasticity")

terrestrial_plastic_summary <- matrix(c(terrestrial_plastic_means_list, terrestrial_plastic_low_ci, terrestrial_plastic_high_ci), 
                                  nrow = 2, ncol = 3, byrow = FALSE, 
                                  dimnames = list(c(terrestrial_plastic_categories), 
                                                  c("Mean", "Low_CI", "High_CI")))
terrestrial_plastic_summary <- data.frame(terrestrial_plastic_summary)

# Preparing Graph - Combined

Terrestrial_Plasticity_rnames <- c("Acclimation", "Developmental Plasticity")

Terrestrial_Plasticity_k <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"], 
                                               Terrestrial_Plasticity_Exploration["Developmental", "Freq"]), 
                                       row.names = Terrestrial_Plasticity_rnames)

Terrestrial_Plasticity_group_no <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                            Terrestrial_Plasticity_Species_Count["Developmental", "Freq"]), 
                                              row.names = Terrestrial_Plasticity_rnames)

Terrestrial_Plasticity_study <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                       Terrestrial_Plasticity_Study_Count["Developmental", "Freq"]), 
                                           row.names = Terrestrial_Plasticity_rnames)

Terrestrial_Plasticity_table <- data.frame(estimate = terrestrial_plastic_summary[,"Mean"], 
                                           lowerCL = terrestrial_plastic_summary[,"Low_CI"], 
                                           upperCL = terrestrial_plastic_summary[,"High_CI"], 
                                           K = Terrestrial_Plasticity_k[,1], 
                                           group_no = Terrestrial_Plasticity_group_no[,1], 
                                           row.names = Terrestrial_Plasticity_rnames)
Terrestrial_Plasticity_table$name <- row.names(Terrestrial_Plasticity_table)

Terrestrial_Plasticity_raw_mean <- c(b_abs_terrestrial_plastic_acc, b_abs_terrestrial_plastic_dev)

Terrestrial_Plasticity_raw_name <- c(replicate(8000, "Acclimation"), 
                                     replicate(8000, "Developmental Plasticity"))

Terrestrial_Plasticity_raw_df <- data.frame("Model" = Terrestrial_Plasticity_raw_name, 
                                            "Effect" = Terrestrial_Plasticity_raw_mean)

# Graph code - Combined

Terrestrial_Plasticity_Order <- c("Developmental Plasticity", "Acclimation")

density_terrestrial_plasticity <- Terrestrial_Plasticity_table %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = Terrestrial_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table)[1], 1)), xmin = max(Terrestrial_Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table)[1], 1)), xmin = min(Terrestrial_Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-0.8, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.01, 1.25)) +
                                  annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Plasticity_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(Terrestrial_Plasticity_table["Developmental Plasticity", "K"], 
                                                                Terrestrial_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                              c(Terrestrial_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                                Terrestrial_Plasticity_table["Acclimation", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(b_abs_terrestrial_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Terrestrial_Plasticity_table$estimate+0.2), y = (seq(1, dim(Terrestrial_Plasticity_table)[1], 1)+0.4)), size = 3.5)
density_terrestrial_plasticity #(400x240)

# Preparing Graph - Part 1

Terrestrial_Plasticity_rnames_1 <- c("Acclimation")

Terrestrial_Plasticity_k_1 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"]), 
                                         row.names = Terrestrial_Plasticity_rnames_1)

Terrestrial_Plasticity_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                row.names = Terrestrial_Plasticity_rnames_1)

Terrestrial_Plasticity_study_1 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                             row.names = Terrestrial_Plasticity_rnames_1)

terrestrial_plastic_summary_Reorder_1 <- terrestrial_plastic_summary[c("Acclimation"), ]

Terrestrial_Plasticity_table_1 <- data.frame(estimate = terrestrial_plastic_summary_Reorder_1[,"Mean"], 
                                             lowerCL = terrestrial_plastic_summary_Reorder_1[,"Low_CI"], 
                                             upperCL = terrestrial_plastic_summary_Reorder_1[,"High_CI"], 
                                             K = Terrestrial_Plasticity_k_1[,1], 
                                             group_no = Terrestrial_Plasticity_group_no_1[,1], 
                                             row.names = Terrestrial_Plasticity_rnames_1)
Terrestrial_Plasticity_table_1$name <- row.names(Terrestrial_Plasticity_table_1)

Terrestrial_Plasticity_raw_mean_1 <- c(b_abs_terrestrial_plastic_acc)

Terrestrial_Plasticity_raw_name_1 <- c(replicate(8000, "Acclimation"))

Terrestrial_Plasticity_raw_df_1 <- data.frame("Model" = Terrestrial_Plasticity_raw_name_1, 
                                              "Effect" = Terrestrial_Plasticity_raw_mean_1)

# Graph code - Part 1

Terrestrial_Plasticity_Order_1 <- c("Acclimation")

density_terrestrial_plasticity_1 <- Terrestrial_Plasticity_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = Terrestrial_Plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.08, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)), xmin = max(Terrestrial_Plasticity_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)), xmin = min(Terrestrial_Plasticity_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                     vjust = c(-2.7))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#2B4E7A")) +
                                    scale_fill_manual(values = c("#2B4E7A")) +
                                    coord_cartesian(xlim = c(-0.01, 1.25)) +
                                    annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(Terrestrial_Plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                c(Terrestrial_Plasticity_table_1["Acclimation", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                   x = rev(Terrestrial_Plasticity_table_1$estimate+0.2), y = (seq(1, dim(Terrestrial_Plasticity_table_1)[1], 1)+0.4)), size = 3.5)
density_terrestrial_plasticity_1 #(400x160)

# Preparing Graph - Part 2

Terrestrial_Plasticity_rnames_2 <- c("Developmental Plasticity")

Terrestrial_Plasticity_k_2 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Developmental", "Freq"]), 
                                         row.names = Terrestrial_Plasticity_rnames_2)

Terrestrial_Plasticity_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Developmental", "Freq"]), 
                                                row.names = Terrestrial_Plasticity_rnames_2)

Terrestrial_Plasticity_study_2 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Developmental", "Freq"]), 
                                             row.names = Terrestrial_Plasticity_rnames_2)

terrestrial_plastic_summary_Reorder_2 <- terrestrial_plastic_summary[c("Developmental Plasticity"), ]

Terrestrial_Plasticity_table_2 <- data.frame(estimate = terrestrial_plastic_summary_Reorder_2[,"Mean"], 
                                             lowerCL = terrestrial_plastic_summary_Reorder_2[,"Low_CI"], 
                                             upperCL = terrestrial_plastic_summary_Reorder_2[,"High_CI"], 
                                             K = Terrestrial_Plasticity_k_2[,1], 
                                             group_no = Terrestrial_Plasticity_group_no_2[,1], 
                                             row.names = Terrestrial_Plasticity_rnames_2)
Terrestrial_Plasticity_table_2$name <- row.names(Terrestrial_Plasticity_table_2)

Terrestrial_Plasticity_raw_mean_2 <- c(b_abs_terrestrial_plastic_dev)

Terrestrial_Plasticity_raw_name_2 <- c(replicate(8000, "Developmental Plasticity"))

Terrestrial_Plasticity_raw_df_2 <- data.frame("Model" = Terrestrial_Plasticity_raw_name_2, 
                                              "Effect" = Terrestrial_Plasticity_raw_mean_2)

# Graph code - Part 2

Terrestrial_Plasticity_Order_2 <- c("Developmental Plasticity")

density_terrestrial_plasticity_2 <- Terrestrial_Plasticity_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = Terrestrial_Plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.08, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)), xmin = max(Terrestrial_Plasticity_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)), xmin = min(Terrestrial_Plasticity_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                     vjust = c(-0.8))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#5D7AA1")) +
                                    scale_fill_manual(values = c("#5D7AA1")) +
                                    coord_cartesian(xlim = c(-0.01, 1.25)) +
                                    annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(Terrestrial_Plasticity_table_2["Developmental Plasticity", "K"]), "~","(", 
                                                                c(Terrestrial_Plasticity_table_2["Developmental Plasticity", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_plastic_dev)-1)*100, 2), nsmall = 2), "%")), 
                                                   x = rev(Terrestrial_Plasticity_table_2$estimate+0.2), y = (seq(1, dim(Terrestrial_Plasticity_table_2)[1], 1)+0.4)), size = 3.5)
density_terrestrial_plasticity_2 #(400x160)

##### Terrestrial Model - Trait Category Meta-regression #####
Terrestrial_Trait_Exploration <- Terrestrial_Data %>% select("Category") %>% table() %>% data.frame()
Terrestrial_Trait_Exploration <- Terrestrial_Trait_Exploration %>% filter(Freq > 10)
rownames(Terrestrial_Trait_Exploration) <- Terrestrial_Trait_Exploration$Category

Terrestrial_Trait_Data <- Terrestrial_Data %>% filter(Category == "Behavioural"| 
                                                      Category == "Biochemical Assay"| 
                                                      Category == "Gene Expression"|
                                                      Category == "Morphology"| 
                                                      Category == "Tolerance")

Terrestrial_Trait_Species_Count <- Terrestrial_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                                   table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                   select("Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Species_Count) <- Terrestrial_Trait_Species_Count$Category

Terrestrial_Trait_Study_Count <- Terrestrial_Trait_Data %>% select("Study_ID", "Category") %>% 
                                 table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                 select("Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Study_Count) <- Terrestrial_Trait_Study_Count$Category

Terrestrial_Trait_Species <- Terrestrial_Trait_Data %>% select("phylo") %>% unique()

Terrestrial_Trait_A <- as.data.frame(A)
Terrestrial_Trait_A <- Terrestrial_Trait_A[c(Terrestrial_Trait_Species$phylo), c(Terrestrial_Trait_Species$phylo)]
Terrestrial_Trait_A <- as.matrix(Terrestrial_Trait_A)

system.time(
  terrestrial_trait <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                                 ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Category, cor = FALSE)),
                                 data = Terrestrial_Trait_Data,
                                 family = gaussian(),
                                 data2 = list(A = Terrestrial_Trait_A), 
                                 chains = 4, 
                                 cores = 4,
                                 iter = 12000,
                                 warmup = 2000,
                                 thin = 5,
                                 prior = priors,
                                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                                 file = "./terrestrial_trait_model",
                                 file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_terrestrial_trait <- as_draws_df(terrestrial_trait, variable = c("b_Intercept", "b_CategoryBiochemicalAssay", 
                                                                   "b_CategoryGeneExpression", "b_CategoryMorphology",
                                                                   "b_CategoryTolerance"))
b_terrestrial_trait <- data.frame("b_Behavioural" = b_terrestrial_trait$b_Intercept, 
                                  "b_BiochemicalAssay" = b_terrestrial_trait$b_CategoryBiochemicalAssay + b_terrestrial_trait$b_Intercept, 
                                  "b_GeneExpression" = b_terrestrial_trait$b_CategoryGeneExpression + b_terrestrial_trait$b_Intercept, 
                                  "b_Morphology" = b_terrestrial_trait$b_CategoryMorphology + b_terrestrial_trait$b_Intercept, 
                                  "b_Tolerance" = b_terrestrial_trait$b_CategoryTolerance + b_terrestrial_trait$b_Intercept)


sd_terrestrial_trait <- as_draws_df(terrestrial_trait, variable = c("sd_obs__Intercept:CategoryBehavioural", "sd_obs__Intercept:CategoryBiochemicalAssay", 
                                                                    "sd_obs__Intercept:CategoryGeneExpression", "sd_obs__Intercept:CategoryMorphology",
                                                                    "sd_obs__Intercept:CategoryTolerance", "sd_phylo__Intercept", 
                                                                    "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_terrestrial_trait <- data.frame("sd_Behavioural" = sd_terrestrial_trait$`sd_obs__Intercept:CategoryBehavioural`,
                                   "sd_BiochemicalAssay" = sd_terrestrial_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`, 
                                   "sd_GeneExpression" = sd_terrestrial_trait$`sd_obs__Intercept:CategoryGeneExpression`, 
                                   "sd_Morphology" = sd_terrestrial_trait$`sd_obs__Intercept:CategoryMorphology`, 
                                   "sd_Tolerance" = sd_terrestrial_trait$`sd_obs__Intercept:CategoryTolerance`,
                                   "sd_phylo__Intercept" = sd_terrestrial_trait$`sd_phylo__Intercept`, 
                                   "sd_Study_ID__Intercept" = sd_terrestrial_trait$`sd_Study_ID__Intercept`, 
                                   "sd_Measurement__Intercept" = sd_terrestrial_trait$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
terrestrial_trait_means <- apply(b_terrestrial_trait, 2, mean)
terrestrial_trait_cis <- apply(b_terrestrial_trait, 2, function(x) HPDinterval(as.mcmc(x)))
terrestrial_trait_pMCMC <- apply(b_terrestrial_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_terrestrial_trait_behavioural <- folded_norm(b_terrestrial_trait$b_Behavioural, sqrt(rowSums(sd_terrestrial_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_trait[, "sd_Behavioural"]^2)))
b_abs_terrestrial_trait_biochem <- folded_norm(b_terrestrial_trait$b_BiochemicalAssay, sqrt(rowSums(sd_terrestrial_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_terrestrial_trait_gene <- folded_norm(b_terrestrial_trait$b_GeneExpression, sqrt(rowSums(sd_terrestrial_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_trait[, "sd_GeneExpression"]^2)))
b_abs_terrestrial_trait_morphology <- folded_norm(b_terrestrial_trait$b_Morphology, sqrt(rowSums(sd_terrestrial_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_trait[, "sd_Morphology"]^2)))
b_abs_terrestrial_trait_tolerance <- folded_norm(b_terrestrial_trait$b_Tolerance, sqrt(rowSums(sd_terrestrial_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_trait[, "sd_Tolerance"]^2)))
mean_abs_b_terrestrial_trait_behavioural <- mean(b_abs_terrestrial_trait_behavioural)
mean_abs_b_terrestrial_trait_biochem <- mean(b_abs_terrestrial_trait_biochem)
mean_abs_b_terrestrial_trait_gene <- mean(b_abs_terrestrial_trait_gene)
mean_abs_b_terrestrial_trait_morphology <- mean(b_abs_terrestrial_trait_morphology)
mean_abs_b_terrestrial_trait_tolerance <- mean(b_abs_terrestrial_trait_tolerance)
ci.abs_terrestrial_behavioural <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_behavioural))
ci.abs_terrestrial_biochem <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_biochem))
ci.abs_terrestrial_gene <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_gene))
ci.abs_terrestrial_morphology <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_morphology))
ci.abs_terrestrial_tolerance <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_tolerance))

# Heterogeneity
terrestrial_i2_trait <- i2(sd_terrestrial_trait, Terrestrial_Trait_Data$Variance_Adjusted) 

# Overall_trait_summary

terrestrial_trait_means_list <- c(mean_abs_b_terrestrial_trait_behavioural, mean_abs_b_terrestrial_trait_biochem, mean_abs_b_terrestrial_trait_gene, 
                                  mean_abs_b_terrestrial_trait_morphology, mean_abs_b_terrestrial_trait_tolerance)
terrestrial_trait_low_ci <- c(ci.abs_terrestrial_behavioural[1], ci.abs_terrestrial_biochem[1], ci.abs_terrestrial_gene[1], 
                              ci.abs_terrestrial_morphology[1], ci.abs_terrestrial_tolerance[1])
terrestrial_trait_high_ci <- c(ci.abs_terrestrial_behavioural[2], ci.abs_terrestrial_biochem[2], ci.abs_terrestrial_gene[2], 
                               ci.abs_terrestrial_morphology[2], ci.abs_terrestrial_tolerance[2])
terrestrial_trait_categories <- c("Behavioural", "Biochemical Assay", "Gene Expression",
                                  "Morphology", "Tolerance")

terrestrial_trait_summary <- matrix(c(terrestrial_trait_means_list, terrestrial_trait_low_ci, terrestrial_trait_high_ci), 
                                    nrow = 5, ncol = 3, byrow = FALSE, 
                                    dimnames = list(c(terrestrial_trait_categories), 
                                                    c("Mean", "Low_CI", "High_CI")))
terrestrial_trait_summary <- data.frame(terrestrial_trait_summary)

# Preparing Graph - Combined

Terrestrial_Trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression",
                              "Morphological", "Tolerance")

Terrestrial_Trait_k <- data.frame("k" = c(Terrestrial_Trait_Exploration["Behavioural", "Freq"], 
                                          Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                          Terrestrial_Trait_Exploration["Gene Expression", "Freq"],
                                          Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                          Terrestrial_Trait_Exploration["Tolerance", "Freq"]), 
                                  row.names = Terrestrial_Trait_rnames)

Terrestrial_Trait_group_no <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Behavioural", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Gene Expression", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Morphology", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Tolerance", "Freq"]), 
                                         row.names = Terrestrial_Trait_rnames)

Terrestrial_Trait_study <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Behavioural", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Gene Expression", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Morphology", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Tolerance", "Freq"]), 
                                      row.names = Terrestrial_Trait_rnames)

terrestrial_trait_summary_Reorder <- terrestrial_trait_summary[c("Behavioural", "Biochemical Assay", "Gene Expression",
                                                                 "Morphology", "Tolerance"), ]

Terrestrial_Trait_table <- data.frame(estimate = terrestrial_trait_summary_Reorder[,"Mean"], 
                                      lowerCL = terrestrial_trait_summary_Reorder[,"Low_CI"], 
                                      upperCL = terrestrial_trait_summary_Reorder[,"High_CI"], 
                                      K = Terrestrial_Trait_k[,1], 
                                      group_no = Terrestrial_Trait_group_no[,1], 
                                      row.names = Terrestrial_Trait_rnames)
Terrestrial_Trait_table$name <- row.names(Terrestrial_Trait_table)

Terrestrial_Trait_raw_mean <- c(b_abs_terrestrial_trait_behavioural, b_abs_terrestrial_trait_biochem, b_abs_terrestrial_trait_gene, 
                                b_abs_terrestrial_trait_morphology, b_abs_terrestrial_trait_tolerance)

Terrestrial_Trait_raw_name <- c(replicate(8000, "Behavioural"), 
                                replicate(8000, "Biochemical Assay"), 
                                replicate(8000, "Gene Expression"),
                                replicate(8000, "Morphological"), 
                                replicate(8000, "Tolerance"))

Terrestrial_Trait_raw_df <- data.frame("Model" = Terrestrial_Trait_raw_name, 
                                       "Effect" = Terrestrial_Trait_raw_mean)

# Graph code - Combined

Terrestrial_Trait_Order <- c("Tolerance", "Morphological",
                             "Gene Expression", "Biochemical Assay", "Behavioural")

density_terrestrial_trait <- Terrestrial_Trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = Terrestrial_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table)[1], 1)), xmin = max(Terrestrial_Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table)[1], 1)), xmin = min(Terrestrial_Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (PRRD)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                              vjust = c(-2.7, -2.7, -0.8, -0.8, -2.7))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             coord_cartesian(xlim = c(-0.01, 1.25)) +
                             annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Trait_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(Terrestrial_Trait_table["Tolerance", "K"], 
                                                           Terrestrial_Trait_table["Morphological", "K"],
                                                           Terrestrial_Trait_table["Gene Expression", "K"], 
                                                           Terrestrial_Trait_table["Biochemical Assay", "K"], 
                                                           Terrestrial_Trait_table["Behavioural", "K"]), "~","(", 
                                                         c(Terrestrial_Trait_table["Tolerance", "group_no"], 
                                                           Terrestrial_Trait_table["Morphological", "group_no"],  
                                                           Terrestrial_Trait_table["Gene Expression", "group_no"], 
                                                           Terrestrial_Trait_table["Biochemical Assay", "group_no"], 
                                                           Terrestrial_Trait_table["Behavioural", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_terrestrial_trait_morphology)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_terrestrial_trait_gene)-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(b_abs_terrestrial_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_terrestrial_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                            x = rev(Terrestrial_Trait_table$estimate+0.2), y = (seq(1, dim(Terrestrial_Trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait #(400x400)

# Preparing Graph - Part 1

Terrestrial_Trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression")

Terrestrial_Trait_k_1 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Behavioural", "Freq"], 
                                            Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                            Terrestrial_Trait_Exploration["Gene Expression", "Freq"]), 
                                    row.names = Terrestrial_Trait_rnames_1)

Terrestrial_Trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Behavioural", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Gene Expression", "Freq"]), 
                                           row.names = Terrestrial_Trait_rnames_1)

Terrestrial_Trait_study_1 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Behavioural", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Gene Expression", "Freq"]), 
                                        row.names = Terrestrial_Trait_rnames_1)

terrestrial_trait_summary_Reorder_1 <- terrestrial_trait_summary[c("Behavioural", "Biochemical Assay", "Gene Expression"), ]

Terrestrial_Trait_table_1 <- data.frame(estimate = terrestrial_trait_summary_Reorder_1[,"Mean"], 
                                        lowerCL = terrestrial_trait_summary_Reorder_1[,"Low_CI"], 
                                        upperCL = terrestrial_trait_summary_Reorder_1[,"High_CI"], 
                                        K = Terrestrial_Trait_k_1[,1], 
                                        group_no = Terrestrial_Trait_group_no_1[,1], 
                                        row.names = Terrestrial_Trait_rnames_1)
Terrestrial_Trait_table_1$name <- row.names(Terrestrial_Trait_table_1)

Terrestrial_Trait_raw_mean_1 <- c(b_abs_terrestrial_trait_behavioural, b_abs_terrestrial_trait_biochem, b_abs_terrestrial_trait_gene)

Terrestrial_Trait_raw_name_1 <- c(replicate(8000, "Behavioural"), 
                                  replicate(8000, "Biochemical Assay"), 
                                  replicate(8000, "Gene Expression"))

Terrestrial_Trait_raw_df_1 <- data.frame("Model" = Terrestrial_Trait_raw_name_1, 
                                         "Effect" = Terrestrial_Trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Trait_Order_1 <- c("Gene Expression", "Biochemical Assay", "Behavioural")

density_terrestrial_trait_1 <- Terrestrial_Trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = Terrestrial_Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_1)[1], 1)), xmin = max(Terrestrial_Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_1)[1], 1)), xmin = min(Terrestrial_Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-0.8, -0.8, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Trait_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Terrestrial_Trait_table_1["Gene Expression", "K"], 
                                                             Terrestrial_Trait_table_1["Biochemical Assay", "K"], 
                                                             Terrestrial_Trait_table_1["Behavioural", "K"]), "~","(", 
                                                           c(Terrestrial_Trait_table_1["Gene Expression", "group_no"], 
                                                             Terrestrial_Trait_table_1["Biochemical Assay", "group_no"], 
                                                             Terrestrial_Trait_table_1["Behavioural", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_trait_gene)-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(b_abs_terrestrial_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_terrestrial_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Terrestrial_Trait_table_1$estimate+0.2), y = (seq(1, dim(Terrestrial_Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_1 #(400x320)

# Preparing Graph - Part 2

Terrestrial_Trait_rnames_2 <- c("Morphological", "Tolerance")

Terrestrial_Trait_k_2 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                            Terrestrial_Trait_Exploration["Tolerance", "Freq"]), 
                                    row.names = Terrestrial_Trait_rnames_2)

Terrestrial_Trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Morphology", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Tolerance", "Freq"]), 
                                           row.names = Terrestrial_Trait_rnames_2)

Terrestrial_Trait_study_2 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Morphology", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Tolerance", "Freq"]), 
                                        row.names = Terrestrial_Trait_rnames_2)

terrestrial_trait_summary_Reorder_2 <- terrestrial_trait_summary[c("Morphology", "Tolerance"), ]

Terrestrial_Trait_table_2 <- data.frame(estimate = terrestrial_trait_summary_Reorder_2[,"Mean"], 
                                        lowerCL = terrestrial_trait_summary_Reorder_2[,"Low_CI"], 
                                        upperCL = terrestrial_trait_summary_Reorder_2[,"High_CI"], 
                                        K = Terrestrial_Trait_k_2[,1], 
                                        group_no = Terrestrial_Trait_group_no_2[,1], 
                                        row.names = Terrestrial_Trait_rnames_2)
Terrestrial_Trait_table_2$name <- row.names(Terrestrial_Trait_table_2)

Terrestrial_Trait_raw_mean_2 <- c(b_abs_terrestrial_trait_morphology, b_abs_terrestrial_trait_tolerance)

Terrestrial_Trait_raw_name_2 <- c(replicate(8000, "Morphological"), 
                                  replicate(8000, "Tolerance"))

Terrestrial_Trait_raw_df_2 <- data.frame("Model" = Terrestrial_Trait_raw_name_2, 
                                         "Effect" = Terrestrial_Trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Trait_Order_2 <- c("Tolerance", "Morphological")

density_terrestrial_trait_2 <- Terrestrial_Trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = Terrestrial_Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_2)[1], 1)), xmin = max(Terrestrial_Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Trait_table_2)[1], 1)), xmin = min(Terrestrial_Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-2.7, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Trait_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Terrestrial_Trait_table_2["Tolerance", "K"], 
                                                             Terrestrial_Trait_table_2["Morphological", "K"]), "~","(", 
                                                           c(Terrestrial_Trait_table_2["Tolerance", "group_no"], 
                                                             Terrestrial_Trait_table_2["Morphological", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_terrestrial_trait_morphology)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Terrestrial_Trait_table_2$estimate+0.2), y = (seq(1, dim(Terrestrial_Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_2 #(400x240)

##### Terrestrial Model - Treatment Meta-regression #####
Terrestrial_Treatment_Exploration <- Terrestrial_Data %>% select("Type") %>% table() %>% data.frame()
Terrestrial_Treatment_Exploration <- Terrestrial_Treatment_Exploration %>% filter(Freq > 10)
rownames(Terrestrial_Treatment_Exploration) <- Terrestrial_Treatment_Exploration$Type

Terrestrial_Treatment_Data <- Terrestrial_Data %>% filter(Type == "Humidity"|
                                                          Type == "Temperature")

Terrestrial_Treatment_Species_Count <- Terrestrial_Treatment_Data %>% select("Scientific_Name", "Type") %>% 
                                       table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                       select("Type") %>% table() %>% data.frame()
rownames(Terrestrial_Treatment_Species_Count) <- Terrestrial_Treatment_Species_Count$Type

Terrestrial_Treatment_Study_Count <- Terrestrial_Treatment_Data %>% select("Study_ID", "Type") %>% 
                                     table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                     select("Type") %>% table() %>% data.frame()
rownames(Terrestrial_Treatment_Study_Count) <- Terrestrial_Treatment_Study_Count$Type

Terrestrial_Treatment_Species <- Terrestrial_Treatment_Data %>% select("phylo") %>% unique()

Terrestrial_Treatment_A <- as.data.frame(A)
Terrestrial_Treatment_A <- Terrestrial_Treatment_A[c(Terrestrial_Treatment_Species$phylo), c(Terrestrial_Treatment_Species$phylo)]
Terrestrial_Treatment_A <- as.matrix(Terrestrial_Treatment_A)

system.time(
  terrestrial_treatment <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                                     ~ Type + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Type, cor = FALSE)),
                                     data = Terrestrial_Treatment_Data,
                                     family = gaussian(),
                                     data2 = list(A = Terrestrial_Treatment_A), 
                                     chains = 4, 
                                     cores = 4,
                                     iter = 12000,
                                     warmup = 2000,
                                     thin = 5,
                                     prior = priors,
                                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                                     file = "./terrestrial_treatment_model",
                                     file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_terrestrial_treatment <- as_draws_df(terrestrial_treatment, variable = c("b_Intercept",
                                                                           "b_TypeTemperature"))
b_terrestrial_treatment <- data.frame("b_Humidity" = b_terrestrial_treatment$b_Intercept, 
                                      "b_Temperature" = b_terrestrial_treatment$b_TypeTemperature + b_terrestrial_treatment$b_Intercept)


sd_terrestrial_treatment <- as_draws_df(terrestrial_treatment, variable = c("sd_obs__Intercept:TypeHumidity", "sd_obs__Intercept:TypeTemperature", 
                                                                            "sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_terrestrial_treatment <- data.frame("sd_Humidity" = sd_terrestrial_treatment$`sd_obs__Intercept:TypeHumidity`, 
                                       "sd_Temperature" = sd_terrestrial_treatment$`sd_obs__Intercept:TypeTemperature`, 
                                       "sd_phylo__Intercept" = sd_terrestrial_treatment$`sd_phylo__Intercept`, 
                                       "sd_Study_ID__Intercept" = sd_terrestrial_treatment$`sd_Study_ID__Intercept`, 
                                       "sd_Measurement__Intercept" = sd_terrestrial_treatment$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
terrestrial_treatment_means <- apply(b_terrestrial_treatment, 2, mean)
terrestrial_treatment_cis <- apply(b_terrestrial_treatment, 2, function(x) HPDinterval(as.mcmc(x)))
terrestrial_treatment_pMCMC <- apply(b_terrestrial_treatment, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_terrestrial_treatment_humidity <- folded_norm(b_terrestrial_treatment$b_Humidity, sqrt(rowSums(sd_terrestrial_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_treatment[, "sd_Humidity"]^2)))
b_abs_terrestrial_treatment_temperature <- folded_norm(b_terrestrial_treatment$b_Temperature, sqrt(rowSums(sd_terrestrial_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_terrestrial_treatment[, "sd_Temperature"]^2)))
mean_abs_b_terrestrial_treatment_humidity <- mean(b_abs_terrestrial_treatment_humidity)
mean_abs_b_terrestrial_treatment_temperature <- mean(b_abs_terrestrial_treatment_temperature)
ci.abs_terrestrial_treatment_humidity <- HPDinterval(as.mcmc(b_abs_terrestrial_treatment_humidity))
ci.abs_terrestrial_treatment_temperature <- HPDinterval(as.mcmc(b_abs_terrestrial_treatment_temperature))

# Heterogeneity
terrestrial_i2_treatment <- i2(sd_terrestrial_treatment, Terrestrial_Treatment_Data$Variance_Adjusted) 

# Overall_trait_summary

terrestrial_treatment_means_list <- c(mean_abs_b_terrestrial_treatment_humidity, mean_abs_b_terrestrial_treatment_temperature)
terrestrial_treatment_low_ci <- c(ci.abs_terrestrial_treatment_humidity[1], ci.abs_terrestrial_treatment_temperature[1])
terrestrial_treatment_high_ci <- c(ci.abs_terrestrial_treatment_humidity[2], ci.abs_terrestrial_treatment_temperature[2])
terrestrial_treatment_categories <- c("Humidity", "Temperature")

terrestrial_treatment_summary <- matrix(c(terrestrial_treatment_means_list, terrestrial_treatment_low_ci, terrestrial_treatment_high_ci), 
                                    nrow = 2, ncol = 3, byrow = FALSE, 
                                    dimnames = list(c(terrestrial_treatment_categories), 
                                                    c("Mean", "Low_CI", "High_CI")))
terrestrial_treatment_summary <- data.frame(terrestrial_treatment_summary)

# Preparing Graph - Combined

Terrestrial_Treatment_rnames <- c("Humidity", "Temperature")

Terrestrial_Treatment_k <- data.frame("k" = c(Terrestrial_Treatment_Exploration["Humidity", "Freq"], 
                                              Terrestrial_Treatment_Exploration["Temperature", "Freq"]), 
                                      row.names = Terrestrial_Treatment_rnames)

Terrestrial_Treatment_group_no <- data.frame("Spp No." = c(Terrestrial_Treatment_Species_Count["Humidity", "Freq"], 
                                                           Terrestrial_Treatment_Species_Count["Temperature", "Freq"]), 
                                             row.names = Terrestrial_Treatment_rnames)

Terrestrial_Treatment_study <- data.frame("Study" = c(Terrestrial_Treatment_Study_Count["Humidity", "Freq"], 
                                                      Terrestrial_Treatment_Study_Count["Temperature", "Freq"]), 
                                          row.names = Terrestrial_Treatment_rnames)

Terrestrial_Treatment_table <- data.frame(estimate = terrestrial_treatment_summary[,"Mean"], 
                                          lowerCL = terrestrial_treatment_summary[,"Low_CI"], 
                                          upperCL = terrestrial_treatment_summary[,"High_CI"], 
                                          K = Terrestrial_Treatment_k[,1], 
                                          group_no = Terrestrial_Treatment_group_no[,1], 
                                          row.names = Terrestrial_Treatment_rnames)
Terrestrial_Treatment_table$name <- row.names(Terrestrial_Treatment_table)

Terrestrial_Treatment_raw_mean <- c(b_abs_terrestrial_treatment_humidity, b_abs_terrestrial_treatment_temperature)

Terrestrial_Treatment_raw_name <- c(replicate(8000, "Humidity"),  
                                    replicate(8000, "Temperature"))

Terrestrial_Treatment_raw_df <- data.frame("Model" = Terrestrial_Treatment_raw_name, 
                                           "Effect" = Terrestrial_Treatment_raw_mean)

# Graph code - Combined

Terrestrial_Treatment_Order <- c("Temperature", "Humidity")

density_terrestrial_treatment <- Terrestrial_Treatment_table %>% mutate(name = fct_relevel(name, Terrestrial_Treatment_Order)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = Terrestrial_Treatment_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Treatment_Order)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table)[1], 1)), xmin = max(Terrestrial_Treatment_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table)[1], 1)), xmin = min(Terrestrial_Treatment_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                                 theme_bw() +
                                 guides(fill = "none", colour = "none") +
                                 labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                  vjust = c(-2.7, -2.7))) +
                                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                 theme(axis.ticks = element_blank()) +
                                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                 scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                 scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                 coord_cartesian(xlim = c(-0.01, 1.25)) +
                                 annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Treatment_table)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(Terrestrial_Treatment_table["Temperature", "K"], 
                                                               Terrestrial_Treatment_table["Humidity", "K"]), "~","(", 
                                                             c(Terrestrial_Treatment_table["Temperature", "group_no"],
                                                               Terrestrial_Treatment_table["Humidity", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_treatment_temperature)-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(b_abs_terrestrial_treatment_humidity)-1)*100, 2), nsmall = 2), "%")), 
                                            x = rev(Terrestrial_Treatment_table$estimate+0.2), y = (seq(1, dim(Terrestrial_Treatment_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_treatment #(400x240)

# Preparing Graph - Part 1

Terrestrial_Treatment_rnames_1 <- c("Humidity")

Terrestrial_Treatment_k_1 <- data.frame("k" = c(Terrestrial_Treatment_Exploration["Humidity", "Freq"]), 
                                        row.names = Terrestrial_Treatment_rnames_1)

Terrestrial_Treatment_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Treatment_Species_Count["Humidity", "Freq"]), 
                                               row.names = Terrestrial_Treatment_rnames_1)

Terrestrial_Treatment_study_1 <- data.frame("Study" = c(Terrestrial_Treatment_Study_Count["Humidity", "Freq"]), 
                                            row.names = Terrestrial_Treatment_rnames_1)

terrestrial_treatment_summary_Reorder_1 <- terrestrial_treatment_summary[c("Humidity"), ]

Terrestrial_Treatment_table_1 <- data.frame(estimate = terrestrial_treatment_summary_Reorder_1[,"Mean"], 
                                            lowerCL = terrestrial_treatment_summary_Reorder_1[,"Low_CI"], 
                                            upperCL = terrestrial_treatment_summary_Reorder_1[,"High_CI"], 
                                            K = Terrestrial_Treatment_k_1[,1], 
                                            group_no = Terrestrial_Treatment_group_no_1[,1], 
                                            row.names = Terrestrial_Treatment_rnames_1)
Terrestrial_Treatment_table_1$name <- row.names(Terrestrial_Treatment_table_1)

Terrestrial_Treatment_raw_mean_1 <- c(b_abs_terrestrial_treatment_humidity)

Terrestrial_Treatment_raw_name_1 <- c(replicate(8000, "Humidity"))

Terrestrial_Treatment_raw_df_1 <- data.frame("Model" = Terrestrial_Treatment_raw_name_1, 
                                             "Effect" = Terrestrial_Treatment_raw_mean_1)

# Graph code - Part 1

Terrestrial_Treatment_Order_1 <- c("Humidity")

density_terrestrial_treatment_1 <- Terrestrial_Treatment_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Treatment_Order_1)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = Terrestrial_Treatment_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Treatment_Order_1)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.06, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)), xmin = max(Terrestrial_Treatment_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)), xmin = min(Terrestrial_Treatment_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#2B4E7A")) +
                                   scale_fill_manual(values = c("#2B4E7A")) +
                                   coord_cartesian(xlim = c(-0.01, 1.25)) +
                                   annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(Terrestrial_Treatment_table_1["Humidity", "K"]), "~","(", 
                                                               c(Terrestrial_Treatment_table_1["Humidity", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_treatment_humidity)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Terrestrial_Treatment_table_1$estimate+0.2), y = (seq(1, dim(Terrestrial_Treatment_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_treatment_1 #(400x160)

# Preparing Graph - Part 2

Terrestrial_Treatment_rnames_2 <- c("Temperature")

Terrestrial_Treatment_k_2 <- data.frame("k" = c(Terrestrial_Treatment_Exploration["Temperature", "Freq"]), 
                                        row.names = Terrestrial_Treatment_rnames_2)

Terrestrial_Treatment_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Treatment_Species_Count["Temperature", "Freq"]), 
                                               row.names = Terrestrial_Treatment_rnames_2)

Terrestrial_Treatment_study_2 <- data.frame("Study" = c(Terrestrial_Treatment_Study_Count["Temperature", "Freq"]), 
                                            row.names = Terrestrial_Treatment_rnames_2)

terrestrial_treatment_summary_Reorder_2 <- terrestrial_treatment_summary[c("Temperature"), ]

Terrestrial_Treatment_table_2 <- data.frame(estimate = terrestrial_treatment_summary_Reorder_2[,"Mean"], 
                                            lowerCL = terrestrial_treatment_summary_Reorder_2[,"Low_CI"], 
                                            upperCL = terrestrial_treatment_summary_Reorder_2[,"High_CI"], 
                                            K = Terrestrial_Treatment_k_2[,1], 
                                            group_no = Terrestrial_Treatment_group_no_2[,1], 
                                            row.names = Terrestrial_Treatment_rnames_2)
Terrestrial_Treatment_table_2$name <- row.names(Terrestrial_Treatment_table_2)

Terrestrial_Treatment_raw_mean_2 <- c(b_abs_terrestrial_treatment_temperature)

Terrestrial_Treatment_raw_name_2 <- c(replicate(8000, "Temperature"))

Terrestrial_Treatment_raw_df_2 <- data.frame("Model" = Terrestrial_Treatment_raw_name_2, 
                                             "Effect" = Terrestrial_Treatment_raw_mean_2)

# Graph code - Part 2

Terrestrial_Treatment_Order_2 <- c("Temperature")

density_terrestrial_treatment_2 <- Terrestrial_Treatment_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Treatment_Order_2)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = Terrestrial_Treatment_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Treatment_Order_2)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)), xmin = max(Terrestrial_Treatment_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)), xmin = min(Terrestrial_Treatment_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1")) +
                                   scale_fill_manual(values = c("#5D7AA1")) +
                                   coord_cartesian(xlim = c(-0.01, 1.25)) +
                                   annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(Terrestrial_Treatment_table_2["Temperature", "K"]), "~","(", 
                                                               c(Terrestrial_Treatment_table_2["Temperature", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_treatment_temperature)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Terrestrial_Treatment_table_2$estimate+0.2), y = (seq(1, dim(Terrestrial_Treatment_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_treatment_2 #(400x160)

##### Aquatic Model #####
Aquatic_Data <- data %>% filter(Ecosystem == "Aquatic")
Aquatic_Species <- Aquatic_Data %>% select("phylo") %>% unique()

Aquatic_A <- as.data.frame(A)
Aquatic_A <- Aquatic_A[c(Aquatic_Species$phylo), c(Aquatic_Species$phylo)]
Aquatic_A <- as.matrix(Aquatic_A)

system.time(
  aquatic <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                       ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|obs),
                       data = Aquatic_Data,
                       family = gaussian(),
                       data2 = list(A = Aquatic_A), 
                       chains = 4, 
                       cores = 4,
                       iter = 12000,
                       warmup = 2000,
                       thin = 5,
                       prior = priors,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       file = "./aquatic_model",
                       file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_aquatic <- as_draws_df(aquatic, variable = "b_Intercept")
b_aquatic <- data.frame(b_aquatic$b_Intercept)

sd_aquatic <- as_draws_df(aquatic, variable = c("sd_Measurement__Intercept", "sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic <- data.frame("sd_Measurement__Intercept" = sd_aquatic$sd_Measurement__Intercept, 
                         "sd_obs__Intercept" = sd_aquatic$sd_obs__Intercept, 
                         "sd_phylo__Intercept" = sd_aquatic$sd_phylo__Intercept, 
                         "sd_Study_ID__Intercept" = sd_aquatic$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic <-  sapply(b_aquatic, mean)
ci_b_aquatic <- as.vector(HPDinterval(as.mcmc(b_aquatic)))
pMCMC_b_aquatic <- 2*(1 - max(table(b_aquatic<0) / nrow(b_aquatic)))

# Absolute magnitude
b_abs_aquatic <- folded_norm(b_aquatic[,1], rowSums(sd_aquatic))
mean_abs_b_aquatic <- mean(b_abs_aquatic)
ci.abs_aquatic <- HPDinterval(as.mcmc(b_abs_aquatic))

# Heterogeneity
aquatic_i2 <- i2(sd_aquatic, Aquatic_Data$Variance_Adjusted) 

##### Aquatic Model - Plasticity Mechanism Meta-regression #####
Aquatic_Plasticity_Exploration <- Aquatic_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Exploration) <- Aquatic_Plasticity_Exploration$Plasticity_Category

Aquatic_Plasticity_Species_Count <- Aquatic_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                    table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                    select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Species_Count) <- Aquatic_Plasticity_Species_Count$Plasticity_Category

Aquatic_Plasticity_Study_Count <- Aquatic_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                  table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                  select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Study_Count) <- Aquatic_Plasticity_Study_Count$Plasticity_Category

system.time(
  aquatic_plastic <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                               ~ Plasticity_Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Plasticity_Category, cor = FALSE)),
                               data = Aquatic_Data,
                               family = gaussian(),
                               data2 = list(A = Aquatic_A), 
                               chains = 4, 
                               cores = 4,
                               iter = 12000,
                               warmup = 2000,
                               thin = 5,
                               prior = priors,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               file = "./aquatic_plastic_model",
                               file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_plastic <- as_draws_df(aquatic_plastic, variable = c("b_Intercept", 
                                                               "b_Plasticity_CategoryDevelopmental", 
                                                               "b_Plasticity_CategoryTransgenerational"))
b_aquatic_plastic <- data.frame("b_Acclimation" = b_aquatic_plastic$b_Intercept, 
                                "b_Developmental" = b_aquatic_plastic$b_Plasticity_CategoryDevelopmental + b_aquatic_plastic$b_Intercept, 
                                "b_Transgenerational" = b_aquatic_plastic$b_Plasticity_CategoryTransgenerational + b_aquatic_plastic$b_Intercept)


sd_aquatic_plastic <- as_draws_df(aquatic_plastic, variable = c("sd_obs__Intercept:Plasticity_CategoryAcclimation",
                                                                "sd_obs__Intercept:Plasticity_CategoryDevelopmental",
                                                                "sd_obs__Intercept:Plasticity_CategoryTransgenerational",
                                                                "sd_phylo__Intercept", "sd_Study_ID__Intercept", 
                                                                "sd_Measurement__Intercept"))
sd_aquatic_plastic <- data.frame("sd_Acclimation" = sd_aquatic_plastic$`sd_obs__Intercept:Plasticity_CategoryAcclimation`,
                                 "sd_Developmental" = sd_aquatic_plastic$`sd_obs__Intercept:Plasticity_CategoryDevelopmental`, 
                                 "sd_Transgenerational" = sd_aquatic_plastic$`sd_obs__Intercept:Plasticity_CategoryTransgenerational`,
                                 "sd_phylo__Intercept" = sd_aquatic_plastic$`sd_phylo__Intercept`, 
                                 "sd_Study_ID__Intercept" = sd_aquatic_plastic$`sd_Study_ID__Intercept`, 
                                 "sd_Measurement__Intercept" = sd_aquatic_plastic$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
aquatic_plastic_means <- apply(b_aquatic_plastic, 2, mean)
aquatic_plastic_cis <- apply(b_aquatic_plastic, 2, function(x) HPDinterval(as.mcmc(x)))
aquatic_plastic_pMCMC <- apply(b_aquatic_plastic, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_aquatic_plastic_acc <- folded_norm(b_aquatic_plastic$b_Acclimation, sqrt(rowSums(sd_aquatic_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_plastic[, "sd_Acclimation"]^2)))
b_abs_aquatic_plastic_dev <- folded_norm(b_aquatic_plastic$b_Developmental, sqrt(rowSums(sd_aquatic_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_plastic[, "sd_Developmental"]^2)))
b_abs_aquatic_plastic_trans <- folded_norm(b_aquatic_plastic$b_Transgenerational, sqrt(rowSums(sd_aquatic_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_plastic[, "sd_Transgenerational"]^2)))
mean_abs_b_aquatic_plastic_acc <- mean(b_abs_aquatic_plastic_acc)
mean_abs_b_aquatic_plastic_dev <- mean(b_abs_aquatic_plastic_dev)
mean_abs_b_aquatic_plastic_trans <- mean(b_abs_aquatic_plastic_trans)
ci.abs_aquatic_acc <- HPDinterval(as.mcmc(b_abs_aquatic_plastic_acc))
ci.abs_aquatic_dev <- HPDinterval(as.mcmc(b_abs_aquatic_plastic_dev))
ci.abs_aquatic_trans <- HPDinterval(as.mcmc(b_abs_aquatic_plastic_trans))

# Heterogeneity
aquatic_i2_plastic <- i2(sd_aquatic_plastic, Aquatic_Data$Variance_Adjusted) 

# Overall_Plasticity_Summary
aquatic_plastic_means_list <- c(mean_abs_b_aquatic_plastic_acc, mean_abs_b_aquatic_plastic_dev, 
                                mean_abs_b_aquatic_plastic_trans)
aquatic_plastic_low_ci <- c(ci.abs_aquatic_acc[1], ci.abs_aquatic_dev[1], 
                            ci.abs_aquatic_trans[1])
aquatic_plastic_high_ci <- c(ci.abs_aquatic_acc[2], ci.abs_aquatic_dev[2], 
                             ci.abs_aquatic_trans[2])
aquatic_plastic_categories <- c("Acclimation", "Developmental Plasticity", 
                                "Transgenerational Effects")

aquatic_plastic_summary <- matrix(c(aquatic_plastic_means_list, aquatic_plastic_low_ci, aquatic_plastic_high_ci), 
                                  nrow = 3, ncol = 3, byrow = FALSE, 
                                  dimnames = list(c(aquatic_plastic_categories), 
                                                  c("Mean", "Low_CI", "High_CI")))
aquatic_plastic_summary <- data.frame(aquatic_plastic_summary)

# Preparing Graph - Combined

Aquatic_Plasticity_rnames <- c("Acclimation", "Developmental Plasticity", "Transgenerational Effects")

Aquatic_Plasticity_k <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"], 
                                           Aquatic_Plasticity_Exploration["Developmental", "Freq"], 
                                           Aquatic_Plasticity_Exploration["Transgenerational", "Freq"]), 
                                   row.names = Aquatic_Plasticity_rnames)

Aquatic_Plasticity_group_no <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                        Aquatic_Plasticity_Species_Count["Developmental", "Freq"], 
                                                        Aquatic_Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                          row.names = Aquatic_Plasticity_rnames)

Aquatic_Plasticity_study <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                   Aquatic_Plasticity_Study_Count["Developmental", "Freq"], 
                                                   Aquatic_Plasticity_Study_Count["Transgenerational", "Freq"]), 
                                       row.names = Aquatic_Plasticity_rnames)

Aquatic_Plasticity_table <- data.frame(estimate = aquatic_plastic_summary[,"Mean"], 
                                       lowerCL = aquatic_plastic_summary[,"Low_CI"], 
                                       upperCL = aquatic_plastic_summary[,"High_CI"], 
                                       K = Aquatic_Plasticity_k[,1], 
                                       group_no = Aquatic_Plasticity_group_no[,1], 
                                       row.names = Aquatic_Plasticity_rnames)
Aquatic_Plasticity_table$name <- row.names(Aquatic_Plasticity_table)

Aquatic_Plasticity_raw_mean <- c(b_abs_aquatic_plastic_acc, b_abs_aquatic_plastic_dev, 
                                 b_abs_aquatic_plastic_trans)

Aquatic_Plasticity_raw_name <- c(replicate(8000, "Acclimation"), 
                                 replicate(8000, "Developmental Plasticity"), 
                                 replicate(8000, "Transgenerational Effects"))

Aquatic_Plasticity_raw_df <- data.frame("Model" = Aquatic_Plasticity_raw_name, 
                                        "Effect" = Aquatic_Plasticity_raw_mean)

# Graph code - Combined

Aquatic_Plasticity_Order <- c("Transgenerational Effects", "Developmental Plasticity", "Acclimation")

density_aquatic_plasticity <- Aquatic_Plasticity_table %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order)) %>%
                              ggplot() +
                              geom_density_ridges(data = Aquatic_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order)), 
                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                              geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table)[1], 1)), xmin = max(Aquatic_Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table)[1], 1)), xmin = min(Aquatic_Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                             size = 1) +
                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                              size = 1, fatten = 2) +
                              theme_bw() +
                              guides(fill = "none", colour = "none") +
                              labs(x = TeX("Effect Size (PRRD)"), y = "") +
                              theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                               vjust = c(-0.8, -0.8, -2.7))) +
                              theme(axis.text.x = element_text(margin = margin(b = 5))) +
                              theme(axis.ticks = element_blank()) +
                              theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                              theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                              scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                              scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                              scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                              coord_cartesian(xlim = c(-0.01, 1.25)) +
                              annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Plasticity_table)[1], 1)+0.4),
                              label= paste("italic(k)==", c(Aquatic_Plasticity_table["Transgenerational Effects", "K"], 
                                                            Aquatic_Plasticity_table["Developmental Plasticity", "K"], 
                                                            Aquatic_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                          c(Aquatic_Plasticity_table["Transgenerational Effects", "group_no"], 
                                                            Aquatic_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                            Aquatic_Plasticity_table["Acclimation", "group_no"]), 
                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                              geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_plastic_trans)-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(b_abs_aquatic_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(b_abs_aquatic_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Aquatic_Plasticity_table$estimate+0.2), y = (seq(1, dim(Aquatic_Plasticity_table)[1], 1)+0.4)), size = 3.5)
density_aquatic_plasticity #(400x320)

# Preparing Graph - Part 1

Aquatic_Plasticity_rnames_1 <- c("Acclimation", "Developmental Plasticity")

Aquatic_Plasticity_k_1 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"], 
                                             Aquatic_Plasticity_Exploration["Developmental", "Freq"]), 
                                     row.names = Aquatic_Plasticity_rnames_1)

Aquatic_Plasticity_group_no_1 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                          Aquatic_Plasticity_Species_Count["Developmental", "Freq"]), 
                                            row.names = Aquatic_Plasticity_rnames_1)

Aquatic_Plasticity_study_1 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                     Aquatic_Plasticity_Study_Count["Developmental", "Freq"]), 
                                         row.names = Aquatic_Plasticity_rnames_1)

aquatic_plastic_summary_Reorder_1 <- aquatic_plastic_summary[c("Acclimation", "Developmental Plasticity"), ]

Aquatic_Plasticity_table_1 <- data.frame(estimate = aquatic_plastic_summary_Reorder_1[,"Mean"], 
                                         lowerCL = aquatic_plastic_summary_Reorder_1[,"Low_CI"], 
                                         upperCL = aquatic_plastic_summary_Reorder_1[,"High_CI"], 
                                         K = Aquatic_Plasticity_k_1[,1], 
                                         group_no = Aquatic_Plasticity_group_no_1[,1], 
                                         row.names = Aquatic_Plasticity_rnames_1)
Aquatic_Plasticity_table_1$name <- row.names(Aquatic_Plasticity_table_1)

Aquatic_Plasticity_raw_mean_1 <- c(b_abs_aquatic_plastic_acc, b_abs_aquatic_plastic_dev)

Aquatic_Plasticity_raw_name_1 <- c(replicate(8000, "Acclimation"), 
                                   replicate(8000, "Developmental Plasticity"))

Aquatic_Plasticity_raw_df_1 <- data.frame("Model" = Aquatic_Plasticity_raw_name_1, 
                                          "Effect" = Aquatic_Plasticity_raw_mean_1)

# Graph code - Part 1

Aquatic_Plasticity_Order_1 <- c("Developmental Plasticity", "Acclimation")

density_aquatic_plasticity_1 <- Aquatic_Plasticity_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_1)) %>%
                                ggplot() +
                                geom_density_ridges(data = Aquatic_Plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_1)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)), xmin = max(Aquatic_Plasticity_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)), xmin = min(Aquatic_Plasticity_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                 vjust = c(-0.8, -2.7))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                coord_cartesian(xlim = c(-0.01, 1.25)) +
                                annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)+0.4),
                                label= paste("italic(k)==", c(Aquatic_Plasticity_table_1["Developmental Plasticity", "K"], 
                                                              Aquatic_Plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                            c(Aquatic_Plasticity_table_1["Developmental Plasticity", "group_no"], 
                                                              Aquatic_Plasticity_table_1["Acclimation", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(b_abs_aquatic_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                               x = rev(Aquatic_Plasticity_table_1$estimate+0.2), y = (seq(1, dim(Aquatic_Plasticity_table_1)[1], 1)+0.4)), size = 3.5)
density_aquatic_plasticity_1 #(400x240)

# Preparing Graph - Part 2

Aquatic_Plasticity_rnames_2 <- c("Transgenerational Effects")

Aquatic_Plasticity_k_2 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Transgenerational", "Freq"]), 
                                     row.names = Aquatic_Plasticity_rnames_2)

Aquatic_Plasticity_group_no_2 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                            row.names = Aquatic_Plasticity_rnames_2)

Aquatic_Plasticity_study_2 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Transgenerational", "Freq"]), 
                                         row.names = Aquatic_Plasticity_rnames_2)

aquatic_plastic_summary_Reorder_2 <- aquatic_plastic_summary[c("Transgenerational Effects"), ]

Aquatic_Plasticity_table_2 <- data.frame(estimate = aquatic_plastic_summary_Reorder_2[,"Mean"], 
                                         lowerCL = aquatic_plastic_summary_Reorder_2[,"Low_CI"], 
                                         upperCL = aquatic_plastic_summary_Reorder_2[,"High_CI"], 
                                         K = Aquatic_Plasticity_k_2[,1], 
                                         group_no = Aquatic_Plasticity_group_no_2[,1], 
                                         row.names = Aquatic_Plasticity_rnames_2)
Aquatic_Plasticity_table_2$name <- row.names(Aquatic_Plasticity_table_2)

Aquatic_Plasticity_raw_mean_2 <- c(b_abs_aquatic_plastic_trans)

Aquatic_Plasticity_raw_name_2 <- c(replicate(8000, "Transgenerational Effects"))

Aquatic_Plasticity_raw_df_2 <- data.frame("Model" = Aquatic_Plasticity_raw_name_2, 
                                          "Effect" = Aquatic_Plasticity_raw_mean_2)

# Graph code - Part 2

Aquatic_Plasticity_Order_2 <- c("Transgenerational Effects")

density_aquatic_plasticity_2 <- Aquatic_Plasticity_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_2)) %>%
                                ggplot() +
                                geom_density_ridges(data = Aquatic_Plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_2)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)), xmin = max(Aquatic_Plasticity_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)), xmin = min(Aquatic_Plasticity_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                 vjust = c(-0.8))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#5D7AA1")) +
                                scale_fill_manual(values = c("#5D7AA1")) +
                                coord_cartesian(xlim = c(-0.01, 1.25)) +
                                annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)+0.4),
                                label= paste("italic(k)==", c(Aquatic_Plasticity_table_2["Transgenerational Effects", "K"]), "~","(", 
                                                            c(Aquatic_Plasticity_table_2["Transgenerational Effects", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_plastic_trans)-1)*100, 2), nsmall = 2), "%")), 
                                               x = rev(Aquatic_Plasticity_table_2$estimate+0.2), y = (seq(1, dim(Aquatic_Plasticity_table_2)[1], 1)+0.4)), size = 3.5)
density_aquatic_plasticity_2 #(400x160)

##### Aquatic Model - Trait Category Meta-regression #####
Aquatic_Trait_Exploration <- Aquatic_Data %>% select("Category") %>% table() %>% data.frame()
Aquatic_Trait_Exploration <- Aquatic_Trait_Exploration %>% filter(Freq > 10)
rownames(Aquatic_Trait_Exploration) <- Aquatic_Trait_Exploration$Category

Aquatic_Trait_Data <- Aquatic_Data %>% filter(Category == "Biochemical Assay"| 
                                              Category == "Life-History Traits"| 
                                              Category == "Morphology"|
                                              Category == "Physiological"| 
                                              Category == "Tolerance")

Aquatic_Trait_Species_Count <- Aquatic_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                               table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                               select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Species_Count) <- Aquatic_Trait_Species_Count$Category

Aquatic_Trait_Study_Count <- Aquatic_Trait_Data %>% select("Study_ID", "Category") %>% 
                             table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                             select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Study_Count) <- Aquatic_Trait_Study_Count$Category

Aquatic_Trait_Species <- Aquatic_Trait_Data %>% select("phylo") %>% unique()

Aquatic_Trait_A <- as.data.frame(A)
Aquatic_Trait_A <- Aquatic_Trait_A[c(Aquatic_Trait_Species$phylo), c(Aquatic_Trait_Species$phylo)]
Aquatic_Trait_A <- as.matrix(Aquatic_Trait_A)

system.time(
  aquatic_trait <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                             ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Category, cor = FALSE)),
                             data = Aquatic_Trait_Data,
                             family = gaussian(),
                             data2 = list(A = Aquatic_Trait_A), 
                             chains = 4, 
                             cores = 4,
                             iter = 12000,
                             warmup = 2000,
                             thin = 5,
                             prior = priors,
                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                             file = "./aquatic_trait_model",
                             file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_trait <- as_draws_df(aquatic_trait, variable = c("b_Intercept", "b_CategoryLifeMHistoryTraits", 
                                                           "b_CategoryMorphology", "b_CategoryPhysiological",
                                                           "b_CategoryTolerance"))
b_aquatic_trait <- data.frame("b_BiochemicalAssay" = b_aquatic_trait$b_Intercept, 
                              "b_LifeMHistoryTraits" = b_aquatic_trait$b_CategoryLifeMHistoryTraits + b_aquatic_trait$b_Intercept, 
                              "b_Morphology" = b_aquatic_trait$b_CategoryMorphology + b_aquatic_trait$b_Intercept, 
                              "b_Physiological" = b_aquatic_trait$b_CategoryPhysiological + b_aquatic_trait$b_Intercept, 
                              "b_Tolerance" = b_aquatic_trait$b_CategoryTolerance + b_aquatic_trait$b_Intercept)


sd_aquatic_trait <- as_draws_df(aquatic_trait, variable = c("sd_obs__Intercept:CategoryBiochemicalAssay", "sd_obs__Intercept:CategoryLife-HistoryTraits", 
                                                            "sd_obs__Intercept:CategoryMorphology", "sd_obs__Intercept:CategoryPhysiological",
                                                            "sd_obs__Intercept:CategoryTolerance", "sd_phylo__Intercept", 
                                                            "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_aquatic_trait <- data.frame("sd_BiochemicalAssay" = sd_aquatic_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`,
                               "sd_Life-HistoryTraits" = sd_aquatic_trait$`sd_obs__Intercept:CategoryLife-HistoryTraits`, 
                               "sd_Morphology" = sd_aquatic_trait$`sd_obs__Intercept:CategoryMorphology`, 
                               "sd_Physiological" = sd_aquatic_trait$`sd_obs__Intercept:CategoryPhysiological`, 
                               "sd_Tolerance" = sd_aquatic_trait$`sd_obs__Intercept:CategoryTolerance`,
                               "sd_phylo__Intercept" = sd_aquatic_trait$`sd_phylo__Intercept`, 
                               "sd_Study_ID__Intercept" = sd_aquatic_trait$`sd_Study_ID__Intercept`, 
                               "sd_Measurement__Intercept" = sd_aquatic_trait$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
aquatic_trait_means <- apply(b_aquatic_trait, 2, mean)
aquatic_trait_cis <- apply(b_aquatic_trait, 2, function(x) HPDinterval(as.mcmc(x)))
aquatic_trait_pMCMC <- apply(b_aquatic_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_aquatic_trait_biochem <- folded_norm(b_aquatic_trait$b_BiochemicalAssay, sqrt(rowSums(sd_aquatic_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_aquatic_trait_life <- folded_norm(b_aquatic_trait$b_LifeMHistoryTraits, sqrt(rowSums(sd_aquatic_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_trait[, "sd_Life.HistoryTraits"]^2)))
b_abs_aquatic_trait_morphology <- folded_norm(b_aquatic_trait$b_Morphology, sqrt(rowSums(sd_aquatic_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_trait[, "sd_Morphology"]^2)))
b_abs_aquatic_trait_physiological <- folded_norm(b_aquatic_trait$b_Physiological, sqrt(rowSums(sd_aquatic_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_trait[, "sd_Physiological"]^2)))
b_abs_aquatic_trait_tolerance <- folded_norm(b_aquatic_trait$b_Tolerance, sqrt(rowSums(sd_aquatic_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_trait[, "sd_Tolerance"]^2)))
mean_abs_b_aquatic_trait_biochem <- mean(b_abs_aquatic_trait_biochem)
mean_abs_b_aquatic_trait_life <- mean(b_abs_aquatic_trait_life)
mean_abs_b_aquatic_trait_morphology <- mean(b_abs_aquatic_trait_morphology)
mean_abs_b_aquatic_trait_physiological <- mean(b_abs_aquatic_trait_physiological)
mean_abs_b_aquatic_trait_tolerance <- mean(b_abs_terrestrial_trait_tolerance)
ci.abs_aquatic_biochem <- HPDinterval(as.mcmc(b_abs_aquatic_trait_biochem))
ci.abs_aquatic_life <- HPDinterval(as.mcmc(b_abs_aquatic_trait_life))
ci.abs_aquatic_morphology <- HPDinterval(as.mcmc(b_abs_aquatic_trait_morphology))
ci.abs_aquatic_physiological <- HPDinterval(as.mcmc(b_abs_aquatic_trait_physiological))
ci.abs_aquatic_tolerance <- HPDinterval(as.mcmc(b_abs_terrestrial_trait_tolerance))

# Heterogeneity
aquatic_i2_trait <- i2(sd_aquatic_trait, Aquatic_Trait_Data$Variance_Adjusted) 

# Overall_trait_summary

aquatic_trait_means_list <- c(mean_abs_b_aquatic_trait_biochem, mean_abs_b_aquatic_trait_life, 
                              mean_abs_b_aquatic_trait_morphology, mean_abs_b_aquatic_trait_physiological, 
                              mean_abs_b_aquatic_trait_tolerance)
aquatic_trait_low_ci <- c(ci.abs_aquatic_biochem[1], ci.abs_aquatic_life[1], 
                          ci.abs_aquatic_morphology[1], ci.abs_aquatic_physiological[1], 
                          ci.abs_aquatic_tolerance[1])
aquatic_trait_high_ci <- c(ci.abs_aquatic_biochem[2], ci.abs_aquatic_life[2], 
                           ci.abs_aquatic_morphology[2], ci.abs_aquatic_physiological[2], 
                           ci.abs_aquatic_tolerance[2])
aquatic_trait_categories <- c("Biochemical Assay", "Life-History Traits", "Morphology",
                              "Physiological", "Tolerance")

aquatic_trait_summary <- matrix(c(aquatic_trait_means_list, aquatic_trait_low_ci, aquatic_trait_high_ci), 
                                  nrow = 5, ncol = 3, byrow = FALSE, 
                                  dimnames = list(c(aquatic_trait_categories), 
                                                  c("Mean", "Low_CI", "High_CI")))
aquatic_trait_summary <- data.frame(aquatic_trait_summary)

# Preparing Graph - Combined

Aquatic_Trait_rnames <- c("Biochemical Assay", "Life-history Traits", "Morphological",
                          "Physiological", "Tolerance")

Aquatic_Trait_k <- data.frame("k" = c(Aquatic_Trait_Exploration["Biochemical Assay", "Freq"], 
                                      Aquatic_Trait_Exploration["Life-History Traits", "Freq"], 
                                      Aquatic_Trait_Exploration["Morphology", "Freq"],
                                      Aquatic_Trait_Exploration["Physiological", "Freq"], 
                                      Aquatic_Trait_Exploration["Tolerance", "Freq"]), 
                              row.names = Aquatic_Trait_rnames)

Aquatic_Trait_group_no <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Life-History Traits", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Morphology", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Physiological", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Tolerance", "Freq"]), 
                                     row.names = Aquatic_Trait_rnames)

Aquatic_Trait_study <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                              Aquatic_Trait_Study_Count["Life-History Traits", "Freq"], 
                                              Aquatic_Trait_Study_Count["Morphology", "Freq"], 
                                              Aquatic_Trait_Study_Count["Physiological", "Freq"], 
                                              Aquatic_Trait_Study_Count["Tolerance", "Freq"]), 
                                  row.names = Aquatic_Trait_rnames)

aquatic_trait_summary_Reorder <- aquatic_trait_summary[c("Biochemical Assay", "Life-History Traits", "Morphology",
                                                         "Physiological", "Tolerance"), ]

Aquatic_Trait_table <- data.frame(estimate = aquatic_trait_summary_Reorder[,"Mean"], 
                                  lowerCL = aquatic_trait_summary_Reorder[,"Low_CI"], 
                                  upperCL = aquatic_trait_summary_Reorder[,"High_CI"], 
                                  K = Aquatic_Trait_k[,1], 
                                  group_no = Aquatic_Trait_group_no[,1], 
                                  row.names = Aquatic_Trait_rnames)
Aquatic_Trait_table$name <- row.names(Aquatic_Trait_table)

Aquatic_Trait_raw_mean <- c(b_abs_aquatic_trait_biochem, b_abs_aquatic_trait_life, b_abs_aquatic_trait_morphology, 
                            b_abs_aquatic_trait_physiological, b_abs_aquatic_trait_tolerance)

Aquatic_Trait_raw_name <- c(replicate(8000, "Biochemical Assay"), 
                            replicate(8000, "Life-history Traits"), 
                            replicate(8000, "Morphological"),
                            replicate(8000, "Physiological"), 
                            replicate(8000, "Tolerance"))

Aquatic_Trait_raw_df <- data.frame("Model" = Aquatic_Trait_raw_name, 
                                   "Effect" = Aquatic_Trait_raw_mean)

# Graph code - Combined

Aquatic_Trait_Order <- c("Tolerance", "Physiological",
                         "Morphological", "Life-history Traits", "Biochemical Assay")

density_aquatic_trait <- Aquatic_Trait_table %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order)) %>%
                         ggplot() +
                         geom_density_ridges(data = Aquatic_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order)), 
                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                         geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table)[1], 1)), xmin = max(Aquatic_Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table)[1], 1)), xmin = min(Aquatic_Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                        size = 1) +
                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                         size = 1, fatten = 2) +
                         theme_bw() +
                         guides(fill = "none", colour = "none") +
                         labs(x = TeX("Effect Size (PRRD)"), y = "") +
                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                          vjust = c(-2.7, -2.7, -2.7, -0.8, -0.8))) +
                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                         theme(axis.ticks = element_blank()) +
                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                         scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                         scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                         coord_cartesian(xlim = c(-0.01, 1.25)) +
                         annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Trait_table)[1], 1)+0.4),
                         label= paste("italic(k)==", c(Aquatic_Trait_table["Tolerance", "K"], 
                                                       Aquatic_Trait_table["Physiological", "K"],
                                                       Aquatic_Trait_table["Morphological", "K"], 
                                                       Aquatic_Trait_table["Life-history Traits", "K"], 
                                                       Aquatic_Trait_table["Biochemical Assay", "K"]), "~","(", 
                                                     c(Aquatic_Trait_table["Tolerance", "group_no"], 
                                                       Aquatic_Trait_table["Physiological", "group_no"],  
                                                       Aquatic_Trait_table["Morphological", "group_no"], 
                                                       Aquatic_Trait_table["Life-history Traits", "group_no"], 
                                                       Aquatic_Trait_table["Biochemical Assay", "group_no"]), 
                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(b_abs_aquatic_trait_physiological)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(b_abs_aquatic_trait_morphology)-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(b_abs_aquatic_trait_life)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(b_abs_aquatic_trait_biochem)-1)*100, 2), nsmall = 2), "%")), 
                                        x = rev(Aquatic_Trait_table$estimate+c(0.2, 0.2, 0.2, -0.2, 0.2)), y = (seq(1, dim(Aquatic_Trait_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait #(400x480)

# Preparing Graph - Part 1

Aquatic_Trait_rnames_1 <- c("Biochemical Assay", "Life-history Traits", "Morphological")

Aquatic_Trait_k_1 <- data.frame("k" = c(Aquatic_Trait_Exploration["Biochemical Assay", "Freq"], 
                                        Aquatic_Trait_Exploration["Life-History Traits", "Freq"], 
                                        Aquatic_Trait_Exploration["Morphology", "Freq"]), 
                                row.names = Aquatic_Trait_rnames_1)

Aquatic_Trait_group_no_1 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                     Aquatic_Trait_Species_Count["Life-History Traits", "Freq"], 
                                                     Aquatic_Trait_Species_Count["Morphology", "Freq"]), 
                                       row.names = Aquatic_Trait_rnames_1)

Aquatic_Trait_study_1 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                Aquatic_Trait_Study_Count["Life-History Traits", "Freq"], 
                                                Aquatic_Trait_Study_Count["Morphology", "Freq"]), 
                                    row.names = Aquatic_Trait_rnames_1)

aquatic_trait_summary_Reorder_1 <- aquatic_trait_summary[c("Biochemical Assay", "Life-History Traits", "Morphology"), ]

Aquatic_Trait_table_1 <- data.frame(estimate = aquatic_trait_summary_Reorder_1[,"Mean"], 
                                    lowerCL = aquatic_trait_summary_Reorder_1[,"Low_CI"], 
                                    upperCL = aquatic_trait_summary_Reorder_1[,"High_CI"], 
                                    K = Aquatic_Trait_k_1[,1], 
                                    group_no = Aquatic_Trait_group_no_1[,1], 
                                    row.names = Aquatic_Trait_rnames_1)
Aquatic_Trait_table_1$name <- row.names(Aquatic_Trait_table_1)

Aquatic_Trait_raw_mean_1 <- c(b_abs_aquatic_trait_biochem, b_abs_aquatic_trait_life, b_abs_aquatic_trait_morphology)

Aquatic_Trait_raw_name_1 <- c(replicate(8000, "Biochemical Assay"), 
                              replicate(8000, "Life-history Traits"), 
                              replicate(8000, "Morphological"))

Aquatic_Trait_raw_df_1 <- data.frame("Model" = Aquatic_Trait_raw_name_1, 
                                     "Effect" = Aquatic_Trait_raw_mean_1)

# Graph code - Part 1

Aquatic_Trait_Order_1 <- c("Morphological", "Life-history Traits", "Biochemical Assay")

density_aquatic_trait_1 <- Aquatic_Trait_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_1)) %>%
                           ggplot() +
                           geom_density_ridges(data = Aquatic_Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_1)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_1)[1], 1)), xmin = max(Aquatic_Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_1)[1], 1)), xmin = min(Aquatic_Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                           size = 1, fatten = 2) +
                           theme_bw() +
                           guides(fill = "none", colour = "none") +
                           labs(x = TeX("Effect Size (PRRD)"), y = "") +
                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                            vjust = c(-2.7, -0.8, -0.8))) +
                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                           theme(axis.ticks = element_blank()) +
                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                           scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                           scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                           coord_cartesian(xlim = c(-0.01, 1.25)) +
                           annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Trait_table_1)[1], 1)+0.4),
                           label= paste("italic(k)==", c(Aquatic_Trait_table_1["Morphological", "K"], 
                                                         Aquatic_Trait_table_1["Life-history Traits", "K"], 
                                                         Aquatic_Trait_table_1["Biochemical Assay", "K"]), "~","(", 
                                                       c(Aquatic_Trait_table_1["Morphological", "group_no"], 
                                                         Aquatic_Trait_table_1["Life-history Traits", "group_no"], 
                                                         Aquatic_Trait_table_1["Biochemical Assay", "group_no"]), 
                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_trait_morphology)-1)*100, 2), nsmall = 2), "%"),
                                                   paste(format(round(mean(exp(b_abs_aquatic_trait_life)-1)*100, 2), nsmall = 2), "%"), 
                                                   paste(format(round(mean(exp(b_abs_aquatic_trait_biochem)-1)*100, 2), nsmall = 2), "%")), 
                                         x = rev(Aquatic_Trait_table_1$estimate+c(0.2, 0.2, 0.2)), y = (seq(1, dim(Aquatic_Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait_1 #(400x320)

# Preparing Graph - Part 2

Aquatic_Trait_rnames_2 <- c("Physiological", "Tolerance")

Aquatic_Trait_k_2 <- data.frame("k" = c(Aquatic_Trait_Exploration["Physiological", "Freq"], 
                                        Aquatic_Trait_Exploration["Tolerance", "Freq"]), 
                                row.names = Aquatic_Trait_rnames_2)

Aquatic_Trait_group_no_2 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Physiological", "Freq"], 
                                                     Aquatic_Trait_Species_Count["Tolerance", "Freq"]), 
                                       row.names = Aquatic_Trait_rnames_2)

Aquatic_Trait_study_2 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Physiological", "Freq"], 
                                                Aquatic_Trait_Study_Count["Tolerance", "Freq"]), 
                                    row.names = Aquatic_Trait_rnames_2)

aquatic_trait_summary_Reorder_2 <- aquatic_trait_summary[c("Physiological", "Tolerance"), ]

Aquatic_Trait_table_2 <- data.frame(estimate = aquatic_trait_summary_Reorder_2[,"Mean"], 
                                    lowerCL = aquatic_trait_summary_Reorder_2[,"Low_CI"], 
                                    upperCL = aquatic_trait_summary_Reorder_2[,"High_CI"], 
                                    K = Aquatic_Trait_k_2[,1], 
                                    group_no = Aquatic_Trait_group_no_2[,1], 
                                    row.names = Aquatic_Trait_rnames_2)
Aquatic_Trait_table_2$name <- row.names(Aquatic_Trait_table_2)

Aquatic_Trait_raw_mean_2 <- c(b_abs_aquatic_trait_physiological, b_abs_aquatic_trait_tolerance)

Aquatic_Trait_raw_name_2 <- c(replicate(8000, "Physiological"), 
                              replicate(8000, "Tolerance"))

Aquatic_Trait_raw_df_2 <- data.frame("Model" = Aquatic_Trait_raw_name_2, 
                                     "Effect" = Aquatic_Trait_raw_mean_2)

# Graph code - Part 2

Aquatic_Trait_Order_2 <- c("Tolerance", "Physiological")

density_aquatic_trait_2 <- Aquatic_Trait_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_2)) %>%
                           ggplot() +
                           geom_density_ridges(data = Aquatic_Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_2)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_2)[1], 1)), xmin = max(Aquatic_Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Trait_table_2)[1], 1)), xmin = min(Aquatic_Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                           size = 1, fatten = 2) +
                           theme_bw() +
                           guides(fill = "none", colour = "none") +
                           labs(x = TeX("Effect Size (PRRD)"), y = "") +
                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                            vjust = c(-2.7, -2.7))) +
                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                           theme(axis.ticks = element_blank()) +
                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                           scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                           scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                           coord_cartesian(xlim = c(-0.01, 1.25)) +
                           annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Trait_table_2)[1], 1)+0.4),
                           label= paste("italic(k)==", c(Aquatic_Trait_table_2["Tolerance", "K"], 
                                                         Aquatic_Trait_table_2["Physiological", "K"]), "~","(", 
                                                       c(Aquatic_Trait_table_2["Tolerance", "group_no"], 
                                                         Aquatic_Trait_table_2["Physiological", "group_no"]), 
                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                   paste(format(round(mean(exp(b_abs_aquatic_trait_physiological)-1)*100, 2), nsmall = 2), "%")), 
                                         x = rev(Aquatic_Trait_table_2$estimate+c(-0.2, 0.2)), y = (seq(1, dim(Aquatic_Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait_2 #(400x240)

##### Aquatic Model - Treatment Meta-regression #####
Aquatic_Treatment_Exploration <- Aquatic_Data %>% select("Type") %>% table() %>% data.frame()
Aquatic_Treatment_Exploration <- Aquatic_Treatment_Exploration %>% filter(Freq > 10)
rownames(Aquatic_Treatment_Exploration) <- Aquatic_Treatment_Exploration$Type

Aquatic_Treatment_Data <- Aquatic_Data %>% filter(Type == "Diet"|
                                                  Type == "Predator"|
                                                  Type == "Salinity"|
                                                  Type == "Temperature"|
                                                  Type == "Water Level")

Aquatic_Treatment_Species_Count <- Aquatic_Treatment_Data %>% select("Scientific_Name", "Type") %>% 
                                   table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                   select("Type") %>% table() %>% data.frame()
rownames(Aquatic_Treatment_Species_Count) <- Aquatic_Treatment_Species_Count$Type

Aquatic_Treatment_Study_Count <- Aquatic_Treatment_Data %>% select("Study_ID", "Type") %>% 
                                 table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                 select("Type") %>% table() %>% data.frame()
rownames(Aquatic_Treatment_Study_Count) <- Aquatic_Treatment_Study_Count$Type

Aquatic_Treatment_Species <- Aquatic_Treatment_Data %>% select("phylo") %>% unique()

Aquatic_Treatment_A <- as.data.frame(A)
Aquatic_Treatment_A <- Aquatic_Treatment_A[c(Aquatic_Treatment_Species$phylo), c(Aquatic_Treatment_Species$phylo)]
Aquatic_Treatment_A <- as.matrix(Aquatic_Treatment_A)

system.time(
  aquatic_treatment <- brms::brm(Effect_Size_Adjusted | se(sqrt(Variance_Adjusted)) 
                                 ~ Type + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|Measurement) + (1|gr(obs, by = Type, cor = FALSE)),
                                 data = Aquatic_Treatment_Data,
                                 family = gaussian(),
                                 data2 = list(A = Aquatic_Treatment_A), 
                                 chains = 4, 
                                 cores = 4,
                                 iter = 12000,
                                 warmup = 2000,
                                 thin = 5,
                                 prior = priors,
                                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                                 file = "./aquatic_treatment_model",
                                 file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_treatment <- as_draws_df(aquatic_treatment, variable = c("b_Intercept",
                                                                   "b_TypePredator", 
                                                                   "b_TypeSalinity", 
                                                                   "b_TypeTemperature", 
                                                                   "b_TypeWaterLevel"))
b_aquatic_treatment <- data.frame("b_Diet" = b_aquatic_treatment$b_Intercept, 
                                  "b_Predator" = b_aquatic_treatment$b_TypePredator + b_aquatic_treatment$b_Intercept, 
                                  "b_Salinity" = b_aquatic_treatment$b_TypeSalinity + b_aquatic_treatment$b_Intercept, 
                                  "b_Temperature" = b_aquatic_treatment$b_TypeTemperature + b_aquatic_treatment$b_Intercept, 
                                  "b_WaterLevel" = b_aquatic_treatment$b_TypeWaterLevel + b_aquatic_treatment$b_Intercept)


sd_aquatic_treatment <- as_draws_df(aquatic_treatment, variable = c("sd_obs__Intercept:TypeDiet", "sd_obs__Intercept:TypePredator", 
                                                                    "sd_obs__Intercept:TypeSalinity", "sd_obs__Intercept:TypeTemperature",
                                                                    "sd_obs__Intercept:TypeWaterLevel",
                                                                    "sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept"))
sd_aquatic_treatment <- data.frame("sd_Diet" = sd_aquatic_treatment$`sd_obs__Intercept:TypeDiet`, 
                                   "sd_Predator" = sd_aquatic_treatment$`sd_obs__Intercept:TypePredator`,
                                   "sd_Salinity" = sd_aquatic_treatment$`sd_obs__Intercept:TypeSalinity`,
                                   "sd_Temperature" = sd_aquatic_treatment$`sd_obs__Intercept:TypeTemperature`,
                                   "sd_WaterLevel" = sd_aquatic_treatment$`sd_obs__Intercept:TypeWaterLevel`, 
                                   "sd_phylo__Intercept" = sd_aquatic_treatment$`sd_phylo__Intercept`, 
                                   "sd_Study_ID__Intercept" = sd_aquatic_treatment$`sd_Study_ID__Intercept`, 
                                   "sd_Measurement__Intercept" = sd_aquatic_treatment$`sd_Measurement__Intercept`)

# Overall estimates
# Signed
aquatic_treatment_means <- apply(b_aquatic_treatment, 2, mean)
aquatic_treatment_cis <- apply(b_aquatic_treatment, 2, function(x) HPDinterval(as.mcmc(x)))
aquatic_treatment_pMCMC <- apply(b_aquatic_treatment, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_aquatic_treatment_diet <- folded_norm(b_aquatic_treatment$b_Diet, sqrt(rowSums(sd_aquatic_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_treatment[, "sd_Diet"]^2)))
b_abs_aquatic_treatment_predator <- folded_norm(b_aquatic_treatment$b_Predator, sqrt(rowSums(sd_aquatic_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_treatment[, "sd_Predator"]^2)))
b_abs_aquatic_treatment_salinity <- folded_norm(b_aquatic_treatment$b_Salinity, sqrt(rowSums(sd_aquatic_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_treatment[, "sd_Salinity"]^2)))
b_abs_aquatic_treatment_temperature <- folded_norm(b_aquatic_treatment$b_Temperature, sqrt(rowSums(sd_aquatic_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_treatment[, "sd_Temperature"]^2)))
b_abs_aquatic_treatment_waterlevel <- folded_norm(b_aquatic_treatment$b_WaterLevel, sqrt(rowSums(sd_aquatic_treatment[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept", "sd_Measurement__Intercept")]^2 + sd_aquatic_treatment[, "sd_WaterLevel"]^2)))
mean_abs_b_aquatic_treatment_diet <- mean(b_abs_aquatic_treatment_diet)
mean_abs_b_aquatic_treatment_predator <- mean(b_abs_aquatic_treatment_predator)
mean_abs_b_aquatic_treatment_salinity <- mean(b_abs_aquatic_treatment_salinity)
mean_abs_b_aquatic_treatment_temperature <- mean(b_abs_aquatic_treatment_temperature)
mean_abs_b_aquatic_treatment_waterlevel <- mean(b_abs_aquatic_treatment_waterlevel)
ci.abs_aquatic_treatment_diet <- HPDinterval(as.mcmc(b_abs_aquatic_treatment_diet))
ci.abs_aquatic_treatment_predator <- HPDinterval(as.mcmc(b_abs_aquatic_treatment_predator))
ci.abs_aquatic_treatment_salinity <- HPDinterval(as.mcmc(b_abs_aquatic_treatment_salinity))
ci.abs_aquatic_treatment_temperature <- HPDinterval(as.mcmc(b_abs_aquatic_treatment_temperature))
ci.abs_aquatic_treatment_waterlevel <- HPDinterval(as.mcmc(b_abs_aquatic_treatment_waterlevel))

# Heterogeneity
aquatic_i2_treatment <- i2(sd_aquatic_treatment, Aquatic_Treatment_Data$Variance_Adjusted) 

# Overall_trait_summary

aquatic_treatment_means_list <- c(mean_abs_b_aquatic_treatment_diet, mean_abs_b_aquatic_treatment_predator,
                                  mean_abs_b_aquatic_treatment_salinity, mean_abs_b_aquatic_treatment_temperature,
                                  mean_abs_b_aquatic_treatment_waterlevel)
aquatic_treatment_low_ci <- c(ci.abs_aquatic_treatment_diet[1], ci.abs_aquatic_treatment_predator[1], 
                              ci.abs_aquatic_treatment_salinity[1], ci.abs_aquatic_treatment_temperature[1], 
                              ci.abs_aquatic_treatment_waterlevel[1])
aquatic_treatment_high_ci <- c(ci.abs_aquatic_treatment_diet[2], ci.abs_aquatic_treatment_predator[2], 
                               ci.abs_aquatic_treatment_salinity[2], ci.abs_aquatic_treatment_temperature[2], 
                               ci.abs_aquatic_treatment_waterlevel[2])
aquatic_treatment_categories <- c("Diet", "Predator", "Salinity", "Temperature", "Water Level")

aquatic_treatment_summary <- matrix(c(aquatic_treatment_means_list, aquatic_treatment_low_ci, aquatic_treatment_high_ci), 
                                    nrow = 5, ncol = 3, byrow = FALSE, 
                                    dimnames = list(c(aquatic_treatment_categories), 
                                                    c("Mean", "Low_CI", "High_CI")))
aquatic_treatment_summary <- data.frame(aquatic_treatment_summary)

# Preparing Graph - Combined

Aquatic_Treatment_rnames <- c("Diet", "Predator", "Salinity", "Temperature", "Water Level")

Aquatic_Treatment_k <- data.frame("k" = c(Aquatic_Treatment_Exploration["Diet", "Freq"],
                                          Aquatic_Treatment_Exploration["Predator", "Freq"],
                                          Aquatic_Treatment_Exploration["Salinity", "Freq"],
                                          Aquatic_Treatment_Exploration["Temperature", "Freq"],
                                          Aquatic_Treatment_Exploration["Water Level", "Freq"]), 
                                  row.names = Aquatic_Treatment_rnames)

Aquatic_Treatment_group_no <- data.frame("Spp No." = c(Aquatic_Treatment_Species_Count["Diet", "Freq"], 
                                                       Aquatic_Treatment_Species_Count["Predator", "Freq"],
                                                       Aquatic_Treatment_Species_Count["Salinity", "Freq"],
                                                       Aquatic_Treatment_Species_Count["Temperature", "Freq"],
                                                       Aquatic_Treatment_Species_Count["Water Level", "Freq"]), 
                                         row.names = Aquatic_Treatment_rnames)

Aquatic_Treatment_study <- data.frame("Study" = c(Aquatic_Treatment_Study_Count["Diet", "Freq"], 
                                                  Aquatic_Treatment_Study_Count["Predator", "Freq"], 
                                                  Aquatic_Treatment_Study_Count["Salinity", "Freq"], 
                                                  Aquatic_Treatment_Study_Count["Temperature", "Freq"], 
                                                  Aquatic_Treatment_Study_Count["Water Level", "Freq"]), 
                                      row.names = Aquatic_Treatment_rnames)

Aquatic_Treatment_table <- data.frame(estimate = aquatic_treatment_summary[,"Mean"], 
                                      lowerCL = aquatic_treatment_summary[,"Low_CI"], 
                                      upperCL = aquatic_treatment_summary[,"High_CI"], 
                                      K = Aquatic_Treatment_k[,1], 
                                      group_no = Aquatic_Treatment_group_no[,1], 
                                      row.names = Aquatic_Treatment_rnames)
Aquatic_Treatment_table$name <- row.names(Aquatic_Treatment_table)

Aquatic_Treatment_raw_mean <- c(b_abs_aquatic_treatment_diet, b_abs_aquatic_treatment_predator, 
                                b_abs_aquatic_treatment_salinity, b_abs_aquatic_treatment_temperature, 
                                b_abs_aquatic_treatment_waterlevel)

Aquatic_Treatment_raw_name <- c(replicate(8000, "Diet"), 
                                replicate(8000, "Predator"),
                                replicate(8000, "Salinity"),
                                replicate(8000, "Temperature"),
                                replicate(8000, "Water Level"))

Aquatic_Treatment_raw_df <- data.frame("Model" = Aquatic_Treatment_raw_name, 
                                       "Effect" = Aquatic_Treatment_raw_mean)

# Graph code - Combined

Aquatic_Treatment_Order <- c("Water Level", "Temperature", "Salinity", "Predator", "Diet")

density_aquatic_treatment <- Aquatic_Treatment_table %>% mutate(name = fct_relevel(name, Aquatic_Treatment_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = Aquatic_Treatment_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Treatment_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table)[1], 1)), xmin = max(Aquatic_Treatment_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table)[1], 1)), xmin = min(Aquatic_Treatment_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Treatment_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                             size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (PRRD)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                              vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             coord_cartesian(xlim = c(-0.01, 1.25)) +
                             annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Treatment_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(Aquatic_Treatment_table["Water Level", "K"], 
                                                           Aquatic_Treatment_table["Temperature", "K"],
                                                           Aquatic_Treatment_table["Salinity", "K"],
                                                           Aquatic_Treatment_table["Predator", "K"],
                                                           Aquatic_Treatment_table["Diet", "K"]), "~","(", 
                                                         c(Aquatic_Treatment_table["Water Level", "group_no"],
                                                           Aquatic_Treatment_table["Temperature", "group_no"],
                                                           Aquatic_Treatment_table["Salinity", "group_no"],
                                                           Aquatic_Treatment_table["Predator", "group_no"],
                                                           Aquatic_Treatment_table["Diet", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_treatment_waterlevel)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_aquatic_treatment_temperature)-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(b_abs_aquatic_treatment_salinity)-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(b_abs_aquatic_treatment_predator)-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(b_abs_aquatic_treatment_diet)-1)*100, 2), nsmall = 2), "%")), 
                                        x = rev(Aquatic_Treatment_table$estimate+0.2), y = (seq(1, dim(Aquatic_Treatment_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_treatment #(400x480)

# Preparing Graph - Part 1

Aquatic_Treatment_rnames_1 <- c("Diet", "Predator", "Salinity")

Aquatic_Treatment_k_1 <- data.frame("k" = c(Aquatic_Treatment_Exploration["Diet", "Freq"],
                                            Aquatic_Treatment_Exploration["Predator", "Freq"],
                                            Aquatic_Treatment_Exploration["Salinity", "Freq"]), 
                                    row.names = Aquatic_Treatment_rnames_1)

Aquatic_Treatment_group_no_1 <- data.frame("Spp No." = c(Aquatic_Treatment_Species_Count["Diet", "Freq"], 
                                                         Aquatic_Treatment_Species_Count["Predator", "Freq"],
                                                         Aquatic_Treatment_Species_Count["Salinity", "Freq"]), 
                                           row.names = Aquatic_Treatment_rnames_1)

Aquatic_Treatment_study_1 <- data.frame("Study" = c(Aquatic_Treatment_Study_Count["Diet", "Freq"], 
                                                    Aquatic_Treatment_Study_Count["Predator", "Freq"], 
                                                    Aquatic_Treatment_Study_Count["Salinity", "Freq"]), 
                                        row.names = Aquatic_Treatment_rnames_1)

aquatic_treatment_summary_Reorder_1 <- aquatic_treatment_summary[c("Diet", "Predator", "Salinity"), ]

Aquatic_Treatment_table_1 <- data.frame(estimate = aquatic_treatment_summary_Reorder_1[,"Mean"], 
                                        lowerCL = aquatic_treatment_summary_Reorder_1[,"Low_CI"], 
                                        upperCL = aquatic_treatment_summary_Reorder_1[,"High_CI"], 
                                        K = Aquatic_Treatment_k_1[,1], 
                                        group_no = Aquatic_Treatment_group_no_1[,1], 
                                        row.names = Aquatic_Treatment_rnames_1)
Aquatic_Treatment_table_1$name <- row.names(Aquatic_Treatment_table_1)

Aquatic_Treatment_raw_mean_1 <- c(b_abs_aquatic_treatment_diet, b_abs_aquatic_treatment_predator, 
                                  b_abs_aquatic_treatment_salinity)

Aquatic_Treatment_raw_name_1 <- c(replicate(8000, "Diet"), 
                                  replicate(8000, "Predator"),
                                  replicate(8000, "Salinity"))

Aquatic_Treatment_raw_df_1 <- data.frame("Model" = Aquatic_Treatment_raw_name_1, 
                                         "Effect" = Aquatic_Treatment_raw_mean_1)

# Graph code - Part 1

Aquatic_Treatment_Order_1 <- c("Salinity", "Predator", "Diet")

density_aquatic_treatment_1 <- Aquatic_Treatment_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Treatment_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = Aquatic_Treatment_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Treatment_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_1)[1], 1)), xmin = max(Aquatic_Treatment_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_1)[1], 1)), xmin = min(Aquatic_Treatment_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Treatment_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-2.7, -2.7, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Treatment_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Aquatic_Treatment_table_1["Salinity", "K"],
                                                             Aquatic_Treatment_table_1["Predator", "K"],
                                                             Aquatic_Treatment_table_1["Diet", "K"]), "~","(", 
                                                           c(Aquatic_Treatment_table_1["Salinity", "group_no"],
                                                             Aquatic_Treatment_table_1["Predator", "group_no"],
                                                             Aquatic_Treatment_table_1["Diet", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_treatment_salinity)-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(b_abs_aquatic_treatment_predator)-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(b_abs_aquatic_treatment_diet)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Aquatic_Treatment_table_1$estimate+0.2), y = (seq(1, dim(Aquatic_Treatment_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_treatment_1 #(400x320)

# Preparing Graph - Part 2

Aquatic_Treatment_rnames_2 <- c("Temperature", "Water Level")

Aquatic_Treatment_k_2 <- data.frame("k" = c(Aquatic_Treatment_Exploration["Temperature", "Freq"],
                                            Aquatic_Treatment_Exploration["Water Level", "Freq"]), 
                                    row.names = Aquatic_Treatment_rnames_2)

Aquatic_Treatment_group_no_2 <- data.frame("Spp No." = c(Aquatic_Treatment_Species_Count["Temperature", "Freq"],
                                                         Aquatic_Treatment_Species_Count["Water Level", "Freq"]), 
                                           row.names = Aquatic_Treatment_rnames_2)

Aquatic_Treatment_study_2 <- data.frame("Study" = c(Aquatic_Treatment_Study_Count["Temperature", "Freq"], 
                                                    Aquatic_Treatment_Study_Count["Water Level", "Freq"]), 
                                        row.names = Aquatic_Treatment_rnames_2)

aquatic_treatment_summary_Reorder_2 <- aquatic_treatment_summary[c("Temperature", "Water Level"), ]

Aquatic_Treatment_table_2 <- data.frame(estimate = aquatic_treatment_summary_Reorder_2[,"Mean"], 
                                        lowerCL = aquatic_treatment_summary_Reorder_2[,"Low_CI"], 
                                        upperCL = aquatic_treatment_summary_Reorder_2[,"High_CI"], 
                                        K = Aquatic_Treatment_k_2[,1], 
                                        group_no = Aquatic_Treatment_group_no_2[,1], 
                                        row.names = Aquatic_Treatment_rnames_2)
Aquatic_Treatment_table_2$name <- row.names(Aquatic_Treatment_table_2)

Aquatic_Treatment_raw_mean_2 <- c(b_abs_aquatic_treatment_temperature, b_abs_aquatic_treatment_waterlevel)

Aquatic_Treatment_raw_name_2 <- c(replicate(8000, "Temperature"),
                                  replicate(8000, "Water Level"))

Aquatic_Treatment_raw_df_2 <- data.frame("Model" = Aquatic_Treatment_raw_name_2, 
                                         "Effect" = Aquatic_Treatment_raw_mean_2)

# Graph code - Part 2

Aquatic_Treatment_Order_2 <- c("Water Level", "Temperature")

density_aquatic_treatment_2 <- Aquatic_Treatment_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Treatment_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = Aquatic_Treatment_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Treatment_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_2)[1], 1)), xmin = max(Aquatic_Treatment_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Treatment_table_2)[1], 1)), xmin = min(Aquatic_Treatment_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Treatment_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-2.7, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Treatment_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Aquatic_Treatment_table_2["Water Level", "K"], 
                                                             Aquatic_Treatment_table_2["Temperature", "K"]), "~","(", 
                                                           c(Aquatic_Treatment_table_2["Water Level", "group_no"],
                                                             Aquatic_Treatment_table_2["Temperature", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_treatment_waterlevel)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_aquatic_treatment_temperature)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Aquatic_Treatment_table_2$estimate+0.2), y = (seq(1, dim(Aquatic_Treatment_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_treatment_2 #(400x240)

##### Temperature Model #####
Temperature_Data <- data %>% filter(Type == "Temperature")
Temperature_Species <- Temperature_Data %>% select("phylo") %>% unique()

Temperature_A <- as.data.frame(A)
Temperature_A <- Temperature_A[c(Temperature_Species$phylo), c(Temperature_Species$phylo)]
Temperature_A <- as.matrix(Temperature_A)

system.time(
  temperature <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                       ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                       data = Temperature_Data,
                       family = gaussian(),
                       data2 = list(A = Temperature_A), 
                       chains = 4, 
                       cores = 4,
                       iter = 12000,
                       warmup = 2000,
                       thin = 5,
                       prior = priors,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       file = "./temperature_model",
                       file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_temperature <- as_draws_df(temperature, variable = "b_Intercept")
b_temperature <- data.frame(b_temperature$b_Intercept)

sd_temperature <- as_draws_df(temperature, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_temperature <- data.frame("sd_obs__Intercept" = sd_temperature$sd_obs__Intercept, 
                             "sd_phylo__Intercept" = sd_temperature$sd_phylo__Intercept, 
                             "sd_Study_ID__Intercept" = sd_temperature$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_temperature <-  sapply(b_temperature, mean)
ci_b_temperature <- as.vector(HPDinterval(as.mcmc(b_temperature)))
pMCMC_b_temperature <- 2*(1 - max(table(b_temperature<0) / nrow(b_temperature)))

# Absolute magnitude
b_abs_temperature <- folded_norm(b_temperature[,1], rowSums(sd_temperature))
mean_abs_b_temperature <- mean(b_abs_temperature)
ci.abs_temperature <- HPDinterval(as.mcmc(b_abs_temperature))

# Heterogeneity
temperature_i2 <- i2(sd_temperature, Temperature_Data$Variance_Type_Adjusted) 

##### Temperature Model - Plasticity Mechanism Meta-regression #####
Temperature_Plasticity_Exploration <- Temperature_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Temperature_Plasticity_Exploration) <- Temperature_Plasticity_Exploration$Plasticity_Category

Temperature_Plasticity_Species_Count <- Temperature_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                        table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                        select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Temperature_Plasticity_Species_Count) <- Temperature_Plasticity_Species_Count$Plasticity_Category

Temperature_Plasticity_Study_Count <- Temperature_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                      table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                      select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Temperature_Plasticity_Study_Count) <- Temperature_Plasticity_Study_Count$Plasticity_Category

system.time(
  temperature_plastic <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                   ~ Plasticity_Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Plasticity_Category, cor = FALSE)),
                                   data = Temperature_Data,
                                   family = gaussian(),
                                   data2 = list(A = Temperature_A), 
                                   chains = 4, 
                                   cores = 4,
                                   iter = 12000,
                                   warmup = 2000,
                                   thin = 5,
                                   prior = priors,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                                   file = "./temperature_plastic_model",
                                   file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_temperature_plastic <- as_draws_df(temperature_plastic, variable = c("b_Intercept", 
                                                                       "b_Plasticity_CategoryDevelopmental", 
                                                                       "b_Plasticity_CategoryTransgenerational"))
b_temperature_plastic <- data.frame("b_Acclimation" = b_temperature_plastic$b_Intercept, 
                                    "b_Developmental" = b_temperature_plastic$b_Plasticity_CategoryDevelopmental + b_temperature_plastic$b_Intercept, 
                                    "b_Transgenerational" = b_temperature_plastic$b_Plasticity_CategoryTransgenerational + b_temperature_plastic$b_Intercept)


sd_temperature_plastic <- as_draws_df(temperature_plastic, variable = c("sd_obs__Intercept:Plasticity_CategoryAcclimation",
                                                                        "sd_obs__Intercept:Plasticity_CategoryDevelopmental",
                                                                        "sd_obs__Intercept:Plasticity_CategoryTransgenerational",
                                                                        "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_temperature_plastic <- data.frame("sd_Acclimation" = sd_temperature_plastic$`sd_obs__Intercept:Plasticity_CategoryAcclimation`,
                                     "sd_Developmental" = sd_temperature_plastic$`sd_obs__Intercept:Plasticity_CategoryDevelopmental`, 
                                     "sd_Transgenerational" = sd_temperature_plastic$`sd_obs__Intercept:Plasticity_CategoryTransgenerational`,
                                     "sd_phylo__Intercept" = sd_temperature_plastic$`sd_phylo__Intercept`, 
                                     "sd_Study_ID__Intercept" = sd_temperature_plastic$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
temperature_plastic_means <- apply(b_temperature_plastic, 2, mean)
temperature_plastic_cis <- apply(b_temperature_plastic, 2, function(x) HPDinterval(as.mcmc(x)))
temperature_plastic_pMCMC <- apply(b_temperature_plastic, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_temperature_plastic_acc <- folded_norm(b_temperature_plastic$b_Acclimation, sqrt(rowSums(sd_temperature_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_plastic[, "sd_Acclimation"]^2)))
b_abs_temperature_plastic_dev <- folded_norm(b_temperature_plastic$b_Developmental, sqrt(rowSums(sd_temperature_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_plastic[, "sd_Developmental"]^2)))
b_abs_temperature_plastic_trans <- folded_norm(b_temperature_plastic$b_Transgenerational, sqrt(rowSums(sd_temperature_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_plastic[, "sd_Transgenerational"]^2)))
mean_abs_b_temperature_plastic_acc <- mean(b_abs_temperature_plastic_acc)
mean_abs_b_temperature_plastic_dev <- mean(b_abs_temperature_plastic_dev)
mean_abs_b_temperature_plastic_trans <- mean(b_abs_temperature_plastic_trans)
ci.abs_temperature_acc <- HPDinterval(as.mcmc(b_abs_temperature_plastic_acc))
ci.abs_temperature_dev <- HPDinterval(as.mcmc(b_abs_temperature_plastic_dev))
ci.abs_temperature_trans <- HPDinterval(as.mcmc(b_abs_temperature_plastic_trans))

# Heterogeneity
temperature_i2_plastic <- i2(sd_temperature_plastic, Temperature_Data$Variance_Type_Adjusted) 

# Overall_Plasticity_Summary
temperature_plastic_means_list <- c(mean_abs_b_temperature_plastic_acc, mean_abs_b_temperature_plastic_dev, 
                                    mean_abs_b_temperature_plastic_trans)
temperature_plastic_low_ci <- c(ci.abs_temperature_acc[1], ci.abs_temperature_dev[1], ci.abs_temperature_trans[1])
temperature_plastic_high_ci <- c(ci.abs_temperature_acc[2], ci.abs_temperature_dev[2], ci.abs_temperature_trans[2])
temperature_plastic_categories <- c("Acclimation", "Developmental Plasticity", "Transgenerational Effects")

temperature_plastic_summary <- matrix(c(temperature_plastic_means_list, temperature_plastic_low_ci, temperature_plastic_high_ci), 
                                      nrow = 3, ncol = 3, byrow = FALSE, 
                                      dimnames = list(c(temperature_plastic_categories), 
                                                      c("Mean", "Low_CI", "High_CI")))
temperature_plastic_summary <- data.frame(temperature_plastic_summary)

# Preparing Graph - Combined

Temperature_Plasticity_rnames <- c("Acclimation", "Developmental Plasticity", "Transgenerational Effects")

Temperature_Plasticity_k <- data.frame("k" = c(Temperature_Plasticity_Exploration["Acclimation", "Freq"], 
                                               Temperature_Plasticity_Exploration["Developmental", "Freq"], 
                                               Temperature_Plasticity_Exploration["Transgenerational", "Freq"]), 
                                       row.names = Temperature_Plasticity_rnames)

Temperature_Plasticity_group_no <- data.frame("Spp No." = c(Temperature_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                            Temperature_Plasticity_Species_Count["Developmental", "Freq"], 
                                                            Temperature_Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                              row.names = Temperature_Plasticity_rnames)

Temperature_Plasticity_study <- data.frame("Study" = c(Temperature_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                       Temperature_Plasticity_Study_Count["Developmental", "Freq"], 
                                                       Temperature_Plasticity_Study_Count["Transgenerational", "Freq"]), 
                                           row.names = Temperature_Plasticity_rnames)

Temperature_Plasticity_table <- data.frame(estimate = temperature_plastic_summary[,"Mean"], 
                                           lowerCL = temperature_plastic_summary[,"Low_CI"], 
                                           upperCL = temperature_plastic_summary[,"High_CI"], 
                                           K = Temperature_Plasticity_k[,1], 
                                           group_no = Temperature_Plasticity_group_no[,1], 
                                           row.names = Temperature_Plasticity_rnames)
Temperature_Plasticity_table$name <- row.names(Temperature_Plasticity_table)

Temperature_Plasticity_raw_mean <- c(b_abs_temperature_plastic_acc, b_abs_temperature_plastic_dev, b_abs_temperature_plastic_trans)

Temperature_Plasticity_raw_name <- c(replicate(8000, "Acclimation"), 
                                     replicate(8000, "Developmental Plasticity"), 
                                     replicate(8000, "Transgenerational Effects"))

Temperature_Plasticity_raw_df <- data.frame("Model" = Temperature_Plasticity_raw_name, 
                                            "Effect" = Temperature_Plasticity_raw_mean)

# Graph code - Combined

Temperature_Plasticity_Order <- c("Transgenerational Effects", "Developmental Plasticity", "Acclimation")

density_temperature_plasticity <- Temperature_Plasticity_table %>% mutate(name = fct_relevel(name, Temperature_Plasticity_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = Temperature_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Temperature_Plasticity_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table)[1], 1)), xmin = max(Temperature_Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table)[1], 1)), xmin = min(Temperature_Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-0.8, -0.8, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.01, 1.25)) +
                                  annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Plasticity_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(Temperature_Plasticity_table["Transgenerational Effects", "K"], 
                                                                Temperature_Plasticity_table["Developmental Plasticity", "K"], 
                                                                Temperature_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                              c(Temperature_Plasticity_table["Transgenerational Effects", "group_no"], 
                                                                Temperature_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                                Temperature_Plasticity_table["Acclimation", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_plastic_trans)-1)*100, 2), nsmall = 2), "%"), 
                                                         paste(format(round(mean(exp(b_abs_temperature_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(b_abs_temperature_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Temperature_Plasticity_table$estimate+0.2), y = (seq(1, dim(Temperature_Plasticity_table)[1], 1)+0.4)), size = 3.5)

density_temperature_plasticity #(400x320)

# Preparing Graph - Part 1

Temperature_Plasticity_rnames_1 <- c("Acclimation", "Developmental Plasticity")

Temperature_Plasticity_k_1 <- data.frame("k" = c(Temperature_Plasticity_Exploration["Acclimation", "Freq"], 
                                                 Temperature_Plasticity_Exploration["Developmental", "Freq"]), 
                                         row.names = Temperature_Plasticity_rnames_1)

Temperature_Plasticity_group_no_1 <- data.frame("Spp No." = c(Temperature_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                              Temperature_Plasticity_Species_Count["Developmental", "Freq"]), 
                                                row.names = Temperature_Plasticity_rnames_1)

Temperature_Plasticity_study_1 <- data.frame("Study" = c(Temperature_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                         Temperature_Plasticity_Study_Count["Developmental", "Freq"]), 
                                             row.names = Temperature_Plasticity_rnames_1)

temperature_plastic_summary_Reorder_1 <- temperature_plastic_summary[c("Acclimation", "Developmental Plasticity"), ]

Temperature_Plasticity_table_1 <- data.frame(estimate = temperature_plastic_summary_Reorder_1[,"Mean"], 
                                           lowerCL = temperature_plastic_summary_Reorder_1[,"Low_CI"], 
                                           upperCL = temperature_plastic_summary_Reorder_1[,"High_CI"], 
                                           K = Temperature_Plasticity_k_1[,1], 
                                           group_no = Temperature_Plasticity_group_no_1[,1], 
                                           row.names = Temperature_Plasticity_rnames_1)
Temperature_Plasticity_table_1$name <- row.names(Temperature_Plasticity_table_1)

Temperature_Plasticity_raw_mean_1 <- c(b_abs_temperature_plastic_acc, b_abs_temperature_plastic_dev)

Temperature_Plasticity_raw_name_1 <- c(replicate(8000, "Acclimation"), 
                                     replicate(8000, "Developmental Plasticity"))

Temperature_Plasticity_raw_df_1 <- data.frame("Model" = Temperature_Plasticity_raw_name_1, 
                                              "Effect" = Temperature_Plasticity_raw_mean_1)

# Graph code - Part 1

Temperature_Plasticity_Order_1 <- c("Developmental Plasticity", "Acclimation")

density_temperature_plasticity_1 <- Temperature_Plasticity_table_1 %>% mutate(name = fct_relevel(name, Temperature_Plasticity_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = Temperature_Plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Temperature_Plasticity_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_1)[1], 1)), xmin = max(Temperature_Plasticity_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_1)[1], 1)), xmin = min(Temperature_Plasticity_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                     vjust = c(-0.8, -2.7))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    coord_cartesian(xlim = c(-0.01, 1.25)) +
                                    annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Plasticity_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(Temperature_Plasticity_table_1["Developmental Plasticity", "K"], 
                                                                  Temperature_Plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                c(Temperature_Plasticity_table_1["Developmental Plasticity", "group_no"], 
                                                                  Temperature_Plasticity_table_1["Acclimation", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                           paste(format(round(mean(exp(b_abs_temperature_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                   x = rev(Temperature_Plasticity_table_1$estimate+0.2), y = (seq(1, dim(Temperature_Plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_temperature_plasticity_1 #(400x240)

# Preparing Graph - Part 2

Temperature_Plasticity_rnames_2 <- c("Transgenerational Effects")

Temperature_Plasticity_k_2 <- data.frame("k" = c(Temperature_Plasticity_Exploration["Transgenerational", "Freq"]), 
                                         row.names = Temperature_Plasticity_rnames_2)

Temperature_Plasticity_group_no_2 <- data.frame("Spp No." = c(Temperature_Plasticity_Species_Count["Transgenerational", "Freq"]), 
                                                row.names = Temperature_Plasticity_rnames_2)

Temperature_Plasticity_study_2 <- data.frame("Study" = c(Temperature_Plasticity_Study_Count["Transgenerational", "Freq"]), 
                                             row.names = Temperature_Plasticity_rnames_2)

temperature_plastic_summary_Reorder_2 <- temperature_plastic_summary[c("Transgenerational Effects"), ]

Temperature_Plasticity_table_2 <- data.frame(estimate = temperature_plastic_summary_Reorder_2[,"Mean"], 
                                           lowerCL = temperature_plastic_summary_Reorder_2[,"Low_CI"], 
                                           upperCL = temperature_plastic_summary_Reorder_2[,"High_CI"], 
                                           K = Temperature_Plasticity_k_2[,1], 
                                           group_no = Temperature_Plasticity_group_no_2[,1], 
                                           row.names = Temperature_Plasticity_rnames_2)
Temperature_Plasticity_table_2$name <- row.names(Temperature_Plasticity_table_2)

Temperature_Plasticity_raw_mean_2 <- c(b_abs_temperature_plastic_trans)

Temperature_Plasticity_raw_name_2 <- c(replicate(8000, "Transgenerational Effects"))

Temperature_Plasticity_raw_df_2 <- data.frame("Model" = Temperature_Plasticity_raw_name_2, 
                                              "Effect" = Temperature_Plasticity_raw_mean_2)

# Graph code - Part 2

Temperature_Plasticity_Order_2 <- c("Transgenerational Effects")

density_temperature_plasticity_2 <- Temperature_Plasticity_table_2 %>% mutate(name = fct_relevel(name, Temperature_Plasticity_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = Temperature_Plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Temperature_Plasticity_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_2)[1], 1)), xmin = max(Temperature_Plasticity_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(Temperature_Plasticity_table_2)[1], 1)), xmin = min(Temperature_Plasticity_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                     vjust = c(-0.8))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#5D7AA1")) +
                                    scale_fill_manual(values = c("#5D7AA1")) +
                                    coord_cartesian(xlim = c(-0.01, 1.25)) +
                                    annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Plasticity_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(Temperature_Plasticity_table_2["Transgenerational Effects", "K"]), "~","(", 
                                                                c(Temperature_Plasticity_table_2["Transgenerational Effects", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_plastic_trans)-1)*100, 2), nsmall = 2), "%")), 
                                                   x = rev(Temperature_Plasticity_table_2$estimate+0.2), y = (seq(1, dim(Temperature_Plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_temperature_plasticity_2 #(400x160)

##### Temperature Model - Trait Category Meta-regression #####
Temperature_Trait_Exploration <- Temperature_Data %>% select("Category") %>% table() %>% data.frame()
Temperature_Trait_Exploration <- Temperature_Trait_Exploration %>% filter(Freq > 10)
rownames(Temperature_Trait_Exploration) <- Temperature_Trait_Exploration$Category

Temperature_Trait_Data <- Temperature_Data %>% filter(Category == "Behavioural"| 
                                                      Category == "Biochemical Assay"| 
                                                      Category == "Morphology"|
                                                      Category == "Tolerance")

Temperature_Trait_Species_Count <- Temperature_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                                   table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                   select("Category") %>% table() %>% data.frame()
rownames(Temperature_Trait_Species_Count) <- Temperature_Trait_Species_Count$Category

Temperature_Trait_Study_Count <- Temperature_Trait_Data %>% select("Study_ID", "Category") %>% 
                                 table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                 select("Category") %>% table() %>% data.frame()
rownames(Temperature_Trait_Study_Count) <- Temperature_Trait_Study_Count$Category

Temperature_Trait_Species <- Temperature_Trait_Data %>% select("phylo") %>% unique()

Temperature_Trait_A <- as.data.frame(A)
Temperature_Trait_A <- Temperature_Trait_A[c(Temperature_Trait_Species$phylo), c(Temperature_Trait_Species$phylo)]
Temperature_Trait_A <- as.matrix(Temperature_Trait_A)

system.time(
  temperature_trait <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                 ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Category, cor = FALSE)),
                                 data = Temperature_Trait_Data,
                                 family = gaussian(),
                                 data2 = list(A = Temperature_Trait_A), 
                                 chains = 4, 
                                 cores = 4,
                                 iter = 12000,
                                 warmup = 2000,
                                 thin = 5,
                                 prior = priors,
                                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                                 file = "./temperature_trait_model",
                                 file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_temperature_trait <- as_draws_df(temperature_trait, variable = c("b_Intercept", "b_CategoryBiochemicalAssay", 
                                                                   "b_CategoryMorphology", "b_CategoryTolerance"))
b_temperature_trait <- data.frame("b_Behavioural" = b_temperature_trait$b_Intercept, 
                                  "b_BiochemicalAssay" = b_temperature_trait$b_CategoryBiochemicalAssay + b_temperature_trait$b_Intercept, 
                                  "b_Morphology" = b_temperature_trait$b_CategoryMorphology + b_temperature_trait$b_Intercept, 
                                  "b_Tolerance" = b_temperature_trait$b_CategoryTolerance + b_temperature_trait$b_Intercept)


sd_temperature_trait <- as_draws_df(temperature_trait, variable = c("sd_obs__Intercept:CategoryBehavioural", "sd_obs__Intercept:CategoryBiochemicalAssay", 
                                                                    "sd_obs__Intercept:CategoryMorphology", "sd_obs__Intercept:CategoryTolerance", 
                                                                    "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_temperature_trait <- data.frame("sd_Behavioural" = sd_temperature_trait$`sd_obs__Intercept:CategoryBehavioural`,
                                   "sd_BiochemicalAssay" = sd_temperature_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`, 
                                   "sd_Morphology" = sd_temperature_trait$`sd_obs__Intercept:CategoryMorphology`, 
                                   "sd_Tolerance" = sd_temperature_trait$`sd_obs__Intercept:CategoryTolerance`,
                                   "sd_phylo__Intercept" = sd_temperature_trait$`sd_phylo__Intercept`, 
                                   "sd_Study_ID__Intercept" = sd_temperature_trait$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
temperature_trait_means <- apply(b_temperature_trait, 2, mean)
temperature_trait_cis <- apply(b_temperature_trait, 2, function(x) HPDinterval(as.mcmc(x)))
temperature_trait_pMCMC <- apply(b_temperature_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_temperature_trait_behavioural <- folded_norm(b_temperature_trait$b_Behavioural, sqrt(rowSums(sd_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_trait[, "sd_Behavioural"]^2)))
b_abs_temperature_trait_biochem <- folded_norm(b_temperature_trait$b_BiochemicalAssay, sqrt(rowSums(sd_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_temperature_trait_morphology <- folded_norm(b_temperature_trait$b_Morphology, sqrt(rowSums(sd_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_trait[, "sd_Morphology"]^2)))
b_abs_temperature_trait_tolerance <- folded_norm(b_temperature_trait$b_Tolerance, sqrt(rowSums(sd_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_trait[, "sd_Tolerance"]^2)))
mean_abs_b_temperature_trait_behavioural <- mean(b_abs_temperature_trait_behavioural)
mean_abs_b_temperature_trait_biochem <- mean(b_abs_temperature_trait_biochem)
mean_abs_b_temperature_trait_morphology <- mean(b_abs_temperature_trait_morphology)
mean_abs_b_temperature_trait_tolerance <- mean(b_abs_temperature_trait_tolerance)
ci.abs_temperature_behavioural <- HPDinterval(as.mcmc(b_abs_temperature_trait_behavioural))
ci.abs_temperature_biochem <- HPDinterval(as.mcmc(b_abs_temperature_trait_biochem))
ci.abs_temperature_morphology <- HPDinterval(as.mcmc(b_abs_temperature_trait_morphology))
ci.abs_temperature_tolerance <- HPDinterval(as.mcmc(b_abs_temperature_trait_tolerance))

# Heterogeneity
temperature_i2_trait <- i2(sd_temperature_trait, Temperature_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

temperature_trait_means_list <- c(mean_abs_b_temperature_trait_behavioural, mean_abs_b_temperature_trait_biochem, 
                                  mean_abs_b_temperature_trait_morphology, mean_abs_b_temperature_trait_tolerance)
temperature_trait_low_ci <- c(ci.abs_temperature_behavioural[1], ci.abs_temperature_biochem[1], 
                              ci.abs_temperature_morphology[1], ci.abs_temperature_tolerance[1])
temperature_trait_high_ci <- c(ci.abs_temperature_behavioural[2], ci.abs_temperature_biochem[2], 
                               ci.abs_temperature_morphology[2], ci.abs_temperature_tolerance[2])
temperature_trait_categories <- c("Behavioural", "Biochemical Assay",
                                  "Morphology", "Tolerance")

temperature_trait_summary <- matrix(c(temperature_trait_means_list, temperature_trait_low_ci, temperature_trait_high_ci), 
                                    nrow = 4, ncol = 3, byrow = FALSE, 
                                    dimnames = list(c(temperature_trait_categories), 
                                                    c("Mean", "Low_CI", "High_CI")))
temperature_trait_summary <- data.frame(temperature_trait_summary)

# Preparing Graph - Combined

Temperature_Trait_rnames <- c("Behavioural", "Biochemical Assay",  
                              "Morphological", "Tolerance")

Temperature_Trait_k <- data.frame("k" = c(Temperature_Trait_Exploration["Behavioural", "Freq"], 
                                          Temperature_Trait_Exploration["Biochemical Assay", "Freq"], 
                                          Temperature_Trait_Exploration["Morphology", "Freq"], 
                                          Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                  row.names = Temperature_Trait_rnames)

Temperature_Trait_group_no <- data.frame("Spp No." = c(Temperature_Trait_Species_Count["Behavioural", "Freq"], 
                                                       Temperature_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Temperature_Trait_Species_Count["Morphology", "Freq"], 
                                                       Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                         row.names = Temperature_Trait_rnames)

Temperature_Trait_study <- data.frame("Study" = c(Temperature_Trait_Study_Count["Behavioural", "Freq"], 
                                                  Temperature_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Temperature_Trait_Study_Count["Morphology", "Freq"], 
                                                  Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                      row.names = Temperature_Trait_rnames)

temperature_trait_summary_Reorder <- temperature_trait_summary[c("Behavioural", "Biochemical Assay",  
                                                                 "Morphology", "Tolerance"), ]

Temperature_Trait_table <- data.frame(estimate = temperature_trait_summary_Reorder[,"Mean"], 
                                      lowerCL = temperature_trait_summary_Reorder[,"Low_CI"], 
                                      upperCL = temperature_trait_summary_Reorder[,"High_CI"], 
                                      K = Temperature_Trait_k[,1], 
                                      group_no = Temperature_Trait_group_no[,1], 
                                      row.names = Temperature_Trait_rnames)
Temperature_Trait_table$name <- row.names(Temperature_Trait_table)

Temperature_Trait_raw_mean <- c(b_abs_temperature_trait_behavioural, b_abs_temperature_trait_biochem, 
                                b_abs_temperature_trait_morphology, b_abs_temperature_trait_tolerance)

Temperature_Trait_raw_name <- c(replicate(8000, "Behavioural"), 
                                replicate(8000, "Biochemical Assay"), 
                                replicate(8000, "Morphological"), 
                                replicate(8000, "Tolerance"))

Temperature_Trait_raw_df <- data.frame("Model" = Temperature_Trait_raw_name, 
                                       "Effect" = Temperature_Trait_raw_mean)

# Graph code - Combined

Temperature_Trait_Order <- c("Tolerance", "Morphological", 
                             "Biochemical Assay", "Behavioural")

density_temperature_trait <- Temperature_Trait_table %>% mutate(name = fct_relevel(name, Temperature_Trait_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = Temperature_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Temperature_Trait_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table)[1], 1)), xmin = max(Temperature_Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table)[1], 1)), xmin = min(Temperature_Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                              vjust = c(-2.7, -2.7, -0.8, -2.7))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                             coord_cartesian(xlim = c(-0.01, 1.25)) +
                             annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Trait_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(Temperature_Trait_table["Tolerance", "K"], 
                                                           Temperature_Trait_table["Morphological", "K"], 
                                                           Temperature_Trait_table["Biochemical Assay", "K"], 
                                                           Temperature_Trait_table["Behavioural", "K"]), "~","(", 
                                                         c(Temperature_Trait_table["Tolerance", "group_no"], 
                                                           Temperature_Trait_table["Morphological", "group_no"], 
                                                           Temperature_Trait_table["Biochemical Assay", "group_no"], 
                                                           Temperature_Trait_table["Behavioural", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_temperature_trait_morphology)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(b_abs_temperature_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                        x = rev(Temperature_Trait_table$estimate+0.2), y = (seq(1, dim(Temperature_Trait_table)[1], 1)+0.4)), size = 3.5)

density_temperature_trait #(400x400)

# Preparing Graph - Part 1

Temperature_Trait_rnames_1 <- c("Behavioural", "Biochemical Assay")

Temperature_Trait_k_1 <- data.frame("k" = c(Temperature_Trait_Exploration["Behavioural", "Freq"], 
                                            Temperature_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                    row.names = Temperature_Trait_rnames_1)

Temperature_Trait_group_no_1 <- data.frame("Spp No." = c(Temperature_Trait_Species_Count["Behavioural", "Freq"], 
                                                         Temperature_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                           row.names = Temperature_Trait_rnames_1)

Temperature_Trait_study_1 <- data.frame("Study" = c(Temperature_Trait_Study_Count["Behavioural", "Freq"], 
                                                    Temperature_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                        row.names = Temperature_Trait_rnames_1)

temperature_trait_summary_Reorder_1 <- temperature_trait_summary[c("Behavioural", "Biochemical Assay"), ]

Temperature_Trait_table_1 <- data.frame(estimate = temperature_trait_summary_Reorder_1[,"Mean"], 
                                        lowerCL = temperature_trait_summary_Reorder_1[,"Low_CI"], 
                                        upperCL = temperature_trait_summary_Reorder_1[,"High_CI"], 
                                        K = Temperature_Trait_k_1[,1], 
                                        group_no = Temperature_Trait_group_no_1[,1], 
                                        row.names = Temperature_Trait_rnames_1)
Temperature_Trait_table_1$name <- row.names(Temperature_Trait_table_1)

Temperature_Trait_raw_mean_1 <- c(b_abs_temperature_trait_behavioural, b_abs_temperature_trait_biochem)

Temperature_Trait_raw_name_1 <- c(replicate(8000, "Behavioural"), 
                                  replicate(8000, "Biochemical Assay"))

Temperature_Trait_raw_df_1 <- data.frame("Model" = Temperature_Trait_raw_name_1, 
                                         "Effect" = Temperature_Trait_raw_mean_1)

# Graph code - Part 1

Temperature_Trait_Order_1 <- c("Biochemical Assay", "Behavioural")

density_temperature_trait_1 <- Temperature_Trait_table_1 %>% mutate(name = fct_relevel(name, Temperature_Trait_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = Temperature_Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Temperature_Trait_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_1)[1], 1)), xmin = max(Temperature_Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_1)[1], 1)), xmin = min(Temperature_Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-0.8, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                               scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Trait_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Temperature_Trait_table_1["Biochemical Assay", "K"], 
                                                             Temperature_Trait_table_1["Behavioural", "K"]), "~","(", 
                                                           c(Temperature_Trait_table_1["Biochemical Assay", "group_no"], 
                                                             Temperature_Trait_table_1["Behavioural", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_temperature_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Temperature_Trait_table_1$estimate+0.2), y = (seq(1, dim(Temperature_Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_temperature_trait_1 #(400x240)

# Preparing Graph - Part 2

Temperature_Trait_rnames_2 <- c("Morphological", "Tolerance")

Temperature_Trait_k_2 <- data.frame("k" = c(Temperature_Trait_Exploration["Morphology", "Freq"], 
                                            Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                    row.names = Temperature_Trait_rnames_2)

Temperature_Trait_group_no_2 <- data.frame("Spp No." = c(Temperature_Trait_Species_Count["Morphology", "Freq"], 
                                                         Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                           row.names = Temperature_Trait_rnames_2)

Temperature_Trait_study_2 <- data.frame("Study" = c(Temperature_Trait_Study_Count["Morphology", "Freq"], 
                                                    Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                        row.names = Temperature_Trait_rnames_2)

temperature_trait_summary_Reorder_2 <- temperature_trait_summary[c("Morphology", "Tolerance"), ]

Temperature_Trait_table_2 <- data.frame(estimate = temperature_trait_summary_Reorder_2[,"Mean"], 
                                      lowerCL = temperature_trait_summary_Reorder_2[,"Low_CI"], 
                                      upperCL = temperature_trait_summary_Reorder_2[,"High_CI"], 
                                      K = Temperature_Trait_k_2[,1], 
                                      group_no = Temperature_Trait_group_no_2[,1], 
                                      row.names = Temperature_Trait_rnames_2)
Temperature_Trait_table_2$name <- row.names(Temperature_Trait_table_2)

Temperature_Trait_raw_mean_2 <- c(b_abs_temperature_trait_morphology, b_abs_temperature_trait_tolerance)

Temperature_Trait_raw_name_2 <- c(replicate(8000, "Morphological"), 
                                  replicate(8000, "Tolerance"))

Temperature_Trait_raw_df_2 <- data.frame("Model" = Temperature_Trait_raw_name_2, 
                                         "Effect" = Temperature_Trait_raw_mean_2)

# Graph code - Part 2

Temperature_Trait_Order_2 <- c("Tolerance", "Morphological")

density_temperature_trait_2 <- Temperature_Trait_table_2 %>% mutate(name = fct_relevel(name, Temperature_Trait_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = Temperature_Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Temperature_Trait_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_2)[1], 1)), xmin = max(Temperature_Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(Temperature_Trait_table_2)[1], 1)), xmin = min(Temperature_Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-2.7, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                               coord_cartesian(xlim = c(-0.01, 1.25)) +
                               annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Trait_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(Temperature_Trait_table_2["Tolerance", "K"], 
                                                             Temperature_Trait_table_2["Morphological", "K"]), "~","(", 
                                                           c(Temperature_Trait_table_2["Tolerance", "group_no"], 
                                                             Temperature_Trait_table_2["Morphological", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_temperature_trait_morphology)-1)*100, 2), nsmall = 2), "%")), 
                                             x = rev(Temperature_Trait_table_2$estimate+0.2), y = (seq(1, dim(Temperature_Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_temperature_trait_2 #(400x240)

##### Temperature Model - Taxonomy Category Meta-regression #####
Temperature_Taxonomy_Exploration <- Temperature_Data %>% select("Class") %>% table() %>% data.frame()
Temperature_Taxonomy_Exploration <- Temperature_Taxonomy_Exploration %>% filter(Freq > 10)
rownames(Temperature_Taxonomy_Exploration) <- Temperature_Taxonomy_Exploration$Class

Temperature_Taxonomy_Data <- Temperature_Data %>% filter(Class == "Actinopteri"| 
                                                         Class == "Insecta")

Temperature_Taxonomy_Species_Count <- Temperature_Taxonomy_Data %>% select("Scientific_Name", "Class") %>% 
                                      table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                      select("Class") %>% table() %>% data.frame()
rownames(Temperature_Taxonomy_Species_Count) <- Temperature_Taxonomy_Species_Count$Class

Temperature_Taxonomy_Study_Count <- Temperature_Taxonomy_Data %>% select("Study_ID", "Class") %>% 
                                    table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                    select("Class") %>% table() %>% data.frame()
rownames(Temperature_Taxonomy_Study_Count) <- Temperature_Taxonomy_Study_Count$Class

Temperature_Taxonomy_Species <- Temperature_Taxonomy_Data %>% select("phylo") %>% unique()

Temperature_Taxonomy_A <- as.data.frame(A)
Temperature_Taxonomy_A <- Temperature_Taxonomy_A[c(Temperature_Taxonomy_Species$phylo), c(Temperature_Taxonomy_Species$phylo)]
Temperature_Taxonomy_A <- as.matrix(Temperature_Taxonomy_A)

system.time(
  temperature_taxonomy <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                    ~ Class + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Class, cor = FALSE)),
                                    data = Temperature_Taxonomy_Data,
                                    family = gaussian(),
                                    data2 = list(A = Temperature_Taxonomy_A), 
                                    chains = 4, 
                                    cores = 4,
                                    iter = 12000,
                                    warmup = 2000,
                                    thin = 5,
                                    prior = priors,
                                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                                    file = "./temperature_taxonomy_model",
                                    file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_temperature_taxonomy <- as_draws_df(temperature_taxonomy, variable = c("b_Intercept", "b_ClassInsecta"))
b_temperature_taxonomy <- data.frame("b_Actinopteri" = b_temperature_taxonomy$b_Intercept, 
                                     "b_Insecta" = b_temperature_taxonomy$b_ClassInsecta + b_temperature_taxonomy$b_Intercept)

sd_temperature_taxonomy <- as_draws_df(temperature_taxonomy, variable = c("sd_obs__Intercept:ClassActinopteri", "sd_obs__Intercept:ClassInsecta", 
                                                                          "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_temperature_taxonomy <- data.frame("sd_Actinopteri" = sd_temperature_taxonomy$`sd_obs__Intercept:ClassActinopteri`,
                                      "sd_Insecta" = sd_temperature_taxonomy$`sd_obs__Intercept:ClassInsecta`, 
                                      "sd_phylo__Intercept" = sd_temperature_taxonomy$`sd_phylo__Intercept`, 
                                      "sd_Study_ID__Intercept" = sd_temperature_taxonomy$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
temperature_taxonomy_means <- apply(b_temperature_taxonomy, 2, mean)
temperature_taxonomy_cis <- apply(b_temperature_taxonomy, 2, function(x) HPDinterval(as.mcmc(x)))
temperature_taxonomy_pMCMC <- apply(b_temperature_taxonomy, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_temperature_taxonomy_actinopteri <- folded_norm(b_temperature_taxonomy$b_Actinopteri, sqrt(rowSums(sd_temperature_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_taxonomy[, "sd_Actinopteri"]^2)))
b_abs_temperature_taxonomy_insecta <- folded_norm(b_temperature_taxonomy$b_Insecta, sqrt(rowSums(sd_temperature_taxonomy[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_temperature_taxonomy[, "sd_Insecta"]^2)))
mean_abs_b_temperature_taxonomy_actinopteri <- mean(b_abs_temperature_taxonomy_actinopteri)
mean_abs_b_temperature_taxonomy_insecta <- mean(b_abs_temperature_taxonomy_insecta)
ci.abs_temperature_taxonomy_actinopteri <- HPDinterval(as.mcmc(b_abs_temperature_taxonomy_actinopteri))
ci.abs_temperature_taxonomy_insecta <- HPDinterval(as.mcmc(b_abs_temperature_taxonomy_insecta))

# Heterogeneity
temperature_i2_taxonomy <- i2(sd_temperature_taxonomy, Temperature_Taxonomy_Data$Variance_Type_Adjusted) 

# Overall_trait_summary
temperature_taxonomy_means_list <- c(mean_abs_b_temperature_taxonomy_actinopteri, mean_abs_b_temperature_taxonomy_insecta)
temperature_taxonomy_low_ci <- c(ci.abs_temperature_taxonomy_actinopteri[1], ci.abs_temperature_taxonomy_insecta[1])
temperature_taxonomy_high_ci <- c(ci.abs_temperature_taxonomy_actinopteri[2], ci.abs_temperature_taxonomy_insecta[2])
temperature_taxonomy_categories <- c("Actinopteri", "Insecta")

temperature_taxonomy_summary <- matrix(c(temperature_taxonomy_means_list, temperature_taxonomy_low_ci, 
                                         temperature_taxonomy_high_ci), 
                                         nrow = 2, ncol = 3, byrow = FALSE, 
                                         dimnames = list(c(temperature_taxonomy_categories), 
                                                         c("Mean", "Low_CI", "High_CI")))
temperature_taxonomy_summary <- data.frame(temperature_taxonomy_summary)

# Preparing Graph - Combined

Temperature_Taxonomy_rnames <- c("Actinopteri", "Insecta")

Temperature_Taxonomy_k <- data.frame("k" = c(Temperature_Taxonomy_Exploration["Actinopteri", "Freq"], 
                                             Temperature_Taxonomy_Exploration["Insecta", "Freq"]), 
                                     row.names = Temperature_Taxonomy_rnames)

Temperature_Taxonomy_group_no <- data.frame("Spp No." = c(Temperature_Taxonomy_Species_Count["Actinopteri", "Freq"], 
                                                          Temperature_Taxonomy_Species_Count["Insecta", "Freq"]), 
                                            row.names = Temperature_Taxonomy_rnames)

Temperature_Taxonomy_study <- data.frame("Study" = c(Temperature_Taxonomy_Study_Count["Actinopteri", "Freq"], 
                                                     Temperature_Taxonomy_Study_Count["Insecta", "Freq"]), 
                                         row.names = Temperature_Taxonomy_rnames)

Temperature_Taxonomy_table <- data.frame(estimate = temperature_taxonomy_summary[,"Mean"], 
                                         lowerCL = temperature_taxonomy_summary[,"Low_CI"], 
                                         upperCL = temperature_taxonomy_summary[,"High_CI"], 
                                         K = Temperature_Taxonomy_k[,1], 
                                         group_no = Temperature_Taxonomy_group_no[,1], 
                                         row.names = Temperature_Taxonomy_rnames)
Temperature_Taxonomy_table$name <- row.names(Temperature_Taxonomy_table)

Temperature_Taxonomy_raw_mean <- c(b_abs_temperature_taxonomy_actinopteri, b_abs_temperature_taxonomy_insecta)

Temperature_Taxonomy_raw_name <- c(replicate(8000, "Actinopteri"), 
                                   replicate(8000, "Insecta"))

Temperature_Taxonomy_raw_df <- data.frame("Model" = Temperature_Taxonomy_raw_name, 
                                          "Effect" = Temperature_Taxonomy_raw_mean)

# Graph code - Combined

Temperature_Taxonomy_Order <- c("Insecta", "Actinopteri")

density_temperature_taxonomy <- Temperature_Taxonomy_table %>% mutate(name = fct_relevel(name, Temperature_Taxonomy_Order)) %>%
                                ggplot() +
                                geom_density_ridges(data = Temperature_Taxonomy_raw_df %>% mutate(Model = fct_relevel(Model, Temperature_Taxonomy_Order)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table)[1], 1)), xmin = max(Temperature_Taxonomy_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table)[1], 1)), xmin = min(Temperature_Taxonomy_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Taxonomy_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                 vjust = c(-2.7, -2.7))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                coord_cartesian(xlim = c(-0.01, 1.25)) +
                                annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Taxonomy_table)[1], 1)+0.4),
                                label= paste("italic(k)==", c(Temperature_Taxonomy_table["Insecta", "K"], 
                                                              Temperature_Taxonomy_table["Actinopteri", "K"]), "~","(", 
                                                            c(Temperature_Taxonomy_table["Insecta", "group_no"], 
                                                              Temperature_Taxonomy_table["Actinopteri", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_taxonomy_insecta)-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(b_abs_temperature_taxonomy_actinopteri)-1)*100, 2), nsmall = 2), "%")), 
                                               x = rev(Temperature_Taxonomy_table$estimate+0.2), y = (seq(1, dim(Temperature_Taxonomy_table)[1], 1)+0.4)), size = 3.5)

density_temperature_taxonomy #(400x240)

# Preparing Graph - Part 1

Temperature_Taxonomy_rnames_1 <- c("Actinopteri")

Temperature_Taxonomy_k_1 <- data.frame("k" = c(Temperature_Taxonomy_Exploration["Actinopteri", "Freq"]), 
                                       row.names = Temperature_Taxonomy_rnames_1)

Temperature_Taxonomy_group_no_1 <- data.frame("Spp No." = c(Temperature_Taxonomy_Species_Count["Actinopteri", "Freq"]), 
                                              row.names = Temperature_Taxonomy_rnames_1)

Temperature_Taxonomy_study_1 <- data.frame("Study" = c(Temperature_Taxonomy_Study_Count["Actinopteri", "Freq"]), 
                                           row.names = Temperature_Taxonomy_rnames_1)

temperature_taxonomy_summary_Reorder_1 <- temperature_taxonomy_summary[c("Actinopteri"), ]

Temperature_Taxonomy_table_1 <- data.frame(estimate = temperature_taxonomy_summary_Reorder_1[,"Mean"], 
                                           lowerCL = temperature_taxonomy_summary_Reorder_1[,"Low_CI"], 
                                           upperCL = temperature_taxonomy_summary_Reorder_1[,"High_CI"], 
                                           K = Temperature_Taxonomy_k_1[,1], 
                                           group_no = Temperature_Taxonomy_group_no_1[,1], 
                                           row.names = Temperature_Taxonomy_rnames_1)
Temperature_Taxonomy_table_1$name <- row.names(Temperature_Taxonomy_table_1)

Temperature_Taxonomy_raw_mean_1 <- c(b_abs_temperature_taxonomy_actinopteri)

Temperature_Taxonomy_raw_name_1 <- c(replicate(8000, "Actinopteri"))

Temperature_Taxonomy_raw_df_1 <- data.frame("Model" = Temperature_Taxonomy_raw_name_1, 
                                            "Effect" = Temperature_Taxonomy_raw_mean_1)

# Graph code - Part 1

Temperature_Taxonomy_Order_1 <- c("Actinopteri")

density_temperature_taxonomy_1 <- Temperature_Taxonomy_table_1 %>% mutate(name = fct_relevel(name, Temperature_Taxonomy_Order_1)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = Temperature_Taxonomy_raw_df_1 %>% mutate(Model = fct_relevel(Model, Temperature_Taxonomy_Order_1)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)), xmin = max(Temperature_Taxonomy_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)), xmin = min(Temperature_Taxonomy_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#2B4E7A")) +
                                  scale_fill_manual(values = c("#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.01, 1.25)) +
                                  annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(Temperature_Taxonomy_table_1["Actinopteri", "K"]), "~","(", 
                                                              c(Temperature_Taxonomy_table_1["Actinopteri", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_taxonomy_actinopteri)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Temperature_Taxonomy_table_1$estimate+0.2), y = (seq(1, dim(Temperature_Taxonomy_table_1)[1], 1)+0.4)), size = 3.5)

density_temperature_taxonomy_1 #(400x160)

# Preparing Graph - Part 2

Temperature_Taxonomy_rnames_2 <- c("Insecta")

Temperature_Taxonomy_k_2 <- data.frame("k" = c(Temperature_Taxonomy_Exploration["Insecta", "Freq"]), 
                                       row.names = Temperature_Taxonomy_rnames_2)

Temperature_Taxonomy_group_no_2 <- data.frame("Spp No." = c(Temperature_Taxonomy_Species_Count["Insecta", "Freq"]), 
                                              row.names = Temperature_Taxonomy_rnames_2)

Temperature_Taxonomy_study_2 <- data.frame("Study" = c(Temperature_Taxonomy_Study_Count["Insecta", "Freq"]), 
                                           row.names = Temperature_Taxonomy_rnames_2)

temperature_taxonomy_summary_Reorder_2 <- temperature_taxonomy_summary[c("Insecta"), ]

Temperature_Taxonomy_table_2 <- data.frame(estimate = temperature_taxonomy_summary_Reorder_2[,"Mean"], 
                                         lowerCL = temperature_taxonomy_summary_Reorder_2[,"Low_CI"], 
                                         upperCL = temperature_taxonomy_summary_Reorder_2[,"High_CI"], 
                                         K = Temperature_Taxonomy_k_2[,1], 
                                         group_no = Temperature_Taxonomy_group_no_2[,1], 
                                         row.names = Temperature_Taxonomy_rnames_2)
Temperature_Taxonomy_table_2$name <- row.names(Temperature_Taxonomy_table_2)

Temperature_Taxonomy_raw_mean_2 <- c(b_abs_temperature_taxonomy_insecta)

Temperature_Taxonomy_raw_name_2 <- c(replicate(8000, "Insecta"))

Temperature_Taxonomy_raw_df_2 <- data.frame("Model" = Temperature_Taxonomy_raw_name_2, 
                                            "Effect" = Temperature_Taxonomy_raw_mean_2)

# Graph code - Part 2

Temperature_Taxonomy_Order_2 <- c("Insecta")

density_temperature_taxonomy_2 <- Temperature_Taxonomy_table_2 %>% mutate(name = fct_relevel(name, Temperature_Taxonomy_Order_2)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = Temperature_Taxonomy_raw_df_2 %>% mutate(Model = fct_relevel(Model, Temperature_Taxonomy_Order_2)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)), xmin = max(Temperature_Taxonomy_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)), xmin = min(Temperature_Taxonomy_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1")) +
                                  scale_fill_manual(values = c("#5D7AA1")) +
                                  coord_cartesian(xlim = c(-0.01, 1.25)) +
                                  annotate('text',  x = 1.25, y = (seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(Temperature_Taxonomy_table_2["Insecta", "K"]), "~","(", 
                                                              c(Temperature_Taxonomy_table_2["Insecta", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_temperature_taxonomy_insecta)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Temperature_Taxonomy_table_2$estimate+0.2), y = (seq(1, dim(Temperature_Taxonomy_table_2)[1], 1)+0.4)), size = 3.5)

density_temperature_taxonomy_2 #(400x160)

##### Terrestrial and Temperature Model #####
Terrestrial_Temperature_Data <- data %>% filter(Ecosystem == "Terrestrial" & Type == "Temperature")
Terrestrial_Temperature_Species <- Terrestrial_Temperature_Data %>% select("phylo") %>% unique()

Terrestrial_Temperature_A <- as.data.frame(A)
Terrestrial_Temperature_A <- Terrestrial_Temperature_A[c(Terrestrial_Temperature_Species$phylo), c(Terrestrial_Temperature_Species$phylo)]
Terrestrial_Temperature_A <- as.matrix(Terrestrial_Temperature_A)

system.time(
  terrestrial_temperature <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                       ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                       data = Terrestrial_Temperature_Data,
                                       family = gaussian(),
                                       data2 = list(A = Terrestrial_Temperature_A), 
                                       chains = 4, 
                                       cores = 4,
                                       iter = 12000,
                                       warmup = 2000,
                                       thin = 5,
                                       prior = priors,
                                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                                       file = "./terrestrial_temperature_model",
                                       file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_terrestrial_temperature <- as_draws_df(terrestrial_temperature, variable = "b_Intercept")
b_terrestrial_temperature <- data.frame(b_terrestrial_temperature$b_Intercept)

sd_terrestrial_temperature <- as_draws_df(terrestrial_temperature, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_terrestrial_temperature <- data.frame("sd_obs__Intercept" = sd_terrestrial_temperature$sd_obs__Intercept, 
                                         "sd_phylo__Intercept" = sd_terrestrial_temperature$sd_phylo__Intercept, 
                                         "sd_Study_ID__Intercept" = sd_terrestrial_temperature$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_terrestrial_temperature <-  sapply(b_terrestrial_temperature, mean)
ci_b_terrestrial_temperature <- as.vector(HPDinterval(as.mcmc(b_terrestrial_temperature)))
pMCMC_b_terrestrial_temperature <- 2*(1 - max(table(b_terrestrial_temperature<0) / nrow(b_terrestrial_temperature)))

# Absolute magnitude
b_abs_terrestrial_temperature <- folded_norm(b_terrestrial_temperature[,1], rowSums(sd_terrestrial_temperature))
mean_abs_b_terrestrial_temperature <- mean(b_abs_terrestrial_temperature)
ci.abs_terrestrial_temperature <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature))

# Heterogeneity
terrestrial_temperature_i2 <- i2(sd_terrestrial_temperature, Terrestrial_Temperature_Data$Variance_Type_Adjusted) 

##### Terrestrial and Temperature Model - Plasticity Mechanism Meta-regression #####
Terrestrial_Temperature_Plasticity_Exploration <- Terrestrial_Temperature_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Temperature_Plasticity_Exploration) <- Terrestrial_Temperature_Plasticity_Exploration$Plasticity_Category

Terrestrial_Temperature_Plasticity_Data <- Terrestrial_Temperature_Data %>% filter(Plasticity_Category == "Acclimation"| 
                                                                                   Plasticity_Category == "Developmental")

Terrestrial_Temperature_Plasticity_Species_Count <- Terrestrial_Temperature_Plasticity_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                                    table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                                    select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Temperature_Plasticity_Species_Count) <- Terrestrial_Temperature_Plasticity_Species_Count$Plasticity_Category

Terrestrial_Temperature_Plasticity_Study_Count <- Terrestrial_Temperature_Plasticity_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                                  table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                                  select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Temperature_Plasticity_Study_Count) <- Terrestrial_Temperature_Plasticity_Study_Count$Plasticity_Category

Terrestrial_Temperature_Plasticity_Species <- Terrestrial_Temperature_Plasticity_Data %>% select("phylo") %>% unique()

Terrestrial_Temperature_Plasticity_A <- as.data.frame(A)
Terrestrial_Temperature_Plasticity_A <- Terrestrial_Temperature_Plasticity_A[c(Terrestrial_Temperature_Plasticity_Species$phylo), c(Terrestrial_Temperature_Plasticity_Species$phylo)]
Terrestrial_Temperature_Plasticity_A <- as.matrix(Terrestrial_Temperature_Plasticity_A)


system.time(
  terrestrial_temperature_plastic <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                               ~ Plasticity_Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Plasticity_Category, cor = FALSE)),
                                               data = Terrestrial_Temperature_Plasticity_Data,
                                               family = gaussian(),
                                               data2 = list(A = Terrestrial_Temperature_Plasticity_A), 
                                               chains = 4, 
                                               cores = 4,
                                               iter = 12000,
                                               warmup = 2000,
                                               thin = 5,
                                               prior = priors,
                                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                                               file = "./terrestrial_temperature_plastic_model",
                                               file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_terrestrial_temperature_plastic <- as_draws_df(terrestrial_temperature_plastic, variable = c("b_Intercept", 
                                                                                               "b_Plasticity_CategoryDevelopmental"))
b_terrestrial_temperature_plastic <- data.frame("b_Acclimation" = b_terrestrial_temperature_plastic$b_Intercept, 
                                                "b_Developmental" = b_terrestrial_temperature_plastic$b_Plasticity_CategoryDevelopmental + b_terrestrial_temperature_plastic$b_Intercept)


sd_terrestrial_temperature_plastic <- as_draws_df(terrestrial_temperature_plastic, variable = c("sd_obs__Intercept:Plasticity_CategoryAcclimation",
                                                                                                "sd_obs__Intercept:Plasticity_CategoryDevelopmental",
                                                                                                "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_terrestrial_temperature_plastic <- data.frame("sd_Acclimation" = sd_terrestrial_temperature_plastic$`sd_obs__Intercept:Plasticity_CategoryAcclimation`,
                                                 "sd_Developmental" = sd_terrestrial_temperature_plastic$`sd_obs__Intercept:Plasticity_CategoryDevelopmental`,
                                                 "sd_phylo__Intercept" = sd_terrestrial_temperature_plastic$`sd_phylo__Intercept`, 
                                                 "sd_Study_ID__Intercept" = sd_terrestrial_temperature_plastic$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
terrestrial_temperature_plastic_means <- apply(b_terrestrial_temperature_plastic, 2, mean)
terrestrial_temperature_plastic_cis <- apply(b_terrestrial_temperature_plastic, 2, function(x) HPDinterval(as.mcmc(x)))
terrestrial_temperature_plastic_pMCMC <- apply(b_terrestrial_temperature_plastic, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_terrestrial_temperature_plastic_acc <- folded_norm(b_terrestrial_temperature_plastic$b_Acclimation, sqrt(rowSums(sd_terrestrial_temperature_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_plastic[, "sd_Acclimation"]^2)))
b_abs_terrestrial_temperature_plastic_dev <- folded_norm(b_terrestrial_temperature_plastic$b_Developmental, sqrt(rowSums(sd_terrestrial_temperature_plastic[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_plastic[, "sd_Developmental"]^2)))
mean_abs_b_terrestrial_temperature_plastic_acc <- mean(b_abs_terrestrial_temperature_plastic_acc)
mean_abs_b_terrestrial_temperature_plastic_dev <- mean(b_abs_terrestrial_temperature_plastic_dev)
ci.abs_terrestrial_temperature_acc <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_plastic_acc))
ci.abs_terrestrial_temperature_dev <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_plastic_dev))

# Heterogeneity
terrestrial_temperature_i2_plastic <- i2(sd_terrestrial_temperature_plastic, Terrestrial_Temperature_Plasticity_Data$Variance_Type_Adjusted) 

# Overall_Plasticity_Summary
terrestrial_temperature_plastic_means_list <- c(mean_abs_b_terrestrial_temperature_plastic_acc, mean_abs_b_terrestrial_temperature_plastic_dev)
terrestrial_temperature_plastic_low_ci <- c(ci.abs_terrestrial_temperature_acc[1], ci.abs_terrestrial_temperature_dev[1])
terrestrial_temperature_plastic_high_ci <- c(ci.abs_terrestrial_temperature_acc[2], ci.abs_terrestrial_temperature_dev[2])
terrestrial_temperature_plastic_categories <- c("Acclimation", "Developmental Plasticity")

terrestrial_temperature_plastic_summary <- matrix(c(terrestrial_temperature_plastic_means_list, terrestrial_temperature_plastic_low_ci, terrestrial_temperature_plastic_high_ci), 
                                                  nrow = 2, ncol = 3, byrow = FALSE, 
                                                  dimnames = list(c(terrestrial_temperature_plastic_categories), 
                                                                  c("Mean", "Low_CI", "High_CI")))
terrestrial_temperature_plastic_summary <- data.frame(terrestrial_temperature_plastic_summary)

# Preparing Graph - Combined

Terrestrial_Temperature_Plasticity_rnames <- c("Acclimation", "Developmental Plasticity")

Terrestrial_Temperature_Plasticity_k <- data.frame("k" = c(Terrestrial_Temperature_Plasticity_Exploration["Acclimation", "Freq"], 
                                                           Terrestrial_Temperature_Plasticity_Exploration["Developmental", "Freq"]), 
                                                   row.names = Terrestrial_Temperature_Plasticity_rnames)

Terrestrial_Temperature_Plasticity_group_no <- data.frame("Spp No." = c(Terrestrial_Temperature_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                                        Terrestrial_Temperature_Plasticity_Species_Count["Developmental", "Freq"]), 
                                                          row.names = Terrestrial_Temperature_Plasticity_rnames)

Terrestrial_Temperature_Plasticity_study <- data.frame("Study" = c(Terrestrial_Temperature_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                                   Terrestrial_Temperature_Plasticity_Study_Count["Developmental", "Freq"]), 
                                                       row.names = Terrestrial_Temperature_Plasticity_rnames)

Terrestrial_Temperature_Plasticity_table <- data.frame(estimate = terrestrial_temperature_plastic_summary[,"Mean"], 
                                                       lowerCL = terrestrial_temperature_plastic_summary[,"Low_CI"], 
                                                       upperCL = terrestrial_temperature_plastic_summary[,"High_CI"], 
                                                       K = Terrestrial_Temperature_Plasticity_k[,1], 
                                                       group_no = Terrestrial_Temperature_Plasticity_group_no[,1], 
                                                       row.names = Terrestrial_Temperature_Plasticity_rnames)
Terrestrial_Temperature_Plasticity_table$name <- row.names(Terrestrial_Temperature_Plasticity_table)

Terrestrial_Temperature_Plasticity_raw_mean <- c(b_abs_terrestrial_temperature_plastic_acc, b_abs_terrestrial_temperature_plastic_dev)

Terrestrial_Temperature_Plasticity_raw_name <- c(replicate(8000, "Acclimation"), 
                                                 replicate(8000, "Developmental Plasticity"))

Terrestrial_Temperature_Plasticity_raw_df <- data.frame("Model" = Terrestrial_Temperature_Plasticity_raw_name, 
                                                        "Effect" = Terrestrial_Temperature_Plasticity_raw_mean)

# Graph code - Combined

Terrestrial_Temperature_Plasticity_Order <- c("Developmental Plasticity", "Acclimation")

density_terrestrial_temperature_plasticity <- Terrestrial_Temperature_Plasticity_table %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Plasticity_Order)) %>%
                                              ggplot() +
                                              geom_density_ridges(data = Terrestrial_Temperature_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Plasticity_Order)), 
                                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                              geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                                 size = 1) +
                                              geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)), xmin = max(Terrestrial_Temperature_Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                             size = 1) +
                                              geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)), xmin = min(Terrestrial_Temperature_Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                             size = 1) +
                                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                                  size = 1, fatten = 2) +
                                              theme_bw() +
                                              guides(fill = "none", colour = "none") +
                                              labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                              theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                               vjust = c(-0.8, -2.7))) +
                                              theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                              theme(axis.ticks = element_blank()) +
                                              theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                              scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                              scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                              coord_cartesian(xlim = c(-0.01, 1.25)) +
                                              annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)+0.4),
                                              label= paste("italic(k)==", c(Terrestrial_Temperature_Plasticity_table["Developmental Plasticity", "K"], 
                                                                            Terrestrial_Temperature_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                                          c(Terrestrial_Temperature_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                                            Terrestrial_Temperature_Plasticity_table["Acclimation", "group_no"]), 
                                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                              geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_plastic_dev)-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(b_abs_terrestrial_temperature_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                         x = rev(Terrestrial_Temperature_Plasticity_table$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_plasticity #(400x240)

# Preparing Graph - Part 1

Terrestrial_Temperature_Plasticity_rnames_1 <- c("Acclimation")

Terrestrial_Temperature_Plasticity_k_1 <- data.frame("k" = c(Terrestrial_Temperature_Plasticity_Exploration["Acclimation", "Freq"]), 
                                                     row.names = Terrestrial_Temperature_Plasticity_rnames_1)

Terrestrial_Temperature_Plasticity_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Temperature_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                            row.names = Terrestrial_Temperature_Plasticity_rnames_1)

Terrestrial_Temperature_Plasticity_study_1 <- data.frame("Study" = c(Terrestrial_Temperature_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                                         row.names = Terrestrial_Temperature_Plasticity_rnames_1)

terrestrial_temperature_plastic_summary_Reorder_1 <- terrestrial_temperature_plastic_summary[c("Acclimation"), ]

Terrestrial_Temperature_Plasticity_table_1 <- data.frame(estimate = terrestrial_temperature_plastic_summary_Reorder_1[,"Mean"], 
                                                       lowerCL = terrestrial_temperature_plastic_summary_Reorder_1[,"Low_CI"], 
                                                       upperCL = terrestrial_temperature_plastic_summary_Reorder_1[,"High_CI"], 
                                                       K = Terrestrial_Temperature_Plasticity_k_1[,1], 
                                                       group_no = Terrestrial_Temperature_Plasticity_group_no_1[,1], 
                                                       row.names = Terrestrial_Temperature_Plasticity_rnames_1)
Terrestrial_Temperature_Plasticity_table_1$name <- row.names(Terrestrial_Temperature_Plasticity_table_1)

Terrestrial_Temperature_Plasticity_raw_mean_1 <- c(b_abs_terrestrial_temperature_plastic_acc)

Terrestrial_Temperature_Plasticity_raw_name_1 <- c(replicate(8000, "Acclimation"))

Terrestrial_Temperature_Plasticity_raw_df_1 <- data.frame("Model" = Terrestrial_Temperature_Plasticity_raw_name_1, 
                                                          "Effect" = Terrestrial_Temperature_Plasticity_raw_mean_1)

# Graph code - Part 1

Terrestrial_Temperature_Plasticity_Order_1 <- c("Acclimation")

density_terrestrial_temperature_plasticity_1 <- Terrestrial_Temperature_Plasticity_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Plasticity_Order_1)) %>%
                                                ggplot() +
                                                geom_density_ridges(data = Terrestrial_Temperature_Plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Plasticity_Order_1)), 
                                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                    scale = 0.08, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                               size = 1) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)), xmin = max(Terrestrial_Temperature_Plasticity_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                               size = 1) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)), xmin = min(Terrestrial_Temperature_Plasticity_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                               size = 1) +
                                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                                size = 1, fatten = 2) +
                                                theme_bw() +
                                                guides(fill = "none", colour = "none") +
                                                labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                                 vjust = c(-2.7))) +
                                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                                theme(axis.ticks = element_blank()) +
                                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                                scale_colour_manual(values = c("#2B4E7A")) +
                                                scale_fill_manual(values = c("#2B4E7A")) +
                                                coord_cartesian(xlim = c(-0.01, 1.25)) +
                                                annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)+0.4),
                                                label= paste("italic(k)==", c(Terrestrial_Temperature_Plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                            c(Terrestrial_Temperature_Plasticity_table_1["Acclimation", "group_no"]), 
                                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_plastic_acc)-1)*100, 2), nsmall = 2), "%")), 
                                                               x = rev(Terrestrial_Temperature_Plasticity_table_1$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_plasticity_1 #(400x160)

# Preparing Graph - Part 2

Terrestrial_Temperature_Plasticity_rnames_2 <- c("Developmental Plasticity")

Terrestrial_Temperature_Plasticity_k_2 <- data.frame("k" = c(Terrestrial_Temperature_Plasticity_Exploration["Developmental", "Freq"]), 
                                                     row.names = Terrestrial_Temperature_Plasticity_rnames_2)

Terrestrial_Temperature_Plasticity_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Temperature_Plasticity_Species_Count["Developmental", "Freq"]), 
                                                            row.names = Terrestrial_Temperature_Plasticity_rnames_2)

Terrestrial_Temperature_Plasticity_study_2 <- data.frame("Study" = c(Terrestrial_Temperature_Plasticity_Study_Count["Developmental", "Freq"]), 
                                                         row.names = Terrestrial_Temperature_Plasticity_rnames_2)

terrestrial_temperature_plastic_summary_Reorder_2 <- terrestrial_temperature_plastic_summary[c("Developmental Plasticity"), ]

Terrestrial_Temperature_Plasticity_table_2 <- data.frame(estimate = terrestrial_temperature_plastic_summary_Reorder_2[,"Mean"], 
                                                         lowerCL = terrestrial_temperature_plastic_summary_Reorder_2[,"Low_CI"], 
                                                         upperCL = terrestrial_temperature_plastic_summary_Reorder_2[,"High_CI"], 
                                                         K = Terrestrial_Temperature_Plasticity_k_2[,1], 
                                                         group_no = Terrestrial_Temperature_Plasticity_group_no_2[,1], 
                                                         row.names = Terrestrial_Temperature_Plasticity_rnames_2)
Terrestrial_Temperature_Plasticity_table_2$name <- row.names(Terrestrial_Temperature_Plasticity_table_2)

Terrestrial_Temperature_Plasticity_raw_mean_2 <- c(b_abs_terrestrial_temperature_plastic_dev)

Terrestrial_Temperature_Plasticity_raw_name_2 <- c(replicate(8000, "Developmental Plasticity"))

Terrestrial_Temperature_Plasticity_raw_df_2 <- data.frame("Model" = Terrestrial_Temperature_Plasticity_raw_name_2, 
                                                          "Effect" = Terrestrial_Temperature_Plasticity_raw_mean_2)

# Graph code - Part 2

Terrestrial_Temperature_Plasticity_Order_2 <- c("Developmental Plasticity")

density_terrestrial_temperature_plasticity_2 <- Terrestrial_Temperature_Plasticity_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Plasticity_Order_2)) %>%
                                                ggplot() +
                                                geom_density_ridges(data = Terrestrial_Temperature_Plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Plasticity_Order_2)), 
                                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                    scale = 0.08, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                               size = 1) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)), xmin = max(Terrestrial_Temperature_Plasticity_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                               size = 1) +
                                                geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)), xmin = min(Terrestrial_Temperature_Plasticity_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                               size = 1) +
                                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                                size = 1, fatten = 2) +
                                                theme_bw() +
                                                guides(fill = "none", colour = "none") +
                                                labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                                 vjust = c(-0.8))) +
                                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                                theme(axis.ticks = element_blank()) +
                                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                                scale_colour_manual(values = c("#5D7AA1")) +
                                                scale_fill_manual(values = c("#5D7AA1")) +
                                                coord_cartesian(xlim = c(-0.01, 1.25)) +
                                                annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)+0.4),
                                                label= paste("italic(k)==", c(Terrestrial_Temperature_Plasticity_table_2["Developmental Plasticity", "K"]), "~","(", 
                                                                            c(Terrestrial_Temperature_Plasticity_table_2["Developmental Plasticity", "group_no"]), 
                                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                                geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_plastic_dev)-1)*100, 2), nsmall = 2), "%")), 
                                                               x = rev(Terrestrial_Temperature_Plasticity_table_2$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_plasticity_2 #(400x160)

##### Terrestrial and Temperature Model - Trait Category Meta-regression #####
Terrestrial_Temperature_Trait_Exploration <- Terrestrial_Temperature_Data %>% select("Category") %>% table() %>% data.frame()
Terrestrial_Temperature_Trait_Exploration <- Terrestrial_Temperature_Trait_Exploration %>% filter(Freq > 10)
rownames(Terrestrial_Temperature_Trait_Exploration) <- Terrestrial_Temperature_Trait_Exploration$Category

Terrestrial_Temperature_Trait_Data <- Terrestrial_Temperature_Data %>% filter(Category == "Behavioural"| 
                                                                              Category == "Biochemical Assay"|
                                                                              Category == "Morphology"| 
                                                                              Category == "Tolerance")

Terrestrial_Temperature_Trait_Species_Count <- Terrestrial_Temperature_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                                               table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                               select("Category") %>% table() %>% data.frame()
rownames(Terrestrial_Temperature_Trait_Species_Count) <- Terrestrial_Temperature_Trait_Species_Count$Category

Terrestrial_Temperature_Trait_Study_Count <- Terrestrial_Temperature_Trait_Data %>% select("Study_ID", "Category") %>% 
                                             table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                             select("Category") %>% table() %>% data.frame()
rownames(Terrestrial_Temperature_Trait_Study_Count) <- Terrestrial_Temperature_Trait_Study_Count$Category

Terrestrial_Temperature_Trait_Species <- Terrestrial_Temperature_Trait_Data %>% select("phylo") %>% unique()

Terrestrial_Temperature_Trait_A <- as.data.frame(A)
Terrestrial_Temperature_Trait_A <- Terrestrial_Temperature_Trait_A[c(Terrestrial_Temperature_Trait_Species$phylo), c(Terrestrial_Temperature_Trait_Species$phylo)]
Terrestrial_Temperature_Trait_A <- as.matrix(Terrestrial_Temperature_Trait_A)

system.time(
  terrestrial_temperature_trait <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                             ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Category, cor = FALSE)),
                                             data = Terrestrial_Temperature_Trait_Data,
                                             family = gaussian(),
                                             data2 = list(A = Terrestrial_Temperature_Trait_A), 
                                             chains = 4, 
                                             cores = 4,
                                             iter = 12000,
                                             warmup = 2000,
                                             thin = 5,
                                             prior = priors,
                                             control = list(adapt_delta = 0.99, max_treedepth = 15),
                                             file = "./terrestrial_temperature_trait_model",
                                             file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_terrestrial_temperature_trait <- as_draws_df(terrestrial_temperature_trait, variable = c("b_Intercept", "b_CategoryBiochemicalAssay", 
                                                                                           "b_CategoryMorphology", "b_CategoryTolerance"))
b_terrestrial_temperature_trait <- data.frame("b_Behavioural" = b_terrestrial_temperature_trait$b_Intercept, 
                                              "b_BiochemicalAssay" = b_terrestrial_temperature_trait$b_CategoryBiochemicalAssay + b_terrestrial_temperature_trait$b_Intercept, 
                                              "b_Morphology" = b_terrestrial_temperature_trait$b_CategoryMorphology + b_terrestrial_temperature_trait$b_Intercept, 
                                              "b_Tolerance" = b_terrestrial_temperature_trait$b_CategoryTolerance + b_terrestrial_temperature_trait$b_Intercept)


sd_terrestrial_temperature_trait <- as_draws_df(terrestrial_temperature_trait, variable = c("sd_obs__Intercept:CategoryBehavioural", "sd_obs__Intercept:CategoryBiochemicalAssay", 
                                                                                            "sd_obs__Intercept:CategoryMorphology", "sd_obs__Intercept:CategoryTolerance", 
                                                                                            "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_terrestrial_temperature_trait <- data.frame("sd_Behavioural" = sd_terrestrial_temperature_trait$`sd_obs__Intercept:CategoryBehavioural`,
                                               "sd_BiochemicalAssay" = sd_terrestrial_temperature_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`, 
                                               "sd_Morphology" = sd_terrestrial_temperature_trait$`sd_obs__Intercept:CategoryMorphology`, 
                                               "sd_Tolerance" = sd_terrestrial_temperature_trait$`sd_obs__Intercept:CategoryTolerance`,
                                               "sd_phylo__Intercept" = sd_terrestrial_temperature_trait$`sd_phylo__Intercept`, 
                                               "sd_Study_ID__Intercept" = sd_terrestrial_temperature_trait$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
terrestrial_temperature_trait_means <- apply(b_terrestrial_temperature_trait, 2, mean)
terrestrial_temperature_trait_cis <- apply(b_terrestrial_temperature_trait, 2, function(x) HPDinterval(as.mcmc(x)))
terrestrial_temperature_trait_pMCMC <- apply(b_terrestrial_temperature_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_terrestrial_temperature_trait_behavioural <- folded_norm(b_terrestrial_temperature_trait$b_Behavioural, sqrt(rowSums(sd_terrestrial_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_trait[, "sd_Behavioural"]^2)))
b_abs_terrestrial_temperature_trait_biochem <- folded_norm(b_terrestrial_temperature_trait$b_BiochemicalAssay, sqrt(rowSums(sd_terrestrial_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_terrestrial_temperature_trait_morphology <- folded_norm(b_terrestrial_temperature_trait$b_Morphology, sqrt(rowSums(sd_terrestrial_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_trait[, "sd_Morphology"]^2)))
b_abs_terrestrial_temperature_trait_tolerance <- folded_norm(b_terrestrial_temperature_trait$b_Tolerance, sqrt(rowSums(sd_terrestrial_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_terrestrial_temperature_trait[, "sd_Tolerance"]^2)))
mean_abs_b_terrestrial_temperature_trait_behavioural <- mean(b_abs_terrestrial_temperature_trait_behavioural)
mean_abs_b_terrestrial_temperature_trait_biochem <- mean(b_abs_terrestrial_temperature_trait_biochem)
mean_abs_b_terrestrial_temperature_trait_morphology <- mean(b_abs_terrestrial_temperature_trait_morphology)
mean_abs_b_terrestrial_temperature_trait_tolerance <- mean(b_abs_terrestrial_temperature_trait_tolerance)
ci.abs_terrestrial_temperature_behavioural <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_trait_behavioural))
ci.abs_terrestrial_temperature_biochem <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_trait_biochem))
ci.abs_terrestrial_temperature_morphology <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_trait_morphology))
ci.abs_terrestrial_temperature_tolerance <- HPDinterval(as.mcmc(b_abs_terrestrial_temperature_trait_tolerance))

# Heterogeneity
terrestrial_temperature_i2_trait <- i2(sd_terrestrial_temperature_trait, Terrestrial_Temperature_Trait_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

terrestrial_temperature_trait_means_list <- c(mean_abs_b_terrestrial_temperature_trait_behavioural, mean_abs_b_terrestrial_temperature_trait_biochem,
                                              mean_abs_b_terrestrial_temperature_trait_morphology, mean_abs_b_terrestrial_temperature_trait_tolerance)
terrestrial_temperature_trait_low_ci <- c(ci.abs_terrestrial_temperature_behavioural[1], ci.abs_terrestrial_temperature_biochem[1], 
                                          ci.abs_terrestrial_temperature_morphology[1], ci.abs_terrestrial_temperature_tolerance[1])
terrestrial_temperature_trait_high_ci <- c(ci.abs_terrestrial_temperature_behavioural[2], ci.abs_terrestrial_temperature_biochem[2], 
                                           ci.abs_terrestrial_temperature_morphology[2], ci.abs_terrestrial_temperature_tolerance[2])
terrestrial_temperature_trait_categories <- c("Behavioural", "Biochemical Assay", 
                                              "Morphology", "Tolerance")

terrestrial_temperature_trait_summary <- matrix(c(terrestrial_temperature_trait_means_list, terrestrial_temperature_trait_low_ci, 
                                                  terrestrial_temperature_trait_high_ci), 
                                         nrow = 4, ncol = 3, byrow = FALSE, 
                                         dimnames = list(c(terrestrial_temperature_trait_categories), 
                                                         c("Mean", "Low_CI", "High_CI")))
terrestrial_temperature_trait_summary <- data.frame(terrestrial_temperature_trait_summary)

# Preparing Graph - Combined

Terrestrial_Temperature_Trait_rnames <- c("Behavioural", "Biochemical Assay",
                                          "Morphological", "Tolerance")

Terrestrial_Temperature_Trait_k <- data.frame("k" = c(Terrestrial_Temperature_Trait_Exploration["Behavioural", "Freq"], 
                                                      Terrestrial_Temperature_Trait_Exploration["Biochemical Assay", "Freq"], 
                                                      Terrestrial_Temperature_Trait_Exploration["Morphology", "Freq"], 
                                                      Terrestrial_Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                              row.names = Terrestrial_Temperature_Trait_rnames)

Terrestrial_Temperature_Trait_group_no <- data.frame("Spp No." = c(Terrestrial_Temperature_Trait_Species_Count["Behavioural", "Freq"], 
                                                                   Terrestrial_Temperature_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                                   Terrestrial_Temperature_Trait_Species_Count["Morphology", "Freq"], 
                                                                   Terrestrial_Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                                     row.names = Terrestrial_Temperature_Trait_rnames)

Terrestrial_Temperature_Trait_study <- data.frame("Study" = c(Terrestrial_Temperature_Trait_Study_Count["Behavioural", "Freq"], 
                                                              Terrestrial_Temperature_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                              Terrestrial_Temperature_Trait_Study_Count["Morphology", "Freq"], 
                                                              Terrestrial_Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                                  row.names = Terrestrial_Temperature_Trait_rnames)

terrestrial_temperature_trait_summary_Reorder <- terrestrial_temperature_trait_summary[c("Behavioural", "Biochemical Assay",
                                                                                         "Morphology", "Tolerance"), ]

Terrestrial_Temperature_Trait_table <- data.frame(estimate = terrestrial_temperature_trait_summary_Reorder[,"Mean"], 
                                                  lowerCL = terrestrial_temperature_trait_summary_Reorder[,"Low_CI"], 
                                                  upperCL = terrestrial_temperature_trait_summary_Reorder[,"High_CI"], 
                                                  K = Terrestrial_Temperature_Trait_k[,1], 
                                                  group_no = Terrestrial_Temperature_Trait_group_no[,1], 
                                                  row.names = Terrestrial_Temperature_Trait_rnames)
Terrestrial_Temperature_Trait_table$name <- row.names(Terrestrial_Temperature_Trait_table)

Terrestrial_Temperature_Trait_raw_mean <- c(b_abs_terrestrial_temperature_trait_behavioural, b_abs_terrestrial_temperature_trait_biochem, 
                                            b_abs_terrestrial_temperature_trait_morphology, b_abs_terrestrial_temperature_trait_tolerance)

Terrestrial_Temperature_Trait_raw_name <- c(replicate(8000, "Behavioural"), 
                                            replicate(8000, "Biochemical Assay"), 
                                            replicate(8000, "Morphological"), 
                                            replicate(8000, "Tolerance"))

Terrestrial_Temperature_Trait_raw_df <- data.frame("Model" = Terrestrial_Temperature_Trait_raw_name, 
                                                   "Effect" = Terrestrial_Temperature_Trait_raw_mean)

# Graph code - Combined

Terrestrial_Temperature_Trait_Order <- c("Tolerance", "Morphological",
                                         "Biochemical Assay", "Behavioural")

density_terrestrial_temperature_trait <- Terrestrial_Temperature_Trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Trait_Order)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = Terrestrial_Temperature_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Trait_Order)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)), xmin = max(Terrestrial_Temperature_Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                        size = 1) +
                                         geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)), xmin = min(Terrestrial_Temperature_Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                          vjust = c(-2.7, -2.7, -0.8, -2.7))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                         scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                         coord_cartesian(xlim = c(-0.01, 1.25)) +
                                         annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(Terrestrial_Temperature_Trait_table["Tolerance", "K"], 
                                                                       Terrestrial_Temperature_Trait_table["Morphological", "K"],
                                                                       Terrestrial_Temperature_Trait_table["Biochemical Assay", "K"], 
                                                                       Terrestrial_Temperature_Trait_table["Behavioural", "K"]), "~","(", 
                                                                     c(Terrestrial_Temperature_Trait_table["Tolerance", "group_no"], 
                                                                       Terrestrial_Temperature_Trait_table["Morphological", "group_no"],  
                                                                       Terrestrial_Temperature_Trait_table["Biochemical Assay", "group_no"], 
                                                                       Terrestrial_Temperature_Trait_table["Behavioural", "group_no"]), 
                                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_morphology)-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                                    x = rev(Terrestrial_Temperature_Trait_table$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_trait #(400x400)

# Preparing Graph - Part 1

Terrestrial_Temperature_Trait_rnames_1 <- c("Behavioural", "Biochemical Assay")

Terrestrial_Temperature_Trait_k_1 <- data.frame("k" = c(Terrestrial_Temperature_Trait_Exploration["Behavioural", "Freq"], 
                                                        Terrestrial_Temperature_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                                row.names = Terrestrial_Temperature_Trait_rnames_1)

Terrestrial_Temperature_Trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Temperature_Trait_Species_Count["Behavioural", "Freq"], 
                                                                     Terrestrial_Temperature_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                                      row.names = Terrestrial_Temperature_Trait_rnames_1)

Terrestrial_Temperature_Trait_study_1 <- data.frame("Study" = c(Terrestrial_Temperature_Trait_Study_Count["Behavioural", "Freq"], 
                                                                Terrestrial_Temperature_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                                    row.names = Terrestrial_Temperature_Trait_rnames_1)

terrestrial_temperature_trait_summary_Reorder_1 <- terrestrial_temperature_trait_summary[c("Behavioural", "Biochemical Assay"), ]

Terrestrial_Temperature_Trait_table_1 <- data.frame(estimate = terrestrial_temperature_trait_summary_Reorder_1[,"Mean"], 
                                                    lowerCL = terrestrial_temperature_trait_summary_Reorder_1[,"Low_CI"], 
                                                    upperCL = terrestrial_temperature_trait_summary_Reorder_1[,"High_CI"], 
                                                    K = Terrestrial_Temperature_Trait_k_1[,1], 
                                                    group_no = Terrestrial_Temperature_Trait_group_no_1[,1], 
                                                    row.names = Terrestrial_Temperature_Trait_rnames_1)
Terrestrial_Temperature_Trait_table_1$name <- row.names(Terrestrial_Temperature_Trait_table_1)

Terrestrial_Temperature_Trait_raw_mean_1 <- c(b_abs_terrestrial_temperature_trait_behavioural, b_abs_terrestrial_temperature_trait_biochem)

Terrestrial_Temperature_Trait_raw_name_1 <- c(replicate(8000, "Behavioural"), 
                                              replicate(8000, "Biochemical Assay"))

Terrestrial_Temperature_Trait_raw_df_1 <- data.frame("Model" = Terrestrial_Temperature_Trait_raw_name_1, 
                                                     "Effect" = Terrestrial_Temperature_Trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Temperature_Trait_Order_1 <- c("Biochemical Assay", "Behavioural")

density_terrestrial_temperature_trait_1 <- Terrestrial_Temperature_Trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Trait_Order_1)) %>%
                                           ggplot() +
                                           geom_density_ridges(data = Terrestrial_Temperature_Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Trait_Order_1)), 
                                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                          size = 1) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)), xmin = max(Terrestrial_Temperature_Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                          size = 1) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)), xmin = min(Terrestrial_Temperature_Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                          size = 1) +
                                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                           size = 1, fatten = 2) +
                                           theme_bw() +
                                           guides(fill = "none", colour = "none") +
                                           labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                            vjust = c(-0.8, -2.7))) +
                                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                           theme(axis.ticks = element_blank()) +
                                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                           scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                           scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                           coord_cartesian(xlim = c(-0.01, 1.25)) +
                                           annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)+0.4),
                                           label= paste("italic(k)==", c(Terrestrial_Temperature_Trait_table_1["Biochemical Assay", "K"], 
                                                                         Terrestrial_Temperature_Trait_table_1["Behavioural", "K"]), "~","(", 
                                                                       c(Terrestrial_Temperature_Trait_table_1["Biochemical Assay", "group_no"], 
                                                                         Terrestrial_Temperature_Trait_table_1["Behavioural", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_behavioural)-1)*100, 2), nsmall = 2), "%")), 
                                                         x = rev(Terrestrial_Temperature_Trait_table_1$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_trait_1 #(400x240)

# Preparing Graph - Part 2

Terrestrial_Temperature_Trait_rnames_2 <- c("Morphological", "Tolerance")

Terrestrial_Temperature_Trait_k_2 <- data.frame("k" = c(Terrestrial_Temperature_Trait_Exploration["Morphology", "Freq"], 
                                                        Terrestrial_Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                                row.names = Terrestrial_Temperature_Trait_rnames_2)

Terrestrial_Temperature_Trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Temperature_Trait_Species_Count["Morphology", "Freq"], 
                                                                     Terrestrial_Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                                       row.names = Terrestrial_Temperature_Trait_rnames_2)

Terrestrial_Temperature_Trait_study_2 <- data.frame("Study" = c(Terrestrial_Temperature_Trait_Study_Count["Morphology", "Freq"], 
                                                                Terrestrial_Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                                    row.names = Terrestrial_Temperature_Trait_rnames_2)

terrestrial_temperature_trait_summary_Reorder_2 <- terrestrial_temperature_trait_summary[c("Morphology", "Tolerance"), ]

Terrestrial_Temperature_Trait_table_2 <- data.frame(estimate = terrestrial_temperature_trait_summary_Reorder_2[,"Mean"], 
                                                    lowerCL = terrestrial_temperature_trait_summary_Reorder_2[,"Low_CI"], 
                                                    upperCL = terrestrial_temperature_trait_summary_Reorder_2[,"High_CI"], 
                                                    K = Terrestrial_Temperature_Trait_k_2[,1], 
                                                    group_no = Terrestrial_Temperature_Trait_group_no_2[,1], 
                                                    row.names = Terrestrial_Temperature_Trait_rnames_2)
Terrestrial_Temperature_Trait_table_2$name <- row.names(Terrestrial_Temperature_Trait_table_2)

Terrestrial_Temperature_Trait_raw_mean_2 <- c(b_abs_terrestrial_temperature_trait_morphology, b_abs_terrestrial_temperature_trait_tolerance)

Terrestrial_Temperature_Trait_raw_name_2 <- c(replicate(8000, "Morphological"), 
                                              replicate(8000, "Tolerance"))

Terrestrial_Temperature_Trait_raw_df_2 <- data.frame("Model" = Terrestrial_Temperature_Trait_raw_name_2, 
                                                     "Effect" = Terrestrial_Temperature_Trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Temperature_Trait_Order_2 <- c("Tolerance", "Morphological")

density_terrestrial_temperature_trait_2 <- Terrestrial_Temperature_Trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Temperature_Trait_Order_2)) %>%
                                           ggplot() +
                                           geom_density_ridges(data = Terrestrial_Temperature_Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Temperature_Trait_Order_2)), 
                                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                          size = 1) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)), xmin = max(Terrestrial_Temperature_Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                          size = 1) +
                                           geom_linerange(aes(y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)), xmin = min(Terrestrial_Temperature_Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                          size = 1) +
                                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                           size = 1, fatten = 2) +
                                           theme_bw() +
                                           guides(fill = "none", colour = "none") +
                                           labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                            vjust = c(-2.7, -2.7))) +
                                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                           theme(axis.ticks = element_blank()) +
                                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                           scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                           scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                           coord_cartesian(xlim = c(-0.01, 1.25)) +
                                           annotate('text',  x = 1.25, y = (seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)+0.4),
                                           label= paste("italic(k)==", c(Terrestrial_Temperature_Trait_table_2["Tolerance", "K"], 
                                                                         Terrestrial_Temperature_Trait_table_2["Morphological", "K"]), "~","(", 
                                                                       c(Terrestrial_Temperature_Trait_table_2["Tolerance", "group_no"], 
                                                                         Terrestrial_Temperature_Trait_table_2["Morphological", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(b_abs_terrestrial_temperature_trait_morphology)-1)*100, 2), nsmall = 2), "%")), 
                                                         x = rev(Terrestrial_Temperature_Trait_table_2$estimate+0.2), y = (seq(1, dim(Terrestrial_Temperature_Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_temperature_trait_2 #(400x240)

##### Aquatic and Temperature Model #####
Aquatic_Temperature_Data <- data %>% filter(Ecosystem == "Aquatic" & Type == "Temperature")
Aquatic_Temperature_Species <- Aquatic_Temperature_Data %>% select("phylo") %>% unique()

Aquatic_Temperature_A <- as.data.frame(A)
Aquatic_Temperature_A <- Aquatic_Temperature_A[c(Aquatic_Temperature_Species$phylo), c(Aquatic_Temperature_Species$phylo)]
Aquatic_Temperature_A <- as.matrix(Aquatic_Temperature_A)

system.time(
  aquatic_temperature <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                   ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                   data = Aquatic_Temperature_Data,
                                   family = gaussian(),
                                   data2 = list(A = Aquatic_Temperature_A), 
                                   chains = 4, 
                                   cores = 4,
                                   iter = 12000,
                                   warmup = 2000,
                                   thin = 5,
                                   prior = priors,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                                   file = "./aquatic_temperature_model",
                                   file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_aquatic_temperature <- as_draws_df(aquatic_temperature, variable = "b_Intercept")
b_aquatic_temperature <- data.frame(b_aquatic_temperature$b_Intercept)

sd_aquatic_temperature <- as_draws_df(aquatic_temperature, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_temperature <- data.frame("sd_obs__Intercept" = sd_aquatic_temperature$sd_obs__Intercept, 
                                     "sd_phylo__Intercept" = sd_aquatic_temperature$sd_phylo__Intercept, 
                                     "sd_Study_ID__Intercept" = sd_aquatic_temperature$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic_temperature <-  sapply(b_aquatic_temperature, mean)
ci_b_aquatic_temperature <- as.vector(HPDinterval(as.mcmc(b_aquatic_temperature)))
pMCMC_b_aquatic_temperature <- 2*(1 - max(table(b_aquatic_temperature<0) / nrow(b_aquatic_temperature)))

# Absolute magnitude
b_abs_aquatic_temperature <- folded_norm(b_aquatic_temperature[,1], rowSums(sd_aquatic_temperature))
mean_abs_b_aquatic_temperature <- mean(b_abs_aquatic_temperature)
ci.abs_aquatic_temperature <- HPDinterval(as.mcmc(b_abs_aquatic_temperature))

# Heterogeneity
aquatic_temperature_i2 <- i2(sd_aquatic_temperature, Aquatic_Temperature_Data$Variance_Type_Adjusted) 

##### Aquatic and Temperature Model - Plasticity Mechanism Meta-regression #####
Aquatic_Temperature_Plasticity_Exploration <- Aquatic_Temperature_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Temperature_Plasticity_Exploration) <- Aquatic_Temperature_Plasticity_Exploration$Plasticity_Category

Aquatic_Temperature_Plasticity_Data <- Aquatic_Temperature_Data %>% filter(Plasticity_Category == "Acclimation")

Aquatic_Temperature_Plasticity_Species_Count <- Aquatic_Temperature_Plasticity_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                                table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                                select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Temperature_Plasticity_Species_Count) <- Aquatic_Temperature_Plasticity_Species_Count$Plasticity_Category

Aquatic_Temperature_Plasticity_Study_Count <- Aquatic_Temperature_Plasticity_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                              table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                              select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Temperature_Plasticity_Study_Count) <- Aquatic_Temperature_Plasticity_Study_Count$Plasticity_Category

Aquatic_Temperature_Plasticity_Species <- Aquatic_Temperature_Plasticity_Data %>% select("phylo") %>% unique()

Aquatic_Temperature_Plasticity_A <- as.data.frame(A)
Aquatic_Temperature_Plasticity_A <- Aquatic_Temperature_Plasticity_A[c(Aquatic_Temperature_Plasticity_Species$phylo), c(Aquatic_Temperature_Plasticity_Species$phylo)]
Aquatic_Temperature_Plasticity_A <- as.matrix(Aquatic_Temperature_Plasticity_A)


system.time(
  aquatic_temperature_plastic <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                           ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                           data = Aquatic_Temperature_Plasticity_Data,
                                           family = gaussian(),
                                           data2 = list(A = Aquatic_Temperature_Plasticity_A), 
                                           chains = 4, 
                                           cores = 4,
                                           iter = 12000,
                                           warmup = 2000,
                                           thin = 5,
                                           prior = priors,
                                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                                           file = "./aquatic_temperature_plastic_model",
                                           file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_temperature_plastic <- as_draws_df(aquatic_temperature_plastic, variable = "b_Intercept")
b_aquatic_temperature_plastic <- data.frame(b_aquatic_temperature_plastic$b_Intercept)

sd_aquatic_temperature_plastic <- as_draws_df(aquatic_temperature_plastic, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_temperature_plastic <- data.frame("sd_obs__Intercept" = sd_aquatic_temperature_plastic$sd_obs__Intercept, 
                                             "sd_phylo__Intercept" = sd_aquatic_temperature_plastic$sd_phylo__Intercept, 
                                             "sd_Study_ID__Intercept" = sd_aquatic_temperature_plastic$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic_temperature_plastic <-  sapply(b_aquatic_temperature_plastic, mean)
ci_b_aquatic_temperature_plastic <- as.vector(HPDinterval(as.mcmc(b_aquatic_temperature_plastic)))
pMCMC_b_aquatic_temperature_plastic <- 2*(1 - max(table(b_aquatic_temperature_plastic<0) / nrow(b_aquatic_temperature_plastic)))

# Absolute magnitude
b_abs_aquatic_temperature_plastic <- folded_norm(b_aquatic_temperature_plastic[,1], rowSums(sd_aquatic_temperature_plastic))
mean_abs_b_aquatic_temperature_plastic <- mean(b_abs_aquatic_temperature_plastic)
ci.abs_aquatic_temperature_plastic <- HPDinterval(as.mcmc(b_abs_aquatic_temperature_plastic))

# Heterogeneity
aquatic_temperature_plastic_i2 <- i2(sd_aquatic_temperature_plastic, Aquatic_Temperature_Plasticity_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

aquatic_temperature_plasticity_means_list <- c(mean_abs_b_aquatic_temperature_plastic)
aquatic_temperature_plasticity_low_ci <- c(ci.abs_aquatic_temperature_plastic[1])
aquatic_temperature_plasticity_high_ci <- c(ci.abs_aquatic_temperature_plastic[2])
aquatic_temperature_plasticity_categories <- c("Acclimation")

aquatic_temperature_plasticity_summary <- matrix(c(aquatic_temperature_plasticity_means_list, aquatic_temperature_plasticity_low_ci, 
                                                   aquatic_temperature_plasticity_high_ci), 
                                                 nrow = 1, ncol = 3, byrow = FALSE, 
                                                 dimnames = list(c(aquatic_temperature_plasticity_categories), 
                                                                 c("Mean", "Low_CI", "High_CI")))
aquatic_temperature_plasticity_summary <- data.frame(aquatic_temperature_plasticity_summary)

# Preparing Graph

Aquatic_Temperature_Plasticity_rnames <- c("Acclimation")

Aquatic_Temperature_Plasticity_k <- data.frame("k" = c(Aquatic_Temperature_Plasticity_Exploration["Acclimation", "Freq"]), 
                                               row.names = Aquatic_Temperature_Plasticity_rnames)

Aquatic_Temperature_Plasticity_group_no <- data.frame("Spp No." = c(Aquatic_Temperature_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                      row.names = Aquatic_Temperature_Plasticity_rnames)

Aquatic_Temperature_Plasticity_study <- data.frame("Study" = c(Aquatic_Temperature_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                                   row.names = Aquatic_Temperature_Plasticity_rnames)

Aquatic_Temperature_Plasticity_table <- data.frame(estimate = aquatic_temperature_plasticity_summary[,"Mean"], 
                                                   lowerCL = aquatic_temperature_plasticity_summary[,"Low_CI"], 
                                                   upperCL = aquatic_temperature_plasticity_summary[,"High_CI"], 
                                                   K = Aquatic_Temperature_Plasticity_k[,1], 
                                                   group_no = Aquatic_Temperature_Plasticity_group_no[,1], 
                                                   row.names = Aquatic_Temperature_Plasticity_rnames)
Aquatic_Temperature_Plasticity_table$name <- row.names(Aquatic_Temperature_Plasticity_table)

Aquatic_Temperature_Plasticity_raw_mean <- c(b_abs_aquatic_temperature_plastic)

Aquatic_Temperature_Plasticity_raw_name <- c(replicate(8000, "Acclimation"))

Aquatic_Temperature_Plasticity_raw_df <- data.frame("Model" = Aquatic_Temperature_Plasticity_raw_name, 
                                                    "Effect" = Aquatic_Temperature_Plasticity_raw_mean)

# Graph code

Aquatic_Temperature_Plasticity_Order <- c("Acclimation")

density_aquatic_temperature_plasticity <- Aquatic_Temperature_Plasticity_table %>% mutate(name = fct_relevel(name, Aquatic_Temperature_Plasticity_Order)) %>%
                                          ggplot() +
                                          geom_density_ridges(data = Aquatic_Temperature_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Temperature_Plasticity_Order)), 
                                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.06, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                          geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)), xmin = max(Aquatic_Temperature_Plasticity_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)), xmin = min(Aquatic_Temperature_Plasticity_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                         size = 1) +
                                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                          size = 1, fatten = 2) +
                                          theme_bw() +
                                          guides(fill = "none", colour = "none") +
                                          labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                          theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                           vjust = c(-2.7))) +
                                          theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                          theme(axis.ticks = element_blank()) +
                                          theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                          scale_colour_manual(values = c("#2B4E7A")) +
                                          scale_fill_manual(values = c("#2B4E7A")) +
                                          coord_cartesian(xlim = c(-0.01, 1.25)) +
                                          annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)+0.4),
                                                   label= paste("italic(k)==", c(Aquatic_Temperature_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                                               c(Aquatic_Temperature_Plasticity_table["Acclimation", "group_no"]), 
                                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                          geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_temperature_plastic)-1)*100, 2), nsmall = 2), "%")), 
                                                         x = rev(Aquatic_Temperature_Plasticity_table$estimate+0.2), y = (seq(1, dim(Aquatic_Temperature_Plasticity_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_temperature_plasticity #(400x160)

##### Aquatic and Temperature Model - Trait Category Meta-regression #####
Aquatic_Temperature_Trait_Exploration <- Aquatic_Temperature_Data %>% select("Category") %>% table() %>% data.frame()
Aquatic_Temperature_Trait_Exploration <- Aquatic_Temperature_Trait_Exploration %>% filter(Freq > 10)
rownames(Aquatic_Temperature_Trait_Exploration) <- Aquatic_Temperature_Trait_Exploration$Category

Aquatic_Temperature_Trait_Data <- Aquatic_Temperature_Data %>% filter(Category == "Biochemical Assay"|
                                                                      Category == "Tolerance")

Aquatic_Temperature_Trait_Species_Count <- Aquatic_Temperature_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                                           table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                           select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Temperature_Trait_Species_Count) <- Aquatic_Temperature_Trait_Species_Count$Category

Aquatic_Temperature_Trait_Study_Count <- Aquatic_Temperature_Trait_Data %>% select("Study_ID", "Category") %>% 
                                         table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                         select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Temperature_Trait_Study_Count) <- Aquatic_Temperature_Trait_Study_Count$Category

Aquatic_Temperature_Trait_Species <- Aquatic_Temperature_Trait_Data %>% select("phylo") %>% unique()

Aquatic_Temperature_Trait_A <- as.data.frame(A)
Aquatic_Temperature_Trait_A <- Aquatic_Temperature_Trait_A[c(Aquatic_Temperature_Trait_Species$phylo), c(Aquatic_Temperature_Trait_Species$phylo)]
Aquatic_Temperature_Trait_A <- as.matrix(Aquatic_Temperature_Trait_A)

system.time(
  aquatic_temperature_trait <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                         ~ Category + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|gr(obs, by = Category, cor = FALSE)),
                                         data = Aquatic_Temperature_Trait_Data,
                                         family = gaussian(),
                                         data2 = list(A = Aquatic_Temperature_Trait_A), 
                                         chains = 4, 
                                         cores = 4,
                                         iter = 12000,
                                         warmup = 2000,
                                         thin = 5,
                                         prior = priors,
                                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                                         file = "./aquatic_temperature_trait_model",
                                         file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_temperature_trait <- as_draws_df(aquatic_temperature_trait, variable = c("b_Intercept", "b_CategoryTolerance"))
b_aquatic_temperature_trait <- data.frame("b_BiochemicalAssay" = b_aquatic_temperature_trait$b_Intercept, 
                                          "b_Tolerance" = b_aquatic_temperature_trait$b_CategoryTolerance + b_aquatic_temperature_trait$b_Intercept)


sd_aquatic_temperature_trait <- as_draws_df(aquatic_temperature_trait, variable = c("sd_obs__Intercept:CategoryBiochemicalAssay", "sd_obs__Intercept:CategoryTolerance", 
                                                                                    "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_temperature_trait <- data.frame("sd_BiochemicalAssay" = sd_aquatic_temperature_trait$`sd_obs__Intercept:CategoryBiochemicalAssay`, 
                                           "sd_Tolerance" = sd_aquatic_temperature_trait$`sd_obs__Intercept:CategoryTolerance`,
                                           "sd_phylo__Intercept" = sd_aquatic_temperature_trait$`sd_phylo__Intercept`, 
                                           "sd_Study_ID__Intercept" = sd_aquatic_temperature_trait$`sd_Study_ID__Intercept`)

# Overall estimates
# Signed
aquatic_temperature_trait_means <- apply(b_aquatic_temperature_trait, 2, mean)
aquatic_temperature_trait_cis <- apply(b_aquatic_temperature_trait, 2, function(x) HPDinterval(as.mcmc(x)))
aquatic_temperature_trait_pMCMC <- apply(b_aquatic_temperature_trait, 2, function(x) 2*(1 - max(table(x<0) / length(x))))

# Absolute magnitude - Check sd numbers based on what random effects you have added.
b_abs_aquatic_temperature_trait_biochem <- folded_norm(b_aquatic_temperature_trait$b_BiochemicalAssay, sqrt(rowSums(sd_aquatic_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_aquatic_temperature_trait[, "sd_BiochemicalAssay"]^2)))
b_abs_aquatic_temperature_trait_tolerance <- folded_norm(b_aquatic_temperature_trait$b_Tolerance, sqrt(rowSums(sd_aquatic_temperature_trait[,c("sd_phylo__Intercept", "sd_Study_ID__Intercept")]^2 + sd_aquatic_temperature_trait[, "sd_Tolerance"]^2)))
mean_abs_b_aquatic_temperature_trait_biochem <- mean(b_abs_aquatic_temperature_trait_biochem)
mean_abs_b_aquatic_temperature_trait_tolerance <- mean(b_abs_aquatic_temperature_trait_tolerance)
ci.abs_aquatic_temperature_biochem <- HPDinterval(as.mcmc(b_abs_aquatic_temperature_trait_biochem))
ci.abs_aquatic_temperature_tolerance <- HPDinterval(as.mcmc(b_abs_aquatic_temperature_trait_tolerance))

# Heterogeneity
aquatic_temperature_i2_trait <- i2(sd_aquatic_temperature_trait, Aquatic_Temperature_Trait_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

aquatic_temperature_trait_means_list <- c(mean_abs_b_aquatic_temperature_trait_biochem, mean_abs_b_aquatic_temperature_trait_tolerance)
aquatic_temperature_trait_low_ci <- c(ci.abs_aquatic_temperature_biochem[1], ci.abs_aquatic_temperature_tolerance[1])
aquatic_temperature_trait_high_ci <- c(ci.abs_aquatic_temperature_biochem[2], ci.abs_aquatic_temperature_tolerance[2])
aquatic_temperature_trait_categories <- c("Biochemical Assay", "Tolerance")

aquatic_temperature_trait_summary <- matrix(c(aquatic_temperature_trait_means_list, aquatic_temperature_trait_low_ci, 
                                              aquatic_temperature_trait_high_ci), 
                                              nrow = 2, ncol = 3, byrow = FALSE, 
                                              dimnames = list(c(aquatic_temperature_trait_categories), 
                                                              c("Mean", "Low_CI", "High_CI")))
aquatic_temperature_trait_summary <- data.frame(aquatic_temperature_trait_summary)

# Preparing Graph - Combined

Aquatic_Temperature_Trait_rnames <- c("Biochemical Assay", "Tolerance")

Aquatic_Temperature_Trait_k <- data.frame("k" = c(Aquatic_Temperature_Trait_Exploration["Biochemical Assay", "Freq"], 
                                                  Aquatic_Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                          row.names = Aquatic_Temperature_Trait_rnames)

Aquatic_Temperature_Trait_group_no <- data.frame("Spp No." = c(Aquatic_Temperature_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                               Aquatic_Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                                 row.names = Aquatic_Temperature_Trait_rnames)

Aquatic_Temperature_Trait_study <- data.frame("Study" = c(Aquatic_Temperature_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                          Aquatic_Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                              row.names = Aquatic_Temperature_Trait_rnames)

Aquatic_Temperature_Trait_table <- data.frame(estimate = aquatic_temperature_trait_summary[,"Mean"], 
                                              lowerCL = aquatic_temperature_trait_summary[,"Low_CI"], 
                                              upperCL = aquatic_temperature_trait_summary[,"High_CI"], 
                                              K = Aquatic_Temperature_Trait_k[,1], 
                                              group_no = Aquatic_Temperature_Trait_group_no[,1], 
                                              row.names = Aquatic_Temperature_Trait_rnames)
Aquatic_Temperature_Trait_table$name <- row.names(Aquatic_Temperature_Trait_table)

Aquatic_Temperature_Trait_raw_mean <- c(b_abs_aquatic_temperature_trait_biochem, b_abs_aquatic_temperature_trait_tolerance)

Aquatic_Temperature_Trait_raw_name <- c(replicate(8000, "Biochemical Assay"), 
                                        replicate(8000, "Tolerance"))

Aquatic_Temperature_Trait_raw_df <- data.frame("Model" = Aquatic_Temperature_Trait_raw_name, 
                                               "Effect" = Aquatic_Temperature_Trait_raw_mean)

# Graph code - Combined

Aquatic_Temperature_Trait_Order <- c("Tolerance", "Biochemical Assay")

density_aquatic_temperature_trait <- Aquatic_Temperature_Trait_table %>% mutate(name = fct_relevel(name, Aquatic_Temperature_Trait_Order)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = Aquatic_Temperature_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Temperature_Trait_Order)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)), xmin = max(Aquatic_Temperature_Trait_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)), xmin = min(Aquatic_Temperature_Trait_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                      vjust = c(-2.7, -0.8))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                     scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                     coord_cartesian(xlim = c(-0.01, 1.25)) +
                                     annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(Aquatic_Temperature_Trait_table["Tolerance", "K"], 
                                                                   Aquatic_Temperature_Trait_table["Biochemical Assay", "K"]), "~","(", 
                                                                 c(Aquatic_Temperature_Trait_table["Tolerance", "group_no"], 
                                                                   Aquatic_Temperature_Trait_table["Biochemical Assay", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(b_abs_aquatic_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%")), 
                                                x = rev(Aquatic_Temperature_Trait_table$estimate+0.2), y = (seq(1, dim(Aquatic_Temperature_Trait_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_temperature_trait #(400x240)

# Preparing Graph - Part 1

Aquatic_Temperature_Trait_rnames_1 <- c("Biochemical Assay")

Aquatic_Temperature_Trait_k_1 <- data.frame("k" = c(Aquatic_Temperature_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                            row.names = Aquatic_Temperature_Trait_rnames_1)

Aquatic_Temperature_Trait_group_no_1 <- data.frame("Spp No." = c(Aquatic_Temperature_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                                   row.names = Aquatic_Temperature_Trait_rnames_1)

Aquatic_Temperature_Trait_study_1 <- data.frame("Study" = c(Aquatic_Temperature_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                                row.names = Aquatic_Temperature_Trait_rnames_1)

aquatic_temperature_trait_summary_Reorder_1 <- aquatic_temperature_trait_summary[c("Biochemical Assay"), ]

Aquatic_Temperature_Trait_table_1 <- data.frame(estimate = aquatic_temperature_trait_summary_Reorder_1[,"Mean"], 
                                                lowerCL = aquatic_temperature_trait_summary_Reorder_1[,"Low_CI"], 
                                                upperCL = aquatic_temperature_trait_summary_Reorder_1[,"High_CI"], 
                                                K = Aquatic_Temperature_Trait_k_1[,1], 
                                                group_no = Aquatic_Temperature_Trait_group_no_1[,1], 
                                                row.names = Aquatic_Temperature_Trait_rnames_1)
Aquatic_Temperature_Trait_table_1$name <- row.names(Aquatic_Temperature_Trait_table_1)

Aquatic_Temperature_Trait_raw_mean_1 <- c(b_abs_aquatic_temperature_trait_biochem)

Aquatic_Temperature_Trait_raw_name_1 <- c(replicate(8000, "Biochemical Assay"))

Aquatic_Temperature_Trait_raw_df_1 <- data.frame("Model" = Aquatic_Temperature_Trait_raw_name_1, 
                                                 "Effect" = Aquatic_Temperature_Trait_raw_mean_1)

# Graph code - Part 1

Aquatic_Temperature_Trait_Order_1 <- c("Biochemical Assay")

density_aquatic_temperature_trait_1 <- Aquatic_Temperature_Trait_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Temperature_Trait_Order_1)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = Aquatic_Temperature_Trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Temperature_Trait_Order_1)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.07, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)), xmin = max(Aquatic_Temperature_Trait_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)), xmin = min(Aquatic_Temperature_Trait_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-0.8))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#2B4E7A")) +
                                       scale_fill_manual(values = c("#2B4E7A")) +
                                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                                       annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(Aquatic_Temperature_Trait_table_1["Biochemical Assay", "K"]), "~","(", 
                                                                   c(Aquatic_Temperature_Trait_table_1["Biochemical Assay", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_temperature_trait_biochem)-1)*100, 2), nsmall = 2), "%")), 
                                                     x = rev(Aquatic_Temperature_Trait_table_1$estimate+0.2), y = (seq(1, dim(Aquatic_Temperature_Trait_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_temperature_trait_1 #(400x160)

# Preparing Graph - Part 2

Aquatic_Temperature_Trait_rnames_2 <- c("Tolerance")

Aquatic_Temperature_Trait_k_2 <- data.frame("k" = c(Aquatic_Temperature_Trait_Exploration["Tolerance", "Freq"]), 
                                            row.names = Aquatic_Temperature_Trait_rnames_2)

Aquatic_Temperature_Trait_group_no_2 <- data.frame("Spp No." = c(Aquatic_Temperature_Trait_Species_Count["Tolerance", "Freq"]), 
                                                   row.names = Aquatic_Temperature_Trait_rnames_2)

Aquatic_Temperature_Trait_study_2 <- data.frame("Study" = c(Aquatic_Temperature_Trait_Study_Count["Tolerance", "Freq"]), 
                                                row.names = Aquatic_Temperature_Trait_rnames_2)

aquatic_temperature_trait_summary_Reorder_2 <- aquatic_temperature_trait_summary[c("Tolerance"), ]

Aquatic_Temperature_Trait_table_2 <- data.frame(estimate = aquatic_temperature_trait_summary_Reorder_2[,"Mean"], 
                                                lowerCL = aquatic_temperature_trait_summary_Reorder_2[,"Low_CI"], 
                                                upperCL = aquatic_temperature_trait_summary_Reorder_2[,"High_CI"], 
                                                K = Aquatic_Temperature_Trait_k_2[,1], 
                                                group_no = Aquatic_Temperature_Trait_group_no_2[,1], 
                                                row.names = Aquatic_Temperature_Trait_rnames_2)
Aquatic_Temperature_Trait_table_2$name <- row.names(Aquatic_Temperature_Trait_table_2)

Aquatic_Temperature_Trait_raw_mean_2 <- c(b_abs_aquatic_temperature_trait_tolerance)

Aquatic_Temperature_Trait_raw_name_2 <- c(replicate(8000, "Tolerance"))

Aquatic_Temperature_Trait_raw_df_2 <- data.frame("Model" = Aquatic_Temperature_Trait_raw_name_2, 
                                                 "Effect" = Aquatic_Temperature_Trait_raw_mean_2)

# Graph code - Part 2

Aquatic_Temperature_Trait_Order_2 <- c("Tolerance")

density_aquatic_temperature_trait_2 <- Aquatic_Temperature_Trait_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Temperature_Trait_Order_2)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = Aquatic_Temperature_Trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Temperature_Trait_Order_2)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.07, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)), xmin = max(Aquatic_Temperature_Trait_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)), xmin = min(Aquatic_Temperature_Trait_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-2.7))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#5D7AA1")) +
                                       scale_fill_manual(values = c("#5D7AA1")) +
                                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                                       annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(Aquatic_Temperature_Trait_table["Tolerance", "K"]), "~","(", 
                                                                   c(Aquatic_Temperature_Trait_table["Tolerance", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_temperature_trait_tolerance)-1)*100, 2), nsmall = 2), "%")), 
                                                     x = rev(Aquatic_Temperature_Trait_table_2$estimate+0.2), y = (seq(1, dim(Aquatic_Temperature_Trait_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_temperature_trait_2 #(400x160)

##### Aquatic and Salinity Model #####
Aquatic_Salinity_Data <- data %>% filter(Ecosystem == "Aquatic" & Type == "Salinity" & Study_ID != 37)
Aquatic_Salinity_Species <- Aquatic_Salinity_Data %>% select("phylo") %>% unique()

Aquatic_Salinity_A <- as.data.frame(A)
Aquatic_Salinity_A <- Aquatic_Salinity_A[c(Aquatic_Salinity_Species$phylo), c(Aquatic_Salinity_Species$phylo)]
Aquatic_Salinity_A <- as.matrix(Aquatic_Salinity_A)

system.time(
  aquatic_salinity <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                data = Aquatic_Salinity_Data,
                                family = gaussian(),
                                data2 = list(A = Aquatic_Salinity_A), 
                                chains = 4, 
                                cores = 4,
                                iter = 12000,
                                warmup = 2000,
                                thin = 5,
                                prior = priors,
                                control = list(adapt_delta = 0.99, max_treedepth = 15),
                                file = "./aquatic_salinity_model",
                                file_refit = "always"))

####-- Bayesian Model/Data Output --#####

# Extracting the posterior distributions
b_aquatic_salinity <- as_draws_df(aquatic_salinity, variable = "b_Intercept")
b_aquatic_salinity <- data.frame(b_aquatic_salinity$b_Intercept)

sd_aquatic_salinity <- as_draws_df(aquatic_salinity, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_salinity <- data.frame("sd_obs__Intercept" = sd_aquatic_salinity$sd_obs__Intercept, 
                                  "sd_phylo__Intercept" = sd_aquatic_salinity$sd_phylo__Intercept, 
                                  "sd_Study_ID__Intercept" = sd_aquatic_salinity$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic_salinity <-  sapply(b_aquatic_salinity, mean)
ci_b_aquatic_salinity <- as.vector(HPDinterval(as.mcmc(b_aquatic_salinity)))
pMCMC_b_aquatic_salinity <- 2*(1 - max(table(b_aquatic_salinity<0) / nrow(b_aquatic_salinity)))

# Absolute magnitude
b_abs_aquatic_salinity <- folded_norm(b_aquatic_salinity[,1], rowSums(sd_aquatic_salinity))
mean_abs_b_aquatic_salinity <- mean(b_abs_aquatic_salinity)
ci.abs_aquatic_salinity <- HPDinterval(as.mcmc(b_abs_aquatic_salinity))

# Heterogeneity
aquatic_salinity_i2 <- i2(sd_aquatic_salinity, Aquatic_Salinity_Data$Variance_Type_Adjusted) 

##### Aquatic and Salinity Model - Plasticity Mechanism Meta-regression #####
Aquatic_Salinity_Plasticity_Exploration <- Aquatic_Salinity_Data %>% select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Salinity_Plasticity_Exploration) <- Aquatic_Salinity_Plasticity_Exploration$Plasticity_Category

Aquatic_Salinity_Plasticity_Data <- Aquatic_Salinity_Data %>% filter(Plasticity_Category == "Acclimation")

Aquatic_Salinity_Plasticity_Species_Count <- Aquatic_Salinity_Plasticity_Data %>% select("Scientific_Name", "Plasticity_Category") %>% 
                                             table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                             select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Salinity_Plasticity_Species_Count) <- Aquatic_Salinity_Plasticity_Species_Count$Plasticity_Category

Aquatic_Salinity_Plasticity_Study_Count <- Aquatic_Salinity_Plasticity_Data %>% select("Study_ID", "Plasticity_Category") %>% 
                                           table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                           select("Plasticity_Category") %>% table() %>% data.frame()
rownames(Aquatic_Salinity_Plasticity_Study_Count) <- Aquatic_Salinity_Plasticity_Study_Count$Plasticity_Category

Aquatic_Salinity_Plasticity_Species <- Aquatic_Salinity_Plasticity_Data %>% select("phylo") %>% unique()

Aquatic_Salinity_Plasticity_A <- as.data.frame(A)
Aquatic_Salinity_Plasticity_A <- Aquatic_Salinity_Plasticity_A[c(Aquatic_Salinity_Plasticity_Species$phylo), c(Aquatic_Salinity_Plasticity_Species$phylo)]
Aquatic_Salinity_Plasticity_A <- as.matrix(Aquatic_Salinity_Plasticity_A)


system.time(
  aquatic_salinity_plastic <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                        ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                        data = Aquatic_Salinity_Plasticity_Data,
                                        family = gaussian(),
                                        data2 = list(A = Aquatic_Salinity_Plasticity_A), 
                                        chains = 4, 
                                        cores = 4,
                                        iter = 12000,
                                        warmup = 2000,
                                        thin = 5,
                                        prior = priors,
                                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                                        file = "./aquatic_salinity_plastic_model",
                                        file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_salinity_plastic <- as_draws_df(aquatic_salinity_plastic, variable = "b_Intercept")
b_aquatic_salinity_plastic <- data.frame(b_aquatic_salinity_plastic$b_Intercept)

sd_aquatic_salinity_plastic <- as_draws_df(aquatic_salinity_plastic, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_salinity_plastic <- data.frame("sd_obs__Intercept" = sd_aquatic_salinity_plastic$sd_obs__Intercept, 
                                          "sd_phylo__Intercept" = sd_aquatic_salinity_plastic$sd_phylo__Intercept, 
                                          "sd_Study_ID__Intercept" = sd_aquatic_salinity_plastic$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic_salinity_plastic <-  sapply(b_aquatic_salinity_plastic, mean)
ci_b_aquatic_salinity_plastic <- as.vector(HPDinterval(as.mcmc(b_aquatic_salinity_plastic)))
pMCMC_b_aquatic_salinity_plastic <- 2*(1 - max(table(b_aquatic_salinity_plastic<0) / nrow(b_aquatic_salinity_plastic)))

# Absolute magnitude
b_abs_aquatic_salinity_plastic <- folded_norm(b_aquatic_salinity_plastic[,1], rowSums(sd_aquatic_salinity_plastic))
mean_abs_b_aquatic_salinity_plastic <- mean(b_abs_aquatic_salinity_plastic)
ci.abs_aquatic_salinity_plastic <- HPDinterval(as.mcmc(b_abs_aquatic_salinity_plastic))

# Heterogeneity
aquatic_salinity_plastic_i2 <- i2(sd_aquatic_salinity_plastic, Aquatic_Salinity_Plasticity_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

aquatic_salinity_plasticity_means_list <- c(mean_abs_b_aquatic_salinity_plastic)
aquatic_salinity_plasticity_low_ci <- c(ci.abs_aquatic_salinity_plastic[1])
aquatic_salinity_plasticity_high_ci <- c(ci.abs_aquatic_salinity_plastic[2])
aquatic_salinity_plasticity_categories <- c("Acclimation")

aquatic_salinity_plasticity_summary <- matrix(c(aquatic_salinity_plasticity_means_list, aquatic_salinity_plasticity_low_ci, 
                                                aquatic_salinity_plasticity_high_ci), 
                                                 nrow = 1, ncol = 3, byrow = FALSE, 
                                                 dimnames = list(c(aquatic_salinity_plasticity_categories), 
                                                                 c("Mean", "Low_CI", "High_CI")))
aquatic_salinity_plasticity_summary <- data.frame(aquatic_salinity_plasticity_summary)

# Preparing Graph
Aquatic_Salinity_Plasticity_rnames <- c("Acclimation")

Aquatic_Salinity_Plasticity_k <- data.frame("k" = c(Aquatic_Salinity_Plasticity_Exploration["Acclimation", "Freq"]), 
                                            row.names = Aquatic_Salinity_Plasticity_rnames)

Aquatic_Salinity_Plasticity_group_no <- data.frame("Spp No." = c(Aquatic_Salinity_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                   row.names = Aquatic_Salinity_Plasticity_rnames)

Aquatic_Salinity_Plasticity_study <- data.frame("Study" = c(Aquatic_Salinity_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                                row.names = Aquatic_Salinity_Plasticity_rnames)

Aquatic_Salinity_Plasticity_table <- data.frame(estimate = aquatic_salinity_plasticity_summary[,"Mean"], 
                                                lowerCL = aquatic_salinity_plasticity_summary[,"Low_CI"], 
                                                upperCL = aquatic_salinity_plasticity_summary[,"High_CI"], 
                                                K = Aquatic_Salinity_Plasticity_k[,1], 
                                                group_no = Aquatic_Salinity_Plasticity_group_no[,1], 
                                                row.names = Aquatic_Salinity_Plasticity_rnames)
Aquatic_Salinity_Plasticity_table$name <- row.names(Aquatic_Salinity_Plasticity_table)

Aquatic_Salinity_Plasticity_raw_mean <- c(b_abs_aquatic_salinity_plastic)

Aquatic_Salinity_Plasticity_raw_name <- c(replicate(8000, "Acclimation"))

Aquatic_Salinity_Plasticity_raw_df <- data.frame("Model" = Aquatic_Salinity_Plasticity_raw_name, 
                                                 "Effect" = Aquatic_Salinity_Plasticity_raw_mean)

# Graph code

Aquatic_Salinity_Plasticity_Order <- c("Acclimation")

density_aquatic_salinity_plasticity <- Aquatic_Salinity_Plasticity_table %>% mutate(name = fct_relevel(name, Aquatic_Salinity_Plasticity_Order)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = Aquatic_Salinity_Plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Salinity_Plasticity_Order)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.06, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Salinity_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Salinity_Plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-2.7))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#2B4E7A")) +
                                       scale_fill_manual(values = c("#2B4E7A")) +
                                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                                       annotate('text',  x = 1.25, y = (seq(1, dim(Aquatic_Salinity_Plasticity_table)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(Aquatic_Salinity_Plasticity_table["Acclimation", "K"]), "~","(", 
                                                                   c(Aquatic_Salinity_Plasticity_table["Acclimation", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_salinity_plastic)-1)*100, 2), nsmall = 2), "%")), 
                                                      x = rev(Aquatic_Salinity_Plasticity_table$estimate+0.2), y = (seq(1, dim(Aquatic_Salinity_Plasticity_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_salinity_plasticity #(400x160)

##### Aquatic and Salinity Model - Trait Category Meta-regression #####
Aquatic_Salinity_Trait_Exploration <- Aquatic_Salinity_Data %>% select("Category") %>% table() %>% data.frame()
Aquatic_Salinity_Trait_Exploration <- Aquatic_Salinity_Trait_Exploration %>% filter(Freq > 10)
rownames(Aquatic_Salinity_Trait_Exploration) <- Aquatic_Salinity_Trait_Exploration$Category

Aquatic_Salinity_Trait_Data <- Aquatic_Salinity_Data %>% filter(Category == "Biochemical Assay")

Aquatic_Salinity_Trait_Species_Count <- Aquatic_Salinity_Trait_Data %>% select("Scientific_Name", "Category") %>% 
                                        table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                        select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Salinity_Trait_Species_Count) <- Aquatic_Salinity_Trait_Species_Count$Category

Aquatic_Salinity_Trait_Study_Count <- Aquatic_Salinity_Trait_Data %>% select("Study_ID", "Category") %>% 
                                      table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                      select("Category") %>% table() %>% data.frame()
rownames(Aquatic_Salinity_Trait_Study_Count) <- Aquatic_Salinity_Trait_Study_Count$Category

Aquatic_Salinity_Trait_Species <- Aquatic_Salinity_Trait_Data %>% select("phylo") %>% unique()

Aquatic_Salinity_Trait_A <- as.data.frame(A)
Aquatic_Salinity_Trait_A <- Aquatic_Salinity_Trait_A[c(Aquatic_Salinity_Trait_Species$phylo), c(Aquatic_Salinity_Trait_Species$phylo)]
Aquatic_Salinity_Trait_A <- as.matrix(Aquatic_Salinity_Trait_A)

system.time(
  aquatic_salinity_trait <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
                                      ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                      data = Aquatic_Salinity_Trait_Data,
                                      family = gaussian(),
                                      data2 = list(A = Aquatic_Salinity_Trait_A), 
                                      chains = 4, 
                                      cores = 4,
                                      iter = 12000,
                                      warmup = 2000,
                                      thin = 5,
                                      prior = priors,
                                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                                      file = "./aquatic_salinity_trait_model",
                                      file_refit = "always"))

####-- Bayesian Model/Data Output --####

# Extracting the posterior distributions
b_aquatic_salinity_trait <- as_draws_df(aquatic_salinity_trait, variable = "b_Intercept")
b_aquatic_salinity_trait <- data.frame(b_aquatic_salinity_trait$b_Intercept)

sd_aquatic_salinity_trait <- as_draws_df(aquatic_salinity_trait, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_aquatic_salinity_trait <- data.frame("sd_obs__Intercept" = sd_aquatic_salinity_trait$sd_obs__Intercept, 
                                        "sd_phylo__Intercept" = sd_aquatic_salinity_trait$sd_phylo__Intercept, 
                                        "sd_Study_ID__Intercept" = sd_aquatic_salinity_trait$sd_Study_ID__Intercept)

# Overall estimates
# Signed
mean_b_aquatic_salinity_trait <-  sapply(b_aquatic_salinity_trait, mean)
ci_b_aquatic_salinity_trait <- as.vector(HPDinterval(as.mcmc(b_aquatic_salinity_trait)))
pMCMC_b_aquatic_salinity_trait <- 2*(1 - max(table(b_aquatic_salinity_trait<0) / nrow(b_aquatic_salinity_trait)))

# Absolute magnitude
b_abs_aquatic_salinity_trait <- folded_norm(b_aquatic_salinity_trait[,1], rowSums(sd_aquatic_salinity_trait))
mean_abs_b_aquatic_salinity_trait <- mean(b_abs_aquatic_salinity_trait)
ci.abs_aquatic_salinity_trait <- HPDinterval(as.mcmc(b_abs_aquatic_salinity_trait))

# Heterogeneity
aquatic_salinity_trait_i2 <- i2(sd_aquatic_salinity_trait, Aquatic_Salinity_Trait_Data$Variance_Type_Adjusted) 

# Overall_trait_summary

aquatic_salinity_trait_means_list <- c(mean_abs_b_aquatic_salinity_trait)
aquatic_salinity_trait_low_ci <- c(ci.abs_aquatic_salinity_trait[1])
aquatic_salinity_trait_high_ci <- c(ci.abs_aquatic_salinity_trait[2])
aquatic_salinity_trait_categories <- c("Biochemical Assay")

aquatic_salinity_trait_summary <- matrix(c(aquatic_salinity_trait_means_list, aquatic_salinity_trait_low_ci, 
                                           aquatic_salinity_trait_high_ci), 
                                           nrow = 1, ncol = 3, byrow = FALSE, 
                                           dimnames = list(c(aquatic_salinity_plasticity_categories), 
                                                           c("Mean", "Low_CI", "High_CI")))
aquatic_salinity_trait_summary <- data.frame(aquatic_salinity_trait_summary)

# Preparing Graph
Aquatic_Salinity_Trait_rnames <- c("Biochemical Assay")

Aquatic_Salinity_Trait_k <- data.frame("k" = c(Aquatic_Salinity_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                               row.names = Aquatic_Salinity_Trait_rnames)

Aquatic_Salinity_Trait_group_no <- data.frame("Spp No." = c(Aquatic_Salinity_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                                            row.names = Aquatic_Salinity_Trait_rnames)

Aquatic_Salinity_Trait_study <- data.frame("Study" = c(Aquatic_Salinity_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                                       row.names = Aquatic_Salinity_Trait_rnames)

Aquatic_Salinity_Trait_table <- data.frame(estimate = aquatic_salinity_trait_summary[,"Mean"], 
                                           lowerCL = aquatic_salinity_trait_summary[,"Low_CI"], 
                                           upperCL = aquatic_salinity_trait_summary[,"High_CI"], 
                                           K = Aquatic_Salinity_Trait_k[,1], 
                                           group_no = Aquatic_Salinity_Trait_group_no[,1], 
                                           row.names = Aquatic_Salinity_Trait_rnames)
Aquatic_Salinity_Trait_table$name <- row.names(Aquatic_Salinity_Trait_table)

Aquatic_Salinity_Trait_raw_mean <- c(b_abs_aquatic_salinity_trait)

Aquatic_Salinity_Trait_raw_name <- c(replicate(8000, "Biochemical Assay"))

Aquatic_Salinity_Trait_raw_df <- data.frame("Model" = Aquatic_Salinity_Trait_raw_name, 
                                            "Effect" = Aquatic_Salinity_Trait_raw_mean)

# Graph code

Aquatic_Salinity_Trait_Order <- c("Biochemical Assay")

density_aquatic_salinity_trait <- Aquatic_Salinity_Trait_table %>% mutate(name = fct_relevel(name, Aquatic_Salinity_Trait_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = Aquatic_Salinity_Trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Salinity_Trait_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.06, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(Aquatic_Salinity_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Aquatic_Salinity_Trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-0.8))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#2B4E7A")) +
                                  scale_fill_manual(values = c("#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.01, 2.25)) +
                                  annotate('text',  x = 2.25, y = (seq(1, dim(Aquatic_Salinity_Trait_table)[1], 1)+0.4),
                                           label= paste("italic(k)==", c(Aquatic_Salinity_Trait_table["Biochemical Assay", "K"]), "~","(", 
                                                                       c(Aquatic_Salinity_Trait_table["Biochemical Assay", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_salinity_trait)-1)*100, 2), nsmall = 2), "%")), 
                                                 x = rev(Aquatic_Salinity_Trait_table$estimate+0.2), y = (seq(1, dim(Aquatic_Salinity_Trait_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_salinity_trait #(400x160)

##### Meta-analytic Models (Intercept Only) - Subset Graph #####

# Preparing Data - Combined

Intercept_rnames <- c("Overall", "Terrestrial", "Aquatic", "Temperature", 
                      "Terrestrial and Temperature", "Aquatic and Temperature", 
                      "Aquatic and Salinity")

Intercept_means_df <- data.frame("Mean" = c(mean_abs_b, 
                                            mean_abs_b_terrestrial, 
                                            mean_abs_b_aquatic,
                                            mean_abs_b_temperature, 
                                            mean_abs_b_terrestrial_temperature, 
                                            mean_abs_b_aquatic_temperature,
                                            mean_abs_b_aquatic_salinity), 
                                 row.names = Intercept_rnames)

Intercept_low_df <- data.frame("Low CI" = c(ci.abs[1], 
                                            ci.abs_terrestrial[1], 
                                            ci.abs_aquatic[1], 
                                            ci.abs_temperature[1], 
                                            ci.abs_terrestrial_temperature[1], 
                                            ci.abs_aquatic_temperature[1], 
                                            ci.abs_aquatic_salinity[1]), 
                               row.names = Intercept_rnames)

Intercept_high_df <- data.frame("High CI" = c(ci.abs[2], 
                                              ci.abs_terrestrial[2], 
                                              ci.abs_aquatic[2], 
                                              ci.abs_temperature[2], 
                                              ci.abs_terrestrial_temperature[2], 
                                              ci.abs_aquatic_temperature[2], 
                                              ci.abs_aquatic_salinity[2]), 
                                row.names = Intercept_rnames)

Intercept_k <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                  length(Terrestrial_Data$Effect_Size_ID), 
                                  length(Aquatic_Data$Effect_Size_ID), 
                                  length(Temperature_Data$Effect_Size_ID), 
                                  length(Terrestrial_Temperature_Data$Effect_Size_ID), 
                                  length(Aquatic_Temperature_Data$Effect_Size_ID), 
                                  length(Aquatic_Salinity_Data$Effect_Size_ID)),
                          row.names = Intercept_rnames)

Intercept_group_no <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                               length(unique(Terrestrial_Data$Scientific_Name)), 
                                               length(unique(Aquatic_Data$Scientific_Name)), 
                                               length(unique(Temperature_Data$Scientific_Name)), 
                                               length(unique(Terrestrial_Temperature_Data$Scientific_Name)), 
                                               length(unique(Aquatic_Temperature_Data$Scientific_Name)), 
                                               length(unique(Aquatic_Salinity_Data$Scientific_Name))), 
                                 row.names = Intercept_rnames)

Intercept_study <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                          length(unique(Terrestrial_Data$Study_ID)), 
                                          length(unique(Aquatic_Data$Study_ID)), 
                                          length(unique(Temperature_Data$Study_ID)), 
                                          length(unique(Terrestrial_Temperature_Data$Study_ID)), 
                                          length(unique(Aquatic_Temperature_Data$Study_ID)), 
                                          length(unique(Aquatic_Salinity_Data$Study_ID))), 
                              row.names = Intercept_rnames)

Intercept_table <- data.frame(estimate = Intercept_means_df[,1], 
                              lowerCL = Intercept_low_df[,1], 
                              upperCL = Intercept_high_df[,1], 
                              K = Intercept_k[,1], 
                              group_no = Intercept_group_no[,1], 
                              row.names = Intercept_rnames)
Intercept_table$name <- row.names(Intercept_table)

Intercept_raw_mean <- c(b_abs, 
                        b_abs_terrestrial, 
                        b_abs_aquatic, 
                        b_abs_temperature, 
                        b_abs_terrestrial_temperature, 
                        b_abs_aquatic_temperature, 
                        b_abs_aquatic_salinity)

Intercept_raw_name <- c(replicate(8000, "Overall"), 
                        replicate(8000, "Terrestrial"), 
                        replicate(8000, "Aquatic"), 
                        replicate(8000, "Temperature"), 
                        replicate(8000, "Terrestrial and Temperature"), 
                        replicate(8000, "Aquatic and Temperature"), 
                        replicate(8000, "Aquatic and Salinity"))

Intercept_raw_df <- data.frame("Model" = Intercept_raw_name, 
                               "Effect" = Intercept_raw_mean)

# Graph Code - Combined

Intercept_Order <- c("Aquatic and Salinity", "Aquatic and Temperature", "Terrestrial and Temperature", 
                     "Temperature", "Aquatic", "Terrestrial", "Overall")

density_intercept <- Intercept_table %>% mutate(name = fct_relevel(name, Intercept_Order)) %>%
                     ggplot() +
                     geom_density_ridges(data = Intercept_raw_df %>% mutate(Model = fct_relevel(Model, Intercept_Order)), 
                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                     geom_linerange(aes(y = rev(seq(1, dim(Intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(Intercept_table)[1], 1)), xmin = max(Intercept_raw_df$Effect)+0.001, xmax = 1.5, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(Intercept_table)[1], 1)), xmin = min(Intercept_raw_df$Effect)-0.001, xmax = -0.2, colour = name),
                                    size = 1) +
                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                     theme_bw() +
                     guides(fill = "none", colour = "none") +
                     labs(x = TeX("Effect Size (PRRD)"), y = "") +
                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                      vjust = c(-0.4, -0.4, -0.4, -2.7, -2.7, -2.7, -2.7))) +
                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                     theme(axis.ticks = element_blank()) +
                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 7)) +
                     scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                     scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                     coord_cartesian(xlim = c(-0.01, 1.25)) +
                     annotate('text',  x = 1.25, y = (seq(1, dim(Intercept_table)[1], 1)+0.4),
                     label= paste("italic(k)==", c(Intercept_table["Aquatic and Salinity", "K"], 
                                                   Intercept_table["Aquatic and Temperature", "K"], 
                                                   Intercept_table["Terrestrial and Temperature", "K"], 
                                                   Intercept_table["Temperature", "K"], 
                                                   Intercept_table["Aquatic", "K"], 
                                                   Intercept_table["Terrestrial", "K"], 
                                                   Intercept_table["Overall", "K"]), "~","(", 
                                                 c(Intercept_table["Aquatic and Salinity", "group_no"], 
                                                   Intercept_table["Aquatic and Temperature", "group_no"], 
                                                   Intercept_table["Terrestrial and Temperature", "group_no"], 
                                                   Intercept_table["Temperature", "group_no"], 
                                                   Intercept_table["Aquatic", "group_no"], 
                                                   Intercept_table["Terrestrial", "group_no"], 
                                                   Intercept_table["Overall", "group_no"]), 
                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                     geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_salinity)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_aquatic_temperature)-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(b_abs_terrestrial_temperature)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_temperature)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_aquatic)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs_terrestrial)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(b_abs)-1)*100, 2), nsmall = 2), "%")), 
                                x = rev(Intercept_table$estimate+0.2), y = (seq(1, dim(Intercept_table)[1], 1)+0.4)), size = 3.5)

density_intercept #(400x640)

# Preparing Data - Part 1

Intercept_rnames_1 <- c("Overall", "Terrestrial", "Aquatic")

Intercept_means_df_1 <- data.frame("Mean" = c(mean_abs_b, 
                                              mean_abs_b_terrestrial, 
                                              mean_abs_b_aquatic), 
                                   row.names = Intercept_rnames_1)

Intercept_low_df_1 <- data.frame("Low CI" = c(ci.abs[1], 
                                              ci.abs_terrestrial[1], 
                                              ci.abs_aquatic[1]), 
                                 row.names = Intercept_rnames_1)

Intercept_high_df_1 <- data.frame("High CI" = c(ci.abs[2], 
                                                ci.abs_terrestrial[2], 
                                                ci.abs_aquatic[2]), 
                                  row.names = Intercept_rnames_1)

Intercept_k_1 <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                    length(Terrestrial_Data$Effect_Size_ID), 
                                    length(Aquatic_Data$Effect_Size_ID)),
                            row.names = Intercept_rnames_1)

Intercept_group_no_1 <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                                 length(unique(Terrestrial_Data$Scientific_Name)), 
                                                 length(unique(Aquatic_Data$Scientific_Name))), 
                                   row.names = Intercept_rnames_1)

Intercept_study_1 <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                            length(unique(Terrestrial_Data$Study_ID)), 
                                            length(unique(Aquatic_Data$Study_ID))), 
                                row.names = Intercept_rnames_1)

Intercept_table_1 <- data.frame(estimate = Intercept_means_df_1[,1], 
                                lowerCL = Intercept_low_df_1[,1], 
                                upperCL = Intercept_high_df_1[,1], 
                                K = Intercept_k_1[,1], 
                                group_no = Intercept_group_no_1[,1], 
                                row.names = Intercept_rnames_1)
Intercept_table_1$name <- row.names(Intercept_table_1)

Intercept_raw_mean_1 <- c(b_abs, 
                          b_abs_terrestrial, 
                          b_abs_aquatic)

Intercept_raw_name_1 <- c(replicate(8000, "Overall"), 
                          replicate(8000, "Terrestrial"), 
                          replicate(8000, "Aquatic"))

Intercept_raw_df_1 <- data.frame("Model" = Intercept_raw_name_1, 
                                 "Effect" = Intercept_raw_mean_1)

# Graph Code - Part 1

Intercept_Order_1 <- c("Aquatic", "Terrestrial", "Overall")

density_intercept_1 <- Intercept_table_1 %>% mutate(name = fct_relevel(name, Intercept_Order_1)) %>%
                       ggplot() +
                       geom_density_ridges(data = Intercept_raw_df_1 %>% mutate(Model = fct_relevel(Model, Intercept_Order_1)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_1)[1], 1)), xmin = max(Intercept_raw_df_1$Effect)+0.001, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_1)[1], 1)), xmin = min(Intercept_raw_df_1$Effect)-0.001, xmax = -0.2, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 7)) +
                       scale_colour_manual(values = c("#2B4E7A", "#1B3D6B", "#0D2A51")) +
                       scale_fill_manual(values = c("#2B4E7A", "#1B3D6B", "#0D2A51")) +
                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                       annotate('text',  x = 1.25, y = (seq(1, dim(Intercept_table_1)[1], 1)+0.4),
                       label= paste("italic(k)==", c(Intercept_table_1["Aquatic", "K"], 
                                                     Intercept_table_1["Terrestrial", "K"], 
                                                     Intercept_table_1["Overall", "K"]), "~","(", 
                                                   c(Intercept_table_1["Aquatic", "group_no"], 
                                                     Intercept_table_1["Terrestrial", "group_no"], 
                                                     Intercept_table_1["Overall", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                        geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic)-1)*100, 2), nsmall = 2), "%"), 
                                               paste(format(round(mean(exp(b_abs_terrestrial)-1)*100, 2), nsmall = 2), "%"), 
                                               paste(format(round(mean(exp(b_abs)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Intercept_table_1$estimate+0.2), y = (seq(1, dim(Intercept_table_1)[1], 1)+0.4)), size = 3.5)

density_intercept_1 #(400x320)

# Preparing Data - Part 2

Intercept_rnames_2 <- c("Temperature", "Terrestrial and Temperature", "Aquatic and Temperature", 
                        "Aquatic and Salinity")

Intercept_means_df_2 <- data.frame("Mean" = c(mean_abs_b_temperature, 
                                              mean_abs_b_terrestrial_temperature, 
                                              mean_abs_b_aquatic_temperature,
                                              mean_abs_b_aquatic_salinity), 
                                   row.names = Intercept_rnames_2)

Intercept_low_df_2 <- data.frame("Low CI" = c(ci.abs_temperature[1], 
                                              ci.abs_terrestrial_temperature[1], 
                                              ci.abs_aquatic_temperature[1], 
                                              ci.abs_aquatic_salinity[1]), 
                                 row.names = Intercept_rnames_2)

Intercept_high_df_2 <- data.frame("High CI" = c(ci.abs_temperature[2], 
                                                ci.abs_terrestrial_temperature[2], 
                                                ci.abs_aquatic_temperature[2], 
                                                ci.abs_aquatic_salinity[2]), 
                                  row.names = Intercept_rnames_2)

Intercept_k_2 <- data.frame("k" = c(length(Temperature_Data$Effect_Size_ID), 
                                    length(Terrestrial_Temperature_Data$Effect_Size_ID), 
                                    length(Aquatic_Temperature_Data$Effect_Size_ID), 
                                    length(Aquatic_Salinity_Data$Effect_Size_ID)),
                            row.names = Intercept_rnames_2)

Intercept_group_no_2 <- data.frame("Spp No." = c(length(unique(Temperature_Data$Scientific_Name)), 
                                                 length(unique(Terrestrial_Temperature_Data$Scientific_Name)), 
                                                 length(unique(Aquatic_Temperature_Data$Scientific_Name)), 
                                                 length(unique(Aquatic_Salinity_Data$Scientific_Name))), 
                                   row.names = Intercept_rnames_2)

Intercept_study_2 <- data.frame("Study" = c(length(unique(Temperature_Data$Study_ID)), 
                                            length(unique(Terrestrial_Temperature_Data$Study_ID)), 
                                            length(unique(Aquatic_Temperature_Data$Study_ID)), 
                                            length(unique(Aquatic_Salinity_Data$Study_ID))), 
                                row.names = Intercept_rnames_2)

Intercept_table_2 <- data.frame(estimate = Intercept_means_df_2[,1], 
                                lowerCL = Intercept_low_df_2[,1], 
                                upperCL = Intercept_high_df_2[,1], 
                                K = Intercept_k_2[,1], 
                                group_no = Intercept_group_no_2[,1], 
                                row.names = Intercept_rnames_2)
Intercept_table_2$name <- row.names(Intercept_table_2)

Intercept_raw_mean_2 <- c(b_abs_temperature, 
                          b_abs_terrestrial_temperature, 
                          b_abs_aquatic_temperature, 
                          b_abs_aquatic_salinity)

Intercept_raw_name_2 <- c(replicate(8000, "Temperature"), 
                          replicate(8000, "Terrestrial and Temperature"), 
                          replicate(8000, "Aquatic and Temperature"), 
                          replicate(8000, "Aquatic and Salinity"))

Intercept_raw_df_2 <- data.frame("Model" = Intercept_raw_name_2, 
                                 "Effect" = Intercept_raw_mean_2)

# Graph Code - Part 2

Intercept_Order_2 <- c("Aquatic and Salinity", "Aquatic and Temperature", "Terrestrial and Temperature", "Temperature")

density_intercept_2 <- Intercept_table_2 %>% mutate(name = fct_relevel(name, Intercept_Order_2)) %>%
                       ggplot() +
                       geom_density_ridges(data = Intercept_raw_df_2 %>% mutate(Model = fct_relevel(Model, Intercept_Order_2)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_2)[1], 1)), xmin = max(Intercept_raw_df_2$Effect)+0.001, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(Intercept_table_2)[1], 1)), xmin = min(Intercept_raw_df_2$Effect)-0.001, xmax = -0.2, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(Intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = expression("Effect Size (PRRD"["S"]*")"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-0.4, -0.4, -0.4, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 7)) +
                       scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D")) +
                       scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D")) +
                       coord_cartesian(xlim = c(-0.01, 1.25)) +
                       annotate('text',  x = 1.25, y = (seq(1, dim(Intercept_table_2)[1], 1)+0.4),
                       label= paste("italic(k)==", c(Intercept_table_2["Aquatic and Salinity", "K"], 
                                                     Intercept_table_2["Aquatic and Temperature", "K"], 
                                                     Intercept_table_2["Terrestrial and Temperature", "K"], 
                                                     Intercept_table_2["Temperature", "K"]), "~","(", 
                                                   c(Intercept_table_2["Aquatic and Salinity", "group_no"], 
                                                     Intercept_table_2["Aquatic and Temperature", "group_no"], 
                                                     Intercept_table_2["Terrestrial and Temperature", "group_no"], 
                                                     Intercept_table_2["Temperature", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(b_abs_aquatic_salinity)-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(b_abs_aquatic_temperature)-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(b_abs_terrestrial_temperature)-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(b_abs_temperature)-1)*100, 2), nsmall = 2), "%")), 
                                     x = rev(Intercept_table_2$estimate+0.2), y = (seq(1, dim(Intercept_table_2)[1], 1)+0.4)), size = 3.5)

density_intercept_2 #(400x400)

##### Supplementary Material Table #####
# Consistency Changes - Studies, Species and Effect Sizes Counts
Pre_Data <- read.csv("./Pre_Data.csv")

Measurement_Studies <- Pre_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Studies) <- Measurement_Studies$Measurement
colnames(Measurement_Studies) <- c("Measurement", "Study")

Measurement_Species <- Pre_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Species) <- Measurement_Species$Measurement
colnames(Measurement_Species) <- c("Measurement", "Species")

Measurement_Effects <- Pre_Data %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Effects) <- Measurement_Effects$Measurement
colnames(Measurement_Effects) <- c("Measurement", "Effect Sizes")

Measurement_Final_Counts <- Measurement_Studies %>% 
                            left_join(Measurement_Species, by = "Measurement") %>% 
                            left_join(Measurement_Effects, by = "Measurement")

write.csv(Measurement_Final_Counts, file = "./Measurement_Final_Counts.csv", row.names = FALSE)

# Category - Studies, Species and Effect Sizes Counts
Measurement_Category_Studies <- data %>% select("Study_ID", "Category") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Studies) <- Measurement_Category_Studies$Category
colnames(Measurement_Category_Studies) <- c("Category", "Study")

Measurement_Category_Species <- data %>% select("Scientific_Name", "Category") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Species) <- Measurement_Category_Species$Category
colnames(Measurement_Category_Species) <- c("Category", "Species")

Measurement_Category_Effects <- data %>% select("Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Effects) <- Measurement_Category_Effects$Category
colnames(Measurement_Category_Effects) <- c("Category", "Effect Sizes")

Measurement_Category_Final_Counts <- Measurement_Category_Studies %>% 
                                     left_join(Measurement_Category_Species, by = "Category") %>% 
                                     left_join(Measurement_Category_Effects, by = "Category")

write.csv(Measurement_Category_Final_Counts, file = "./Measurement_Category_Final_Counts.csv", row.names = FALSE)

# Measurements and their Categories
Behaviour_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Behavioural") %>%
                  select("Measurement") %>% unique ()

Biochemical_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Biochemical Assay") %>%
                    select("Measurement") %>% unique ()

Gene_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Gene Expression") %>%
             select("Measurement") %>% unique ()

Life_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Life-History Traits") %>%
             select("Measurement") %>% unique ()

Morphology_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Morphology") %>%
                   select("Measurement") %>% unique ()

Physiology_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Physiological") %>%
                   select("Measurement") %>% unique ()

Tolerance_Name <- data %>% select("Study_ID", "Category", "Measurement") %>% filter(`Category` == "Tolerance") %>%
                  select("Measurement") %>% unique ()

# Phylogenetic Tree with labels
labelled_tree <- tree

Scientific_Name_Effects <- data %>% select("Scientific_Name") %>% table() %>% data.frame()
rownames(Scientific_Name_Effects) <- Scientific_Name_Effects$Scientific_Name
colnames(Scientific_Name_Effects) <- c("Scientific_Name", "Effect_Sizes")
Scientific_Name_Effects <- Scientific_Name_Effects[c(labelled_tree$tip.label), ]

Scientific_Name_Studies <- data %>% select("Study_ID", "Scientific_Name") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Scientific_Name") %>% table() %>% data.frame()
rownames(Scientific_Name_Studies) <- Scientific_Name_Studies$Scientific_Name
colnames(Scientific_Name_Studies) <- c("Scientific_Name", "Study")
Scientific_Name_Studies <- Scientific_Name_Studies[c(labelled_tree$tip.label), ]

labelled_tree$tip.label <- paste(labelled_tree$tip.label, " ", Scientific_Name_Effects$Effect_Sizes, "(", Scientific_Name_Studies$Study, ")")
node.depth(labelled_tree, method = 2)
plot(labelled_tree, node.color = "#183357")

# Raw Data Tables

Raw_Overall <- data.frame("Overall" = c("MLMA"),
                          "Studies" = c(Intercept_study["Overall", "Study"]), 
                          "Species" = c(Intercept_table["Overall", "group_no"]), 
                          "Effect Sizes" = c(Intercept_table["Overall", "K"]),
                          "Estimate" = c(mean_abs_b),
                          "CI Low" = c(ci.abs[1]), 
                          "CI High" = c(ci.abs[2]))

Raw_Plasticity <- data.frame("Plasticity Category" = c("Acclimation", "Developmental Plasticity", "Transgenerational Effects"), 
                             "Studies" = c(Plasticity_study["Acclimation", "Study"], 
                                           Plasticity_study["Developmental Plasticity", "Study"], 
                                           Plasticity_study["Transgenerational Effects", "Study"]), 
                             "Species" = c(Plasticity_table["Acclimation", "group_no"], 
                                           Plasticity_table["Developmental Plasticity", "group_no"], 
                                           Plasticity_table["Transgenerational Effects", "group_no"]), 
                             "Effect Sizes" = c(Plasticity_table["Acclimation", "K"], 
                                                Plasticity_table["Developmental Plasticity", "K"], 
                                                Plasticity_table["Transgenerational Effects", "K"]), 
                             "Estimate" = c(mean_abs_b_plastic_acc, mean_abs_b_plastic_dev, mean_abs_b_plastic_trans), 
                             "CI Low" = c(ci.abs_acc[1], ci.abs_dev[1], ci.abs_trans[1]), 
                             "CI High" = c(ci.abs_acc[2], ci.abs_dev[2], ci.abs_trans[2]))

Raw_Trait <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-History Traits", 
                                                          "Morphology", "Physiological", "Tolerance"), 
                        "Studies" = c(Trait_study["Behavioural", "Study"], 
                                      Trait_study["Biochemical Assay", "Study"], 
                                      Trait_study["Gene Expression", "Study"], 
                                      Trait_study["Life-history Traits", "Study"], 
                                      Trait_study["Morphological", "Study"], 
                                      Trait_study["Physiological", "Study"], 
                                      Trait_study["Tolerance", "Study"]), 
                        "Species" = c(Trait_table["Behavioural", "group_no"], 
                                      Trait_table["Biochemical Assay", "group_no"], 
                                      Trait_table["Gene Expression", "group_no"], 
                                      Trait_table["Life-history Traits", "group_no"], 
                                      Trait_table["Morphological", "group_no"], 
                                      Trait_table["Physiological", "group_no"], 
                                      Trait_table["Tolerance", "group_no"]), 
                        "Effect Sizes" = c(Trait_table["Behavioural", "K"], 
                                           Trait_table["Biochemical Assay", "K"], 
                                           Trait_table["Gene Expression", "K"], 
                                           Trait_table["Life-history Traits", "K"], 
                                           Trait_table["Morphological", "K"], 
                                           Trait_table["Physiological", "K"], 
                                           Trait_table["Tolerance", "K"]), 
                        "Estimate" = c(mean_abs_b_trait_behavioural, mean_abs_b_trait_biochem, 
                                       mean_abs_b_trait_gene, mean_abs_b_trait_life, 
                                       mean_abs_b_trait_morphology, mean_abs_b_trait_physiological, 
                                       mean_abs_b_trait_tolerance), 
                        "CI Low" = c(ci.abs_behavioural[1], ci.abs_biochem[1], ci.abs_gene[1], 
                                     ci.abs_life[1], ci.abs_morphology[1], ci.abs_physiological[1], 
                                     ci.abs_tolerance[1]), 
                        "CI High" = c(ci.abs_behavioural[2], ci.abs_biochem[2], ci.abs_gene[2], 
                                      ci.abs_life[2], ci.abs_morphology[2], ci.abs_physiological[2], 
                                      ci.abs_tolerance[2]))

Raw_Taxonomy <- data.frame("Taxonomic Class" = c("Actinopteri", "Amphibia", "Branchiopoda", 
                                                 "Gastropoda", "Insecta"), 
                           "Studies" = c(Taxonomy_study["Actinopteri", "Study"], 
                                         Taxonomy_study["Amphibia", "Study"], 
                                         Taxonomy_study["Branchiopoda", "Study"], 
                                         Taxonomy_study["Gastropoda", "Study"], 
                                         Taxonomy_study["Insecta", "Study"]), 
                           "Species" = c(Taxonomy_table["Actinopteri", "group_no"], 
                                         Taxonomy_table["Amphibia", "group_no"], 
                                         Taxonomy_table["Branchiopoda", "group_no"], 
                                         Taxonomy_table["Gastropoda", "group_no"], 
                                         Taxonomy_table["Insecta", "group_no"]), 
                           "Effect Sizes" = c(Taxonomy_table["Actinopteri", "K"], 
                                              Taxonomy_table["Amphibia", "K"], 
                                              Taxonomy_table["Branchiopoda", "K"], 
                                              Taxonomy_table["Gastropoda", "K"], 
                                              Taxonomy_table["Insecta", "K"]), 
                           "Estimate" = c(mean_abs_b_taxonomy_actinopteri, mean_abs_b_taxonomy_amphibia, 
                                          mean_abs_b_taxonomy_branchiopoda, mean_abs_b_taxonomy_gastropoda, 
                                          mean_abs_b_taxonomy_insecta), 
                           "CI Low" = c(ci.abs_taxonomy_actinopteri[1], ci.abs_taxonomy_amphibia[1], 
                                        ci.abs_taxonomy_branchiopoda[1], ci.abs_taxonomy_gastropoda[1], 
                                        ci.abs_taxonomy_insecta[1]), 
                           "CI High" = c(ci.abs_taxonomy_actinopteri[2], ci.abs_taxonomy_amphibia[2], 
                                         ci.abs_taxonomy_branchiopoda[2], ci.abs_taxonomy_gastropoda[2], 
                                         ci.abs_taxonomy_insecta[2]))

Raw_Treatment <- data.frame("Treatment Type" = c("Diet", "Humidity", "Predator", 
                                                 "Salinity", "Temperature", "Water Level"), 
                            "Studies" = c(Treatment_study["Diet", "Study"], 
                                          Treatment_study["Humidity", "Study"], 
                                          Treatment_study["Predator", "Study"], 
                                          Treatment_study["Salinity", "Study"], 
                                          Treatment_study["Temperature", "Study"], 
                                          Treatment_study["Water Level", "Study"]), 
                            "Species" = c(Treatment_table["Diet", "group_no"], 
                                          Treatment_table["Humidity", "group_no"], 
                                          Treatment_table["Predator", "group_no"], 
                                          Treatment_table["Salinity", "group_no"], 
                                          Treatment_table["Temperature", "group_no"],
                                          Treatment_table["Water Level", "group_no"]), 
                            "Effect Sizes" = c(Treatment_table["Diet", "K"], 
                                               Treatment_table["Humidity", "K"], 
                                               Treatment_table["Predator", "K"], 
                                               Treatment_table["Salinity", "K"], 
                                               Treatment_table["Temperature", "K"],
                                               Treatment_table["Water Level", "K"]), 
                            "Estimate" = c(mean_abs_b_treatment_diet, mean_abs_b_treatment_humidity, 
                                           mean_abs_b_treatment_predator, mean_abs_b_treatment_salinity, 
                                           mean_abs_b_treatment_temperature, mean_abs_b_treatment_waterlevel), 
                            "CI Low" = c(ci.abs_treatment_diet[1], ci.abs_treatment_humidity[1], 
                                         ci.abs_treatment_predator[1], ci.abs_treatment_salinity[1], 
                                         ci.abs_treatment_temperature[1], ci.abs_treatment_waterlevel[1]), 
                            "CI High" = c(ci.abs_treatment_diet[2], ci.abs_treatment_humidity[2], 
                                          ci.abs_treatment_predator[2], ci.abs_treatment_salinity[2], 
                                          ci.abs_treatment_temperature[2], ci.abs_treatment_waterlevel[2]))

Raw_Terrestrial <- data.frame("Terrestrial" = c("MLMA"), 
                              "Studies" = c(Intercept_study ["Terrestrial", "Study"]), 
                              "Species" = c(Intercept_table["Terrestrial", "group_no"]), 
                              "Effect Sizes" = c(Intercept_table["Terrestrial", "K"]),
                              "Estimate" = c(mean_abs_b_terrestrial),
                              "CI Low" = c(ci.abs_terrestrial[1]), 
                              "CI High" = c(ci.abs_terrestrial[2]))

Raw_Terrestrial_Plasticity <- data.frame("Plasticity Category" = c("Acclimation", "Developmental Plasticity"), 
                                         "Studies" = c(Terrestrial_Plasticity_study["Acclimation", "Study"], 
                                                       Terrestrial_Plasticity_study["Developmental Plasticity", "Study"]), 
                                         "Species" = c(Terrestrial_Plasticity_table["Acclimation", "group_no"], 
                                                       Terrestrial_Plasticity_table["Developmental Plasticity", "group_no"]), 
                                         "Effect Sizes" = c(Terrestrial_Plasticity_table["Acclimation", "K"], 
                                                            Terrestrial_Plasticity_table["Developmental Plasticity", "K"]), 
                                         "Estimate" = c(mean_abs_b_terrestrial_plastic_acc, mean_abs_b_terrestrial_plastic_dev), 
                                         "CI Low" = c(ci.abs_terrestrial_acc[1], ci.abs_terrestrial_dev[1]), 
                                         "CI High" = c(ci.abs_terrestrial_acc[2], ci.abs_terrestrial_dev[2]))

Raw_Terrestrial_Trait <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", 
                                                                      "Morphology", "Tolerance"), 
                                    "Studies" = c(Terrestrial_Trait_study["Behavioural", "Study"], 
                                                  Terrestrial_Trait_study["Biochemical Assay", "Study"], 
                                                  Terrestrial_Trait_study["Gene Expression", "Study"], 
                                                  Terrestrial_Trait_study["Morphological", "Study"], 
                                                  Terrestrial_Trait_study["Tolerance", "Study"]), 
                                    "Species" = c(Terrestrial_Trait_table["Behavioural", "group_no"], 
                                                  Terrestrial_Trait_table["Biochemical Assay", "group_no"], 
                                                  Terrestrial_Trait_table["Gene Expression", "group_no"], 
                                                  Terrestrial_Trait_table["Morphological", "group_no"], 
                                                  Terrestrial_Trait_table["Tolerance", "group_no"]), 
                                    "Effect Sizes" = c(Terrestrial_Trait_table["Behavioural", "K"], 
                                                       Terrestrial_Trait_table["Biochemical Assay", "K"], 
                                                       Terrestrial_Trait_table["Gene Expression", "K"], 
                                                       Terrestrial_Trait_table["Morphological", "K"], 
                                                       Terrestrial_Trait_table["Tolerance", "K"]), 
                                    "Estimate" = c(mean_abs_b_terrestrial_trait_behavioural, mean_abs_b_terrestrial_trait_biochem, 
                                                   mean_abs_b_terrestrial_trait_gene, mean_abs_b_terrestrial_trait_morphology,
                                                   mean_abs_b_terrestrial_trait_tolerance), 
                                    "CI Low" = c(ci.abs_terrestrial_behavioural[1], ci.abs_terrestrial_biochem[1], ci.abs_terrestrial_gene[1], 
                                                 ci.abs_terrestrial_morphology[1], ci.abs_terrestrial_tolerance[1]), 
                                    "CI High" = c(ci.abs_terrestrial_behavioural[2], ci.abs_terrestrial_biochem[2], ci.abs_terrestrial_gene[2], 
                                                  ci.abs_terrestrial_morphology[2], ci.abs_terrestrial_tolerance[2]))

Raw_Terrestrial_Treatment <- data.frame("Treatment Type" = c("Humidity", "Temperature"), 
                                        "Studies" = c(Terrestrial_Treatment_study["Humidity", "Study"], 
                                                      Terrestrial_Treatment_study["Temperature", "Study"]), 
                                        "Species" = c(Terrestrial_Treatment_table["Humidity", "group_no"], 
                                                      Terrestrial_Treatment_table["Temperature", "group_no"]), 
                                        "Effect Sizes" = c(Terrestrial_Treatment_table["Humidity", "K"], 
                                                           Terrestrial_Treatment_table["Temperature", "K"]), 
                                        "Estimate" = c(mean_abs_b_terrestrial_treatment_humidity, mean_abs_b_terrestrial_treatment_temperature), 
                                        "CI Low" = c(ci.abs_terrestrial_treatment_humidity[1], ci.abs_terrestrial_treatment_temperature[1]), 
                                        "CI High" = c(ci.abs_terrestrial_treatment_humidity[2], ci.abs_terrestrial_treatment_temperature[2]))

Raw_Aquatic <- data.frame("Aquatic" = c("MLMA"), 
                          "Studies" = c(Intercept_study ["Aquatic", "Study"]), 
                          "Species" = c(Intercept_table["Aquatic", "group_no"]), 
                          "Effect Sizes" = c(Intercept_table["Aquatic", "K"]),
                          "Estimate" = c(mean_abs_b_aquatic),
                          "CI Low" = c(ci.abs_aquatic[1]), 
                          "CI High" = c(ci.abs_aquatic[2]))

Raw_Aquatic_Plasticity <- data.frame("Plasticity Category" = c("Acclimation", "Developmental Plasticity", "Transgenerational Effects"), 
                                     "Studies" = c(Aquatic_Plasticity_study["Acclimation", "Study"], 
                                                   Aquatic_Plasticity_study["Developmental Plasticity", "Study"], 
                                                   Aquatic_Plasticity_study["Transgenerational Effects", "Study"]), 
                                     "Species" = c(Aquatic_Plasticity_table["Acclimation", "group_no"], 
                                                   Aquatic_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                   Aquatic_Plasticity_table["Transgenerational Effects", "group_no"]), 
                                     "Effect Sizes" = c(Aquatic_Plasticity_table["Acclimation", "K"], 
                                                        Aquatic_Plasticity_table["Developmental Plasticity", "K"], 
                                                        Aquatic_Plasticity_table["Transgenerational Effects", "K"]), 
                                     "Estimate" = c(mean_abs_b_aquatic_plastic_acc, mean_abs_b_aquatic_plastic_dev, 
                                                    mean_abs_b_aquatic_plastic_trans), 
                                     "CI Low" = c(ci.abs_aquatic_acc[1], ci.abs_aquatic_dev[1], 
                                                  ci.abs_aquatic_trans[1]), 
                                     "CI High" = c(ci.abs_aquatic_acc[2], ci.abs_aquatic_dev[2], 
                                                   ci.abs_aquatic_trans[2]))

Raw_Aquatic_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay", "Life-History Traits", "Morphology", 
                                                                  "Physiological", "Tolerance"), 
                                "Studies" = c(Aquatic_Trait_study["Biochemical Assay", "Study"], 
                                              Aquatic_Trait_study["Life-history Traits", "Study"], 
                                              Aquatic_Trait_study["Morphological", "Study"], 
                                              Aquatic_Trait_study["Physiological", "Study"], 
                                              Aquatic_Trait_study["Tolerance", "Study"]), 
                                "Species" = c(Aquatic_Trait_table["Biochemical Assay", "group_no"], 
                                              Aquatic_Trait_table["Life-history Traits", "group_no"], 
                                              Aquatic_Trait_table["Morphological", "group_no"], 
                                              Aquatic_Trait_table["Physiological", "group_no"], 
                                              Aquatic_Trait_table["Tolerance", "group_no"]), 
                                "Effect Sizes" = c(Aquatic_Trait_table["Biochemical Assay", "K"], 
                                                   Aquatic_Trait_table["Life-history Traits", "K"], 
                                                   Aquatic_Trait_table["Morphological", "K"], 
                                                   Aquatic_Trait_table["Physiological", "K"], 
                                                   Aquatic_Trait_table["Tolerance", "K"]), 
                                "Estimate" = c(mean_abs_b_aquatic_trait_biochem, mean_abs_b_aquatic_trait_life, 
                                               mean_abs_b_aquatic_trait_morphology, mean_abs_b_aquatic_trait_physiological,
                                               mean_abs_b_aquatic_trait_tolerance), 
                                "CI Low" = c(ci.abs_aquatic_biochem[1], ci.abs_aquatic_life[1], ci.abs_aquatic_morphology[1], 
                                             ci.abs_aquatic_physiological[1], ci.abs_aquatic_tolerance[1]), 
                                "CI High" = c(ci.abs_aquatic_biochem[2], ci.abs_aquatic_life[2], ci.abs_aquatic_morphology[2], 
                                              ci.abs_aquatic_physiological[2], ci.abs_aquatic_tolerance[2]))

Raw_Aquatic_Treatment <- data.frame("Treatment Type" = c("Diet", "Predator", "Salinity", "Temperature", "Water Level"), 
                                    "Studies" = c(Aquatic_Treatment_study["Diet", "Study"], 
                                                  Aquatic_Treatment_study["Predator", "Study"], 
                                                  Aquatic_Treatment_study["Salinity", "Study"], 
                                                  Aquatic_Treatment_study["Temperature", "Study"], 
                                                  Aquatic_Treatment_study["Water Level", "Study"]), 
                                    "Species" = c(Aquatic_Treatment_table["Diet", "group_no"], 
                                                  Aquatic_Treatment_table["Predator", "group_no"], 
                                                  Aquatic_Treatment_table["Salinity", "group_no"], 
                                                  Aquatic_Treatment_table["Temperature", "group_no"], 
                                                  Aquatic_Treatment_table["Water Level", "group_no"]), 
                                    "Effect Sizes" = c(Aquatic_Treatment_table["Diet", "K"], 
                                                       Aquatic_Treatment_table["Predator", "K"], 
                                                       Aquatic_Treatment_table["Salinity", "K"], 
                                                       Aquatic_Treatment_table["Temperature", "K"], 
                                                       Aquatic_Treatment_table["Water Level", "K"]), 
                                    "Estimate" = c(mean_abs_b_aquatic_treatment_diet, mean_abs_b_aquatic_treatment_predator, 
                                                   mean_abs_b_aquatic_treatment_salinity, mean_abs_b_aquatic_treatment_temperature, 
                                                   mean_abs_b_aquatic_treatment_waterlevel), 
                                    "CI Low" = c(ci.abs_aquatic_treatment_diet[1], ci.abs_aquatic_treatment_predator[1], 
                                                 ci.abs_aquatic_treatment_salinity[1], ci.abs_aquatic_treatment_temperature[1], 
                                                 ci.abs_aquatic_treatment_waterlevel[1]), 
                                    "CI High" = c(ci.abs_aquatic_treatment_diet[2], ci.abs_aquatic_treatment_predator[2], 
                                                  ci.abs_aquatic_treatment_salinity[2], ci.abs_aquatic_treatment_temperature[2], 
                                                  ci.abs_aquatic_treatment_waterlevel[2]))

Raw_Temperature <- data.frame("Temperature" = c("MLMA"),
                              "Studies" = c(Intercept_study["Temperature", "Study"]), 
                              "Species" = c(Intercept_table["Temperature", "group_no"]), 
                              "Effect Sizes" = c(Intercept_table["Temperature", "K"]),
                              "Estimate" = c(mean_abs_b_temperature),
                              "CI Low" = c(ci.abs_temperature[1]), 
                              "CI High" = c(ci.abs_temperature[2]))

Raw_Temperature_Plasticity <- data.frame("Plasticity Category" = c("Acclimation", "Developmental Plasticity", "Transgenerational Effects"), 
                                         "Studies" = c(Temperature_Plasticity_study["Acclimation", "Study"], 
                                                       Temperature_Plasticity_study["Developmental Plasticity", "Study"], 
                                                       Temperature_Plasticity_study["Transgenerational Effects", "Study"]), 
                                         "Species" = c(Temperature_Plasticity_table["Acclimation", "group_no"], 
                                                       Temperature_Plasticity_table["Developmental Plasticity", "group_no"], 
                                                       Temperature_Plasticity_table["Transgenerational Effects", "group_no"]), 
                                         "Effect Sizes" = c(Temperature_Plasticity_table["Acclimation", "K"], 
                                                            Temperature_Plasticity_table["Developmental Plasticity", "K"], 
                                                            Temperature_Plasticity_table["Transgenerational Effects", "K"]), 
                                         "Estimate" = c(mean_abs_b_temperature_plastic_acc, mean_abs_b_temperature_plastic_dev, mean_abs_b_temperature_plastic_trans), 
                                         "CI Low" = c(ci.abs_temperature_acc[1], ci.abs_temperature_dev[1], ci.abs_temperature_trans[1]), 
                                         "CI High" = c(ci.abs_temperature_acc[2], ci.abs_temperature_dev[2], ci.abs_temperature_trans[2]))

Raw_Temperature_Trait <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", 
                                                                      "Morphology", "Tolerance"), 
                                    "Studies" = c(Temperature_Trait_study["Behavioural", "Study"], 
                                                  Temperature_Trait_study["Biochemical Assay", "Study"], 
                                                  Temperature_Trait_study["Morphological", "Study"], 
                                                  Temperature_Trait_study["Tolerance", "Study"]), 
                                    "Species" = c(Temperature_Trait_table["Behavioural", "group_no"], 
                                                  Temperature_Trait_table["Biochemical Assay", "group_no"], 
                                                  Temperature_Trait_table["Morphological", "group_no"], 
                                                  Temperature_Trait_table["Tolerance", "group_no"]), 
                                    "Effect Sizes" = c(Temperature_Trait_table["Behavioural", "K"], 
                                                       Temperature_Trait_table["Biochemical Assay", "K"], 
                                                       Temperature_Trait_table["Morphological", "K"], 
                                                       Temperature_Trait_table["Tolerance", "K"]), 
                                    "Estimate" = c(mean_abs_b_temperature_trait_behavioural, mean_abs_b_temperature_trait_biochem, 
                                                   mean_abs_b_temperature_trait_morphology, mean_abs_b_temperature_trait_tolerance), 
                                    "CI Low" = c(ci.abs_temperature_behavioural[1], ci.abs_temperature_biochem[1], 
                                                 ci.abs_temperature_morphology[1], ci.abs_temperature_tolerance[1]), 
                                    "CI High" = c(ci.abs_temperature_behavioural[2], ci.abs_temperature_biochem[2], 
                                                  ci.abs_temperature_morphology[2], ci.abs_temperature_tolerance[2]))

Raw_Temperature_Taxonomy <- data.frame("Taxonomic Class" = c("Actinopteri", "Insecta"), 
                                       "Studies" = c(Temperature_Taxonomy_study["Actinopteri", "Study"], 
                                                     Temperature_Taxonomy_study["Insecta", "Study"]), 
                                       "Species" = c(Temperature_Taxonomy_table["Actinopteri", "group_no"], 
                                                     Temperature_Taxonomy_table["Insecta", "group_no"]), 
                                       "Effect Sizes" = c(Temperature_Taxonomy_table["Actinopteri", "K"], 
                                                          Temperature_Taxonomy_table["Insecta", "K"]), 
                                       "Estimate" = c(mean_abs_b_temperature_taxonomy_actinopteri, mean_abs_b_temperature_taxonomy_insecta), 
                                       "CI Low" = c(ci.abs_temperature_taxonomy_actinopteri[1], ci.abs_temperature_taxonomy_insecta[1]), 
                                       "CI High" = c(ci.abs_temperature_taxonomy_actinopteri[2], ci.abs_temperature_taxonomy_insecta[2]))

Raw_Terrestrial_Temperature <- data.frame("Terrestrial and Temperature" = c("MLMA"), 
                                          "Studies" = c(Intercept_study ["Terrestrial and Temperature", "Study"]), 
                                          "Species" = c(Intercept_table["Terrestrial and Temperature", "group_no"]), 
                                          "Effect Sizes" = c(Intercept_table["Terrestrial and Temperature", "K"]),
                                          "Estimate" = c(mean_abs_b_terrestrial_temperature),
                                          "CI Low" = c(ci.abs_terrestrial_temperature[1]), 
                                          "CI High" = c(ci.abs_terrestrial_temperature[2]))

Raw_Terrestrial_Temperature_Plasticity <- data.frame("Plasticity Category" = c("Acclimation", "Developmental Plasticity"), 
                                                     "Studies" = c(Terrestrial_Temperature_Plasticity_study["Acclimation", "Study"], 
                                                                   Terrestrial_Temperature_Plasticity_study["Developmental Plasticity", "Study"]), 
                                                     "Species" = c(Terrestrial_Temperature_Plasticity_table["Acclimation", "group_no"], 
                                                                   Terrestrial_Temperature_Plasticity_table["Developmental Plasticity", "group_no"]), 
                                                     "Effect Sizes" = c(Terrestrial_Temperature_Plasticity_table["Acclimation", "K"], 
                                                                        Terrestrial_Temperature_Plasticity_table["Developmental Plasticity", "K"]), 
                                                     "Estimate" = c(mean_abs_b_terrestrial_temperature_plastic_acc, 
                                                                    mean_abs_b_terrestrial_temperature_plastic_dev), 
                                                     "CI Low" = c(ci.abs_terrestrial_temperature_acc[1], 
                                                                  ci.abs_terrestrial_temperature_dev[1]), 
                                                     "CI High" = c(ci.abs_terrestrial_temperature_acc[2], 
                                                                   ci.abs_terrestrial_temperature_dev[2]))

Raw_Terrestrial_Temperature_Trait <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay",  
                                                                                  "Morphology", "Tolerance"), 
                                                "Studies" = c(Terrestrial_Temperature_Trait_study["Behavioural", "Study"], 
                                                              Terrestrial_Temperature_Trait_study["Biochemical Assay", "Study"], 
                                                              Terrestrial_Temperature_Trait_study["Morphological", "Study"], 
                                                              Terrestrial_Temperature_Trait_study["Tolerance", "Study"]), 
                                                "Species" = c(Terrestrial_Temperature_Trait_table["Behavioural", "group_no"], 
                                                              Terrestrial_Temperature_Trait_table["Biochemical Assay", "group_no"], 
                                                              Terrestrial_Temperature_Trait_table["Morphological", "group_no"], 
                                                              Terrestrial_Temperature_Trait_table["Tolerance", "group_no"]), 
                                                "Effect Sizes" = c(Terrestrial_Temperature_Trait_table["Behavioural", "K"], 
                                                                   Terrestrial_Temperature_Trait_table["Biochemical Assay", "K"], 
                                                                   Terrestrial_Temperature_Trait_table["Morphological", "K"], 
                                                                   Terrestrial_Temperature_Trait_table["Tolerance", "K"]), 
                                                "Estimate" = c(mean_abs_b_terrestrial_temperature_trait_behavioural, mean_abs_b_terrestrial_temperature_trait_biochem, 
                                                               mean_abs_b_terrestrial_temperature_trait_morphology, mean_abs_b_terrestrial_temperature_trait_tolerance), 
                                                "CI Low" = c(ci.abs_terrestrial_temperature_behavioural[1], ci.abs_terrestrial_temperature_biochem[1],  
                                                             ci.abs_terrestrial_temperature_morphology[1], ci.abs_terrestrial_temperature_tolerance[1]), 
                                                "CI High" = c(ci.abs_terrestrial_temperature_behavioural[2], ci.abs_terrestrial_temperature_biochem[2],  
                                                              ci.abs_terrestrial_temperature_morphology[2], ci.abs_terrestrial_temperature_tolerance[2]))

Raw_Aquatic_Temperature <- data.frame("Aquatic and Temperature" = c("MLMA"), 
                                      "Studies" = c(Intercept_study ["Aquatic and Temperature", "Study"]), 
                                      "Species" = c(Intercept_table["Aquatic and Temperature", "group_no"]), 
                                      "Effect Sizes" = c(Intercept_table["Aquatic and Temperature", "K"]),
                                      "Estimate" = c(mean_abs_b_aquatic_temperature),
                                      "CI Low" = c(ci.abs_aquatic_temperature[1]), 
                                      "CI High" = c(ci.abs_aquatic_temperature[2]))

Raw_Aquatic_Temperature_Plasticity <- data.frame("Plasticity Category" = c("Acclimation"), 
                                                 "Studies" = c(Aquatic_Temperature_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                                 "Species" = c(Aquatic_Temperature_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                 "Effect Sizes" = c(length(Aquatic_Temperature_Plasticity_Data$Effect_Size_ID)), 
                                                 "Estimate" = c(mean_abs_b_aquatic_temperature_plastic), 
                                                 "CI Low" = c(ci.abs_aquatic_temperature_plastic[1]), 
                                                 "CI High" = c(ci.abs_aquatic_temperature_plastic[2]))

Raw_Aquatic_Temperature_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay", "Tolerance"), 
                                            "Studies" = c(Aquatic_Temperature_Trait_study["Biochemical Assay", "Study"], 
                                                          Aquatic_Temperature_Trait_study["Tolerance", "Study"]), 
                                            "Species" = c(Aquatic_Temperature_Trait_table["Biochemical Assay", "group_no"], 
                                                          Aquatic_Temperature_Trait_table["Tolerance", "group_no"]), 
                                            "Effect Sizes" = c(Aquatic_Temperature_Trait_table["Biochemical Assay", "K"], 
                                                               Aquatic_Temperature_Trait_table["Tolerance", "K"]), 
                                            "Estimate" = c(mean_abs_b_aquatic_temperature_trait_biochem, 
                                                           mean_abs_b_aquatic_temperature_trait_tolerance), 
                                            "CI Low" = c(ci.abs_aquatic_temperature_biochem[1], 
                                                         ci.abs_aquatic_temperature_tolerance[1]), 
                                            "CI High" = c(ci.abs_aquatic_temperature_biochem[2], 
                                                          ci.abs_aquatic_temperature_tolerance[2]))

Raw_Aquatic_Salinity <- data.frame("Aquatic and Salinity" = c("MLMA"), 
                                   "Studies" = c(Intercept_study["Aquatic and Salinity", "Study"]), 
                                   "Species" = c(Intercept_table["Aquatic and Salinity", "group_no"]), 
                                   "Effect Sizes" = c(Intercept_table["Aquatic and Salinity", "K"]),
                                   "Estimate" = c(mean_abs_b_aquatic_salinity),
                                   "CI Low" = c(ci.abs_aquatic_salinity[1]), 
                                   "CI High" = c(ci.abs_aquatic_salinity[2]))

Raw_Aquatic_Salinity_Plasticity <- data.frame("Plasticity Category" = c("Acclimation"), 
                                              "Studies" = c(Aquatic_Salinity_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                              "Species" = c(Aquatic_Salinity_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                              "Effect Sizes" = c(length(Aquatic_Salinity_Plasticity_Data$Effect_Size_ID)), 
                                              "Estimate" = c(mean_abs_b_aquatic_salinity_plastic), 
                                              "CI Low" = c(ci.abs_aquatic_salinity_plastic[1]), 
                                              "CI High" = c(ci.abs_aquatic_salinity_plastic[2]))

Raw_Aquatic_Salinity_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay"), 
                                         "Studies" = c(Aquatic_Salinity_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                         "Species" = c(Aquatic_Salinity_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                         "Effect Sizes" = c(length(Aquatic_Salinity_Trait_Data$Effect_Size_ID)), 
                                         "Estimate" = c(mean_abs_b_aquatic_salinity_trait), 
                                         "CI Low" = c(ci.abs_aquatic_salinity_trait[1]), 
                                         "CI High" = c(ci.abs_aquatic_salinity_trait[2]))

write.csv(Raw_Overall, file = "./Raw_Overall.csv", row.names = FALSE)
write.csv(Raw_Plasticity, file = "./Raw_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Trait, file = "./Raw_Trait.csv", row.names = FALSE)
write.csv(Raw_Taxonomy, file = "./Raw_Taxonomy.csv", row.names = FALSE)
write.csv(Raw_Treatment, file = "./Raw_Treatment.csv", row.names = FALSE)
write.csv(Raw_Terrestrial, file = "./Raw_Terrestrial.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Plasticity, file = "./Raw_Terrestrial_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Trait, file = "./Raw_Terrestrial_Trait.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Treatment, file = "./Raw_Terrestrial_Treatment.csv", row.names = FALSE)
write.csv(Raw_Aquatic, file = "./Raw_Aquatic.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Plasticity, file = "./Raw_Aquatic_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Trait, file = "./Raw_Aquatic_Trait.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Treatment, file = "./Raw_Aquatic_Treatment.csv", row.names = FALSE)
write.csv(Raw_Temperature, file = "./Raw_Temperature.csv", row.names = FALSE)
write.csv(Raw_Temperature_Plasticity, file = "./Raw_Temperature_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Temperature_Trait, file = "./Raw_Temperature_Trait.csv", row.names = FALSE)
write.csv(Raw_Temperature_Taxonomy, file = "./Raw_Temperature_Taxonomy.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Temperature, file = "./Raw_Terrestrial_Temperature.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Temperature_Plasticity, file = "./Raw_Terrestrial_Temperature_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Temperature_Trait, file = "./Raw_Terrestrial_Temperature_Trait.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Temperature, file = "./Raw_Aquatic_Temperature.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Temperature_Plasticity, file = "./Raw_Aquatic_Temperature_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Temperature_Trait, file = "./Raw_Aquatic_Temperature_Trait.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Salinity, file = "./Raw_Aquatic_Salinity.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Salinity_Plasticity, file = "./Raw_Aquatic_Salinity_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Salinity_Trait, file = "./Raw_Aquatic_Salinity_Trait.csv", row.names = FALSE)


# Heterogeneity Tables
Heterogeneity_Subsets <- data.frame("Models" = c("Overall", "Terrestrial", "Aquatic", "Temperature", 
                                                 "Terrestrial and Temperature", "Aquatic and Temperature", "Aquatic and Salinity"), 
                                    "Measurement" = c(overall_i2["sd_Measurement__Intercept", 1], terrestrial_i2["sd_Measurement__Intercept", 1], 
                                                      aquatic_i2["sd_Measurement__Intercept", 1], NA, NA, NA, NA),
                                    "Observational" = c(overall_i2["sd_obs__Intercept", 1], terrestrial_i2["sd_obs__Intercept", 1], 
                                                        aquatic_i2["sd_obs__Intercept", 1], temperature_i2["sd_obs__Intercept", 1], 
                                                        terrestrial_temperature_i2["sd_obs__Intercept", 1], aquatic_temperature_i2["sd_obs__Intercept", 1], 
                                                        aquatic_salinity_i2["sd_obs__Intercept", 1]), 
                                    "Phylogenetic Relatedness" = c(overall_i2["sd_phylo__Intercept", 1], terrestrial_i2["sd_phylo__Intercept", 1], 
                                                                   aquatic_i2["sd_phylo__Intercept", 1], temperature_i2["sd_phylo__Intercept", 1], 
                                                                   terrestrial_temperature_i2["sd_phylo__Intercept", 1], aquatic_temperature_i2["sd_phylo__Intercept", 1], 
                                                                   aquatic_salinity_i2["sd_phylo__Intercept", 1]), 
                                    "Study" = c(overall_i2["sd_Study_ID__Intercept", 1], terrestrial_i2["sd_Study_ID__Intercept", 1], 
                                                aquatic_i2["sd_Study_ID__Intercept", 1], temperature_i2["sd_Study_ID__Intercept", 1], 
                                                terrestrial_temperature_i2["sd_Study_ID__Intercept", 1], aquatic_temperature_i2["sd_Study_ID__Intercept", 1], 
                                                aquatic_salinity_i2["sd_Study_ID__Intercept", 1]), 
                                    "Total" = c(overall_i2["i2_total", 1], terrestrial_i2["i2_total", 1], 
                                                aquatic_i2["i2_total", 1], temperature_i2["i2_total", 1], 
                                                terrestrial_temperature_i2["i2_total", 1], aquatic_temperature_i2["i2_total", 1], 
                                                aquatic_salinity_i2["i2_total", 1]))

Heterogeneity_Plasticity <- data.frame("Models" = c("Overall", "Terrestrial", "Aquatic", "Temperature", 
                                                    "Terrestrial and Temperature", "Aquatic and Temperature", "Aquatic and Salinity"), 
                                       "Acclimation" = c(overall_i2_plastic["sd_Acclimation", 1], terrestrial_i2_plastic["sd_Acclimation", 1], 
                                                         aquatic_i2_plastic["sd_Acclimation", 1], temperature_i2_plastic["sd_Acclimation", 1], 
                                                         terrestrial_temperature_i2_plastic["sd_Acclimation", 1], aquatic_temperature_plastic_i2["sd_obs__Intercept", 1], 
                                                         aquatic_salinity_plastic_i2["sd_obs__Intercept", 1]),
                                       "Developmental Plasticity" = c(overall_i2_plastic["sd_Developmental", 1], terrestrial_i2_plastic["sd_Developmental", 1], 
                                                                      aquatic_i2_plastic["sd_Developmental", 1], temperature_i2_plastic["sd_Developmental", 1], 
                                                                      terrestrial_temperature_i2_plastic["sd_Developmental", 1], NA, NA), 
                                       "Transgenerational Effects" = c(overall_i2_plastic["sd_Transgenerational", 1], NA, 
                                                                       aquatic_i2_plastic["sd_Transgenerational", 1], temperature_i2_plastic["sd_Transgenerational", 1], 
                                                                       NA, NA, NA), 
                                       "Measurement" = c(overall_i2_plastic["sd_Measurement__Intercept", 1], terrestrial_i2_plastic["sd_Measurement__Intercept", 1], 
                                                         aquatic_i2_plastic["sd_Measurement__Intercept", 1], NA, NA, NA, NA),
                                       "Phylogenetic Relatedness" = c(overall_i2_plastic["sd_phylo__Intercept", 1], terrestrial_i2_plastic["sd_phylo__Intercept", 1], 
                                                                      aquatic_i2_plastic["sd_phylo__Intercept", 1], temperature_i2_plastic["sd_phylo__Intercept", 1], 
                                                                      terrestrial_temperature_i2_plastic["sd_phylo__Intercept", 1], aquatic_temperature_plastic_i2["sd_phylo__Intercept", 1], 
                                                                      aquatic_salinity_plastic_i2["sd_phylo__Intercept", 1]), 
                                       "Study" = c(overall_i2_plastic["sd_Study_ID__Intercept", 1], terrestrial_i2_plastic["sd_Study_ID__Intercept", 1], 
                                                   aquatic_i2_plastic["sd_Study_ID__Intercept", 1], temperature_i2_plastic["sd_Study_ID__Intercept", 1], 
                                                   terrestrial_temperature_i2_plastic["sd_Study_ID__Intercept", 1], aquatic_temperature_plastic_i2["sd_Study_ID__Intercept", 1], 
                                                   aquatic_salinity_plastic_i2["sd_Study_ID__Intercept", 1]), 
                                       "Total" = c(overall_i2_plastic["i2_total", 1], terrestrial_i2_plastic["i2_total", 1], 
                                                   aquatic_i2_plastic["i2_total", 1], temperature_i2_plastic["i2_total", 1], 
                                                   terrestrial_temperature_i2_plastic["i2_total", 1], aquatic_temperature_plastic_i2["i2_total", 1], 
                                                   aquatic_salinity_plastic_i2["i2_total", 1]))

Heterogeneity_Trait <- data.frame("Models" = c("Overall", "Terrestrial", "Aquatic", "Temperature", 
                                               "Terrestrial and Temperature", "Aquatic and Temperature", "Aquatic and Salinity"), 
                                 "Behavioural" = c(overall_i2_trait["sd_Behavioural", 1], terrestrial_i2_trait["sd_Behavioural", 1], 
                                                   NA, temperature_i2_trait["sd_Behavioural", 1], 
                                                   terrestrial_temperature_i2_trait["sd_Behavioural", 1], NA, NA),
                                 "Biochemical Assay" = c(overall_i2_trait["sd_BiochemicalAssay", 1], terrestrial_i2_trait["sd_BiochemicalAssay", 1], 
                                                         aquatic_i2_trait["sd_BiochemicalAssay", 1], temperature_i2_trait["sd_BiochemicalAssay", 1], 
                                                         terrestrial_temperature_i2_trait["sd_BiochemicalAssay", 1], aquatic_temperature_i2_trait["sd_BiochemicalAssay", 1], 
                                                         aquatic_salinity_trait_i2["sd_obs__Intercept", 1]), 
                                 "Gene Expression" = c(overall_i2_trait["sd_GeneExpression", 1], terrestrial_i2_trait["sd_GeneExpression", 1], 
                                                       NA, NA, NA, NA, NA), 
                                 "Life-History Traits" = c(overall_i2_trait["sd_LifeHistory", 1], NA, 
                                                           aquatic_i2_trait["sd_Life.HistoryTraits", 1], NA, NA, NA, NA),
                                 "Morphology" = c(overall_i2_trait["sd_Morphology", 1], terrestrial_i2_trait["sd_Morphology", 1], 
                                                  aquatic_i2_trait["sd_Morphology", 1], temperature_i2_trait["sd_Morphology", 1], 
                                                  terrestrial_temperature_i2_trait["sd_Morphology", 1], NA, NA), 
                                 "Physiological" = c(overall_i2_trait["sd_Physiological", 1], NA, 
                                                     aquatic_i2_trait["sd_Physiological", 1], NA, NA, NA, NA),
                                 "Tolerance" = c(overall_i2_trait["sd_Tolerance", 1], terrestrial_i2_trait["sd_Tolerance", 1], 
                                                 aquatic_i2_trait["sd_Tolerance", 1], temperature_i2_trait["sd_Tolerance", 1], 
                                                 terrestrial_temperature_i2_trait["sd_Tolerance", 1], aquatic_temperature_i2_trait["sd_Tolerance", 1], 
                                                 NA),
                                 "Measurement" = c(overall_i2_trait["sd_Measurement__Intercept", 1], terrestrial_i2_trait["sd_Measurement__Intercept", 1], 
                                                   aquatic_i2_trait["sd_Measurement__Intercept", 1], NA, NA, NA, NA),
                                 "Phylogenetic Relatedness" = c(overall_i2_trait["sd_phylo__Intercept", 1], terrestrial_i2_trait["sd_phylo__Intercept", 1], 
                                                                aquatic_i2_trait["sd_phylo__Intercept", 1], temperature_i2_trait["sd_phylo__Intercept", 1], 
                                                                terrestrial_temperature_i2_trait["sd_phylo__Intercept", 1], aquatic_temperature_i2_trait["sd_phylo__Intercept", 1], 
                                                                aquatic_salinity_trait_i2["sd_phylo__Intercept", 1]), 
                                 "Study" = c(overall_i2_trait["sd_Study_ID__Intercept", 1], terrestrial_i2_trait["sd_Study_ID__Intercept", 1], 
                                             aquatic_i2_trait["sd_Study_ID__Intercept", 1], temperature_i2_trait["sd_Study_ID__Intercept", 1], 
                                             terrestrial_temperature_i2_trait["sd_Study_ID__Intercept", 1], aquatic_temperature_i2_trait["sd_Study_ID__Intercept", 1], 
                                             aquatic_salinity_trait_i2["sd_Study_ID__Intercept", 1]), 
                                 "Total" = c(overall_i2_trait["i2_total", 1], terrestrial_i2_trait["i2_total", 1], 
                                             aquatic_i2_trait["i2_total", 1], temperature_i2_trait["i2_total", 1], 
                                             terrestrial_temperature_i2_trait["i2_total", 1], aquatic_temperature_i2_trait["i2_total", 1], 
                                             aquatic_salinity_trait_i2["i2_total", 1]))

Heterogeneity_Tax <- data.frame("Models" = c("Overall","Temperature"),
                                "Actinopteri" = c(overall_i2_taxonomy["sd_Actinopteri", 1], temperature_i2_taxonomy["sd_Actinopteri", 1]),
                                "Amphibia" = c(overall_i2_taxonomy["sd_Amphibia", 1], NA),
                                "Branchiopoda" = c(overall_i2_taxonomy["sd_Branchiopoda", 1], NA),
                                "Gastropoda" = c(overall_i2_taxonomy["sd_Gastropoda", 1], NA),
                                "Insecta" = c(overall_i2_taxonomy["sd_Insecta", 1], temperature_i2_taxonomy["sd_Insecta", 1]),
                                "Measurement" = c(overall_i2_taxonomy["sd_Measurement__Intercept", 1], NA),
                                "Phylogenetic Relatedness" = c(overall_i2_taxonomy["sd_phylo__Intercept", 1], temperature_i2_taxonomy["sd_phylo__Intercept", 1]), 
                                "Study" = c(overall_i2_taxonomy["sd_Study_ID__Intercept", 1], temperature_i2_taxonomy["sd_Study_ID__Intercept", 1]), 
                                "Total" = c(overall_i2_taxonomy["i2_total", 1], temperature_i2_taxonomy["i2_total", 1]))

Heterogeneity_Treatment <- data.frame("Models" = c("Overall","Terrestrial", "Aquatic"),
                                      "Diet" = c(overall_i2_treatment["sd_Diet", 1], NA, 
                                                 aquatic_i2_treatment["sd_Diet", 1]),
                                      "Humidity" = c(overall_i2_treatment["sd_Humidity", 1], terrestrial_i2_treatment["sd_Humidity", 1], 
                                                     NA),
                                      "Predator" = c(overall_i2_treatment["sd_Predator", 1], NA, 
                                                     aquatic_i2_treatment["sd_Predator", 1]),
                                      "Salinity" = c(overall_i2_treatment["sd_Salinity", 1], NA, 
                                                     aquatic_i2_treatment["sd_Salinity", 1]),
                                      "Temperature" = c(overall_i2_treatment["sd_Temperature", 1], terrestrial_i2_treatment["sd_Temperature", 1], 
                                                        aquatic_i2_treatment["sd_Temperature", 1]),
                                      "Water Level" = c(overall_i2_treatment["sd_WaterLevel", 1], NA, 
                                                        aquatic_i2_treatment["sd_WaterLevel", 1]),
                                      "Measurement" = c(overall_i2_treatment["sd_Measurement__Intercept", 1], terrestrial_i2_treatment["sd_Measurement__Intercept", 1], 
                                                        aquatic_i2_treatment["sd_Measurement__Intercept", 1]),
                                      "Phylogenetic Relatedness" = c(overall_i2_treatment["sd_phylo__Intercept", 1], terrestrial_i2_treatment["sd_phylo__Intercept", 1], 
                                                                     aquatic_i2_treatment["sd_phylo__Intercept", 1]), 
                                      "Study" = c(overall_i2_treatment["sd_Study_ID__Intercept", 1], terrestrial_i2_treatment["sd_Study_ID__Intercept", 1], 
                                                  aquatic_i2_treatment["sd_Study_ID__Intercept", 1]), 
                                      "Total" = c(overall_i2_treatment["i2_total", 1], terrestrial_i2_treatment["i2_total", 1], 
                                                  aquatic_i2_treatment["i2_total", 1]))

write.csv(Heterogeneity_Subsets, file = "./Heterogeneity_Subsets.csv", row.names = FALSE)
write.csv(Heterogeneity_Plasticity, file = "./Heterogeneity_Plasticity.csv", row.names = FALSE)
write.csv(Heterogeneity_Trait, file = "./Heterogeneity_Trait.csv", row.names = FALSE)
write.csv(Heterogeneity_Tax, file = "./Heterogeneity_Tax.csv", row.names = FALSE)
write.csv(Heterogeneity_Treatment, file = "./Heterogeneity_Treatment.csv", row.names = FALSE)
