rm(list = ls())
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
if (!require("metaAidR")) devtools::install_gitgub("daniel1noble/metaAidR", force = TRUE)
if (!require("orchaRd")) devtools::install_gitgub("daniel1noble/orchaRd", force = TRUE)
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
               robumeta, ggpmisc, ggpubr)

source("./4_Laboratory_Plasticity/3_Data_Analysis/1_R_code/func.R")

# Importing Data Set
data <- read.csv("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/data/Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

Temp_Subset_Data <- data %>% filter(!is.na(Effect_Size_Type_Adjusted))
Temp_Subset_Data <- Temp_Subset_Data %>% filter(`Type` == "Temperature")
Temp_Subset_Data$obs <- 1:nrow(Temp_Subset_Data)
Temp_Species <- Temp_Subset_Data %>% select("phylo") %>% unique()

Sal_Subset_Data <- data %>% filter(!is.na(Effect_Size_Type_Adjusted))
Sal_Subset_Data <- Sal_Subset_Data %>% filter(`Type` == "Salinity")
Sal_Subset_Data$obs <- 1:nrow(Sal_Subset_Data)
Sal_Species <- Sal_Subset_Data %>% select("phylo") %>% unique()

# Phylogenetic covariance matrix
tree <- ape::read.tree("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/phylogeny/tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

Temp_A <- as.data.frame(A)
Temp_A <- Temp_A[c(Temp_Species$phylo), c(Temp_Species$phylo)]
Temp_A <- as.matrix(Temp_A)

Temp_A_cor <- as.data.frame(A_cor)
Temp_A_cor <- Temp_A_cor[c(Temp_Species$phylo), c(Temp_Species$phylo)]
Temp_A_cor <- as.matrix(Temp_A_cor)

Sal_A <- as.data.frame(A)
Sal_A <- Sal_A[c(Sal_Species$phylo), c(Sal_Species$phylo)]
Sal_A <- as.matrix(Sal_A)

Sal_A_cor <- as.data.frame(A_cor)
Sal_A_cor <- Sal_A_cor[c(Sal_Species$phylo), c(Sal_Species$phylo)]
Sal_A_cor <- as.matrix(Sal_A_cor)

priors <-  prior(student_t(3, 0, 20), class = "sd")

##### Publication Bias - Overall #####
# Importing Model
Model <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/overall_model.rds")

b_overall <- as_draws_df(Model, variable = "b_Intercept")
b_overall <- data.frame(b_overall$b_Intercept)

sd <- as_draws_df(Model, variable = c("sd_Measurement__Intercept", "sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd <- data.frame("sd_Measurement__Intercept" = sd$sd_Measurement__Intercept, 
                 "sd_obs__Intercept" = sd$sd_obs__Intercept, 
                 "sd_phylo__Intercept" = sd$sd_phylo__Intercept, 
                 "sd_Study_ID__Intercept" = sd$sd_Study_ID__Intercept)

b_abs <- folded_norm(b_overall[,1], rowSums(sd))
mean_abs_b <- mean(b_abs)
mode_abs_b <- posterior.mode(as.mcmc(b_abs))
ci.abs <- HPDinterval(as.mcmc(b_abs))
overall_i2 <- i2(sd, data$Variance_Adjusted) 

# Residuals
Effect_Sizes <- data %>% select(c("obs", "Effect_Size_Adjusted"))

Predict_Model <- data.frame(predict(Model))
Predict_Model$obs <- 1:nrow(Predict_Model)
Residuals <- Predict_Model %>% 
             left_join(Effect_Sizes, by = "obs")
Residuals <- Residuals %>% mutate(Residuals = `Effect_Size_Adjusted` - `Estimate`)
Residuals <- Residuals[, c("Residuals", "Est.Error")]

# Funnel Plot
Funnel_Plot <- funnel(Residuals$Residuals, Residuals$Est.Error, yaxis = "seinv",
                      ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                      pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("top", legend = c("0.05 < p ≤ 1.00", "0 < p ≤ 0.05", "Studies"), 
       pch = c(22, 22, 21, 21), pt.bg = c("#FFFFFF","#D3DDEB", "#183357"), box.lwd = 2)

# Publication Year Graph
Graph_Data <- data
Graph_Data <- Graph_Data %>% mutate(n_category = ifelse(n_.P1T1. <= 15, "15", 
                                                 ifelse(n_.P1T1. > 15 & n_.P1T1. <= 30, "30", 
                                                 ifelse(n_.P1T1. > 30 & n_.P1T1. <= 45, "45", "> 45"))))


Publication_Graph <- ggplot(Graph_Data, aes(x = Publication_Year, y = abs(Effect_Size_Adjusted))) + 
                     geom_point(aes(x = Publication_Year, y = abs(Effect_Size_Adjusted), 
                                    size = fct_relevel(n_category, c("15", "30", "45", "> 45"))), 
                                    shape = 21, fill = "#4292c6", alpha = 0.5) + 
                     labs(x = "Publication Year", y = "Effect Size (PRRD)", 
                          size = "Sample Size") +
                     theme_bw() +
                     theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                     theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                     theme(legend.position = "bottom", legend.direction = "horizontal") + 
                     geom_hline(yintercept = mode_abs_b[1], lty = 2) + 
                     geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                     stat_poly_eq(formula = y ~ x, 
                     aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                     parse = TRUE)

Publication_Graph #(750x500)

# Time-lag Bias 
run <- FALSE
system.time( #  5ish minutes
  if(run){
    Year_Precision_rma <- metafor::rma.mv(abs(Effect_Size_Adjusted), V = Variance_Adjusted, test = "t", dfs = "contain",
                                          mods = ~ Publication_Year_Z + Precision - 1,
                                          random = list(~1|phylo, ~1|Study_ID, ~1|Measurement, ~1|obs), 
                                          R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                          control=list(rel.tol=1e-9))
    saveRDS(Year_Precision_rma, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_rma")
  } else {
    Year_Precision_rma <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_rma")})

Year_Precision_rma_Estimates <- data.frame(estimate = Year_Precision_rma$b, ci.lb = Year_Precision_rma$ci.lb, 
                                    ci.ub = Year_Precision_rma$ci.ub)
Year_Precision_rma_i2 <- data.frame(round(orchaRd::i2_ml(Year_Precision_rma), 2))

# cooks Distance (metafor)
run <- FALSE
system.time( #  5ish minutes
  if(run){
    Cooks_Model <- metafor::rma.mv(Effect_Size_Adjusted, V = Variance_Adjusted, test = "t", dfs = "contain",
                                   random = list(~1|phylo, ~1|Study_ID, ~1|Measurement, ~1|obs), 
                                   R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                   control=list(rel.tol=1e-9))
    saveRDS(Cooks_Model, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model")
  } else {
    Cooks_Model <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model")})

Cooks_Model_Estimates <- data.frame(estimate = Cooks_Model$b, ci.lb = Cooks_Model$ci.lb, 
                                    ci.ub = Cooks_Model$ci.ub)
Cooks_Model_i2 <- data.frame(round(orchaRd::i2_ml(Cooks_Model), 2))

#run <- FALSE
#system.time( # hasn't been ran yet
#  if(run){
#    Overall_Cooks <- cooks.distance(Cooks_Model)
#    saveRDS(Overall_Cooks, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks")
#  } else {
#    Overall_Cooks <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks")})

#dev.off()
#Cooks_Plot <- plot(Overall_Cooks, type = "o", pch = 21, xlab = "Observed Outcome", 
#                   ylab = "Cook's Distance", bg = "#183357")
#box(lwd = 2)

# Untransformed Model
system.time(  # 35ish minutes
  Untransformed <- brms::brm(Effect_Size | se(sqrt(Variance)) 
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
                       file = "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/untransformed_model",
                       file_refit = "on_change"))

b_untransformed <- as_draws_df(Untransformed, variable = "b_Intercept")
b_untransformed <- data.frame(b_untransformed$b_Intercept)

sd_untransformed <- as_draws_df(Untransformed, variable = c("sd_Measurement__Intercept", "sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_untransformed <- data.frame("sd_Measurement__Intercept" = sd_untransformed$sd_Measurement__Intercept, 
                               "sd_obs__Intercept" = sd_untransformed$sd_obs__Intercept, 
                               "sd_phylo__Intercept" = sd_untransformed$sd_phylo__Intercept, 
                               "sd_Study_ID__Intercept" = sd_untransformed$sd_Study_ID__Intercept)

mean_b_untransformed <-  sapply(b_untransformed, mean)
ci_b_untransformed <- as.vector(HPDinterval(as.mcmc(b_untransformed)))
pMCMC_b_untransformed <- 2*(1 - max(table(b_untransformed<0) / nrow(b_untransformed)))

b_abs_untransformed <- folded_norm(b_untransformed[,1], rowSums(sd_untransformed))
mean_abs_b_untransformed <- mean(b_abs_untransformed)
mode_abs_b_untransformed <- posterior.mode(as.mcmc(b_abs_untransformed))
ci.abs_untransformed <- HPDinterval(as.mcmc(b_abs_untransformed))

overall_i2_untransformed <- i2(sd_untransformed, data$Variance) 

##### Publication Bias - Temperature #####
# Importing Model
Model_Temp <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/temp_model.rds")

b_overall_temp <- as_draws_df(Model_Temp, variable = "b_Intercept")
b_overall_temp <- data.frame(b_overall_temp$b_Intercept)

sd_temp <- as_draws_df(Model_Temp, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_temp <- data.frame("sd_obs__Intercept" = sd_temp$sd_obs__Intercept, 
                      "sd_phylo__Intercept" = sd_temp$sd_phylo__Intercept, 
                      "sd_Study_ID__Intercept" = sd_temp$sd_Study_ID__Intercept)

b_abs_temp <- folded_norm(b_overall_temp[,1], rowSums(sd_temp))
mean_abs_b_temp <- mean(b_abs_temp)
mode_abs_b_temp <- posterior.mode(as.mcmc(b_abs_temp))
ci.abs_temp <- HPDinterval(as.mcmc(b_abs_temp))
overall_temp_i2 <- i2(sd_temp, Temp_Subset_Data$Variance_Type_Adjusted) 


# Residuals
Effect_Sizes_Temp <- Temp_Subset_Data %>% select(c("obs", "Effect_Size_Type_Adjusted"))

Predict_Model_Temp <- data.frame(predict(Model_Temp))
Predict_Model_Temp$obs <- 1:nrow(Predict_Model_Temp)
Residuals_Temp <- Predict_Model_Temp %>% 
                  left_join(Effect_Sizes_Temp, by = "obs")
Residuals_Temp <- Residuals_Temp %>% mutate(Residuals_Temp = `Effect_Size_Type_Adjusted` - `Estimate`)
Residuals_Temp <- Residuals_Temp[, c("Residuals_Temp", "Est.Error")]

# Funnel Plot
Funnel_Plot_Temp <- funnel(Residuals_Temp$Residuals_Temp, Residuals_Temp$Est.Error, yaxis = "seinv",
                           ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                           pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
     legend("top", legend = c("0.05 < p ≤ 1.00", "0 < p ≤ 0.05", "Studies"), 
     pch = c(22, 22, 21, 21), pt.bg = c("#FFFFFF","#D3DDEB", "#183357"), box.lwd = 2)

# Publication Year Graph
Graph_Data_Temp <- Temp_Subset_Data
Graph_Data_Temp <- Graph_Data_Temp %>% mutate(n_category = ifelse(n_.P1T1. <= 15, "15", 
                                                           ifelse(n_.P1T1. > 15 & n_.P1T1. <= 30, "30", 
                                                           ifelse(n_.P1T1. > 30 & n_.P1T1. <= 45, "45", "> 45"))))


Publication_Graph_Temp <- ggplot(Graph_Data_Temp, aes(x = Publication_Year, y = abs(Effect_Size_Type_Adjusted))) + 
                          geom_point(aes(x = Publication_Year, y = abs(Effect_Size_Type_Adjusted), 
                          size = fct_relevel(n_category, c("15", "30", "45", "> 45"))), 
                          shape = 21, fill = "#4292c6", alpha = 0.5) + 
                          labs(x = "Publication Year", y = expression("Effect Size (PRRD"["S"]*")"), 
                               size = "Sample Size") +
                          theme_bw() +
                          theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                          theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                          theme(legend.position = "bottom", legend.direction = "horizontal") + 
                          geom_hline(yintercept = mode_abs_b_temp[1], lty = 2) + 
                          geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                          stat_poly_eq(formula = y ~ x, 
                                       aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                       parse = TRUE) #+
                          #coord_cartesian(xlim = c(0, 25), 
                          #                ylim = c(-5, 5))

Publication_Graph_Temp #(750x500)

# Time-lag Bias
run <- FALSE
system.time( #  5ish minutes
  if(run){
    Year_Precision_Temp_rma <- metafor::rma.mv(abs(Effect_Size_Type_Adjusted), V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                               mods = ~ Publication_Year_Z + Precision - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                               R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Year_Precision_Temp_rma, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_Temp_rma")
  } else {
    Year_Precision_Temp_rma <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_Temp_rma")})

Year_Precision_Temp_rma_Estimates <- data.frame(estimate = Year_Precision_Temp_rma$b, ci.lb = Year_Precision_Temp_rma$ci.lb, 
                                           ci.ub = Year_Precision_Temp_rma$ci.ub)
Year_Precision_Temp_rma_i2 <- data.frame(round(orchaRd::i2_ml(Year_Precision_Temp_rma), 2))

# Cooks Distance (Metafor)
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Cooks_Model_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                        R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                        control=list(rel.tol=1e-9))
    saveRDS(Cooks_Model_Temp, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model_temp")
  } else {
    Cooks_Model_Temp <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model_temp")})

Cooks_Model_Estimates_Temp <- data.frame(estimate = Cooks_Model_Temp$b, ci.lb = Cooks_Model_Temp$ci.lb, 
                                         ci.ub = Cooks_Model_Temp$ci.ub)
Cooks_Model_i2_Temp <- data.frame(round(orchaRd::i2_ml(Cooks_Model_Temp), 2))

#run <- FALSE
#system.time( # hasn't been ran yet
#  if(run){
#    Overall_Cooks_Temp <- cooks.distance(Cooks_Model_Temp)
#    saveRDS(Overall_Cooks_Temp, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks_Temp")
#  } else {
#    Overall_Cooks_Temp <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks_Temp")})

#dev.off()
#Cooks_Plot_Temp <- plot(Overall_Cooks_Temp, type = "o", pch = 21, xlab = "Observed Outcome", 
#                   ylab = "Cook's Distance", bg = "#183357")
#box(lwd = 2)

# Untransformed Model
system.time(  # 3ish minutes
  Untransformed_Temp <- brms::brm(Effect_Size_Type | se(sqrt(Variance_Type)) 
                                  ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                  data = Temp_Subset_Data,
                                  family = gaussian(),
                                  data2 = list(A = Temp_A), 
                                  chains = 4, 
                                  cores = 4,
                                  iter = 12000,
                                  warmup = 2000,
                                  thin = 5,
                                  prior = priors,
                                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                                  file = "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/untransformed_temp_model",
                                  file_refit = "on_change"))

b_untransformed_temp <- as_draws_df(Untransformed_Temp, variable = "b_Intercept")
b_untransformed_temp <- data.frame(b_untransformed_temp$b_Intercept)

sd_untransformed_temp <- as_draws_df(Untransformed_Temp, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_untransformed_temp <- data.frame("sd_obs__Intercept" = sd_untransformed_temp$sd_obs__Intercept, 
                                    "sd_phylo__Intercept" = sd_untransformed_temp$sd_phylo__Intercept, 
                                    "sd_Study_ID__Intercept" = sd_untransformed_temp$sd_Study_ID__Intercept)

mean_b_untransformed_temp <-  sapply(b_untransformed_temp, mean)
ci_b_untransformed_temp <- as.vector(HPDinterval(as.mcmc(b_untransformed_temp)))
pMCMC_b_untransformed_temp <- 2*(1 - max(table(b_untransformed_temp<0) / nrow(b_untransformed_temp)))

b_abs_untransformed_temp <- folded_norm(b_untransformed_temp[,1], rowSums(sd_untransformed_temp))
mean_abs_b_untransformed_temp <- mean(b_abs_untransformed_temp)
mode_abs_b_untransformed_temp <- posterior.mode(as.mcmc(b_abs_untransformed_temp))
ci.abs_untransformed_temp <- HPDinterval(as.mcmc(b_abs_untransformed_temp))

overall_i2_untransformed_temp <- i2(sd_untransformed_temp, Temp_Subset_Data$Variance_Type) 

##### Publication Bias - Salinity #####
# Importing Model
Model_Sal <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/sal_model.rds")

b_overall_sal <- as_draws_df(Model_Sal, variable = "b_Intercept")
b_overall_sal <- data.frame(b_overall_sal$b_Intercept)

sd_sal <- as_draws_df(Model_Sal, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_sal <- data.frame("sd_obs__Intercept" = sd_sal$sd_obs__Intercept, 
                      "sd_phylo__Intercept" = sd_sal$sd_phylo__Intercept, 
                      "sd_Study_ID__Intercept" = sd_sal$sd_Study_ID__Intercept)

b_abs_sal <- folded_norm(b_overall_sal[,1], rowSums(sd_sal))
mean_abs_b_sal <- mean(b_abs_sal)
mode_abs_b_sal <- posterior.mode(as.mcmc(b_abs_sal))
ci.abs_sal <- HPDinterval(as.mcmc(b_abs_sal))
overall_sal_i2 <- i2(sd_sal, Sal_Subset_Data$Variance_Type_Adjusted) 


# Residuals
Effect_Sizes_Sal <- Sal_Subset_Data %>% select(c("obs", "Effect_Size_Type_Adjusted"))

Predict_Model_Sal <- data.frame(predict(Model_Sal))
Predict_Model_Sal$obs <- 1:nrow(Predict_Model_Sal)
Residuals_Sal <- Predict_Model_Sal %>% 
                 left_join(Effect_Sizes_Sal, by = "obs")
Residuals_Sal <- Residuals_Sal %>% mutate(Residuals_Sal = `Effect_Size_Type_Adjusted` - `Estimate`)
Residuals_Sal <- Residuals_Sal[, c("Residuals_Sal", "Est.Error")]

# Funnel Plot
Funnel_Plot_Sal <- funnel(Residuals_Sal$Residuals_Sal, Residuals_Sal$Est.Error, yaxis = "seinv",
                          ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                          pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("top", legend = c("0.05 < p ≤ 1.00", "0 < p ≤ 0.05", "Studies"), 
       pch = c(22, 22, 21, 21), pt.bg = c("#FFFFFF","#D3DDEB", "#183357"), box.lwd = 2)

# Publication Year Graph
Graph_Data_Sal <- Sal_Subset_Data
Graph_Data_Sal <- Graph_Data_Sal %>% mutate(n_category = ifelse(n_.P1T1. <= 10, "10", 
                                                         ifelse(n_.P1T1. > 10 & n_.P1T1. <= 20, "20", 
                                                         ifelse(n_.P1T1. > 20 & n_.P1T1. <= 30, "> 20", "> 30"))))


Publication_Graph_Sal <- ggplot(Graph_Data_Sal, aes(x = Publication_Year, y = abs(Effect_Size_Type_Adjusted))) + 
                         geom_point(aes(x = Publication_Year, y = abs(Effect_Size_Type_Adjusted), 
                                        size = fct_relevel(n_category, c("10", "20", "> 20"))), 
                                        shape = 21, fill = "#4292c6", alpha = 0.5) + 
                         labs(x = "Publication Year", y = expression("Effect Size (PRRD"["S"]*")"), 
                              size = "Sample Size") +
                         theme_bw() +
                         theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                         theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                         theme(legend.position = "bottom", legend.direction = "horizontal") + 
                         geom_hline(yintercept = mode_abs_b_sal[1], lty = 2) + 
                         geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                         stat_poly_eq(formula = y ~ x, 
                                      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                      parse = TRUE) #+
                         #coord_cartesian(xlim = c(0, 25), 
                         #                ylim = c(-5, 5))

Publication_Graph_Sal #(750x500)

# Time-lag Bias
run <- FALSE
system.time( #  5ish minutes
  if(run){
    Year_Precision_Sal_rma <- metafor::rma.mv(abs(Effect_Size_Type_Adjusted), V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                              mods = ~ Publication_Year_Z + Precision - 1,
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                              R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(Year_Precision_Sal_rma, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_Sal_rma")
  } else {
    Year_Precision_Sal_rma <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Year_Precision_Sal_rma")})

Year_Precision_Sal_rma_Estimates <- data.frame(estimate = Year_Precision_Sal_rma$b, ci.lb = Year_Precision_Sal_rma$ci.lb, 
                                                ci.ub = Year_Precision_Sal_rma$ci.ub)
Year_Precision_Sal_rma_i2 <- data.frame(round(orchaRd::i2_ml(Year_Precision_Sal_rma), 2))

# Cooks Distance (Metafor)
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Cooks_Model_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                        R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                        control=list(rel.tol=1e-9))
    saveRDS(Cooks_Model_Sal, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model_sal")
  } else {
    Cooks_Model_Sal <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Cooks_model_sal")})

Cooks_Model_Estimates_Sal <- data.frame(estimate = Cooks_Model_Sal$b, ci.lb = Cooks_Model_Sal$ci.lb, 
                                         ci.ub = Cooks_Model_Sal$ci.ub)
Cooks_Model_i2_Sal <- data.frame(round(orchaRd::i2_ml(Cooks_Model_Sal), 2))

#run <- FALSE
#system.time( # hasn't been ran yet
#  if(run){
#    Overall_Cooks_Sal <- cooks.distance(Cooks_Model_Sal)
#    saveRDS(Overall_Cooks_Sal, "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks_Sal")
#  } else {
#    Overall_Cooks_Sal <- readRDS("./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/Overall_Cooks_Sal")})

#dev.off()
#Cooks_Plot_Sal <- plot(Overall_Cooks_Sal, type = "o", pch = 21, xlab = "Observed Outcome", 
#                   ylab = "Cook's Distance", bg = "#183357")
#box(lwd = 2)

# Untransformed Model
system.time(  # 3ish minutes
  Untransformed_Sal <- brms::brm(Effect_Size_Type | se(sqrt(Variance_Type)) 
                                 ~ 1 + (1|gr(phylo, cov = A)) + (1|Study_ID) + (1|obs),
                                 data = Sal_Subset_Data,
                                 family = gaussian(),
                                 data2 = list(A = Sal_A), 
                                 chains = 4, 
                                 cores = 4,
                                 iter = 12000,
                                 warmup = 2000,
                                 thin = 5,
                                 prior = priors,
                                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                                 file = "./4_Laboratory_Plasticity/3_Data_Analysis/2_Output/models/untransformed_sal_model",
                                 file_refit = "on_change"))

b_untransformed_sal <- as_draws_df(Untransformed_Sal, variable = "b_Intercept")
b_untransformed_sal <- data.frame(b_untransformed_sal$b_Intercept)

sd_untransformed_sal <- as_draws_df(Untransformed_Sal, variable = c("sd_obs__Intercept", "sd_phylo__Intercept", "sd_Study_ID__Intercept"))
sd_untransformed_sal <- data.frame("sd_obs__Intercept" = sd_untransformed_sal$sd_obs__Intercept, 
                                    "sd_phylo__Intercept" = sd_untransformed_sal$sd_phylo__Intercept, 
                                    "sd_Study_ID__Intercept" = sd_untransformed_sal$sd_Study_ID__Intercept)

mean_b_untransformed_sal <-  sapply(b_untransformed_sal, mean)
ci_b_untransformed_sal <- as.vector(HPDinterval(as.mcmc(b_untransformed_sal)))
pMCMC_b_untransformed_sal <- 2*(1 - max(table(b_untransformed_sal<0) / nrow(b_untransformed_sal)))

b_abs_untransformed_sal <- folded_norm(b_untransformed_sal[,1], rowSums(sd_untransformed_sal))
mean_abs_b_untransformed_sal <- mean(b_abs_untransformed_sal)
mode_abs_b_untransformed_sal <- posterior.mode(as.mcmc(b_abs_untransformed_sal))
ci.abs_untransformed_sal <- HPDinterval(as.mcmc(b_abs_untransformed_sal))

overall_i2_untransformed_sal <- i2(sd_untransformed_sal, Sal_Subset_Data$Variance_Type) 
