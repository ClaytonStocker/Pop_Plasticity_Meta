rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, robumeta, cli)

##### Model Selection #####
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

# All Possible Random Effects
run <- TRUE
system.time(
if(run){
  All <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                       ~1|Scientific_Name, ~1|Shared_Animal_Number, ~1|Measurement), 
                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                         control=list(rel.tol=1e-9))
  saveRDS(All, "./All.rds")
} else {
  All <- readRDS("./All.rds")})

All_i2 <- data.frame(round(orchaRd::i2_ml(All), 2))
All_aic <- fitstats(All)

# No Scientific Name Random Effect
run <- TRUE
system.time(
  if(run){
    No_Species <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                  R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                  control=list(rel.tol=1e-9))
    saveRDS(No_Species, "./No_Species.rds")
  } else {
    No_Species <- readRDS("./No_Species.rds")})

No_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Species), 2))
No_Species_aic <- fitstats(No_Species)

# No Shared Animal Number Random Effect
run <- TRUE
system.time(
  if(run){
    No_Animal <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                               ~1|Scientific_Name, ~1|Measurement), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(No_Animal, "./No_Animal.rds")
  } else {
    No_Animal <- readRDS("./No_Animal.rds")})

No_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal), 2))
No_Animal_aic <- fitstats(No_Animal)

# No Measurement Random Effect
run <- TRUE
system.time(
  if(run){
    No_Measurement <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                    ~1|Scientific_Name, ~1|Shared_Animal_Number), 
                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Measurement, "./No_Measurement.rds")
  } else {
    No_Measurement <- readRDS("./No_Measurement.rds")})

No_Measurement_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement), 2))
No_Measurement_aic <- fitstats(No_Measurement)

# No Scientific Name or Shared Animal Number Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Animal <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                       ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(No_Species_Animal, "./No_Species_Animal.rds")
  } else {
    No_Species_Animal <- readRDS("./No_Species_Animal.rds")})

No_Species_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Animal), 2))
No_Species_Animal_aic <- fitstats(No_Species_Animal)

# No Scientific Name or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Measurement <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                            ~1|Shared_Animal_Number), 
                                              R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(No_Species_Measurement, "./No_Species_Measurement.rds")
  } else {
    No_Species_Measurement <- readRDS("./No_Species_Measurement.rds")})

No_Species_Measurement_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Measurement), 2))
No_Species_Measurement_aic <- fitstats(No_Species_Measurement)

# No Shared Animal Number or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Animal_Measurement <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                           ~1|Scientific_Name), 
                                             R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Measurement, "./No_Animal_Measurement.rds")
  } else {
    No_Animal_Measurement <- readRDS("./No_Animal_Measurement.rds")})

No_Animal_Measurement_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Measurement), 2))
No_Animal_Measurement_aic <- fitstats(No_Animal_Measurement)

# The Base Model
run <- TRUE
system.time(
  if(run){
    Base <- metafor::rma.mv(Effect_Size_Adjusted ~ 1, V = Variance_Adjusted, test = "t", dfs = "contain",
                            random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                            R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                            control=list(rel.tol=1e-9))
    saveRDS(Base, "./Base.rds")
  } else {
    Base <- readRDS("./Base.rds")})

Base_i2 <- data.frame(round(orchaRd::i2_ml(Base), 2))
Base_aic <- fitstats(Base)

# AICc Summary
AICc <- data.frame("Models" = c("All", "No_Species", "No_Animal", "No_Measurement", "No_Species_Animal", 
                                "No_Species_Measurement", "No_Animal_Measurement", "Base"), 
                   "AICc" = c(All_aic[5], No_Species_aic[5], No_Animal_aic[5], No_Measurement_aic[5], 
                              No_Species_Animal_aic[5], No_Species_Measurement_aic[5], 
                              No_Animal_Measurement_aic[5], Base_aic[5]))

# Heterogeneity Summary
Heterogeneity <- data.frame("Random Effects" = c("Animal", "Measurement", "Obs", "Phylo", "Species", "Study", "Total"), 
                            "All" = c(All_i2["I2_Shared_Animal_Number", 1], All_i2["I2_Measurement", 1], 
                                      All_i2["I2_obs", 1], All_i2["I2_phylo", 1], All_i2["I2_Scientific_Name", 1], 
                                      All_i2["I2_Study_ID", 1], All_i2["I2_Total", 1]), 
                            "No_Species" = c(No_Species_i2["I2_Shared_Animal_Number", 1], No_Species_i2["I2_Measurement", 1], 
                                             No_Species_i2["I2_obs", 1], No_Species_i2["I2_phylo", 1], NA, 
                                             No_Species_i2["I2_Study_ID", 1], No_Species_i2["I2_Total", 1]), 
                            "No_Animal" = c(NA, No_Animal_i2["I2_Measurement", 1], 
                                            No_Animal_i2["I2_obs", 1], No_Animal_i2["I2_phylo", 1], No_Animal_i2["I2_Scientific_Name", 1], 
                                            No_Animal_i2["I2_Study_ID", 1], No_Animal_i2["I2_Total", 1]), 
                            "No_Measurement" = c(No_Measurement_i2["I2_Shared_Animal_Number", 1], NA, 
                                                 No_Measurement_i2["I2_obs", 1], No_Measurement_i2["I2_phylo", 1], No_Measurement_i2["I2_Scientific_Name", 1], 
                                                 No_Measurement_i2["I2_Study_ID", 1], No_Measurement_i2["I2_Total", 1]), 
                            "No_Species_Animal" = c(NA, No_Species_Animal_i2["I2_Measurement", 1], 
                                                    No_Species_Animal_i2["I2_obs", 1], No_Species_Animal_i2["I2_phylo", 1], NA, 
                                                    No_Species_Animal_i2["I2_Study_ID", 1], No_Species_Animal_i2["I2_Total", 1]), 
                            "No_Species_Measurement" = c(No_Species_Measurement_i2["I2_Shared_Animal_Number", 1], NA, 
                                                         No_Species_Measurement_i2["I2_obs", 1], No_Species_Measurement_i2["I2_phylo", 1], NA, 
                                                         No_Species_Measurement_i2["I2_Study_ID", 1], No_Species_Measurement_i2["I2_Total", 1]), 
                            "No_Animal_Measurement" = c(NA, NA, No_Animal_Measurement_i2["I2_obs", 1], 
                                                        No_Animal_Measurement_i2["I2_phylo", 1], No_Animal_Measurement_i2["I2_Scientific_Name", 1], 
                                                        No_Animal_Measurement_i2["I2_Study_ID", 1], No_Animal_Measurement_i2["I2_Total", 1]),
                            "Base" = c(NA, NA, 
                                       Base_i2["I2_obs", 1], Base_i2["I2_phylo", 1], NA, 
                                       Base_i2["I2_Study_ID", 1], Base_i2["I2_Total", 1]))

##### Model Selection - Temperature Data #####
# Phylogenetic covariance matrix
Temp_Subset_Data <- data %>% filter(!is.na(Effect_Size_Type_Adjusted))
Temp_Subset_Data <- Temp_Subset_Data %>% filter(`Type` == "Temperature")
Temp_Species <- Temp_Subset_Data %>% select("phylo") %>% unique()

Temp_A_cor <- as.data.frame(A_cor)
Temp_A_cor <- Temp_A_cor[c(Temp_Species$phylo), c(Temp_Species$phylo)]
Temp_A_cor <- as.matrix(Temp_A_cor)

# All Possible Random Effects
run <- TRUE
system.time(
  if(run){
    All_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                              ~1|Scientific_Name, ~1|Shared_Animal_Number, ~1|Measurement), 
                                R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                control=list(rel.tol=1e-9))
    saveRDS(All_Temp, "./All_Temp.rds")
  } else {
    All_Temp <- readRDS("./All_Temp.rds")})

All_Temp_i2 <- data.frame(round(orchaRd::i2_ml(All_Temp), 2))
All_Temp_aic <- fitstats(All_Temp)

# No Scientific Name Random Effect
run <- TRUE
system.time(
  if(run){
    No_Species_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                       R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                       control=list(rel.tol=1e-9))
    saveRDS(No_Species_Temp, "./No_Species_Temp.rds")
  } else {
    No_Species_Temp <- readRDS("./No_Species_Temp.rds")})

No_Species_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Temp), 2))
No_Species_Temp_aic <- fitstats(No_Species_Temp)

# No Shared Animal Number Random Effect
run <- TRUE
system.time(
  if(run){
    No_Animal_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                    ~1|Scientific_Name, ~1|Measurement), 
                                      R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Temp, "./No_Animal_Temp.rds")
  } else {
    No_Animal_Temp <- readRDS("./No_Animal_Temp.rds")})

No_Animal_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Temp), 2))
No_Animal_Temp_aic <- fitstats(No_Animal_Temp)

# No Measurement Random Effect
run <- TRUE
system.time(
  if(run){
    No_Measurement_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                         ~1|Scientific_Name, ~1|Shared_Animal_Number), 
                                           R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Temp, "./No_Measurement_Temp.rds")
  } else {
    No_Measurement_Temp <- readRDS("./No_Measurement_Temp.rds")})

No_Measurement_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Temp), 2))
No_Measurement_Temp_aic <- fitstats(No_Measurement_Temp)

# No Scientific Name or Shared Animal Number Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Animal_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                            ~1|Measurement), 
                                              R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(No_Species_Animal_Temp, "./No_Species_Animal_Temp.rds")
  } else {
    No_Species_Animal_Temp <- readRDS("./No_Species_Animal_Temp.rds")})

No_Species_Animal_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Animal_Temp), 2))
No_Species_Animal_Temp_aic <- fitstats(No_Species_Animal_Temp)

# No Scientific Name or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Measurement_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                                 ~1|Shared_Animal_Number), 
                                                   R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(No_Species_Measurement_Temp, "./No_Species_Measurement_Temp.rds")
  } else {
    No_Species_Measurement_Temp <- readRDS("./No_Species_Measurement_Temp.rds")})

No_Species_Measurement_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Measurement_Temp), 2))
No_Species_Measurement_Temp_aic <- fitstats(No_Species_Measurement_Temp)

# No Shared Animal Number or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Animal_Measurement_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                                ~1|Scientific_Name), 
                                                  R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Measurement_Temp, "./No_Animal_Measurement_Temp.rds")
  } else {
    No_Animal_Measurement_Temp <- readRDS("./No_Animal_Measurement_Temp.rds")})

No_Animal_Measurement_Temp_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Measurement_Temp), 2))
No_Animal_Measurement_Temp_aic <- fitstats(No_Animal_Measurement_Temp)

# The Base Model
run <- TRUE
system.time(
  if(run){
    Base_Temp <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                 R = list(phylo=Temp_A_cor), data = Temp_Subset_Data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(Base_Temp, "./Base_Temp.rds")
  } else {
    Base_Temp <- readRDS("./Base_Temp.rds")})

Base_Temp_i2 <- data.frame(round(orchaRd::i2_ml(Base_Temp), 2))
Base_Temp_aic <- fitstats(Base_Temp)

# AICc Summary
AICc_Temp <- data.frame("Models" = c("All", "No_Species", "No_Animal", "No_Measurement", "No_Species_Animal", 
                                     "No_Species_Measurement", "No_Animal_Measurement", "Base"), 
                        "AICc" = c(All_Temp_aic[5], No_Species_Temp_aic[5], No_Animal_Temp_aic[5], No_Measurement_Temp_aic[5], 
                                   No_Species_Animal_Temp_aic[5], No_Species_Measurement_Temp_aic[5], 
                                   No_Animal_Measurement_Temp_aic[5], Base_Temp_aic[5]))

# Heterogeneity Summary
Heterogeneity_Temp <- data.frame("Random Effects" = c("Animal", "Measurement", "Obs", "Phylo", "Species", "Study", "Total"), 
                                 "All" = c(All_Temp_i2["I2_Shared_Animal_Number", 1], All_Temp_i2["I2_Measurement", 1], 
                                           All_Temp_i2["I2_obs", 1], All_Temp_i2["I2_phylo", 1], All_Temp_i2["I2_Scientific_Name", 1], 
                                           All_Temp_i2["I2_Study_ID", 1], All_Temp_i2["I2_Total", 1]), 
                                 "No_Species" = c(No_Species_Temp_i2["I2_Shared_Animal_Number", 1], No_Species_Temp_i2["I2_Measurement", 1], 
                                                  No_Species_Temp_i2["I2_obs", 1], No_Species_Temp_i2["I2_phylo", 1], NA, 
                                                  No_Species_Temp_i2["I2_Study_ID", 1], No_Species_Temp_i2["I2_Total", 1]), 
                                 "No_Animal" = c(NA, No_Animal_Temp_i2["I2_Measurement", 1], 
                                                 No_Animal_Temp_i2["I2_obs", 1], No_Animal_Temp_i2["I2_phylo", 1], No_Animal_Temp_i2["I2_Scientific_Name", 1], 
                                                 No_Animal_Temp_i2["I2_Study_ID", 1], No_Animal_Temp_i2["I2_Total", 1]), 
                                 "No_Measurement" = c(No_Measurement_Temp_i2["I2_Shared_Animal_Number", 1], NA, 
                                                      No_Measurement_Temp_i2["I2_obs", 1], No_Measurement_Temp_i2["I2_phylo", 1], No_Measurement_Temp_i2["I2_Scientific_Name", 1], 
                                                      No_Measurement_Temp_i2["I2_Study_ID", 1], No_Measurement_Temp_i2["I2_Total", 1]), 
                                 "No_Species_Animal" = c(NA, No_Species_Animal_Temp_i2["I2_Measurement", 1], 
                                                         No_Species_Animal_Temp_i2["I2_obs", 1], No_Species_Animal_Temp_i2["I2_phylo", 1], NA, 
                                                         No_Species_Animal_Temp_i2["I2_Study_ID", 1], No_Species_Animal_Temp_i2["I2_Total", 1]), 
                                 "No_Species_Measurement" = c(No_Species_Measurement_Temp_i2["I2_Shared_Animal_Number", 1], NA, 
                                                              No_Species_Measurement_Temp_i2["I2_obs", 1], No_Species_Measurement_Temp_i2["I2_phylo", 1], NA, 
                                                              No_Species_Measurement_Temp_i2["I2_Study_ID", 1], No_Species_Measurement_Temp_i2["I2_Total", 1]), 
                                 "No_Animal_Measurement" = c(NA, NA, No_Animal_Measurement_Temp_i2["I2_obs", 1], 
                                                             No_Animal_Measurement_Temp_i2["I2_phylo", 1], No_Animal_Measurement_Temp_i2["I2_Scientific_Name", 1], 
                                                             No_Animal_Measurement_Temp_i2["I2_Study_ID", 1], No_Animal_Measurement_Temp_i2["I2_Total", 1]),
                                 "Base" = c(NA, NA, 
                                            Base_Temp_i2["I2_obs", 1], Base_Temp_i2["I2_phylo", 1], NA, 
                                            Base_Temp_i2["I2_Study_ID", 1], Base_Temp_i2["I2_Total", 1]))

##### Model Selection - Salinity Data #####
# Phylogenetic covariance matrix
Sal_Subset_Data <- data %>% filter(!is.na(Effect_Size_Type_Adjusted))
Sal_Subset_Data <- Sal_Subset_Data %>% filter(`Type` == "Salinity")
Sal_Species <- Sal_Subset_Data %>% select("phylo") %>% unique()

Sal_A_cor <- as.data.frame(A_cor)
Sal_A_cor <- Sal_A_cor[c(Sal_Species$phylo), c(Sal_Species$phylo)]
Sal_A_cor <- as.matrix(Sal_A_cor)

# All Possible Random Effects
run <- TRUE
system.time(
  if(run){
    All_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                             ~1|Scientific_Name, ~1|Shared_Animal_Number, ~1|Measurement), 
                               R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                               control=list(rel.tol=1e-9))
    saveRDS(All_Sal, "./All_Sal.rds")
  } else {
    All_Sal <- readRDS("./All_Sal.rds")})

All_Sal_i2 <- data.frame(round(orchaRd::i2_ml(All_Sal), 2))
All_Sal_aic <- fitstats(All_Sal)

# No Scientific Name Random Effect
run <- TRUE
system.time(
  if(run){
    No_Species_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                    ~1|Shared_Animal_Number, ~1|Measurement), 
                                      R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Species_Sal, "./No_Species_Sal.rds")
  } else {
    No_Species_Sal <- readRDS("./No_Species_Sal.rds")})

No_Species_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Sal), 2))
No_Species_Sal_aic <- fitstats(No_Species_Sal)

# No Shared Animal Number Random Effect
run <- TRUE
system.time(
  if(run){
    No_Animal_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                    ~1|Scientific_Name, ~1|Measurement), 
                                      R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Sal, "./No_Animal_Sal.rds")
  } else {
    No_Animal_Sal <- readRDS("./No_Animal_Sal.rds")})

No_Animal_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Sal), 2))
No_Animal_Sal_aic <- fitstats(No_Animal_Sal)

# No Measurement Random Effect
run <- TRUE
system.time(
  if(run){
    No_Measurement_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                          random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                        ~1|Scientific_Name, ~1|Shared_Animal_Number), 
                                          R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                          control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Sal, "./No_Measurement_Sal.rds")
  } else {
    No_Measurement_Sal <- readRDS("./No_Measurement_Sal.rds")})

No_Measurement_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Sal), 2))
No_Measurement_Sal_aic <- fitstats(No_Measurement_Sal)

# No Scientific Name or Shared Animal Number Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Animal_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                           ~1|Measurement), 
                                             R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(No_Species_Animal_Sal, "./No_Species_Animal_Sal.rds")
  } else {
    No_Species_Animal_Sal <- readRDS("./No_Species_Animal_Sal.rds")})

No_Species_Animal_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Animal_Sal), 2))
No_Species_Animal_Sal_aic <- fitstats(No_Species_Animal_Sal)

# No Scientific Name or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Species_Measurement_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                                ~1|Shared_Animal_Number), 
                                                  R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(No_Species_Measurement_Sal, "./No_Species_Measurement_Sal.rds")
  } else {
    No_Species_Measurement_Sal <- readRDS("./No_Species_Measurement_Sal.rds")})

No_Species_Measurement_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_Measurement_Sal), 2))
No_Species_Measurement_Sal_aic <- fitstats(No_Species_Measurement_Sal)

# No Shared Animal Number or Measurement Random Effects
run <- TRUE
system.time(
  if(run){
    No_Animal_Measurement_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                               ~1|Scientific_Name), 
                                                 R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Measurement_Sal, "./No_Animal_Measurement_Sal.rds")
  } else {
    No_Animal_Measurement_Sal <- readRDS("./No_Animal_Measurement_Sal.rds")})

No_Animal_Measurement_Sal_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Measurement_Sal), 2))
No_Animal_Measurement_Sal_aic <- fitstats(No_Animal_Measurement_Sal)

# The Base Model
run <- TRUE
system.time(
  if(run){
    Base_Sal <- metafor::rma.mv(Effect_Size_Type_Adjusted ~ 1, V = Variance_Type_Adjusted, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                                 R = list(phylo=Sal_A_cor), data = Sal_Subset_Data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(Base_Sal, "./Base_Sal.rds")
  } else {
    Base_Sal <- readRDS("./Base_Sal.rds")})

Base_Sal_i2 <- data.frame(round(orchaRd::i2_ml(Base_Sal), 2))
Base_Sal_aic <- fitstats(Base_Sal)

# AICc Summary
AICc_Sal <- data.frame("Models" = c("All", "No_Species", "No_Animal", "No_Measurement", "No_Species_Animal", 
                                     "No_Species_Measurement", "No_Animal_Measurement", "Base"), 
                        "AICc" = c(All_Sal_aic[5], No_Species_Sal_aic[5], No_Animal_Sal_aic[5], No_Measurement_Sal_aic[5], 
                                   No_Species_Animal_Sal_aic[5], No_Species_Measurement_Sal_aic[5], 
                                   No_Animal_Measurement_Sal_aic[5], Base_Sal_aic[5]))

# Heterogeneity Summary
Heterogeneity_Sal <- data.frame("Random Effects" = c("Animal", "Measurement", "Obs", "Phylo", "Species", "Study", "Total"), 
                                 "All" = c(All_Sal_i2["I2_Shared_Animal_Number", 1], All_Sal_i2["I2_Measurement", 1], 
                                           All_Sal_i2["I2_obs", 1], All_Sal_i2["I2_phylo", 1], All_Sal_i2["I2_Scientific_Name", 1], 
                                           All_Sal_i2["I2_Study_ID", 1], All_Sal_i2["I2_Total", 1]), 
                                 "No_Species" = c(No_Species_Sal_i2["I2_Shared_Animal_Number", 1], No_Species_Sal_i2["I2_Measurement", 1], 
                                                  No_Species_Sal_i2["I2_obs", 1], No_Species_Sal_i2["I2_phylo", 1], NA, 
                                                  No_Species_Sal_i2["I2_Study_ID", 1], No_Species_Sal_i2["I2_Total", 1]), 
                                 "No_Animal" = c(NA, No_Animal_Sal_i2["I2_Measurement", 1], 
                                                 No_Animal_Sal_i2["I2_obs", 1], No_Animal_Sal_i2["I2_phylo", 1], No_Animal_Sal_i2["I2_Scientific_Name", 1], 
                                                 No_Animal_Sal_i2["I2_Study_ID", 1], No_Animal_Sal_i2["I2_Total", 1]), 
                                 "No_Measurement" = c(No_Measurement_Sal_i2["I2_Shared_Animal_Number", 1], NA, 
                                                      No_Measurement_Sal_i2["I2_obs", 1], No_Measurement_Sal_i2["I2_phylo", 1], No_Measurement_Sal_i2["I2_Scientific_Name", 1], 
                                                      No_Measurement_Sal_i2["I2_Study_ID", 1], No_Measurement_Sal_i2["I2_Total", 1]), 
                                 "No_Species_Animal" = c(NA, No_Species_Animal_Sal_i2["I2_Measurement", 1], 
                                                         No_Species_Animal_Sal_i2["I2_obs", 1], No_Species_Animal_Sal_i2["I2_phylo", 1], NA, 
                                                         No_Species_Animal_Sal_i2["I2_Study_ID", 1], No_Species_Animal_Sal_i2["I2_Total", 1]), 
                                 "No_Species_Measurement" = c(No_Species_Measurement_Sal_i2["I2_Shared_Animal_Number", 1], NA, 
                                                              No_Species_Measurement_Sal_i2["I2_obs", 1], No_Species_Measurement_Sal_i2["I2_phylo", 1], NA, 
                                                              No_Species_Measurement_Sal_i2["I2_Study_ID", 1], No_Species_Measurement_Sal_i2["I2_Total", 1]), 
                                 "No_Animal_Measurement" = c(NA, NA, No_Animal_Measurement_Sal_i2["I2_obs", 1], 
                                                             No_Animal_Measurement_Sal_i2["I2_phylo", 1], No_Animal_Measurement_Sal_i2["I2_Scientific_Name", 1], 
                                                             No_Animal_Measurement_Sal_i2["I2_Study_ID", 1], No_Animal_Measurement_Sal_i2["I2_Total", 1]),
                                 "Base" = c(NA, NA, 
                                            Base_Sal_i2["I2_obs", 1], Base_Sal_i2["I2_phylo", 1], NA, 
                                            Base_Sal_i2["I2_Study_ID", 1], Base_Sal_i2["I2_Total", 1]))

##### MCMC Diagnostic - Overall Data Set #####

# Overall BRMS Model
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

predict(overall)

# Bayesian Model/Data Output

# Extracting the posterior distributions
b_overall <- as_draws_df(overall, variable = "b_Intercept")

# Graphing autocorrelation of MCMC chains
auto_plots <- bayesplot::mcmc_acf(b_overall, lags = 10) +
              theme_bw() +
              theme(strip.text.y = element_blank()) +
              theme(strip.text.x.top = element_blank()) +
              theme(axis.text.y = element_text(size = 10, colour = "black")) +
              theme(axis.text.x = element_text(size = 10, colour = "black")) +
              theme(axis.title.y = element_text(size = 11, colour = "black", margin = margin(r = 10))) +
              theme(axis.title.x = element_text(size = 11, colour = "black", margin = margin(t = 5)))

auto_plots # 400x800

##### MCMC Diagnostics - Temperature Data

Temp_A <- as.data.frame(A)
Temp_A <- Temp_A[c(Temp_Species$phylo), c(Temp_Species$phylo)]
Temp_A <- as.matrix(Temp_A)

# Temperature BRMS Model
priors <-  prior(student_t(3, 0, 20), class = "sd")

system.time(
  temp <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
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
                    file = "./temp_model",
                    file_refit = "always"))

# Bayesian Model/Data Output

# Extracting the posterior distributions
b_temp <- as_draws_df(temp, variable = "b_Intercept")

# Graphing autocorrelation of MCMC chains
auto_plots_temp <- bayesplot::mcmc_acf(b_temp, lags = 10) +
                   theme_bw() +
                   theme(strip.text.y = element_blank()) +
                   theme(strip.text.x.top = element_blank()) +
                   theme(axis.text.y = element_text(size = 10, colour = "black")) +
                   theme(axis.text.x = element_text(size = 10, colour = "black")) +
                   theme(axis.title.y = element_text(size = 11, colour = "black", margin = margin(r = 10))) +
                   theme(axis.title.x = element_text(size = 11, colour = "black", margin = margin(t = 5)))

auto_plots_temp # 400x800

##### MCMC Diagnostics - Salinity Data

Sal_A <- as.data.frame(A)
Sal_A <- Sal_A[c(Sal_Species$phylo), c(Sal_Species$phylo)]
Sal_A <- as.matrix(Sal_A)

# Salinity BRMS Model
priors <-  prior(student_t(3, 0, 20), class = "sd")

system.time(
  sal <- brms::brm(Effect_Size_Type_Adjusted | se(sqrt(Variance_Type_Adjusted)) 
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
                   file = "./sal_model",
                   file_refit = "always"))

# Bayesian Model/Data Output

# Extracting the posterior distributions
b_sal <- as_draws_df(sal, variable = "b_Intercept")

# Graphing autocorrelation of MCMC chains
auto_plots_sal <- bayesplot::mcmc_acf(b_sal, lags = 10) +
                  theme_bw() +
                  theme(strip.text.y = element_blank()) +
                  theme(strip.text.x.top = element_blank()) +
                  theme(axis.text.y = element_text(size = 10, colour = "black")) +
                  theme(axis.text.x = element_text(size = 10, colour = "black")) +
                  theme(axis.title.y = element_text(size = 11, colour = "black", margin = margin(r = 10))) +
                  theme(axis.title.x = element_text(size = 11, colour = "black", margin = margin(t = 5)))

auto_plots_sal # 400x800
