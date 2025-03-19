rm(list = ls())
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
if (!require("metaAidR")) devtools::install_github("daniel1noble/metaAidR", force = TRUE)
if (!require("orchaRd")) devtools::install_github("daniel1noble/orchaRd", force = TRUE)
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, robumeta)

##### Model Selection #####
# Importing Data Set
data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Complex_Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

# Phylogenetic Covariance Matrix
tree <- ape::read.tree("./3.Data_Analysis/2.Outputs/Phylogeny/Complex_tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

# Variance Matrix (Shared Control)
VCV <- make_VCV_matrix(data, V = "v_InRR", cluster = "Shared_Control_Number")

# All Possible Random Effects
run <- FALSE
system.time( #  1ish minutes
  if(run){
    All <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(All, "./3.Data_Analysis/2.Outputs/Models/Complex_All_Randoms.rds")
  } else {
    All <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_All_Randoms.rds")})

All_i2 <- data.frame(round(orchaRd::i2_ml(All), 2))
All_aic <- fitstats(All)

# No Measurement Random Effect
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Measurement <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                    ~1|Scientific_Name, ~1|Shared_Animal_Number), 
                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Measurement, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement.rds")
  } else {
    No_Measurement <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement.rds")})

No_Measurement_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement), 2))
No_Measurement_aic <- fitstats(No_Measurement)

# No Shared Animal Number Random Effect
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                               ~1|Scientific_Name, ~1|Measurement), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(No_Animal, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Animal.rds")
  } else {
    No_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Animal.rds")})

No_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal), 2))
No_Animal_aic <- fitstats(No_Animal)

# No Species Random Effect
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                  R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                  control=list(rel.tol=1e-9))
    saveRDS(No_Species, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Species.rds")
  } else {
    No_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Species.rds")})

No_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Species), 2))
No_Species_aic <- fitstats(No_Species)

# No Measurement or Shared Animal Number Random Effects
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Measurement_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name), 
                                             R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Animal, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement_Animal.rds")
  } else {
    No_Measurement_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement_Animal.rds")})

No_Measurement_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Animal), 2))
No_Measurement_Animal_aic <- fitstats(No_Measurement_Animal)

# No Measurement or Species Random Effects
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Measurement_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Animal_Number), 
                                              R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Species, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement_Species.rds")
  } else {
    No_Measurement_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Measurement_Species.rds")})

No_Measurement_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Species), 2))
No_Measurement_Species_aic <- fitstats(No_Measurement_Species)

# No Shared Animal or Species Random Effects
run <- FALSE
system.time( #  1ish minutes
  if(run){
    No_Animal_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Species, "./3.Data_Analysis/2.Outputs/Models/Complex_No_Animal_Species.rds")
  } else {
    No_Animal_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_No_Animal_Species.rds")})

No_Animal_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Species), 2))
No_Animal_Species_aic <- fitstats(No_Animal_Species)

# The Base Model
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Base <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                            random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                            R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                            control=list(rel.tol=1e-9))
    saveRDS(Base, "./3.Data_Analysis/2.Outputs/Models/Complex_Base.rds")
  } else {
    Base <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Base.rds")})

Base_i2 <- data.frame(round(orchaRd::i2_ml(Base), 2))
Base_aic <- fitstats(Base)

# AICc Summary
AICc <- data.frame("Models" = c("All", "No_Measurement", "No_Animal", "No_Species", 
                                "No_Measurement_Animal", "No_Measurement_Species", 
                                "No_Animal_Species", "Base"), 
                   "AICc" = c(All_aic[5], No_Measurement_aic[5], No_Animal_aic[5], No_Species_aic[5], 
                              No_Measurement_Animal_aic[5], No_Measurement_Species_aic[5], 
                              No_Animal_Species_aic[5], Base_aic[5]))

# Heterogeneity Summary
Heterogeneity <- data.frame("Random Effect" = c("Phylogeny", "Scientific_Name", "Shared_Animal_Number", "Study_ID", "Measurement", "Obs", "Total"), 
                            "All" = c(All_i2["I2_phylo", 1], All_i2["I2_Scientific_Name", 1], All_i2["I2_Shared_Animal_Number", 1],  
                                      All_i2["I2_Study_ID", 1], All_i2["I2_Measurement", 1], All_i2["I2_obs", 1], All_i2["I2_Total", 1]), 
                            "No_Measurement" = c(No_Measurement_i2["I2_phylo", 1], No_Measurement_i2["I2_Scientific_Name", 1], No_Measurement_i2["I2_Shared_Animal_Number", 1], 
                                                 No_Measurement_i2["I2_Study_ID", 1], NA, 
                                                 No_Measurement_i2["I2_obs", 1], No_Measurement_i2["I2_Total", 1]),
                            "No_Animal" = c(No_Animal_i2["I2_phylo", 1], No_Animal_i2["I2_Scientific_Name", 1], NA, 
                                            No_Animal_i2["I2_Study_ID", 1], No_Animal_i2["I2_Measurement", 1], 
                                            No_Animal_i2["I2_obs", 1], No_Animal_i2["I2_Total", 1]),
                            "No_Species" = c(No_Species_i2["I2_phylo", 1], NA, No_Species_i2["I2_Shared_Animal_Number", 1], 
                                             No_Species_i2["I2_Study_ID", 1], No_Species_i2["I2_Measurement", 1], 
                                             No_Species_i2["I2_obs", 1], No_Species_i2["I2_Total", 1]),
                            "No_Measurement_Animal" = c(No_Measurement_Animal_i2["I2_phylo", 1], No_Measurement_Animal_i2["I2_Scientific_Name", 1], 
                                                        NA, No_Measurement_Animal_i2["I2_Study_ID", 1], NA, 
                                                        No_Measurement_Animal_i2["I2_obs", 1], No_Measurement_Animal_i2["I2_Total", 1]), 
                            "No_Measurement_Species" = c(No_Measurement_Species_i2["I2_phylo", 1], NA, 
                                                         No_Measurement_Species_i2["I2_Shared_Animal_Number", 1],  
                                                         No_Measurement_Species_i2["I2_Study_ID", 1], NA, 
                                                         No_Measurement_Species_i2["I2_obs", 1], No_Measurement_Species_i2["I2_Total", 1]), 
                            "No_Animal_Species" = c(No_Animal_Species_i2["I2_phylo", 1], NA, 
                                                    NA, No_Animal_Species_i2["I2_Study_ID", 1], No_Animal_Species_i2["I2_Measurement", 1], 
                                                    No_Animal_Species_i2["I2_obs", 1], No_Animal_Species_i2["I2_Total", 1]),
                            "Base" = c(Base_i2["I2_phylo", 1], NA, NA,  
                                       Base_i2["I2_Study_ID", 1], NA, Base_i2["I2_obs", 1], Base_i2["I2_Total", 1]))
