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

# Importing Data Set
data <- read.csv("./Complex_Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

# Phylogenetic covariance matrix
tree <- ape::read.tree("./Complex_tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

#### Calculate Effect sizes #####

data <- data  %>% 
  mutate(PRRD = PRRD(t1 = T1_constant, t2 = T2_constant, 
                     t1_c = Mean_T1_C_Add, t2_c = Mean_T2_C_Add, t1_f = Mean_T1_F_Add, t2_f = Mean_T2_F_Add, 
                     sd_t1_c = SD_Final_T1_C_Add, sd_t2_c= SD_Final_T2_C_Add, sd_t1_f = SD_Final_T1_F_Add, sd_t2_f = SD_Final_T2_F_Add, 
                     n_t1_c =  n_T1_C, n_t2_c = n_T2_C, n_t1_f = n_T1_F, n_t2_f = n_T2_F, type = 'ef'),
         v_PRRD = PRRD(t1 = T1_constant, t2 = T2_constant, 
                       t1_c = Mean_T1_C_Add, t2_c = Mean_T2_C_Add, t1_f = Mean_T1_F_Add, t2_f = Mean_T2_F_Add, 
                       sd_t1_c = SD_Final_T1_C_Add, sd_t2_c= SD_Final_T2_C_Add, sd_t1_f = SD_Final_T1_F_Add, sd_t2_f = SD_Final_T2_F_Add, 
                       n_t1_c =  n_T1_C, n_t2_c = n_T2_C, n_t1_f = n_T1_F, n_t2_f = n_T2_F, type = 'v'))

# Check data with extreme precision. Mix up of SD/SE?
  checks <- data %>% filter(1/sqrt(v_PRRD) > 100)
  write.csv(checks, file = "checks.csv")

# Variance Matrix (Shared Control)
VCV <- make_VCV_matrix(data, V = "v_PRRD", cluster = "Shared_Control_Number")

#VCV_Untransformed <- make_VCV_matrix(data, V = "v_InRR_Untransformed", cluster = "Shared_Control_Number")

##### Publication Bias #####

# Model of the Residuals from the Overall Model. 
Model <- readRDS("./output/models/Complex_Overall_Model.rds")
Residuals <- rstandard(Model)
Residuals[c("slab", "digits")] = NULL
Residuals_df <- do.call("rbind", Residuals)
Residuals_df2 <- t(Residuals_df)
colnames(Residuals_df2) <- c("yi", "vi", "z")

Residuals_Model <- rma(yi = yi, vi = vi, data = Residuals_df2)

# Publication Year Graph
Graph_Data <- data
Graph_Data <- Graph_Data %>% mutate(n_category = ifelse(n_T1_C <= 25, "25", 
                                                 ifelse(n_T1_C > 25 & n_T1_C <= 50, "50", 
                                                 ifelse(n_T1_C > 50 & n_T1_C <= 75, "75", "> 75")))) %>% filter(!1/sqrt(v_PRRD) > 100)


Publication_Graph <- ggplot(Graph_Data, aes(x = Year, y = PRRD)) + 
                     geom_point(aes(x = Year, y = PRRD, 
                                    size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                shape = 21, fill = "#4292c6", alpha = 0.5) + 
                     labs(x = "Publication Year", y = expression("Effect Size (PRRD"["S"]*")"), 
                          size = "Sample Size") +
                     theme_bw() +
                     theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                     theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                     theme(legend.position = "bottom", legend.direction = "horizontal") + 
                     geom_hline(yintercept = Model$b, lty = 2) + 
                     geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                     stat_poly_eq(formula = y ~ x, 
                     aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                         parse = TRUE) #+
                     #coord_cartesian(xlim = c(0, 25), 
                     #                ylim = c(-5, 5))
Publication_Graph

# Funnel Plot
par(mfrow = c(1,2))
Funnel_Plot <- funnel(Residuals_Model, yaxis = "seinv",
                      ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                      pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

# Time-lag Bias
run <- TRUE

# Variance Matrix (Shared Control)
VCV <- make_VCV_matrix(Graph_Data, V = "v_PRRD", cluster = "Shared_Control_Number")

Graph_Data$Precision <- 1/sqrt(Graph_Data$v_PRRD) # 1/SE

system.time(
  if(run){
    TL_Model <- metafor::rma.mv(PRRD, V = VCV, test = "t", dfs = "contain",
                                mods = ~ Year_Z + Precision,
                                random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                              ~1|Shared_Animal_Number, ~1|Measurement), 
                                R = list(phylo=A_cor), data = Graph_Data, method = "REML", sparse = TRUE, 
                                control=list(rel.tol=1e-9))
    saveRDS(TL_Model, "./Complex_TL_Model.rds")
  } else {
    TL_Model <- readRDS("./Complex_TL_Model.rds")})

##### Sensitivity Analysis #####

run <- TRUE
system.time(
  if(run){
    Cooks_Overall_Model <- metafor::rma.mv(PRRD ~ 1, V = v_PRRD, test = "t", dfs = "contain",
                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                     R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                     control=list(rel.tol=1e-9))
    saveRDS(Cooks_Overall_Model, "./Cooks_Complex_Overall_Model.rds")
  } else {
    Cooks_Overall_Model <- readRDS("./Cooks_Complex_Overall_Model.rds")})

# Cooks Distance
run <- TRUE
system.time( #  9.5 Days
  if(run){
    Overall_Cooks <- cooks.distance(Cooks_Overall_Model)
    saveRDS(Overall_Cooks, "./Complex_Overall_Cooks.rds")
  } else {
    Overall_Cooks <- readRDS("./Complex_Overall_Cooks.rds")})

dev.off()
Cooks_Plot <- plot(Overall_Cooks, type = "o", pch = 21, xlab = "Observed Outcome", 
                   ylab = "Cook's Distance", bg = "#183357")
box(lwd = 2)

# Untransformed
Model_Estimates <- data.frame(estimate = Model$b, ci.lb = Model$ci.lb, ci.ub = Model$ci.ub)
Model_i2 <- data.frame(round(orchaRd::i2_ml(Model), 2))

run <- TRUE
system.time(
  if(run){
    Untransformed_Model <- metafor::rma.mv(InRR_Untransformed ~ 1, V = VCV_Untransformed, test = "t", dfs = "contain",
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Untransformed_Model, "./Complex_Untransformed_Model.rds")
  } else {
    Untransformed_Model <- readRDS("./Complex_Untransformed_Model.rds")})

Untransformed_Model_Estimates <- data.frame(estimate = Untransformed_Model$b, ci.lb = Untransformed_Model$ci.lb, ci.ub = Untransformed_Model$ci.ub)
Untransformed_Model_i2 <- data.frame(round(orchaRd::i2_ml(Untransformed_Model), 2))
