rm(list = ls())
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
if (!require("metaAidR")) devtools::install_gitgub("daniel1noble/metaAidR", force = TRUE)
if (!require("orchaRd")) devtools::install_gitgub("daniel1noble/orchaRd", force = TRUE)
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
               robumeta, ggpmisc, ggridges, ggbeeswarm, gridExtra)

# Importing Data Set
data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Complex_Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

# Phylogenetic covariance matrix
tree <- ape::read.tree("./3.Data_Analysis/2.Outputs/Phylogeny/Complex_tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

# Variance Matrix (Shared Control)
VCV <- make_VCV_matrix(data, V = "v_InRR", cluster = "Shared_Control_Number")

##### Overall Model #####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Overall_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = VCV, test = "t", dfs = "contain",
                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                     R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                     control=list(rel.tol=1e-9))
    saveRDS(Overall_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Overall_Model.rds")
  } else {
    Overall_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Overall_Model.rds")})

Overall_Model_rob <- robust(Overall_Model, cluster = data$Study_ID, adjust = TRUE)

Overall_Model_Estimates <- data.frame(estimate = Overall_Model$b, ci.lb = Overall_Model$ci.lb, 
                                      ci.ub = Overall_Model$ci.ub)
Overall_Model_i2 <- data.frame(round(orchaRd::i2_ml(Overall_Model), 2))

#### Overall Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = VCV, test = "t", dfs = "contain",
                                       mods = ~ Fluctuation_Magnitude - 1,
                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                       R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                       control=list(rel.tol=1e-9))
    saveRDS(Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Amplitude_Model.rds")
  } else {
    Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Amplitude_Model.rds")})

Amplitude_Model_rob <- robust(Amplitude_Model, cluster = data$Study_ID, adjust = TRUE)

Amplitude_Model_Estimates <- data.frame(estimate = Amplitude_Model$b, ci.lb = Amplitude_Model$ci.lb, 
                                        ci.ub = Amplitude_Model$ci.ub)
Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Amplitude_Model), 2))

Plot_Data <- data
Plot_Data <- Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                               ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                               ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

Amplitude_Plot <- ggplot(Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                  geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                 size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                 shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                  labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                  size = "Sample Size", title = "Overall Analysis") +
                  theme_bw() +
                  theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                  theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                  theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                  #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                  geom_hline(yintercept = Overall_Model_Estimates$estimate, lty = 2) + 
                  geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                  stat_poly_eq(formula = y ~ x, 
                               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                               parse = TRUE) +
                  coord_cartesian(xlim = c(0, 25), 
                                  ylim = c(-0.25, 0.25))

Amplitude_Plot #(330x330)

#### Overall Model - Type of Fluctuation Meta-Regression ####
Fluctuation_Data <- data %>% filter(!is.na(Fluctuation_Category))

Fluctuation_Exploration <- Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Exploration) <- Fluctuation_Exploration$Fluctuation_Category

Fluctuation_Data <- Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")

Fluctuation_Species_Count <- Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% 
                             table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                             select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Species_Count) <- Fluctuation_Species_Count$Fluctuation_Category

Fluctuation_Study_Count <- Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% 
                           table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                           select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Study_Count) <- Fluctuation_Study_Count$Fluctuation_Category

Fluctuation_Species <- Fluctuation_Data %>% select("phylo") %>% unique()

Fluctuation_A_cor <- as.data.frame(A_cor)
Fluctuation_A_cor <- Fluctuation_A_cor[c(Fluctuation_Species$phylo), c(Fluctuation_Species$phylo)]
Fluctuation_A_cor <- as.matrix(Fluctuation_A_cor)

Fluctuation_VCV <- make_VCV_matrix(Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Fluctuation_VCV, test = "t", dfs = "contain",
                                         mods = ~ Fluctuation_Category - 1,
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                         R = list(phylo=Fluctuation_A_cor), data = Fluctuation_Data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Fluctuation_Model.rds")
  } else {
    Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Fluctuation_Model.rds")})

Fluctuation_Model_rob <- robust(Fluctuation_Model, cluster = Fluctuation_Data$Study_ID, adjust = TRUE)

Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Fluctuation_Model$b), 21, 100),
                                          estimate = Fluctuation_Model$b, ci.lb = Fluctuation_Model$ci.lb, 
                                          ci.ub = Fluctuation_Model$ci.ub)
rownames(Fluctuation_Model_Estimates) <- Fluctuation_Model_Estimates$Category
Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Fluctuation_Model), 2))

# Preparing Graph - Combined

fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

fluctuation_k <- data.frame("k" = c(Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                    Fluctuation_Exploration["Alternating", "Freq"], 
                                    Fluctuation_Exploration["Stepwise", "Freq"]), 
                            row.names = fluctuation_rnames)

fluctuation_group_no <- data.frame("Spp No." = c(Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                 Fluctuation_Species_Count["Alternating", "Freq"], 
                                                 Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                   row.names = fluctuation_rnames)

fluctuation_study <- data.frame("Study" = c(Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                            Fluctuation_Study_Count["Alternating", "Freq"], 
                                            Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                row.names = fluctuation_rnames)


Fluctuation_Model_Estimates_Reorder <- Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

fluctuation_table <- data.frame(estimate = Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                lowerCL = Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                upperCL = Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                K = fluctuation_k[,1], 
                                group_no = fluctuation_group_no[,1], 
                                row.names = fluctuation_rnames)
fluctuation_table$name <- row.names(fluctuation_table)

fluctuation_raw_mean <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                          select("InRR_Transformed"))), 
                          unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                          select("InRR_Transformed"))), 
                          unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                          select("InRR_Transformed"))))

fluctuation_raw_name <- c(replicate(80, "Sinusoidal (Sine Curve)"), 
                          replicate(54, "Alternating"), 
                          replicate(48, "Stepwise"))

fluctuation_raw_df <- data.frame("Model" = fluctuation_raw_name, 
                                 "Effect" = fluctuation_raw_mean)

# Graph code - Combined

Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_fluctuation <- fluctuation_table %>% mutate(name = fct_relevel(name, Fluctuation_Order)) %>%
                       ggplot() +
                       geom_density_ridges(data = fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Fluctuation_Order)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table)[1], 1)), xmin = min(fluctuation_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table)[1], 1)), xmin = max(fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                           size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-2.7, -2.7, -0.8))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                       scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                       coord_cartesian(xlim = c(-0.5, 0.5)) +
                       annotate('text',  x = 0.5, y = (seq(1, dim(fluctuation_table)[1], 1)+0.4),
                       label= paste("italic(k)==", c(fluctuation_table["Stepwise", "K"], 
                                                     fluctuation_table["Alternating", "K"], 
                                                     fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                   c(fluctuation_table["Stepwise", "group_no"], 
                                                     fluctuation_table["Alternating", "group_no"], 
                                                     fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                      x = -0.4, y = (seq(1, dim(fluctuation_table)[1], 1)+0.4)), size = 3.5)
density_fluctuation #(400x320)

# Preparing Graph - Part 1

fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

fluctuation_k_1 <- data.frame("k" = c(Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                      Fluctuation_Exploration["Alternating", "Freq"]), 
                              row.names = fluctuation_rnames_1)

fluctuation_group_no_1 <- data.frame("Spp No." = c(Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                   Fluctuation_Species_Count["Alternating", "Freq"]), 
                                     row.names = fluctuation_rnames_1)

fluctuation_study_1 <- data.frame("Study" = c(Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                              Fluctuation_Study_Count["Alternating", "Freq"]), 
                                  row.names = fluctuation_rnames_1)


Fluctuation_Model_Estimates_Reorder_1 <- Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

fluctuation_table_1 <- data.frame(estimate = Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                  lowerCL = Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                  upperCL = Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                  K = fluctuation_k_1[,1], 
                                  group_no = fluctuation_group_no_1[,1], 
                                  row.names = fluctuation_rnames_1)
fluctuation_table_1$name <- row.names(fluctuation_table_1)

fluctuation_raw_mean_1 <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                          select("InRR_Transformed"))), 
                            unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                          select("InRR_Transformed"))))

fluctuation_raw_name_1 <- c(replicate(80, "Sinusoidal (Sine Curve)"), 
                            replicate(54, "Alternating"))

fluctuation_raw_df_1 <- data.frame("Model" = fluctuation_raw_name_1, 
                                   "Effect" = fluctuation_raw_mean_1)

# Graph code - Part 1

Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_fluctuation_1 <- fluctuation_table_1 %>% mutate(name = fct_relevel(name, Fluctuation_Order_1)) %>%
                         ggplot() +
                         geom_density_ridges(data = fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Fluctuation_Order_1)), 
                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_1)[1], 1)), xmin = min(fluctuation_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_1)[1], 1)), xmin = max(fluctuation_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                        size = 1) +
                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                         size = 1, fatten = 2) +
                         theme_bw() +
                         guides(fill = "none", colour = "none") +
                         labs(x = TeX("Effect Size (PRRD)"), y = "") +
                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                               vjust = c(-2.7, -0.8))) +
                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                         theme(axis.ticks = element_blank()) +
                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                         scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                         scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                         coord_cartesian(xlim = c(-0.5, 0.5)) +
                         annotate('text',  x = 0.5, y = (seq(1, dim(fluctuation_table_1)[1], 1)+0.4),
                         label= paste("italic(k)==", c(fluctuation_table_1["Alternating", "K"], 
                                                       fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                     c(fluctuation_table_1["Alternating", "group_no"], 
                                                       fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                               paste(format(round(mean(exp(Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                   x = -0.4, y = (seq(1, dim(fluctuation_table_1)[1], 1)+0.4)), size = 3.5)
density_fluctuation_1 #(400x240)

# Preparing Graph - Part 2

fluctuation_rnames_2 <- c("Stepwise")

fluctuation_k_2 <- data.frame("k" = c(Fluctuation_Exploration["Stepwise", "Freq"]), 
                              row.names = fluctuation_rnames_2)

fluctuation_group_no_2 <- data.frame("Spp No." = c(Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                     row.names = fluctuation_rnames_2)

fluctuation_study_2 <- data.frame("Study" = c(Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                  row.names = fluctuation_rnames_2)


Fluctuation_Model_Estimates_Reorder_2 <- Fluctuation_Model_Estimates[c("Stepwise"), ]

fluctuation_table_2 <- data.frame(estimate = Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                  lowerCL = Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                  upperCL = Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                  K = fluctuation_k_2[,1], 
                                  group_no = fluctuation_group_no_2[,1], 
                                  row.names = fluctuation_rnames_2)
fluctuation_table_2$name <- row.names(fluctuation_table_2)

fluctuation_raw_mean_2 <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                          select("InRR_Transformed"))))

fluctuation_raw_name_2 <- c(replicate(48, "Stepwise"))

fluctuation_raw_df_2 <- data.frame("Model" = fluctuation_raw_name_2, 
                                   "Effect" = fluctuation_raw_mean_2)

# Graph code - Part 2

Fluctuation_Order_2 <- c("Stepwise")

density_fluctuation_2 <- fluctuation_table_2 %>% mutate(name = fct_relevel(name, Fluctuation_Order_2)) %>%
                         ggplot() +
                         geom_density_ridges(data = fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Fluctuation_Order_2)), 
                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_2)[1], 1)), xmin = min(fluctuation_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_2)[1], 1)), xmin = max(fluctuation_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                        size = 1) +
                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                         coord_cartesian(xlim = c(-0.5, 0.5)) +
                         annotate('text',  x = 0.5, y = (seq(1, dim(fluctuation_table_2)[1], 1)+0.4),
                         label= paste("italic(k)==", c(fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                     c(fluctuation_table_2["Stepwise", "group_no"]), 
                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                    x = -0.4, y = (seq(1, dim(fluctuation_table_2)[1], 1)+0.4)), size = 3.5)
density_fluctuation_2 #(400x160)

##### Overall Model - Trait Meta-Regression #####
Trait_Exploration <- data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Exploration) <- Trait_Exploration$Trait_Category

Trait_Data <- data %>% filter(Trait_Category != "Behavioural" &
                              Trait_Category != "Gene Expression" &
                              Trait_Category != "Population")

Trait_Species_Count <- Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Species_Count) <- Trait_Species_Count$Trait_Category

Trait_Study_Count <- Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
                     filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Study_Count) <- Trait_Study_Count$Trait_Category

Trait_Species <- Trait_Data %>% select("phylo") %>% unique()

Trait_A_cor <- as.data.frame(A_cor)
Trait_A_cor <- Trait_A_cor[c(Trait_Species$phylo), c(Trait_Species$phylo)]
Trait_A_cor <- as.matrix(Trait_A_cor)

Trait_VCV <- make_VCV_matrix(Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Trait_VCV, test = "t", dfs = "contain",
                                   mods = ~ Trait_Category - 1,
                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                   R = list(phylo=Trait_A_cor), data = Trait_Data, method = "REML", sparse = TRUE, 
                                   control=list(rel.tol=1e-9))
    saveRDS(Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Trait_Model.rds")
  } else {
    Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Trait_Model.rds")})

Trait_Model_rob <- robust(Trait_Model, cluster = Trait_Data$Study_ID, adjust = TRUE)

Trait_Model_Estimates <- data.frame(Category = substr(row.names(Trait_Model$b), 15, 100),
                                    estimate = Trait_Model$b, 
                                    ci.lb = Trait_Model$ci.lb, 
                                    ci.ub = Trait_Model$ci.ub)
rownames(Trait_Model_Estimates) <- Trait_Model_Estimates$Category
Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Trait_Model), 2))

# Preparing Graph - Combined

trait_rnames <- c("Biochemical Assay", "Life-history Traits", 
                  "Morphology", "Physiological")

trait_k <- data.frame("k" = c(Trait_Exploration["Biochemical Assay", "Freq"], 
                              Trait_Exploration["Life-History Traits", "Freq"], 
                              Trait_Exploration["Morphology", "Freq"], 
                              Trait_Exploration["Physiological", "Freq"]), 
                      row.names = trait_rnames)

trait_group_no <- data.frame("Spp No." = c(Trait_Species_Count["Biochemical Assay", "Freq"], 
                                           Trait_Species_Count["Life-History Traits", "Freq"],
                                           Trait_Species_Count["Morphology", "Freq"],
                                           Trait_Species_Count["Physiological", "Freq"]), 
                             row.names = trait_rnames)

trait_study <- data.frame("Study" = c(Trait_Study_Count["Biochemical Assay", "Freq"], 
                                      Trait_Study_Count["Life-History Traits", "Freq"],
                                      Trait_Study_Count["Morphology", "Freq"],
                                      Trait_Study_Count["Physiological", "Freq"]), 
                          row.names = trait_rnames)

trait_table <- data.frame(estimate = Trait_Model_Estimates[,"estimate"], 
                          lowerCL = Trait_Model_Estimates[,"ci.lb"], 
                          upperCL = Trait_Model_Estimates[,"ci.ub"], 
                          K = trait_k[,1], 
                          group_no = trait_group_no[,1], 
                          row.names = trait_rnames)
trait_table$name <- row.names(trait_table)

trait_raw_mean <- c(unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                    select("InRR_Transformed"))), 
                    unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                    select("InRR_Transformed"))), 
                    unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                    select("InRR_Transformed"))),
                    unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                    select("InRR_Transformed"))))

trait_raw_name <- c(replicate(32, "Biochemical Assay"), 
                    replicate(68, "Life-history Traits"), 
                    replicate(54, "Morphology"),
                    replicate(41, "Physiological"))

trait_raw_df <- data.frame("Model" = trait_raw_name, 
                           "Effect" = trait_raw_mean)

# Graph code - Combined

Trait_Order <- c("Physiological", "Morphology", "Life-history Traits", "Biochemical Assay")

density_trait <- trait_table %>% mutate(name = fct_relevel(name, Trait_Order)) %>%
                 ggplot() +
                 geom_density_ridges(data = trait_raw_df %>% mutate(Model = fct_relevel(Model, Trait_Order)), 
                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                 geom_linerange(aes(y = rev(seq(1, dim(trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(trait_table)[1], 1)), xmin = min(trait_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(trait_table)[1], 1)), xmin = max(trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                size = 1) +
                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                 theme_bw() +
                 guides(fill = "none", colour = "none") +
                 labs(x = TeX("Effect Size (PRRD)"), y = "") +
                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                   vjust = c(-2.7, -2.7, -0.8, -0.8))) +
                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                 theme(axis.ticks = element_blank()) +
                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                 scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                 scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                 annotate('text',  x = 0.5, y = (seq(1, dim(trait_table)[1], 1)+0.4),
                 label= paste("italic(k)==", c(trait_table["Physiological", "K"], 
                                               trait_table["Morphology", "K"], 
                                               trait_table["Life-history Traits", "K"],
                                               trait_table["Biochemical Assay", "K"]), "~","(", 
                                             c(trait_table["Physiological", "group_no"], 
                                               trait_table["Morphology", "group_no"], 
                                               trait_table["Life-history Traits", "group_no"],
                                               trait_table["Biochemical Assay", "group_no"]), 
                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                 geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                        paste(format(round(mean(exp(Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                x = -0.4, y = (seq(1, dim(trait_table)[1], 1)+0.4)), size = 3.5)

density_trait #(400x400)

# Preparing Graph - Part 1

trait_rnames_1 <- c("Biochemical Assay", "Life-history Traits")

trait_k_1 <- data.frame("k" = c(Trait_Exploration["Biochemical Assay", "Freq"], 
                                Trait_Exploration["Life-History Traits", "Freq"]), 
                        row.names = trait_rnames_1)

trait_group_no_1 <- data.frame("Spp No." = c(Trait_Species_Count["Biochemical Assay", "Freq"], 
                                             Trait_Species_Count["Life-History Traits", "Freq"]), 
                               row.names = trait_rnames_1)

trait_study_1 <- data.frame("Study" = c(Trait_Study_Count["Biochemical Assay", "Freq"], 
                                        Trait_Study_Count["Life-History Traits", "Freq"]), 
                            row.names = trait_rnames_1)

Trait_Model_Estimates_Reorder_1 <- Trait_Model_Estimates[c("Biochemical Assay", "Life-History Traits"), ]

trait_table_1 <- data.frame(estimate = Trait_Model_Estimates_Reorder_1[,"estimate"], 
                            lowerCL = Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                            upperCL = Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                            K = trait_k_1[,1], 
                            group_no = trait_group_no_1[,1], 
                            row.names = trait_rnames_1)
trait_table_1$name <- row.names(trait_table_1)

trait_raw_mean_1 <- c(unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                    select("InRR_Transformed"))), 
                      unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                    select("InRR_Transformed"))))

trait_raw_name_1 <- c(replicate(32, "Biochemical Assay"), 
                      replicate(68, "Life-history Traits"))

trait_raw_df_1 <- data.frame("Model" = trait_raw_name_1, 
                            "Effect" = trait_raw_mean_1)

# Graph code - Part 1

Trait_Order_1 <- c("Life-history Traits", "Biochemical Assay")

density_trait_1 <- trait_table_1 %>% mutate(name = fct_relevel(name, Trait_Order_1)) %>%
                   ggplot() +
                   geom_density_ridges(data = trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Trait_Order_1)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_1)[1], 1)), xmin = min(trait_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_1)[1], 1)), xmin = max(trait_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                   size = 1, fatten = 2) +
                   theme_bw() +
                   guides(fill = "none", colour = "none") +
                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                    vjust = c(-0.8, -0.8))) +
                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                   theme(axis.ticks = element_blank()) +
                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                   annotate('text',  x = 0.5, y = (seq(1, dim(trait_table_1)[1], 1)+0.4),
                   label= paste("italic(k)==", c(trait_table_1["Life-history Traits", "K"],
                                                 trait_table_1["Biochemical Assay", "K"]), "~","(", 
                                               c(trait_table_1["Life-history Traits", "group_no"],
                                                 trait_table_1["Biochemical Assay", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                           paste(format(round(mean(exp(Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                               x = -0.4, y = (seq(1, dim(trait_table_1)[1], 1)+0.4)), size = 3.5)

density_trait_1 #(400x240)

# Preparing Graph - Part 2

trait_rnames_2 <- c("Morphology", "Physiological")

trait_k_2 <- data.frame("k" = c(Trait_Exploration["Morphology", "Freq"], 
                                Trait_Exploration["Physiological", "Freq"]), 
                        row.names = trait_rnames_2)

trait_group_no_2 <- data.frame("Spp No." = c(Trait_Species_Count["Morphology", "Freq"],
                                             Trait_Species_Count["Physiological", "Freq"]), 
                               row.names = trait_rnames_2)

trait_study_2 <- data.frame("Study" = c(Trait_Study_Count["Morphology", "Freq"],
                                        Trait_Study_Count["Physiological", "Freq"]), 
                            row.names = trait_rnames_2)

Trait_Model_Estimates_Reorder_2 <- Trait_Model_Estimates[c("Morphology", "Physiological"), ]

trait_table_2 <- data.frame(estimate = Trait_Model_Estimates_Reorder_2[,"estimate"], 
                            lowerCL = Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                            upperCL = Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                            K = trait_k_2[,1], 
                            group_no = trait_group_no_2[,1], 
                            row.names = trait_rnames_2)
trait_table_2$name <- row.names(trait_table_2)

trait_raw_mean_2 <- c(unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                    select("InRR_Transformed"))),
                      unlist(unname(Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                    select("InRR_Transformed"))))

trait_raw_name_2 <- c(replicate(54, "Morphology"),
                      replicate(41, "Physiological"))

trait_raw_df_2 <- data.frame("Model" = trait_raw_name_2, 
                             "Effect" = trait_raw_mean_2)

# Graph code - Part 2

Trait_Order_2 <- c("Physiological", "Morphology")

density_trait_2 <- trait_table_2 %>% mutate(name = fct_relevel(name, Trait_Order_2)) %>%
                   ggplot() +
                   geom_density_ridges(data = trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Trait_Order_2)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_2)[1], 1)), xmin = min(trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(trait_table_2)[1], 1)), xmin = max(trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                   annotate('text',  x = 0.5, y = (seq(1, dim(trait_table_2)[1], 1)+0.4),
                   label= paste("italic(k)==", c(trait_table["Physiological", "K"], 
                                                 trait_table["Morphology", "K"]), "~","(", 
                                               c(trait_table["Physiological", "group_no"], 
                                                 trait_table["Morphology", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                           paste(format(round(mean(exp(Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                               x = -0.4, y = (seq(1, dim(trait_table_2)[1], 1)+0.4)), size = 3.5)

density_trait_2 #(400x240)

##### Overall Model - Class Meta-Regression #####
Class_Exploration <- data %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Exploration) <- Class_Exploration$Class

Class_Data <- data %>% filter(Class != "Actinopteri" &
                              Class != "Amphibia" &
                              Class != "Anthozoa" &
                              Class != "Branchiopoda" &
                              Class != "Clitellata" &
                              Class != "Gastropoda" &
                              Class != "Holothuroidea")

Class_Species_Count <- Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Species_Count) <- Class_Species_Count$Class

Class_Study_Count <- Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>% 
                     filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Study_Count) <- Class_Study_Count$Class

Class_Species <- Class_Data %>% select("phylo") %>% unique()

Class_A_cor <- as.data.frame(A_cor)
Class_A_cor <- Class_A_cor[c(Class_Species$phylo), c(Class_Species$phylo)]
Class_A_cor <- as.matrix(Class_A_cor)

Class_VCV <- make_VCV_matrix(Class_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Class_Model <- metafor::rma.mv(InRR_Transformed, V = Class_VCV, test = "t", dfs = "contain",
                                   mods = ~ Class - 1,
                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                   R = list(phylo=Class_A_cor), data = Class_Data, method = "REML", sparse = TRUE, 
                                   control=list(rel.tol=1e-9))
    saveRDS(Class_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Class_Model.rds")
  } else {
    Class_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Class_Model.rds")})

Class_Model_rob <- robust(Class_Model, cluster = Class_Data$Study_ID, adjust = TRUE)

Class_Model_Estimates <- data.frame(Class = substr(row.names(Class_Model$b), 6, 100),
                                    estimate = Class_Model$b, 
                                    ci.lb = Class_Model$ci.lb, 
                                    ci.ub = Class_Model$ci.ub)
rownames(Class_Model_Estimates) <- Class_Model_Estimates$Class
Class_Model_i2 <- data.frame(round(orchaRd::i2_ml(Class_Model), 2))

# Preparing Graph - Combined

class_rnames <- c("Arachnida", "Insecta", "Malacostraca")

class_k <- data.frame("k" = c(Class_Exploration["Arachnida", "Freq"], 
                              Class_Exploration["Insecta", "Freq"],
                              Class_Exploration["Malacostraca", "Freq"]), 
                      row.names = class_rnames)

class_group_no <- data.frame("Spp No." = c(Class_Species_Count["Arachnida", "Freq"],
                                           Class_Species_Count["Insecta", "Freq"],
                                           Class_Species_Count["Malacostraca", "Freq"]), 
                             row.names = class_rnames)

class_study <- data.frame("Study" = c(Class_Study_Count["Arachnida", "Freq"],
                                      Class_Study_Count["Insecta", "Freq"],
                                      Class_Study_Count["Malacostraca", "Freq"]), 
                          row.names = class_rnames)

class_table <- data.frame(estimate = Class_Model_Estimates[,"estimate"], 
                          lowerCL = Class_Model_Estimates[,"ci.lb"], 
                          upperCL = Class_Model_Estimates[,"ci.ub"], 
                          K = class_k[,1], 
                          group_no = class_group_no[,1], 
                          row.names = class_rnames)
class_table$name <- row.names(class_table)

class_raw_mean <- c(unlist(unname(Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                    select("InRR_Transformed"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Insecta") %>% 
                                    select("InRR_Transformed"))),
                    unlist(unname(Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                    select("InRR_Transformed"))))

class_raw_name <- c(replicate(11, "Arachnida"), 
                    replicate(108, "Insecta"),
                    replicate(11, "Malacostraca"))

class_raw_df <- data.frame("Model" = class_raw_name, 
                           "Effect" = class_raw_mean)

# Graph code - Combined

Class_Order <- c("Malacostraca", "Insecta", "Arachnida")

density_class <- class_table %>% mutate(name = fct_relevel(name, Class_Order)) %>%
                 ggplot() +
                 geom_density_ridges(data = class_raw_df %>% mutate(Model = fct_relevel(Model, Class_Order)), 
                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                 geom_linerange(aes(y = rev(seq(1, dim(class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(class_table)[1], 1)), xmin = min(class_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                size = 1) +
                 geom_linerange(aes(y = rev(seq(1, dim(class_table)[1], 1)), xmin = max(class_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                size = 1) +
                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                 scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                 scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                 annotate('text',  x = 0.5, y = (seq(1, dim(class_table)[1], 1)+0.4),
                 label= paste("italic(k)==", c(class_table["Malacostraca", "K"], 
                                               class_table["Insecta", "K"], 
                                               class_table["Arachnida", "K"]), "~","(", 
                                             c(class_table["Malacostraca", "group_no"], 
                                               class_table["Insecta", "group_no"], 
                                               class_table["Arachnida", "group_no"]), 
                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                 geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                        paste(format(round(mean(exp(Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                        paste(format(round(mean(exp(Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                x = -0.4, y = (seq(1, dim(class_table)[1], 1)+0.4)), size = 3.5)

density_class #(400x320)

# Preparing Graph - Part 1

class_rnames_1 <- c("Arachnida", "Insecta")

class_k_1 <- data.frame("k" = c(Class_Exploration["Arachnida", "Freq"], 
                                Class_Exploration["Insecta", "Freq"]), 
                        row.names = class_rnames_1)

class_group_no_1 <- data.frame("Spp No." = c(Class_Species_Count["Arachnida", "Freq"],
                                             Class_Species_Count["Insecta", "Freq"]), 
                               row.names = class_rnames_1)

class_study_1 <- data.frame("Study" = c(Class_Study_Count["Arachnida", "Freq"],
                                        Class_Study_Count["Insecta", "Freq"]), 
                            row.names = class_rnames_1)

Class_Model_Estimates_Reorder_1 <- Class_Model_Estimates[c("Arachnida", "Insecta"), ]

class_table_1 <- data.frame(estimate = Class_Model_Estimates_Reorder_1[,"estimate"], 
                            lowerCL = Class_Model_Estimates_Reorder_1[,"ci.lb"], 
                            upperCL = Class_Model_Estimates_Reorder_1[,"ci.ub"], 
                            K = class_k_1[,1], 
                            group_no = class_group_no_1[,1], 
                            row.names = class_rnames_1)
class_table_1$name <- row.names(class_table_1)

class_raw_mean_1 <- c(unlist(unname(Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                    select("InRR_Transformed"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Insecta") %>% 
                                    select("InRR_Transformed"))))

class_raw_name_1 <- c(replicate(11, "Arachnida"), 
                      replicate(108, "Insecta"))

class_raw_df_1 <- data.frame("Model" = class_raw_name_1, 
                             "Effect" = class_raw_mean_1)

# Graph code - Part 1

Class_Order_1 <- c("Insecta", "Arachnida")

density_class_1 <- class_table_1 %>% mutate(name = fct_relevel(name, Class_Order_1)) %>%
                   ggplot() +
                   geom_density_ridges(data = class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Class_Order_1)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_1)[1], 1)), xmin = min(class_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_1)[1], 1)), xmin = max(class_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                   scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                   scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                   annotate('text',  x = 0.5, y = (seq(1, dim(class_table_1)[1], 1)+0.4),
                   label= paste("italic(k)==", c(class_table_1["Insecta", "K"], 
                                                 class_table_1["Arachnida", "K"]), "~","(", 
                                               c(class_table_1["Insecta", "group_no"], 
                                                 class_table_1["Arachnida", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                           paste(format(round(mean(exp(Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                               x = -0.4, y = (seq(1, dim(class_table_1)[1], 1)+0.4)), size = 3.5)

density_class_1 #(400x240)

# Preparing Graph - Part 2

class_rnames_2 <- c("Malacostraca")

class_k_2 <- data.frame("k" = c(Class_Exploration["Malacostraca", "Freq"]), 
                        row.names = class_rnames_2)

class_group_no_2 <- data.frame("Spp No." = c(Class_Species_Count["Malacostraca", "Freq"]), 
                               row.names = class_rnames_2)

class_study_2 <- data.frame("Study" = c(Class_Study_Count["Malacostraca", "Freq"]), 
                            row.names = class_rnames_2)

Class_Model_Estimates_Reorder_2 <- Class_Model_Estimates[c("Malacostraca"), ]

class_table_2 <- data.frame(estimate = Class_Model_Estimates_Reorder_2[,"estimate"], 
                            lowerCL = Class_Model_Estimates_Reorder_2[,"ci.lb"], 
                            upperCL = Class_Model_Estimates_Reorder_2[,"ci.ub"], 
                            K = class_k_2[,1], 
                            group_no = class_group_no_2[,1], 
                            row.names = class_rnames_2)
class_table_2$name <- row.names(class_table_2)

class_raw_mean_2 <- c(unlist(unname(Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                    select("InRR_Transformed"))))

class_raw_name_2 <- c(replicate(11, "Malacostraca"))

class_raw_df_2 <- data.frame("Model" = class_raw_name_2, 
                             "Effect" = class_raw_mean_2)

# Graph code - Part 2

Class_Order_2 <- c("Malacostraca")

density_class_2 <- class_table_2 %>% mutate(name = fct_relevel(name, Class_Order_2)) %>%
                   ggplot() +
                   geom_density_ridges(data = class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Class_Order_2)), 
                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                           scale = 0.015, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_2)[1], 1)), xmin = min(class_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                  size = 1) +
                   geom_linerange(aes(y = rev(seq(1, dim(class_table_2)[1], 1)), xmin = max(class_raw_df_2$Effect)+0.01, xmax = 1.5, colour = name),
                                  size = 1) +
                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                   annotate('text',  x = 0.5, y = (seq(1, dim(class_table_2)[1], 1)+0.4),
                   label= paste("italic(k)==", c(class_table_2["Malacostraca", "K"]), "~","(", 
                                               c(class_table_2["Malacostraca", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                    geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                               x = -0.4, y = (seq(1, dim(class_table_2)[1], 1)+0.4)), size = 3.5)

density_class_2 #(400x160)

##### Overall Model - Specific Trait Meta-Regression #####
Specific_Trait_Exploration <- data %>% select("Measurement") %>% table() %>% data.frame()
Specific_Trait_Exploration <- Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Specific_Trait_Exploration) <- Specific_Trait_Exploration$Measurement

Specific_Trait_Data <- data %>% filter(Measurement == "Development Time"| 
                                       Measurement == "Length"|
                                       Measurement == "Mass"|
                                       Measurement == "Metabolic Rate")

Specific_Trait_Species_Count <- Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Specific_Trait_Species_Count) <- Specific_Trait_Species_Count$Measurement

Specific_Trait_Study_Count <- Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                              filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Specific_Trait_Study_Count) <- Specific_Trait_Study_Count$Measurement

Specific_Trait_Species <- Specific_Trait_Data %>% select("phylo") %>% unique()

Specific_Trait_A_cor <- as.data.frame(A_cor)
Specific_Trait_A_cor <- Specific_Trait_A_cor[c(Specific_Trait_Species$phylo), c(Specific_Trait_Species$phylo)]
Specific_Trait_A_cor <- as.matrix(Specific_Trait_A_cor)

Specific_Trait_VCV <- make_VCV_matrix(Specific_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Specific_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Specific_Trait_VCV, test = "t", dfs = "contain",
                                            mods = ~ Measurement - 1,
                                            random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                          ~1|Shared_Animal_Number), 
                                            R = list(phylo=Specific_Trait_A_cor), data = Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                            control=list(rel.tol=1e-9))
    saveRDS(Specific_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Specific_Trait_Model.rds")
  } else {
    Specific_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Specific_Trait_Model.rds")})

Specific_Trait_Model_rob <- robust(Specific_Trait_Model, cluster = Specific_Trait_Data$Study_ID, adjust = TRUE)

Specific_Trait_Model_Estimates <- data.frame(Trait = substr(row.names(Specific_Trait_Model$b), 12, 100),
                                             estimate = Specific_Trait_Model$b, ci.lb = Specific_Trait_Model$ci.lb, 
                                             ci.ub = Specific_Trait_Model$ci.ub)
rownames(Specific_Trait_Model_Estimates) <- Specific_Trait_Model_Estimates$Trait
Specific_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Specific_Trait_Model), 2))

# Preparing Graph - Combined

specific_trait_rnames <- c("Development Time", "Length", "Mass", "Metabolic Rate")

specific_trait_k <- data.frame("k" = c(Specific_Trait_Exploration["Development Time", "Freq"], 
                                       Specific_Trait_Exploration["Length", "Freq"], 
                                       Specific_Trait_Exploration["Mass", "Freq"], 
                                       Specific_Trait_Exploration["Metabolic Rate", "Freq"]), 
                               row.names = specific_trait_rnames)

specific_trait_group_no <- data.frame("Spp No." = c(Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                    Specific_Trait_Species_Count["Length", "Freq"], 
                                                    Specific_Trait_Species_Count["Mass", "Freq"], 
                                                    Specific_Trait_Species_Count["Metabolic Rate", "Freq"]), 
                                      row.names = specific_trait_rnames)

specific_trait_study <- data.frame("Study" = c(Specific_Trait_Study_Count["Development Time", "Freq"], 
                                               Specific_Trait_Study_Count["Length", "Freq"], 
                                               Specific_Trait_Study_Count["Mass", "Freq"], 
                                               Specific_Trait_Study_Count["Metabolic Rate", "Freq"]), 
                                   row.names = specific_trait_rnames)

Specific_Trait_Model_Estimates_Reorder <- Specific_Trait_Model_Estimates[c("Development Time", "Length", 
                                                                           "Mass", "Metabolic Rate"), ]

specific_trait_table <- data.frame(estimate = Specific_Trait_Model_Estimates_Reorder[,"estimate"], 
                                   lowerCL = Specific_Trait_Model_Estimates_Reorder[,"ci.lb"], 
                                   upperCL = Specific_Trait_Model_Estimates_Reorder[,"ci.ub"], 
                                   K = specific_trait_k[,1], 
                                   group_no = specific_trait_group_no[,1], 
                                   row.names = specific_trait_rnames)
specific_trait_table$name <- row.names(specific_trait_table)

specific_trait_raw_mean <- c(unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                             select("InRR_Transformed"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                             select("InRR_Transformed"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                             select("InRR_Transformed"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                             select("InRR_Transformed"))))

specific_trait_raw_name <- c(replicate(46, "Development Time"), 
                             replicate(14, "Length"), 
                             replicate(25, "Mass"), 
                             replicate(12, "Metabolic Rate"))

specific_trait_raw_df <- data.frame("Model" = specific_trait_raw_name, 
                                    "Effect" = specific_trait_raw_mean)

# Graph code - Combined

Specific_Trait_Order <- c("Metabolic Rate", "Mass", "Length", "Development Time")

density_specific_trait <- specific_trait_table %>% mutate(name = fct_relevel(name, Specific_Trait_Order)) %>%
                          ggplot() +
                          geom_density_ridges(data = specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Specific_Trait_Order)), 
                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                          geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                             size = 1) +
                          geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table)[1], 1)), xmin = min(specific_trait_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                         size = 1) +
                          geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table)[1], 1)), xmin = max(specific_trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                         size = 1) +
                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                              size = 1, fatten = 2) +
                          theme_bw() +
                          guides(fill = "none", colour = "none") +
                          labs(x = TeX("Effect Size (PRRD)"), y = "") +
                          theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                   vjust = c(-0.8, -2.7, -2.7, -0.8))) +
                          theme(axis.text.x = element_text(margin = margin(b = 5))) +
                          theme(axis.ticks = element_blank()) +
                          theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                          theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                          scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                          scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                          scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                          coord_cartesian(xlim = c(-0.5, 0.5)) +
                          annotate('text',  x = 0.5, y = (seq(1, dim(specific_trait_table)[1], 1)+0.4),
                          label= paste("italic(k)==", c(specific_trait_table["Metabolic Rate", "K"], 
                                                        specific_trait_table["Mass", "K"], 
                                                        specific_trait_table["Length", "K"], 
                                                        specific_trait_table["Development Time", "K"]), "~","(", 
                                                      c(specific_trait_table["Metabolic Rate", "group_no"], 
                                                        specific_trait_table["Mass", "group_no"], 
                                                        specific_trait_table["Length", "group_no"], 
                                                        specific_trait_table["Development Time", "group_no"]), 
                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                          geom_label(aes(label=c(paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                 paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                 paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                 paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                         x = -0.4, y = (seq(1, dim(specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_specific_trait #(400x400)

# Preparing Graph - Part 1

specific_trait_rnames_1 <- c("Development Time", "Length")

specific_trait_k_1 <- data.frame("k" = c(Specific_Trait_Exploration["Development Time", "Freq"], 
                                         Specific_Trait_Exploration["Length", "Freq"]), 
                                 row.names = specific_trait_rnames_1)

specific_trait_group_no_1 <- data.frame("Spp No." = c(Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                      Specific_Trait_Species_Count["Length", "Freq"]), 
                                        row.names = specific_trait_rnames_1)

specific_trait_study_1 <- data.frame("Study" = c(Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                 Specific_Trait_Study_Count["Length", "Freq"]), 
                                     row.names = specific_trait_rnames_1)

Specific_Trait_Model_Estimates_Reorder_1 <- Specific_Trait_Model_Estimates[c("Development Time", "Length"), ]

specific_trait_table_1 <- data.frame(estimate = Specific_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                     lowerCL = Specific_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                     upperCL = Specific_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                     K = specific_trait_k_1[,1], 
                                     group_no = specific_trait_group_no_1[,1], 
                                     row.names = specific_trait_rnames_1)
specific_trait_table_1$name <- row.names(specific_trait_table_1)

specific_trait_raw_mean_1 <- c(unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                             select("InRR_Transformed"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                             select("InRR_Transformed"))))

specific_trait_raw_name_1 <- c(replicate(46, "Development Time"), 
                               replicate(14, "Length"))

specific_trait_raw_df_1 <- data.frame("Model" = specific_trait_raw_name_1, 
                                      "Effect" = specific_trait_raw_mean_1)

# Graph code - Part 1

Specific_Trait_Order_1 <- c("Length", "Development Time")

density_specific_trait_1 <- specific_trait_table_1 %>% mutate(name = fct_relevel(name, Specific_Trait_Order_1)) %>%
                            ggplot() +
                            geom_density_ridges(data = specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Specific_Trait_Order_1)), 
                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                           size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_1)[1], 1)), xmin = min(specific_trait_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                           size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_1)[1], 1)), xmin = max(specific_trait_raw_df_1$Effect)+0.01, xmax = 1.5, colour = name),
                                           size = 1) +
                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                            size = 1, fatten = 2) +
                            theme_bw() +
                            guides(fill = "none", colour = "none") +
                            labs(x = TeX("Effect Size (PRRD)"), y = "") +
                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                             vjust = c(-2.7, -0.8))) +
                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                            theme(axis.ticks = element_blank()) +
                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                            scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                            scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                            coord_cartesian(xlim = c(-0.5, 0.5)) +
                            annotate('text',  x = 0.5, y = (seq(1, dim(specific_trait_table_1)[1], 1)+0.4),
                            label= paste("italic(k)==", c(specific_trait_table_1["Length", "K"], 
                                                          specific_trait_table_1["Development Time", "K"]), "~","(", 
                                                        c(specific_trait_table_1["Length", "group_no"], 
                                                          specific_trait_table_1["Development Time", "group_no"]), 
                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.4, y = (seq(1, dim(specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_specific_trait_1 #(400x240)

# Preparing Graph - Part 2

specific_trait_rnames_2 <- c("Mass", "Metabolic Rate")

specific_trait_k_2 <- data.frame("k" = c(Specific_Trait_Exploration["Mass", "Freq"], 
                                         Specific_Trait_Exploration["Metabolic Rate", "Freq"]), 
                                 row.names = specific_trait_rnames_2)

specific_trait_group_no_2 <- data.frame("Spp No." = c(Specific_Trait_Species_Count["Mass", "Freq"], 
                                                      Specific_Trait_Species_Count["Metabolic Rate", "Freq"]), 
                                        row.names = specific_trait_rnames_2)

specific_trait_study_2 <- data.frame("Study" = c(Specific_Trait_Study_Count["Mass", "Freq"], 
                                                 Specific_Trait_Study_Count["Metabolic Rate", "Freq"]), 
                                     row.names = specific_trait_rnames_2)

Specific_Trait_Model_Estimates_Reorder_2 <- Specific_Trait_Model_Estimates[c("Mass", "Metabolic Rate"), ]

specific_trait_table_2 <- data.frame(estimate = Specific_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                     lowerCL = Specific_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                     upperCL = Specific_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                     K = specific_trait_k_2[,1], 
                                     group_no = specific_trait_group_no_2[,1], 
                                     row.names = specific_trait_rnames_2)
specific_trait_table_2$name <- row.names(specific_trait_table_2)

specific_trait_raw_mean_2 <- c(unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                             select("InRR_Transformed"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                             select("InRR_Transformed"))))

specific_trait_raw_name_2 <- c(replicate(25, "Mass"), 
                               replicate(12, "Metabolic Rate"))

specific_trait_raw_df_2 <- data.frame("Model" = specific_trait_raw_name_2, 
                                      "Effect" = specific_trait_raw_mean_2)

# Graph code - Part 2

Specific_Trait_Order_2 <- c("Metabolic Rate", "Mass")

density_specific_trait_2 <- specific_trait_table_2 %>% mutate(name = fct_relevel(name, Specific_Trait_Order_2)) %>%
                            ggplot() +
                            geom_density_ridges(data = specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Specific_Trait_Order_2)), 
                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                           size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_2)[1], 1)), xmin = min(specific_trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                           size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_2)[1], 1)), xmin = max(specific_trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                           size = 1) +
                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                            scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                            scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                            coord_cartesian(xlim = c(-0.5, 0.5)) +
                            annotate('text',  x = 0.5, y = (seq(1, dim(specific_trait_table_2)[1], 1)+0.4),
                            label= paste("italic(k)==", c(specific_trait_table_2["Metabolic Rate", "K"], 
                                                          specific_trait_table_2["Mass", "K"]), "~","(", 
                                                        c(specific_trait_table_2["Metabolic Rate", "group_no"], 
                                                          specific_trait_table_2["Mass", "group_no"]), 
                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.4, y = (seq(1, dim(specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_specific_trait_2 #(400x240)

##### Individual-Level Trait Subset Model #####
Individual_Subset_Data <- data %>% filter(Trait_Category != "Population")
Individual_Species <- Individual_Subset_Data %>% select("phylo") %>% unique()

Individual_A_cor <- as.data.frame(A_cor)
Individual_A_cor <- Individual_A_cor[c(Individual_Species$phylo), c(Individual_Species$phylo)]
Individual_A_cor <- as.matrix(Individual_A_cor)

Individual_VCV <- make_VCV_matrix(Individual_Subset_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Individual_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Individual_VCV, test = "t", dfs = "contain",
                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                      ~1|Shared_Animal_Number, ~1|Measurement), 
                                        R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE,
                                        control=list(rel.tol=1e-9))
    saveRDS(Individual_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Model.rds")
  } else {
    Individual_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Model.rds")})

Individual_Model_rob <- robust(Individual_Model, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)

Individual_Model_Estimates <- data.frame(estimate = Individual_Model$b, ci.lb = Individual_Model$ci.lb, ci.ub = Individual_Model$ci.ub)
Individual_Model_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Model), 2))

#### Individual-Level Subset Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Individual_Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = Individual_VCV, test = "t", dfs = "contain",
                                                  mods = ~ Fluctuation_Magnitude - 1,
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                                  R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(Individual_Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Amplitude_Model.rds")
  } else {
    Individual_Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Amplitude_Model.rds")})

Individual_Amplitude_Model_rob <- robust(Individual_Amplitude_Model, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)

Individual_Amplitude_Model_Estimates <- data.frame(estimate = Individual_Amplitude_Model$b, ci.lb = Individual_Amplitude_Model$ci.lb, 
                                                   ci.ub = Individual_Amplitude_Model$ci.ub)
Individual_Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Amplitude_Model), 2))

# Preparing Graph

Individual_Plot_Data <- Individual_Subset_Data
Individual_Plot_Data <- Individual_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                     ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                     ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Individual_Amplitude_Plot <- ggplot(Individual_Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                             geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                            size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                            shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                             labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                                  size = "Sample Size", title = "Individual-level Traits") +
                             theme_bw() +
                             theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                             theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                             theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                             #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                             geom_hline(yintercept = Individual_Model_Estimates$estimate, lty = 2) + 
                             geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                             stat_poly_eq(formula = y ~ x, 
                                          aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                          parse = TRUE) +
                             coord_cartesian(xlim = c(0, 25), 
                                             ylim = c(-0.25, 0.25))

Individual_Amplitude_Plot #(330x330)

#### Individual-Level Subset Model - Type of Fluctuation Meta-Regression ####
Individual_Fluctuation_Data <- Individual_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Individual_Fluctuation_Exploration <- Individual_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Exploration) <- Individual_Fluctuation_Exploration$Fluctuation_Category

Individual_Fluctuation_Data <- Individual_Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")

Individual_Fluctuation_Species_Count <- Individual_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                        filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Species_Count) <- Individual_Fluctuation_Species_Count$Fluctuation_Category

Individual_Fluctuation_Study_Count <- Individual_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                      filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Study_Count) <- Individual_Fluctuation_Study_Count$Fluctuation_Category

Individual_Fluctuation_Species <- Individual_Fluctuation_Data %>% select("phylo") %>% unique()

Individual_Fluctuation_A_cor <- as.data.frame(A_cor)
Individual_Fluctuation_A_cor <- Individual_Fluctuation_A_cor[c(Individual_Fluctuation_Species$phylo), c(Individual_Fluctuation_Species$phylo)]
Individual_Fluctuation_A_cor <- as.matrix(Individual_Fluctuation_A_cor)

Individual_Fluctuation_VCV <- make_VCV_matrix(Individual_Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Individual_Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Individual_Fluctuation_VCV, test = "t", dfs = "contain",
                                                    mods = ~ Fluctuation_Category - 1,
                                                    random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                  ~1|Shared_Animal_Number, ~1|Measurement), 
                                                    R = list(phylo=Individual_Fluctuation_A_cor), data = Individual_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                    control=list(rel.tol=1e-9))
    saveRDS(Individual_Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Fluctuation_Model.rds")
  } else {
    Individual_Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Fluctuation_Model.rds")})

Individual_Fluctuation_Model_rob <- robust(Individual_Fluctuation_Model, cluster = Individual_Fluctuation_Data$Study_ID, adjust = TRUE)

Individual_Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Individual_Fluctuation_Model$b), 21, 100),
                                                     estimate = Individual_Fluctuation_Model$b, 
                                                     ci.lb = Individual_Fluctuation_Model$ci.lb, 
                                                     ci.ub = Individual_Fluctuation_Model$ci.ub)
rownames(Individual_Fluctuation_Model_Estimates) <- Individual_Fluctuation_Model_Estimates$Category
Individual_Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Fluctuation_Model), 2))

# Preparing Graph - Combined

individual_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

individual_fluctuation_k <- data.frame("k" = c(Individual_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                               Individual_Fluctuation_Exploration["Alternating", "Freq"], 
                                               Individual_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                       row.names = individual_fluctuation_rnames)

individual_fluctuation_group_no <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                            Individual_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                            Individual_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                              row.names = individual_fluctuation_rnames)

individual_fluctuation_study <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                       Individual_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                       Individual_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                           row.names = individual_fluctuation_rnames)

Individual_Fluctuation_Model_Estimates_Reorder <- Individual_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

individual_fluctuation_table <- data.frame(estimate = Individual_Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                           lowerCL = Individual_Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                           upperCL = Individual_Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                           K = individual_fluctuation_k[,1], 
                                           group_no = individual_fluctuation_group_no[,1], 
                                           row.names = individual_fluctuation_rnames)
individual_fluctuation_table$name <- row.names(individual_fluctuation_table)

individual_fluctuation_raw_mean <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InRR_Transformed"))))

individual_fluctuation_raw_name <- c(replicate(74, "Sinusoidal (Sine Curve)"), 
                                     replicate(53, "Alternating"), 
                                     replicate(47, "Stepwise"))

individual_fluctuation_raw_df <- data.frame("Model" = individual_fluctuation_raw_name, 
                                            "Effect" = individual_fluctuation_raw_mean)

# Graph code - Combined

Individual_Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_individual_fluctuation <- individual_fluctuation_table %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = individual_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)), xmin = min(individual_fluctuation_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)), xmin = max(individual_fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                        vjust = c(-2.7, -2.7, -0.8))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.5, 0.5)) +
                                  annotate('text',  x = 0.5, y = (seq(1, dim(individual_fluctuation_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(individual_fluctuation_table["Stepwise", "K"], 
                                                                individual_fluctuation_table["Alternating", "K"], 
                                                                individual_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                              c(individual_fluctuation_table["Stepwise", "group_no"], 
                                                                individual_fluctuation_table["Alternating", "group_no"], 
                                                                individual_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                         paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(individual_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation #(400x320)

# Preparing Graph - Part 1

individual_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

individual_fluctuation_k_1 <- data.frame("k" = c(Individual_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                 Individual_Fluctuation_Exploration["Alternating", "Freq"]), 
                                         row.names = individual_fluctuation_rnames_1)

individual_fluctuation_group_no_1 <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                              Individual_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                row.names = individual_fluctuation_rnames_1)

individual_fluctuation_study_1 <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                         Individual_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                             row.names = individual_fluctuation_rnames_1)

Individual_Fluctuation_Model_Estimates_Reorder_1 <- Individual_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

individual_fluctuation_table_1 <- data.frame(estimate = Individual_Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Individual_Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Individual_Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                             K = individual_fluctuation_k_1[,1], 
                                             group_no = individual_fluctuation_group_no_1[,1], 
                                             row.names = individual_fluctuation_rnames_1)
individual_fluctuation_table_1$name <- row.names(individual_fluctuation_table_1)

individual_fluctuation_raw_mean_1 <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InRR_Transformed"))), 
                                       unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InRR_Transformed"))))

individual_fluctuation_raw_name_1 <- c(replicate(74, "Sinusoidal (Sine Curve)"), 
                                       replicate(53, "Alternating"))

individual_fluctuation_raw_df_1 <- data.frame("Model" = individual_fluctuation_raw_name_1, 
                                              "Effect" = individual_fluctuation_raw_mean_1)

# Graph code - Part 1

Individual_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_individual_fluctuation_1 <- individual_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = individual_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)), xmin = min(individual_fluctuation_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)), xmin = max(individual_fluctuation_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                     vjust = c(-2.7, -0.8))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(individual_fluctuation_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(individual_fluctuation_table_1["Alternating", "K"], 
                                                                  individual_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                c(individual_fluctuation_table_1["Alternating", "group_no"], 
                                                                  individual_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(individual_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation_1 #(400x240)

# Preparing Graph - part 2

individual_fluctuation_rnames_2 <- c("Stepwise")

individual_fluctuation_k_2 <- data.frame("k" = c(Individual_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                         row.names = individual_fluctuation_rnames_2)

individual_fluctuation_group_no_2 <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                row.names = individual_fluctuation_rnames_2)

individual_fluctuation_study_2 <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                             row.names = individual_fluctuation_rnames_2)

Individual_Fluctuation_Model_Estimates_Reorder_2 <- Individual_Fluctuation_Model_Estimates[c("Stepwise"), ]

individual_fluctuation_table_2 <- data.frame(estimate = Individual_Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Individual_Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Individual_Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                             K = individual_fluctuation_k_2[,1], 
                                             group_no = individual_fluctuation_group_no_2[,1], 
                                             row.names = individual_fluctuation_rnames_2)
individual_fluctuation_table_2$name <- row.names(individual_fluctuation_table_2)

individual_fluctuation_raw_mean_2 <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InRR_Transformed"))))

individual_fluctuation_raw_name_2 <- c(replicate(47, "Stepwise"))

individual_fluctuation_raw_df_2 <- data.frame("Model" = individual_fluctuation_raw_name_2, 
                                              "Effect" = individual_fluctuation_raw_mean_2)

# Graph code - Part 2

Individual_Fluctuation_Order_2 <- c("Stepwise")

density_individual_fluctuation_2 <- individual_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = individual_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)), xmin = min(individual_fluctuation_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)), xmin = max(individual_fluctuation_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(individual_fluctuation_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(individual_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                c(individual_fluctuation_table_2["Stepwise", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(individual_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation_2 #(400x160)

##### Individual-Level Subset Model - Class Meta-Regression #####
Individual_Class_Exploration <- Individual_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Exploration) <- Individual_Class_Exploration$Class

Individual_Class_Data <- Individual_Subset_Data %>% filter(Class == "Arachnida"|
                                                           Class == "Insecta")

Individual_Class_Species_Count <- Individual_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>%
                                  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Species_Count) <- Individual_Class_Species_Count$Class

Individual_Class_Study_Count <- Individual_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>%
                                filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Study_Count) <- Individual_Class_Study_Count$Class

Individual_Class_Species <- Individual_Class_Data %>% select("phylo") %>% unique()

Individual_Class_A_cor <- as.data.frame(A_cor)
Individual_Class_A_cor <- Individual_Class_A_cor[c(Individual_Class_Species$phylo), c(Individual_Class_Species$phylo)]
Individual_Class_A_cor <- as.matrix(Individual_Class_A_cor)

Individual_Class_VCV <- make_VCV_matrix(Individual_Class_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Individual_Class_Model <- metafor::rma.mv(InRR_Transformed, V = Individual_Class_VCV, test = "t", dfs = "contain",
                                              mods = ~ Class - 1,
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                            ~1|Shared_Animal_Number, ~1|Measurement), 
                                              R = list(phylo=Individual_Class_A_cor), data = Individual_Class_Data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(Individual_Class_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Class_Model.rds")
  } else {
    Individual_Class_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Individual_Class_Model.rds")})

Individual_Class_Model_rob <- robust(Individual_Class_Model, cluster = Individual_Class_Data$Study_ID, adjust = TRUE)

Individual_Class_Model_Estimates <- data.frame(Class = substr(row.names(Individual_Class_Model$b), 6, 100),
                                               estimate = Individual_Class_Model$b, ci.lb = Individual_Class_Model$ci.lb, 
                                               ci.ub = Individual_Class_Model$ci.ub)
rownames(Individual_Class_Model_Estimates) <- Individual_Class_Model_Estimates$Class
Individual_Class_Model_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Class_Model), 2))

# Preparing Graph - Combined

individual_class_rnames <- c("Arachnida", "Insecta")

individual_class_k <- data.frame("k" = c(Individual_Class_Exploration["Arachnida", "Freq"], 
                                         Individual_Class_Exploration["Insecta", "Freq"]), 
                                 row.names = individual_class_rnames)

individual_class_group_no <- data.frame("Spp No." = c(Individual_Class_Species_Count["Arachnida", "Freq"],
                                                      Individual_Class_Species_Count["Insecta", "Freq"]), 
                                        row.names = individual_class_rnames)

individual_class_study <- data.frame("Study" = c(Individual_Class_Study_Count["Arachnida", "Freq"],
                                                 Individual_Class_Study_Count["Insecta", "Freq"]), 
                                     row.names = individual_class_rnames)

individual_class_table <- data.frame(estimate = Individual_Class_Model_Estimates[,"estimate"], 
                                     lowerCL = Individual_Class_Model_Estimates[,"ci.lb"], 
                                     upperCL = Individual_Class_Model_Estimates[,"ci.ub"], 
                                     K = individual_class_k[,1], 
                                     group_no = individual_class_group_no[,1], 
                                     row.names = individual_class_rnames)
individual_class_table$name <- row.names(individual_class_table)

individual_class_raw_mean <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InRR_Transformed"))), 
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InRR_Transformed"))))

individual_class_raw_name <- c(replicate(11, "Arachnida"), 
                               replicate(104, "Insecta"))

individual_class_raw_df <- data.frame("Model" = individual_class_raw_name, 
                                      "Effect" = individual_class_raw_mean)

# Graph code - Combined

Individual_Class_Order <- c("Insecta", "Arachnida")

density_individual_class <- individual_class_table %>% mutate(name = fct_relevel(name, Individual_Class_Order)) %>%
                            ggplot() +
                            geom_density_ridges(data = individual_class_raw_df %>% mutate(Model = fct_relevel(Model, Individual_Class_Order)), 
                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                            geom_linerange(aes(y = rev(seq(1, dim(individual_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(individual_class_table)[1], 1)), xmin = min(individual_class_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                           size = 1) +
                            geom_linerange(aes(y = rev(seq(1, dim(individual_class_table)[1], 1)), xmin = max(individual_class_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                           size = 1) +
                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                            coord_cartesian(xlim = c(-0.5, 0.5)) +
                            annotate('text',  x = 0.5, y = (seq(1, dim(individual_class_table)[1], 1)+0.4),
                            label= paste("italic(k)==", c(individual_class_table["Insecta", "K"], 
                                                          individual_class_table["Arachnida", "K"]), "~","(", 
                                                        c(individual_class_table["Insecta", "group_no"], 
                                                          individual_class_table["Arachnida", "group_no"]), 
                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                   paste(format(round(mean(exp(Individual_Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(individual_class_table)[1], 1)+0.4)), size = 3.5)

density_individual_class #(400x240)

# Preparing Graph - Part 1

individual_class_rnames_1 <- c("Arachnida")

individual_class_k_1 <- data.frame("k" = c(Individual_Class_Exploration["Arachnida", "Freq"]), 
                                   row.names = individual_class_rnames_1)

individual_class_group_no_1 <- data.frame("Spp No." = c(Individual_Class_Species_Count["Arachnida", "Freq"]), 
                                          row.names = individual_class_rnames_1)

individual_class_study_1 <- data.frame("Study" = c(Individual_Class_Study_Count["Arachnida", "Freq"]), 
                                       row.names = individual_class_rnames_1)

Individual_Class_Model_Estimates_Reorder_1 <- Individual_Class_Model_Estimates[c("Arachnida"), ]

individual_class_table_1 <- data.frame(estimate = Individual_Class_Model_Estimates_Reorder_1[,"estimate"], 
                                     lowerCL = Individual_Class_Model_Estimates_Reorder_1[,"ci.lb"], 
                                     upperCL = Individual_Class_Model_Estimates_Reorder_1[,"ci.ub"], 
                                     K = individual_class_k_1[,1], 
                                     group_no = individual_class_group_no_1[,1], 
                                     row.names = individual_class_rnames_1)
individual_class_table_1$name <- row.names(individual_class_table_1)

individual_class_raw_mean_1 <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InRR_Transformed"))))

individual_class_raw_name_1 <- c(replicate(11, "Arachnida"))

individual_class_raw_df_1 <- data.frame("Model" = individual_class_raw_name_1, 
                                        "Effect" = individual_class_raw_mean_1)

# Graph code - Part 1

Individual_Class_Order_1 <- c("Arachnida")

density_individual_class_1 <- individual_class_table_1 %>% mutate(name = fct_relevel(name, Individual_Class_Order_1)) %>%
                              ggplot() +
                              geom_density_ridges(data = individual_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Individual_Class_Order_1)), 
                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_1)[1], 1)), xmin = min(individual_class_raw_df_1$Effect)-0.06, xmax = -1.5, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_1)[1], 1)), xmin = max(individual_class_raw_df_1$Effect)+0.06, xmax = 1.5, colour = name),
                                             size = 1) +
                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                              coord_cartesian(xlim = c(-0.5, 0.5)) +
                              annotate('text',  x = 0.5, y = (seq(1, dim(individual_class_table_1)[1], 1)+0.4),
                              label= paste("italic(k)==", c(individual_class_table_1["Arachnida", "K"]), "~","(", 
                                                          c(individual_class_table_1["Arachnida", "group_no"]), 
                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(individual_class_table_1)[1], 1)+0.4)), size = 3.5)

density_individual_class_1 #(400x160)

# Preparing Graph - Part 2

individual_class_rnames_2 <- c("Insecta")

individual_class_k_2 <- data.frame("k" = c(Individual_Class_Exploration["Insecta", "Freq"]), 
                                   row.names = individual_class_rnames_2)

individual_class_group_no_2 <- data.frame("Spp No." = c(Individual_Class_Species_Count["Insecta", "Freq"]), 
                                          row.names = individual_class_rnames_2)

individual_class_study_2 <- data.frame("Study" = c(Individual_Class_Study_Count["Insecta", "Freq"]), 
                                       row.names = individual_class_rnames_2)

Individual_Class_Model_Estimates_Reorder_2 <- Individual_Class_Model_Estimates[c("Insecta"), ]

individual_class_table_2 <- data.frame(estimate = Individual_Class_Model_Estimates_Reorder_2[,"estimate"], 
                                     lowerCL = Individual_Class_Model_Estimates_Reorder_2[,"ci.lb"], 
                                     upperCL = Individual_Class_Model_Estimates_Reorder_2[,"ci.ub"], 
                                     K = individual_class_k_2[,1], 
                                     group_no = individual_class_group_no_2[,1], 
                                     row.names = individual_class_rnames_2)
individual_class_table_2$name <- row.names(individual_class_table_2)

individual_class_raw_mean_2 <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InRR_Transformed"))))

individual_class_raw_name_2 <- c(replicate(104, "Insecta"))

individual_class_raw_df_2 <- data.frame("Model" = individual_class_raw_name_2, 
                                        "Effect" = individual_class_raw_mean_2)

# Graph code - Part 2

Individual_Class_Order_2 <- c("Insecta")

density_individual_class_2 <- individual_class_table_2 %>% mutate(name = fct_relevel(name, Individual_Class_Order_2)) %>%
                              ggplot() +
                              geom_density_ridges(data = individual_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Individual_Class_Order_2)), 
                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_2)[1], 1)), xmin = min(individual_class_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_2)[1], 1)), xmin = max(individual_class_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                             size = 1) +
                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                              coord_cartesian(xlim = c(-0.5, 0.5)) +
                              annotate('text',  x = 0.5, y = (seq(1, dim(individual_class_table_2)[1], 1)+0.4),
                              label= paste("italic(k)==", c(individual_class_table_2["Insecta", "K"]), "~","(", 
                                                          c(individual_class_table_2["Insecta", "group_no"]), 
                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                          x = -0.4, y = (seq(1, dim(individual_class_table_2)[1], 1)+0.4)), size = 3.5)

density_individual_class_2 #(400x160)

##### Aquatic Subset Model #####
Aquatic_Subset_Data <- Individual_Subset_Data %>% filter(Ecosystem == "Aquatic")
Aquatic_Species <- Aquatic_Subset_Data %>% select("phylo") %>% unique()

Aquatic_A_cor <- as.data.frame(A_cor)
Aquatic_A_cor <- Aquatic_A_cor[c(Aquatic_Species$phylo), c(Aquatic_Species$phylo)]
Aquatic_A_cor <- as.matrix(Aquatic_A_cor)

Aquatic_VCV <- make_VCV_matrix(Aquatic_Subset_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Aquatic_VCV, test = "t", dfs = "contain",
                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                     R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                     control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Model.rds")
  } else {
    Aquatic_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Model.rds")})

Aquatic_Model_rob <- robust(Aquatic_Model, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Model_Estimates <- data.frame(estimate = Aquatic_Model$b, ci.lb = Aquatic_Model$ci.lb, ci.ub = Aquatic_Model$ci.ub)
Aquatic_Model_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Model), 2))

#### Aquatic Subset Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = Aquatic_VCV, test = "t", dfs = "contain",
                                               mods = ~ Fluctuation_Magnitude - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Amplitude_Model.rds")
  } else {
    Aquatic_Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Amplitude_Model.rds")})

Aquatic_Amplitude_Model_rob <- robust(Aquatic_Amplitude_Model, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Amplitude_Model_Estimates <- data.frame(estimate = Aquatic_Amplitude_Model$b, ci.lb = Aquatic_Amplitude_Model$ci.lb, 
                                                ci.ub = Aquatic_Amplitude_Model$ci.ub)
Aquatic_Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Amplitude_Model), 2))

# Graph Preparing

Aquatic_Plot_Data <- Aquatic_Subset_Data
Aquatic_Plot_Data <- Aquatic_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                               ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                               ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Aquatic_Amplitude_Plot <- ggplot(Aquatic_Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                          geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                         size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                         shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                          labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                               size = "Sample Size", title = "Aquatic Organisms") +
                          theme_bw() +
                          theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                          theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                          theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                          #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                          geom_hline(yintercept = Aquatic_Model_Estimates$estimate, lty = 2) + 
                          geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                          stat_poly_eq(formula = y ~ x, 
                                       aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                       parse = TRUE) +
                          coord_cartesian(xlim = c(0, 25), 
                                          ylim = c(-0.25, 0.25))

Aquatic_Amplitude_Plot #(330x330)

#### Aquatic Subset Model - Type of Fluctuation Meta-Regression ####
Aquatic_Fluctuation_Data <- Aquatic_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Aquatic_Fluctuation_Exploration <- Aquatic_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Exploration) <- Aquatic_Fluctuation_Exploration$Fluctuation_Category

Aquatic_Fluctuation_Data <- Aquatic_Fluctuation_Data %>% filter(Fluctuation_Category != "Stepwise" &
                                                                Fluctuation_Category != "Stochastic")

Aquatic_Fluctuation_Species_Count <- Aquatic_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>% 
                                     filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Species_Count) <- Aquatic_Fluctuation_Species_Count$Fluctuation_Category

Aquatic_Fluctuation_Study_Count <- Aquatic_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>% 
                                   filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Study_Count) <- Aquatic_Fluctuation_Study_Count$Fluctuation_Category

Aquatic_Fluctuation_Species <- Aquatic_Fluctuation_Data %>% select("phylo") %>% unique()

Aquatic_Fluctuation_A_cor <- as.data.frame(A_cor)
Aquatic_Fluctuation_A_cor <- Aquatic_Fluctuation_A_cor[c(Aquatic_Fluctuation_Species$phylo), c(Aquatic_Fluctuation_Species$phylo)]
Aquatic_Fluctuation_A_cor <- as.matrix(Aquatic_Fluctuation_A_cor)

Aquatic_Fluctuation_VCV <- make_VCV_matrix(Aquatic_Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Aquatic_Fluctuation_VCV, test = "t", dfs = "contain",
                                                 mods = ~ Fluctuation_Category - 1,
                                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                               ~1|Shared_Animal_Number, ~1|Measurement), 
                                                 R = list(phylo=Aquatic_Fluctuation_A_cor), data = Aquatic_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Fluctuation_Model.rds")
  } else {
    Aquatic_Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Fluctuation_Model.rds")})

Aquatic_Fluctuation_Model_rob <- robust(Aquatic_Fluctuation_Model, cluster = Aquatic_Fluctuation_Data$Study_ID, adjust = TRUE)

Aquatic_Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Aquatic_Fluctuation_Model$b), 21, 100),
                                                  estimate = Aquatic_Fluctuation_Model$b, 
                                                  ci.lb = Aquatic_Fluctuation_Model$ci.lb, 
                                                  ci.ub = Aquatic_Fluctuation_Model$ci.ub)
rownames(Aquatic_Fluctuation_Model_Estimates) <- Aquatic_Fluctuation_Model_Estimates$Category
Aquatic_Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Fluctuation_Model), 2))

# Preparing Graph - Combined

aquatic_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating")

aquatic_fluctuation_k <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                            Aquatic_Fluctuation_Exploration["Alternating", "Freq"]), 
                                    row.names = aquatic_fluctuation_rnames)

aquatic_fluctuation_group_no <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                         Aquatic_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                           row.names = aquatic_fluctuation_rnames)

aquatic_fluctuation_study <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                    Aquatic_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                        row.names = aquatic_fluctuation_rnames)

Aquatic_Fluctuation_Model_Estimates_Reorder <- Aquatic_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

aquatic_fluctuation_table <- data.frame(estimate = Aquatic_Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                        lowerCL = Aquatic_Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                        upperCL = Aquatic_Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                        K = aquatic_fluctuation_k[,1], 
                                        group_no = aquatic_fluctuation_group_no[,1], 
                                        row.names = aquatic_fluctuation_rnames)
aquatic_fluctuation_table$name <- row.names(aquatic_fluctuation_table)

aquatic_fluctuation_raw_mean <- c(unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                  select("InRR_Transformed"))), 
                                  unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                  select("InRR_Transformed"))))

aquatic_fluctuation_raw_name <- c(replicate(39, "Sinusoidal (Sine Curve)"), 
                                  replicate(11, "Alternating"))

aquatic_fluctuation_raw_df <- data.frame("Model" = aquatic_fluctuation_raw_name, 
                                         "Effect" = aquatic_fluctuation_raw_mean)

# Graph code - Combined

Aquatic_Fluctuation_Order <- c("Alternating", "Sinusoidal (Sine Curve)")

density_aquatic_fluctuation <- aquatic_fluctuation_table %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order)) %>%
                               ggplot() +
                               geom_density_ridges(data = aquatic_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)), xmin = min(aquatic_fluctuation_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)), xmin = max(aquatic_fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                     vjust = c(-2.7, -0.8))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                               scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_fluctuation_table)[1], 1)+0.4),
                               label= paste("italic(k)==", c(aquatic_fluctuation_table["Alternating", "K"], 
                                                             aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                           c(aquatic_fluctuation_table["Alternating", "group_no"], 
                                                             aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                      paste(format(round(mean(exp(Aquatic_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.4, y = (seq(1, dim(aquatic_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation #(400x240)

# Preparing Graph - Part 1

aquatic_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)")

aquatic_fluctuation_k_1 <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"]), 
                                      row.names = aquatic_fluctuation_rnames_1)

aquatic_fluctuation_group_no_1 <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                             row.names = aquatic_fluctuation_rnames_1)

aquatic_fluctuation_study_1 <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                          row.names = aquatic_fluctuation_rnames_1)

Aquatic_Fluctuation_Model_Estimates_Reorder_1 <- Aquatic_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)"), ]

aquatic_fluctuation_table_1 <- data.frame(estimate = Aquatic_Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                          lowerCL = Aquatic_Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                          upperCL = Aquatic_Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                          K = aquatic_fluctuation_k_1[,1], 
                                          group_no = aquatic_fluctuation_group_no_1[,1], 
                                          row.names = aquatic_fluctuation_rnames_1)
aquatic_fluctuation_table_1$name <- row.names(aquatic_fluctuation_table_1)

aquatic_fluctuation_raw_mean_1 <- c(unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                  select("InRR_Transformed"))))

aquatic_fluctuation_raw_name_1 <- c(replicate(39, "Sinusoidal (Sine Curve)"))

aquatic_fluctuation_raw_df_1 <- data.frame("Model" = aquatic_fluctuation_raw_name_1, 
                                           "Effect" = aquatic_fluctuation_raw_mean_1)

# Graph code - Part 1

Aquatic_Fluctuation_Order_1 <- c("Sinusoidal (Sine Curve)")

density_aquatic_fluctuation_1 <- aquatic_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order_1)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = aquatic_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order_1)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.025, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)), xmin = min(aquatic_fluctuation_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)), xmin = max(aquatic_fluctuation_raw_df_1$Effect)+0.01, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 scale_colour_manual(values = c("#2B4E7A")) +
                                 scale_fill_manual(values = c("#2B4E7A")) +
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_fluctuation_table_1)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(aquatic_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                             c(aquatic_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(aquatic_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation_1 #(400x160)

# Preparing Graph - Part 2

aquatic_fluctuation_rnames_2 <- c("Alternating")

aquatic_fluctuation_k_2 <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Alternating", "Freq"]), 
                                      row.names = aquatic_fluctuation_rnames_2)

aquatic_fluctuation_group_no_2 <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                             row.names = aquatic_fluctuation_rnames_2)

aquatic_fluctuation_study_2 <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                          row.names = aquatic_fluctuation_rnames_2)

Aquatic_Fluctuation_Model_Estimates_Reorder_2 <- Aquatic_Fluctuation_Model_Estimates[c("Alternating"), ]

aquatic_fluctuation_table_2 <- data.frame(estimate = Aquatic_Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                          lowerCL = Aquatic_Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                          upperCL = Aquatic_Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                          K = aquatic_fluctuation_k_2[,1], 
                                          group_no = aquatic_fluctuation_group_no_2[,1], 
                                          row.names = aquatic_fluctuation_rnames_2)
aquatic_fluctuation_table_2$name <- row.names(aquatic_fluctuation_table_2)

aquatic_fluctuation_raw_mean_2 <- c(unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                  select("InRR_Transformed"))))

aquatic_fluctuation_raw_name_2 <- c(replicate(11, "Alternating"))

aquatic_fluctuation_raw_df_2 <- data.frame("Model" = aquatic_fluctuation_raw_name_2, 
                                           "Effect" = aquatic_fluctuation_raw_mean_2)

# Graph code - Part 2

Aquatic_Fluctuation_Order_2 <- c("Alternating")

density_aquatic_fluctuation_2 <- aquatic_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order_2)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = aquatic_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order_2)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)), xmin = min(aquatic_fluctuation_raw_df_2$Effect)-0.06, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)), xmin = max(aquatic_fluctuation_raw_df_2$Effect)+0.06, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_fluctuation_table_2)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(aquatic_fluctuation_table_2["Alternating", "K"]), "~","(", 
                                                             c(aquatic_fluctuation_table_2["Alternating", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(aquatic_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation_2 #(400x160)

##### Aquatic Subset Model - Trait Meta-Regression #####
Aquatic_Trait_Exploration <- Aquatic_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Exploration) <- Aquatic_Trait_Exploration$Trait_Category

Aquatic_Trait_Data <- Aquatic_Subset_Data %>% filter(Trait_Category != "Behavioural" &
                                                     Trait_Category != "Biochemical Assay" &
                                                     Trait_Category != "Morphology")

Aquatic_Trait_Species_Count <- Aquatic_Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
                               filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Species_Count) <- Aquatic_Trait_Species_Count$Trait_Category

Aquatic_Trait_Study_Count <- Aquatic_Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
                             filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Study_Count) <- Aquatic_Trait_Study_Count$Trait_Category

Aquatic_Trait_Species <- Aquatic_Trait_Data %>% select("phylo") %>% unique()

Aquatic_Trait_A_cor <- as.data.frame(A_cor)
Aquatic_Trait_A_cor <- Aquatic_Trait_A_cor[c(Aquatic_Trait_Species$phylo), c(Aquatic_Trait_Species$phylo)]
Aquatic_Trait_A_cor <- as.matrix(Aquatic_Trait_A_cor)

Aquatic_Trait_VCV <- make_VCV_matrix(Aquatic_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Aquatic_Trait_VCV, test = "t", dfs = "contain",
                                           mods = ~ Trait_Category - 1,
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                                           R = list(phylo=Aquatic_Trait_A_cor), data = Aquatic_Trait_Data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Trait_Model.rds")
  } else {
    Aquatic_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Trait_Model.rds")})

Aquatic_Trait_Model_rob <- robust(Aquatic_Trait_Model, cluster = Aquatic_Trait_Data$Study_ID, adjust = TRUE)

Aquatic_Trait_Model_Estimates <- data.frame(Category = substr(row.names(Aquatic_Trait_Model$b), 15, 100),
                                            estimate = Aquatic_Trait_Model$b, ci.lb = Aquatic_Trait_Model$ci.lb, 
                                            ci.ub = Aquatic_Trait_Model$ci.ub)
rownames(Aquatic_Trait_Model_Estimates) <- Aquatic_Trait_Model_Estimates$Category
Aquatic_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Trait_Model), 2))

# Preparing Graph - Combined

aquatic_trait_rnames <- c("Life-history Traits", "Physiological")

aquatic_trait_k <- data.frame("k" = c(Aquatic_Trait_Exploration["Life-History Traits", "Freq"], 
                                      Aquatic_Trait_Exploration["Physiological", "Freq"]), 
                              row.names = aquatic_trait_rnames)

aquatic_trait_group_no <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Life-History Traits", "Freq"],
                                                   Aquatic_Trait_Species_Count["Physiological", "Freq"]), 
                                     row.names = aquatic_trait_rnames)

aquatic_trait_study <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Life-History Traits", "Freq"],
                                              Aquatic_Trait_Study_Count["Physiological", "Freq"]), 
                                  row.names = aquatic_trait_rnames)

aquatic_trait_table <- data.frame(estimate = Aquatic_Trait_Model_Estimates[,"estimate"], 
                                  lowerCL = Aquatic_Trait_Model_Estimates[,"ci.lb"], 
                                  upperCL = Aquatic_Trait_Model_Estimates[,"ci.ub"], 
                                  K = aquatic_trait_k[,1], 
                                  group_no = aquatic_trait_group_no[,1], 
                                  row.names = aquatic_trait_rnames)
aquatic_trait_table$name <- row.names(aquatic_trait_table)

aquatic_trait_raw_mean <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                            select("InRR_Transformed"))), 
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                            select("InRR_Transformed"))))

aquatic_trait_raw_name <- c(replicate(15, "Life-history Traits"), 
                            replicate(23, "Physiological"))

aquatic_trait_raw_df <- data.frame("Model" = aquatic_trait_raw_name, 
                                   "Effect" = aquatic_trait_raw_mean)

# Graph code - Combined

Aquatic_Trait_Order <- c("Physiological", "Life-history Traits")

density_aquatic_trait <- aquatic_trait_table %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order)) %>%
                         ggplot() +
                         geom_density_ridges(data = aquatic_trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order)), 
                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                         geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table)[1], 1)), xmin = min(aquatic_trait_raw_df$Effect)-0.04, xmax = -1.5, colour = name),
                                        size = 1) +
                         geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table)[1], 1)), xmin = max(aquatic_trait_raw_df$Effect)+0.04, xmax = 1.5, colour = name),
                                        size = 1) +
                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                             size = 1, fatten = 2) +
                         theme_bw() +
                         guides(fill = "none", colour = "none") +
                         labs(x = TeX("Effect Size (PRRD)"), y = "") +
                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                               vjust = c(-2.7, -0.8))) +
                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                         theme(axis.ticks = element_blank()) +
                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                         scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                         scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                         coord_cartesian(xlim = c(-0.5, 0.5)) +
                         annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_trait_table)[1], 1)+0.4),
                         label= paste("italic(k)==", c(aquatic_trait_table["Physiological", "K"], 
                                                       aquatic_trait_table["Life-history Traits", "K"]), "~","(", 
                                                     c(aquatic_trait_table["Physiological", "group_no"], 
                                                       aquatic_trait_table["Life-history Traits", "group_no"]), 
                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Aquatic_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.4, y = (seq(1, dim(aquatic_trait_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait #(400x240)

# Preparing Graph - Part 1

aquatic_trait_rnames_1 <- c("Life-history Traits")

aquatic_trait_k_1 <- data.frame("k" = c(Aquatic_Trait_Exploration["Life-History Traits", "Freq"]), 
                                row.names = aquatic_trait_rnames_1)

aquatic_trait_group_no_1 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Life-History Traits", "Freq"]), 
                                       row.names = aquatic_trait_rnames_1)

aquatic_trait_study_1 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Life-History Traits", "Freq"]), 
                                    row.names = aquatic_trait_rnames_1)

Aquatic_Trait_Model_Estimates_Reorder_1 <- Aquatic_Trait_Model_Estimates[c("Life-History Traits"), ]

aquatic_trait_table_1 <- data.frame(estimate = Aquatic_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                    lowerCL = Aquatic_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                    upperCL = Aquatic_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                    K = aquatic_trait_k_1[,1], 
                                    group_no = aquatic_trait_group_no_1[,1], 
                                    row.names = aquatic_trait_rnames_1)
aquatic_trait_table_1$name <- row.names(aquatic_trait_table_1)

aquatic_trait_raw_mean_1 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                            select("InRR_Transformed"))))

aquatic_trait_raw_name_1 <- c(replicate(15, "Life-history Traits"))

aquatic_trait_raw_df_1 <- data.frame("Model" = aquatic_trait_raw_name_1, 
                                     "Effect" = aquatic_trait_raw_mean_1)

# Graph code - Part 1

Aquatic_Trait_Order_1 <- c("Life-history Traits")

density_aquatic_trait_1 <- aquatic_trait_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_1)) %>%
                           ggplot() +
                           geom_density_ridges(data = aquatic_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_1)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)), xmin = min(aquatic_trait_raw_df_1$Effect)-0.07, xmax = -1.5, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)), xmin = max(aquatic_trait_raw_df_1$Effect)+0.07, xmax = 1.5, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                           scale_colour_manual(values = c("#2B4E7A")) +
                           scale_fill_manual(values = c("#2B4E7A")) +
                           coord_cartesian(xlim = c(-0.5, 0.5)) +
                           annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_trait_table_1)[1], 1)+0.4),
                           label= paste("italic(k)==", c(aquatic_trait_table["Life-history Traits", "K"]), "~","(", 
                                                       c(aquatic_trait_table["Life-history Traits", "group_no"]), 
                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                       x = -0.4, y = (seq(1, dim(aquatic_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait_1 #(400x160)

# Preparing Graph - Part 2

aquatic_trait_rnames_2 <- c("Physiological")

aquatic_trait_k_2 <- data.frame("k" = c(Aquatic_Trait_Exploration["Physiological", "Freq"]), 
                                row.names = aquatic_trait_rnames_2)

aquatic_trait_group_no_2 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Physiological", "Freq"]), 
                                       row.names = aquatic_trait_rnames_2)

aquatic_trait_study_2 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Physiological", "Freq"]), 
                                    row.names = aquatic_trait_rnames_2)

Aquatic_Trait_Model_Estimates_Reorder_2 <- Aquatic_Trait_Model_Estimates[c("Physiological"), ]

aquatic_trait_table_2 <- data.frame(estimate = Aquatic_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                  lowerCL = Aquatic_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                  upperCL = Aquatic_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                  K = aquatic_trait_k_2[,1], 
                                  group_no = aquatic_trait_group_no_2[,1], 
                                  row.names = aquatic_trait_rnames_2)
aquatic_trait_table_2$name <- row.names(aquatic_trait_table_2)

aquatic_trait_raw_mean_2 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                            select("InRR_Transformed"))))

aquatic_trait_raw_name_2 <- c(replicate(23, "Physiological"))

aquatic_trait_raw_df_2 <- data.frame("Model" = aquatic_trait_raw_name_2, 
                                     "Effect" = aquatic_trait_raw_mean_2)

# Graph code - Part 2

Aquatic_Trait_Order_2 <- c("Physiological")

density_aquatic_trait_2 <- aquatic_trait_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_2)) %>%
                           ggplot() +
                           geom_density_ridges(data = aquatic_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_2)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)), xmin = min(aquatic_trait_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                          size = 1) +
                           geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)), xmin = max(aquatic_trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                           coord_cartesian(xlim = c(-0.5, 0.5)) +
                           annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_trait_table_2)[1], 1)+0.4),
                           label= paste("italic(k)==", c(aquatic_trait_table_2["Physiological", "K"]), "~","(", 
                                                       c(aquatic_trait_table_2["Physiological", "group_no"]), 
                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                       x = -0.4, y = (seq(1, dim(aquatic_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_trait_2 #(400x160)

##### Aquatic Subset Model - Plasticity Mechanism Meta-Regression #####
Aquatic_Plasticity_Exploration <- Aquatic_Subset_Data %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Exploration) <- Aquatic_Plasticity_Exploration$Plasticity_Mechanism

Aquatic_Plasticity_Species_Count <- Aquatic_Subset_Data %>% select("Scientific_Name", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Species_Count) <- Aquatic_Plasticity_Species_Count$Plasticity_Mechanism

Aquatic_Plasticity_Study_Count <- Aquatic_Subset_Data %>% select("Study_ID", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Study_Count) <- Aquatic_Plasticity_Study_Count$Plasticity_Mechanism

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Plasticity_Model <- metafor::rma.mv(InRR_Transformed, V = Aquatic_VCV, test = "t", dfs = "contain",
                                                mods = ~ Plasticity_Mechanism - 1,
                                                random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                              ~1|Shared_Animal_Number, ~1|Measurement), 
                                                R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                                control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Plasticity_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Plasticity_Model.rds")
  } else {
    Aquatic_Plasticity_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Aquatic_Plasticity_Model.rds")})

Aquatic_Plasticity_Model_rob <- robust(Aquatic_Plasticity_Model, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Plasticity_Model_Estimates <- data.frame(Plasticity_Mechanism = substr(row.names(Aquatic_Plasticity_Model$b), 21, 100),
                                                 estimate = Aquatic_Plasticity_Model$b, ci.lb = Aquatic_Plasticity_Model$ci.lb, 
                                                 ci.ub = Aquatic_Plasticity_Model$ci.ub)
rownames(Aquatic_Plasticity_Model_Estimates) <- Aquatic_Plasticity_Model_Estimates$Plasticity_Mechanism
Aquatic_Plasticity_Model_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Plasticity_Model), 2))

# Preparing Graph - Combined

aquatic_plasticity_rnames <- c("Acclimation", "Development")

aquatic_plasticity_k <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"], 
                                           Aquatic_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                   row.names = aquatic_plasticity_rnames)

aquatic_plasticity_group_no <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                        Aquatic_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                          row.names = aquatic_plasticity_rnames)

aquatic_plasticity_study <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                   Aquatic_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                       row.names = aquatic_plasticity_rnames)

aquatic_plasticity_table <- data.frame(estimate = Aquatic_Plasticity_Model_Estimates[,"estimate"], 
                                       lowerCL = Aquatic_Plasticity_Model_Estimates[,"ci.lb"], 
                                       upperCL = Aquatic_Plasticity_Model_Estimates[,"ci.ub"], 
                                       K = aquatic_plasticity_k[,1], 
                                       group_no = aquatic_plasticity_group_no[,1], 
                                       row.names = aquatic_plasticity_rnames)
aquatic_plasticity_table$name <- row.names(aquatic_plasticity_table)

aquatic_plasticity_raw_mean <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                 select("InRR_Transformed"))), 
                                 unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                 select("InRR_Transformed"))))

aquatic_plasticity_raw_name <- c(replicate(28, "Acclimation"), 
                                 replicate(23, "Development"))

aquatic_plasticity_raw_df <- data.frame("Model" = aquatic_plasticity_raw_name, 
                                        "Effect" = aquatic_plasticity_raw_mean)

# Graph code - Combined

Aquatic_Plasticity_Order <- c("Development", "Acclimation")

density_aquatic_plasticity <- aquatic_plasticity_table %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order)) %>%
                              ggplot() +
                              geom_density_ridges(data = aquatic_plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order)), 
                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                              geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)), xmin = min(aquatic_plasticity_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                             size = 1) +
                              geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)), xmin = max(aquatic_plasticity_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                             size = 1) +
                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                              coord_cartesian(xlim = c(-0.5, 0.5)) +
                              annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_plasticity_table)[1], 1)+0.4),
                              label= paste("italic(k)==", c(aquatic_plasticity_table["Development", "K"], 
                                                            aquatic_plasticity_table["Acclimation", "K"]), "~","(", 
                                                          c(aquatic_plasticity_table["Development", "group_no"], 
                                                            aquatic_plasticity_table["Acclimation", "group_no"]), 
                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                              geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Aquatic_Plasticity_Model_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(aquatic_plasticity_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity #(400x240)

# Preparing Graph - Part 1

aquatic_plasticity_rnames_1 <- c("Acclimation")

aquatic_plasticity_k_1 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"]), 
                                     row.names = aquatic_plasticity_rnames_1)

aquatic_plasticity_group_no_1 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                            row.names = aquatic_plasticity_rnames_1)

aquatic_plasticity_study_1 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                         row.names = aquatic_plasticity_rnames_1)

Aquatic_Plasticity_Model_Estimates_Reorder_1 <- Aquatic_Plasticity_Model_Estimates[c("Acclimation"), ]

aquatic_plasticity_table_1 <- data.frame(estimate = Aquatic_Plasticity_Model_Estimates_Reorder_1[,"estimate"], 
                                         lowerCL = Aquatic_Plasticity_Model_Estimates_Reorder_1[,"ci.lb"], 
                                         upperCL = Aquatic_Plasticity_Model_Estimates_Reorder_1[,"ci.ub"], 
                                         K = aquatic_plasticity_k_1[,1], 
                                         group_no = aquatic_plasticity_group_no_1[,1], 
                                         row.names = aquatic_plasticity_rnames_1)
aquatic_plasticity_table_1$name <- row.names(aquatic_plasticity_table_1)

aquatic_plasticity_raw_mean_1 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                 select("InRR_Transformed"))))

aquatic_plasticity_raw_name_1 <- c(replicate(28, "Acclimation"))

aquatic_plasticity_raw_df_1 <- data.frame("Model" = aquatic_plasticity_raw_name_1, 
                                          "Effect" = aquatic_plasticity_raw_mean_1)

# Graph code - Part 1

Aquatic_Plasticity_Order_1 <- c("Acclimation")

density_aquatic_plasticity_1 <- aquatic_plasticity_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_1)) %>%
                                ggplot() +
                                geom_density_ridges(data = aquatic_plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_1)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)), xmin = min(aquatic_plasticity_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)), xmin = max(aquatic_plasticity_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                coord_cartesian(xlim = c(-0.5, 0.5)) +
                                annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_plasticity_table_1)[1], 1)+0.4),
                                label= paste("italic(k)==", c(aquatic_plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                            c(aquatic_plasticity_table_1["Acclimation", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(aquatic_plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity_1 #(400x160)

# Preparing Graph - Part 2

aquatic_plasticity_rnames_2 <- c("Development")

aquatic_plasticity_k_2 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                     row.names = aquatic_plasticity_rnames_2)

aquatic_plasticity_group_no_2 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                            row.names = aquatic_plasticity_rnames_2)

aquatic_plasticity_study_2 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                         row.names = aquatic_plasticity_rnames_2)

Aquatic_Plasticity_Model_Estimates_Reorder_2 <- Aquatic_Plasticity_Model_Estimates[c("Development"), ]

aquatic_plasticity_table_2 <- data.frame(estimate = Aquatic_Plasticity_Model_Estimates_Reorder_2[,"estimate"], 
                                       lowerCL = Aquatic_Plasticity_Model_Estimates_Reorder_2[,"ci.lb"], 
                                       upperCL = Aquatic_Plasticity_Model_Estimates_Reorder_2[,"ci.ub"], 
                                       K = aquatic_plasticity_k_2[,1], 
                                       group_no = aquatic_plasticity_group_no_2[,1], 
                                       row.names = aquatic_plasticity_rnames_2)
aquatic_plasticity_table_2$name <- row.names(aquatic_plasticity_table_2)

aquatic_plasticity_raw_mean_2 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                 select("InRR_Transformed"))))

aquatic_plasticity_raw_name_2 <- c(replicate(23, "Development"))

aquatic_plasticity_raw_df_2 <- data.frame("Model" = aquatic_plasticity_raw_name_2, 
                                          "Effect" = aquatic_plasticity_raw_mean_2)

# Graph code - Part 2

Aquatic_Plasticity_Order_2 <- c("Development")

density_aquatic_plasticity_2 <- aquatic_plasticity_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_2)) %>%
                                ggplot() +
                                geom_density_ridges(data = aquatic_plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_2)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)), xmin = min(aquatic_plasticity_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)), xmin = max(aquatic_plasticity_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                coord_cartesian(xlim = c(-0.5, 0.5)) +
                                annotate('text',  x = 0.5, y = (seq(1, dim(aquatic_plasticity_table_2)[1], 1)+0.4),
                                label= paste("italic(k)==", c(aquatic_plasticity_table_2["Development", "K"]), "~","(", 
                                                            c(aquatic_plasticity_table_2["Development", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(aquatic_plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity_2 #(400x160)

##### Terrestrial Subset Model #####
Terrestrial_Subset_Data <- Individual_Subset_Data %>% filter(Ecosystem == "Terrestrial")
Terrestrial_Species <- Terrestrial_Subset_Data %>% select("phylo") %>% unique()

Terrestrial_A_cor <- as.data.frame(A_cor)
Terrestrial_A_cor <- Terrestrial_A_cor[c(Terrestrial_Species$phylo), c(Terrestrial_Species$phylo)]
Terrestrial_A_cor <- as.matrix(Terrestrial_A_cor)

Terrestrial_VCV <- make_VCV_matrix(Terrestrial_Subset_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Terrestrial_VCV, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                         R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Model.rds")
  } else {
    Terrestrial_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Model.rds")})

Terrestrial_Model_rob <- robust(Terrestrial_Model, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Model_Estimates <- data.frame(estimate = Terrestrial_Model$b, ci.lb = Terrestrial_Model$ci.lb, ci.ub = Terrestrial_Model$ci.ub)
Terrestrial_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Model), 2))

#### Terrestrial Subset Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = Terrestrial_VCV, test = "t", dfs = "contain",
                                                   mods = ~ Fluctuation_Magnitude - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Amplitude_Model.rds")
  } else {
    Terrestrial_Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Amplitude_Model.rds")})

Terrestrial_Amplitude_Model_rob <- robust(Terrestrial_Amplitude_Model, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Amplitude_Model_Estimates <- data.frame(estimate = Terrestrial_Amplitude_Model$b, ci.lb = Terrestrial_Amplitude_Model$ci.lb, 
                                                    ci.ub = Terrestrial_Amplitude_Model$ci.ub)
Terrestrial_Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Amplitude_Model), 2))

# Graph Preparing

Terrestrial_Plot_Data <- Terrestrial_Subset_Data
Terrestrial_Plot_Data <- Terrestrial_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                       ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                       ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Terrestrial_Amplitude_Plot <- ggplot(Terrestrial_Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                              geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                             size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                              labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                                   size = "Sample Size", title = "Terrestrial Organisms") +
                              theme_bw() +
                              theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                              theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                              theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                              #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                              geom_hline(yintercept = Terrestrial_Model_Estimates$estimate, lty = 2) + 
                              geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                              stat_poly_eq(formula = y ~ x, 
                                           aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                           parse = TRUE) +
                              coord_cartesian(xlim = c(0, 25), 
                                              ylim = c(-0.25, 0.25))

Terrestrial_Amplitude_Plot #(330x330)

#### Terrestrial Subset Model - Type of Fluctuation Meta-Regression ####
Terrestrial_Fluctuation_Data <- Terrestrial_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Terrestrial_Fluctuation_Exploration <- Terrestrial_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Exploration) <- Terrestrial_Fluctuation_Exploration$Fluctuation_Category

Terrestrial_Fluctuation_Data <- Terrestrial_Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")

Terrestrial_Fluctuation_Species_Count <- Terrestrial_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                         filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Species_Count) <- Terrestrial_Fluctuation_Species_Count$Fluctuation_Category

Terrestrial_Fluctuation_Study_Count <- Terrestrial_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                       filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Study_Count) <- Terrestrial_Fluctuation_Study_Count$Fluctuation_Category

Terrestrial_Fluctuation_Species <- Terrestrial_Fluctuation_Data %>% select("phylo") %>% unique()

Terrestrial_Fluctuation_A_cor <- as.data.frame(A_cor)
Terrestrial_Fluctuation_A_cor <- Terrestrial_Fluctuation_A_cor[c(Terrestrial_Fluctuation_Species$phylo), c(Terrestrial_Fluctuation_Species$phylo)]
Terrestrial_Fluctuation_A_cor <- as.matrix(Terrestrial_Fluctuation_A_cor)

Terrestrial_Fluctuation_VCV <- make_VCV_matrix(Terrestrial_Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Terrestrial_Fluctuation_VCV, test = "t", dfs = "contain",
                                                     mods = ~ Fluctuation_Category - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Terrestrial_Fluctuation_A_cor), data = Terrestrial_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Fluctuation_Model.rds")
  } else {
    Terrestrial_Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Fluctuation_Model.rds")})

Terrestrial_Fluctuation_Model_rob <- robust(Terrestrial_Fluctuation_Model, cluster = Terrestrial_Fluctuation_Data$Study_ID, adjust = TRUE)

Terrestrial_Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Terrestrial_Fluctuation_Model$b), 21, 100),
                                                      estimate = Terrestrial_Fluctuation_Model$b, ci.lb = Terrestrial_Fluctuation_Model$ci.lb, 
                                                      ci.ub = Terrestrial_Fluctuation_Model$ci.ub)
rownames(Terrestrial_Fluctuation_Model_Estimates) <- Terrestrial_Fluctuation_Model_Estimates$Category
Terrestrial_Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Fluctuation_Model), 2))

# Preparing Graph - Combined

terrestrial_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

terrestrial_fluctuation_k <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                Terrestrial_Fluctuation_Exploration["Alternating", "Freq"], 
                                                Terrestrial_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                        row.names = terrestrial_fluctuation_rnames)

terrestrial_fluctuation_group_no <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                             Terrestrial_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                             Terrestrial_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                               row.names = terrestrial_fluctuation_rnames)

terrestrial_fluctuation_study <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                        Terrestrial_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                        Terrestrial_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                            row.names = terrestrial_fluctuation_rnames)

Terrestrial_Fluctuation_Model_Estimates_Reorder <- Terrestrial_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

terrestrial_fluctuation_table <- data.frame(estimate = Terrestrial_Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                            lowerCL = Terrestrial_Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                            upperCL = Terrestrial_Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                            K = terrestrial_fluctuation_k[,1], 
                                            group_no = terrestrial_fluctuation_group_no[,1], 
                                            row.names = terrestrial_fluctuation_rnames)
terrestrial_fluctuation_table$name <- row.names(terrestrial_fluctuation_table)

terrestrial_fluctuation_raw_mean <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InRR_Transformed"))), 
                                      unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InRR_Transformed"))), 
                                      unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InRR_Transformed"))))

terrestrial_fluctuation_raw_name <- c(replicate(35, "Sinusoidal (Sine Curve)"), 
                                      replicate(42, "Alternating"), 
                                      replicate(46, "Stepwise"))

terrestrial_fluctuation_raw_df <- data.frame("Model" = terrestrial_fluctuation_raw_name, 
                                             "Effect" = terrestrial_fluctuation_raw_mean)

# Graph code - Combined

Terrestrial_Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_terrestrial_fluctuation <- terrestrial_fluctuation_table %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = terrestrial_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)), xmin = min(terrestrial_fluctuation_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)), xmin = max(terrestrial_fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7, -0.8))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                                   annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_fluctuation_table)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(terrestrial_fluctuation_table["Stepwise", "K"], 
                                                                 terrestrial_fluctuation_table["Alternating", "K"], 
                                                                 terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                               c(terrestrial_fluctuation_table["Stepwise", "group_no"], 
                                                                 terrestrial_fluctuation_table["Alternating", "group_no"], 
                                                                 terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.4, y = (seq(1, dim(terrestrial_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation #(400x320)

# Preparing Graph - Part 1

terrestrial_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

terrestrial_fluctuation_k_1 <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                Terrestrial_Fluctuation_Exploration["Alternating", "Freq"]), 
                                          row.names = terrestrial_fluctuation_rnames_1)

terrestrial_fluctuation_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                               Terrestrial_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                 row.names = terrestrial_fluctuation_rnames_1)

terrestrial_fluctuation_study_1 <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                          Terrestrial_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                              row.names = terrestrial_fluctuation_rnames_1)

Terrestrial_Fluctuation_Model_Estimates_Reorder_1 <- Terrestrial_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

terrestrial_fluctuation_table_1 <- data.frame(estimate = Terrestrial_Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                              lowerCL = Terrestrial_Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                              upperCL = Terrestrial_Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                              K = terrestrial_fluctuation_k_1[,1], 
                                              group_no = terrestrial_fluctuation_group_no_1[,1], 
                                              row.names = terrestrial_fluctuation_rnames_1)
terrestrial_fluctuation_table_1$name <- row.names(terrestrial_fluctuation_table_1)

terrestrial_fluctuation_raw_mean_1 <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InRR_Transformed"))), 
                                        unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InRR_Transformed"))))

terrestrial_fluctuation_raw_name_1 <- c(replicate(35, "Sinusoidal (Sine Curve)"), 
                                        replicate(42, "Alternating"))

terrestrial_fluctuation_raw_df_1 <- data.frame("Model" = terrestrial_fluctuation_raw_name_1, 
                                               "Effect" = terrestrial_fluctuation_raw_mean_1)

# Graph code - Part 1

Terrestrial_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_terrestrial_fluctuation_1 <- terrestrial_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order_1)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = terrestrial_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order_1)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)), xmin = min(terrestrial_fluctuation_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)), xmin = max(terrestrial_fluctuation_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                      vjust = c(-2.7, -0.8))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                     scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                                     annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(terrestrial_fluctuation_table_1["Alternating", "K"], 
                                                                   terrestrial_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                 c(terrestrial_fluctuation_table_1["Alternating", "group_no"], 
                                                                   terrestrial_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation_1 #(400x240)

# Preparing Graph - Part 2

terrestrial_fluctuation_rnames_2 <- c("Stepwise")

terrestrial_fluctuation_k_2 <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                          row.names = terrestrial_fluctuation_rnames_2)

terrestrial_fluctuation_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                 row.names = terrestrial_fluctuation_rnames_2)

terrestrial_fluctuation_study_2 <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                              row.names = terrestrial_fluctuation_rnames_2)

Terrestrial_Fluctuation_Model_Estimates_Reorder_2 <- Terrestrial_Fluctuation_Model_Estimates[c("Stepwise"), ]

terrestrial_fluctuation_table_2 <- data.frame(estimate = Terrestrial_Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                              lowerCL = Terrestrial_Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                              upperCL = Terrestrial_Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                              K = terrestrial_fluctuation_k_2[,1], 
                                              group_no = terrestrial_fluctuation_group_no_2[,1], 
                                              row.names = terrestrial_fluctuation_rnames_2)
terrestrial_fluctuation_table_2$name <- row.names(terrestrial_fluctuation_table_2)

terrestrial_fluctuation_raw_mean_2 <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InRR_Transformed"))))

terrestrial_fluctuation_raw_name_2 <- c(replicate(46, "Stepwise"))

terrestrial_fluctuation_raw_df_2 <- data.frame("Model" = terrestrial_fluctuation_raw_name_2, 
                                               "Effect" = terrestrial_fluctuation_raw_mean_2)

# Graph code - Part 2

Terrestrial_Fluctuation_Order_2 <- c("Stepwise")

density_terrestrial_fluctuation_2 <- terrestrial_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order_2)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = terrestrial_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order_2)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)), xmin = min(terrestrial_fluctuation_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)), xmin = max(terrestrial_fluctuation_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                                     annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(terrestrial_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                 c(terrestrial_fluctuation_table_2["Stepwise", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation_2 #(400x160)

##### Terrestrial Subset Model - Trait Meta-Regression #####
Terrestrial_Trait_Exploration <- Terrestrial_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Exploration) <- Terrestrial_Trait_Exploration$Trait_Category

Terrestrial_Trait_Data <- Terrestrial_Subset_Data %>% filter(Trait_Category != "Behavioural" &
                                                             Trait_Category != "Gene Expression")

Terrestrial_Trait_Species_Count <- Terrestrial_Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Species_Count) <- Terrestrial_Trait_Species_Count$Trait_Category

Terrestrial_Trait_Study_Count <- Terrestrial_Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Study_Count) <- Terrestrial_Trait_Study_Count$Trait_Category

Terrestrial_Trait_Species <- Terrestrial_Trait_Data %>% select("phylo") %>% unique()

Terrestrial_Trait_A_cor <- as.data.frame(A_cor)
Terrestrial_Trait_A_cor <- Terrestrial_Trait_A_cor[c(Terrestrial_Trait_Species$phylo), c(Terrestrial_Trait_Species$phylo)]
Terrestrial_Trait_A_cor <- as.matrix(Terrestrial_Trait_A_cor)

Terrestrial_Trait_VCV <- make_VCV_matrix(Terrestrial_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Terrestrial_Trait_VCV, test = "t", dfs = "contain",
                                               mods = ~ Trait_Category - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Terrestrial_Trait_A_cor), data = Terrestrial_Trait_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Trait_Model.rds")
  } else {
    Terrestrial_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Trait_Model.rds")})

Terrestrial_Trait_Model_rob <- robust(Terrestrial_Trait_Model, cluster = Terrestrial_Trait_Data$Study_ID, adjust = TRUE)

Terrestrial_Trait_Model_Estimates <- data.frame(Category = substr(row.names(Terrestrial_Trait_Model$b), 15, 100),
                                                estimate = Terrestrial_Trait_Model$b, ci.lb = Terrestrial_Trait_Model$ci.lb, 
                                                ci.ub = Terrestrial_Trait_Model$ci.ub)
rownames(Terrestrial_Trait_Model_Estimates) <- Terrestrial_Trait_Model_Estimates$Category
Terrestrial_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Trait_Model), 2))

# Preparing Graph - Combined

terrestrial_trait_rnames <- c("Biochemical Assay", "Life-history Traits", "Morphology", "Physiological")

terrestrial_trait_k <- data.frame("k" = c(Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                          Terrestrial_Trait_Exploration["Life-History Traits", "Freq"], 
                                          Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                          Terrestrial_Trait_Exploration["Physiological", "Freq"]), 
                                  row.names = terrestrial_trait_rnames)

terrestrial_trait_group_no <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Life-History Traits", "Freq"],
                                                       Terrestrial_Trait_Species_Count["Morphology", "Freq"],
                                                       Terrestrial_Trait_Species_Count["Physiological", "Freq"]), 
                                         row.names = terrestrial_trait_rnames)

terrestrial_trait_study <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Life-History Traits", "Freq"],
                                                  Terrestrial_Trait_Study_Count["Morphology", "Freq"],
                                                  Terrestrial_Trait_Study_Count["Physiological", "Freq"]), 
                                      row.names = terrestrial_trait_rnames)

terrestrial_trait_table <- data.frame(estimate = Terrestrial_Trait_Model_Estimates[,"estimate"], 
                                      lowerCL = Terrestrial_Trait_Model_Estimates[,"ci.lb"], 
                                      upperCL = Terrestrial_Trait_Model_Estimates[,"ci.ub"], 
                                      K = terrestrial_trait_k[,1], 
                                      group_no = terrestrial_trait_group_no[,1], 
                                      row.names = terrestrial_trait_rnames)
terrestrial_trait_table$name <- row.names(terrestrial_trait_table)

terrestrial_trait_raw_mean <- c(unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                select("InRR_Transformed"))),
                                unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InRR_Transformed"))))

terrestrial_trait_raw_name <- c(replicate(28, "Biochemical Assay"), 
                                replicate(53, "Life-history Traits"), 
                                replicate(46, "Morphology"),
                                replicate(18, "Physiological"))

terrestrial_trait_raw_df <- data.frame("Model" = terrestrial_trait_raw_name, 
                                       "Effect" = terrestrial_trait_raw_mean)

# Graph code - Combined

Terrestrial_Trait_Order <- c("Physiological", "Morphology", "Life-history Traits", "Biochemical Assay")

density_terrestrial_trait <- terrestrial_trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = terrestrial_trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)), xmin = min(terrestrial_trait_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)), xmin = max(terrestrial_trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (PRRD)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                              vjust = c(-2.7, -2.7, -0.8, -0.8))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                             coord_cartesian(xlim = c(-0.5, 0.5)) +
                             annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_trait_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(terrestrial_trait_table["Physiological", "K"], 
                                                           terrestrial_trait_table["Morphology", "K"], 
                                                           terrestrial_trait_table["Life-history Traits", "K"],
                                                           terrestrial_trait_table["Biochemical Assay", "K"]), "~","(", 
                                                         c(terrestrial_trait_table["Physiological", "group_no"], 
                                                           terrestrial_trait_table["Morphology", "group_no"], 
                                                           terrestrial_trait_table["Life-history Traits", "group_no"],
                                                           terrestrial_trait_table["Biochemical Assay", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(terrestrial_trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait #(400x400)

# Preparing Graph - Part 1

terrestrial_trait_rnames_1 <- c("Biochemical Assay", "Life-history Traits")

terrestrial_trait_k_1 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                            Terrestrial_Trait_Exploration["Life-History Traits", "Freq"]), 
                                    row.names = terrestrial_trait_rnames_1)

terrestrial_trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Life-History Traits", "Freq"]), 
                                           row.names = terrestrial_trait_rnames_1)

terrestrial_trait_study_1 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Life-History Traits", "Freq"]), 
                                        row.names = terrestrial_trait_rnames_1)

Terrestrial_Trait_Model_Estimates_Reorder_1 <- Terrestrial_Trait_Model_Estimates[c("Biochemical Assay", "Life-History Traits"), ]

terrestrial_trait_table_1 <- data.frame(estimate = Terrestrial_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                      lowerCL = Terrestrial_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                      upperCL = Terrestrial_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                      K = terrestrial_trait_k_1[,1], 
                                      group_no = terrestrial_trait_group_no_1[,1], 
                                      row.names = terrestrial_trait_rnames_1)
terrestrial_trait_table_1$name <- row.names(terrestrial_trait_table_1)

terrestrial_trait_raw_mean_1 <- c(unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InRR_Transformed"))), 
                                  unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InRR_Transformed"))))

terrestrial_trait_raw_name_1 <- c(replicate(28, "Biochemical Assay"), 
                                  replicate(53, "Life-history Traits"))

terrestrial_trait_raw_df_1 <- data.frame("Model" = terrestrial_trait_raw_name_1, 
                                         "Effect" = terrestrial_trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Trait_Order_1 <- c("Life-history Traits", "Biochemical Assay")

density_terrestrial_trait_1 <- terrestrial_trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = terrestrial_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)), xmin = min(terrestrial_trait_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)), xmin = max(terrestrial_trait_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-0.8, -0.8))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                               scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_trait_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(terrestrial_trait_table_1["Life-history Traits", "K"],
                                                             terrestrial_trait_table_1["Biochemical Assay", "K"]), "~","(", 
                                                           c(terrestrial_trait_table_1["Life-history Traits", "group_no"],
                                                             terrestrial_trait_table_1["Biochemical Assay", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(terrestrial_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_1 #(400x240)

# Preparing Graph - Part 2

terrestrial_trait_rnames_2 <- c("Morphology", "Physiological")

terrestrial_trait_k_2 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                            Terrestrial_Trait_Exploration["Physiological", "Freq"]), 
                                    row.names = terrestrial_trait_rnames_2)

terrestrial_trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Morphology", "Freq"],
                                                         Terrestrial_Trait_Species_Count["Physiological", "Freq"]), 
                                           row.names = terrestrial_trait_rnames_2)

terrestrial_trait_study_2 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Morphology", "Freq"],
                                                    Terrestrial_Trait_Study_Count["Physiological", "Freq"]), 
                                        row.names = terrestrial_trait_rnames_2)

Terrestrial_Trait_Model_Estimates_Reorder_2 <- Terrestrial_Trait_Model_Estimates[c("Morphology", "Physiological"), ]

terrestrial_trait_table_2 <- data.frame(estimate = Terrestrial_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Terrestrial_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Terrestrial_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                        K = terrestrial_trait_k_2[,1], 
                                        group_no = terrestrial_trait_group_no_2[,1], 
                                        row.names = terrestrial_trait_rnames_2)
terrestrial_trait_table_2$name <- row.names(terrestrial_trait_table_2)

terrestrial_trait_raw_mean_2 <- c(unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                select("InRR_Transformed"))),
                                  unlist(unname(Terrestrial_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InRR_Transformed"))))

terrestrial_trait_raw_name_2 <- c(replicate(46, "Morphology"),
                                  replicate(18, "Physiological"))

terrestrial_trait_raw_df_2 <- data.frame("Model" = terrestrial_trait_raw_name_2, 
                                         "Effect" = terrestrial_trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Trait_Order_2 <- c("Physiological", "Morphology")

density_terrestrial_trait_2 <- terrestrial_trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = terrestrial_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)), xmin = min(terrestrial_trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)), xmin = max(terrestrial_trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_trait_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(terrestrial_trait_table_2["Physiological", "K"], 
                                                             terrestrial_trait_table_2["Morphology", "K"]), "~","(", 
                                                           c(terrestrial_trait_table_2["Physiological", "group_no"], 
                                                             terrestrial_trait_table_2["Morphology", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Terrestrial_Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(terrestrial_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_2 #(400x240)

##### Terrestrial Subset Model - Plasticity Mechanism Meta-Regression #####
Terrestrial_Plasticity_Exploration <- Terrestrial_Subset_Data %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Exploration) <- Terrestrial_Plasticity_Exploration$Plasticity_Mechanism

Terrestrial_Plasticity_Species_Count <- Terrestrial_Subset_Data %>% select("Scientific_Name", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Species_Count) <- Terrestrial_Plasticity_Species_Count$Plasticity_Mechanism

Terrestrial_Plasticity_Study_Count <- Terrestrial_Subset_Data %>% select("Study_ID", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Study_Count) <- Terrestrial_Plasticity_Study_Count$Plasticity_Mechanism

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Plasticity_Model <- metafor::rma.mv(InRR_Transformed, V = Terrestrial_VCV, test = "t", dfs = "contain",
                                                    mods = ~ Plasticity_Mechanism - 1,
                                                    random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                  ~1|Shared_Animal_Number, ~1|Measurement), 
                                                    R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                                    control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Plasticity_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Plasticity_Model.rds")
  } else {
    Terrestrial_Plasticity_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Plasticity_Model.rds")})

Terrestrial_Plasticity_Model_rob <- robust(Terrestrial_Plasticity_Model, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Plasticity_Model_Estimates <- data.frame(Plasticity_Mechanism = substr(row.names(Terrestrial_Plasticity_Model$b), 21, 100),
                                                     estimate = Terrestrial_Plasticity_Model$b, ci.lb = Terrestrial_Plasticity_Model$ci.lb, 
                                                     ci.ub = Terrestrial_Plasticity_Model$ci.ub)
rownames(Terrestrial_Plasticity_Model_Estimates) <- Terrestrial_Plasticity_Model_Estimates$Plasticity_Mechanism
Terrestrial_Plasticity_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Plasticity_Model), 2))

# Preparing Graph - Combined

terrestrial_plasticity_rnames <- c("Acclimation", "Development")

terrestrial_plasticity_k <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"], 
                                               Terrestrial_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                       row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_group_no <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                            Terrestrial_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                              row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_study <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                       Terrestrial_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                           row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_table <- data.frame(estimate = Terrestrial_Plasticity_Model_Estimates[,"estimate"], 
                                           lowerCL = Terrestrial_Plasticity_Model_Estimates[,"ci.lb"], 
                                           upperCL = Terrestrial_Plasticity_Model_Estimates[,"ci.ub"], 
                                           K = terrestrial_plasticity_k[,1], 
                                           group_no = terrestrial_plasticity_group_no[,1], 
                                           row.names = terrestrial_plasticity_rnames)
terrestrial_plasticity_table$name <- row.names(terrestrial_plasticity_table)

terrestrial_plasticity_raw_mean <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                     select("InRR_Transformed"))))

terrestrial_plasticity_raw_name <- c(replicate(37, "Acclimation"), 
                                     replicate(115, "Development"))

terrestrial_plasticity_raw_df <- data.frame("Model" = terrestrial_plasticity_raw_name, 
                                            "Effect" = terrestrial_plasticity_raw_mean)

# Graph code - Combined

Terrestrial_Plasticity_Order <- c("Development", "Acclimation")

density_terrestrial_plasticity <- terrestrial_plasticity_table %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = terrestrial_plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)), xmin = min(terrestrial_plasticity_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)), xmin = max(terrestrial_plasticity_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                  coord_cartesian(xlim = c(-0.5, 0.5)) +
                                  annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_plasticity_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(terrestrial_plasticity_table["Development", "K"], 
                                                                terrestrial_plasticity_table["Acclimation", "K"]), "~","(", 
                                                              c(terrestrial_plasticity_table["Development", "group_no"], 
                                                                terrestrial_plasticity_table["Acclimation", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Terrestrial_Plasticity_Model_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(terrestrial_plasticity_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity #(400x240)

# Preparing Graph - Part 1

terrestrial_plasticity_rnames_1 <- c("Acclimation")

terrestrial_plasticity_k_1 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"]), 
                                         row.names = terrestrial_plasticity_rnames_1)

terrestrial_plasticity_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                row.names = terrestrial_plasticity_rnames_1)

terrestrial_plasticity_study_1 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                             row.names = terrestrial_plasticity_rnames_1)

Terrestrial_Plasticity_Model_Estimates_Reorder_1 <- Terrestrial_Plasticity_Model_Estimates[c("Acclimation"), ]

terrestrial_plasticity_table_1 <- data.frame(estimate = Terrestrial_Plasticity_Model_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Terrestrial_Plasticity_Model_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Terrestrial_Plasticity_Model_Estimates_Reorder_1[,"ci.ub"], 
                                             K = terrestrial_plasticity_k_1[,1], 
                                             group_no = terrestrial_plasticity_group_no_1[,1], 
                                             row.names = terrestrial_plasticity_rnames_1)
terrestrial_plasticity_table_1$name <- row.names(terrestrial_plasticity_table_1)

terrestrial_plasticity_raw_mean_1 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                     select("InRR_Transformed"))))

terrestrial_plasticity_raw_name_1 <- c(replicate(37, "Acclimation"))

terrestrial_plasticity_raw_df_1 <- data.frame("Model" = terrestrial_plasticity_raw_name_1, 
                                              "Effect" = terrestrial_plasticity_raw_mean_1)

# Graph code - Part 1

Terrestrial_Plasticity_Order_1 <- c("Acclimation")

density_terrestrial_plasticity_1 <- terrestrial_plasticity_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = terrestrial_plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)), xmin = min(terrestrial_plasticity_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)), xmin = max(terrestrial_plasticity_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_plasticity_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(terrestrial_plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                c(terrestrial_plasticity_table_1["Acclimation", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(terrestrial_plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity_1 #(400x160)

# Preparing Graph - Part 2

terrestrial_plasticity_rnames_2 <- c("Development")

terrestrial_plasticity_k_2 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                         row.names = terrestrial_plasticity_rnames_2)

terrestrial_plasticity_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                                row.names = terrestrial_plasticity_rnames_2)

terrestrial_plasticity_study_2 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                             row.names = terrestrial_plasticity_rnames_2)

Terrestrial_Plasticity_Model_Estimates_Reorder_2 <- Terrestrial_Plasticity_Model_Estimates[c("Development"), ]

terrestrial_plasticity_table_2 <- data.frame(estimate = Terrestrial_Plasticity_Model_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Terrestrial_Plasticity_Model_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Terrestrial_Plasticity_Model_Estimates_Reorder_2[,"ci.ub"], 
                                             K = terrestrial_plasticity_k_2[,1], 
                                             group_no = terrestrial_plasticity_group_no_2[,1], 
                                             row.names = terrestrial_plasticity_rnames_2)
terrestrial_plasticity_table_2$name <- row.names(terrestrial_plasticity_table_2)

terrestrial_plasticity_raw_mean_2 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                     select("InRR_Transformed"))))

terrestrial_plasticity_raw_name_2 <- c(replicate(115, "Development"))

terrestrial_plasticity_raw_df_2 <- data.frame("Model" = terrestrial_plasticity_raw_name_2, 
                                              "Effect" = terrestrial_plasticity_raw_mean_2)

# Graph code - Part 2

Terrestrial_Plasticity_Order_2 <- c("Development")

density_terrestrial_plasticity_2 <- terrestrial_plasticity_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = terrestrial_plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)), xmin = min(terrestrial_plasticity_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)), xmin = max(terrestrial_plasticity_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_plasticity_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(terrestrial_plasticity_table_2["Development", "K"]), "~","(", 
                                                                c(terrestrial_plasticity_table_2["Development", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(terrestrial_plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity_2 #(400x160)

##### Terrestrial Subset Model - Specific Trait Meta-Regression #####
Terrestrial_Specific_Trait_Exploration <- Terrestrial_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Terrestrial_Specific_Trait_Exploration <- Terrestrial_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Terrestrial_Specific_Trait_Exploration) <- Terrestrial_Specific_Trait_Exploration$Measurement

Terrestrial_Specific_Trait_Data <- Terrestrial_Subset_Data %>% filter(Measurement == "Development Time"| 
                                                                      Measurement == "Length"|
                                                                      Measurement == "Mass")

Terrestrial_Specific_Trait_Species_Count <- Terrestrial_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Terrestrial_Specific_Trait_Species_Count) <- Terrestrial_Specific_Trait_Species_Count$Measurement

Terrestrial_Specific_Trait_Study_Count <- Terrestrial_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Terrestrial_Specific_Trait_Study_Count) <- Terrestrial_Specific_Trait_Study_Count$Measurement

Terrestrial_Specific_Trait_Species <- Terrestrial_Specific_Trait_Data %>% select("phylo") %>% unique()

Terrestrial_Specific_Trait_A_cor <- as.data.frame(A_cor)
Terrestrial_Specific_Trait_A_cor <- Terrestrial_Specific_Trait_A_cor[c(Terrestrial_Specific_Trait_Species$phylo), c(Terrestrial_Specific_Trait_Species$phylo)]
Terrestrial_Specific_Trait_A_cor <- as.matrix(Terrestrial_Specific_Trait_A_cor)

Terrestrial_Specific_Trait_VCV <- make_VCV_matrix(Terrestrial_Specific_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Terrestrial_Specific_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Terrestrial_Specific_Trait_VCV, test = "t", dfs = "contain",
                                                        mods = ~ Measurement - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number), 
                                                        R = list(phylo=Terrestrial_Specific_Trait_A_cor), data = Terrestrial_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Specific_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Specific_Trait_Model.rds")
  } else {
    Terrestrial_Specific_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Terrestrial_Specific_Trait_Model.rds")})

Terrestrial_Specific_Trait_Model_rob <- robust(Terrestrial_Specific_Trait_Model, cluster = Terrestrial_Specific_Trait_Data$Study_ID, adjust = TRUE)

Terrestrial_Specific_Trait_Model_Estimates <- data.frame(Trait = substr(row.names(Terrestrial_Specific_Trait_Model$b), 12, 100),
                                                         estimate = Terrestrial_Specific_Trait_Model$b, ci.lb = Terrestrial_Specific_Trait_Model$ci.lb, 
                                                         ci.ub = Terrestrial_Specific_Trait_Model$ci.ub)
rownames(Terrestrial_Specific_Trait_Model_Estimates) <- Terrestrial_Specific_Trait_Model_Estimates$Trait
Terrestrial_Specific_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Specific_Trait_Model), 2))

# Preparing Graph - Combined

terrestrial_specific_trait_rnames <- c("Development Time", "Length", "Mass")

terrestrial_specific_trait_k <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Length", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Mass", "Freq"]), 
                                           row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_group_no <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Length", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Mass", "Freq"]), 
                                                  row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_study <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Length", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Mass", "Freq"]), 
                                               row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_table <- data.frame(estimate = Terrestrial_Specific_Trait_Model_Estimates[,"estimate"], 
                                               lowerCL = Terrestrial_Specific_Trait_Model_Estimates[,"ci.lb"], 
                                               upperCL = Terrestrial_Specific_Trait_Model_Estimates[,"ci.ub"], 
                                               K = terrestrial_specific_trait_k[,1], 
                                               group_no = terrestrial_specific_trait_group_no[,1], 
                                               row.names = terrestrial_specific_trait_rnames)
terrestrial_specific_trait_table$name <- row.names(terrestrial_specific_trait_table)

terrestrial_specific_trait_raw_mean <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                         select("InRR_Transformed"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                         select("InRR_Transformed"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                         select("InRR_Transformed"))))

terrestrial_specific_trait_raw_name <- c(replicate(40, "Development Time"), 
                                         replicate(11, "Length"), 
                                         replicate(23, "Mass"))

terrestrial_specific_trait_raw_df <- data.frame("Model" = terrestrial_specific_trait_raw_name, 
                                                "Effect" = terrestrial_specific_trait_raw_mean)

# Graph code - Combined

Terrestrial_Specific_Trait_Order <- c("Mass", "Length", "Development Time")

density_terrestrial_specific_trait <- terrestrial_specific_trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = terrestrial_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                      geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)), xmin = min(terrestrial_specific_trait_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                                     size = 1) +
                                      geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)), xmin = max(terrestrial_specific_trait_raw_df$Effect)+0.01, xmax = 1.5, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                          size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                            vjust = c(-2.7, -2.7, -0.8))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                      coord_cartesian(xlim = c(-0.5, 0.5)) +
                                      annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_specific_trait_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(terrestrial_specific_trait_table["Mass", "K"],
                                                                    terrestrial_specific_trait_table["Length", "K"],
                                                                    terrestrial_specific_trait_table["Development Time", "K"]), "~","(", 
                                                                  c(terrestrial_specific_trait_table["Mass", "group_no"],
                                                                    terrestrial_specific_trait_table["Length", "group_no"], 
                                                                    terrestrial_specific_trait_table["Development Time", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.4, y = (seq(1, dim(terrestrial_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait #(400x320)

# Preparing Graph - Part 1

terrestrial_specific_trait_rnames_1 <- c("Development Time", "Length")

terrestrial_specific_trait_k_1 <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Length", "Freq"]), 
                                             row.names = terrestrial_specific_trait_rnames_1)

terrestrial_specific_trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Length", "Freq"]), 
                                                    row.names = terrestrial_specific_trait_rnames_1)

terrestrial_specific_trait_study_1 <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Length", "Freq"]), 
                                                 row.names = terrestrial_specific_trait_rnames_1)

Terrestrial_Specific_Trait_Model_Estimates_Reorder_1 <- Terrestrial_Specific_Trait_Model_Estimates[c("Development Time", "Length"), ]

terrestrial_specific_trait_table_1 <- data.frame(estimate = Terrestrial_Specific_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                                 lowerCL = Terrestrial_Specific_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                                 upperCL = Terrestrial_Specific_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                                 K = terrestrial_specific_trait_k_1[,1], 
                                                 group_no = terrestrial_specific_trait_group_no_1[,1], 
                                                 row.names = terrestrial_specific_trait_rnames_1)
terrestrial_specific_trait_table_1$name <- row.names(terrestrial_specific_trait_table_1)

terrestrial_specific_trait_raw_mean_1 <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                         select("InRR_Transformed"))), 
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                         select("InRR_Transformed"))))

terrestrial_specific_trait_raw_name_1 <- c(replicate(40, "Development Time"), 
                                           replicate(11, "Length"))

terrestrial_specific_trait_raw_df_1 <- data.frame("Model" = terrestrial_specific_trait_raw_name_1, 
                                                  "Effect" = terrestrial_specific_trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Specific_Trait_Order_1 <- c("Length", "Development Time")

density_terrestrial_specific_trait_1 <- terrestrial_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = terrestrial_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)), xmin = min(terrestrial_specific_trait_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)), xmin = max(terrestrial_specific_trait_raw_df_1$Effect)+0.01, xmax = 1.5, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                         vjust = c(-2.7, -0.8))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                        scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                        coord_cartesian(xlim = c(-0.5, 0.5)) +
                                        annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(terrestrial_specific_trait_table_1["Length", "K"],
                                                                      terrestrial_specific_trait_table_1["Development Time", "K"]), "~","(", 
                                                                    c(terrestrial_specific_trait_table_1["Length", "group_no"], 
                                                                      terrestrial_specific_trait_table_1["Development Time", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.4, y = (seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait_1 #(400x240)

# Preparing Graph - Part 2

terrestrial_specific_trait_rnames_2 <- c("Mass")

terrestrial_specific_trait_k_2 <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Mass", "Freq"]), 
                                             row.names = terrestrial_specific_trait_rnames_2)

terrestrial_specific_trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Mass", "Freq"]), 
                                                    row.names = terrestrial_specific_trait_rnames_2)

terrestrial_specific_trait_study_2 <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Mass", "Freq"]), 
                                                 row.names = terrestrial_specific_trait_rnames_2)

Terrestrial_Specific_Trait_Model_Estimates_Reorder_2 <- Terrestrial_Specific_Trait_Model_Estimates[c("Mass"), ]

terrestrial_specific_trait_table_2 <- data.frame(estimate = Terrestrial_Specific_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                               lowerCL = Terrestrial_Specific_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                               upperCL = Terrestrial_Specific_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                               K = terrestrial_specific_trait_k_2[,1], 
                                               group_no = terrestrial_specific_trait_group_no_2[,1], 
                                               row.names = terrestrial_specific_trait_rnames_2)
terrestrial_specific_trait_table_2$name <- row.names(terrestrial_specific_trait_table_2)

terrestrial_specific_trait_raw_mean_2 <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                         select("InRR_Transformed"))))

terrestrial_specific_trait_raw_name_2 <- c(replicate(23, "Mass"))

terrestrial_specific_trait_raw_df_2 <- data.frame("Model" = terrestrial_specific_trait_raw_name_2, 
                                                  "Effect" = terrestrial_specific_trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Specific_Trait_Order_2 <- c("Mass")

density_terrestrial_specific_trait_2 <- terrestrial_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = terrestrial_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.02, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)), xmin = min(terrestrial_specific_trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)), xmin = max(terrestrial_specific_trait_raw_df_2$Effect)+0.01, xmax = 1.5, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                        coord_cartesian(xlim = c(-0.5, 0.5)) +
                                        annotate('text',  x = 0.5, y = (seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(terrestrial_specific_trait_table["Mass", "K"]), "~","(", 
                                                                    c(terrestrial_specific_trait_table["Mass", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.4, y = (seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait_2 #(400x160)

##### Acclimation Subset Model #####
Acclimation_Subset_Data <- Individual_Subset_Data %>% filter(Plasticity_Mechanism == "Acclimation")
Acclimation_Species <- Acclimation_Subset_Data %>% select("phylo") %>% unique()

Acclimation_A_cor <- as.data.frame(A_cor)
Acclimation_A_cor <- Acclimation_A_cor[c(Acclimation_Species$phylo), c(Acclimation_Species$phylo)]
Acclimation_A_cor <- as.matrix(Acclimation_A_cor)

Acclimation_VCV <- make_VCV_matrix(Acclimation_Subset_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Acclimation_VCV, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                         R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Model.rds")
  } else {
    Acclimation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Model.rds")})

Acclimation_Model_rob <- robust(Acclimation_Model, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Model_Estimates <- data.frame(estimate = Acclimation_Model$b, ci.lb = Acclimation_Model$ci.lb, ci.ub = Acclimation_Model$ci.ub)
Acclimation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Model), 2))

#### Acclimation Subset Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_VCV, test = "t", dfs = "contain",
                                                   mods = ~ Fluctuation_Magnitude - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Amplitude_Model.rds")
  } else {
    Acclimation_Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Amplitude_Model.rds")})

Acclimation_Amplitude_Model_rob <- robust(Acclimation_Amplitude_Model, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Amplitude_Model_Estimates <- data.frame(estimate = Acclimation_Amplitude_Model$b, ci.lb = Acclimation_Amplitude_Model$ci.lb, 
                                                    ci.ub = Acclimation_Amplitude_Model$ci.ub)
Acclimation_Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Amplitude_Model), 2))

# Graph Preparing

Acclimation_Plot_Data <- Acclimation_Subset_Data
Acclimation_Plot_Data <- Acclimation_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                       ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                       ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Acclimation_Amplitude_Plot <- ggplot(Acclimation_Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                              geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                             size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                              labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                                   size = "Sample Size", title = "Acclimation Treatments") +
                              theme_bw() +
                              theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                              theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                              theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                              #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                              geom_hline(yintercept = Acclimation_Model_Estimates$estimate, lty = 2) + 
                              geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                              stat_poly_eq(formula = y ~ x, 
                                           aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                           parse = TRUE) +
                              coord_cartesian(xlim = c(0, 25), 
                                              ylim = c(-0.25, 0.25))

Acclimation_Amplitude_Plot #(330x330)

#### Acclimation Subset Model - Exposure Time Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Exposure_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_VCV, test = "t", dfs = "contain",
                                                  mods = ~ Acclimation_Exposure_Time - 1,
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                                  R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Exposure_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Exposure_Model.rds")
  } else {
    Acclimation_Exposure_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Exposure_Model.rds")})

Acclimation_Exposure_Model_rob <- robust(Acclimation_Exposure_Model, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Exposure_Model_Estimates <- data.frame(estimate = Acclimation_Exposure_Model$b, ci.lb = Acclimation_Exposure_Model$ci.lb, 
                                                   ci.ub = Acclimation_Exposure_Model$ci.ub)
Acclimation_Exposure_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Exposure_Model), 2))

# Graph Preparing

Acclimation_Exposure_Plot_Data <- Acclimation_Subset_Data
Acclimation_Exposure_Plot_Data <- Acclimation_Exposure_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                                         ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                                         ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Acclimation_Exposure_Amplitude_Plot <- ggplot(Acclimation_Exposure_Plot_Data, aes(x = Acclimation_Exposure_Time, y = InRR_Transformed)) + 
                                       geom_point(aes(x = Acclimation_Exposure_Time, y = InRR_Transformed, 
                                                      size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                                      shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                       labs(x = "Exposure Time (Days)", y = "Effect Size (PRRD)", 
                                            size = "Sample Size", title = "Acclimation Treatments") +
                                       theme_bw() +
                                       theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                       theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                       theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                       geom_hline(yintercept = Acclimation_Model_Estimates$estimate, lty = 2) + 
                                       geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                       stat_poly_eq(formula = y ~ x, 
                                                    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                                        parse = TRUE) +
                                       coord_cartesian(xlim = c(0, 80), 
                                                       ylim = c(c(-0.25, 0.25)))

Acclimation_Exposure_Amplitude_Plot #(400x400)

#### Acclimation Subset Model - Fluctuation Frequency Meta-Regression ####
Acclimation_Frequency_Data <- Acclimation_Subset_Data %>% filter(!is.na(Number_Of_Fluctuations))
Acclimation_Frequency_Species <- Acclimation_Frequency_Data %>% select("phylo") %>% unique()

Acclimation_Frequency_A_cor <- as.data.frame(A_cor)
Acclimation_Frequency_A_cor <- Acclimation_Frequency_A_cor[c(Acclimation_Frequency_Species$phylo), c(Acclimation_Frequency_Species$phylo)]
Acclimation_Frequency_A_cor <- as.matrix(Acclimation_Frequency_A_cor)

Acclimation_Frequency_VCV <- make_VCV_matrix(Acclimation_Frequency_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Frequency_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_Frequency_VCV, test = "t", dfs = "contain",
                                                   mods = ~ Number_Of_Fluctuations - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Acclimation_Frequency_A_cor), data = Acclimation_Frequency_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Frequency_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Frequency_Model.rds")
  } else {
    Acclimation_Frequency_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Frequency_Model.rds")})

Acclimation_Frequency_Model_rob <- robust(Acclimation_Frequency_Model, cluster = Acclimation_Frequency_Data$Study_ID, adjust = TRUE)

Acclimation_Frequency_Model_Estimates <- data.frame(estimate = Acclimation_Frequency_Model$b, ci.lb = Acclimation_Frequency_Model$ci.lb, 
                                                    ci.ub = Acclimation_Frequency_Model$ci.ub)
Acclimation_Frequency_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Frequency_Model), 2))

# Graph Preparing

Acclimation_Frequency_Plot_Data <- Acclimation_Frequency_Data
Acclimation_Frequency_Plot_Data <- Acclimation_Frequency_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                                           ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                                           ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Acclimation_Frequency_Plot <- ggplot(Acclimation_Frequency_Plot_Data, aes(x = Number_Of_Fluctuations, y = InRR_Transformed)) + 
                              geom_point(aes(x = Number_Of_Fluctuations, y = InRR_Transformed, 
                                             size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5) + 
                              labs(x = "Number of Fluctuations", y = "Effect Size (PRRD)", 
                                   size = "Sample Size", title = "Acclimation Treatments") +
                              theme_bw() +
                              theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                              theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                              theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                              theme(legend.position = "bottom", legend.direction = "horizontal") + 
                              geom_hline(yintercept = Acclimation_Model_Estimates$estimate, lty = 2) + 
                              geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                              stat_poly_eq(formula = y ~ x, 
                                           aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                               parse = TRUE) +
                              coord_cartesian(xlim = c(0, 80),
                                              ylim = c(-0.25, 0.25))

Acclimation_Frequency_Plot #(400x400)

#### Acclimation Subset Model - Type of Fluctuation Meta-Regression ####
Acclimation_Fluctuation_Data <- Acclimation_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Acclimation_Fluctuation_Exploration <- Acclimation_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Exploration) <- Acclimation_Fluctuation_Exploration$Fluctuation_Category

Acclimation_Fluctuation_Data <- Acclimation_Fluctuation_Data %>% filter(Fluctuation_Category != "Alternating" &
                                                                        Fluctuation_Category != "Stochastic")

Acclimation_Fluctuation_Species_Count <- Acclimation_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Species_Count) <- Acclimation_Fluctuation_Species_Count$Fluctuation_Category

Acclimation_Fluctuation_Study_Count <- Acclimation_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Study_Count) <- Acclimation_Fluctuation_Study_Count$Fluctuation_Category

Acclimation_Fluctuation_Species <- Acclimation_Fluctuation_Data %>% select("phylo") %>% unique()

Acclimation_Fluctuation_A_cor <- as.data.frame(A_cor)
Acclimation_Fluctuation_A_cor <- Acclimation_Fluctuation_A_cor[c(Acclimation_Fluctuation_Species$phylo), c(Acclimation_Fluctuation_Species$phylo)]
Acclimation_Fluctuation_A_cor <- as.matrix(Acclimation_Fluctuation_A_cor)

Acclimation_Fluctuation_VCV <- make_VCV_matrix(Acclimation_Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_Fluctuation_VCV, test = "t", dfs = "contain",
                                                     mods = ~ Fluctuation_Category - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Acclimation_Fluctuation_A_cor), data = Acclimation_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Fluctuation_Model.rds")
  } else {
    Acclimation_Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Fluctuation_Model.rds")})

Acclimation_Fluctuation_Model_rob <- robust(Acclimation_Fluctuation_Model, cluster = Acclimation_Fluctuation_Data$Study_ID, adjust = TRUE)

Acclimation_Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Acclimation_Fluctuation_Model$b), 21, 100),
                                                      estimate = Acclimation_Fluctuation_Model$b, ci.lb = Acclimation_Fluctuation_Model$ci.lb, 
                                                      ci.ub = Acclimation_Fluctuation_Model$ci.ub)
rownames(Acclimation_Fluctuation_Model_Estimates) <- Acclimation_Fluctuation_Model_Estimates$Category
Acclimation_Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Fluctuation_Model), 2))

# Preparing Graph - Combined

acclimation_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Stepwise")

acclimation_fluctuation_k <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                Acclimation_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                        row.names = acclimation_fluctuation_rnames)

acclimation_fluctuation_group_no <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                             Acclimation_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                               row.names = acclimation_fluctuation_rnames)

acclimation_fluctuation_study <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                        Acclimation_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                            row.names = acclimation_fluctuation_rnames)

Acclimation_Fluctuation_Model_Estimates_Reorder <- Acclimation_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Stepwise"), ]

acclimation_fluctuation_table <- data.frame(estimate = Acclimation_Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                            lowerCL = Acclimation_Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                            upperCL = Acclimation_Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                            K = acclimation_fluctuation_k[,1], 
                                            group_no = acclimation_fluctuation_group_no[,1], 
                                            row.names = acclimation_fluctuation_rnames)
acclimation_fluctuation_table$name <- row.names(acclimation_fluctuation_table)

acclimation_fluctuation_raw_mean <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InRR_Transformed"))), 
                                      unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InRR_Transformed"))))

acclimation_fluctuation_raw_name <- c(replicate(25, "Sinusoidal (Sine Curve)"), 
                                      replicate(30, "Stepwise"))

acclimation_fluctuation_raw_df <- data.frame("Model" = acclimation_fluctuation_raw_name, 
                                             "Effect" = acclimation_fluctuation_raw_mean)

# Graph code - Combined

Acclimation_Fluctuation_Order <- c("Stepwise", "Sinusoidal (Sine Curve)")

density_acclimation_fluctuation <- acclimation_fluctuation_table %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)), xmin = min(acclimation_fluctuation_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)), xmin = max(acclimation_fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7, -0.8))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                   coord_cartesian(xlim = c(-0.5, 0.5)) +
                                   annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_fluctuation_table)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_fluctuation_table["Stepwise", "K"], 
                                                                 acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                               c(acclimation_fluctuation_table["Stepwise", "group_no"], 
                                                                 acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.4, y = (seq(1, dim(acclimation_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation #(400x240)

# Preparing Graph - Part 1

acclimation_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)")

acclimation_fluctuation_k_1 <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"]), 
                                          row.names = acclimation_fluctuation_rnames_1)

acclimation_fluctuation_group_no_1 <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                                 row.names = acclimation_fluctuation_rnames_1)

acclimation_fluctuation_study_1 <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                              row.names = acclimation_fluctuation_rnames_1)

Acclimation_Fluctuation_Model_Estimates_Reorder_1 <- Acclimation_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)"), ]

acclimation_fluctuation_table_1 <- data.frame(estimate = Acclimation_Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                              lowerCL = Acclimation_Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                              upperCL = Acclimation_Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                              K = acclimation_fluctuation_k_1[,1], 
                                              group_no = acclimation_fluctuation_group_no_1[,1], 
                                              row.names = acclimation_fluctuation_rnames_1)
acclimation_fluctuation_table_1$name <- row.names(acclimation_fluctuation_table_1)

acclimation_fluctuation_raw_mean_1 <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InRR_Transformed"))))

acclimation_fluctuation_raw_name_1 <- c(replicate(25, "Sinusoidal (Sine Curve)"))

acclimation_fluctuation_raw_df_1 <- data.frame("Model" = acclimation_fluctuation_raw_name_1, 
                                               "Effect" = acclimation_fluctuation_raw_mean_1)

# Graph code - Part 1

Acclimation_Fluctuation_Order_1 <- c("Sinusoidal (Sine Curve)")

density_acclimation_fluctuation_1 <- acclimation_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order_1)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = acclimation_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order_1)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)), xmin = min(acclimation_fluctuation_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)), xmin = max(acclimation_fluctuation_raw_df_1$Effect)+0.03, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                     scale_colour_manual(values = c("#2B4E7A")) +
                                     scale_fill_manual(values = c("#2B4E7A")) +
                                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                                     annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_fluctuation_table_1)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(acclimation_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                 c(acclimation_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(acclimation_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation_1 #(400x160)

# Preparing Graph - Part 2

acclimation_fluctuation_rnames_2 <- c("Stepwise")

acclimation_fluctuation_k_2 <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                          row.names = acclimation_fluctuation_rnames_2)

acclimation_fluctuation_group_no_2 <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                row.names = acclimation_fluctuation_rnames_2)

acclimation_fluctuation_study_2 <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                              row.names = acclimation_fluctuation_rnames_2)

Acclimation_Fluctuation_Model_Estimates_Reorder_2 <- Acclimation_Fluctuation_Model_Estimates[c("Stepwise"), ]

acclimation_fluctuation_table_2 <- data.frame(estimate = Acclimation_Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                            lowerCL = Acclimation_Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                            upperCL = Acclimation_Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                            K = acclimation_fluctuation_k_2[,1], 
                                            group_no = acclimation_fluctuation_group_no_2[,1], 
                                            row.names = acclimation_fluctuation_rnames_2)
acclimation_fluctuation_table_2$name <- row.names(acclimation_fluctuation_table_2)

acclimation_fluctuation_raw_mean_2 <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InRR_Transformed"))))

acclimation_fluctuation_raw_name_2 <- c(replicate(30, "Stepwise"))

acclimation_fluctuation_raw_df_2 <- data.frame("Model" = acclimation_fluctuation_raw_name_2, 
                                               "Effect" = acclimation_fluctuation_raw_mean_2)

# Graph code - Part 2

Acclimation_Fluctuation_Order_2 <- c("Stepwise")

density_acclimation_fluctuation_2 <- acclimation_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order_2)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = acclimation_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order_2)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)), xmin = min(acclimation_fluctuation_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)), xmin = max(acclimation_fluctuation_raw_df_2$Effect)+0.03, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                                     annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_fluctuation_table_2)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(acclimation_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                 c(acclimation_fluctuation_table_2["Stepwise", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(acclimation_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation_2 #(400x160)

##### Acclimation Subset Model - Trait Meta-Regression #####
Acclimation_Trait_Exploration <- Acclimation_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Exploration) <- Acclimation_Trait_Exploration$Trait_Category

Acclimation_Trait_Data <- Acclimation_Subset_Data %>% filter(Trait_Category != "Behavioural" &
                                                             Trait_Category != "Gene Expression" &
                                                             Trait_Category != "Life-History Traits")

Acclimation_Trait_Species_Count <- Acclimation_Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Species_Count) <- Acclimation_Trait_Species_Count$Trait_Category

Acclimation_Trait_Study_Count <- Acclimation_Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Study_Count) <- Acclimation_Trait_Study_Count$Trait_Category

Acclimation_Trait_Species <- Acclimation_Trait_Data %>% select("phylo") %>% unique()

Acclimation_Trait_A_cor <- as.data.frame(A_cor)
Acclimation_Trait_A_cor <- Acclimation_Trait_A_cor[c(Acclimation_Trait_Species$phylo), c(Acclimation_Trait_Species$phylo)]
Acclimation_Trait_A_cor <- as.matrix(Acclimation_Trait_A_cor)

Acclimation_Trait_VCV <- make_VCV_matrix(Acclimation_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_Trait_VCV, test = "t", dfs = "contain",
                                               mods = ~ Trait_Category - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Acclimation_Trait_A_cor), data = Acclimation_Trait_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Trait_Model.rds")
  } else {
    Acclimation_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Trait_Model.rds")})

Acclimation_Trait_Model_rob <- robust(Acclimation_Trait_Model, cluster = Acclimation_Trait_Data$Study_ID, adjust = TRUE)

Acclimation_Trait_Model_Estimates <- data.frame(Category = substr(row.names(Acclimation_Trait_Model$b), 15, 100),
                                                estimate = Acclimation_Trait_Model$b, ci.lb = Acclimation_Trait_Model$ci.lb, 
                                                ci.ub = Acclimation_Trait_Model$ci.ub)
rownames(Acclimation_Trait_Model_Estimates) <- Acclimation_Trait_Model_Estimates$Category
Acclimation_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Trait_Model), 2))

# Preparing Graph - Combined

acclimation_trait_rnames <- c("Biochemical Assay", "Physiological")

acclimation_trait_k <- data.frame("k" = c(Acclimation_Trait_Exploration["Biochemical Assay", "Freq"],
                                          Acclimation_Trait_Exploration["Physiological", "Freq"]), 
                                  row.names = acclimation_trait_rnames)

acclimation_trait_group_no <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Acclimation_Trait_Species_Count["Physiological", "Freq"]), 
                                         row.names = acclimation_trait_rnames)

acclimation_trait_study <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Acclimation_Trait_Study_Count["Physiological", "Freq"]), 
                                      row.names = acclimation_trait_rnames)

acclimation_trait_table <- data.frame(estimate = Acclimation_Trait_Model_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Trait_Model_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Trait_Model_Estimates[,"ci.ub"], 
                                      K = acclimation_trait_k[,1], 
                                      group_no = acclimation_trait_group_no[,1], 
                                      row.names = acclimation_trait_rnames)
acclimation_trait_table$name <- row.names(acclimation_trait_table)

acclimation_trait_raw_mean <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InRR_Transformed"))))

acclimation_trait_raw_name <- c(replicate(27, "Biochemical Assay"), 
                                replicate(32, "Physiological"))

acclimation_trait_raw_df <- data.frame("Model" = acclimation_trait_raw_name, 
                                       "Effect" = acclimation_trait_raw_mean)

# Graph code - Combined

Acclimation_Trait_Order <- c("Physiological", "Biochemical Assay")

density_acclimation_trait <- acclimation_trait_table %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = acclimation_trait_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table)[1], 1)), xmin = min(acclimation_trait_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table)[1], 1)), xmin = max(acclimation_trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (PRRD)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                              vjust = c(-2.7, -0.8))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                             scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                             coord_cartesian(xlim = c(-0.5, 0.5)) +
                             annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_trait_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(acclimation_trait_table["Physiological", "K"],
                                                           acclimation_trait_table["Biochemical Assay", "K"]), "~","(", 
                                                         c(acclimation_trait_table["Physiological", "group_no"],
                                                           acclimation_trait_table["Biochemical Assay", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(Acclimation_Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(acclimation_trait_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait #(400x240)

# Preparing Graph

acclimation_trait_rnames_1 <- c("Biochemical Assay")

acclimation_trait_k_1 <- data.frame("k" = c(Acclimation_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                    row.names = acclimation_trait_rnames_1)

acclimation_trait_group_no_1 <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                           row.names = acclimation_trait_rnames_1)

acclimation_trait_study_1 <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                        row.names = acclimation_trait_rnames_1)

Acclimation_Trait_Model_Estimates_Reorder_1 <- Acclimation_Trait_Model_Estimates[c("Biochemical Assay"), ]

acclimation_trait_table_1 <- data.frame(estimate = Acclimation_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                      lowerCL = Acclimation_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                      upperCL = Acclimation_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                      K = acclimation_trait_k_1[,1], 
                                      group_no = acclimation_trait_group_no_1[,1], 
                                      row.names = acclimation_trait_rnames_1)
acclimation_trait_table_1$name <- row.names(acclimation_trait_table_1)

acclimation_trait_raw_mean_1 <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InRR_Transformed"))))

acclimation_trait_raw_name_1 <- c(replicate(27, "Biochemical Assay"))

acclimation_trait_raw_df_1 <- data.frame("Model" = acclimation_trait_raw_name_1, 
                                         "Effect" = acclimation_trait_raw_mean_1)

# Graph code

Acclimation_Trait_Order_1 <- c("Biochemical Assay")

density_acclimation_trait_1 <- acclimation_trait_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = acclimation_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)), xmin = min(acclimation_trait_raw_df_1$Effect)-0.03, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)), xmin = max(acclimation_trait_raw_df_1$Effect)+0.03, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               scale_colour_manual(values = c("#2B4E7A")) +
                               scale_fill_manual(values = c("#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_trait_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(acclimation_trait_table_1["Biochemical Assay", "K"]), "~","(", 
                                                           c(acclimation_trait_table_1["Biochemical Assay", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(acclimation_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait_1 #(400x160)

# Preparing Graph

acclimation_trait_rnames_2 <- c("Physiological")

acclimation_trait_k_2 <- data.frame("k" = c(Acclimation_Trait_Exploration["Physiological", "Freq"]), 
                                    row.names = acclimation_trait_rnames_2)

acclimation_trait_group_no_2 <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Physiological", "Freq"]), 
                                           row.names = acclimation_trait_rnames_2)

acclimation_trait_study_2 <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Physiological", "Freq"]), 
                                        row.names = acclimation_trait_rnames_2)

Acclimation_Trait_Model_Estimates_Reorder_2 <- Acclimation_Trait_Model_Estimates[c("Physiological"), ]

acclimation_trait_table_2 <- data.frame(estimate = Acclimation_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                      lowerCL = Acclimation_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                      upperCL = Acclimation_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                      K = acclimation_trait_k_2[,1], 
                                      group_no = acclimation_trait_group_no_2[,1], 
                                      row.names = acclimation_trait_rnames_2)
acclimation_trait_table_2$name <- row.names(acclimation_trait_table_2)

acclimation_trait_raw_mean_2 <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InRR_Transformed"))))

acclimation_trait_raw_name_2 <- c(replicate(32, "Physiological"))

acclimation_trait_raw_df_2 <- data.frame("Model" = acclimation_trait_raw_name_2, 
                                         "Effect" = acclimation_trait_raw_mean_2)

# Graph code

Acclimation_Trait_Order_2 <- c("Physiological")

density_acclimation_trait_2 <- acclimation_trait_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = acclimation_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)), xmin = min(acclimation_trait_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)), xmin = max(acclimation_trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_trait_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(acclimation_trait_table_2["Physiological", "K"]), "~","(", 
                                                           c(acclimation_trait_table_2["Physiological", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(acclimation_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait_2 #(400x160)

##### Acclimation Subset Model - Life-History Stage Meta-Regression #####
Acclimation_Stage_Exploration <- Acclimation_Subset_Data %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Exploration) <- Acclimation_Stage_Exploration$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Data <- Acclimation_Subset_Data %>% filter(Acclimation_Life.History_Stage_Category != "Embryo")

Acclimation_Stage_Species_Count <- Acclimation_Stage_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Species_Count) <- Acclimation_Stage_Species_Count$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Study_Count <- Acclimation_Stage_Data %>% select("Study_ID", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Study_Count) <- Acclimation_Stage_Study_Count$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Species <- Acclimation_Stage_Data %>% select("phylo") %>% unique()

Acclimation_Stage_A_cor <- as.data.frame(A_cor)
Acclimation_Stage_A_cor <- Acclimation_Stage_A_cor[c(Acclimation_Stage_Species$phylo), c(Acclimation_Stage_Species$phylo)]
Acclimation_Stage_A_cor <- as.matrix(Acclimation_Stage_A_cor)

Acclimation_Stage_VCV <- make_VCV_matrix(Acclimation_Stage_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Stage_Model <- metafor::rma.mv(InRR_Transformed, V = Acclimation_Stage_VCV, test = "t", dfs = "contain",
                                               mods = ~ Acclimation_Life.History_Stage_Category - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Acclimation_Stage_A_cor), data = Acclimation_Stage_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Stage_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Stage_Model.rds")
  } else {
    Acclimation_Stage_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Stage_Model.rds")})

Acclimation_Stage_Model_rob <- robust(Acclimation_Stage_Model, cluster = Acclimation_Stage_Data$Study_ID, adjust = TRUE)

Acclimation_Stage_Model_Estimates <- data.frame(Stage = substr(row.names(Acclimation_Stage_Model$b), 40, 100),
                                                estimate = Acclimation_Stage_Model$b, ci.lb = Acclimation_Stage_Model$ci.lb, 
                                                ci.ub = Acclimation_Stage_Model$ci.ub)
rownames(Acclimation_Stage_Model_Estimates) <- Acclimation_Stage_Model_Estimates$Stage
Acclimation_Stage_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Stage_Model), 2))

# Preparing Graph - Combined

acclimation_stage_rnames <- c("Adult", "Juvenile", "Larva")

acclimation_stage_k <- data.frame("k" = c(Acclimation_Stage_Exploration["Adult", "Freq"], 
                                          Acclimation_Stage_Exploration["Juvenile", "Freq"],
                                          Acclimation_Stage_Exploration["Larvae", "Freq"]), 
                                  row.names = acclimation_stage_rnames)

acclimation_stage_group_no <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Adult", "Freq"], 
                                                       Acclimation_Stage_Species_Count["Juvenile", "Freq"],
                                                       Acclimation_Stage_Species_Count["Larvae", "Freq"]), 
                                         row.names = acclimation_stage_rnames)

acclimation_stage_study <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Adult", "Freq"], 
                                                  Acclimation_Stage_Study_Count["Juvenile", "Freq"],
                                                  Acclimation_Stage_Study_Count["Larvae", "Freq"]), 
                                      row.names = acclimation_stage_rnames)

acclimation_stage_table <- data.frame(estimate = Acclimation_Stage_Model_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Stage_Model_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Stage_Model_Estimates[,"ci.ub"], 
                                      K = acclimation_stage_k[,1], 
                                      group_no = acclimation_stage_group_no[,1], 
                                      row.names = acclimation_stage_rnames)
acclimation_stage_table$name <- row.names(acclimation_stage_table)

acclimation_stage_raw_mean <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Adult") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Juvenile") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Larvae") %>% 
                                                select("InRR_Transformed"))))

acclimation_stage_raw_name <- c(replicate(16, "Adult"), 
                                replicate(15, "Juvenile"), 
                                replicate(30, "Larva"))

acclimation_stage_raw_df <- data.frame("Model" = acclimation_stage_raw_name, 
                                       "Effect" = acclimation_stage_raw_mean)

# Graph code - Combined

Acclimation_Stage_Order <- c("Larva", "Juvenile", "Adult")

density_acclimation_stage <- acclimation_stage_table %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = acclimation_stage_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table)[1], 1)), xmin = min(acclimation_stage_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table)[1], 1)), xmin = max(acclimation_stage_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                             coord_cartesian(xlim = c(-0.5, 0.5)) +
                             annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_stage_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(acclimation_stage_table["Larva", "K"],
                                                           acclimation_stage_table["Juvenile", "K"],
                                                           acclimation_stage_table["Adult", "K"]), "~","(", 
                                                         c(acclimation_stage_table["Larva", "group_no"],
                                                           acclimation_stage_table["Juvenile", "group_no"],
                                                           acclimation_stage_table["Adult", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                            geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                   paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                   paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Adult", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(acclimation_stage_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage #(400x320)

# Preparing Graph - Part 1

acclimation_stage_rnames_1 <- c("Adult", "Juvenile")

acclimation_stage_k_1 <- data.frame("k" = c(Acclimation_Stage_Exploration["Adult", "Freq"], 
                                            Acclimation_Stage_Exploration["Juvenile", "Freq"]), 
                                    row.names = acclimation_stage_rnames_1)

acclimation_stage_group_no_1 <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Adult", "Freq"], 
                                                         Acclimation_Stage_Species_Count["Juvenile", "Freq"]), 
                                           row.names = acclimation_stage_rnames_1)

acclimation_stage_study_1 <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Adult", "Freq"], 
                                                    Acclimation_Stage_Study_Count["Juvenile", "Freq"]), 
                                        row.names = acclimation_stage_rnames_1)

Acclimation_Stage_Model_Estimates_Reorder_1 <- Acclimation_Stage_Model_Estimates[c("Adult", "Juvenile"), ]

acclimation_stage_table_1 <- data.frame(estimate = Acclimation_Stage_Model_Estimates_Reorder_1[,"estimate"], 
                                        lowerCL = Acclimation_Stage_Model_Estimates_Reorder_1[,"ci.lb"], 
                                        upperCL = Acclimation_Stage_Model_Estimates_Reorder_1[,"ci.ub"], 
                                        K = acclimation_stage_k_1[,1], 
                                        group_no = acclimation_stage_group_no_1[,1], 
                                        row.names = acclimation_stage_rnames_1)
acclimation_stage_table_1$name <- row.names(acclimation_stage_table_1)

acclimation_stage_raw_mean_1 <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Adult") %>% 
                                                select("InRR_Transformed"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Juvenile") %>% 
                                                select("InRR_Transformed"))))

acclimation_stage_raw_name_1 <- c(replicate(16, "Adult"), 
                                  replicate(15, "Juvenile"))

acclimation_stage_raw_df_1 <- data.frame("Model" = acclimation_stage_raw_name_1, 
                                         "Effect" = acclimation_stage_raw_mean_1)

# Graph code - Part 1

Acclimation_Stage_Order_1 <- c("Juvenile", "Adult")

density_acclimation_stage_1 <- acclimation_stage_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = acclimation_stage_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)), xmin = min(acclimation_stage_raw_df_1$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)), xmin = max(acclimation_stage_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                               scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_stage_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(acclimation_stage_table_1["Juvenile", "K"],
                                                             acclimation_stage_table_1["Adult", "K"]), "~","(", 
                                                           c(acclimation_stage_table_1["Juvenile", "group_no"],
                                                             acclimation_stage_table_1["Adult", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Adult", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(acclimation_stage_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage_1 #(400x240)

# Preparing Graph - Part 2

acclimation_stage_rnames_2 <- c("Larva")

acclimation_stage_k_2 <- data.frame("k" = c(Acclimation_Stage_Exploration["Larvae", "Freq"]), 
                                    row.names = acclimation_stage_rnames_2)

acclimation_stage_group_no_2 <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Larvae", "Freq"]), 
                                           row.names = acclimation_stage_rnames_2)

acclimation_stage_study_2 <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Larvae", "Freq"]), 
                                        row.names = acclimation_stage_rnames_2)

Acclimation_Stage_Model_Estimates_Reorder_2 <- Acclimation_Stage_Model_Estimates[c("Larva"), ]

acclimation_stage_table_2 <- data.frame(estimate = Acclimation_Stage_Model_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Acclimation_Stage_Model_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Acclimation_Stage_Model_Estimates_Reorder_2[,"ci.ub"], 
                                        K = acclimation_stage_k_2[,1], 
                                        group_no = acclimation_stage_group_no_2[,1], 
                                        row.names = acclimation_stage_rnames_2)
acclimation_stage_table_2$name <- row.names(acclimation_stage_table_2)

acclimation_stage_raw_mean_2 <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Larvae") %>% 
                                                select("InRR_Transformed"))))

acclimation_stage_raw_name_2 <- c(replicate(30, "Larva"))

acclimation_stage_raw_df_2 <- data.frame("Model" = acclimation_stage_raw_name_2, 
                                         "Effect" = acclimation_stage_raw_mean_2)

# Graph code - Part 2

Acclimation_Stage_Order_2 <- c("Larva")

density_acclimation_stage_2 <- acclimation_stage_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = acclimation_stage_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)), xmin = min(acclimation_stage_raw_df_2$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)), xmin = max(acclimation_stage_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_stage_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(acclimation_stage_table_2["Larva", "K"]), "~","(", 
                                                           c(acclimation_stage_table_2["Larva", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.4, y = (seq(1, dim(acclimation_stage_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage_2 #(400x160)

##### Acclimation Subset Model - Class Meta-Regression #####
Acclimation_Class_Exploration <- Acclimation_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Exploration) <- Acclimation_Class_Exploration$Class

Acclimation_Class_Data <- Acclimation_Subset_Data %>% filter(Class != "Clitellata" &
                                                             Class != "Gastropoda" &
                                                             Class != "Holothuroidea" &
                                                             Class != "Malacostraca")

Acclimation_Class_Species_Count <- Acclimation_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Species_Count) <- Acclimation_Class_Species_Count$Class

Acclimation_Class_Study_Count <- Acclimation_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Study_Count) <- Acclimation_Class_Study_Count$Class

Acclimation_Class_Species <- Acclimation_Class_Data %>% select("phylo") %>% unique()

Acclimation_Class_A_cor <- as.data.frame(A_cor)
Acclimation_Class_A_cor <- Acclimation_Class_A_cor[c(Acclimation_Class_Species$phylo), c(Acclimation_Class_Species$phylo)]
Acclimation_Class_A_cor <- as.matrix(Acclimation_Class_A_cor)

Acclimation_Class_VCV <- make_VCV_matrix(Acclimation_Class_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Class_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Acclimation_Class_VCV, test = "t", dfs = "contain",
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Acclimation_Class_A_cor), data = Acclimation_Class_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Class_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Class_Model.rds")
  } else {
    Acclimation_Class_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Class_Model.rds")})

Acclimation_Class_Model_rob <- robust(Acclimation_Class_Model, cluster = Acclimation_Class_Data$Study_ID, adjust = TRUE)

Acclimation_Class_Model_Estimates <- data.frame(Class = substr(row.names(Acclimation_Class_Model$b), 6, 100),
                                                estimate = Acclimation_Class_Model$b, 
                                                ci.lb = Acclimation_Class_Model$ci.lb, 
                                                ci.ub = Acclimation_Class_Model$ci.ub)
rownames(Acclimation_Class_Model_Estimates) <- c("Insecta")
Acclimation_Class_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Class_Model), 2))

# Preparing Graph

acclimation_class_rnames <- c("Insecta")

acclimation_class_k <- data.frame("k" = c(Acclimation_Class_Exploration["Insecta", "Freq"]), 
                                  row.names = acclimation_class_rnames)

acclimation_class_group_no <- data.frame("Spp No." = c(Acclimation_Class_Species_Count["Insecta", "Freq"]), 
                                         row.names = acclimation_class_rnames)

acclimation_class_study <- data.frame("Study" = c(Acclimation_Class_Study_Count["Insecta", "Freq"]), 
                                      row.names = acclimation_class_rnames)

acclimation_class_table <- data.frame(estimate = Acclimation_Class_Model_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Class_Model_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Class_Model_Estimates[,"ci.ub"], 
                                      K = acclimation_class_k[,1], 
                                      group_no = acclimation_class_group_no[,1], 
                                      row.names = acclimation_class_rnames)
acclimation_class_table$name <- row.names(acclimation_class_table)

acclimation_class_raw_mean <- c(unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                select("InRR_Transformed"))))

acclimation_class_raw_name <- c(replicate(31, "Insecta"))

acclimation_class_raw_df <- data.frame("Model" = acclimation_class_raw_name, 
                                       "Effect" = acclimation_class_raw_mean)

# Graph code

Acclimation_Class_Order <- c("Insecta")

density_acclimation_class <- acclimation_class_table %>% mutate(name = fct_relevel(name, Acclimation_Class_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = acclimation_class_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Class_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table)[1], 1)), xmin = min(acclimation_class_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                            size = 1) +
                             geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table)[1], 1)), xmin = max(acclimation_class_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                             coord_cartesian(xlim = c(-0.5, 0.5)) +
                             annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_class_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(acclimation_class_table["Insecta", "K"]), "~","(", 
                                                         c(acclimation_class_table["Insecta", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.4, y = (seq(1, dim(acclimation_class_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_class #(400x160)

##### Acclimation Subset Model - Specific Trait Meta-Regression #####
Acclimation_Specific_Trait_Exploration <- Acclimation_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Acclimation_Specific_Trait_Exploration <- Acclimation_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Acclimation_Specific_Trait_Exploration) <- Acclimation_Specific_Trait_Exploration$Measurement

Acclimation_Specific_Trait_Data <- Acclimation_Subset_Data %>% filter(Measurement == "Metabolic Rate")

Acclimation_Specific_Trait_Species_Count <- Acclimation_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Acclimation_Specific_Trait_Species_Count) <- Acclimation_Specific_Trait_Species_Count$Measurement

Acclimation_Specific_Trait_Study_Count <- Acclimation_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Acclimation_Specific_Trait_Study_Count) <- Acclimation_Specific_Trait_Study_Count$Measurement

Acclimation_Specific_Trait_Species <- Acclimation_Specific_Trait_Data %>% select("phylo") %>% unique()

Acclimation_Specific_Trait_A_cor <- as.data.frame(A_cor)
Acclimation_Specific_Trait_A_cor <- Acclimation_Specific_Trait_A_cor[c(Acclimation_Specific_Trait_Species$phylo), c(Acclimation_Specific_Trait_Species$phylo)]
Acclimation_Specific_Trait_A_cor <- as.matrix(Acclimation_Specific_Trait_A_cor)

Acclimation_Specific_Trait_VCV <- make_VCV_matrix(Acclimation_Specific_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Specific_Trait_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Acclimation_Specific_Trait_VCV, test = "t", dfs = "contain",
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number), 
                                                        R = list(phylo=Acclimation_Specific_Trait_A_cor), data = Acclimation_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Specific_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Specific_Trait_Model.rds")
  } else {
    Acclimation_Specific_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Acclimation_Specific_Trait_Model.rds")})

Acclimation_Specific_Trait_Model_rob <- robust(Acclimation_Specific_Trait_Model, cluster = Acclimation_Specific_Trait_Data$Study_ID, adjust = TRUE)

Acclimation_Specific_Trait_Model_Estimates <- data.frame(Trait = substr(row.names(Acclimation_Specific_Trait_Model$b), 12, 100),
                                                         estimate = Acclimation_Specific_Trait_Model$b, ci.lb = Acclimation_Specific_Trait_Model$ci.lb, 
                                                         ci.ub = Acclimation_Specific_Trait_Model$ci.ub)
rownames(Acclimation_Specific_Trait_Model_Estimates) <- c("Metabolic Rate")
Acclimation_Specific_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Specific_Trait_Model), 2))

# Preparing Graph

acclimation_specific_trait_rnames <- c("Metabolic Rate")

acclimation_specific_trait_k <- data.frame("k" = c(Acclimation_Specific_Trait_Exploration["Metabolic Rate", "Freq"]), 
                                           row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_group_no <- data.frame("Spp No." = c(Acclimation_Specific_Trait_Species_Count["Metabolic Rate", "Freq"]), 
                                                  row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_study <- data.frame("Study" = c(Acclimation_Specific_Trait_Study_Count["Metabolic Rate", "Freq"]), 
                                               row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_table <- data.frame(estimate = Acclimation_Specific_Trait_Model_Estimates[,"estimate"], 
                                               lowerCL = Acclimation_Specific_Trait_Model_Estimates[,"ci.lb"], 
                                               upperCL = Acclimation_Specific_Trait_Model_Estimates[,"ci.ub"], 
                                               K = acclimation_specific_trait_k[,1], 
                                               group_no = acclimation_specific_trait_group_no[,1], 
                                               row.names = acclimation_specific_trait_rnames)
acclimation_specific_trait_table$name <- row.names(acclimation_specific_trait_table)

acclimation_specific_trait_raw_mean <- c(unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                         select("InRR_Transformed"))))

acclimation_specific_trait_raw_name <- c(replicate(12, "Metabolic Rate"))

acclimation_specific_trait_raw_df <- data.frame("Model" = acclimation_specific_trait_raw_name, 
                                                "Effect" = acclimation_specific_trait_raw_mean)

# Graph code

Acclimation_Specific_Trait_Order <- c("Metabolic Rate")

density_acclimation_specific_trait <- acclimation_specific_trait_table %>% mutate(name = fct_relevel(name, Acclimation_Specific_Trait_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = acclimation_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Specific_Trait_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                      geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)), xmin = min(acclimation_specific_trait_raw_df$Effect)-0.03, xmax = -1.5, colour = name),
                                                     size = 1) +
                                      geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)), xmin = max(acclimation_specific_trait_raw_df$Effect)+0.03, xmax = 1.5, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                      scale_colour_manual(values = c("#2B4E7A")) +
                                      scale_fill_manual(values = c("#2B4E7A")) +
                                      coord_cartesian(xlim = c(-0.5, 0.5)) +
                                      annotate('text',  x = 0.5, y = (seq(1, dim(acclimation_specific_trait_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(acclimation_specific_trait_table["Metabolic Rate", "K"]), "~","(", 
                                                                  c(acclimation_specific_trait_table["Metabolic Rate", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.4, y = (seq(1, dim(acclimation_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_specific_trait #(400x160)

##### Developmental Subset Model #####
Developmental_Subset_Data <- Individual_Subset_Data %>% filter(Plasticity_Mechanism == "Developmental Plasticity")
Developmental_Species <- Developmental_Subset_Data %>% select("phylo") %>% unique()

Developmental_A_cor <- as.data.frame(A_cor)
Developmental_A_cor <- Developmental_A_cor[c(Developmental_Species$phylo), c(Developmental_Species$phylo)]
Developmental_A_cor <- as.matrix(Developmental_A_cor)

Developmental_VCV <- make_VCV_matrix(Developmental_Subset_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Model <- metafor::rma.mv(InRR_Transformed ~ 1, V = Developmental_VCV, test = "t", dfs = "contain",
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                                           R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Developmental_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Model.rds")
  } else {
    Developmental_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Model.rds")})

Developmental_Model_rob <- robust(Developmental_Model, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Model_Estimates <- data.frame(estimate = Developmental_Model$b, ci.lb = Developmental_Model$ci.lb, ci.ub = Developmental_Model$ci.ub)
Developmental_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Model), 2))

#### Developmental Subset Model - Fluctuation Amplitude Meta-Regression ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Amplitude_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_VCV, test = "t", dfs = "contain",
                                                     mods = ~ Fluctuation_Magnitude - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Developmental_Amplitude_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Amplitude_Model.rds")
  } else {
    Developmental_Amplitude_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Amplitude_Model.rds")})

Developmental_Amplitude_Model_rob <- robust(Developmental_Amplitude_Model, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Amplitude_Model_Estimates <- data.frame(estimate = Developmental_Amplitude_Model$b, ci.lb = Developmental_Amplitude_Model$ci.lb, 
                                                      ci.ub = Developmental_Amplitude_Model$ci.ub)
Developmental_Amplitude_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Amplitude_Model), 2))

# Graph Preparing

Developmental_Plot_Data <- Developmental_Subset_Data
Developmental_Plot_Data <- Developmental_Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                                           ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                           ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))

# Graph Code

Developmental_Amplitude_Plot <- ggplot(Developmental_Plot_Data, aes(x = Fluctuation_Magnitude, y = InRR_Transformed)) + 
                                geom_point(aes(x = Fluctuation_Magnitude, y = InRR_Transformed, 
                                               size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                                               shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
                                labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (PRRD)", 
                                     size = "Sample Size", title = "Developmental Treatments") +
                                theme_bw() +
                                theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                #theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                geom_hline(yintercept = Developmental_Model_Estimates$estimate, lty = 2) + 
                                geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                stat_poly_eq(formula = y ~ x, 
                                             aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                                 parse = TRUE) +
                                coord_cartesian(xlim = c(0, 25), 
                                                ylim = c(-0.25, 0.25))

Developmental_Amplitude_Plot #(330x330)

#### Developmental Subset Model - Type of Fluctuation Meta-Regression ####
Developmental_Fluctuation_Data <- Developmental_Subset_Data %>% filter(!is.na(Fluctuation_Category))
Developmental_Fluctuation_Exploration <- Developmental_Subset_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Exploration) <- Developmental_Fluctuation_Exploration$Fluctuation_Category

Developmental_Fluctuation_Data <- Developmental_Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")

Developmental_Fluctuation_Species_Count <- Developmental_Subset_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Species_Count) <- Developmental_Fluctuation_Species_Count$Fluctuation_Category

Developmental_Fluctuation_Study_Count <- Developmental_Subset_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Study_Count) <- Developmental_Fluctuation_Study_Count$Fluctuation_Category

Developmental_Fluctuation_Species <- Developmental_Fluctuation_Data %>% select("phylo") %>% unique()

Developmental_Fluctuation_A_cor <- as.data.frame(A_cor)
Developmental_Fluctuation_A_cor <- Developmental_Fluctuation_A_cor[c(Developmental_Fluctuation_Species$phylo), c(Developmental_Fluctuation_Species$phylo)]
Developmental_Fluctuation_A_cor <- as.matrix(Developmental_Fluctuation_A_cor)

Developmental_Fluctuation_VCV <- make_VCV_matrix(Developmental_Fluctuation_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Fluctuation_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_Fluctuation_VCV, test = "t", dfs = "contain",
                                                       mods = ~ Fluctuation_Category - 1,
                                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                                       R = list(phylo=Developmental_Fluctuation_A_cor), data = Developmental_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                       control=list(rel.tol=1e-9))
    saveRDS(Developmental_Fluctuation_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Fluctuation_Model.rds")
  } else {
    Developmental_Fluctuation_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Fluctuation_Model.rds")})

Developmental_Fluctuation_Model_rob <- robust(Developmental_Fluctuation_Model, cluster = Developmental_Fluctuation_Data$Study_ID, adjust = TRUE)

Developmental_Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Developmental_Fluctuation_Model$b), 21, 100),
                                                        estimate = Developmental_Fluctuation_Model$b, ci.lb = Developmental_Fluctuation_Model$ci.lb, 
                                                        ci.ub = Developmental_Fluctuation_Model$ci.ub)
rownames(Developmental_Fluctuation_Model_Estimates) <- Developmental_Fluctuation_Model_Estimates$Category
Developmental_Fluctuation_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Fluctuation_Model), 2))

# Preparing Graph - Combined

developmental_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

developmental_fluctuation_k <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                  Developmental_Fluctuation_Exploration["Alternating", "Freq"], 
                                                  Developmental_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                          row.names = developmental_fluctuation_rnames)

developmental_fluctuation_group_no <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                               Developmental_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                               Developmental_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                 row.names = developmental_fluctuation_rnames)

developmental_fluctuation_study <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                          Developmental_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                          Developmental_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                              row.names = developmental_fluctuation_rnames)

Developmental_Fluctuation_Model_Estimates_Reorder <- Developmental_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

developmental_fluctuation_table <- data.frame(estimate = Developmental_Fluctuation_Model_Estimates_Reorder[,"estimate"], 
                                              lowerCL = Developmental_Fluctuation_Model_Estimates_Reorder[,"ci.lb"], 
                                              upperCL = Developmental_Fluctuation_Model_Estimates_Reorder[,"ci.ub"], 
                                              K = developmental_fluctuation_k[,1], 
                                              group_no = developmental_fluctuation_group_no[,1], 
                                              row.names = developmental_fluctuation_rnames)
developmental_fluctuation_table$name <- row.names(developmental_fluctuation_table)

developmental_fluctuation_raw_mean <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                        select("InRR_Transformed"))), 
                                        unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                        select("InRR_Transformed"))), 
                                        unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                        select("InRR_Transformed"))))

developmental_fluctuation_raw_name <- c(replicate(49, "Sinusoidal (Sine Curve)"), 
                                        replicate(50, "Alternating"), 
                                        replicate(17, "Stepwise"))

developmental_fluctuation_raw_df <- data.frame("Model" = developmental_fluctuation_raw_name, 
                                               "Effect" = developmental_fluctuation_raw_mean)

# Graph code - Combined

Developmental_Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_developmental_fluctuation <- developmental_fluctuation_table %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = developmental_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)), xmin = min(developmental_fluctuation_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)), xmin = max(developmental_fluctuation_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                      vjust = c(-2.7, -2.7, -0.8))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                     scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                                     annotate('text',  x = 0.5, y = (seq(1, dim(developmental_fluctuation_table)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(developmental_fluctuation_table["Stepwise", "K"], 
                                                                   developmental_fluctuation_table["Alternating", "K"], 
                                                                   developmental_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                 c(developmental_fluctuation_table["Stepwise", "group_no"], 
                                                                   developmental_fluctuation_table["Alternating", "group_no"], 
                                                                   developmental_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.4, y = (seq(1, dim(developmental_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation #(400x320)

# Preparing Graph - Part 1

developmental_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

developmental_fluctuation_k_1 <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                    Developmental_Fluctuation_Exploration["Alternating", "Freq"]), 
                                            row.names = developmental_fluctuation_rnames_1)

developmental_fluctuation_group_no_1 <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                                 Developmental_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                   row.names = developmental_fluctuation_rnames_1)

developmental_fluctuation_study_1 <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                            Developmental_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                                row.names = developmental_fluctuation_rnames_1)

Developmental_Fluctuation_Model_Estimates_Reorder_1 <- Developmental_Fluctuation_Model_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

developmental_fluctuation_table_1 <- data.frame(estimate = Developmental_Fluctuation_Model_Estimates_Reorder_1[,"estimate"], 
                                                lowerCL = Developmental_Fluctuation_Model_Estimates_Reorder_1[,"ci.lb"], 
                                                upperCL = Developmental_Fluctuation_Model_Estimates_Reorder_1[,"ci.ub"], 
                                                K = developmental_fluctuation_k_1[,1], 
                                                group_no = developmental_fluctuation_group_no_1[,1], 
                                                row.names = developmental_fluctuation_rnames_1)
developmental_fluctuation_table_1$name <- row.names(developmental_fluctuation_table_1)

developmental_fluctuation_raw_mean_1 <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                        select("InRR_Transformed"))), 
                                          unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                        select("InRR_Transformed"))))

developmental_fluctuation_raw_name_1 <- c(replicate(49, "Sinusoidal (Sine Curve)"), 
                                          replicate(50, "Alternating"))

developmental_fluctuation_raw_df_1 <- data.frame("Model" = developmental_fluctuation_raw_name_1, 
                                                 "Effect" = developmental_fluctuation_raw_mean_1)

# Graph code - Part 1

Developmental_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_developmental_fluctuation_1 <- developmental_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order_1)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = developmental_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order_1)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)), xmin = min(developmental_fluctuation_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)), xmin = max(developmental_fluctuation_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-2.7, -0.8))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                       scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                       coord_cartesian(xlim = c(-0.5, 0.5)) +
                                       annotate('text',  x = 0.5, y = (seq(1, dim(developmental_fluctuation_table_1)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(developmental_fluctuation_table_1["Alternating", "K"], 
                                                                     developmental_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                   c(developmental_fluctuation_table_1["Alternating", "group_no"], 
                                                                     developmental_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.4, y = (seq(1, dim(developmental_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation_1 #(400x240)

# Preparing Graph - Part 2

developmental_fluctuation_rnames_2 <- c("Stepwise")

developmental_fluctuation_k_2 <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                            row.names = developmental_fluctuation_rnames_2)

developmental_fluctuation_group_no_2 <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                   row.names = developmental_fluctuation_rnames_2)

developmental_fluctuation_study_2 <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                                row.names = developmental_fluctuation_rnames_2)

Developmental_Fluctuation_Model_Estimates_Reorder_2 <- Developmental_Fluctuation_Model_Estimates[c("Stepwise"), ]

developmental_fluctuation_table_2 <- data.frame(estimate = Developmental_Fluctuation_Model_Estimates_Reorder_2[,"estimate"], 
                                                lowerCL = Developmental_Fluctuation_Model_Estimates_Reorder_2[,"ci.lb"], 
                                                upperCL = Developmental_Fluctuation_Model_Estimates_Reorder_2[,"ci.ub"], 
                                                K = developmental_fluctuation_k_2[,1], 
                                                group_no = developmental_fluctuation_group_no_2[,1], 
                                                row.names = developmental_fluctuation_rnames_2)
developmental_fluctuation_table_2$name <- row.names(developmental_fluctuation_table_2)

developmental_fluctuation_raw_mean_2 <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                        select("InRR_Transformed"))))

developmental_fluctuation_raw_name_2 <- c(replicate(17, "Stepwise"))

developmental_fluctuation_raw_df_2 <- data.frame("Model" = developmental_fluctuation_raw_name_2, 
                                                 "Effect" = developmental_fluctuation_raw_mean_2)

# Graph code - Part 2

Developmental_Fluctuation_Order_2 <- c("Stepwise")

density_developmental_fluctuation_2 <- developmental_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order_2)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = developmental_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order_2)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                               scale = 0.025, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)), xmin = min(developmental_fluctuation_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)), xmin = max(developmental_fluctuation_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                       coord_cartesian(xlim = c(-0.5, 0.5)) +
                                       annotate('text',  x = 0.5, y = (seq(1, dim(developmental_fluctuation_table_2)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(developmental_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                   c(developmental_fluctuation_table_2["Stepwise", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.4, y = (seq(1, dim(developmental_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation_2 #(400x160)

##### Developmental Subset Model - Trait Meta-Regression #####
Developmental_Trait_Exploration <- Developmental_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Exploration) <- Developmental_Trait_Exploration$Trait_Category

Developmental_Trait_Data <- Developmental_Subset_Data %>% filter(Trait_Category != "Behavioural" &
                                                                 Trait_Category != "Biochemical Assay" &
                                                                 Trait_Category != "Gene Expression" &
                                                                 Trait_Category != "Physiological")

Developmental_Trait_Species_Count <- Developmental_Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Species_Count) <- Developmental_Trait_Species_Count$Trait_Category

Developmental_Trait_Study_Count <- Developmental_Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Study_Count) <- Developmental_Trait_Study_Count$Trait_Category

Developmental_Trait_Species <- Developmental_Trait_Data %>% select("phylo") %>% unique()

Developmental_Trait_A_cor <- as.data.frame(A_cor)
Developmental_Trait_A_cor <- Developmental_Trait_A_cor[c(Developmental_Trait_Species$phylo), c(Developmental_Trait_Species$phylo)]
Developmental_Trait_A_cor <- as.matrix(Developmental_Trait_A_cor)

Developmental_Trait_VCV <- make_VCV_matrix(Developmental_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_Trait_VCV, test = "t", dfs = "contain",
                                                 mods = ~ Trait_Category - 1,
                                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                               ~1|Shared_Animal_Number, ~1|Measurement), 
                                                 R = list(phylo=Developmental_Trait_A_cor), data = Developmental_Trait_Data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
    saveRDS(Developmental_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Trait_Model.rds")
  } else {
    Developmental_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Trait_Model.rds")})

Developmental_Trait_Model_rob <- robust(Developmental_Trait_Model, cluster = Developmental_Trait_Data$Study_ID, adjust = TRUE)

Developmental_Trait_Model_Estimates <- data.frame(Category = substr(row.names(Developmental_Trait_Model$b), 15, 100),
                                                  estimate = Developmental_Trait_Model$b, 
                                                  ci.lb = Developmental_Trait_Model$ci.lb, 
                                                  ci.ub = Developmental_Trait_Model$ci.ub)
rownames(Developmental_Trait_Model_Estimates) <- Developmental_Trait_Model_Estimates$Category
Developmental_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Trait_Model), 2))

# Preparing Graph - Combined

developmental_trait_rnames <- c("Life-history Traits", "Morphology")

developmental_trait_k <- data.frame("k" = c(Developmental_Trait_Exploration["Life-History Traits", "Freq"], 
                                            Developmental_Trait_Exploration["Morphology", "Freq"]), 
                                    row.names = developmental_trait_rnames)

developmental_trait_group_no <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Life-History Traits", "Freq"],
                                                         Developmental_Trait_Species_Count["Morphology", "Freq"]), 
                                           row.names = developmental_trait_rnames)

developmental_trait_study <- data.frame("Study" = c(Developmental_Trait_Study_Count["Life-History Traits", "Freq"],
                                                    Developmental_Trait_Study_Count["Morphology", "Freq"]), 
                                        row.names = developmental_trait_rnames)

developmental_trait_table <- data.frame(estimate = Developmental_Trait_Model_Estimates[,"estimate"], 
                                        lowerCL = Developmental_Trait_Model_Estimates[,"ci.lb"], 
                                        upperCL = Developmental_Trait_Model_Estimates[,"ci.ub"], 
                                        K = developmental_trait_k[,1], 
                                        group_no = developmental_trait_group_no[,1], 
                                        row.names = developmental_trait_rnames)
developmental_trait_table$name <- row.names(developmental_trait_table)

developmental_trait_raw_mean <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                  select("InRR_Transformed"))), 
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                  select("InRR_Transformed"))))

developmental_trait_raw_name <- c(replicate(65, "Life-history Traits"), 
                                  replicate(54, "Morphology"))

developmental_trait_raw_df <- data.frame("Model" = developmental_trait_raw_name, 
                                         "Effect" = developmental_trait_raw_mean)

# Graph code - Combined

Developmental_Trait_Order <- c("Morphology", "Life-history Traits")

density_developmental_trait <- developmental_trait_table %>% mutate(name = fct_relevel(name, Developmental_Trait_Order)) %>%
                               ggplot() +
                               geom_density_ridges(data = developmental_trait_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table)[1], 1)), xmin = min(developmental_trait_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table)[1], 1)), xmin = max(developmental_trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (PRRD)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                vjust = c(-2.7, -0.8))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                               scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(developmental_trait_table)[1], 1)+0.4),
                               label= paste("italic(k)==", c(developmental_trait_table["Morphology", "K"], 
                                                             developmental_trait_table["Life-history Traits", "K"]), "~","(", 
                                                           c(developmental_trait_table["Morphology", "group_no"], 
                                                             developmental_trait_table["Life-history Traits", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                      paste(format(round(mean(exp(Developmental_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.4, y = (seq(1, dim(developmental_trait_table)[1], 1)+0.4)), size = 3.5)

density_developmental_trait #(400x240)

# Preparing Graph - Part 1

developmental_trait_rnames_1 <- c("Life-history Traits")

developmental_trait_k_1 <- data.frame("k" = c(Developmental_Trait_Exploration["Life-History Traits", "Freq"]), 
                                      row.names = developmental_trait_rnames_1)

developmental_trait_group_no_1 <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Life-History Traits", "Freq"]), 
                                             row.names = developmental_trait_rnames_1)

developmental_trait_study_1 <- data.frame("Study" = c(Developmental_Trait_Study_Count["Life-History Traits", "Freq"]), 
                                          row.names = developmental_trait_rnames_1)

Developmental_Trait_Model_Estimates_Reorder_1 <- Developmental_Trait_Model_Estimates[c("Life-History Traits"), ]

developmental_trait_table_1 <- data.frame(estimate = Developmental_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                          lowerCL = Developmental_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                          upperCL = Developmental_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                          K = developmental_trait_k_1[,1], 
                                          group_no = developmental_trait_group_no_1[,1], 
                                          row.names = developmental_trait_rnames_1)
developmental_trait_table_1$name <- row.names(developmental_trait_table_1)

developmental_trait_raw_mean_1 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                  select("InRR_Transformed"))))

developmental_trait_raw_name_1 <- c(replicate(65, "Life-history Traits"))

developmental_trait_raw_df_1 <- data.frame("Model" = developmental_trait_raw_name_1, 
                                           "Effect" = developmental_trait_raw_mean_1)

# Graph code - Part 1

Developmental_Trait_Order_1 <- c("Life-history Traits")

density_developmental_trait_1 <- developmental_trait_table_1 %>% mutate(name = fct_relevel(name, Developmental_Trait_Order_1)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = developmental_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order_1)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)), xmin = min(developmental_trait_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)), xmin = max(developmental_trait_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 scale_colour_manual(values = c("#2B4E7A")) +
                                 scale_fill_manual(values = c("#2B4E7A")) +
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(developmental_trait_table_1)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(developmental_trait_table_1["Life-history Traits", "K"]), "~","(", 
                                                             c(developmental_trait_table_1["Life-history Traits", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(developmental_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_trait_1 #(400x160)

# Preparing Graph - Part 2

developmental_trait_rnames_2 <- c("Morphology")

developmental_trait_k_2 <- data.frame("k" = c(Developmental_Trait_Exploration["Morphology", "Freq"]), 
                                      row.names = developmental_trait_rnames_2)

developmental_trait_group_no_2 <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Morphology", "Freq"]), 
                                             row.names = developmental_trait_rnames_2)

developmental_trait_study_2 <- data.frame("Study" = c(Developmental_Trait_Study_Count["Morphology", "Freq"]), 
                                          row.names = developmental_trait_rnames_2)

Developmental_Trait_Model_Estimates_Reorder_2 <- Developmental_Trait_Model_Estimates[c("Morphology"), ]

developmental_trait_table_2 <- data.frame(estimate = Developmental_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Developmental_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Developmental_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                        K = developmental_trait_k_2[,1], 
                                        group_no = developmental_trait_group_no_2[,1], 
                                        row.names = developmental_trait_rnames_2)
developmental_trait_table_2$name <- row.names(developmental_trait_table_2)

developmental_trait_raw_mean_2 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                  select("InRR_Transformed"))))

developmental_trait_raw_name_2 <- c(replicate(54, "Morphology"))

developmental_trait_raw_df_2 <- data.frame("Model" = developmental_trait_raw_name_2, 
                                           "Effect" = developmental_trait_raw_mean_2) 

# Graph code - Part 2

Developmental_Trait_Order_2 <- c("Morphology")

density_developmental_trait_2 <- developmental_trait_table_2 %>% mutate(name = fct_relevel(name, Developmental_Trait_Order_2)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = developmental_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order_2)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.02, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)), xmin = min(developmental_trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)), xmin = max(developmental_trait_raw_df_2$Effect)+0.01, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(developmental_trait_table_2)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(developmental_trait_table_2["Morphology", "K"]), "~","(", 
                                                             c(developmental_trait_table_2["Morphology", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.4, y = (seq(1, dim(developmental_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_trait_2 #(400x160)

##### Developmental Subset Model - Exposure Time Meta-Regression #####
Developmental_Exposure_Exploration <- Developmental_Subset_Data %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Exploration) <- Developmental_Exposure_Exploration$Developmental_Exposure_Time_Category

Developmental_Exposure_Species_Count <- Developmental_Subset_Data %>% select("Scientific_Name", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Species_Count) <- Developmental_Exposure_Species_Count$Developmental_Exposure_Time_Category

Developmental_Exposure_Study_Count <- Developmental_Subset_Data %>% select("Study_ID", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Study_Count) <- Developmental_Exposure_Study_Count$Developmental_Exposure_Time_Category

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Exposure_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_VCV, test = "t", dfs = "contain",
                                                    mods = ~ Developmental_Exposure_Time_Category - 1,
                                                    random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                  ~1|Shared_Animal_Number, ~1|Measurement), 
                                                    R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                                    control=list(rel.tol=1e-9))
    saveRDS(Developmental_Exposure_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Exposure_Model.rds")
  } else {
    Developmental_Exposure_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Exposure_Model.rds")})

Developmental_Exposure_Model_rob <- robust(Developmental_Exposure_Model, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Exposure_Model_Estimates <- data.frame(Category = substr(row.names(Developmental_Exposure_Model$b), 37, 100),
                                                     estimate = Developmental_Exposure_Model$b, ci.lb = Developmental_Exposure_Model$ci.lb, 
                                                     ci.ub = Developmental_Exposure_Model$ci.ub)
rownames(Developmental_Exposure_Model_Estimates) <- Developmental_Exposure_Model_Estimates$Category
Developmental_Exposure_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Exposure_Model), 2))

# Preparing Graph - Combined

developmental_exposure_rnames <- c("Embryo", "Juvenile", "Larva")

developmental_exposure_k <- data.frame("k" = c(Developmental_Exposure_Exploration["Embryo", "Freq"], 
                                               Developmental_Exposure_Exploration["Juvenile", "Freq"], 
                                               Developmental_Exposure_Exploration["Larvae", "Freq"]), 
                                       row.names = developmental_exposure_rnames)

developmental_exposure_group_no <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Embryo", "Freq"], 
                                                            Developmental_Exposure_Species_Count["Juvenile", "Freq"], 
                                                            Developmental_Exposure_Species_Count["Larvae", "Freq"]), 
                                              row.names = developmental_exposure_rnames)

developmental_exposure_study <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Embryo", "Freq"], 
                                                       Developmental_Exposure_Study_Count["Juvenile", "Freq"], 
                                                       Developmental_Exposure_Study_Count["Larvae", "Freq"]), 
                                           row.names = developmental_exposure_rnames)

developmental_exposure_table <- data.frame(estimate = Developmental_Exposure_Model_Estimates[,"estimate"], 
                                           lowerCL = Developmental_Exposure_Model_Estimates[,"ci.lb"], 
                                           upperCL = Developmental_Exposure_Model_Estimates[,"ci.ub"], 
                                           K = developmental_exposure_k[,1], 
                                           group_no = developmental_exposure_group_no[,1], 
                                           row.names = developmental_exposure_rnames)
developmental_exposure_table$name <- row.names(developmental_exposure_table)

developmental_exposure_raw_mean <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Embryo") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Juvenile") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Larvae") %>% 
                                                     select("InRR_Transformed"))))

developmental_exposure_raw_name <- c(replicate(70, "Embryo"), 
                                     replicate(13, "Juvenile"), 
                                     replicate(55, "Larva"))

developmental_exposure_raw_df <- data.frame("Model" = developmental_exposure_raw_name, 
                                            "Effect" = developmental_exposure_raw_mean)

# Graph code - Combined

Developmental_Exposure_Order <- c("Larva", "Juvenile", "Embryo")

density_developmental_exposure <- developmental_exposure_table %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = developmental_exposure_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table)[1], 1)), xmin = min(developmental_exposure_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table)[1], 1)), xmin = max(developmental_exposure_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                  scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-0.5, 0.5)) +
                                  annotate('text',  x = 0.5, y = (seq(1, dim(developmental_exposure_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(developmental_exposure_table["Larva", "K"], 
                                                                developmental_exposure_table["Juvenile", "K"], 
                                                                developmental_exposure_table["Embryo", "K"]), "~","(", 
                                                              c(developmental_exposure_table["Larva", "group_no"], 
                                                                developmental_exposure_table["Juvenile", "group_no"], 
                                                                developmental_exposure_table["Embryo", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                         paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.4, y = (seq(1, dim(developmental_exposure_table)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure #(400x320)

# Preparing Graph - Part 1

developmental_exposure_rnames_1 <- c("Embryo", "Juvenile")

developmental_exposure_k_1 <- data.frame("k" = c(Developmental_Exposure_Exploration["Embryo", "Freq"], 
                                                 Developmental_Exposure_Exploration["Juvenile", "Freq"]), 
                                         row.names = developmental_exposure_rnames_1)

developmental_exposure_group_no_1 <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Embryo", "Freq"], 
                                                              Developmental_Exposure_Species_Count["Juvenile", "Freq"]), 
                                                row.names = developmental_exposure_rnames_1)

developmental_exposure_study_1 <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Embryo", "Freq"], 
                                                         Developmental_Exposure_Study_Count["Juvenile", "Freq"]), 
                                             row.names = developmental_exposure_rnames_1)

Developmental_Exposure_Model_Estimates_Reorder_1 <- Developmental_Exposure_Model_Estimates[c("Embryo", "Juvenile"), ]

developmental_exposure_table_1 <- data.frame(estimate = Developmental_Exposure_Model_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Developmental_Exposure_Model_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Developmental_Exposure_Model_Estimates_Reorder_1[,"ci.ub"], 
                                             K = developmental_exposure_k_1[,1], 
                                             group_no = developmental_exposure_group_no_1[,1], 
                                             row.names = developmental_exposure_rnames_1)
developmental_exposure_table_1$name <- row.names(developmental_exposure_table_1)

developmental_exposure_raw_mean_1 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Embryo") %>% 
                                                     select("InRR_Transformed"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Juvenile") %>% 
                                                     select("InRR_Transformed"))))

developmental_exposure_raw_name_1 <- c(replicate(70, "Embryo"), 
                                       replicate(13, "Juvenile"))

developmental_exposure_raw_df_1 <- data.frame("Model" = developmental_exposure_raw_name_1, 
                                              "Effect" = developmental_exposure_raw_mean_1)

# Graph code - Part 1

Developmental_Exposure_Order_1 <- c("Juvenile", "Embryo")

density_developmental_exposure_1 <- developmental_exposure_table_1 %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = developmental_exposure_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)), xmin = min(developmental_exposure_raw_df_1$Effect)-0.04, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)), xmin = max(developmental_exposure_raw_df_1$Effect)+0.04, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                    scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(developmental_exposure_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(developmental_exposure_table_1["Juvenile", "K"], 
                                                                  developmental_exposure_table_1["Embryo", "K"]), "~","(", 
                                                                c(developmental_exposure_table_1["Juvenile", "group_no"], 
                                                                  developmental_exposure_table_1["Embryo", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(developmental_exposure_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure_1 #(400x240)

# Preparing Graph - Part 2

developmental_exposure_rnames_2 <- c("Larva")

developmental_exposure_k_2 <- data.frame("k" = c(Developmental_Exposure_Exploration["Larvae", "Freq"]), 
                                         row.names = developmental_exposure_rnames_2)

developmental_exposure_group_no_2 <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Larvae", "Freq"]), 
                                                row.names = developmental_exposure_rnames_2)

developmental_exposure_study_2 <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Larvae", "Freq"]), 
                                             row.names = developmental_exposure_rnames_2)

Developmental_Exposure_Model_Estimates_Reorder_2 <- Developmental_Exposure_Model_Estimates[c("Larva"), ]

developmental_exposure_table_2 <- data.frame(estimate = Developmental_Exposure_Model_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Developmental_Exposure_Model_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Developmental_Exposure_Model_Estimates_Reorder_2[,"ci.ub"], 
                                             K = developmental_exposure_k_2[,1], 
                                             group_no = developmental_exposure_group_no_2[,1], 
                                             row.names = developmental_exposure_rnames_2)
developmental_exposure_table_2$name <- row.names(developmental_exposure_table_2)

developmental_exposure_raw_mean_2 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Larvae") %>% 
                                                     select("InRR_Transformed"))))

developmental_exposure_raw_name_2 <- c(replicate(55, "Larva"))

developmental_exposure_raw_df_2 <- data.frame("Model" = developmental_exposure_raw_name_2, 
                                              "Effect" = developmental_exposure_raw_mean_2)

# Graph code - Part 2

Developmental_Exposure_Order_2 <- c("Larva")

density_developmental_exposure_2 <- developmental_exposure_table_2 %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = developmental_exposure_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.04, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)), xmin = min(developmental_exposure_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)), xmin = max(developmental_exposure_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                    coord_cartesian(xlim = c(-0.5, 0.5)) +
                                    annotate('text',  x = 0.5, y = (seq(1, dim(developmental_exposure_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(developmental_exposure_table_2["Larva", "K"]), "~","(", 
                                                                c(developmental_exposure_table_2["Larva", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.4, y = (seq(1, dim(developmental_exposure_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure_2 #(400x160)

##### Developmental Subset Model - Class Meta-Regression #####
Developmental_Class_Exploration <- Developmental_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Exploration) <- Developmental_Class_Exploration$Class

Developmental_Class_Data <- Developmental_Subset_Data %>% filter(Class != "Actinopteri" &
                                                                 Class != "Amphibia" &
                                                                 Class != "Anthozoa" &
                                                                 Class != "Branchiopoda")

Developmental_Class_Species_Count <- Developmental_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Species_Count) <- Developmental_Class_Species_Count$Class

Developmental_Class_Study_Count <- Developmental_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>%
  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Study_Count) <- Developmental_Class_Study_Count$Class

Developmental_Class_Species <- Developmental_Class_Data %>% select("phylo") %>% unique()

Developmental_Class_A_cor <- as.data.frame(A_cor)
Developmental_Class_A_cor <- Developmental_Class_A_cor[c(Developmental_Class_Species$phylo), c(Developmental_Class_Species$phylo)]
Developmental_Class_A_cor <- as.matrix(Developmental_Class_A_cor)

Developmental_Class_VCV <- make_VCV_matrix(Developmental_Class_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Class_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_Class_VCV, test = "t", dfs = "contain",
                                                 mods = ~ Class - 1,
                                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                               ~1|Shared_Animal_Number, ~1|Measurement), 
                                                 R = list(phylo=Developmental_Class_A_cor), data = Developmental_Class_Data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
    saveRDS(Developmental_Class_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Class_Model.rds")
  } else {
    Developmental_Class_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Class_Model.rds")})

Developmental_Class_Model_rob <- robust(Developmental_Class_Model, cluster = Developmental_Class_Data$Study_ID, adjust = TRUE)

Developmental_Class_Model_Estimates <- data.frame(Class = substr(row.names(Developmental_Class_Model$b), 6, 100),
                                                  estimate = Developmental_Class_Model$b, 
                                                  ci.lb = Developmental_Class_Model$ci.lb, 
                                                  ci.ub = Developmental_Class_Model$ci.ub)
rownames(Developmental_Class_Model_Estimates) <- Developmental_Class_Model_Estimates$Class
Developmental_Class_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Class_Model), 2))

# Preparing Graph - Combined

developmental_class_rnames <- c("Arachnida", "Insecta")

developmental_class_k <- data.frame("k" = c(Developmental_Class_Exploration["Arachnida", "Freq"], 
                                            Developmental_Class_Exploration["Insecta", "Freq"]), 
                                    row.names = developmental_class_rnames)

developmental_class_group_no <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Arachnida", "Freq"],
                                                         Developmental_Class_Species_Count["Insecta", "Freq"]), 
                                           row.names = developmental_class_rnames)

developmental_class_study <- data.frame("Study" = c(Developmental_Class_Study_Count["Arachnida", "Freq"],
                                                    Developmental_Class_Study_Count["Insecta", "Freq"]), 
                                        row.names = developmental_class_rnames)

developmental_class_table <- data.frame(estimate = Developmental_Class_Model_Estimates[,"estimate"], 
                                        lowerCL = Developmental_Class_Model_Estimates[,"ci.lb"], 
                                        upperCL = Developmental_Class_Model_Estimates[,"ci.ub"], 
                                        K = developmental_class_k[,1], 
                                        group_no = developmental_class_group_no[,1], 
                                        row.names = developmental_class_rnames)
developmental_class_table$name <- row.names(developmental_class_table)

developmental_class_raw_mean <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                                  select("InRR_Transformed"))), 
                                  unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                  select("InRR_Transformed"))))

developmental_class_raw_name <- c(replicate(11, "Arachnida"), 
                                  replicate(73, "Insecta"))

developmental_class_raw_df <- data.frame("Model" = developmental_class_raw_name, 
                                         "Effect" = developmental_class_raw_mean)

# Graph code - Combined

Developmental_Class_Order <- c("Insecta", "Arachnida")

density_developmental_class <- developmental_class_table %>% mutate(name = fct_relevel(name, Developmental_Class_Order)) %>%
                               ggplot() +
                               geom_density_ridges(data = developmental_class_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table)[1], 1)), xmin = min(developmental_class_raw_df$Effect)-0.02, xmax = -1.5, colour = name),
                                              size = 1) +
                               geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table)[1], 1)), xmin = max(developmental_class_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                               coord_cartesian(xlim = c(-0.5, 0.5)) +
                               annotate('text',  x = 0.5, y = (seq(1, dim(developmental_class_table)[1], 1)+0.4),
                               label= paste("italic(k)==", c(developmental_class_table["Insecta", "K"], 
                                                             developmental_class_table["Arachnida", "K"]), "~","(", 
                                                           c(developmental_class_table["Insecta", "group_no"],
                                                             developmental_class_table["Arachnida", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                      paste(format(round(mean(exp(Developmental_Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.4, y = (seq(1, dim(developmental_class_table)[1], 1)+0.4)), size = 3.5)

density_developmental_class #(400x240)

# Preparing Graph - Part 1

developmental_class_rnames_1 <- c("Arachnida")

developmental_class_k_1 <- data.frame("k" = c(Developmental_Class_Exploration["Arachnida", "Freq"]), 
                                      row.names = developmental_class_rnames_1)

developmental_class_group_no_1 <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Arachnida", "Freq"]), 
                                             row.names = developmental_class_rnames_1)

developmental_class_study_1 <- data.frame("Study" = c(Developmental_Class_Study_Count["Arachnida", "Freq"]), 
                                          row.names = developmental_class_rnames_1)

Developmental_Class_Model_Estimates_Reorder_1 <- Developmental_Class_Model_Estimates[c("Arachnida"), ]

developmental_class_table_1 <- data.frame(estimate = Developmental_Class_Model_Estimates_Reorder_1[,"estimate"], 
                                        lowerCL = Developmental_Class_Model_Estimates_Reorder_1[,"ci.lb"], 
                                        upperCL = Developmental_Class_Model_Estimates_Reorder_1[,"ci.ub"], 
                                        K = developmental_class_k_1[,1], 
                                        group_no = developmental_class_group_no_1[,1], 
                                        row.names = developmental_class_rnames_1)
developmental_class_table_1$name <- row.names(developmental_class_table_1)

developmental_class_raw_mean_1 <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                                  select("InRR_Transformed"))))

developmental_class_raw_name_1 <- c(replicate(11, "Arachnida"))

developmental_class_raw_df_1 <- data.frame("Model" = developmental_class_raw_name_1, 
                                           "Effect" = developmental_class_raw_mean_1)

# Graph code - Part 1

Developmental_Class_Order_1 <- c("Arachnida")

density_developmental_class_1 <- developmental_class_table_1 %>% mutate(name = fct_relevel(name, Developmental_Class_Order_1)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = developmental_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order_1)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.05, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_1)[1], 1)), xmin = min(developmental_class_raw_df_1$Effect)-0.05, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_1)[1], 1)), xmin = max(developmental_class_raw_df_1$Effect)+0.06, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(developmental_class_table_1)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(developmental_class_table_1["Arachnida", "K"]), "~","(", 
                                                             c(developmental_class_table_1["Arachnida", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(developmental_class_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_class_1 #(400x160)

# Preparing Graph - Part 2

developmental_class_rnames_2 <- c("Insecta")

developmental_class_k_2 <- data.frame("k" = c(Developmental_Class_Exploration["Insecta", "Freq"]), 
                                      row.names = developmental_class_rnames_2)

developmental_class_group_no_2 <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Insecta", "Freq"]), 
                                             row.names = developmental_class_rnames_2)

developmental_class_study_2 <- data.frame("Study" = c(Developmental_Class_Study_Count["Insecta", "Freq"]), 
                                          row.names = developmental_class_rnames_2)

Developmental_Class_Model_Estimates_Reorder_2 <- Developmental_Class_Model_Estimates[c("Insecta"), ]

developmental_class_table_2 <- data.frame(estimate = Developmental_Class_Model_Estimates_Reorder_2[,"estimate"], 
                                          lowerCL = Developmental_Class_Model_Estimates_Reorder_2[,"ci.lb"], 
                                          upperCL = Developmental_Class_Model_Estimates_Reorder_2[,"ci.ub"], 
                                          K = developmental_class_k_2[,1], 
                                          group_no = developmental_class_group_no_2[,1], 
                                          row.names = developmental_class_rnames_2)
developmental_class_table_2$name <- row.names(developmental_class_table_2)

developmental_class_raw_mean_2 <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                  select("InRR_Transformed"))))

developmental_class_raw_name_2 <- c(replicate(73, "Insecta"))

developmental_class_raw_df_2 <- data.frame("Model" = developmental_class_raw_name_2, 
                                           "Effect" = developmental_class_raw_mean_2)

# Graph code - Part 2

Developmental_Class_Order_2 <- c("Insecta")

density_developmental_class_2 <- developmental_class_table_2 %>% mutate(name = fct_relevel(name, Developmental_Class_Order_2)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = developmental_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order_2)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                         scale = 0.03, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_2)[1], 1)), xmin = min(developmental_class_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                size = 1) +
                                 geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_2)[1], 1)), xmin = max(developmental_class_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                 coord_cartesian(xlim = c(-0.5, 0.5)) +
                                 annotate('text',  x = 0.5, y = (seq(1, dim(developmental_class_table_2)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(developmental_class_table_2["Insecta", "K"]), "~","(", 
                                                             c(developmental_class_table_2["Insecta", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.4, y = (seq(1, dim(developmental_class_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_class_2 #(400x160)

##### Developmental Subset Model - Specific Trait Meta-Regression #####
Developmental_Specific_Trait_Exploration <- Developmental_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Developmental_Specific_Trait_Exploration <- Developmental_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Developmental_Specific_Trait_Exploration) <- Developmental_Specific_Trait_Exploration$Measurement

Developmental_Specific_Trait_Data <- Developmental_Subset_Data %>% filter(Measurement == "Development Time"| 
                                                                          Measurement == "Length"|
                                                                          Measurement == "Mass")

Developmental_Specific_Trait_Species_Count <- Developmental_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Developmental_Specific_Trait_Species_Count) <- Developmental_Specific_Trait_Species_Count$Measurement

Developmental_Specific_Trait_Study_Count <- Developmental_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Developmental_Specific_Trait_Study_Count) <- Developmental_Specific_Trait_Study_Count$Measurement

Developmental_Specific_Trait_Species <- Developmental_Specific_Trait_Data %>% select("phylo") %>% unique()

Developmental_Specific_Trait_A_cor <- as.data.frame(A_cor)
Developmental_Specific_Trait_A_cor <- Developmental_Specific_Trait_A_cor[c(Developmental_Specific_Trait_Species$phylo), c(Developmental_Specific_Trait_Species$phylo)]
Developmental_Specific_Trait_A_cor <- as.matrix(Developmental_Specific_Trait_A_cor)

Developmental_Specific_Trait_VCV <- make_VCV_matrix(Developmental_Specific_Trait_Data, V = "v_InRR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Developmental_Specific_Trait_Model <- metafor::rma.mv(InRR_Transformed, V = Developmental_Specific_Trait_VCV, test = "t", dfs = "contain",
                                                          mods = ~ Measurement - 1,
                                                          random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                        ~1|Shared_Animal_Number), 
                                                          R = list(phylo=Developmental_Specific_Trait_A_cor), data = Developmental_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                          control=list(rel.tol=1e-9))
    saveRDS(Developmental_Specific_Trait_Model, "./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Specific_Trait_Model.rds")
  } else {
    Developmental_Specific_Trait_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Complex_Developmental_Specific_Trait_Model.rds")})

Developmental_Specific_Trait_Model_rob <- robust(Developmental_Specific_Trait_Model, cluster = Developmental_Specific_Trait_Data$Study_ID, adjust = TRUE)

Developmental_Specific_Trait_Model_Estimates <- data.frame(Trait = substr(row.names(Developmental_Specific_Trait_Model$b), 12, 100),
                                                           estimate = Developmental_Specific_Trait_Model$b, ci.lb = Developmental_Specific_Trait_Model$ci.lb, 
                                                           ci.ub = Developmental_Specific_Trait_Model$ci.ub)
rownames(Developmental_Specific_Trait_Model_Estimates) <- Developmental_Specific_Trait_Model_Estimates$Trait
Developmental_Specific_Trait_Model_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Specific_Trait_Model), 2))

# Preparing Graph - Combined

developmental_specific_trait_rnames <- c("Development Time", "Length", "Mass")

developmental_specific_trait_k <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Length", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Mass", "Freq"]), 
                                             row.names = developmental_specific_trait_rnames)

developmental_specific_trait_group_no <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Length", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Mass", "Freq"]), 
                                                    row.names = developmental_specific_trait_rnames)

developmental_specific_trait_study <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Length", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Mass", "Freq"]), 
                                                 row.names = developmental_specific_trait_rnames)

developmental_specific_trait_table <- data.frame(estimate = Developmental_Specific_Trait_Model_Estimates[,"estimate"], 
                                                 lowerCL = Developmental_Specific_Trait_Model_Estimates[,"ci.lb"], 
                                                 upperCL = Developmental_Specific_Trait_Model_Estimates[,"ci.ub"], 
                                                 K = developmental_specific_trait_k[,1], 
                                                 group_no = developmental_specific_trait_group_no[,1], 
                                                 row.names = developmental_specific_trait_rnames)
developmental_specific_trait_table$name <- row.names(developmental_specific_trait_table)

developmental_specific_trait_raw_mean <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                           select("InRR_Transformed"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                           select("InRR_Transformed"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                           select("InRR_Transformed"))))

developmental_specific_trait_raw_name <- c(replicate(46, "Development Time"), 
                                           replicate(14, "Length"), 
                                           replicate(25, "Mass"))

developmental_specific_trait_raw_df <- data.frame("Model" = developmental_specific_trait_raw_name, 
                                                  "Effect" = developmental_specific_trait_raw_mean)

# Graph code - Combined

Developmental_Specific_Trait_Order <- c("Mass", "Length", "Development Time")

density_developmental_specific_trait <- developmental_specific_trait_table %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = developmental_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)), xmin = min(developmental_specific_trait_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)), xmin = max(developmental_specific_trait_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                         vjust = c(-2.7, -2.7, -0.8))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                        scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                        coord_cartesian(xlim = c(-0.5, 0.5)) +
                                        annotate('text',  x = 0.5, y = (seq(1, dim(developmental_specific_trait_table)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(developmental_specific_trait_table["Mass", "K"],
                                                                      developmental_specific_trait_table["Length", "K"],
                                                                      developmental_specific_trait_table["Development Time", "K"]), "~","(", 
                                                                    c(developmental_specific_trait_table["Mass", "group_no"],
                                                                      developmental_specific_trait_table["Length", "group_no"], 
                                                                      developmental_specific_trait_table["Development Time", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                       x = -0.4, y = (seq(1, dim(developmental_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait #(400x360)

# Preparing Graph - Part 1

developmental_specific_trait_rnames_1 <- c("Development Time", "Length")

developmental_specific_trait_k_1 <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Length", "Freq"]), 
                                               row.names = developmental_specific_trait_rnames_1)

developmental_specific_trait_group_no_1 <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                    Developmental_Specific_Trait_Species_Count["Length", "Freq"]), 
                                                      row.names = developmental_specific_trait_rnames_1)

developmental_specific_trait_study_1 <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                               Developmental_Specific_Trait_Study_Count["Length", "Freq"]), 
                                                   row.names = developmental_specific_trait_rnames_1)

Developmental_Specific_Trait_Model_Estimates_Reorder_1 <- Developmental_Specific_Trait_Model_Estimates[c("Development Time", "Length"), ]

developmental_specific_trait_table_1 <- data.frame(estimate = Developmental_Specific_Trait_Model_Estimates_Reorder_1[,"estimate"], 
                                                 lowerCL = Developmental_Specific_Trait_Model_Estimates_Reorder_1[,"ci.lb"], 
                                                 upperCL = Developmental_Specific_Trait_Model_Estimates_Reorder_1[,"ci.ub"], 
                                                 K = developmental_specific_trait_k_1[,1], 
                                                 group_no = developmental_specific_trait_group_no_1[,1], 
                                                 row.names = developmental_specific_trait_rnames_1)
developmental_specific_trait_table_1$name <- row.names(developmental_specific_trait_table_1)

developmental_specific_trait_raw_mean_1 <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                           select("InRR_Transformed"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                           select("InRR_Transformed"))))

developmental_specific_trait_raw_name_1 <- c(replicate(46, "Development Time"), 
                                             replicate(14, "Length"))

developmental_specific_trait_raw_df_1 <- data.frame("Model" = developmental_specific_trait_raw_name_1, 
                                                    "Effect" = developmental_specific_trait_raw_mean_1)

# Graph code - Part 1

Developmental_Specific_Trait_Order_1 <- c("Length", "Development Time")

density_developmental_specific_trait_1 <- developmental_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order_1)) %>%
                                          ggplot() +
                                          geom_density_ridges(data = developmental_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order_1)), 
                                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)), xmin = min(developmental_specific_trait_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)), xmin = max(developmental_specific_trait_raw_df_1$Effect)+0.01, xmax = 1.5, colour = name),
                                                         size = 1) +
                                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                          size = 1, fatten = 2) +
                                          theme_bw() +
                                          guides(fill = "none", colour = "none") +
                                          labs(x = TeX("Effect Size (PRRD)"), y = "") +
                                          theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                           vjust = c(-2.7, -0.8))) +
                                          theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                          theme(axis.ticks = element_blank()) +
                                          theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                          scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                          scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                          coord_cartesian(xlim = c(-0.5, 0.5)) +
                                          annotate('text',  x = 0.5, y = (seq(1, dim(developmental_specific_trait_table_1)[1], 1)+0.4),
                                          label= paste("italic(k)==", c(developmental_specific_trait_table_1["Length", "K"],
                                                                        developmental_specific_trait_table_1["Development Time", "K"]), "~","(", 
                                                                      c(developmental_specific_trait_table_1["Length", "group_no"], 
                                                                        developmental_specific_trait_table_1["Development Time", "group_no"]), 
                                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                          geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.4, y = (seq(1, dim(developmental_specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait_1 #(400x240)

# Preparing Graph - Part 2

developmental_specific_trait_rnames_2 <- c("Mass")

developmental_specific_trait_k_2 <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Mass", "Freq"]), 
                                               row.names = developmental_specific_trait_rnames_2)

developmental_specific_trait_group_no_2 <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Mass", "Freq"]), 
                                                      row.names = developmental_specific_trait_rnames_2)

developmental_specific_trait_study_2 <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Mass", "Freq"]), 
                                                   row.names = developmental_specific_trait_rnames_2)

Developmental_Specific_Trait_Model_Estimates_Reorder_2 <- Developmental_Specific_Trait_Model_Estimates[c("Mass"), ]

developmental_specific_trait_table_2 <- data.frame(estimate = Developmental_Specific_Trait_Model_Estimates_Reorder_2[,"estimate"], 
                                                 lowerCL = Developmental_Specific_Trait_Model_Estimates_Reorder_2[,"ci.lb"], 
                                                 upperCL = Developmental_Specific_Trait_Model_Estimates_Reorder_2[,"ci.ub"], 
                                                 K = developmental_specific_trait_k_2[,1], 
                                                 group_no = developmental_specific_trait_group_no_2[,1], 
                                                 row.names = developmental_specific_trait_rnames_2)
developmental_specific_trait_table_2$name <- row.names(developmental_specific_trait_table_2)

developmental_specific_trait_raw_mean_2 <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                           select("InRR_Transformed"))))

developmental_specific_trait_raw_name_2 <- c(replicate(25, "Mass"))

developmental_specific_trait_raw_df_2 <- data.frame("Model" = developmental_specific_trait_raw_name_2, 
                                                    "Effect" = developmental_specific_trait_raw_mean_2)

# Graph code - Part 2

Developmental_Specific_Trait_Order_2 <- c("Mass")

density_developmental_specific_trait_2 <- developmental_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order_2)) %>%
                                          ggplot() +
                                          geom_density_ridges(data = developmental_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order_2)), 
                                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                  scale = 0.02, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)), xmin = min(developmental_specific_trait_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                                         size = 1) +
                                          geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)), xmin = max(developmental_specific_trait_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                                         size = 1) +
                                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                                          coord_cartesian(xlim = c(-0.5, 0.5)) +
                                          annotate('text',  x = 0.5, y = (seq(1, dim(developmental_specific_trait_table_2)[1], 1)+0.4),
                                          label= paste("italic(k)==", c(developmental_specific_trait_table_2["Mass", "K"]), "~","(", 
                                                                      c(developmental_specific_trait_table_2["Mass", "group_no"]), 
                                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                          geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.4, y = (seq(1, dim(developmental_specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait_2 #(400x160)

##### Meta-analytic Models (Intercept Only) - Subset Graph #####

# Preparing Data - Combined

intercept_rnames <- c("Overall", "Individual-level Traits", "Aquatic", 
                      "Terrestrial", "Acclimation", "Developmental")

intercept_means_df <- data.frame("Mean" = c(Overall_Model_Estimates$estimate, 
                                            Individual_Model_Estimates$estimate, 
                                            Aquatic_Model_Estimates$estimate, 
                                            Terrestrial_Model_Estimates$estimate, 
                                            Acclimation_Model_Estimates$estimate, 
                                            Developmental_Model_Estimates$estimate), 
                                 row.names = intercept_rnames)

intercept_low_df <- data.frame("Low CI" = c(Overall_Model_Estimates$ci.lb, 
                                            Individual_Model_Estimates$ci.lb, 
                                            Aquatic_Model_Estimates$ci.lb, 
                                            Terrestrial_Model_Estimates$ci.lb, 
                                            Acclimation_Model_Estimates$ci.lb, 
                                            Developmental_Model_Estimates$ci.lb), 
                               row.names = intercept_rnames)

intercept_high_df <- data.frame("High CI" = c(Overall_Model_Estimates$ci.ub, 
                                              Individual_Model_Estimates$ci.ub, 
                                              Aquatic_Model_Estimates$ci.ub, 
                                              Terrestrial_Model_Estimates$ci.ub, 
                                              Acclimation_Model_Estimates$ci.ub, 
                                              Developmental_Model_Estimates$ci.ub), 
                                row.names = intercept_rnames)

intercept_k <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                  length(Individual_Subset_Data$Effect_Size_ID), 
                                  length(Aquatic_Subset_Data$Effect_Size_ID), 
                                  length(Terrestrial_Subset_Data$Effect_Size_ID), 
                                  length(Acclimation_Subset_Data$Effect_Size_ID),
                                  length(Developmental_Subset_Data$Effect_Size_ID)), 
                          row.names = intercept_rnames)

intercept_group_no <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                               length(unique(Individual_Subset_Data$Scientific_Name)), 
                                               length(unique(Aquatic_Subset_Data$Scientific_Name)), 
                                               length(unique(Terrestrial_Subset_Data$Scientific_Name)), 
                                               length(unique(Acclimation_Subset_Data$Scientific_Name)),
                                               length(unique(Developmental_Subset_Data$Scientific_Name))),
                                 row.names = intercept_rnames)

intercept_study <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                          length(unique(Individual_Subset_Data$Study_ID)), 
                                          length(unique(Aquatic_Subset_Data$Study_ID)), 
                                          length(unique(Terrestrial_Subset_Data$Study_ID)), 
                                          length(unique(Acclimation_Subset_Data$Study_ID)),
                                          length(unique(Developmental_Subset_Data$Study_ID))),
                              row.names = intercept_rnames)


intercept_table <- data.frame(estimate = intercept_means_df[,1], 
                              lowerCL = intercept_low_df[,1], 
                              upperCL = intercept_high_df[,1], 
                              K = intercept_k[,1], 
                              group_no = intercept_group_no[,1], 
                              row.names = intercept_rnames)
intercept_table$name <- row.names(intercept_table)

intercept_raw_mean <- c(data$InRR_Transformed, 
                        Individual_Subset_Data$InRR_Transformed, 
                        Aquatic_Subset_Data$InRR_Transformed, 
                        Terrestrial_Subset_Data$InRR_Transformed, 
                        Acclimation_Subset_Data$InRR_Transformed, 
                        Developmental_Subset_Data$InRR_Transformed)

intercept_raw_name <- c(replicate(212, "Overall"), 
                        replicate(203, "Individual-level Traits"), 
                        replicate(51, "Aquatic"), 
                        replicate(152, "Terrestrial"), 
                        replicate(65, "Acclimation"), 
                        replicate(138, "Developmental"))

intercept_raw_df <- data.frame("Model" = intercept_raw_name, 
                               "Effect" = intercept_raw_mean)
# Graph code - Combined

Intercept_Order <- c("Developmental", "Acclimation", "Terrestrial", "Aquatic", "Individual-level Traits", "Overall")

density_intercept <- intercept_table %>% mutate(name = fct_relevel(name, Intercept_Order)) %>%
                     ggplot() +
                     geom_density_ridges(data = intercept_raw_df %>% mutate(Model = fct_relevel(Model, Intercept_Order)), 
                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                     geom_linerange(aes(y = rev(seq(1, dim(intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(intercept_table)[1], 1)), xmin = min(intercept_raw_df$Effect)-0.01, xmax = -1.5, colour = name),
                                    size = 1) +
                     geom_linerange(aes(y = rev(seq(1, dim(intercept_table)[1], 1)), xmin = max(intercept_raw_df$Effect)+0.02, xmax = 1.5, colour = name),
                                    size = 1) +
                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                     theme_bw() +
                     guides(fill = "none", colour = "none") +
                     labs(x = TeX("Effect Size (PRRD)"), y = "") +
                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                           vjust = c(-2.7, -2.7, -2.7, -2.7, -0.8, -2.7))) +
                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                     theme(axis.ticks = element_blank()) +
                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                     #theme(axis.text.y = element_text())) +
                     scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                     scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                     coord_cartesian(xlim = c(-0.5, 0.5)) +
                     annotate('text',  x = 0.5, y = (seq(1, dim(intercept_table)[1], 1)+0.4),
                     label= paste("italic(k)==", c(intercept_table["Developmental", "K"], 
                                                   intercept_table["Acclimation", "K"], 
                                                   intercept_table["Terrestrial", "K"], 
                                                   intercept_table["Aquatic", "K"], 
                                                   intercept_table["Individual-level Traits", "K"], 
                                                   intercept_table["Overall", "K"]), "~","(", 
                                                 c(intercept_table["Developmental", "group_no"], 
                                                   intercept_table["Acclimation", "group_no"], 
                                                   intercept_table["Terrestrial", "group_no"], 
                                                   intercept_table["Aquatic", "group_no"], 
                                                   intercept_table["Individual-level Traits", "group_no"], 
                                                   intercept_table["Overall", "group_no"]), 
                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Acclimation_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Terrestrial_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Aquatic_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Individual_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Overall_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                               x = -0.4, y = (seq(1, dim(intercept_table)[1], 1)+0.4)), size = 3.5)
density_intercept #(400x560)

# Preparing Data - Part 1

intercept_rnames_1 <- c("Overall", "Individual-level Traits", "Aquatic")

intercept_means_df_1 <- data.frame("Mean" = c(Overall_Model_Estimates$estimate, 
                                              Individual_Model_Estimates$estimate, 
                                              Aquatic_Model_Estimates$estimate), 
                                   row.names = intercept_rnames_1)

intercept_low_df_1 <- data.frame("Low CI" = c(Overall_Model_Estimates$ci.lb, 
                                              Individual_Model_Estimates$ci.lb, 
                                              Aquatic_Model_Estimates$ci.lb), 
                                 row.names = intercept_rnames_1)

intercept_high_df_1 <- data.frame("High CI" = c(Overall_Model_Estimates$ci.ub, 
                                                Individual_Model_Estimates$ci.ub, 
                                                Aquatic_Model_Estimates$ci.ub), 
                                  row.names = intercept_rnames_1)

intercept_k_1 <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                    length(Individual_Subset_Data$Effect_Size_ID), 
                                    length(Aquatic_Subset_Data$Effect_Size_ID)), 
                            row.names = intercept_rnames_1)

intercept_group_no_1 <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                                 length(unique(Individual_Subset_Data$Scientific_Name)), 
                                                 length(unique(Aquatic_Subset_Data$Scientific_Name))),
                                   row.names = intercept_rnames_1)

intercept_study_1 <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                            length(unique(Individual_Subset_Data$Study_ID)), 
                                            length(unique(Aquatic_Subset_Data$Study_ID))),
                                row.names = intercept_rnames_1)


intercept_table_1 <- data.frame(estimate = intercept_means_df_1[,1], 
                                lowerCL = intercept_low_df_1[,1], 
                                upperCL = intercept_high_df_1[,1], 
                                K = intercept_k_1[,1], 
                                group_no = intercept_group_no_1[,1], 
                                row.names = intercept_rnames_1)
intercept_table_1$name <- row.names(intercept_table_1)

intercept_raw_mean_1 <- c(data$InRR_Transformed, 
                          Individual_Subset_Data$InRR_Transformed, 
                          Aquatic_Subset_Data$InRR_Transformed)

intercept_raw_name_1 <- c(replicate(212, "Overall"), 
                          replicate(203, "Individual-level Traits"), 
                          replicate(51, "Aquatic"))

intercept_raw_df_1 <- data.frame("Model" = intercept_raw_name_1, 
                                 "Effect" = intercept_raw_mean_1)
# Graph code - Part 1

Intercept_Order_1 <- c("Aquatic", "Individual-level Traits", "Overall")

density_intercept_1 <- intercept_table_1 %>% mutate(name = fct_relevel(name, Intercept_Order_1)) %>%
                       ggplot() +
                       geom_density_ridges(data = intercept_raw_df_1 %>% mutate(Model = fct_relevel(Model, Intercept_Order_1)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_1)[1], 1)), xmin = min(intercept_raw_df_1$Effect)-0.01, xmax = -1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_1)[1], 1)), xmin = max(intercept_raw_df_1$Effect)+0.02, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (PRRD)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                        vjust = c(-2.7, -0.8, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       #theme(axis.text.y = element_text())) +
                       scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                       scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                       coord_cartesian(xlim = c(-0.5, 0.5)) +
                       annotate('text',  x = 0.5, y = (seq(1, dim(intercept_table_1)[1], 1)+0.4),
                       label= paste("italic(k)==", c(intercept_table_1["Aquatic", "K"], 
                                                     intercept_table_1["Individual-level Traits", "K"], 
                                                     intercept_table_1["Overall", "K"]), "~","(", 
                                                   c(intercept_table_1["Aquatic", "group_no"], 
                                                     intercept_table_1["Individual-level Traits", "group_no"], 
                                                     intercept_table_1["Overall", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Individual_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Overall_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                                  x = -0.4, y = (seq(1, dim(intercept_table_1)[1], 1)+0.4)), size = 3.5)
density_intercept_1 #(400x320)

# Preparing Data - Part 2

intercept_rnames_2 <- c("Terrestrial", "Acclimation", "Developmental")

intercept_means_df_2 <- data.frame("Mean" = c(Terrestrial_Model_Estimates$estimate, 
                                              Acclimation_Model_Estimates$estimate, 
                                              Developmental_Model_Estimates$estimate), 
                                   row.names = intercept_rnames_2)

intercept_low_df_2 <- data.frame("Low CI" = c(Terrestrial_Model_Estimates$ci.lb, 
                                              Acclimation_Model_Estimates$ci.lb, 
                                              Developmental_Model_Estimates$ci.lb), 
                                 row.names = intercept_rnames_2)

intercept_high_df_2 <- data.frame("High CI" = c(Terrestrial_Model_Estimates$ci.ub, 
                                                Acclimation_Model_Estimates$ci.ub, 
                                                Developmental_Model_Estimates$ci.ub), 
                                  row.names = intercept_rnames_2)

intercept_k_2 <- data.frame("k" = c(length(Terrestrial_Subset_Data$Effect_Size_ID), 
                                    length(Acclimation_Subset_Data$Effect_Size_ID),
                                    length(Developmental_Subset_Data$Effect_Size_ID)), 
                            row.names = intercept_rnames_2)

intercept_group_no_2 <- data.frame("Spp No." = c(length(unique(Terrestrial_Subset_Data$Scientific_Name)), 
                                                 length(unique(Acclimation_Subset_Data$Scientific_Name)),
                                                 length(unique(Developmental_Subset_Data$Scientific_Name))),
                                   row.names = intercept_rnames_2)

intercept_study_2 <- data.frame("Study" = c(length(unique(Terrestrial_Subset_Data$Study_ID)), 
                                            length(unique(Acclimation_Subset_Data$Study_ID)),
                                            length(unique(Developmental_Subset_Data$Study_ID))),
                                row.names = intercept_rnames_2)


intercept_table_2 <- data.frame(estimate = intercept_means_df_2[,1], 
                                lowerCL = intercept_low_df_2[,1], 
                                upperCL = intercept_high_df_2[,1], 
                                K = intercept_k_2[,1], 
                                group_no = intercept_group_no_2[,1], 
                                row.names = intercept_rnames_2)
intercept_table_2$name <- row.names(intercept_table_2)

intercept_raw_mean_2 <- c(Terrestrial_Subset_Data$InRR_Transformed, 
                          Acclimation_Subset_Data$InRR_Transformed, 
                          Developmental_Subset_Data$InRR_Transformed)

intercept_raw_name_2 <- c(replicate(152, "Terrestrial"), 
                          replicate(65, "Acclimation"), 
                          replicate(138, "Developmental"))

intercept_raw_df_2 <- data.frame("Model" = intercept_raw_name_2, 
                                 "Effect" = intercept_raw_mean_2)
# Graph code - Part 2

Intercept_Order_2 <- c("Developmental", "Acclimation", "Terrestrial")

density_intercept_2 <- intercept_table_2 %>% mutate(name = fct_relevel(name, Intercept_Order_2)) %>%
                       ggplot() +
                       geom_density_ridges(data = intercept_raw_df_2 %>% mutate(Model = fct_relevel(Model, Intercept_Order_2)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_2)[1], 1)), xmin = min(intercept_raw_df_2$Effect)-0.01, xmax = -1.5, colour = name),
                                      size = 1) +
                       geom_linerange(aes(y = rev(seq(1, dim(intercept_table_2)[1], 1)), xmin = max(intercept_raw_df_2$Effect)+0.02, xmax = 1.5, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
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
                       #theme(axis.text.y = element_text())) +
                       scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       coord_cartesian(xlim = c(-0.5, 0.5)) +
                       annotate('text',  x = 0.5, y = (seq(1, dim(intercept_table_2)[1], 1)+0.4),
                       label= paste("italic(k)==", c(intercept_table_2["Developmental", "K"], 
                                                     intercept_table_2["Acclimation", "K"], 
                                                     intercept_table_2["Terrestrial", "K"]), "~","(", 
                                                   c(intercept_table_2["Developmental", "group_no"], 
                                                     intercept_table_2["Acclimation", "group_no"], 
                                                     intercept_table_2["Terrestrial", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Acclimation_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Terrestrial_Model_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                                  x = -0.4, y = (seq(1, dim(intercept_table_2)[1], 1)+0.4)), size = 3.5)
density_intercept_2 #(400x320)

##### Supplementary Material Tables #####

# Consistency Changes - Studies, Species and Effect Sizes Counts

Pre_Data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Pre_Data.csv")
Complex_Pre_Data <- Pre_Data %>% filter(Complex_Design == "Yes")

Individual_Pre_Data <- Complex_Pre_Data %>% filter(Trait_Category != "Population")
Developmental_Pre_Data <- Individual_Pre_Data %>% filter(Plasticity_Mechanism == "Developmental Plasticity") %>% 
                          select("Study_ID", "Species_ID", "Treatment_ID", "Trait_ID", "Developmental_Exposure_Time") %>% 
                          dplyr::rename("Treatment_ID_T1" = `Treatment_ID`, 
                                        "Developmental_Exposure_Time_Original" = `Developmental_Exposure_Time`)
Developmental_Pre_Data <- Developmental_Subset_Data %>% 
                          left_join(Developmental_Pre_Data, by = c("Study_ID", "Species_ID", 
                                                                   "Treatment_ID_T1", "Trait_ID"))

Developmental_Exposure_Time_Studies <- Developmental_Pre_Data %>% select("Study_ID", "Developmental_Exposure_Time_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Original") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Studies) <- Developmental_Exposure_Time_Studies$Developmental_Exposure_Time_Original
colnames(Developmental_Exposure_Time_Studies) <- c("Developmental_Exposure_Time_Original", "Study")

Developmental_Exposure_Time_Species <- Developmental_Pre_Data %>% select("Scientific_Name", "Developmental_Exposure_Time_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Original") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Species) <- Developmental_Exposure_Time_Species$Developmental_Exposure_Time_Original
colnames(Developmental_Exposure_Time_Species) <- c("Developmental_Exposure_Time_Original", "Species")

Developmental_Exposure_Time_Effects <- Developmental_Pre_Data %>% select("Developmental_Exposure_Time_Original") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Effects) <- Developmental_Exposure_Time_Effects$Developmental_Exposure_Time_Original
colnames(Developmental_Exposure_Time_Effects) <- c("Developmental_Exposure_Time_Original", "Effect Sizes")

Developmental_Exposure_Time_Final_Counts <- Developmental_Exposure_Time_Studies %>% 
  left_join(Developmental_Exposure_Time_Species, by = "Developmental_Exposure_Time_Original") %>% 
  left_join(Developmental_Exposure_Time_Effects, by = "Developmental_Exposure_Time_Original")

write.csv(Developmental_Exposure_Time_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Developmental_Exposure_Time_Final_Counts.csv", row.names = FALSE)

Acclimation_Pre_Data <- Individual_Pre_Data %>% filter(Plasticity_Mechanism == "Acclimation") %>% 
  select("Study_ID", "Species_ID", "Treatment_ID", "Trait_ID", "Acclimation_Life.History_Stage") %>% 
  dplyr::rename("Treatment_ID_T1" = `Treatment_ID`, 
                "Acclimation_Life.History_Stage_Original" = `Acclimation_Life.History_Stage`)
Acclimation_Pre_Data <- Acclimation_Subset_Data %>% 
  left_join(Acclimation_Pre_Data, by = c("Study_ID", "Species_ID", 
                                           "Treatment_ID_T1", "Trait_ID"))

Acclimation_LH_Stage_Studies <- Acclimation_Pre_Data %>% select("Study_ID", "Acclimation_Life.History_Stage_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Original") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Studies) <- Acclimation_LH_Stage_Studies$Acclimation_Life.History_Stage_Original
colnames(Acclimation_LH_Stage_Studies) <- c("Acclimation_Life.History_Stage_Original", "Study")

Acclimation_LH_Stage_Species <- Acclimation_Pre_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Original") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Species) <- Acclimation_LH_Stage_Species$Acclimation_Life.History_Stage_Original
colnames(Acclimation_LH_Stage_Species) <- c("Acclimation_Life.History_Stage_Original", "Species")

Acclimation_LH_Stage_Effects <- Acclimation_Pre_Data %>% select("Acclimation_Life.History_Stage_Original") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Effects) <- Acclimation_LH_Stage_Effects$Acclimation_Life.History_Stage_Original
colnames(Acclimation_LH_Stage_Effects) <- c("Acclimation_Life.History_Stage_Original", "Effect Sizes")

Acclimation_LH_Stage_Final_Counts <- Acclimation_LH_Stage_Studies %>% 
  left_join(Acclimation_LH_Stage_Species, by = "Acclimation_Life.History_Stage_Original") %>% 
  left_join(Acclimation_LH_Stage_Effects, by = "Acclimation_Life.History_Stage_Original")

write.csv(Acclimation_LH_Stage_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Acclimation_LH_Stage_Final_Counts.csv", row.names = FALSE)

Measurements_Pre_Data <- Complex_Pre_Data %>%
  select("Study_ID", "Species_ID", "Treatment_ID", "Trait_ID", "Measurement") %>% 
  dplyr::rename("Treatment_ID_T1" = `Treatment_ID`, 
                "Measurement_Original" = `Measurement`)
Measurements_Pre_Data <- data %>% 
  left_join(Measurements_Pre_Data, by = c("Study_ID", "Species_ID", 
                                          "Treatment_ID_T1", "Trait_ID"))

Measurement_Studies <- Measurements_Pre_Data %>% select("Study_ID", "Measurement_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement_Original") %>% table() %>% data.frame()
rownames(Measurement_Studies) <- Measurement_Studies$Measurement_Original
colnames(Measurement_Studies) <- c("Measurement_Original", "Study")

Measurement_Species <- Measurements_Pre_Data %>% select("Scientific_Name", "Measurement_Original") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Measurement_Original") %>% table() %>% data.frame()
rownames(Measurement_Species) <- Measurement_Species$Measurement_Original
colnames(Measurement_Species) <- c("Measurement_Original", "Species")

Measurement_Effects <- Measurements_Pre_Data %>% select("Measurement_Original") %>% table() %>% data.frame()
rownames(Measurement_Effects) <- Measurement_Effects$Measurement_Original
colnames(Measurement_Effects) <- c("Measurement_Original", "Effect Sizes")

Measurement_Final_Counts <- Measurement_Studies %>% 
  left_join(Measurement_Species, by = "Measurement_Original") %>% 
  left_join(Measurement_Effects, by = "Measurement_Original")

write.csv(Measurement_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Measurement_Final_Counts.csv", row.names = FALSE)

# Category - Studies, Species and Effect Sizes

Developmental_Exposure_Time_Category_Studies <- Developmental_Subset_Data %>% select("Study_ID", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Studies) <- Developmental_Exposure_Time_Category_Studies$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Studies) <- c("Developmental_Exposure_Time_Category", "Study")

Developmental_Exposure_Time_Category_Species <- Developmental_Subset_Data %>% select("Scientific_Name", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Species) <- Developmental_Exposure_Time_Category_Species$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Species) <- c("Developmental_Exposure_Time_Category", "Species")

Developmental_Exposure_Time_Category_Effects <- Developmental_Subset_Data %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Effects) <- Developmental_Exposure_Time_Category_Effects$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Effects) <- c("Developmental_Exposure_Time_Category", "Effect Sizes")

Developmental_Exposure_Time_Category_Final_Counts <- Developmental_Exposure_Time_Category_Studies %>% 
  left_join(Developmental_Exposure_Time_Category_Species, by = "Developmental_Exposure_Time_Category") %>% 
  left_join(Developmental_Exposure_Time_Category_Effects, by = "Developmental_Exposure_Time_Category")

write.csv(Developmental_Exposure_Time_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Developmental_Exposure_Time_Category_Final_Counts.csv", row.names = FALSE)

Acclimation_LH_Stage_Category_Studies <- Acclimation_Subset_Data %>% select("Study_ID", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Studies) <- Acclimation_LH_Stage_Category_Studies$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Studies) <- c("Acclimation_Life.History_Stage_Category", "Study")

Acclimation_LH_Stage_Category_Species <- Acclimation_Subset_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Species) <- Acclimation_LH_Stage_Category_Species$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Species) <- c("Acclimation_Life.History_Stage_Category", "Species")

Acclimation_LH_Stage_Category_Effects <- Acclimation_Subset_Data %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Effects) <- Acclimation_LH_Stage_Category_Effects$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Effects) <- c("Acclimation_Life.History_Stage_Category", "Effect Sizes")

Acclimation_LH_Stage_Category_Final_Counts <- Acclimation_LH_Stage_Category_Studies %>% 
  left_join(Acclimation_LH_Stage_Category_Species, by = "Acclimation_Life.History_Stage_Category") %>% 
  left_join(Acclimation_LH_Stage_Category_Effects, by = "Acclimation_Life.History_Stage_Category")

write.csv(Acclimation_LH_Stage_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Acclimation_LH_Stage_Category_Final_Counts.csv", row.names = FALSE)

Measurement_Category_Studies <- data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Studies) <- Measurement_Category_Studies$Trait_Category
colnames(Measurement_Category_Studies) <- c("Trait_Category", "Study")

Measurement_Category_Species <- data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
  filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Species) <- Measurement_Category_Species$Trait_Category
colnames(Measurement_Category_Species) <- c("Trait_Category", "Species")

Measurement_Category_Effects <- data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Effects) <- Measurement_Category_Effects$Trait_Category
colnames(Measurement_Category_Effects) <- c("Trait_Category", "Effect Sizes")

Measurement_Category_Final_Counts <- Measurement_Category_Studies %>% 
  left_join(Measurement_Category_Species, by = "Trait_Category") %>% 
  left_join(Measurement_Category_Effects, by = "Trait_Category")

write.csv(Measurement_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Measurement_Category_Final_Counts.csv", row.names = FALSE)

# Phylogenetic Tree with labels

dev.off()
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

# Raw Data Table

Raw_Overall <- data.frame("Overall" = c("MLMA"),
                          "Studies" = c(intercept_study["Overall", "Study"]), 
                          "Species" = c(intercept_table["Overall", "group_no"]), 
                          "Effect Sizes" = c(intercept_table["Overall", "K"]),
                          "Estimate" = c(Overall_Model$b[1]),
                          "CI Low" = c(Overall_Model$ci.lb), 
                          "CI High" = c(Overall_Model$ci.ub), 
                          "df" = c(Overall_Model$ddf[[1]]), 
                          "p-value" = c(Overall_Model$pval))

Raw_Amplitude <- data.frame("Overall" = c("Fluctuation Amplitude"),
                            "Studies" = c(length(unique(data$Study_ID))), 
                            "Species" = c(length(unique(data$Scientific_Name))), 
                            "Effect Sizes" = c(length(data$Effect_Size_ID)),
                            "Estimate" = c(Amplitude_Model$b[1]),
                            "CI Low" = c(Amplitude_Model$ci.lb), 
                            "CI High" = c(Amplitude_Model$ci.ub), 
                            "df" = c(Amplitude_Model$ddf[[1]]), 
                            "p-value" = c(Amplitude_Model$pval))

Raw_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                          "Stepwise"),
                                   "Studies" = c(fluctuation_study["Sinusoidal (Sine Curve)", "Study"], fluctuation_study["Alternating", "Study"], 
                                                 fluctuation_study["Stepwise", "Study"]), 
                                   "Species" = c(fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], fluctuation_table["Alternating", "group_no"], 
                                                 fluctuation_table["Stepwise", "group_no"]), 
                                   "Effect Sizes" = c(fluctuation_table["Sinusoidal (Sine Curve)", "K"], fluctuation_table["Alternating", "K"], 
                                                      fluctuation_table["Stepwise", "K"]),
                                   "Estimate" = c(Fluctuation_Model$b[[2]], Fluctuation_Model$b[[1]],
                                                  Fluctuation_Model$b[[3]]),
                                   "CI Low" = c(Fluctuation_Model$ci.lb[2], Fluctuation_Model$ci.lb[1], 
                                                Fluctuation_Model$ci.lb[3]), 
                                   "CI High" = c(Fluctuation_Model$ci.ub[2], Fluctuation_Model$ci.ub[1], 
                                                 Fluctuation_Model$ci.ub[3]), 
                                   "df" = c(Fluctuation_Model$ddf[[2]], Fluctuation_Model$ddf[[1]], 
                                            Fluctuation_Model$ddf[[3]]), 
                                   "p-value" = c(Fluctuation_Model$pval[2], Fluctuation_Model$pval[1], 
                                                 Fluctuation_Model$pval[3]))

Raw_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay", "Life-history Traits",  
                                                          "Morphology", "Physiological"),
                        "Studies" = c(trait_study["Biochemical Assay", "Study"], trait_study["Life-history Traits", "Study"],
                                      trait_study["Morphology", "Study"], trait_study["Physiological", "Study"]), 
                        "Species" = c(trait_table["Biochemical Assay", "group_no"], trait_table["Life-history Traits", "group_no"],
                                      trait_table["Morphology", "group_no"], trait_table["Physiological", "group_no"]), 
                        "Effect Sizes" = c(trait_table["Biochemical Assay", "K"], trait_table["Life-history Traits", "K"],
                                           trait_table["Morphology", "K"], trait_table["Physiological", "K"]),
                        "Estimate" = c(Trait_Model$b[[1]], Trait_Model$b[[2]], Trait_Model$b[[3]], Trait_Model$b[[4]]),
                        "CI Low" = c(Trait_Model$ci.lb[1], Trait_Model$ci.lb[2], Trait_Model$ci.lb[3], Trait_Model$ci.lb[4]), 
                        "CI High" = c(Trait_Model$ci.ub[1], Trait_Model$ci.ub[2], Trait_Model$ci.ub[3], Trait_Model$ci.ub[4]), 
                        "df" = c(Trait_Model$ddf[[1]], Trait_Model$ddf[[2]], Trait_Model$ddf[[3]], Trait_Model$ddf[[4]]), 
                        "p-value" = c(Trait_Model$pval[1], Trait_Model$pval[2], Trait_Model$pval[3], Trait_Model$pval[4]))

Raw_Class <- data.frame("Taxonomic Class" = c("Arachnida", "Insecta", "Malacostraca"),
                        "Studies" = c(class_study["Arachnida", "Study"], class_study["Insecta", "Study"], 
                                      class_study["Malacostraca", "Study"]), 
                        "Species" = c(class_table["Arachnida", "group_no"], class_table["Insecta", "group_no"], 
                                      class_table["Malacostraca", "group_no"]), 
                        "Effect Sizes" = c(class_table["Arachnida", "K"], class_table["Insecta", "K"], 
                                           class_table["Malacostraca", "K"]),
                        "Estimate" = c(Class_Model$b[[1]], Class_Model$b[[2]], Class_Model$b[[3]]),
                        "CI Low" = c(Class_Model$ci.lb[1], Class_Model$ci.lb[2], Class_Model$ci.lb[3]), 
                        "CI High" = c(Class_Model$ci.ub[1], Class_Model$ci.ub[2], Class_Model$ci.ub[3]), 
                        "df" = c(Class_Model$ddf[[1]], Class_Model$ddf[[2]], Class_Model$ddf[[3]]), 
                        "p-value" = c(Class_Model$pval[1], Class_Model$pval[2], Class_Model$pval[3]))

Raw_Specific_Trait <- data.frame("Specific Phenotypic Traits" = c("Development Time", "Length", "Mass", "Metabolic Rate"),
                                 "Studies" = c(specific_trait_study["Development Time", "Study"], specific_trait_study["Length", "Study"], 
                                               specific_trait_study["Mass", "Study"], specific_trait_study["Metabolic Rate", "Study"]), 
                                 "Species" = c(specific_trait_table["Development Time", "group_no"], specific_trait_table["Length", "group_no"], 
                                               specific_trait_table["Mass", "group_no"], specific_trait_table["Metabolic Rate", "group_no"]), 
                                 "Effect Sizes" = c(specific_trait_table["Development Time", "K"], specific_trait_table["Length", "K"], 
                                                    specific_trait_table["Mass", "K"], specific_trait_table["Metabolic Rate", "K"]),
                                 "Estimate" = c(Specific_Trait_Model$b[[1]], Specific_Trait_Model$b[[2]], 
                                                Specific_Trait_Model$b[[3]], Specific_Trait_Model$b[[4]]),
                                 "CI Low" = c(Specific_Trait_Model$ci.lb[1], Specific_Trait_Model$ci.lb[2], 
                                              Specific_Trait_Model$ci.lb[3], Specific_Trait_Model$ci.lb[4]), 
                                 "CI High" = c(Specific_Trait_Model$ci.ub[1], Specific_Trait_Model$ci.ub[2], 
                                               Specific_Trait_Model$ci.ub[3], Specific_Trait_Model$ci.ub[4]), 
                                 "df" = c(Specific_Trait_Model$ddf[[1]], Specific_Trait_Model$ddf[[2]], 
                                          Specific_Trait_Model$ddf[[3]], Specific_Trait_Model$ddf[[4]]), 
                                 "p-value" = c(Specific_Trait_Model$pval[1], Specific_Trait_Model$pval[2], 
                                               Specific_Trait_Model$pval[3], Specific_Trait_Model$pval[4]))

Raw_Individual <- data.frame("Individual-level Traits" = c("MLMA"),
                             "Studies" = c(intercept_study["Individual-level Traits", "Study"]), 
                             "Species" = c(intercept_table["Individual-level Traits", "group_no"]), 
                             "Effect Sizes" = c(intercept_table["Individual-level Traits", "K"]),
                             "Estimate" = c(Individual_Model$b[1]),
                             "CI Low" = c(Individual_Model$ci.lb), 
                             "CI High" = c(Individual_Model$ci.ub), 
                             "df" = c(Individual_Model$ddf[[1]]), 
                             "p-value" = c(Individual_Model$pval))

Raw_Individual_Amplitude <- data.frame("Individual-level Traits" = c("Fluctuation Amplitude"),
                                       "Studies" = c(length(unique(Individual_Subset_Data$Study_ID))), 
                                       "Species" = c(length(unique(Individual_Subset_Data$Scientific_Name))), 
                                       "Effect Sizes" = c(length(Individual_Subset_Data$Effect_Size_ID)),
                                       "Estimate" = c(Individual_Amplitude_Model$b[1]),
                                       "CI Low" = c(Individual_Amplitude_Model$ci.lb), 
                                       "CI High" = c(Individual_Amplitude_Model$ci.ub), 
                                       "df" = c(Individual_Amplitude_Model$ddf[[1]]), 
                                       "p-value" = c(Individual_Amplitude_Model$pval))

Raw_Individual_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"),
                                              "Studies" = c(individual_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], individual_fluctuation_study["Alternating", "Study"], 
                                                            individual_fluctuation_study["Stepwise", "Study"]), 
                                              "Species" = c(individual_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], individual_fluctuation_table["Alternating", "group_no"], 
                                                            individual_fluctuation_table["Stepwise", "group_no"]), 
                                              "Effect Sizes" = c(individual_fluctuation_table["Sinusoidal (Sine Curve)", "K"], individual_fluctuation_table["Alternating", "K"], 
                                                                 individual_fluctuation_table["Stepwise", "K"]),
                                              "Estimate" = c(Individual_Fluctuation_Model$b[[2]], Individual_Fluctuation_Model$b[[1]],
                                                             Individual_Fluctuation_Model$b[[3]]),
                                              "CI Low" = c(Individual_Fluctuation_Model$ci.lb[2], Individual_Fluctuation_Model$ci.lb[1], 
                                                           Individual_Fluctuation_Model$ci.lb[3]), 
                                              "CI High" = c(Individual_Fluctuation_Model$ci.ub[2], Individual_Fluctuation_Model$ci.ub[1], 
                                                            Individual_Fluctuation_Model$ci.ub[3]), 
                                              "df" = c(Individual_Fluctuation_Model$ddf[[2]], Individual_Fluctuation_Model$ddf[[1]], 
                                                       Individual_Fluctuation_Model$ddf[[3]]), 
                                              "p-value" = c(Individual_Fluctuation_Model$pval[2], Individual_Fluctuation_Model$pval[1], 
                                                            Individual_Fluctuation_Model$pval[3]))

Raw_Individual_Class <- data.frame("Taxonomic Class" = c("Arachnida", "Insecta"),
                                   "Studies" = c(individual_class_study["Arachnida", "Study"], individual_class_study["Insecta", "Study"]), 
                                   "Species" = c(individual_class_table["Arachnida", "group_no"], individual_class_table["Insecta", "group_no"]), 
                                   "Effect Sizes" = c(individual_class_table["Arachnida", "K"], individual_class_table["Insecta", "K"]),
                                   "Estimate" = c(Individual_Class_Model$b[[1]], Individual_Class_Model$b[[2]]),
                                   "CI Low" = c(Individual_Class_Model$ci.lb[1], Individual_Class_Model$ci.lb[2]), 
                                   "CI High" = c(Individual_Class_Model$ci.ub[1], Individual_Class_Model$ci.ub[2]), 
                                   "df" = c(Individual_Class_Model$ddf[[1]], Individual_Class_Model$ddf[[2]]), 
                                   "p-value" = c(Individual_Class_Model$pval[1], Individual_Class_Model$pval[2]))

Raw_Aquatic <- data.frame("Aquatic Organisms" = c("MLMA"),
                          "Studies" = c(intercept_study["Aquatic", "Study"]), 
                          "Species" = c(intercept_table["Aquatic", "group_no"]), 
                          "Effect Sizes" = c(intercept_table["Aquatic", "K"]),
                          "Estimate" = c(Aquatic_Model$b[1]),
                          "CI Low" = c(Aquatic_Model$ci.lb), 
                          "CI High" = c(Aquatic_Model$ci.ub), 
                          "df" = c(Aquatic_Model$ddf[[1]]), 
                          "p-value" = c(Aquatic_Model$pval))

Raw_Aquatic_Amplitude <- data.frame("Aquatic Organisms" = c("Fluctuation Amplitude"),
                                    "Studies" = c(length(unique(Aquatic_Subset_Data$Study_ID))), 
                                    "Species" = c(length(unique(Aquatic_Subset_Data$Scientific_Name))), 
                                    "Effect Sizes" = c(length(Aquatic_Subset_Data$Effect_Size_ID)),
                                    "Estimate" = c(Aquatic_Amplitude_Model$b[1]),
                                    "CI Low" = c(Aquatic_Amplitude_Model$ci.lb), 
                                    "CI High" = c(Aquatic_Amplitude_Model$ci.ub), 
                                    "df" = c(Aquatic_Amplitude_Model$ddf[[1]]), 
                                    "p-value" = c(Aquatic_Amplitude_Model$pval))

Raw_Aquatic_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating"),
                                           "Studies" = c(aquatic_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], aquatic_fluctuation_study["Alternating", "Study"]), 
                                           "Species" = c(aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], aquatic_fluctuation_table["Alternating", "group_no"]), 
                                           "Effect Sizes" = c(aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "K"], aquatic_fluctuation_table["Alternating", "K"]),
                                           "Estimate" = c(Aquatic_Fluctuation_Model$b[[2]], Aquatic_Fluctuation_Model$b[[1]]),
                                           "CI Low" = c(Aquatic_Fluctuation_Model$ci.lb[2], Aquatic_Fluctuation_Model$ci.lb[1]), 
                                           "CI High" = c(Aquatic_Fluctuation_Model$ci.ub[2], Aquatic_Fluctuation_Model$ci.ub[1]), 
                                           "df" = c(Aquatic_Fluctuation_Model$ddf[[2]], Aquatic_Fluctuation_Model$ddf[[1]]), 
                                           "p-value" = c(Aquatic_Fluctuation_Model$pval[2], Aquatic_Fluctuation_Model$pval[1]))

Raw_Aquatic_Trait <- data.frame("Phenotypic Trait Categories" = c("Life-history Traits", "Physiological"),
                                "Studies" = c(aquatic_trait_study["Life-history Traits", "Study"], aquatic_trait_study["Physiological", "Study"]), 
                                "Species" = c(aquatic_trait_table["Life-history Traits", "group_no"], aquatic_trait_table["Physiological", "group_no"]), 
                                "Effect Sizes" = c(aquatic_trait_table["Life-history Traits", "K"], aquatic_trait_table["Physiological", "K"]),
                                "Estimate" = c(Aquatic_Trait_Model$b[[1]], Aquatic_Trait_Model$b[[2]]),
                                "CI Low" = c(Aquatic_Trait_Model$ci.lb[1], Aquatic_Trait_Model$ci.lb[2]), 
                                "CI High" = c(Aquatic_Trait_Model$ci.ub[1], Aquatic_Trait_Model$ci.ub[2]), 
                                "df" = c(Aquatic_Trait_Model$ddf[[1]], Aquatic_Trait_Model$ddf[[2]]), 
                                "p-value" = c(Aquatic_Trait_Model$pval[1], Aquatic_Trait_Model$pval[2]))

Raw_Aquatic_Plasticity <- data.frame("Exposure Type" = c("Acclimation", "Developmental"),
                                     "Studies" = c(aquatic_plasticity_study["Acclimation", "Study"], aquatic_plasticity_study["Development", "Study"]), 
                                     "Species" = c(aquatic_plasticity_table["Acclimation", "group_no"], aquatic_plasticity_table["Development", "group_no"]), 
                                     "Effect Sizes" = c(aquatic_plasticity_table["Acclimation", "K"], aquatic_plasticity_table["Development", "K"]),
                                     "Estimate" = c(Aquatic_Plasticity_Model$b[[1]], Aquatic_Plasticity_Model$b[[2]]),
                                     "CI Low" = c(Aquatic_Plasticity_Model$ci.lb[1], Aquatic_Plasticity_Model$ci.lb[2]), 
                                     "CI High" = c(Aquatic_Plasticity_Model$ci.ub[1], Aquatic_Plasticity_Model$ci.ub[2]), 
                                     "df" = c(Aquatic_Plasticity_Model$ddf[[1]], Aquatic_Plasticity_Model$ddf[[2]]), 
                                     "p-value" = c(Aquatic_Plasticity_Model$pval[1], Aquatic_Plasticity_Model$pval[2]))

Raw_Terrestrial <- data.frame("Terrestrial Organisms" = c("MLMA"),
                              "Studies" = c(intercept_study["Terrestrial", "Study"]), 
                              "Species" = c(intercept_table["Terrestrial", "group_no"]), 
                              "Effect Sizes" = c(intercept_table["Terrestrial", "K"]),
                              "Estimate" = c(Terrestrial_Model$b[1]),
                              "CI Low" = c(Terrestrial_Model$ci.lb), 
                              "CI High" = c(Terrestrial_Model$ci.ub), 
                              "df" = c(Terrestrial_Model$ddf[[1]]), 
                              "p-value" = c(Terrestrial_Model$pval))

Raw_Terrestrial_Amplitude <- data.frame("Terrestrial Organisms" = c("Fluctuation Amplitude"),
                                        "Studies" = c(length(unique(Terrestrial_Subset_Data$Study_ID))), 
                                        "Species" = c(length(unique(Terrestrial_Subset_Data$Scientific_Name))), 
                                        "Effect Sizes" = c(length(Terrestrial_Subset_Data$Effect_Size_ID)),
                                        "Estimate" = c(Terrestrial_Amplitude_Model$b[1]),
                                        "CI Low" = c(Terrestrial_Amplitude_Model$ci.lb), 
                                        "CI High" = c(Terrestrial_Amplitude_Model$ci.ub), 
                                        "df" = c(Terrestrial_Amplitude_Model$ddf[[1]]), 
                                        "p-value" = c(Terrestrial_Amplitude_Model$pval))

Raw_Terrestrial_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                      "Stepwise"),
                                               "Studies" = c(terrestrial_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], terrestrial_fluctuation_study["Alternating", "Study"], 
                                                             terrestrial_fluctuation_study["Stepwise", "Study"]), 
                                               "Species" = c(terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], terrestrial_fluctuation_table["Alternating", "group_no"], 
                                                             terrestrial_fluctuation_table["Stepwise", "group_no"]), 
                                               "Effect Sizes" = c(terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "K"], terrestrial_fluctuation_table["Alternating", "K"], 
                                                                  terrestrial_fluctuation_table["Stepwise", "K"]),
                                               "Estimate" = c(Terrestrial_Fluctuation_Model$b[[2]], Terrestrial_Fluctuation_Model$b[[1]],
                                                              Terrestrial_Fluctuation_Model$b[[3]]),
                                               "CI Low" = c(Terrestrial_Fluctuation_Model$ci.lb[2], Terrestrial_Fluctuation_Model$ci.lb[1], 
                                                            Terrestrial_Fluctuation_Model$ci.lb[3]), 
                                               "CI High" = c(Terrestrial_Fluctuation_Model$ci.ub[2], Terrestrial_Fluctuation_Model$ci.ub[1], 
                                                             Terrestrial_Fluctuation_Model$ci.ub[3]), 
                                               "df" = c(Terrestrial_Fluctuation_Model$ddf[[2]], Terrestrial_Fluctuation_Model$ddf[[1]], 
                                                        Terrestrial_Fluctuation_Model$ddf[[3]]), 
                                               "p-value" = c(Terrestrial_Fluctuation_Model$pval[2], Terrestrial_Fluctuation_Model$pval[1], 
                                                             Terrestrial_Fluctuation_Model$pval[3]))

Raw_Terrestrial_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay", "Life-history Traits",  
                                                                      "Morphology", "Physiological"),
                                    "Studies" = c(terrestrial_trait_study["Biochemical Assay", "Study"], terrestrial_trait_study["Life-history Traits", "Study"],
                                                  terrestrial_trait_study["Morphology", "Study"], terrestrial_trait_study["Physiological", "Study"]), 
                                    "Species" = c(terrestrial_trait_table["Biochemical Assay", "group_no"], terrestrial_trait_table["Life-history Traits", "group_no"],
                                                  terrestrial_trait_table["Morphology", "group_no"], terrestrial_trait_table["Physiological", "group_no"]), 
                                    "Effect Sizes" = c(terrestrial_trait_table["Biochemical Assay", "K"], terrestrial_trait_table["Life-history Traits", "K"],
                                                       terrestrial_trait_table["Morphology", "K"], terrestrial_trait_table["Physiological", "K"]),
                                    "Estimate" = c(Terrestrial_Trait_Model$b[[1]], Terrestrial_Trait_Model$b[[2]], Terrestrial_Trait_Model$b[[3]], Terrestrial_Trait_Model$b[[4]]),
                                    "CI Low" = c(Terrestrial_Trait_Model$ci.lb[1], Terrestrial_Trait_Model$ci.lb[2], Terrestrial_Trait_Model$ci.lb[3], Terrestrial_Trait_Model$ci.lb[4]), 
                                    "CI High" = c(Terrestrial_Trait_Model$ci.ub[1], Terrestrial_Trait_Model$ci.ub[2], Terrestrial_Trait_Model$ci.ub[3], Terrestrial_Trait_Model$ci.ub[4]), 
                                    "df" = c(Terrestrial_Trait_Model$ddf[[1]], Terrestrial_Trait_Model$ddf[[2]], Terrestrial_Trait_Model$ddf[[3]], Terrestrial_Trait_Model$ddf[[4]]), 
                                    "p-value" = c(Terrestrial_Trait_Model$pval[1], Terrestrial_Trait_Model$pval[2], Terrestrial_Trait_Model$pval[3], Terrestrial_Trait_Model$pval[4]))

Raw_Terrestrial_Plasticity <- data.frame("Exposure Type" = c("Acclimation", "Developmental"),
                                         "Studies" = c(terrestrial_plasticity_study["Acclimation", "Study"], terrestrial_plasticity_study["Development", "Study"]), 
                                         "Species" = c(terrestrial_plasticity_table["Acclimation", "group_no"], terrestrial_plasticity_table["Development", "group_no"]), 
                                         "Effect Sizes" = c(terrestrial_plasticity_table["Acclimation", "K"], terrestrial_plasticity_table["Development", "K"]),
                                         "Estimate" = c(Terrestrial_Plasticity_Model$b[[1]], Terrestrial_Plasticity_Model$b[[2]]),
                                         "CI Low" = c(Terrestrial_Plasticity_Model$ci.lb[1], Terrestrial_Plasticity_Model$ci.lb[2]), 
                                         "CI High" = c(Terrestrial_Plasticity_Model$ci.ub[1], Terrestrial_Plasticity_Model$ci.ub[2]), 
                                         "df" = c(Terrestrial_Plasticity_Model$ddf[[1]], Terrestrial_Plasticity_Model$ddf[[2]]), 
                                         "p-value" = c(Terrestrial_Plasticity_Model$pval[1], Terrestrial_Plasticity_Model$pval[2]))

Raw_Terrestrial_Specific_Trait <- data.frame("Specific Phenotypic Traits" = c("Development Time", "Length", "Mass"),
                                             "Studies" = c(terrestrial_specific_trait_study["Development Time", "Study"], terrestrial_specific_trait_study["Length", "Study"], 
                                                           terrestrial_specific_trait_study["Mass", "Study"]), 
                                             "Species" = c(terrestrial_specific_trait_table["Development Time", "group_no"], terrestrial_specific_trait_table["Length", "group_no"], 
                                                           terrestrial_specific_trait_table["Mass", "group_no"]), 
                                             "Effect Sizes" = c(terrestrial_specific_trait_table["Development Time", "K"], terrestrial_specific_trait_table["Length", "K"], 
                                                                terrestrial_specific_trait_table["Mass", "K"]),
                                             "Estimate" = c(Terrestrial_Specific_Trait_Model$b[[1]], Terrestrial_Specific_Trait_Model$b[[2]], Terrestrial_Specific_Trait_Model$b[[3]]),
                                             "CI Low" = c(Terrestrial_Specific_Trait_Model$ci.lb[1], Terrestrial_Specific_Trait_Model$ci.lb[2], Terrestrial_Specific_Trait_Model$ci.lb[3]), 
                                             "CI High" = c(Terrestrial_Specific_Trait_Model$ci.ub[1], Terrestrial_Specific_Trait_Model$ci.ub[2], Terrestrial_Specific_Trait_Model$ci.ub[3]), 
                                             "df" = c(Terrestrial_Specific_Trait_Model$ddf[[1]], Terrestrial_Specific_Trait_Model$ddf[[2]], Terrestrial_Specific_Trait_Model$ddf[[3]]), 
                                             "p-value" = c(Terrestrial_Specific_Trait_Model$pval[1], Terrestrial_Specific_Trait_Model$pval[2], Terrestrial_Specific_Trait_Model$pval[3]))

Raw_Acclimation <- data.frame("Acclimation Exposure" = c("MLMA"),
                              "Studies" = c(intercept_study["Acclimation", "Study"]), 
                              "Species" = c(intercept_table["Acclimation", "group_no"]), 
                              "Effect Sizes" = c(intercept_table["Acclimation", "K"]),
                              "Estimate" = c(Acclimation_Model$b[1]),
                              "CI Low" = c(Acclimation_Model$ci.lb), 
                              "CI High" = c(Acclimation_Model$ci.ub), 
                              "df" = c(Acclimation_Model$ddf[[1]]), 
                              "p-value" = c(Acclimation_Model$pval))

Raw_Acclimation_Amplitude <- data.frame("Acclimation Exposure" = c("Fluctuation Amplitude"),
                                        "Studies" = c(length(unique(Acclimation_Subset_Data$Study_ID))), 
                                        "Species" = c(length(unique(Acclimation_Subset_Data$Scientific_Name))), 
                                        "Effect Sizes" = c(length(Acclimation_Subset_Data$Effect_Size_ID)),
                                        "Estimate" = c(Acclimation_Amplitude_Model$b[1]),
                                        "CI Low" = c(Acclimation_Amplitude_Model$ci.lb), 
                                        "CI High" = c(Acclimation_Amplitude_Model$ci.ub), 
                                        "df" = c(Acclimation_Amplitude_Model$ddf[[1]]), 
                                        "p-value" = c(Acclimation_Amplitude_Model$pval))

Raw_Acclimation_Exposure <- data.frame("Acclimation Exposure" = c("Exposure Time"),
                                       "Studies" = c(length(unique(Acclimation_Subset_Data$Study_ID))), 
                                       "Species" = c(length(unique(Acclimation_Subset_Data$Scientific_Name))), 
                                       "Effect Sizes" = c(length(Acclimation_Subset_Data$Effect_Size_ID)),
                                       "Estimate" = c(Acclimation_Exposure_Model$b[1]),
                                       "CI Low" = c(Acclimation_Exposure_Model$ci.lb), 
                                       "CI High" = c(Acclimation_Exposure_Model$ci.ub), 
                                       "df" = c(Acclimation_Exposure_Model$ddf[[1]]), 
                                       "p-value" = c(Acclimation_Exposure_Model$pval))

Raw_Acclimation_Frequency <- data.frame("Acclimation Exposure" = c("Number of Fluctuations"),
                                        "Studies" = c(length(unique(Acclimation_Frequency_Data$Study_ID))), 
                                        "Species" = c(length(unique(Acclimation_Frequency_Data$Scientific_Name))), 
                                        "Effect Sizes" = c(length(Acclimation_Frequency_Data$Effect_Size_ID)),
                                        "Estimate" = c(Acclimation_Frequency_Model$b[1]),
                                        "CI Low" = c(Acclimation_Frequency_Model$ci.lb), 
                                        "CI High" = c(Acclimation_Frequency_Model$ci.ub), 
                                        "df" = c(Acclimation_Frequency_Model$ddf[[1]]), 
                                        "p-value" = c(Acclimation_Frequency_Model$pval))

Raw_Acclimation_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Stepwise"),
                                               "Studies" = c(acclimation_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], acclimation_fluctuation_study["Stepwise", "Study"]), 
                                               "Species" = c(acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], acclimation_fluctuation_table["Stepwise", "group_no"]), 
                                               "Effect Sizes" = c(acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "K"], acclimation_fluctuation_table["Stepwise", "K"]),
                                               "Estimate" = c(Acclimation_Fluctuation_Model$b[[2]], Acclimation_Fluctuation_Model$b[[1]]),
                                               "CI Low" = c(Acclimation_Fluctuation_Model$ci.lb[2], Acclimation_Fluctuation_Model$ci.lb[1]), 
                                               "CI High" = c(Acclimation_Fluctuation_Model$ci.ub[2], Acclimation_Fluctuation_Model$ci.ub[1]), 
                                               "df" = c(Acclimation_Fluctuation_Model$ddf[[2]], Acclimation_Fluctuation_Model$ddf[[1]]), 
                                               "p-value" = c(Acclimation_Fluctuation_Model$pval[2], Acclimation_Fluctuation_Model$pval[1]))

Raw_Acclimation_Trait <- data.frame("Phenotypic Trait Categories" = c("Biochemical Assay", "Physiological"),
                                    "Studies" = c(acclimation_trait_study["Biochemical Assay", "Study"], acclimation_trait_study["Physiological", "Study"]), 
                                    "Species" = c(acclimation_trait_table["Biochemical Assay", "group_no"], acclimation_trait_table["Physiological", "group_no"]), 
                                    "Effect Sizes" = c(acclimation_trait_table["Biochemical Assay", "K"], acclimation_trait_table["Physiological", "K"]),
                                    "Estimate" = c(Acclimation_Trait_Model$b[[1]], Acclimation_Trait_Model$b[[2]]),
                                    "CI Low" = c(Acclimation_Trait_Model$ci.lb[1], Acclimation_Trait_Model$ci.lb[2]), 
                                    "CI High" = c(Acclimation_Trait_Model$ci.ub[1], Acclimation_Trait_Model$ci.ub[2]), 
                                    "df" = c(Acclimation_Trait_Model$ddf[[1]], Acclimation_Trait_Model$ddf[[2]]), 
                                    "p-value" = c(Acclimation_Trait_Model$pval[1], Acclimation_Trait_Model$pval[2]))

Raw_Acclimation_Stage <- data.frame("Life-history Stages" = c("Adult", "Juvenile", "Larva"),
                                    "Studies" = c(acclimation_stage_study["Adult", "Study"], acclimation_stage_study["Juvenile", "Study"], acclimation_stage_study["Larva", "Study"]), 
                                    "Species" = c(acclimation_stage_table["Adult", "group_no"], acclimation_stage_table["Juvenile", "group_no"], acclimation_stage_table["Larva", "group_no"]), 
                                    "Effect Sizes" = c(acclimation_stage_table["Adult", "K"], acclimation_stage_table["Juvenile", "K"], acclimation_stage_table["Larva", "K"]),
                                    "Estimate" = c(Acclimation_Stage_Model$b[[1]], Acclimation_Stage_Model$b[[2]], Acclimation_Stage_Model$b[[3]]),
                                    "CI Low" = c(Acclimation_Stage_Model$ci.lb[1], Acclimation_Stage_Model$ci.lb[2], Acclimation_Stage_Model$ci.lb[3]), 
                                    "CI High" = c(Acclimation_Stage_Model$ci.ub[1], Acclimation_Stage_Model$ci.ub[2], Acclimation_Stage_Model$ci.ub[3]), 
                                    "df" = c(Acclimation_Stage_Model$ddf[[1]], Acclimation_Stage_Model$ddf[[2]], Acclimation_Stage_Model$ddf[[3]]), 
                                    "p-value" = c(Acclimation_Stage_Model$pval[1], Acclimation_Stage_Model$pval[2], Acclimation_Stage_Model$pval[3]))

Raw_Acclimation_Class <- data.frame("Taxonomic Class" = c("Insecta"),
                                    "Studies" = c(acclimation_class_study["Insecta", "Study"]), 
                                    "Species" = c(acclimation_class_table["Insecta", "group_no"]), 
                                    "Effect Sizes" = c(acclimation_class_table["Insecta", "K"]),
                                    "Estimate" = c(Acclimation_Class_Model$b[[1]]),
                                    "CI Low" = c(Acclimation_Class_Model$ci.lb[1]), 
                                    "CI High" = c(Acclimation_Class_Model$ci.ub[1]), 
                                    "df" = c(Acclimation_Class_Model$ddf[[1]]), 
                                    "p-value" = c(Acclimation_Class_Model$pval[1]))

Raw_Acclimation_Specific_Trait <- data.frame("Specific Phenotypic Traits" = c("Metabolic Rate"),
                                             "Studies" = c(acclimation_specific_trait_study["Metabolic Rate", "Study"]), 
                                             "Species" = c(acclimation_specific_trait_table["Metabolic Rate", "group_no"]), 
                                             "Effect Sizes" = c(acclimation_specific_trait_table["Metabolic Rate", "K"]),
                                             "Estimate" = c(Acclimation_Specific_Trait_Model$b[[1]]),
                                             "CI Low" = c(Acclimation_Specific_Trait_Model$ci.lb[1]), 
                                             "CI High" = c(Acclimation_Specific_Trait_Model$ci.ub[1]), 
                                             "df" = c(Acclimation_Specific_Trait_Model$ddf[[1]]), 
                                             "p-value" = c(Acclimation_Specific_Trait_Model$pval[1]))

Raw_Developmental <- data.frame("Developmental Exposure" = c("MLMA"),
                                "Studies" = c(intercept_study["Developmental", "Study"]), 
                                "Species" = c(intercept_table["Developmental", "group_no"]), 
                                "Effect Sizes" = c(intercept_table["Developmental", "K"]),
                                "Estimate" = c(Developmental_Model$b[1]),
                                "CI Low" = c(Developmental_Model$ci.lb), 
                                "CI High" = c(Developmental_Model$ci.ub), 
                                "df" = c(Developmental_Model$ddf[[1]]), 
                                "p-value" = c(Developmental_Model$pval))

Raw_Developmental_Amplitude <- data.frame("Developmental Exposure" = c("Fluctuation Amplitude"),
                                          "Studies" = c(length(unique(Developmental_Subset_Data$Study_ID))), 
                                          "Species" = c(length(unique(Developmental_Subset_Data$Scientific_Name))), 
                                          "Effect Sizes" = c(length(Developmental_Subset_Data$Effect_Size_ID)),
                                          "Estimate" = c(Developmental_Amplitude_Model$b[1]),
                                          "CI Low" = c(Developmental_Amplitude_Model$ci.lb), 
                                          "CI High" = c(Developmental_Amplitude_Model$ci.ub), 
                                          "df" = c(Developmental_Amplitude_Model$ddf[[1]]), 
                                          "p-value" = c(Developmental_Amplitude_Model$pval))

Raw_Developmental_Fluctuation_Type <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                        "Stepwise"),
                                                 "Studies" = c(developmental_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], developmental_fluctuation_study["Alternating", "Study"], 
                                                               developmental_fluctuation_study["Stepwise", "Study"]), 
                                                 "Species" = c(developmental_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], developmental_fluctuation_table["Alternating", "group_no"], 
                                                               developmental_fluctuation_table["Stepwise", "group_no"]), 
                                                 "Effect Sizes" = c(developmental_fluctuation_table["Sinusoidal (Sine Curve)", "K"], developmental_fluctuation_table["Alternating", "K"], 
                                                                    developmental_fluctuation_table["Stepwise", "K"]),
                                                 "Estimate" = c(Developmental_Fluctuation_Model$b[[2]], Developmental_Fluctuation_Model$b[[1]],
                                                                Developmental_Fluctuation_Model$b[[3]]),
                                                 "CI Low" = c(Developmental_Fluctuation_Model$ci.lb[2], Developmental_Fluctuation_Model$ci.lb[1], 
                                                              Developmental_Fluctuation_Model$ci.lb[3]), 
                                                 "CI High" = c(Developmental_Fluctuation_Model$ci.ub[2], Developmental_Fluctuation_Model$ci.ub[1], 
                                                               Developmental_Fluctuation_Model$ci.ub[3]), 
                                                 "df" = c(Developmental_Fluctuation_Model$ddf[[2]], Developmental_Fluctuation_Model$ddf[[1]], 
                                                          Developmental_Fluctuation_Model$ddf[[3]]), 
                                                 "p-value" = c(Developmental_Fluctuation_Model$pval[2], Developmental_Fluctuation_Model$pval[1], 
                                                               Developmental_Fluctuation_Model$pval[3]))

Raw_Developmental_Trait <- data.frame("Phenotypic Trait Categories" = c("Life-history Traits", "Morphology"),
                                      "Studies" = c(developmental_trait_study["Life-history Traits", "Study"], developmental_trait_study["Morphology", "Study"]), 
                                      "Species" = c(developmental_trait_table["Life-history Traits", "group_no"], developmental_trait_table["Morphology", "group_no"]), 
                                      "Effect Sizes" = c(developmental_trait_table["Life-history Traits", "K"], developmental_trait_table["Morphology", "K"]),
                                      "Estimate" = c(Developmental_Trait_Model$b[[1]], Developmental_Trait_Model$b[[2]]),
                                      "CI Low" = c(Developmental_Trait_Model$ci.lb[1], Developmental_Trait_Model$ci.lb[2]), 
                                      "CI High" = c(Developmental_Trait_Model$ci.ub[1], Developmental_Trait_Model$ci.ub[2]), 
                                      "df" = c(Developmental_Trait_Model$ddf[[1]], Developmental_Trait_Model$ddf[[2]]), 
                                      "p-value" = c(Developmental_Trait_Model$pval[1], Developmental_Trait_Model$pval[2]))

Raw_Developmental_Exposure <- data.frame("Exposure Time" = c("Embryo", "Juvenile", "Larva"),
                                         "Studies" = c(developmental_exposure_study["Embryo", "Study"], developmental_exposure_study["Juvenile", "Study"], 
                                                       developmental_exposure_study["Larva", "Study"]), 
                                         "Species" = c(developmental_exposure_table["Embryo", "group_no"], developmental_exposure_table["Juvenile", "group_no"], 
                                                       developmental_exposure_table["Larva", "group_no"]), 
                                         "Effect Sizes" = c(developmental_exposure_table["Embryo", "K"], developmental_exposure_table["Juvenile", "K"], 
                                                            developmental_exposure_table["Larva", "K"]),
                                         "Estimate" = c(Developmental_Exposure_Model$b[[1]], Developmental_Exposure_Model$b[[2]],
                                                        Developmental_Exposure_Model$b[[3]]),
                                         "CI Low" = c(Developmental_Exposure_Model$ci.lb[1], Developmental_Exposure_Model$ci.lb[2], 
                                                      Developmental_Exposure_Model$ci.lb[3]), 
                                         "CI High" = c(Developmental_Exposure_Model$ci.ub[1], Developmental_Exposure_Model$ci.ub[2], 
                                                       Developmental_Exposure_Model$ci.ub[3]), 
                                         "df" = c(Developmental_Exposure_Model$ddf[[1]], Developmental_Exposure_Model$ddf[[2]], 
                                                  Developmental_Exposure_Model$ddf[[3]]), 
                                         "p-value" = c(Developmental_Exposure_Model$pval[1], Developmental_Exposure_Model$pval[2], 
                                                       Developmental_Exposure_Model$pval[3]))

Raw_Developmental_Class <- data.frame("Taxonomic Class" = c("Arachnida", "Insecta"),
                                      "Studies" = c(developmental_class_study["Arachnida", "Study"], developmental_class_study["Insecta", "Study"]), 
                                      "Species" = c(developmental_class_table["Arachnida", "group_no"], developmental_class_table["Insecta", "group_no"]), 
                                      "Effect Sizes" = c(developmental_class_table["Arachnida", "K"], developmental_class_table["Insecta", "K"]),
                                      "Estimate" = c(Developmental_Class_Model$b[[1]], Developmental_Class_Model$b[[2]]),
                                      "CI Low" = c(Developmental_Class_Model$ci.lb[1], Developmental_Class_Model$ci.lb[2]), 
                                      "CI High" = c(Developmental_Class_Model$ci.ub[1], Developmental_Class_Model$ci.ub[2]), 
                                      "df" = c(Developmental_Class_Model$ddf[[1]], Developmental_Class_Model$ddf[[2]]), 
                                      "p-value" = c(Developmental_Class_Model$pval[1], Developmental_Class_Model$pval[2]))

Raw_Developmental_Specific_Trait <- data.frame("Specific Phenotypic Traits" = c("Development Time", "Length", "Mass"),
                                               "Studies" = c(developmental_specific_trait_study["Development Time", "Study"], developmental_specific_trait_study["Length", "Study"], 
                                                             developmental_specific_trait_study["Mass", "Study"]), 
                                               "Species" = c(developmental_specific_trait_table["Development Time", "group_no"], developmental_specific_trait_table["Length", "group_no"], 
                                                             developmental_specific_trait_table["Mass", "group_no"]), 
                                               "Effect Sizes" = c(developmental_specific_trait_table["Development Time", "K"], developmental_specific_trait_table["Length", "K"], 
                                                                  developmental_specific_trait_table["Mass", "K"]),
                                               "Estimate" = c(Developmental_Specific_Trait_Model$b[[1]], Developmental_Specific_Trait_Model$b[[2]], Developmental_Specific_Trait_Model$b[[3]]),
                                               "CI Low" = c(Developmental_Specific_Trait_Model$ci.lb[1], Developmental_Specific_Trait_Model$ci.lb[2], Developmental_Specific_Trait_Model$ci.lb[3]), 
                                               "CI High" = c(Developmental_Specific_Trait_Model$ci.ub[1], Developmental_Specific_Trait_Model$ci.ub[2], Developmental_Specific_Trait_Model$ci.ub[3]), 
                                               "df" = c(Developmental_Specific_Trait_Model$ddf[[1]], Developmental_Specific_Trait_Model$ddf[[2]], Developmental_Specific_Trait_Model$ddf[[3]]), 
                                               "p-value" = c(Developmental_Specific_Trait_Model$pval[1], Developmental_Specific_Trait_Model$pval[2], Developmental_Specific_Trait_Model$pval[3]))

write.csv(Raw_Overall, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Overall.csv", row.names = FALSE)
write.csv(Raw_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Trait.csv", row.names = FALSE)
write.csv(Raw_Class, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Class.csv", row.names = FALSE)
write.csv(Raw_Specific_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Specific_Trait.csv", row.names = FALSE)
write.csv(Raw_Individual, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Individual.csv", row.names = FALSE)
write.csv(Raw_Individual_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Individual_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Individual_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Individual_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Individual_Class, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Individual_Class.csv", row.names = FALSE)
write.csv(Raw_Aquatic, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Aquatic.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Aquatic_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Aquatic_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Aquatic_Trait.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Plasticity, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Aquatic_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Terrestrial, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial_Trait.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Plasticity, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial_Plasticity.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Specific_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Terrestrial_Specific_Trait.csv", row.names = FALSE)
write.csv(Raw_Acclimation, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Exposure, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Exposure.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Frequency, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Frequency.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Trait.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Stage, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Stage.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Class, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Class.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Specific_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Acclimation_Specific_Trait.csv", row.names = FALSE)
write.csv(Raw_Developmental, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental.csv", row.names = FALSE)
write.csv(Raw_Developmental_Amplitude, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Amplitude.csv", row.names = FALSE)
write.csv(Raw_Developmental_Fluctuation_Type, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Fluctuation_Type.csv", row.names = FALSE)
write.csv(Raw_Developmental_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Trait.csv", row.names = FALSE)
write.csv(Raw_Developmental_Exposure, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Exposure.csv", row.names = FALSE)
write.csv(Raw_Developmental_Class, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Class.csv", row.names = FALSE)
write.csv(Raw_Developmental_Specific_Trait, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Raw_Developmental_Specific_Trait.csv", row.names = FALSE)

# Heterogeneity Table

Heterogeneity_Overall <- data.frame("Models" = c("Overall", "Fluctuation Amplitude", "Fluctuation Type", 
                                                 "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                    "Shared Animal" = c(Overall_Model_i2[6, 1], Amplitude_Model_i2[6, 1], Fluctuation_Model_i2[6, 1], 
                                                        Trait_Model_i2[6, 1], Specific_Trait_Model_i2[6, 1], Class_Model_i2[6, 1]), 
                                    "Measurement" = c(Overall_Model_i2[7, 1], Amplitude_Model_i2[7, 1], Fluctuation_Model_i2[7, 1], 
                                                      Trait_Model_i2[7, 1], NA, Class_Model_i2[7, 1]),
                                    "Observational" = c(Overall_Model_i2[4, 1], Amplitude_Model_i2[4, 1], Fluctuation_Model_i2[4, 1], 
                                                        Trait_Model_i2[4, 1], Specific_Trait_Model_i2[4, 1], Class_Model_i2[4, 1]), 
                                    "Phylogenetic Relatedness" = c(Overall_Model_i2[2, 1], Amplitude_Model_i2[2, 1], Fluctuation_Model_i2[2, 1], 
                                                                   Trait_Model_i2[2, 1], Specific_Trait_Model_i2[2, 1], Class_Model_i2[2, 1]), 
                                    "Species" = c(Overall_Model_i2[5, 1], Amplitude_Model_i2[5, 1], Fluctuation_Model_i2[5, 1], 
                                                  Trait_Model_i2[5, 1], Specific_Trait_Model_i2[5, 1], Class_Model_i2[5, 1]),
                                    "Study" = c(Overall_Model_i2[3, 1], Amplitude_Model_i2[3, 1], Fluctuation_Model_i2[3, 1], 
                                                Trait_Model_i2[3, 1], Specific_Trait_Model_i2[3, 1], Class_Model_i2[3, 1]), 
                                    "Total" = c(Overall_Model_i2[1, 1], Amplitude_Model_i2[1, 1], Fluctuation_Model_i2[1, 1], 
                                                Trait_Model_i2[1, 1], Specific_Trait_Model_i2[1, 1], Class_Model_i2[1, 1]))

Heterogeneity_Individual <- data.frame("Models" = c("Individual-level Traits", "Fluctuation Amplitude", 
                                                    "Fluctuation Type", "Taxonomic Class"), 
                                       "Shared Animal" = c(Individual_Model_i2[6, 1], Individual_Amplitude_Model_i2[6, 1], 
                                                           Individual_Fluctuation_Model_i2[6, 1], Individual_Class_Model_i2[6, 1]),
                                       "Measurement" = c(Individual_Model_i2[7, 1], Individual_Amplitude_Model_i2[7, 1], 
                                                         Individual_Fluctuation_Model_i2[7, 1], Individual_Class_Model_i2[7, 1]),
                                       "Observational" = c(Individual_Model_i2[4, 1], Individual_Amplitude_Model_i2[4, 1], 
                                                           Individual_Fluctuation_Model_i2[4, 1], Individual_Class_Model_i2[4, 1]), 
                                       "Phylogenetic Relatedness" = c(Individual_Model_i2[2, 1], Individual_Amplitude_Model_i2[2, 1], 
                                                                      Individual_Fluctuation_Model_i2[2, 1], Individual_Class_Model_i2[2, 1]), 
                                       "Species" = c(Individual_Model_i2[5, 1], Individual_Amplitude_Model_i2[5, 1], 
                                                     Individual_Fluctuation_Model_i2[5, 1], Individual_Class_Model_i2[5, 1]),
                                       "Study" = c(Individual_Model_i2[3, 1], Individual_Amplitude_Model_i2[3, 1], 
                                                   Individual_Fluctuation_Model_i2[3, 1], Individual_Class_Model_i2[3, 1]), 
                                       "Total" = c(Individual_Model_i2[1, 1], Individual_Amplitude_Model_i2[1, 1], 
                                                   Individual_Fluctuation_Model_i2[1, 1], Individual_Class_Model_i2[1, 1]))

Heterogeneity_Aquatic <- data.frame("Models" = c("Aquatic Organisms", "Exposure Type", "Fluctuation Amplitude", 
                                                 "Fluctuation Type", "Phenotypic Trait Categories"), 
                                    "Shared Animal" = c(Aquatic_Model_i2[6, 1], Aquatic_Plasticity_Model_i2[6, 1], Aquatic_Amplitude_Model_i2[6, 1], 
                                                        Aquatic_Fluctuation_Model_i2[6, 1], Aquatic_Trait_Model_i2[6, 1]),
                                    "Measurement" = c(Aquatic_Model_i2[7, 1], Aquatic_Plasticity_Model_i2[7, 1], Aquatic_Amplitude_Model_i2[7, 1], 
                                                      Aquatic_Fluctuation_Model_i2[7, 1], Aquatic_Trait_Model_i2[7, 1]),
                                    "Observational" = c(Aquatic_Model_i2[4, 1], Aquatic_Plasticity_Model_i2[4, 1], Aquatic_Amplitude_Model_i2[4, 1], 
                                                        Aquatic_Fluctuation_Model_i2[4, 1], Aquatic_Trait_Model_i2[4, 1]), 
                                    "Phylogenetic Relatedness" = c(Aquatic_Model_i2[2, 1], Aquatic_Plasticity_Model_i2[2, 1], Aquatic_Amplitude_Model_i2[2, 1], 
                                                                   Aquatic_Fluctuation_Model_i2[2, 1], Aquatic_Trait_Model_i2[2, 1]), 
                                    "Species" = c(Aquatic_Model_i2[5, 1], Aquatic_Plasticity_Model_i2[5, 1], Aquatic_Amplitude_Model_i2[5, 1], 
                                                  Aquatic_Fluctuation_Model_i2[5, 1], Aquatic_Trait_Model_i2[5, 1]),
                                    "Study" = c(Aquatic_Model_i2[3, 1], Aquatic_Plasticity_Model_i2[3, 1], Aquatic_Amplitude_Model_i2[3, 1], 
                                                Aquatic_Fluctuation_Model_i2[3, 1], Aquatic_Trait_Model_i2[3, 1]), 
                                    "Total" = c(Aquatic_Model_i2[1, 1], Aquatic_Plasticity_Model_i2[1, 1], Aquatic_Amplitude_Model_i2[1, 1], 
                                                Aquatic_Fluctuation_Model_i2[1, 1], Aquatic_Trait_Model_i2[1, 1]))

Heterogeneity_Terrestrial <- data.frame("Models" = c("Terrestrial Organisms", "Exposure Type", "Fluctuation Amplitude", 
                                                     "Fluctuation Type", "Phenotypic Trait Categories", "Specific Phenotypic Traits"), 
                                        "Shared Animal" = c(Terrestrial_Model_i2[6, 1], Terrestrial_Plasticity_Model_i2[6, 1], Terrestrial_Amplitude_Model_i2[6, 1], 
                                                            Terrestrial_Fluctuation_Model_i2[6, 1], Terrestrial_Trait_Model_i2[6, 1], Terrestrial_Specific_Trait_Model_i2[6, 1]),
                                        "Measurement" = c(Terrestrial_Model_i2[7, 1], Terrestrial_Plasticity_Model_i2[7, 1], Terrestrial_Amplitude_Model_i2[7, 1], 
                                                          Terrestrial_Fluctuation_Model_i2[7, 1], Terrestrial_Trait_Model_i2[7, 1], NA),
                                        "Observational" = c(Terrestrial_Model_i2[4, 1], Terrestrial_Plasticity_Model_i2[4, 1], Terrestrial_Amplitude_Model_i2[4, 1], 
                                                            Terrestrial_Fluctuation_Model_i2[4, 1], Terrestrial_Trait_Model_i2[4, 1], Terrestrial_Specific_Trait_Model_i2[4, 1]), 
                                        "Phylogenetic Relatedness" = c(Terrestrial_Model_i2[2, 1], Terrestrial_Plasticity_Model_i2[2, 1], Terrestrial_Amplitude_Model_i2[2, 1], 
                                                                       Terrestrial_Fluctuation_Model_i2[2, 1], Terrestrial_Trait_Model_i2[2, 1], Terrestrial_Specific_Trait_Model_i2[2, 1]), 
                                        "Species" = c(Terrestrial_Model_i2[5, 1], Terrestrial_Plasticity_Model_i2[5, 1], Terrestrial_Amplitude_Model_i2[5, 1], 
                                                      Terrestrial_Fluctuation_Model_i2[5, 1], Terrestrial_Trait_Model_i2[5, 1], Terrestrial_Specific_Trait_Model_i2[5, 1]),
                                        "Study" = c(Terrestrial_Model_i2[3, 1], Terrestrial_Plasticity_Model_i2[3, 1], Terrestrial_Amplitude_Model_i2[3, 1], 
                                                    Terrestrial_Fluctuation_Model_i2[3, 1], Terrestrial_Trait_Model_i2[3, 1], Terrestrial_Specific_Trait_Model_i2[3, 1]), 
                                        "Total" = c(Terrestrial_Model_i2[1, 1], Terrestrial_Plasticity_Model_i2[1, 1], Terrestrial_Amplitude_Model_i2[1, 1], 
                                                    Terrestrial_Fluctuation_Model_i2[1, 1], Terrestrial_Trait_Model_i2[1, 1], Terrestrial_Specific_Trait_Model_i2[1, 1]))

Heterogeneity_Acclimation <- data.frame("Models" = c("Acclimation", "Exposure Time", "Fluctuation Amplitude", "Fluctuation Type", "Life-history Stage", 
                                                     "Number of Fluctuations", "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                        "Shared Animal" = c(Acclimation_Model_i2[6, 1], Acclimation_Exposure_Model_i2[6, 1], Acclimation_Amplitude_Model_i2[6, 1], Acclimation_Fluctuation_Model_i2[6, 1], Acclimation_Stage_Model_i2[6, 1], 
                                                            Acclimation_Frequency_Model_i2[6, 1], Acclimation_Trait_Model_i2[6, 1], Acclimation_Specific_Trait_Model_i2[6, 1], Acclimation_Class_Model_i2[6, 1]),
                                        "Measurement" = c(Acclimation_Model_i2[7, 1], Acclimation_Exposure_Model_i2[7, 1], Acclimation_Amplitude_Model_i2[7, 1], Acclimation_Fluctuation_Model_i2[7, 1], Acclimation_Stage_Model_i2[7, 1], 
                                                          Acclimation_Frequency_Model_i2[7, 1], Acclimation_Trait_Model_i2[7, 1], NA, Acclimation_Class_Model_i2[7, 1]),
                                        "Observational" = c(Acclimation_Model_i2[4, 1], Acclimation_Exposure_Model_i2[4, 1], Acclimation_Amplitude_Model_i2[4, 1], Acclimation_Fluctuation_Model_i2[4, 1], Acclimation_Stage_Model_i2[4, 1], 
                                                            Acclimation_Frequency_Model_i2[4, 1], Acclimation_Trait_Model_i2[4, 1], Acclimation_Specific_Trait_Model_i2[4, 1], Acclimation_Class_Model_i2[4, 1]), 
                                        "Phylogenetic Relatedness" = c(Acclimation_Model_i2[2, 1], Acclimation_Exposure_Model_i2[2, 1], Acclimation_Amplitude_Model_i2[2, 1], Acclimation_Fluctuation_Model_i2[2, 1], Acclimation_Stage_Model_i2[2, 1], 
                                                                       Acclimation_Frequency_Model_i2[2, 1], Acclimation_Trait_Model_i2[2, 1], Acclimation_Specific_Trait_Model_i2[2, 1], Acclimation_Class_Model_i2[2, 1]), 
                                        "Species" = c(Acclimation_Model_i2[5, 1], Acclimation_Exposure_Model_i2[5, 1], Acclimation_Amplitude_Model_i2[5, 1], Acclimation_Fluctuation_Model_i2[5, 1], Acclimation_Stage_Model_i2[5, 1], 
                                                      Acclimation_Frequency_Model_i2[5, 1], Acclimation_Trait_Model_i2[5, 1], Acclimation_Specific_Trait_Model_i2[5, 1], Acclimation_Class_Model_i2[5, 1]),
                                        "Study" = c(Acclimation_Model_i2[3, 1], Acclimation_Exposure_Model_i2[3, 1], Acclimation_Amplitude_Model_i2[3, 1], Acclimation_Fluctuation_Model_i2[3, 1], Acclimation_Stage_Model_i2[3, 1], 
                                                    Acclimation_Frequency_Model_i2[3, 1], Acclimation_Trait_Model_i2[3, 1], Acclimation_Specific_Trait_Model_i2[3, 1], Acclimation_Class_Model_i2[3, 1]), 
                                        "Total" = c(Acclimation_Model_i2[1, 1], Acclimation_Exposure_Model_i2[1, 1], Acclimation_Amplitude_Model_i2[1, 1], Acclimation_Fluctuation_Model_i2[1, 1], Acclimation_Stage_Model_i2[1, 1], 
                                                    Acclimation_Frequency_Model_i2[1, 1], Acclimation_Trait_Model_i2[1, 1], Acclimation_Specific_Trait_Model_i2[1, 1], Acclimation_Class_Model_i2[1, 1]))

Heterogeneity_Developmental <- data.frame("Models" = c("Developmental", "Exposure Time", "Fluctuation Amplitude", "Fluctuation Type", 
                                                       "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                          "Shared Animal" = c(Developmental_Model_i2[6, 1], Developmental_Exposure_Model_i2[6, 1], Developmental_Amplitude_Model_i2[6, 1], Developmental_Fluctuation_Model_i2[6, 1], 
                                                              Developmental_Trait_Model_i2[6, 1], Developmental_Specific_Trait_Model_i2[6, 1], Developmental_Class_Model_i2[6, 1]),
                                          "Measurement" = c(Developmental_Model_i2[7, 1], Developmental_Exposure_Model_i2[7, 1], Developmental_Amplitude_Model_i2[7, 1], Developmental_Fluctuation_Model_i2[7, 1], 
                                                            Developmental_Trait_Model_i2[7, 1], NA, Developmental_Class_Model_i2[7, 1]),
                                          "Observational" = c(Developmental_Model_i2[4, 1], Developmental_Exposure_Model_i2[4, 1], Developmental_Amplitude_Model_i2[4, 1], Developmental_Fluctuation_Model_i2[4, 1], 
                                                              Developmental_Trait_Model_i2[4, 1], Developmental_Specific_Trait_Model_i2[4, 1], Developmental_Class_Model_i2[4, 1]), 
                                          "Phylogenetic Relatedness" = c(Developmental_Model_i2[2, 1], Developmental_Exposure_Model_i2[2, 1], Developmental_Amplitude_Model_i2[2, 1], Developmental_Fluctuation_Model_i2[2, 1], 
                                                                         Developmental_Trait_Model_i2[2, 1], Developmental_Specific_Trait_Model_i2[2, 1], Developmental_Class_Model_i2[2, 1]), 
                                          "Species" = c(Developmental_Model_i2[5, 1], Developmental_Exposure_Model_i2[5, 1], Developmental_Amplitude_Model_i2[5, 1], Developmental_Fluctuation_Model_i2[5, 1], 
                                                        Developmental_Trait_Model_i2[5, 1], Developmental_Specific_Trait_Model_i2[5, 1], Developmental_Class_Model_i2[5, 1]),
                                          "Study" = c(Developmental_Model_i2[3, 1], Developmental_Exposure_Model_i2[3, 1], Developmental_Amplitude_Model_i2[3, 1], Developmental_Fluctuation_Model_i2[3, 1], 
                                                      Developmental_Trait_Model_i2[3, 1], Developmental_Specific_Trait_Model_i2[3, 1], Developmental_Class_Model_i2[3, 1]), 
                                          "Total" = c(Developmental_Model_i2[1, 1], Developmental_Exposure_Model_i2[1, 1], Developmental_Amplitude_Model_i2[1, 1], Developmental_Fluctuation_Model_i2[1, 1], 
                                                      Developmental_Trait_Model_i2[1, 1], Developmental_Specific_Trait_Model_i2[1, 1], Developmental_Class_Model_i2[1, 1]))

write.csv(Heterogeneity_Overall, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Overall.csv", row.names = FALSE)
write.csv(Heterogeneity_Individual, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Individual.csv", row.names = FALSE)
write.csv(Heterogeneity_Aquatic, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Aquatic.csv", row.names = FALSE)
write.csv(Heterogeneity_Terrestrial, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Terrestrial.csv", row.names = FALSE)
write.csv(Heterogeneity_Acclimation, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Acclimation.csv", row.names = FALSE)
write.csv(Heterogeneity_Developmental, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Complex_Heterogeneity_Developmental.csv", row.names = FALSE)
