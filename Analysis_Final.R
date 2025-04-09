##### Setup #####
# Clean working space and load packages
  rm(list = ls())
  if (!require("pacman")) install.packages("pacman")
  if (!require("devtools")) install.packages("devtools")
  if (!require("metaAidR")) devtools::install_github("daniel1noble/metaAidR", force = TRUE)
  if (!require("orchaRd")) remotes::install_github("daniel1noble/orchaRd", dependencies = TRUE, force = TRUE)
  pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
                 tidyr, ggplot2, rotl, DescTools, stringr, ape, 
                 emmeans, patchwork, latex2exp, metafor, brms, 
                 flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
                 robumeta, ggpmisc, ggridges, ggbeeswarm, gridExtra, janitor, rphylopic)
  
  source("func.R")

# Importing Data Set
                    data <- read.csv("./Complex_Final_Data.csv")
                data$obs <- 1:nrow(data)
    data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
              data$phylo <- data$Scientific_Name
        data$vert_invert <- ifelse(data$Phylum == "Chordata" , "Vertebrate", "Invertebrate")
        
# Calculate Effect sizes 
        
        data <- data  %>% 
          mutate(PRRD = PRRD(t1 = T1_constant, t2 = T2_constant, 
                             t1_c = Mean_Transformed_T1_C, t2_c = Mean_Transformed_T2_C, t1_f = Mean_Transformed_T1_F, t2_f = Mean_Transformed_T2_F, 
                             sd_t1_c = SD_Final_Transformed_T1_C, sd_t2_c= SD_Final_Transformed_T2_C , sd_t1_f = SD_Final_Transformed_T1_F, sd_t2_f = SD_Final_Transformed_T2_F, 
                             n_t1_c =  n_T1_C, n_t2_c = n_T2_C, n_t1_f = n_T1_F, n_t2_f = n_T2_F, type = 'ef'),
                 v_PRRD = PRRD(t1 = T1_constant, t2 = T2_constant, 
                               t1_c = Mean_Transformed_T1_C, t2_c = Mean_Transformed_T2_C, t1_f = Mean_Transformed_T1_F, t2_f = Mean_Transformed_T2_F, 
                               sd_t1_c = SD_Final_Transformed_T1_C, sd_t2_c=  SD_Final_Transformed_T2_C, sd_t1_f = SD_Final_Transformed_T1_F, sd_t2_f = SD_Final_Transformed_T2_F, 
                               n_t1_c =  n_T1_C, n_t2_c = n_T2_C, n_t1_f = n_T1_F, n_t2_f = n_T2_F, type = 'v'))

# Phylogenetic covariance matrix
            tree <- ape::read.tree("./Complex_tree")
             phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
               A <- ape::vcv.phylo(phy)
    row.names(A) <- colnames(A) <- row.names(A)
           A_cor <- ape::vcv.phylo(phy, corr = TRUE)

# Periods used in different studies 
        sum_period <- data %>% 
                      group_by(Fluctuation_Unit)  %>% 
                      summarise(n = length(unique(Study_ID)), per = n/44*100) 

# Variance Matrix (Shared Control)
  VCV <- make_VCV_matrix(data, V = "v_PRRD", cluster = "Shared_Control_Number")
  
  
##### Overall Model #####
  run <- TRUE
  system.time(
    if(run){
      Overall_Model <- metafor::rma.mv(PRRD ~ 1, V = VCV, test = "t", 
                                        random = list(~1|phylo, 
                                                     ~1|Study_ID, 
                                                     ~1|obs, 
                                                     ~1|Scientific_Name, 
                                                     ~1|Shared_Animal_Number, 
                                                     ~1|Measurement), 
                                       R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                       control=list(rel.tol=1e-9))
      saveRDS(Overall_Model, "./output/models/Complex_Overall_Model.rds")
    } else {
      Overall_Model <- readRDS("./output/models/Complex_Overall_Model.rds")
    })
  
  # Check robustness of results to non-independence
  Overall_Model_rob <- robust(Overall_Model, cluster = data$Study_ID, adjust = TRUE)
  
  # Extract Estimates
  Overall_Model_Estimates <- data.frame(estimate = Overall_Model$b, 
                                           ci.lb = Overall_Model$ci.lb, 
                                           ci.ub = Overall_Model$ci.ub)
  # Heterogeneity 
        Overall_Model_i2 <- data.frame(round(orchaRd::i2_ml(Overall_Model), 2))

##### Individual-Level Trait Subset Model #####
        
  # Subset out individual level trait data
        Individual_Subset_Data <- data %>% filter(Trait_Category != "Population")
        Individual_Species <- Individual_Subset_Data %>% select("phylo") %>% unique()
  
  # Prune the phylogeny
        Individual_A_cor <- as.data.frame(A_cor)
        Individual_A_cor <- Individual_A_cor[c(Individual_Species$phylo), c(Individual_Species$phylo)]
        Individual_A_cor <- as.matrix(Individual_A_cor)
  
  # Create VCV      
        Individual_VCV <- make_VCV_matrix(Individual_Subset_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
  
  # Fit model
        run <- TRUE
        system.time(
          if(run){
            Individual_Model <- metafor::rma.mv(PRRD ~ 1, V = Individual_VCV, test = "t", 
                                                random = list(~1|phylo, 
                                                              ~1|Study_ID, 
                                                              ~1|obs, 
                                                              ~1|Scientific_Name, 
                                                              ~1|Shared_Animal_Number, 
                                                              ~1|Measurement), 
                                                R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE,
                                                control=list(rel.tol=1e-9))
            saveRDS(Individual_Model, "./output/models/Complex_Individual_Model.rds")
          } else {
            Individual_Model <- readRDS("./output/models/Complex_Individual_Model.rds")
          })
        
        # Check robustness
        Individual_Model_rob <- robust(Individual_Model, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)
        
        # Extract estimates
        Individual_Model_Estimates <- data.frame(estimate = Individual_Model$b, 
                                                 ci.lb = Individual_Model$ci.lb, 
                                                 ci.ub = Individual_Model$ci.ub)

##### Figure 2 #####
        my_theme <- function() {list( theme_classic() ,theme(axis.text.y = element_text(size = 16), 
                                                             axis.text.x = element_text(margin = margin(b = 5), size = 16), 
                                                             axis.ticks = element_blank(),
                                                             axis.title = element_text(size = 18),
                                                             legend.title = element_text(size = 15),
                                                             legend.text = element_text(size = 15), 
                                                             legend.position = "top",
                                                             plot.tag = element_text(size = 16, face = "italic")))
        }
  
        trunk.size = 1
        density_orchard_overall <- orchard_plot(Overall_Model, group = "Study_ID", mod = "1", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size) + ylim(-0.18, 0.18) + my_theme() + 
          annotate('text',  x =1+0.1, y = 0.18,
                   label= paste("italic(k)==", dim(data)[1], "~","(", length(unique(data$Study_ID)), ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label= paste(format(round(mean(exp(Overall_Model_Estimates[1, "estimate"])-1)*100, 2), nsmall = 2), "%"),
                   x = 1+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Intrcpt" = "Overall")) + scale_fill_manual(values = "gray") +  scale_colour_manual(values = "black")
        
        
        indivdual_orchard_overall <- orchard_plot(Individual_Model, group = "Study_ID", mod = "1", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size) + ylim(-0.18, 0.18) + my_theme() + 
          annotate('text',  x =1+0.1, y = 0.18,
                   label= paste("italic(k)==", dim(Individual_Subset_Data)[1], "~","(", length(unique(Individual_Subset_Data$Study_ID)), ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label= paste(format(round(mean(exp(Individual_Model_Estimates[1, "estimate"])-1)*100, 2), nsmall = 2), "%"),
                   x = 1+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Intrcpt" = "")) + scale_fill_manual(values = "gray") + scale_colour_manual(values = "black")
        
        size = 24
        position = "topleft"
        fig2 <- (density_orchard_overall + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))  | indivdual_orchard_overall + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) ) + plot_annotation(tag_levels = "a", tag_suffix = ")")
        
        ggsave(filename = "./output/figs/fig2.png", , width = 12, height =  5.8)
        
##### Overall Model - Trait Meta-Regression #####
        
        # Have a look at the data
        Trait_Exploration <- data %>% select("Trait_Category") %>% table() %>% data.frame()
        rownames(Trait_Exploration) <- Trait_Exploration$Trait_Category
        
        # Exclude some categories with low numbers of effects
        Trait_Data <- data %>% filter(Trait_Category != "Behavioural" &
                                        Trait_Category != "Gene Expression" &
                                        Trait_Category != "Population")
        
        # How many species?
        Trait_Species_Count <- Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
          filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
        rownames(Trait_Species_Count) <- Trait_Species_Count$Trait_Category
        
        # How many studies
        Trait_Study_Count <- Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
          filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
        rownames(Trait_Study_Count) <- Trait_Study_Count$Trait_Category
        
        # Phylo matrix
        Trait_Species <- Trait_Data %>% select("phylo") %>% unique()
        Trait_A_cor <- as.data.frame(A_cor)
        Trait_A_cor <- Trait_A_cor[c(Trait_Species$phylo), c(Trait_Species$phylo)]
        Trait_A_cor <- as.matrix(Trait_A_cor)
        
        # VCV matrix
        Trait_VCV <- make_VCV_matrix(Trait_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
        
        run <- TRUE
        system.time(
          if(run){
            Trait_Model <- metafor::rma.mv(PRRD, V = Trait_VCV, test = "t", 
                                           mods = ~ Trait_Category - 1,
                                           random = list(~1|phylo, 
                                                         ~1|Study_ID, 
                                                         ~1|obs, 
                                                         ~1|Scientific_Name, 
                                                         ~1|Shared_Animal_Number, 
                                                         ~1|Measurement), 
                                           R = list(phylo=Trait_A_cor), data = Trait_Data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
            saveRDS(Trait_Model, "./output/models/Complex_Trait_Model.rds")
          } else {
            Trait_Model <- readRDS("./output/models/Complex_Trait_Model.rds")
          })
        
        # Check robustness
        Trait_Model_rob <- robust(Trait_Model, cluster = Trait_Data$Study_ID, adjust = TRUE)
        
        # Extract estimates
        Trait_Model_Estimates <- data.frame(Category = substr(row.names(Trait_Model$b), 15, 100),
                                            estimate = Trait_Model$b, 
                                            ci.lb = Trait_Model$ci.lb, 
                                            ci.ub = Trait_Model$ci.ub,
                                            df = Trait_Model$ddf,
                                            pval = Trait_Model$pval)
        rownames(Trait_Model_Estimates) <- Trait_Model_Estimates$Category
        
##### Overall Model - Specific Trait Meta-Regression #####
        
        # Check data
        Specific_Trait_Exploration <- data %>% select("Measurement") %>% table() %>% data.frame()
        Specific_Trait_Exploration <- Specific_Trait_Exploration %>% filter(Freq > 10)
        rownames(Specific_Trait_Exploration) <- Specific_Trait_Exploration$Measurement
        
        # Filter to categories with sufficient data
        Specific_Trait_Data <- data %>% filter(Measurement == "Development Time"| 
                                                 Measurement == "Length"|
                                                 Measurement == "Mass"|
                                                 Measurement == "Metabolic Rate")
        # How many species?
        Specific_Trait_Species_Count <- Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% 
          table() %>% data.frame() %>% 
          filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
        rownames(Specific_Trait_Species_Count) <- Specific_Trait_Species_Count$Measurement
        
        # How many studies
        Specific_Trait_Study_Count <- Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
          filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
        rownames(Specific_Trait_Study_Count) <- Specific_Trait_Study_Count$Measurement
        
        # Prune Phylogeny
        Specific_Trait_Species <- Specific_Trait_Data %>% select("phylo") %>% unique()
        Specific_Trait_A_cor <- as.data.frame(A_cor)
        Specific_Trait_A_cor <- Specific_Trait_A_cor[c(Specific_Trait_Species$phylo), c(Specific_Trait_Species$phylo)]
        Specific_Trait_A_cor <- as.matrix(Specific_Trait_A_cor)
        
        # Build VCV
        Specific_Trait_VCV <- make_VCV_matrix(Specific_Trait_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
        
        run <- TRUE
        system.time(
          if(run){
            Specific_Trait_Model <- metafor::rma.mv(PRRD, V = Specific_Trait_VCV, test = "t", 
                                                    mods = ~ Measurement - 1,
                                                    random = list(~1|phylo, 
                                                                  ~1|Study_ID, 
                                                                  ~1|obs, 
                                                                  ~1|Scientific_Name, 
                                                                  ~1|Shared_Animal_Number), 
                                                    R = list(phylo=Specific_Trait_A_cor), data = Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                    control=list(rel.tol=1e-9))
            saveRDS(Specific_Trait_Model, "./output/models/Complex_Specific_Trait_Model.rds")
          } else {
            Specific_Trait_Model <- readRDS("./output/models/Complex_Specific_Trait_Model.rds")
          })
        
        # Check robustness
        Specific_Trait_Model_rob <- robust(Specific_Trait_Model, cluster = Specific_Trait_Data$Study_ID, adjust = TRUE)
        
        # Extract model estimates
        Specific_Trait_Model_Estimates <- data.frame(Trait = substr(row.names(Specific_Trait_Model$b), 12, 100),
                                                     estimate = Specific_Trait_Model$b, ci.lb = Specific_Trait_Model$ci.lb, 
                                                     ci.ub = Specific_Trait_Model$ci.ub)
        rownames(Specific_Trait_Model_Estimates) <- Specific_Trait_Model_Estimates$Trait
        
##### Figure 3  #####
        # Preparing Graph 
        
        trait_rnames <- c("Biochemical Assay", "Life-history Traits", 
                          "Morphological", "Physiological")
        
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
                                  df = Trait_Model_Estimates[,"df"],
                                  p = Trait_Model_Estimates[,"pval"],
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
                            replicate(54, "Morphological"),
                            replicate(41, "Physiological"))
        
        trait_raw_df <- data.frame("Model" = trait_raw_name, 
                                   "Effect" = trait_raw_mean)
        
        
        # Preparing Graph 
        
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
        
        trunk.size = 1
        branch.size = 1.5
        density_trait_orchard <- orchard_plot(Trait_Model, group = "Study_ID", mod = "Trait_Category", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.18, 0.18) + 
          my_theme() + 
          annotate('text',  x = c(1,2,3,4)+0.25, y = 0.18, label = 
                     paste("italic(k)==", c(trait_table["Biochemical Assay", "K"], 
                                            trait_table["Life-history Traits", "K"],
                                            trait_table["Morphological", "K"],
                                            trait_table["Physiological", "K"]), "~","(", 
                           c(trait_table["Biochemical Assay", "group_no"],
                             trait_table["Life-history Traits", "group_no"],
                             trait_table["Morphological", "group_no"],
                             trait_table["Physiological", "group_no"]), ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label=c(
            paste(format(round(mean(exp(Trait_Model_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
            paste(format(round(mean(exp(Trait_Model_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"),
            paste(format(round(mean(exp(Trait_Model_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"),
            paste(format(round(mean(exp(Trait_Model_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
            x = c(1,2,3,4)+0.25, y = -0.10, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") +  scale_colour_manual(values = c("black", "black", "black", "black"))
        
        
        density_specific_trait_orchard <- orchard_plot(Specific_Trait_Model, group = "Study_ID", mod = "Measurement", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.12, 0.12) + 
          my_theme() + 
          annotate('text',  x = c(1,2,3,4)+0.25, y = 0.12, label= paste("italic(k)==", c(specific_trait_table["Development Time", "K"],
                                                                                        specific_trait_table["Length", "K"],
                                                                                        specific_trait_table["Mass", "K"],
                                                                                        specific_trait_table["Metabolic Rate", "K"]), "~","(", 
                                                                       c(specific_trait_table["Development Time", "group_no"],
                                                                         specific_trait_table["Length", "group_no"],
                                                                         specific_trait_table["Mass", "group_no"],
                                                                         specific_trait_table["Metabolic Rate", "group_no"]), 
                                                                       ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label=c(paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%"), paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                   paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                   paste(format(round(mean(exp(Specific_Trait_Model_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                   x = c(1,2,3,4)+0.25, y = -0.08, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + annotate('text',  x = c(1,3)+0.25, y = -0.025, label = "*", size = 10) +  scale_colour_manual(values = c("black", "black", "black", "black"))
        
        size = 24
        position = "topleft"
        fig3 <- (density_trait_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) | density_specific_trait_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))) + plot_annotation(tag_levels = "a", tag_suffix = ")") 
        
        ggsave(filename = "./output/figs/fig3.png", fig3, width = 15, height =  7.4125)
        
        
##### Overall Model - Invertebrate/Vertebrate Meta-Regression #####
        
        # Lets have a look at data in each category
        vert_invert_Exploration <- data %>% select("vert_invert") %>% table() %>% data.frame()
        rownames(vert_invert_Exploration) <- vert_invert_Exploration$vert_invert
        
        # Run model
        run <- TRUE
        system.time(
          if(run){
            vert_invert_Model <- metafor::rma.mv(PRRD ~ vert_invert-1, V = VCV, test = "t", 
                                                 random = list(~1|phylo, 
                                                               ~1|Study_ID, 
                                                               ~1|obs, 
                                                               ~1|Scientific_Name, 
                                                               ~1|Shared_Animal_Number, 
                                                               ~1|Measurement), 
                                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
            saveRDS(vert_invert_Model, "./output/models/Complex_vert_invert_Model.rds")
          } else {
            vert_invert_Model <- readRDS("./output/models/Complex_vert_invert_Model.rds")
          })
        
        # Check robustness
        vert_invert_Model_rob <- robust(vert_invert_Model, cluster = data$Study_ID, adjust = TRUE)
        
        # Extract estimates
        vert_invert_Model_Estimates <- data.frame(vert_invert = substr(row.names(vert_invert_Model_rob$b), 6, 100),
                                                  estimate = vert_invert_Model_rob$b, 
                                                  ci.lb = vert_invert_Model_rob$ci.lb, 
                                                  ci.ub = vert_invert_Model_rob$ci.ub,
                                                  df = vert_invert_Model$ddf,
                                                  pval = vert_invert_Model$pval)
        rownames(vert_invert_Model_Estimates) <- NULL
        
##### Overall Model - Habitat Meta-Regression #####
        
        run <- TRUE
        system.time(
          if(run){
            habitat_Model <- metafor::rma.mv(PRRD ~ Ecosystem-1, V = VCV, test = "t", 
                                             random = list(~1|phylo, 
                                                           ~1|Study_ID, 
                                                           ~1|obs, 
                                                           ~1|Scientific_Name, 
                                                           ~1|Shared_Animal_Number, 
                                                           ~1|Measurement), 
                                             R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
            saveRDS(habitat_Model, "./output/models/Complex_habitat_Model.rds")
          } else {
            habitat_Model <- readRDS("./output/models/Complex_habitat_Model.rds")
          })
        
        # Check robustness
        habitat_Model_rob <- robust(habitat_Model, cluster = data$Study_ID, adjust = TRUE)
        
        # Extract estimates
        habitat_Model_Estimates <- data.frame(habitat = substr(row.names(habitat_Model_rob$b), 6, 100),
                                              estimate = habitat_Model_rob$b, 
                                              ci.lb = habitat_Model_rob$ci.lb, 
                                              ci.ub = habitat_Model_rob$ci.ub,
                                              df = habitat_Model$ddf,
                                              pval = habitat_Model$pval)
        rownames(habitat_Model_Estimates) <- NULL
        
##### Figure 4 #####
        
        invert_vert_table <- data  %>% group_by(vert_invert) %>% summarise(group_no = n_distinct(Study_ID), spp = n_distinct(phylo), k = n()) %>% cbind(vert_invert_Model_Estimates[,-1])
        
        habitat_table <- data  %>% group_by(Ecosystem) %>% summarise(group_no = n_distinct(Study_ID), spp = n_distinct(phylo), k = n()) %>% cbind(habitat_Model_Estimates[,-1])
        
        trunk.size = 1
        branch.size = 1.5
        col <- c("#0871B9", "#1BB908")
        
        # Get a single image UUID for a species
        uuid_terrs <- get_uuid(name = "Anolis_sagrei")
        # Get the image for that UUID
        terrest <- get_phylopic(uuid = uuid_terrs)
        
        # Get a single image UUID for a species
        uuid_aq <- get_uuid(name = "Oryzias_latipes")
        # Get the image for that UUID
        aq <- get_phylopic(uuid = uuid_aq)
        
        
        density_habitat_orchard <- orchard_plot(habitat_Model, group = "Study_ID", mod = "Ecosystem", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.2, 0.2) + 
          my_theme() + 
          annotate('text',  x = c(1,2)+0.1, y = 0.18, 
                   label= paste("italic(k)==", 
                                c(habitat_table[1, "k"], 
                                  habitat_table[2, "k"]), "~","(", 
                                c(habitat_table[1, "group_no"], 
                                  habitat_table[2, "group_no"]),
                                ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label=c(paste(format(round(mean(exp(habitat_table[1, "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                   paste(format(round(mean(exp(habitat_table[2, "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                   x = c(1,2)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_fill_manual(values = col)  +  scale_colour_manual(values = c("black", "black")) + add_phylopic(img = terrest, x = 2.4, y = -0.15, height = 0.35) + add_phylopic(img = aq, x = 1.4, y = -0.15, height = 0.15)
        
        
        # Get a single image UUID for a species
        uuid_inv <- get_uuid(name = "Aedes_aegypti")
        # Get the image for that UUID
        invert <- get_phylopic(uuid = uuid_inv)
        
        # Get a single image UUID for a species
        uuid_vert <- get_uuid(name = "Chrysemys_picta")
        # Get the image for that UUID
        vert <- get_phylopic(uuid = uuid_vert)
        vert <- rotate_phylopic(img = vert, angle = 45)
        
        
        density_vert_invert_orchard <- orchard_plot(vert_invert_Model, group = "Study_ID", mod = "vert_invert", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.2, 0.2) + 
          my_theme() + 
          annotate('text',  x = c(1,2)+0.1, y = 0.18, 
                   label= paste("italic(k)==", 
                                c(invert_vert_table[1, "k"], 
                                  invert_vert_table[2, "k"]), "~","(", 
                                c(invert_vert_table[1, "group_no"], 
                                  invert_vert_table[2, "group_no"]),
                                ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label=c(paste(format(round(mean(exp(vert_invert_Model_Estimates[1, "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                   paste(format(round(mean(exp(vert_invert_Model_Estimates[2, "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                   x = c(1,2)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + add_phylopic(img = vert, x = 2.4, y = -0.15, height = 0.45) + add_phylopic(img = invert, x = 1.4, y = -0.15, height = 0.25) +  scale_colour_manual(values = c("black", "black"))
        
        size = 24
        position = "topleft"
        fig4 <- (density_habitat_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) | density_vert_invert_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))) + plot_annotation(tag_levels = "a", tag_suffix = ")") 
        
        ggsave(filename = "./output/figs/fig4.png", fig4, width =14, height =  6)
        
        
        
##### Overall Model - Fluctuation Amplitude Meta-Regression ####
        run <- TRUE
        system.time(
          if(run){
            Amplitude_Model <- metafor::rma.mv(PRRD, V = VCV, test = "t", 
                                               mods = ~ Fluctuation_Magnitude,
                                               random = list(~1|phylo, 
                                                             ~1|Study_ID, 
                                                             ~1|obs, 
                                                             ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, 
                                                             ~1|Measurement), 
                                               R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
            saveRDS(Amplitude_Model, "./output/models/Complex_Amplitude_Model.rds")
          } else {
            Amplitude_Model <- readRDS("./output/models/Complex_Amplitude_Model.rds")
          })
        
    # Check robustness of results to non-independence
        Amplitude_Model_rob <- robust(Amplitude_Model, cluster = data$Study_ID, adjust = TRUE)
    
    # Extract estimates   
        Amplitude_Model_Estimates <- data.frame(estimate = Amplitude_Model$b, 
                                                   ci.lb = Amplitude_Model$ci.lb, 
                                                   ci.ub = Amplitude_Model$ci.ub)
    
##### Overall Model - Type of Fluctuation Meta-Regression ####     
       # Filter missing data
          Fluctuation_Data <- data %>% filter(!is.na(Fluctuation_Category))
        
       # Have a look at breakdown
          Fluctuation_Exploration <- Fluctuation_Data %>% 
                                      select("Fluctuation_Category") %>% 
                                      table() %>% data.frame()
          rownames(Fluctuation_Exploration) <- Fluctuation_Exploration$Fluctuation_Category
          
        # Exclude Stochastic as only k = 3
          Fluctuation_Data <- Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")
        
        # How many species?
        Fluctuation_Species_Count <- Fluctuation_Data %>% 
                                      select("Scientific_Name", "Fluctuation_Category") %>% 
                                        table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                                        select("Fluctuation_Category") %>% table() %>% data.frame()
        rownames(Fluctuation_Species_Count) <- Fluctuation_Species_Count$Fluctuation_Category
        
        # How many studies?
          Fluctuation_Study_Count <- Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% 
            table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
            select("Fluctuation_Category") %>% table() %>% data.frame()
          rownames(Fluctuation_Study_Count) <- Fluctuation_Study_Count$Fluctuation_Category
        
        # How many species
        Fluctuation_Species <- Fluctuation_Data %>% select("phylo") %>% unique()
        Fluctuation_A_cor <- as.data.frame(A_cor)
        Fluctuation_A_cor <- Fluctuation_A_cor[c(Fluctuation_Species$phylo), c(Fluctuation_Species$phylo)]
        Fluctuation_A_cor <- as.matrix(Fluctuation_A_cor)
        
        # Create VCV for the fluctuation category
        Fluctuation_VCV <- make_VCV_matrix(Fluctuation_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
        
        run <- TRUE
        system.time(
          if(run){
            Fluctuation_Model <- metafor::rma.mv(PRRD, V = Fluctuation_VCV, test = "t", 
                                                 mods = ~ Fluctuation_Category-1,
                                                 random = list(~1|phylo, 
                                                               ~1|Study_ID, 
                                                               ~1|obs, 
                                                               ~1|Scientific_Name, 
                                                               ~1|Shared_Animal_Number, 
                                                               ~1|Measurement), 
                                                 R = list(phylo=Fluctuation_A_cor), data = Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
            saveRDS(Fluctuation_Model, "./output/models/Complex_Fluctuation_Model.rds")
          } else {
            Fluctuation_Model <- readRDS("./output/models/Complex_Fluctuation_Model.rds")
          })
        
        # Check robustness of results
                Fluctuation_Model_rob <- robust(Fluctuation_Model, cluster = Fluctuation_Data$Study_ID, adjust = TRUE)
        
        # Extract estimates
                Fluctuation_Model_Estimates <- data.frame(Category = substr(row.names(Fluctuation_Model$b), 21, 100),
                                                  estimate = Fluctuation_Model$b, ci.lb = Fluctuation_Model$ci.lb, 
                                                  ci.ub = Fluctuation_Model$ci.ub)
      rownames(Fluctuation_Model_Estimates) <- Fluctuation_Model_Estimates$Category   
        
##### Figure 5 ####
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
      
      trunk.size = 1
      branch.size = 1.5
      density_fluctuation_orchard <- orchard_plot(Fluctuation_Model, group = "Study_ID", mod = "Fluctuation_Category", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.18, 0.20)  + 
        my_theme() + 
        annotate('text',  x = c(1,2,3)+0.25, y = 0.2,
                 label= paste("italic(k)==", c(fluctuation_table["Alternating", "K"],
                                               fluctuation_table["Sinusoidal (Sine Curve)", "K"], 
                                               fluctuation_table["Stepwise", "K"]), "~","(", 
                              c(fluctuation_table["Alternating", "group_no"],
                                fluctuation_table["Sinusoidal (Sine Curve)","group_no"], 
                                fluctuation_table["Stepwise", "group_no"]), 
                              ")"), parse = TRUE, hjust = "right", size = 6) +
        annotate('text', label=c(paste(format(round(mean(exp(Fluctuation_Model_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                 paste(format(round(mean(exp(Fluctuation_Model_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                 paste(format(round(mean(exp(Fluctuation_Model_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                 x = c(1,2,3)+0.25, y = -0.1, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") +  scale_colour_manual(values = c("black", "black", "black"))
      
      ggsave(filename = "./output/figs/fig5.png", density_fluctuation_orchard, width = 8.185185, height =  6.975309)
      
##--------------------------------------------##      

##### Overall model - Plasticity_Mechanism Meta-Regression #####
      
      # Filter out population as this doesn't jive with inidviudal-level plasticity
      plasticity_mec_data  <- data %>%  filter(Trait_Category != "Population")
      
      # Prune phylogeny
      Individual_Species <- Individual_Subset_Data %>% select("phylo") %>% unique()
      Individual_A_cor <- as.data.frame(A_cor)
      Individual_A_cor <- Individual_A_cor[c(Individual_Species$phylo), c(Individual_Species$phylo)]
      Individual_A_cor <- as.matrix(Individual_A_cor)
      
      # Create VCV Matrix
      Individual_VCV <- make_VCV_matrix(Individual_Subset_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
      
      # Run model
      run <- TRUE
      system.time(
        if(run){
          PlasticityMechanism_Model <- metafor::rma.mv(PRRD, V = Individual_VCV, test = "t", 
                                                       mods = ~ Plasticity_Mechanism - 1,
                                                       random = list(~1|phylo, 
                                                                     ~1|Study_ID, 
                                                                     ~1|obs, 
                                                                     ~1|Scientific_Name, 
                                                                     ~1|Shared_Animal_Number), 
                                                       R = list(phylo=Individual_A_cor), data = plasticity_mec_data, method = "REML", sparse = TRUE, 
                                                       control=list(rel.tol=1e-9))
          saveRDS(PlasticityMechanism_Model, "./output/models/Complex_PlasticityMechanism_Model.rds")
        } else {
          PlasticityMechanism_Model <- readRDS("./output/models/Complex_PlasticityMechanism_Model.rds")
        })
      
      # Check robustness
      PlasticityMechanism_Model_rob <- robust(PlasticityMechanism_Model, cluster = data$Study_ID, adjust = TRUE)
      
      # Extract effects
      PlasticityMechanism_Model_Estimates <- data.frame(estimate = PlasticityMechanism_Model$b, 
                                                        ci.lb = PlasticityMechanism_Model$ci.lb, 
                                                        ci.ub = PlasticityMechanism_Model$ci.ub,
                                                        df = PlasticityMechanism_Model$ddf,
                                                        pval = PlasticityMechanism_Model$pval)
      
      # Summarise data
      plasticity_mechanism_dat <- plasticity_mec_data %>% group_by(Plasticity_Mechanism) %>% 
                                  summarise(group_no = length(unique(Study_ID)), spp = length(unique(phylo)), k = n())  %>%
                                  cbind(PlasticityMechanism_Model_Estimates) 
      rownames(plasticity_mechanism_dat) <- plasticity_mechanism_dat$Plasticity_Mechanism
      
##### Figure 6 #####
      
      plasticity_mechanism_dat <- plasticity_mec_data %>% group_by(Plasticity_Mechanism) %>% summarise(group_no = length(unique(Study_ID)), spp = length(unique(phylo)), k = n())  %>% cbind(PlasticityMechanism_Model_Estimates) 
      rownames(plasticity_mechanism_dat) <- plasticity_mechanism_dat$Plasticity_Mechanism
      
      trunk.size = 1
      branch.size = 1.5
      density_plasticiyMechanism_orchard <- orchard_plot(PlasticityMechanism_Model, group = "Study_ID", mod = "Plasticity_Mechanism", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = trunk.size, branch.size = branch.size) + ylim(-0.18, 0.20) + 
        my_theme() + 
        annotate('text',  x = c(1,2)+0.25, y = 0.20,
                 label= paste("italic(k)==", c(plasticity_mechanism_dat["Acclimation", "k"],
                                               plasticity_mechanism_dat["Developmental Plasticity", "k"]), "~","(", 
                              c(plasticity_mechanism_dat["Acclimation", "group_no"],
                                plasticity_mechanism_dat["Developmental Plasticity", "group_no"]), ")"), parse = TRUE, hjust = "right", size = 6) +
        annotate('text', 
                 label=c(paste(format(round(mean(exp(plasticity_mechanism_dat["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                         paste(format(round(mean(exp(plasticity_mechanism_dat["Developmental Plasticity", "estimate"])-1)*100, 2), nsmall = 2), "%")), x = c(1,2)+0.25, y = -0.1, size = 6) + 
        geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Developmental Plasticity" = "Development")) + scale_colour_manual(values = c("black", "black", "black")) + scale_fill_manual(values=c("#453781FF", "#287D8EFF"))
      
      ggsave(filename = "./output/figs/fig6.png", density_plasticiyMechanism_orchard, width = 7, height =  5)

##### Individual-Level Subset Model - Fluctuation Amplitude Meta-Regression ####
      run <- TRUE
      system.time(
        if(run){
          Individual_Amplitude_Model <- metafor::rma.mv(PRRD, V = Individual_VCV, test = "t", 
                                                        mods = ~ Fluctuation_Magnitude,
                                                        random = list(~1|phylo, 
                                                                      ~1|Study_ID, 
                                                                      ~1|obs, 
                                                                      ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number, 
                                                                      ~1|Measurement), 
                                                        R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, 
                                                        method = "REML", sparse = TRUE, control=list(rel.tol=1e-9))
          saveRDS(Individual_Amplitude_Model, "./output/models/Complex_Individual_Amplitude_Model.rds")
        } else {
          Individual_Amplitude_Model <- readRDS("./output/models/Complex_Individual_Amplitude_Model.rds")
        })
      
      # Check robustness
      Individual_Amplitude_Model_rob <- robust(Individual_Amplitude_Model, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)
      # Extract model estimates
      Individual_Amplitude_Model_Estimates <- data.frame(estimate = Individual_Amplitude_Model$b, 
                                                         ci.lb = Individual_Amplitude_Model$ci.lb, 
                                                         ci.ub = Individual_Amplitude_Model$ci.ub)
   
##### Individual-Level Subset Model - Type of Fluctuation Meta-Regression ####
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
      
      Individual_Fluctuation_VCV <- make_VCV_matrix(Individual_Fluctuation_Data, V = "v_PRRD", cluster = "Shared_Control_Number")
      
      run <- TRUE
      system.time(
        if(run){
          Individual_Fluctuation_Model <- metafor::rma.mv(PRRD, V = Individual_Fluctuation_VCV, test = "t", dfs = "contain",
                                                          mods = ~ Fluctuation_Category - 1,
                                                          random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                        ~1|Shared_Animal_Number, ~1|Measurement), 
                                                          R = list(phylo=Individual_Fluctuation_A_cor), data = Individual_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                          control=list(rel.tol=1e-9))
          saveRDS(Individual_Fluctuation_Model, "./output/models/Complex_Individual_Fluctuation_Model.rds")
        } else {
          Individual_Fluctuation_Model <- readRDS("./output/models/Complex_Individual_Fluctuation_Model.rds")})
      
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
      
##### Figure 7 ####
      # Plot the fluctuation relationship for overall data set
      Plot_Data <- data
      Plot_Data <- Plot_Data %>% mutate(n_category = ifelse(n_T1_C <= 10, "10", 
                                                            ifelse(n_T1_C > 10 & n_T1_C <= 20, "20", 
                                                                   ifelse(n_T1_C > 20 & n_T1_C <= 30, "30", "> 30"))))
      
      Amplitude_Plot <- ggplot(Plot_Data, aes(x = Fluctuation_Magnitude, y = PRRD)) + 
        geom_point(aes(x = Fluctuation_Magnitude, y = PRRD, 
                       size = fct_relevel(n_category, c("10", "20", "30", "> 30"))), 
                   shape = 21, fill = "#4292c6", alpha = 0.5, show.legend = FALSE) + 
        labs(x = "Fluctuation Amplitude (\u00B0C)", y = expression("Effect Size (PRRD"["S"]*")"), 
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
        labs(x = "Fluctuation Amplitude (\u00B0C)", y = expression("Effect Size (PRRD"["S"]*")"), 
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

      size = 16
      position = "topleft"
      t <-  function() {theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))}
      
      fig7 <- (Amplitude_Plot + t() | Individual_Amplitude_Plot + t()) + plot_annotation(tag_levels = "a", tag_suffix = ")")
      ggsave(fig7, filename= "./output/figs/fig7.png", height = 5, width = 10)
      
##### Individual-Level Subset Model - Invertebrate/Vertebrate Meta-Regression
      # Fit model
      run <- TRUE
      system.time(
        if(run){
          Individual_Taxa <- metafor::rma.mv(PRRD ~ vert_invert-1, V = Individual_VCV, test = "t", 
                                              random = list(~1|phylo, 
                                                            ~1|Study_ID, 
                                                            ~1|obs, 
                                                            ~1|Scientific_Name, 
                                                            ~1|Shared_Animal_Number, 
                                                            ~1|Measurement), 
                                              R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE,
                                              control=list(rel.tol=1e-9))
          saveRDS(Individual_Taxa, "./output/models/Complex_Individual_Taxa.rds")
        } else {
          Individual_Taxa <- readRDS("./output/models/Complex_Individual_Taxa.rds")
        })
      
      # Check robustness
      Individual_Taxa_rob <- robust(Individual_Taxa, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)
      
      # Extract estimates
      Individual_Taxa_Estimates <- data.frame(estimate = Individual_Taxa$b, 
                                               ci.lb = Individual_Taxa$ci.lb, 
                                               ci.ub = Individual_Taxa$ci.ub)
      
##### Supplementary Material Results #####
      
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
      
    # Model Results Tables
      
                     Raw_Overall <- table_results(Overall_Model, study_name = "Study_ID", species_name = "Scientific_Name")
                       Raw_Trait <- table_results(Trait_Model, group = "Trait_Category", study_name = "Study_ID", species_name = "Scientific_Name")
              Raw_Specific_Trait <- table_results(Specific_Trait_Model, group = "Measurement", study_name = "Study_ID", species_name = "Scientific_Name")
                 Raw_Vert_Invert <- table_results(vert_invert_Model, group = "vert_invert", study_name = "Study_ID", species_name = "Scientific_Name")
                     Raw_Habitat <- table_results(habitat_Model, group = "Ecosystem", study_name = "Study_ID", species_name = "Scientific_Name")
                   Raw_Amplitude <- table_results(Amplitude_Model,  study_name = "Study_ID", species_name = "Scientific_Name") # Note intercept and slope (row 2)
            Raw_Fluctuation_Type <- table_results(Fluctuation_Model, group = "Fluctuation_Category", study_name = "Study_ID", species_name = "Scientific_Name")
                  Raw_Individual <- table_results(Individual_Model, study_name = "Study_ID", species_name = "Scientific_Name")
        Raw_Individual_Amplitude <- table_results(Individual_Amplitude_Model,  study_name = "Study_ID", species_name = "Scientific_Name") # Note intercept and slope (row 2)
 Raw_Individual_Fluctuation_Type <- table_results(Individual_Fluctuation_Model, group = "Fluctuation_Category", study_name = "Study_ID", species_name = "Scientific_Name")
 Raw_Individual_Taxa <- table_results(Individual_Taxa, group = "vert_invert", study_name = "Study_ID", species_name = "Scientific_Name")
 
 # Publication bias
  data$Precision <- 1/sqrt(data$v_PRRD)
 
  if(rerun){
    Overall_PubBias <- metafor::rma.mv(PRRD ~ 1 + Year_Z + Precision, V = VCV, test = "t", 
                                     random = list(~1|phylo, 
                                                   ~1|Study_ID, 
                                                   ~1|obs, 
                                                   ~1|Scientific_Name, 
                                                   ~1|Shared_Animal_Number, 
                                                   ~1|Measurement), 
                                     R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                     control=list(rel.tol=1e-9))
      saveRDS(Overall_PubBias, "./output/models/Overall_PubBias")
  } else {
    Overall_PubBias <- readRDS("./output/models/Overall_PubBias")
  }
  
  bubble_plot(Overall_PubBias, mod = "Year_Z", group = "Study_ID", ylab = "PRRDs", xlab = "Publication Year (scaled)")
  
 # Calculate influence diagnostics. Takes a long time so avoid re-running.
 rerun = FALSE
 if(rerun){
   inf <- cooks.distance(Overall_Model)
   saveRDS(inf,"./output/models/inf.rds")
 } else {
   inf <- readRDS("./output/models/inf.rds")
 }
 
 # Check if there are influential (CD => 1). Appears to be one effect that is problematic. 
    plot(inf, type = "o", ylab = "Cook's Distance", xlab = "Effect")
    which(inf >=0.8) # row 119 appears potentially problematic. 
    
 # Refit model without this outlier, do we get different results? Answer is generally, no
    run <- TRUE
    system.time(
      if(run){
        data2 <- data[-119,]
        VCV2 <- VCV[-119, -119]
        
        Overall_Model_outlier <- metafor::rma.mv(PRRD ~ 1, V = VCV2, test = "t", 
                                         random = list(~1|phylo, 
                                                       ~1|Study_ID, 
                                                       ~1|obs, 
                                                       ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, 
                                                       ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data2, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
        saveRDS(Overall_Model_outlier, "./output/models/Complex_Overall_Model_outlier.rds")
      } else {
        Overall_Model_outlier <- readRDS("./output/models/Complex_Overall_Model_outlier.rds")
      })
    
 
 # Write tables for supp
  write.csv(Raw_Overall, file = "./output/tables/Raw_Overall.csv")
  write.csv(Raw_Trait, file = "./output/tables/Raw_Trait.csv")
  write.csv(Raw_Specific_Trait, file = "./output/tables/Raw_Specific_Trait.csv")
  write.csv(Raw_Vert_Invert, file = "./output/tables/Raw_Vert_Invert.csv")
  write.csv(Raw_Habitat, file = "./output/tables/Raw_Habitat.csv")
  write.csv(Raw_Amplitude, file = "./output/tables/Raw_Amplitude.csv")
  write.csv(Raw_Fluctuation_Type, file = "./output/tables/Raw_Fluctuation_Type.csv")
  write.csv(Raw_Individual, file = "./output/tables/Raw_Individual.csv")
  write.csv(Raw_Individual_Fluctuation_Type, file = "./output/tables/Raw_Individual_Fluctuation_Type.csv")
  write.csv(Raw_Individual_Amplitude, file = "./output/tables/Raw_Individual_Amplitude.csv")
  write.csv(Raw_Individual_Taxa, file = "./output/tables/Raw_Individual_Taxa.csv")
  
 # Heterogeneity Table
      write.csv(t(Overall_Model_i2), file = "./output/tables/Complex_Heterogeneity_Overall.csv", row.names = FALSE)
      