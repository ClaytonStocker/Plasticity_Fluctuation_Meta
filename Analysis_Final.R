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
                 robumeta, ggpmisc, ggridges, ggbeeswarm, gridExtra, janitor)
  
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
                                                             legend.title = element_text(size = 16),
                                                             legend.text = element_text(size = 16), 
                                                             legend.position = "top",
                                                             plot.tag = element_text(size = 16, face = "italic")))
        }
        
        density_orchard_overall <- orchard_plot(Overall_Model, group = "Study_ID", mod = "1", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + my_theme() + 
          annotate('text',  x =1+0.1, y = 0.18,
                   label= paste("italic(k)==", dim(data)[1], "~","(", length(unique(data$Study_ID)), ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label= paste(format(round(mean(exp(Overall_Model_Estimates[1, "estimate"])-1)*100, 2), nsmall = 2), "%"),
                   x = 1+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Intrcpt" = "Overall")) 
        
        
        indivdual_orchard_overall <- orchard_plot(Individual_Model, group = "Study_ID", mod = "1", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + my_theme() + 
          annotate('text',  x =1+0.1, y = 0.18,
                   label= paste("italic(k)==", dim(Individual_Subset_Data)[1], "~","(", length(unique(Individual_Subset_Data$Study_ID)), ")"), parse = TRUE, hjust = "right", size = 6) +
          annotate('text', label= paste(format(round(mean(exp(Individual_Model_Estimates[1, "estimate"])-1)*100, 2), nsmall = 2), "%"),
                   x = 1+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Intrcpt" = "")) 
        
        size = 24
        position = "topleft"
        fig2 <- (density_orchard_overall + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))  | indivdual_orchard_overall + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) ) + plot_annotation(tag_levels = "a", tag_suffix = ")")
        
        ggsave(filename = "./output/figs/fig2.png", , width = 11.2, height =  5.8)
        
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
        density_trait_orchard <- orchard_plot(Trait_Model, group = "Study_ID", mod = "Trait_Category", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + 
          my_theme() + 
          annotate('text',  x = c(1,2,3,4)+0.1, y = 0.18, label = 
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
            x = c(1,2,3,4)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80")
        
        
        density_specific_trait_orchard <- orchard_plot(Specific_Trait_Model, group = "Study_ID", mod = "Measurement", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.25, 0.25) + 
          my_theme() + 
          annotate('text',  x = c(1,2,3,4)+0.1, y = 0.22, label= paste("italic(k)==", c(specific_trait_table["Development Time", "K"],
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
                   x = c(1,2,3,4)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80")
        
        size = 24
        position = "topleft"
        fig3 <- (density_trait_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) | density_specific_trait_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))) + plot_annotation(tag_levels = "a", tag_suffix = ")") 
        
        ggsave(filename = "./output/figs/fig3.png", fig3, width = 13.7125, height =  7.4125)
        
        
        
        
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
        
        
        density_habitat_orchard <- orchard_plot(habitat_Model, group = "Study_ID", mod = "Ecosystem", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + 
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
                   x = c(1,2)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80")
        
        
        density_vert_invert_orchard <- orchard_plot(vert_invert_Model, group = "Study_ID", mod = "vert_invert", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + 
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
                   x = c(1,2)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80")
        
        size = 24
        position = "topleft"
        fig4 <- (density_habitat_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic")) | density_vert_invert_orchard + theme(plot.tag.position = position, plot.tag = element_text(size = size, face = "italic"))) + plot_annotation(tag_levels = "a", tag_suffix = ")") 
        
        ggsave(filename = "./output/figs/fig4.png", fig4, width = 11.9125, height =  8.049383)
        
        
        
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
      
    
      density_fluctuation_orchard <- orchard_plot(Fluctuation_Model, group = "Study_ID", mod = "Fluctuation_Category", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + 
        my_theme() + 
        annotate('text',  x = c(1,2,3)+0.1, y = 0.18,
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
                 x = c(1,2,3)+0.1, y = -0.15, size = 6) + geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80")
      
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
      
      density_plasticiyMechanism_orchard <- orchard_plot(PlasticityMechanism_Model, group = "Study_ID", mod = "Plasticity_Mechanism", xlab = TeX(" Effect Size ($PRRD_{S}$)"), angle = 45, k = FALSE, g = FALSE, trunk.size = 2) + ylim(-0.2, 0.2) + 
        my_theme() + 
        annotate('text',  x = c(1,2)+0.1, y = 0.18,
                 label= paste("italic(k)==", c(plasticity_mechanism_dat["Acclimation", "k"],
                                               plasticity_mechanism_dat["Developmental Plasticity", "k"]), "~","(", 
                              c(plasticity_mechanism_dat["Acclimation", "group_no"],
                                plasticity_mechanism_dat["Developmental Plasticity", "group_no"]), ")"), parse = TRUE, hjust = "right", size = 6) +
        annotate('text', 
                 label=c(paste(format(round(mean(exp(plasticity_mechanism_dat["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                         paste(format(round(mean(exp(plasticity_mechanism_dat["Developmental Plasticity", "estimate"])-1)*100, 2), nsmall = 2), "%")), x = c(1,2)+0.1, y = -0.15, size = 6) + 
        geom_hline(yintercept =  c(-0.2, -0.1, 0.1, 0.2), linetype = "dashed", colour = "gray80") + scale_x_discrete(labels = c("Developmental Plasticity" = "Development"))
      
      ggsave(filename = "./output/figs/fig6.png", density_plasticiyMechanism_orchard, width = 8.825, height =  7.200)

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
      
##### Supplementary Material Tables #####
      
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
                                                        terrestrial_trait_study["Morphological", "Study"], terrestrial_trait_study["Physiological", "Study"]), 
                                          "Species" = c(terrestrial_trait_table["Biochemical Assay", "group_no"], terrestrial_trait_table["Life-history Traits", "group_no"],
                                                        terrestrial_trait_table["Morphological", "group_no"], terrestrial_trait_table["Physiological", "group_no"]), 
                                          "Effect Sizes" = c(terrestrial_trait_table["Biochemical Assay", "K"], terrestrial_trait_table["Life-history Traits", "K"],
                                                             terrestrial_trait_table["Morphological", "K"], terrestrial_trait_table["Physiological", "K"]),
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
                                            "Studies" = c(developmental_trait_study["Life-history Traits", "Study"], developmental_trait_study["Morphological", "Study"]), 
                                            "Species" = c(developmental_trait_table["Life-history Traits", "group_no"], developmental_trait_table["Morphological", "group_no"]), 
                                            "Effect Sizes" = c(developmental_trait_table["Life-history Traits", "K"], developmental_trait_table["Morphological", "K"]),
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
      
      write.csv(Raw_Overall, file = "./Complex_Raw_Overall.csv", row.names = FALSE)
      write.csv(Raw_Amplitude, file = "./Complex_Raw_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Fluctuation_Type, file = "./Complex_Raw_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Trait, file = "./Complex_Raw_Trait.csv", row.names = FALSE)
      write.csv(Raw_Class, file = "./Complex_Raw_Class.csv", row.names = FALSE)
      write.csv(Raw_Specific_Trait, file = "./Complex_Raw_Specific_Trait.csv", row.names = FALSE)
      write.csv(Raw_Individual, file = "./Complex_Raw_Individual.csv", row.names = FALSE)
      write.csv(Raw_Individual_Amplitude, file = "./Complex_Raw_Individual_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Individual_Fluctuation_Type, file = "./Complex_Raw_Individual_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Individual_Class, file = "./Complex_Raw_Individual_Class.csv", row.names = FALSE)
      write.csv(Raw_Aquatic, file = "./Complex_Raw_Aquatic.csv", row.names = FALSE)
      write.csv(Raw_Aquatic_Amplitude, file = "./Complex_Raw_Aquatic_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Aquatic_Fluctuation_Type, file = "./Complex_Raw_Aquatic_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Aquatic_Trait, file = "./Complex_Raw_Aquatic_Trait.csv", row.names = FALSE)
      write.csv(Raw_Aquatic_Plasticity, file = "./Complex_Raw_Aquatic_Plasticity.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial, file = "./Complex_Raw_Terrestrial.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial_Amplitude, file = "./Complex_Raw_Terrestrial_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial_Fluctuation_Type, file = "./Complex_Raw_Terrestrial_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial_Trait, file = "./Complex_Raw_Terrestrial_Trait.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial_Plasticity, file = "./Complex_Raw_Terrestrial_Plasticity.csv", row.names = FALSE)
      write.csv(Raw_Terrestrial_Specific_Trait, file = "./Complex_Raw_Terrestrial_Specific_Trait.csv", row.names = FALSE)
      write.csv(Raw_Acclimation, file = "./Complex_Raw_Acclimation.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Amplitude, file = "./Complex_Raw_Acclimation_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Exposure, file = "./Complex_Raw_Acclimation_Exposure.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Frequency, file = "./Complex_Raw_Acclimation_Frequency.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Fluctuation_Type, file = "./Complex_Raw_Acclimation_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Trait, file = "./Complex_Raw_Acclimation_Trait.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Stage, file = "./Complex_Raw_Acclimation_Stage.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Class, file = "./Complex_Raw_Acclimation_Class.csv", row.names = FALSE)
      write.csv(Raw_Acclimation_Specific_Trait, file = "./Complex_Raw_Acclimation_Specific_Trait.csv", row.names = FALSE)
      write.csv(Raw_Developmental, file = "./Complex_Raw_Developmental.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Amplitude, file = "./Complex_Raw_Developmental_Amplitude.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Fluctuation_Type, file = "./Complex_Raw_Developmental_Fluctuation_Type.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Trait, file = "./Complex_Raw_Developmental_Trait.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Exposure, file = "./Complex_Raw_Developmental_Exposure.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Class, file = "./Complex_Raw_Developmental_Class.csv", row.names = FALSE)
      write.csv(Raw_Developmental_Specific_Trait, file = "./Complex_Raw_Developmental_Specific_Trait.csv", row.names = FALSE)
      
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
      
      write.csv(Heterogeneity_Overall, file = "./Complex_Heterogeneity_Overall.csv", row.names = FALSE)
      