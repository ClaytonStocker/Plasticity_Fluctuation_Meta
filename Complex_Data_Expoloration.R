rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, tidyr, ggplot2, rotl, DescTools, stringr, ape)

# Importing Final Data Set
Complex_Data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Complex_Final_Data.csv")

# Testing for Negative Results (Means and variances)

Negative_Test <- Complex_Data %>% filter(Mean_Transformed_T1_C < 0 |
                                         Mean_Transformed_T1_F < 0 |
                                         Mean_Transformed_T2_C < 0 |
                                         Mean_Transformed_T2_F < 0 |
                                         SD_Final_Transformed_T1_C < 0 |
                                         SD_Final_Transformed_T1_F < 0 |
                                         SD_Final_Transformed_T2_C < 0 |
                                         SD_Final_Transformed_T2_F < 0) %>% nrow() 

NA_Zero_Test <- Complex_Data %>% filter(InRR_Transformed == 0|
                                        is.na(InRR_Transformed))

Identical_SD_Test <- Complex_Data %>% filter(SD_Final_Transformed_T1_C == SD_Final_Transformed_T1_F |
                                             SD_Final_Transformed_T2_C == SD_Final_Transformed_T2_F)

# Mean and Variance Relationship (T1_C)
ggplot(Complex_Data, aes(x = log(c(Mean_Transformed_T1_C)), y = log(SD_Final_Transformed_T1_C))) + facet_wrap(vars(Trait_Category)) +  
  geom_point() + labs(y = "log(SD)", x = "log(mean)") + theme_bw()# + ylim(min = -10, max = 10) 

Graph_Points <- Complex_Data %>% select("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T1", 
                                        "Trait_ID", "Trait_Category", "Mean_Transformed_T1_C", 
                                        "SD_Final_Transformed_T1_C") %>% 
                                 filter(`Trait_Category` == "Physiological")
Graph_Points$logmean <- log(Graph_Points$Mean_Transformed_T1_C)
Graph_Points$logsd <- log(Graph_Points$SD_Final_Transformed_T1_C)

# Mean and variance Relationship (T1_F)
ggplot(Complex_Data, aes(x = log(c(Mean_Transformed_T1_F)), y = log(SD_Final_Transformed_T1_F))) + facet_wrap(vars(Trait_Category)) +  
  geom_point() + labs(y = "log(SD)", x = "log(mean)") + theme_bw()# + ylim(min = -10, max = 10) 

Graph_Points2 <- Complex_Data %>% select("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T1", 
                                          "Trait_ID", "Trait_Category", "Mean_Transformed_T1_F", 
                                          "SD_Final_Transformed_T1_F") %>% 
                                  filter(`Trait_Category` == "Physiological")
Graph_Points2$logmean <- log(Graph_Points2$Mean_Transformed_T1_F)
Graph_Points2$logsd <- log(Graph_Points2$SD_Final_Transformed_T1_F)

# Mean and Variance Relationship (T2_C)
ggplot(Complex_Data, aes(x = log(c(Mean_Transformed_T2_C)), y = log(SD_Final_Transformed_T2_C))) + facet_wrap(vars(Trait_Category)) +  
  geom_point() + labs(y = "log(SD)", x = "log(mean)") + theme_bw()# + ylim(min = -10, max = 10) 

Graph_Points3 <- Complex_Data %>% select("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T2", 
                                        "Trait_ID", "Trait_Category", "Mean_Transformed_T2_C", 
                                        "SD_Final_Transformed_T2_C") %>% 
                                  filter(`Trait_Category` == "Physiological")
Graph_Points3$logmean <- log(Graph_Points3$Mean_Transformed_T2_C)
Graph_Points3$logsd <- log(Graph_Points3$SD_Final_Transformed_T2_C)

# Mean and variance Relationship (T2_F)
ggplot(Complex_Data, aes(x = log(c(Mean_Transformed_T2_F)), y = log(SD_Final_Transformed_T2_F))) + facet_wrap(vars(Trait_Category)) +  
  geom_point() + labs(y = "log(SD)", x = "log(mean)") + theme_bw()# + ylim(min = -10, max = 10) 

Graph_Points4 <- Complex_Data %>% select("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T2", 
                                         "Trait_ID", "Trait_Category", "Mean_Transformed_T2_F", 
                                         "SD_Final_Transformed_T2_F") %>% 
                                  filter(`Trait_Category` == "Physiological")
Graph_Points4$logmean <- log(Graph_Points4$Mean_Transformed_T2_F)
Graph_Points4$logsd <- log(Graph_Points4$SD_Final_Transformed_T2_F)

##### Sanity Checks #####

Development_Exposure_Category_Check <- Complex_Data %>% select("Developmental_Exposure_Time_Category") %>% 
  filter(!Developmental_Exposure_Time_Category %in% c("Embryo", "Larvae", 
                                                      "Pupae", "Juvenile")) %>% 
  na.omit() %>% nrow()

Acclimation_Stage__Category_Check <- Complex_Data %>% select("Acclimation_Life.History_Stage_Category") %>% 
  filter(!`Acclimation_Life.History_Stage_Category` %in%
           c("Embryo", "Larvae", "Pupae", "Juvenile", "Adult")) %>% 
  na.omit() %>% nrow()

Behavioural_Check <- Complex_Data %>% filter(Trait_Category == "Behavioural") %>% select("Measurement") %>%
  filter(!Measurement %in% c("Anlge of Expoloration", "Exploration Speed", 
                             "Food Consumption", "Latency", "Righting Time", 
                             "Tortuosity")) %>% nrow()

Biochemical_Check <- Complex_Data %>% filter(Trait_Category == "Biochemical Assay") %>% select("Measurement") %>% 
  filter(!Measurement %in% c("AChE Activity", "Alkadienes", "Alkanes", "Alkenes", 
                             "Ash Content", "Bleaching Index", "Calcium", 
                             "Catalase Activity", "Chlorine", "Cortisol", 
                             "CytP450 Activity", "Enzyme Activity", "Erythrocytes", 
                             "ETS Activity", "Fat Content", "Free Sugar Content", 
                             "Glucose", "Glycogen", "Haemocytes", "Haemoglobin", 
                             "Haemolymph Phenoloxidase Activity", "Hexokinase Activity", 
                             "Lipid Content", "MDA Levels", "Monounsaturated Fatty Acids", 
                             "NKA", "PO Activity", "Polyunsaturated Fatty Acids", 
                             "Potassium", "Protein Content", "Pyruvate Kinase Activity", 
                             "Saturated Fatty Acids", "SOD Activity", "Sodium", 
                             "Triglyceride")) %>% nrow()

Gene_Expression_Check <- Complex_Data %>% filter(Trait_Category == "Gene Expression") %>% select("Measurement") %>% 
  filter(!Measurement %in% c("agrp", "cck8", "COX", "crh", "CS", "Cyp19A1", 
                             "Dmrt1", "Foxl2", "gh", "ghrelin", "grp", "hsp70", 
                             "hsp90", "igf1", "igf1a", "igf2a", "Kdm6b", "LDH", 
                             "mch1", "mch2", "npy", "orexin", "Sox4", "Sox9")) %>% nrow()

LH_Traits_Check <- Complex_Data %>% filter(Trait_Category == "Life-History Traits") %>% select("Measurement") %>% 
  filter(!Measurement %in% c("Brood Duration", "Brood Size", "Budding Rate", 
                             "Calcification Rate", "Cannibalism", "Development Time", 
                             "Egg Production", "Fecundity", "Gonotrophic Cycles", 
                             "Individual Growth Rate", "Longevity", 
                             "Ovarian Development Grade", "Reproductive Rate", 
                             "Retained Unfertilized Eggs", "Spermatophore Production")) %>% nrow()

Morphology_Check <- Complex_Data %>% filter(Trait_Category == "Morphology") %>% select("Measurement") %>% 
  filter(!Measurement %in% c("Abdomen Length", "Background brightness of eyespot", 
                             "Body Condition", "Carapace Length", "Carapace Width", 
                             "Collar Volume", "Femur Length", "Forelimb Length", 
                             "Forewing Length", "Forewing Spot Area", "Head Depth", 
                             "Head Length", "Head Width", "Hindlimb Length", 
                             "Hindwing Band", "Hindwing Band Area", "Hindwing Spot Area", 
                             "Length", "Lip Volume", "Mass", "Plastron Colouration", 
                             "Plastron Length", "Relative eyespot size", "Tail Length", 
                             "Thorax Size", "Tibia Size", "Total Boutons", "Width", 
                             "Wing Area", "Wing Aspect", "Wing Centriod", "Wing Length", 
                             "Wing Load", "Wing Size", "Yolk Dry Mass", 
                             "Yolk Free Mass", "Yolk Sac Volume")) %>% nrow()

Physiological_Check <- Complex_Data %>% filter(Trait_Category == "Physiological") %>% select("Measurement") %>% 
  filter(!Measurement %in% c("Apparent Digestability Coefficient", 
                             "Byssal Thread Regeneration", "Byssal Thread Strength", 
                             "Cellular Energy Allocation", "Desiccation Tolerance", 
                             "Energy Content", "Erythrocytic Haemoglobin", 
                             "Erythrocytic Volume", "Excretion Rate", "Exuvia Energy", 
                             "Faeces Energy", "Feeding Efficiency", "Filtration Rate", 
                             "Flight Distance", "Flight Velocity", "Haematocrit", 
                             "Heart Rate", "Immune Defense", "Locomotor Performance", 
                             "Maximum Buccal Pressure", "Maximum Excitation Pressure", 
                             "Metabolic Rate", "Minimum Opercular Pressure", 
                             "Muscular Strength", "Osmolality", "PHA Response", 
                             "Plasma Protein Levels", "Red Cell Size", 
                             "Starvation Tolerance", "Stroke Volume", "Swelling", 
                             "Travel Duration", "Ventilatory Flow", 
                             "Ventilatory Rate", "Water Content")) %>% nrow()

Population_Check <- Complex_Data %>% filter(Trait_Category == "Population") %>% select("Measurement") %>%
  filter(!Measurement %in% c("Gender", "Growth Rate", "Hatchling Abnormalities", 
                             "Mortality", "Sex Reversal", "Survival")) %>% nrow()

Development_Exposure_Check <- Complex_Data %>% select("Developmental_Exposure_Time") %>%
  filter(!Developmental_Exposure_Time %in% c("Egg Incubation", "Embryo Development", 
                                             "Final Larvae Period", "Larvae Period", 
                                             "Stage 0 - 21", "Stage 0 - 24", 
                                             "Third Instar Larvae Period", 
                                             "Until Gosner Stage 35-37", "Pupae Period", 
                                             "Juvenile Period", "Neonate Period", 
                                             "Nymphal Period")) %>% na.omit() %>% nrow()

Acclimation_Stage_check <- Complex_Data %>% select("Acclimation_Life.History_Stage") %>% 
  filter(!`Acclimation_Life.History_Stage` %in% c("Egg", "Embryo", "Cocoon", "Fourth Instar", 
                                                  "Larvae", "Tadpole", "Pupae", "Juvenile", 
                                                  "Adult")) %>% na.omit() %>% nrow()


