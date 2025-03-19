rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, tidyr, ggplot2, rotl, DescTools, stringr, ape)

##### Importing Raw Data - Initial Search #####
Data <- read_excel("./1.Data_Extraction/3.Finalised_Data_Spreadsheets/Complex_Data_Spreadsheet.xlsx", sheet = "Data")

# Adding Animal codes for Second Treatment
Treatment_Fix_T1 <- Data[, c("Study_ID", "Species_ID", "Treatment_ID_T1", "Trait_ID")]
Treatment_Fix_T2 <- Data[, c("Study_ID", "Species_ID", "Treatment_ID_T2", "Trait_ID")]

Treatment_Fix_Data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Complex_Data.csv")
Treatment_Fix_Data_T1 <- Treatment_Fix_Data %>% select("Study_ID", "Species_ID", "Treatment_ID", "Trait_ID", "Shared_Control")
Treatment_Fix_Data_T1 <- Treatment_Fix_Data_T1 %>% dplyr::rename("Treatment_ID_T1" = `Treatment_ID`,
                                                                 "Shared_Control_T1" = `Shared_Control`)
Treatment_Fix_T1 <- Treatment_Fix_T1 %>% left_join(Treatment_Fix_Data_T1, by = c("Study_ID", "Species_ID", "Treatment_ID_T1", "Trait_ID"))

Treatment_Fix_Data_T2 <- Treatment_Fix_Data %>% select("Study_ID", "Species_ID", "Treatment_ID", "Trait_ID", 
                                                       "Animal_Overlap_T1", "Animal_Overlap_T2", 
                                                       "Animal_Overlap", "Shared_Control")
Treatment_Fix_Data_T2 <- Treatment_Fix_Data_T2 %>% dplyr::rename("Treatment_ID_T2" = `Treatment_ID`, 
                                                                 "Animal_Overlap_T2_Constant" = `Animal_Overlap_T1`, 
                                                                 "Animal_Overlap_T2_Fluctuation" = `Animal_Overlap_T2`, 
                                                                 "Animal_Overlap_Trait_T2" = `Animal_Overlap`, 
                                                                 "Shared_Control_T2" = `Shared_Control`)
Treatment_Fix_T2 <- Treatment_Fix_T2 %>% left_join(Treatment_Fix_Data_T2, by = c("Study_ID", "Species_ID", "Treatment_ID_T2", "Trait_ID"))
  
Data <- Data %>% 
        left_join(Treatment_Fix_T1, by = c("Study_ID", "Species_ID", "Treatment_ID_T1", "Trait_ID")) %>%
        left_join(Treatment_Fix_T2, by = c("Study_ID", "Species_ID", "Treatment_ID_T2", "Trait_ID"))

# Re-ordering Data Spreadsheet
Complex_Data <- Data[, c("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T1", "Treatment_ID_T2", 
                         "Trait_ID", "First_Author", "Title", "Year", "Journal", 
                         "Journal_Impact_Factor-2021", "Kingdom", "Phylum", "Class", "Order", "Family", "Scientific_Name", "Ecosystem", 
                         "Plasticity_Mechanism", "Developmental_Exposure_Time_Category", "Developmental_Exposure_Time", 
                         "Acclimation_Exposure_Time", "Exposure_Units", "T1_constant", "T1_fluctuation", "T2_constant", "T2_fluctuation",
                         "Fluctuation_Magnitude", "Fluctuation_Category", "Fluctuation_Period", "Fluctuation_Unit", "Number_Of_Fluctuations", 
                         "Acclimation_Life-History_Stage", "Acclimation_Life-History_Stage_Category", "Trait_Category", "Measurement", "Trait_Unit", 
                         "Sex", "Performance_Curve", "Complex_Design",
                         "Species_Overlap", "Animal_Overlap_T1_Constant", "Animal_Overlap_T1_Fluctuation", "Animal_Overlap_Trait_T1", 
                         "Animal_Overlap_T2_Constant", "Animal_Overlap_T2_Fluctuation", "Animal_Overlap_Trait_T2", 
                         "Shared_Control_T1", "Shared_Control_T2",
                         "n_T1_C", "Mean_T1_C", "SD_Final_T1_C", "n_T1_F", "Mean_T1_F", "SD_Final_T1_F", 
                         "n_T2_C", "Mean_T2_C", "SD_Final_T2_C", "n_T2_F", "Mean_T2_F", "SD_Final_T2_F", 
                         "Percentage_Transformation_T1", "Proportion_Transformation_T1", "In_Transformation_T1", 
                         "Percentage_Transformation_T2", "Proportion_Transformation_T2", "In_Transformation_T2", 
                         "Mean_Transformed_T1_C", "SD_Final_Transformed_T1_C", "Mean_Transformed_T1_F", "SD_Final_Transformed_T1_F",
                         "Mean_Transformed_T2_C", "SD_Final_Transformed_T2_C", "Mean_Transformed_T2_F", "SD_Final_Transformed_T2_F")]

##### Adding Shared Animal and Shared Control Columns #####
# Shared Animal Code
Complex_Data <- Complex_Data %>% mutate(Animal_Code = paste(Study_ID, Species_Overlap, 
                                                            Animal_Overlap_T1_Constant, Animal_Overlap_T1_Fluctuation,
                                                            Animal_Overlap_T2_Constant, Animal_Overlap_T2_Fluctuation, 
                                                            `Animal_Overlap_Trait_T1`, `Animal_Overlap_Trait_T2`, sep = ":"))
Shared_Animal <- Complex_Data %>% select("Animal_Code") %>% unique()
Shared_Animal <- Shared_Animal %>% mutate(Shared_Animal_Number = c(1:length(Animal_Code)))
Complex_Data <- Complex_Data %>% left_join(Shared_Animal, by = "Animal_Code")

# Shared Control
Complex_Data <- Complex_Data %>% mutate(Shared_Control_Code = paste(Study_ID, Species_Overlap, 
                                                                    Shared_Control_T1, Shared_Control_T2, 
                                                                    Trait_ID, sep = ":"))

Shared_Control <- Complex_Data %>% select("Shared_Control_Code") %>% unique()
Shared_Control <- Shared_Control %>% mutate(Shared_Control_Number = c(1:length(Shared_Control_Code)))
Complex_Data <- Complex_Data %>% left_join(Shared_Control, by = "Shared_Control_Code")

##### Effect Size Calculations #####

# Identifying the Higher Temperature
Complex_Data <- Complex_Data %>% mutate(High = ifelse(T1_constant > T2_constant,"T1", "T2")) 

# Log Response Ratio
Complex_Data <- Complex_Data %>% mutate(InRR_Transformed = ifelse(High == "T2", 
                                                                 (log(Mean_Transformed_T2_F/Mean_Transformed_T1_F) - 
                                                                  log(Mean_Transformed_T2_C/Mean_Transformed_T1_C))/
                                                                     (T2_constant - T1_constant), 
                                                                 (log(Mean_Transformed_T1_F/Mean_Transformed_T2_F) - 
                                                                  log(Mean_Transformed_T1_C/Mean_Transformed_T2_C))/
                                                                     (T1_constant - T2_constant))) %>%
                                 mutate(v_InRR = ifelse(High == "T2", 
                                                        (1/(T2_constant-T1_constant))^2 +
                                                        (SD_Final_Transformed_T1_C^2/(n_T1_C * Mean_Transformed_T1_C^2)) + 
                                                        (SD_Final_Transformed_T1_F^2/(n_T1_F * Mean_Transformed_T1_F^2)) +
                                                        (SD_Final_Transformed_T2_C^2/(n_T2_C * Mean_Transformed_T2_C^2)) +
                                                        (SD_Final_Transformed_T2_F^2/(n_T2_F * Mean_Transformed_T2_F^2)), 
                                                        (1/(T1_constant-T2_constant))^2 +
                                                        (SD_Final_Transformed_T1_C^2/(n_T1_C * Mean_Transformed_T1_C^2)) + 
                                                        (SD_Final_Transformed_T1_F^2/(n_T1_F * Mean_Transformed_T1_F^2)) +
                                                        (SD_Final_Transformed_T2_C^2/(n_T2_C * Mean_Transformed_T2_C^2)) +
                                                        (SD_Final_Transformed_T2_F^2/(n_T2_F * Mean_Transformed_T2_F^2))))

# Log Response Ratio (Untransformed)
Complex_Data <- Complex_Data %>% 
                mutate(Mean_T1_C_Add = Mean_T1_C+ 0.5) %>%
                mutate(SD_Final_T1_C_Add = SD_Final_T1_C+ 0.5) %>% 
                mutate(Mean_T1_F_Add = Mean_T1_F + 0.5) %>% 
                mutate(SD_Final_T1_F_Add = SD_Final_T1_F + 0.5) %>%
                mutate(Mean_T2_C_Add = Mean_T2_C+ 0.5) %>%
                mutate(SD_Final_T2_C_Add = SD_Final_T2_C+ 0.5) %>% 
                mutate(Mean_T2_F_Add = Mean_T2_F + 0.5) %>% 
                mutate(SD_Final_T2_F_Add = SD_Final_T2_F + 0.5)

Complex_Data <- Complex_Data %>% mutate(InRR_Untransformed = ifelse(High == "T2", 
                                                                   (log(Mean_T2_F_Add/Mean_T1_F_Add) - 
                                                                    log(Mean_T2_C_Add/Mean_T1_C_Add))/
                                                                       (T2_constant - T1_constant), 
                                                                   (log(Mean_T1_F_Add/Mean_T2_F_Add) - 
                                                                    log(Mean_T1_C_Add/Mean_T2_C_Add))/
                                                                       (T1_constant - T2_constant))) %>%
                                 mutate(v_InRR_Untransformed = ifelse(High == "T2", 
                                                                     (1/(T2_constant-T1_constant))^2 +
                                                                     (SD_Final_T1_C_Add^2/(n_T1_C * Mean_T1_C_Add^2)) + 
                                                                     (SD_Final_T1_F_Add^2/(n_T1_F * Mean_T1_F_Add^2)) +
                                                                     (SD_Final_T2_C_Add^2/(n_T2_C * Mean_T2_C_Add^2)) +
                                                                     (SD_Final_T2_F_Add^2/(n_T2_F * Mean_T2_F_Add^2)), 
                                                                     (1/(T1_constant-T2_constant))^2 +
                                                                     (SD_Final_T1_C_Add^2/(n_T1_C * Mean_T1_C_Add^2)) + 
                                                                     (SD_Final_T1_F_Add^2/(n_T1_F * Mean_T1_F_Add^2)) +
                                                                     (SD_Final_T2_C_Add^2/(n_T2_C * Mean_T2_C_Add^2)) +
                                                                     (SD_Final_T2_F_Add^2/(n_T2_F * Mean_T2_F_Add^2))))

##### Transformations for Sensitivity Analysis #####
Complex_Data <- Complex_Data %>%
                mutate(Year_Z = scale(Year)) %>%
                mutate(Precision = scale(1/v_InRR))

##### Changing Column Order and Creating CSV #####
Complex_Data <- Complex_Data[, c("Effect_Size_ID", "Study_ID", "Species_ID", "Treatment_ID_T1", "Treatment_ID_T2", 
                                 "Trait_ID", "First_Author", "Title", "Year", "Journal", 
                                 "Journal_Impact_Factor-2021", "Kingdom", "Phylum", "Class", "Order", "Family", "Scientific_Name", "Ecosystem", 
                                 "Plasticity_Mechanism", "Developmental_Exposure_Time_Category", "Developmental_Exposure_Time", 
                                 "Acclimation_Exposure_Time", "Exposure_Units", "T1_constant", "T1_fluctuation", "T2_constant", "T2_fluctuation", "High", 
                                 "Fluctuation_Magnitude", "Fluctuation_Category", "Fluctuation_Period", "Fluctuation_Unit", "Number_Of_Fluctuations", 
                                 "Acclimation_Life-History_Stage", "Acclimation_Life-History_Stage_Category", "Trait_Category", "Measurement", "Trait_Unit", 
                                 "Sex", "Performance_Curve", "Complex_Design",
                                 "Species_Overlap", "Animal_Overlap_T1_Constant", "Animal_Overlap_T1_Fluctuation", "Animal_Overlap_Trait_T1", 
                                 "Animal_Overlap_T2_Constant", "Animal_Overlap_T2_Fluctuation", "Animal_Overlap_Trait_T2", "Animal_Code", "Shared_Animal_Number",  
                                 "Shared_Control_T1", "Shared_Control_T2", "Shared_Control_Code", "Shared_Control_Number", 
                                 "n_T1_C", "Mean_T1_C", "SD_Final_T1_C", "n_T1_F", "Mean_T1_F", "SD_Final_T1_F", 
                                 "n_T2_C", "Mean_T2_C", "SD_Final_T2_C", "n_T2_F", "Mean_T2_F", "SD_Final_T2_F", 
                                 "Percentage_Transformation_T1", "Proportion_Transformation_T1", "In_Transformation_T1", 
                                 "Percentage_Transformation_T2", "Proportion_Transformation_T2", "In_Transformation_T2", 
                                 "Mean_Transformed_T1_C", "SD_Final_Transformed_T1_C", "Mean_Transformed_T1_F", "SD_Final_Transformed_T1_F",
                                 "Mean_Transformed_T2_C", "SD_Final_Transformed_T2_C", "Mean_Transformed_T2_F", "SD_Final_Transformed_T2_F", 
                                 "Mean_T1_C_Add", "SD_Final_T1_C_Add", "Mean_T1_F_Add", "SD_Final_T1_F_Add", 
                                 "Mean_T2_C_Add", "SD_Final_T2_C_Add", "Mean_T2_F_Add", "SD_Final_T2_F_Add", 
                                 "InRR_Transformed", "v_InRR", "InRR_Untransformed", "v_InRR_Untransformed", "Year_Z", "Precision")]

write.csv(Complex_Data, file = "./3.Data_Analysis/2.Outputs/Data/Complex_Final_Data.csv", row.names = FALSE)

##### Creating Phylogeneitc Tree #####
# Build the tree:
Myspecies <- str_sort(as.character(unique(Complex_Data$Scientific_Name)))
Taxa <- tnrs_match_names(names = Myspecies, include_suppressed = TRUE)
tree <- tol_induced_subtree(ott_ids = Taxa[["ott_id"]], label_format = "name")

# Resolving polytomies at random
binary.tree <- multi2di(tree, random=T)
plot(binary.tree, node.color = "#183357")

# Exporting Tree
write.tree(binary.tree, "./3.Data_Analysis/2.Outputs/Phylogeny/Complex_tree")
