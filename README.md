# Temperature variability does not impact the capacity for phenotypic plasticity in ectotherms – a meta-analysis

This repository contains data and R code for the associated manuscript.  

Clayton W. Stocker, Stephanie M. Bamford, Miki Jahn, Geoffrey P. F. Mazué, Amanda K. Pettersen, Daniel Ritchie, Alexander Rubin, Daniel W.A. Noble, Frank Seebacher. (2025) Temperature variability does not impact the capacity for phenotypic plasticity in ectotherms – a meta-analysis.

Please email clayton.stocker@sydney.edu.au for any inquiries. 

## Analyses

All code for analyses reproducing the results and figures is found in `Analysis_Final.R` script. The code also generates supplementary tables and results.

## Data
The data is the `Complex_Final_Data.csv` file. This file contains the data used in the meta-analysis. The meta-data that describes the columns are as follows:

| Column Name | Description |
------------- | ------------
Effect_Size_ID	| Unique identifiers for individual effect sizes.
Study_ID	| Unique identifiers for each paper.
Species_ID	| Identifier for each species that a study investigates.
Treatment_ID_T1	| Identifier for the studies’ first constant and fluctuating temperature treatment pair (fluctuating mean matching constant temperature). 
Treatment_ID_T2	| Identifier for the studies’ second constant and fluctuating temperature treatment pair (fluctuating mean matching constant temperature). 
Trait_ID	| Identifier for each phenotypic trait a study measures.
First_Author | Initials and surname for the first author of the study.
Title	| Title of the study.
Year	| Year of publication.
Journal	| Name of the journal the study was published in.
Journal_Impact_Factor-2021| 	Journal impact factor as of 2021 (most recent records at time of data collection).
Kingdom	| Kingdom for each species that a study investigates.
Phylum	| Phylum for each species that a study investigates.
Class	| Class for each species that a study investigates.
Order	| Order for each species that a study investigates.
Family	| Family for each species that a study investigates.
Scientific_Name	| Binomial nomenclature for species.
Ecosystem	| Ecosystem the species is naturally observed (Aquatic and Terrestrial). Amphibians considered aquatic.
Plasticity_Mechanism	| Mechanism of plasticity the treatments were exposed to (Acclimation or Developmental Plasticity).
Developmental_Exposure_Time_Category	| Categorisation of Developmental Exposure Time.
Developmental_Exposure_Time	| The period of exposure for treatments imposed during development. 
Acclimation_Exposure_Time	| The duration of exposure for acclimation treatments. 
Exposure_Units	|  Units of Acclimation Exposure Time (Days).
T1_constant	| Temperature of the first constant temperature treatment. 
T1_fluctuation	| Mean temperature of the first fluctuating temperature treatment.
T2_constant	| Temperature of the second constant temperature treatment.
T2_fluctuation	| Mean temperature of the second fluctuating temperature treatment.
High	| Whether the first or second constant and fluctuating temperature treatment pair is around the higher temperature (T1 or T2). 
Fluctuation_Magnitude	| Amplitude of the two fluctuating temperature treatments.
Fluctuation_Category	| Type of fluctuations imposed (Sinusoidal, Alternating, Stepwise, Stochastic). 
Fluctuation_Period	| Period of one fluctuation oscillation.
Fluctuation_Unit	| Units of Fluctuation Period (Days).
Number_Of_Fluctuations	| Acclimation Exposure Time/Fluctuation Period for acclimation treatments. 
Acclimation_Life-History_Stage	| Life-history stage of organisms for acclimation treatments. 
Acclimation_Life-History_Stage_Category	| Categorisation of Acclimation Life-history Stages. 
Trait_Category	| Categorisation of Trait Measurements.
Measurement	| Phenotypic traits measured following treatment exposure. 
Trait_Unit	| Units for measurements. 
Sex	| Sex of the organisms being investigated (Both, Female or Male). NA = sex not specified.
Performance_Curve	| Whether a performance curve was recorded in the study (Yes, No). 
Complex_Design	| Whether a comparison between constant and fluctuating treatments was made at multiple temperatures (Yes). 
Species_Overlap	| Identifier for effect sizes that use the same species within a study
Animal_Overlap_T1_Constant	| Identifier for effect sizes where the constant treatment was conducted on the same animals within a study (first constant and fluctuating temperature treatment pair). 
Animal_Overlap_T1_Fluctuation	| Identifier for effect sizes where the fluctuating treatment was conducted on the same animals within a study (first constant and fluctuating temperature treatment pair). 
Animal_Overlap_Trait_T1	| Identifier for effect sizes that measure phenotypic traits on the same animal within a study (first constant and fluctuating temperature treatment pair). 
Animal_Overlap_T2_Constant	| Identifier for effect sizes where the constant treatment was conducted on the same animals within a study (second constant and fluctuating temperature treatment pair). 
Animal_Overlap_T2_Fluctuation	| Identifier for effect sizes where the fluctuating treatment was conducted on the same animals within a study (second constant and fluctuating temperature treatment pair). 
Animal_Overlap_Trait_T2	| Identifier for effect sizes that measure phenotypic traits on the same animal within a study (second constant and fluctuating temperature treatment pair). 
Animal_Code	| Study ID: Species Overlap: Animal Overlap T1 Constant: Animal Overlap T1 Fluctuation: Animal Overlap T2 | Constant: Animal Overlap T2 Fluctuation: Animal Overlap Trait T1; Animal Overlap Trait T2
Shared_Animal_Number| 	Unique identifier for shared animal codes across effect sizes. 
Shared_Control_T1	| Identifier for effect sizes that use the same constant temperature treatment within a study (first constant and fluctuating temperature treatment pair). 
Shared_Control_T2	| Identifier for effect sizes that use the same constant temperature treatment within a study (second constant and fluctuating temperature treatment pair). 
Shared_Control_Code	| Study ID: Species Overlap: Shared Control T1: Shared Control T2: Trait ID. 
Shared_Control_Number| 	Unique identifier for shared control codes across effect sizes. 
n_T1_C	| Sample size of the constant treatment (first constant and fluctuating temperature treatment pair). 
Mean_T1_C	| Mean response of the constant treatment (first constant and fluctuating temperature treatment pair).
SD_Final_T1_C	| Standard deviation of the constant treatment (first constant and fluctuating temperature treatment pair).
n_T1_F	| Sample size of the fluctuating treatment (first constant and fluctuating temperature treatment pair). 
Mean_T1_F	| Mean response of the fluctuating treatment (first constant and fluctuating temperature treatment pair).
SD_Final_T1_F	| Standard deviation of the fluctuating treatment (first constant and fluctuating temperature treatment pair).
n_T2_C	| Sample size of the constant treatment (second constant and fluctuating temperature treatment pair). 
Mean_T2_C	| Mean response of the constant treatment (second constant and fluctuating temperature treatment pair).
SD_Final_T2_C	| Standard deviation of the constant treatment (second constant and fluctuating temperature treatment pair).
n_T2_F	| Sample size of the fluctuating treatment (second constant and fluctuating temperature treatment pair). 
Mean_T2_F	| Mean response of the fluctuating treatment (second constant and fluctuating temperature treatment pair).
SD_Final_T2_F	| Standard deviation of the fluctuating treatment (second constant and fluctuating temperature treatment pair).
Percentage_Transformation_T1	| Whether the recorded means for the first constant and fluctuating temperature treatment pair were recorded as a percentage (Yes, No).
Proportion_Transformation_T1	| Whether the recorded means for the first constant and fluctuating temperature treatment pair were recorded as a proportion (Yes, No). 
In_Transformation_T1	| Whether the recorded means for the first constant and fluctuating temperature pair were recorded as the natural log (Yes, No). 
Percentage_Transformation_T2	| Whether the recorded means for the second constant and fluctuating temperature treatment pair were recorded as a percentage (Yes, No).
Proportion_Transformation_T2| 	Whether the recorded means for the second constant and fluctuating temperature treatment pair were recorded as a proportion (Yes, No). 
In_Transformation_T2	| Whether the recorded means for the second constant and fluctuating temperature pair were recorded as the natural log (Yes, No). 
Mean_Transformed_T1_C	| Constant treatment means transformed for percentages, proportions, or natural logs, and with a constant of 0.5 added (first constant and fluctuating temperature treatment pair). 
SD_Final_Transformed_T1_C	| Constant treatment standard deviations transformed for percentages, proportions, or natural logs, and with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
Mean_Transformed_T1_F	| Fluctuating treatment means transformed for percentages, proportions, or natural logs, and with a constant of 0.5 added (first constant and fluctuating temperature treatment pair). 
SD_Final_Transformed_T1_F | 	Fluctuating treatment standard deviations transformed for percentages, proportions, or natural logs, and with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
Mean_Transformed_T2_C	|  Constant treatment means transformed for percentages, proportions, or natural logs, and with a constant of 0.5 added (second constant and fluctuating temperature treatment pair). 
SD_Final_Transformed_T2_C | Constant treatment standard deviations transformed for percentages, proportions, or natural logs, and with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
Mean_Transformed_T2_F	| Fluctuating treatment means transformed for percentages, proportions, or natural logs, and with a constant of 0.5 added (second constant and fluctuating temperature treatment pair). 
SD_Final_Transformed_T2_F | 	Fluctuating treatment standard deviations transformed for percentages, proportions, or natural logs, and with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
Mean_T1_C_Add	|  Untransformed constant treatment means with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
SD_Final_T1_C_Add	| Untransformed constant standard deviations with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
Mean_T1_F_Add	| Untransformed fluctuating treatment means with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
SD_Final_T1_F_Add	| Untransformed fluctuating treatment standard deviations with a constant 0.5 added (first constant and fluctuating temperature treatment pair). 
Mean_T2_C_Add| Untransformed constant treatment means with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
SD_Final_T2_C_Add| Untransformed constant standard deviations with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
Mean_T2_F_Add| Untransformed fluctuating treatment means with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
SD_Final_T2_F_Add| Untransformed fluctuating treatment standard deviations with a constant 0.5 added (second constant and fluctuating temperature treatment pair). 
InRR| Plasticity response ratio difference. Represents the difference in plastic responses in constant and fluctuating temperatures, standardised to a one-degree change in treatment temperature. 
InRR_Transformed| Plasticity response ratio difference transformed to account for reciprocal transformation. 
v_InRR| Plasticity response ratio difference sampling variance.
InRR_Untransformed| Plasticity response ratio difference calculated from untransformed data. 
v_InRR_Untransformed| Plasticity response ratio difference sampling variance calculated from untransformed data. 
Year_Z| Z-transformed year of publication.
Precision| Z-transformed inverse of plasticity response ratio difference sampling variance.
