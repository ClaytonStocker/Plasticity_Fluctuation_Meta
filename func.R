#### ------------------------------------- ####
##     Functions
#### ------------------------------------- ####

#| @title PRRD
#| @description Calculate the PRRD
#| @param t1 Temperature 1
#| @param t2 Temperature 2
#| @param t1_c Mean for temperature 1 for constant group
#| @param t2_c Mean for temperature 2 for constant group
#| @param t1_f Mean for temperature 1 for fluctuating group
#| @param t2_f Mean for temperature 2 for fluctuating group
#| @param sd_t1_c Standard deviation for temperature 1 for constant group
#| @param sd_t2_c Standard deviation for temperature 2 for constant group
#| @param sd_t1_f Standard deviation for temperature 1 for fluctuating group
#| @param sd_t2_f Standard deviation for temperature 2 for fluctuating group
#| @param n_t1_c Sample size for temperature 1 for constant group
#| @param n_t2_c Sample size for temperature 2 for constant group
#| @param n_t1_f Sample size for temperature 1 for fluctuating group		
#| @param n_t2_f Sample size for temperature 2 for fluctuating group
#| @param type Type of output: 'ef' for effect size, 'v' for sampling variance
#| @return PRRD or its sampling variance
#| @examples
#| PRRD(t1 = 20, t2 = 30, t1_c = 5, t2_c = 10, t1_f = 5, t2_f = 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be 0
#| PRRD(t1 = 20, t2 = 30, t1_c = 10, t2_c = 5, t1_f = 5, t2_f = 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be negative because c is smaller than f
#| PRRD(t1 = 20, t2 = 30, t1_c = 5, t2_c = 10, t1_f = 10, t2_f = 5, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'ef') # Should be negative because c is bigger than f
#| PRRD(20, 30, 5, 10, 5, 10, 0.1, 0.2, 0.1, 0.2, 5, 10, 5, 10, type = 'v')
#| PRRD(20, 30, 5, 10, 5, 10, 3, 6, 2, 8, 2, 3, 2, 3, type = 'v')  # Should be bigger than l57
#|  check with single row, need to check direction etc
#| PRRD(t2 = 23.7, t1 = 17.7, 
#|   	 t1_c = 33.23927393, 
#| 	 t2_c = 17.81023102, 
#| 	 t1_f = 26.63861386, 
#| 	 t2_f = 16.82013201, 
#| 	 sd_t1_c = 2.25026431, 
#| 	 sd_t2_c = 1.666842873, 
#| 	 sd_t1_f = 2.833685746, 
#| 	 sd_t2_f = 2.833685746, 
#| 	 n_t1_c = 50, n_t2_c = 50, n_t1_f = 50, n_t2_f = 50, type = 'v')  

PRRD <- function(t1, t2, t1_c, t2_c, t1_f, t2_f, sd_t1_c, sd_t2_c, sd_t1_f, sd_t2_f, 
                 n_t1_c, n_t2_c, n_t1_f, n_t2_f, type = c('ef', 'v')) {
  
  type <- match.arg(type)
  
  # Check for invalid values
  if (any(c(t1_c, t2_c, t1_f, t2_f) <= 0)) {
    stop("All mean values must be positive to compute log response ratios.")
  }
  
  inv_temp_diff2 <- 1 / (t2 - t1)^2  # Precompute for efficiency
  
  # Effect size calculation
  if (type == 'ef') { 
    prrd <- ((log(t2_f) - log(t1_f)) - (log(t2_c) - log(t1_c))) / (t2 - t1)
    return(prrd)
  }
  
  # Sampling variance calculation
  if (type == 'v') { 
    v_prrd <- inv_temp_diff2 * (
      (sd_t1_c^2 / (n_t1_c * t1_c^2)) + 
      (sd_t2_c^2 / (n_t2_c * t2_c^2)) + 
      (sd_t1_f^2 / (n_t1_f * t1_f^2)) + 
      (sd_t2_f^2 / (n_t2_f * t2_f^2))
    )
    return(v_prrd)
  }
}

round_df <- function(df, digits = 2) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) round(col, digits) else col
  })
  return(df)
}


# Build table with results for supp
table_results <- function(model, group = NULL, study_name = "Study_ID", species_name = "Scientific_Name"){
  
  dat = model$data
  
  if(is.null(group)){
    tab <- dat %>% summarise(study = length(unique(!!sym(study_name))),
                             species = length(unique(!!sym(species_name))),
                             k = n())
  } else {
    tab <- dat %>% group_by(!!sym(group)) %>% summarise(study = length(unique(!!sym(study_name))),
                                                     species = length(unique(!!sym(species_name))),
                                                     k = n())
  }
  
  modelvals <- round_df(data.frame(
    "Estimate" = as.vector(model$b),
    "CI Low" = as.vector(model$ci.lb), 
    "CI High" = as.vector(model$ci.ub), 
    "df" = as.vector(model$ddf), 
    "p-value" = as.vector(model$pval), 
    check.names = FALSE, row.names = NULL), 3)
  
  return(cbind(tab, modelvals))
}
