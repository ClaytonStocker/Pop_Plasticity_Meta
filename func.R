##---------------------------------------------------------------------------##
## Functions for effect size calculations, data analysis and data exploration
##---------------------------------------------------------------------------##

#' @title sdp
#' @description Calculate the pooled standard deviation
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param sd3 Standard deviation of group 3
#' @param sd4 Standard deviation of group 4
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4

sdp <- function(sd1, n1, sd2, n2, sd3, n3, sd4, n4){
  sdp = sqrt( (((n1-1)*sd1^2) + ((n2-1)*sd2^2) + ((n3-1)*sd3^2) + ((n4-1)*sd4^2)) / (n1 + n2 + n3 + n4) )
  return(sdp)
}

#' @title gI
#' @description Caluclate the interaction based standardised mean difference with small sample correction when ignoring nuisance heterogeneity. 
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param sdp Pooled standard deviation
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' 
gI <- function(x1, x2, x3, x4, n1, n2, n3, n4, sdp){
  df = sum(c(n1, n2, n3, n4)) - 4
  
  J = 1 - (3 / ((4*df) - 1)) 
  
  gI = (((x1 - x2) - (x3 - x4)) / sdp) * J
  return(gI)
}


#' @title gI_cor
#' @description Caluclate the interaction based standardised mean difference with small sample correction while accounting for nuisance heterogeneity. 
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param sdp Pooled standard deviation
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity.
#' 
gI_cor <- function(x1, x2, x3, x4, n1, n2, n3, n4, t1, t2, sdp){
  df = sum(c(n1, n2, n3, n4)) - 4
  
  J = 1 - (3 / ((4*df) - 1)) 
  
  gI_cor = (((x1 - x2) - (x3 - x4)) / (sdp*(t2-t1))) * J
  return(gI_cor)
}

#' @title v_gI_cor
#' @description Calculate the sampling variance for the interaction based standardised mean difference with small sample correction when accounting for nuisance heterogeneity. 
#' @param gI_cor Interaction based standardised mean difference with small sample correction and nuisance heterogenity aaccounted for
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity.

v_gI_cor <- function(gI_cor, n1, n2, n3, n4, t1, t2){
  v_gI_cor = (((gI_cor^2) / (2*sum(c(n1, n2, n3, n4)))) + (1/n1) + (1/n2) + (1/n3) + (1/n4))*((1/t2-t1)^2)
  return(v_gI_cor)
}

#' @title v_gI
#' @description Calculate the sampling variance for the interaction based standardised mean difference with small sample correction when ignoring nuisance heterogeneity. 
#' @param gI Interaction based standardised mean difference with small sample correction
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4

v_gI <- function(gI, n1, n2, n3, n4){
  v_gI = ((gI^2) / (2*sum(c(n1, n2, n3, n4)))) + (1/n1) + (1/n2) + (1/n3) + (1/n4)
  return(v_gI)
}


#' @title SMD_int
#' @description Caluclate the interaction based standardised mean difference with small sample correction and associated sampling variance when correcting for nuisance heterogenity 
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param sd3 Standard deviation of group 3
#' @param sd4 Standard deviation of group 4
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity. Only needed when cor = TRUE.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity. Only needed when cor = TRUE
#' @param cor Logical. Whether to correct for nuisance heterogeneity (if possible) ('TRUE') or not ('FALSE'). Defaults to TRUE.
#' 
SMD_int <- function(x1, x2, x3, x4, sd1, sd2, sd3, sd4, n1, n2, n3, n4, t1, t2, cor = TRUE){
  if(cor){
         SDp <- sdp(sd1, n1, sd2, n2, sd3, n3, sd4, n4)
      gI_cor <- gI_cor(x1, x2, x3, x4, n1, n2, n3, n4, t1, t2, SDp)
    v_gI_cor <- v_gI_cor(gI_cor, n1, n2, n3, n4, t1, t2)
    
    return(data.frame(gI_cor = gI_cor, v_gI_cor = v_gI_cor))
  
    } else {
     SDp <- sdp(sd1, n1, sd2, n2, sd3, n3, sd4, n4)
      gI <- gI(x1, x2, x3, x4, n1, n2, n3, n4, SDp)
    v_gI <- v_gI(gI, n1, n2, n3, n4)
    return(data.frame(gI = gI, v_gI = v_gI))
  }
}

#' @title lnRR_i
#' @description Calculate the interaction based log response ratio when correcting for nuisance heterogenity or not
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity. Only needed when cor = TRUE.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity. Only needed when cor = TRUE
#' @param cor Logical. Whether to correct for nuisance heterogeneity (if possible) ('TRUE') or not ('FALSE'). Defaults to TRUE.
#' 
lnRR_i <- function(x1, x2, x3, x4, t1, t2, cor = TRUE){
  if(cor){
    rr_i_cor <- ((log(x2) - log(x1)) - (log(x4) - log(x3))) / (t2-t1)
    return(rr_i_cor)
  } else{
    rr_i <- (log(x2) - log(x1)) - (log(x4) - log(x3))
    return(rr_i)
  }
   
}

#' @title v_lnRR_i
#' @description Calculate the interaction based log response ratio sampling variance when correcting for nuisance heterogenity or not
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param sd3 Standard deviation of group 3
#' @param sd4 Standard deviation of group 4
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity. Only needed when cor = TRUE.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity. Only needed when cor = TRUE
#' @param cor Logical. Whether to correct for nuisance heterogeneity (if possible) ('TRUE') or not ('FALSE'). Defaults to TRUE.
#' 
v_lnRR_i <- function(x1, x2, x3, x4, sd1, sd2, sd3, sd4, n1, n2, n3, n4, t1, t2, cor = TRUE){
  if(cor){
    v_rr_i_cor <- (1 / (t2-t1))^2 + (sd1^2 / (x1^2*n1)) + (sd2^2 / (x2^2*n2)) + 
                    (sd3^2 / (x3^2*n3)) + (sd4^2 / (x4^2*n4))
    return(v_rr_i_cor)
  } else{
    v_rr_i <- (sd1^2 / (x1^2*n1)) + (sd2^2 / (x2^2*n2)) + (sd3^2 / (x3^2*n3)) + (sd4^2 / (x4^2*n4))
    return(v_rr_i)
  }
  
}

#' @title lnRR_int
#' @description Calculate the interaction based log response ratio and associated sampling variance when correcting for nuisance heterogenity or not
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param x3 Mean of group 3
#' @param x4 Mean of group 4
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param sd3 Standard deviation of group 3
#' @param sd4 Standard deviation of group 4
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param n3 Sample size of group 3
#' @param n4 Sample size of group 4
#' @param t1 Nuisance variable for group 1. For example, temperature or salinity. Only needed when cor = TRUE.
#' @param t2 Nuisance variable for group 2. For example, temperature or salinity. Only needed when cor = TRUE
#' @param cor Logical. Whether to correct for nuisance heterogeneity (if possible) ('TRUE') or not ('FALSE'). Defaults to TRUE.
#' 
lnRR_int <- function(x1, x2, x3, x4, sd1, sd2, sd3, sd4, n1, n2, n3, n4, t1, t2, cor = TRUE){
              lnrr_i <- lnRR_i(x1, x2, x3, x4, t1, t2, cor = cor)
            v_lnrr_i <- v_lnRR_i(x1, x2, x3, x4, sd1, sd2, sd3, sd4, n1, n2, n3, n4, t1, t2, cor = cor)
            
            if(cor){
              return(data.frame(lnrr_i_cor = lnrr_i, v_lnrr_i_cor = v_lnrr_i))
            } else{
              return(data.frame(lnrr_i = lnrr_i, v_lnrr_i = v_lnrr_i))
            }
}

#' @title tree_checks
#' @description Checks that the names of species in the dataset and tree tips match exactly. Identifies which taxa are not in data or tree etc. When no species exist in the data that are not found in the tree it also allows one to prune the tree.
#' @param data The dataframe with data on species and all other variables
#' @param tree A phylogenetic tree 
#' @param dataCol Column in the dataframe that hosts the species name data as a character string
#' @param type Whether to only implement 'checks' that compare species names in tree and data or whther to 'prune' the tree and return the pruned tree. Defaults to 'checks' first because pruning is only available when all species in data are also in tree.
#'
tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
  type = match.arg(type)
  # How many unique species exist in data and tree
  Numbers <- matrix(nrow = 2, ncol = 1)
  Numbers[1,1] <- length(unique(data[,dataCol])) 
  Numbers[2,1] <- length(tree$tip.label) 
  rownames(Numbers)<- c("Species in data:", "Species in tree:")
  # Missing species or species not spelt correct      
  species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
  species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )
  if(type == "checks"){
    return(list(SpeciesNumbers = data.frame(Numbers), 
                Species_InTree_But_NotData=species_list1, 
                Species_InData_But_NotTree=species_list2))
  }
  if(type == "prune"){
    if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
    return(ape::drop.tip(tree, species_list1))
  }
}


#' @title folded_norm
#' @description Converts a mean on the normal scale to a folded normal distribution (absolute)
#' @param mu The mean assuming normality
#' @param sigma The standard deviation

folded_norm <-  function(mu, sigma){
  dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)
}


#' @title mean_ci
#' @description Calculates the mean and 95% credible interval from an MCMC chain
#' @param mcmc MCMC chain of some statistic of interest

mean_ci <- function(mcmc){
      m <- mean(mcmc)      
     ci <- coda::HPDinterval(as.mcmc(mcmc))
  return(format(round(data.frame(Est. = m, `95% LCI` = ci[1], `95% UCI` = ci[2], check.names = FALSE), digits = 2), nsmall = 2))
}

#' @title i2
#' @description Calculates I2 for all random effects in a multilevel meta-analysis model
#' @param sd_mcmc MCMC chain of standard deviation estimates for all random effects from brms model
#' @param sd_mcmc Vector of sampling variances for effect sizes

i2 <- function(sd_mcmc, sv){
  # i2 for all random effects
  tmp <- apply(sd_mcmc, 2, function(x) (x^2 / rowSums(sd_mcmc^2))*100)
  val <- do.call(rbind, apply(tmp, 2, function(x) mean_ci(x)))
  
  # total i2
  sv <- sum(1 / sv) * (length(sv) - 1) / (sum(1 / sv)^2 - sum((1 / sv)^2))
  
  i2_tot <- format(round((rowSums(sd_mcmc^2) / (rowSums(sd_mcmc^2) + sv))*100, 2), 2)
  i2_tot <- mean_ci(as.numeric(i2_tot))
  
  # Clean up and return
     val <- rbind(val, i2_tot)
  row.names(val)[nrow(val)] <- "i2_total"
  
  return(val)
}