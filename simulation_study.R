


library(parallelMap)

trial_design <- readRDS("repStudy_trial_design.rds")
dat_complete <- readRDS(url("http://sz.hochschule-stralsund.de/jclassif/repStudy_dat_complete.rds"))
dat_incomplete <- readRDS(url("http://sz.hochschule-stralsund.de/jclassif/repStudy_dat_incomplete.rds"))

# names of the 108 data situations
td_unique <- stringr::str_remove_all(paste(trial_design[1,-c(1,5)], collapse = "_"),"\\.")
for(rn_td in as.numeric(rownames(unique(trial_design[,-c(1,5)])))[-1]){
  td_unique <- c(td_unique, stringr::str_remove_all(paste(trial_design[rn_td,-c(1,5)], collapse = "_"),"\\."))
}

par_fun <- function(n) {
  
  # Set specifications
  nO = trial_design[n,"nO"]; nC = trial_design[n, "nC"]; nV = trial_design[n, "nV"]; symm = trial_design[n, "symm"]
  fac_prop <- trial_design[n, "fac_prop"]; fac_lev <- trial_design[n, "fac_lev"]; overlap <- trial_design[n, "overlap"]
  missing_type <- as.character(trial_design[n, "missing_type"]); perc_miss <- trial_design[n, "perc_miss"]
  nstarts <- trial_design[n, "nstarts"]; N <- trial_design[n, "N"]
  
  #objects for evaluation:
  eval_objects <- c("ignore", "kpod_in", "kpod_ex", "mice_pooling")
  
  #evaluate computation time for imputation and clustering
  eval_time <- matrix(rep(NA, length(eval_objects)*N), nrow=N); colnames(eval_time) <- eval_objects
  
  #evaluate total sum of within cluster distances of origin data and the different cluster partitions
  eval_total_within_dists <- matrix(rep(NA, length(eval_objects)*N), nrow=N); colnames(eval_total_within_dists) <- eval_objects
  #evaluate the distance between kpres_sil$kp_obj and the different centers 
  eval_centers <- matrix(rep(NA, length(eval_objects)*N), nrow=N); colnames(eval_centers) <- eval_objects
  #evaluate the distance between origin data and the different data with imputed values
  eval_imp_values <- matrix(rep(NA, length(eval_objects)*N), nrow=N); colnames(eval_imp_values) <- eval_objects
  #evaluate the cluster partitions of origin data and with imputed values via rand index
  eval_rand <- matrix(rep(NA, length(eval_objects)*N), nrow=N); colnames(eval_rand) <- eval_objects
  
  duration <- NULL
  
  for(i in 1:N){
    
    ### reproducibility
    set.seed(i)
    
    ### simulated data
    m <- which(stringr::str_remove_all(paste0(trial_design[n,-c(1,5)], collapse="_"),"\\.") == td_unique)
    origin <- dat_complete[[m]][[i]][[1]]

    ### simulated missings
    data_NA <- dat_incomplete[[n]][[i]]
    
    ### Imputation and cluster
    # ignore missings
    time_1 <- Sys.time()
    kpres_ignore <- kproto_kPOD(x = data_NA, k = dat_complete[[m]][[i]][[2]], lambda = dat_complete[[m]][[i]][[3]], nstart = 3, verbose = FALSE, na.rm = "no_imp")
    duration_time <- Sys.time() - time_1
    
    # impute with via k-POD (intern)
    time_1 <- Sys.time()
    kpres_kpod_intern <- kproto_kPOD(x = data_NA, k = dat_complete[[m]][[i]][[2]], lambda = dat_complete[[m]][[i]][[3]], nstart = 3, verbose = FALSE, na.rm = "impute")
    duration_time <- c(duration_time, Sys.time() - time_1)
    
    # impute with via k-POD (extern)
    time_1 <- Sys.time()
    kpres_kpod_extern <- kproto_kPOD_ex(x = data_NA, k = dat_complete[[m]][[i]][[2]], lambda = dat_complete[[m]][[i]][[3]], nstart = 3, iters = 100)
    duration_time <- c(duration_time, Sys.time() - time_1)
    
    # impute via MICE, then cluster with kproto and pooling via Gionis et al.
    time_1 <- Sys.time()
    kpm_try <- try(suppressWarnings(kproto_mice(data_NA = data_NA, n_imp = 5,  k = dat_complete[[m]][[i]][[2]], lambda = dat_complete[[m]][[i]][[3]], nstart = 3, iters = 100)),
                   silent = TRUE)
    duration_time <- c(duration_time, Sys.time() - time_1)
    if(inherits(kpm_try, "try-error")){
      kpres_mice <- NA
    }else{
      kpres_mice <- kpm_try
    }
    
    duration <- rbind(duration, duration_time)
    
    ### evaluate
    # - (Goal: Minimization) total sum of within cluster distances (distance origin data and kpres_XXX$centers)
    # - (Goal: Minimization) distance between cluster centers of origin data and calculated centers
    # - (Goal: Minimization) distance between origin data and data with imputed values
    # - (Goal: Maximization) rand index (\in [0,1], where 1 means the two clustering outcomes match identicaly)
    
    for(v in 1:length(eval_objects)){
      kpres_out <- switch(eval_objects[v],
                          "ignore" = kpres_ignore, "kpod_in" = kpres_kpod_intern, "kpod_ex" = kpres_kpod_extern, 
                          "mice_pooling" = kpres_mice)
      if(!is.na(unlist(kpres_out)[[1]])){
        eval_total_within_dists[i,v] <- eval_within_dist(origin = origin[,-c(1,2)], kpres = kpres_out, 
                                                         lambda = dat_complete[[m]][[i]][[3]], k_opt = dat_complete[[m]][[i]][[2]])
        eval_centers[i,v] <- eval_dist_protos(protos = kpres_out$centers, protos_origin = dat_complete[[m]][[i]][[4]], 
                                              lambda = dat_complete[[m]][[i]][[3]], k_opt = dat_complete[[m]][[i]][[2]])
        eval_imp_values[i,v] <- eval_dist_iv(origin = origin[,-c(1,2)], x_iv = kpres_out$data, lambda = dat_complete[[m]][[i]][[3]])
        eval_rand[i,v] <- fossil::adj.rand.index(kpres_out$cluster, origin$kpres)
      }
    }
    
    
  }
  
  eval_twd <- apply(eval_total_within_dists, 1, rank); eval_twd <- apply(eval_twd, 1, mean)
  eval_c <- apply(eval_centers, 1, rank); eval_c <- apply(eval_c, 1, mean)
  eval_iv <- apply(eval_imp_values, 1, rank); eval_iv <- apply(eval_iv, 1, mean)
  eval_r <- 5 - apply(eval_rand, 1, rank); eval_r <- apply(eval_r, 1, mean)
  
  results_eval <- c(list(list(trial_design[n,], duration, rbind(eval_twd, eval_c, eval_iv,eval_r),
                              eval_total_within_dists, eval_centers, eval_imp_values, eval_rand)))

  return(results_eval)
}

# set the number of cores to be used
n_cores <- 4
parallelStartSocket(n_cores, logging = TRUE, load.balancing = TRUE)
parallelLibrary("clustMixType", "ggplot2", "tibble", "mice", "rpart")
parallelSource("functions_clustering.R", "functions_evaluation.R")
parallelExport("trial_design", "dat_complete", "dat_incomplete", "td_unique")
output <- parallelLapply(1:nrow(trial_design), fun = par_fun)
parallelStop()

results_all <- lapply(output, function(x) x[[1]])


