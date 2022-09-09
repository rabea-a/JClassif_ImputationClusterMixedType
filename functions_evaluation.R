
#########################################################################################

# calculate "total sum of within cluster distance" of origin data with different cluster partitions
eval_within_dist <- function(origin, kpres, lambda, k_opt){
  numvars <- sapply(origin, is.numeric)
  catvars <- sapply(origin, is.factor)
  
  nrows <- nrow(origin)
  dists <- matrix(NA, nrow = nrows, ncol = k_opt)
  for(i in 1:k_opt){
    d1 <- (origin[, numvars, drop = FALSE] - matrix(rep(as.numeric(kpres$centers[i, numvars, drop = FALSE]), nrows), nrow = nrows, byrow=T))^2
    d1 <- rowSums(d1, na.rm = TRUE)
    
    d2 <- sapply(which(catvars), function(j) return(origin[,j] != rep(kpres$centers[i,j], nrows)))
    d2[is.na(d2)] <- FALSE
    d2 <- lambda * rowSums(d2)
    
    dists[,i] <- d1 + d2
  }
  
  min.dists     <- apply(cbind(kpres$cluster, dists), 1, function(z) z[z[1]+1])
  within        <- as.numeric(by(min.dists, kpres$cluster, sum))
  return(sum(within))
}



#########################################################################################

# calculate the distance between kpres_origin$centers and centers
eval_dist_protos <- function(protos, protos_origin, lambda, k_opt){
  numvars <- sapply(protos, is.numeric)
  catvars <- sapply(protos, is.factor)
  
  nrows <- nrow(protos_origin)
  dists <- matrix(NA, nrow = nrows, ncol = k_opt)
  for(i in 1:k_opt){
    d1 <- (protos_origin[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow = nrows, byrow=T))^2
    d1 <- rowSums(d1, na.rm = TRUE)
    
    d2 <- sapply(which(catvars), function(j) return(protos_origin[,j] != rep(protos[i,j], nrows)))
    d2[is.na(d2)] <- FALSE
    d2 <- lambda * rowSums(d2)
    
    dists[,i] <- d1 + d2
  }
  
  if(nrow(protos_origin) == 4){
    sum_dists <- apply(gtools::permutations(4,4), 1, function(z) sum(dists[1,z[1]], dists[2,z[2]], dists[3,z[3]], dists[4,z[4]]))
  }else{
    sum_dists <- apply(gtools::permutations(2,2), 1, function(z) sum(dists[1,z[1]], dists[2,z[2]]))
  }
  
  return(min(sum_dists))
}



#########################################################################################

# calculate the distance between origin and x with imputed values
eval_dist_iv <- function(origin, x_iv, lambda){
  numvars <- sapply(origin, is.numeric)
  catvars <- sapply(origin, is.factor)
  
  d1 <- (origin[,numvars, drop = FALSE] - x_iv[, numvars, drop = FALSE])^2
  d1 <- rowSums(d1, na.rm = TRUE)
  
  d2 <- sapply(which(catvars), function(j) return(origin[,j] != x_iv[,j]))
  d2 <- lambda * rowSums(d2)
  
  dists <- sum(d1 + d2)
  
  return(dists)
}