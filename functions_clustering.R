

kproto_kPOD <- function(x, k, lambda = NULL, iter.max = 100, nstart = 1, na.rm = "yes", keep.data = TRUE, verbose = TRUE, ...){
  
  # enable input of tibbles
  if(is_tibble(x) == TRUE){x <- as.data.frame(x)}
  
  # initial error checks
  if(!is.data.frame(x)) stop("x should be a data frame!")
  if(ncol(x) < 2) stop("For clustering x should contain at least two variables!")
  if(iter.max < 1 | nstart < 1) stop("iter.max and nstart must not be specified < 1!")
  if(!is.null(lambda)){
    if(any(lambda < 0)) stop("lambda must be specified >= 0!")
    if(!any(lambda > 0)) stop("lambda must be specified > 0 for at least one variable!")
  }
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  if(!anynum) stop("\n No numeric variables in x! Try using kmodes() from package klaR...\n\n")
  if(!anyfact) stop("\n No factor variables in x! Try using kmeans()...\n\n")
  
  # treatment of missings
  NAcount <- apply(x, 2, function(z) sum(is.na(z)))
  if(verbose){
    cat("# NAs in variables:\n")
    print(NAcount)
  }
  if(any(NAcount == nrow(x))) stop(paste("Variable(s) have only NAs please remove them:",names(NAcount)[NAcount == nrow(x)],"!"))
  if(na.rm == "yes") {
    miss <- apply(x, 1, function(z) any(is.na(z)))
    if(verbose){
      cat(sum(miss), "observation(s) with NAs.\n")
      if(sum(miss) > 0) message("Observations with NAs are removed.\n")
      cat("\n")
    } 
    x <- x[!miss,]
  } # remove missings
  
  if(na.rm != "yes"){
    allNAs <- apply(x,1,function(z) all(is.na(z)))
    if(sum(allNAs) > 0){
      if(verbose) cat(sum(allNAs), "observation(s) where all variables NA.\n")
      warning("No meaningful cluster assignment possible for observations where all variables NA.\n")
      if(verbose) cat("\n")
      
    }
  }
  
  if(na.rm == "imp_internal"){origin <- x}
  
  if(nrow(x) == 1) stop("Only one observation clustering not meaningful.")
  
  k_input <- k # store input k for nstart > 1 as clusters can be merged 
  
  # initialize prototypes
  if(!is.data.frame(k)){
    if (length(k) == 1){
      if(as.integer(k) != k){k <- as.integer(k); warning(paste("k has been set to", k,"!"))}
      if(sum(complete.cases(x)) < k) stop("Data frame has less complete observations than clusters!")
      ids <- sample(row.names(x[complete.cases(x),]), k)
      protos <- x[ids,]
    }
    if (length(k) > 1){
      if(nrow(x) < length(k)) stop("Data frame has less observations than clusters!")
      ids <- k
      k <- length(ids)
      if(length(unique(ids)) != length(ids)) stop("If k is specified as a vector it should contain different indices!")
      if(any(ids<1)|any(ids>nrow(x))) stop("If k is specified as a vector all elements must be valid indices of x!")
      #check for integer
      protos <- x[ids,]
      if(any(!complete.cases(protos))) stop("Choose initial prototypes without missing values!")
    }
    rm(ids)
  }
  if(is.data.frame(k)){
    if(nrow(x) < nrow(k)) stop("Data frame has less observations than clusters!")
    if(length(names(k)) != length(names(x))) stop("k and x have different numbers of columns!")
    if(any(names(k) != names(x))) stop("k and x have different column names!")
    if(anynum) {if( any(sapply(k, is.numeric) != numvars)) stop("Numeric variables of k and x do not match!")}
    if(anyfact) {if( any(sapply(k, is.factor) != catvars)) stop("Factor variables of k and x do not match!")}
    protos <- k
    if(any(!complete.cases(protos))) stop("Prototypes with missing values. Choose initial prototypes without missing values!")
    k <- nrow(protos)
  }
  if(k < 1) stop("Number of clusters k must not be smaller than 1!")
  
  # automatic calculation of lambda
  if(length(lambda) > 1) {if(length(lambda) != sum(c(numvars,catvars))) stop("If lambda is a vector, its length should be the sum of numeric and factor variables in the data frame!")}
  if(is.null(lambda)){
    if(anynum & anyfact){
      vnum <- mean(sapply(x[,numvars, drop = FALSE], var, na.rm = TRUE))
      vcat <- mean(sapply(x[,catvars, drop = FALSE], function(z) return(1-sum((table(z)/sum(!is.na(z)))^2))))
      if (vnum == 0){
        if(verbose) warning("All numerical variables have zero variance.")
        anynum <- FALSE
      } 
      if (vcat == 0){
        if(verbose) warning("All categorical variables have zero variance.")
        anyfact <- FALSE
      } 
      if(anynum & anyfact){
        lambda <- vnum/vcat
        if(verbose) cat("Estimated lambda:", lambda, "\n\n")
      }else{
        lambda <- 1
      }
    }
  }
  
  if(na.rm == "yes"){
    miss <- apply(x, 1, function(z) any(is.na(z)))
    if(verbose){
      cat(sum(miss), "observation(s) with NAs.\n")
      if(sum(miss) > 0) message("Observations with NAs are removed.\n")
      cat("\n")
    }
    x <- x[!miss,]
  }
  
  # initialize clusters
  clusters  <- numeric(nrow(x)) 
  tot.dists <- NULL
  moved   <- NULL
  iter <- 1
  
  # check for any equal prototypes and reduce cluster number in case of occurence
  if(k > 1){
    keep.protos <- rep(TRUE,k)
    for(l in 1:(k-1)){
      for(m in (l+1):k){
        d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
        d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # wtd simple matching for categorics 
        if((d1+d2) == 0) keep.protos[m] <- FALSE 
      }
    }
    if(!all(keep.protos)){
      protos <- protos[keep.protos,]
      k <- sum(keep.protos)
      if(verbose) message("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")      
    }
  }
  
  # special case only one cluster
  if(k == 1){clusters <- rep(1, nrow(x)); size  <- table(clusters); iter <- iter.max} # REM: named vector size is needed later...
  
  # start iterations for standard case (i.e. k > 1)
  while(iter < iter.max){
    
    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
      if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
      d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
      dists[,i] <- d1 + d2
    }
    
    # assign clusters 
    old.clusters  <- clusters
    clusters      <- apply(dists, 1, function(z) {a <- which(z == min(z)); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
    size          <- table(clusters)  
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
    
    # ...check for empty clusters and eventually reduce number of prototypes    
    if (length(size) < k){
      k <- length(size)
      protos <- protos[1:length(size),]  
      if(verbose) cat("Empty clusters occur. Cluster number reduced to:", k, "\n\n")
    }
    
    # trace
    tot.dists <- c(tot.dists, sum(tot.within))      
    moved <- c(moved, sum(clusters != old.clusters))
    
    # compute new prototypes
    remids <- as.integer(names(size))
    for(i in remids){
      some_vals <- sapply(x[clusters == i, , drop = FALSE], function(z) !all(is.na(z))) # only update variables if not all values are NA
      if(any(some_vals & numvars)){
        protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], mean, na.rm = TRUE)
      }
      if(any(some_vals & catvars)){
        protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      }
    }
    
    # update missing values, if na.rm = "imp_internal"
    if(na.rm == "imp_internal"){
      x_old <- x
      x <- origin
      
      if(any(is.na(x[,numvars]))){
        x[,numvars][is.na(x[,numvars])] <- protos[clusters,numvars][is.na(x[,numvars])]
      }
      if(any(is.na(x[,catvars]))){
        x[,catvars][is.na(x[,catvars])] <- protos[clusters,catvars][is.na(x[,catvars])]
      }
      
      # observation with all values NA won't be imputed
      x[allNAs,] <- NA
    }
    
    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}
    
    # check for any equal prototypes and reduce cluster number in case of occurence
    if(iter == (iter.max-1)){ # REM: for last iteration equal prototypes are allowed. otherwise less prototypes than assigned clusters.
      keep.protos <- rep(TRUE,k)
      for(l in 1:(k-1)){
        for(m in (l+1):k){
          d1 <- sum((protos[l,numvars, drop = FALSE]-protos[m,numvars, drop = FALSE])^2) # euclidean for numerics
          d2 <- sum(protos[l,catvars, drop = FALSE] != protos[m,catvars, drop = FALSE]) # wtd simple matching for categorics 
          if((d1+d2) == 0) keep.protos[m] <- FALSE 
        }
      }
      if(!all(keep.protos)){
        protos <- protos[keep.protos,]
        k <- sum(keep.protos)
        if(verbose) cat("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")      
      }
    }
    
    # add stopping rules depending on imputed values
    if(moved[length(moved)] ==  0){
      if(na.rm != "imp_internal"){
        break
      }else{
        if(sum(!(x == x_old), na.rm = TRUE) == 0){
          break
        }
      }
    }
    
    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}
    
    iter <- iter+1
  }
  
  
  # Final update of prototypes and dists
  if(iter == iter.max){ # otherwise there have been no moves anymore and prototypes correspond to cluster assignments 
    # compute new prototypes
    remids <- as.integer(names(size))
    for(i in remids){
      some_vals <- sapply(x[clusters == i, , drop = FALSE], function(z) !all(is.na(z))) # only update variables if not all values are NA
      if(any(some_vals & numvars)){
        protos[which(remids == i), some_vals & numvars] <- sapply(x[clusters == i, some_vals & numvars, drop = FALSE], mean, na.rm = TRUE)
      }
      if(any(some_vals & catvars)){
        protos[which(remids == i), some_vals & catvars] <- sapply(x[clusters == i, some_vals & catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
      }
    }
    
    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      d1 <- (x[,numvars, drop = FALSE] - matrix(rep(as.numeric(protos[i, numvars, drop = FALSE]), nrows), nrow=nrows, byrow=T))^2
      if(length(lambda) == 1) d1 <- rowSums(d1, na.rm = TRUE)
      if(length(lambda) > 1) d1 <- as.matrix(d1) %*% lambda[numvars]
      d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
      d2[is.na(d2)] <- FALSE
      if(length(lambda) == 1) d2 <- lambda * rowSums(d2)
      if(length(lambda) > 1) d2 <- as.matrix(d2) %*% lambda[catvars]
      dists[,i] <- d1 + d2
    }
    
    size          <- table(clusters)  
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
  }
  
  # observations with all NA are not assigned to a cluster  
  if(na.rm != "yes"){
    if(sum(allNAs) > 0){
      clusters[allNAs] <- NA
      dists[allNAs,] <- NA
    }
  }
  
  names(clusters) <- row.names(dists) <- row.names(x)
  rownames(protos) <- NULL
  # create result: 
  res <- list(cluster = clusters,  
              centers = protos, 
              lambda = lambda, 
              size = size,
              withinss = within,
              tot.withinss = tot.within,   
              dists = dists, 
              iter = iter, 
              trace = list(tot.dists = tot.dists, moved = moved))
  
  if(keep.data) res$data = x
  
  # loop: if nstart > 1:
  if(nstart > 1)
    for(j in 2:nstart){
      # start function kproto_kPOD with the origin data in case of imputated data
      if(na.rm == "imp_internal"){
        x <- origin
      }
      res.new <- kproto_kPOD(x=x, k=k_input, lambda = lambda,  iter.max = iter.max, nstart=1, verbose=verbose, na.rm = na.rm)
      if(res.new$tot.withinss < res$tot.withinss) res <- res.new
    }  
  
  if(na.rm == "imp_onestep"){
    x <- res$data
    if(any(is.na(x[,numvars]))){
      x[,numvars][is.na(x[,numvars])] <- protos[clusters,numvars][is.na(x[,numvars])]
    }
    if(any(is.na(x[,catvars]))){
      x[,catvars][is.na(x[,catvars])] <- protos[clusters,catvars][is.na(x[,catvars])]
    }
    res$data <- x
  }
  
  class(res) <- "kproto"
  return(res)
}



#########################################################################################

kproto_kPOD_ex <- function(x, k, lambda = NULL, nstart = 1, iters = 100){
  x_iv <- x
  numvars <- sapply(x, is.numeric)
  catvars <- sapply(x, is.factor)
  
  iter <- 1
  while(iter < iters){
    res_no <- kproto_kPOD(x = x_iv, k = k, lambda = lambda, nstart = nstart, na.rm = "no", verbose = FALSE, iter.max = 100)
    
    # update imputed values
    x_old <- x_iv
    x_iv <- x
    if(any(is.na(x[,numvars]))){
      x_iv[,numvars][is.na(x_iv[,numvars])] <- res_no$centers[res_no$cluster,numvars][is.na(x_iv[,numvars])]
    }
    if(any(is.na(x[,catvars]))){
      x_iv[,catvars][is.na(x_iv[,catvars])] <- res_no$centers[res_no$cluster,catvars][is.na(x_iv[,catvars])]
    }
    
    if(iter > 1){
      if(sum(!(x_iv == x_old), na.rm = TRUE) == 0){
        break
      }
    }else{
      comp <- rowSums(!(x_iv == x_old))
      comp[!is.na(comp)] <- FALSE; comp[is.na(comp)] <- TRUE
      if(sum(comp) == 0) break
    }
    
    iter <- iter + 1
  }
  
  res_no$data <- x_iv
  return(res_no)
}



#########################################################################################

kproto_mice <- function(data_NA, n_imp,  k, lambda, nstart, iters){
  
  # MICE und je imputed dataset excecute kproto
  mice_model <- mice(data = data_NA, m = n_imp, print = FALSE)
  cluster_mice <- NULL; data_all_imp <- NULL
  for(n in 1:n_imp){
    data_mice <- complete(data = mice_model, action = n)
    data_all_imp <- rbind(data_all_imp, data_mice)
    
    kpres_mice <- kproto(x = data_mice, k = k, lambda = lambda, iters = iters, nstart = nstart, verbose = FALSE, na.rm = FALSE)
    cluster_mice <- cbind(cluster_mice, kpres_mice$cluster)
  }
  
  #pooling via Gionis et al "Clustering Aggregation"
  d_fct <- function(b1, b2){sum(b1 != b2)/length(b1)}
  dists_myd <- proxy::dist(cluster_mice, method = d_fct)
  x <- hclust(d = dists_myd, method = "average")
  kpres_pooled <- cutree(x, k = k)
  
  #calculate prototypes of pooled cluster
  clusters <- rep(kpres_pooled, n_imp)
  
  #compute new prototypes
  protos <- data_all_imp[unique(clusters),]
  numvars <- sapply(data_all_imp, is.numeric); catvars <- sapply(data_all_imp, is.factor)
  for(i in unique(clusters)){
    protos[i, numvars] <- sapply(data_all_imp[clusters == i, numvars, drop = FALSE], mean, na.rm = TRUE)
    protos[i, catvars] <- sapply(data_all_imp[clusters == i, catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
  }
  
  #impute values
  data_imp <- data_NA
  if(any(is.na(data_NA[,numvars]))){
    data_imp[,numvars][is.na(data_imp[,numvars])] <- protos[kpres_pooled, numvars][is.na(data_imp[,numvars])]
  }
  if(any(is.na(data_NA[,catvars]))){
    data_imp[,catvars][is.na(data_imp[,catvars])] <- protos[kpres_pooled, catvars][is.na(data_imp[,catvars])]
  }
  
  protos_2 <- data_imp[unique(kpres_pooled),]
  for(i in unique(kpres_pooled)){
    protos_2[i, numvars] <- sapply(data_imp[kpres_pooled == i, numvars, drop = FALSE], mean, na.rm = TRUE)
    protos_2[i, catvars] <- sapply(data_imp[kpres_pooled == i, catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
  }
  
  return(list("centers" = protos_2,
              "cluster" = kpres_pooled, 
              "data" = data_imp))
}