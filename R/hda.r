#############################
### Definition of method hda
hda <- function (x, ...) 
    UseMethod("hda")

################################################
### formuala interface for hda
hda.formula <- function(formula, data = NULL, ...)
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- hda(x, grouping, ...)
    res$terms <- Terms
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

##############################################
### print function
print.hda <- function(x, ...){
  cat("\nHeteroscedastic discriminant analysis\n\n")
  cat("Call:\n")
  print(x$hda.call)

  cat("\nDimension of reduced space: \n", x$reduced.dimension, "\n\n")
  
  cat("\nClass means in reduced space: \n")
  print(x$class.dist$new.classmeans)
  cat("\n\n")

  cat("\nClass covariances in reduced space: \n")
  print(x$class.dist$new.classcovs)
  cat("\n\n")
  
  invisible(x)
  }

######################
### visualize loadings
showloadings <- function(object, comps = 1:object$reduced.dimension, ...){
  if (max(comps) > nrow(object$hda.loadings)) 
    stop("Component ids have to be <= dimension of object$hda.loadings")
  vnames <- rownames(object$hda.loadings)
  d <- length(comps)
  if(d > 2){
    op <- par(mfrow=c(d,d))
    for(i in comps){
        for(j in comps){
            plot(object$hda.loadings[,c(j,i)], type="n", ...)
            for(k in 1:nrow(object$hda.loadings)) 
                text(object$hda.loadings[k,j], object$hda.loadings[k,i], vnames[k]) 
        }
    }
  }
  if(d==2){
    op <- par(mfrow=c(1,1))
    plot(object$hda.loadings[,c(comps[1],comps[2])], type="n", ...)
    for(k in 1:nrow(object$hda.loadings)) 
        text(object$hda.loadings[k,comps[1]], object$hda.loadings[k,comps[2]], vnames[k]) 
  }    
  if(d==1){
    op <- par(mfrow=c(1,1))
    plot(object$hda.loadings[,comps], type="n", xlab="Variable index", ylab=colnames(object$hda.loadings)[comps], ...)
    for(k in 1:nrow(object$hda.loadings)) 
        text(k, object$hda.loadings[k,comps], vnames[k]) 
  }     
  par(op)
}


######################
### plot function to visualize score
plot.hda <- function(x, comps = 1:x$reduced.dimension, col = x$grouping, ...){
  if (max(comps) > nrow(x$hda.loadings)) 
    stop("Component ids have to be <= dimension of object$hda.loadings!")
  if (length(comps) > 1) 
    plot(as.data.frame(x$hda.scores[,comps]), col = col, ...)
  if (length(comps) == 1) 
    plot(x$hda.scores[,comps], col = col, ylab = paste("comp",comps,sep=" "), ...)
}

###########################################################################
### predict function for easy transformation of new data with a given model
predict.hda <- function(object, newdata, alldims = FALSE, task = c("dr", "c"), ...){
  if (class(object) != "hda") 
    stop("Object must be of class hda!")
  task <- match.arg(task)
  if(is.data.frame(newdata)) 
    newdata <- as.matrix(newdata)
  if(!is.matrix(newdata)) 
    stop("Newdata must be of type matrix or data frame!")
  if(dim(object$hda.loadings)[2] != dim(newdata)[2]) 
    stop("Newdata must be of same dimension as the explanatory input variables data set!")
  if(task == "dr"){
    new.transformed.data <- newdata %*% object$hda.loadings
    if(!alldims) 
        new.transformed.data <- new.transformed.data[,1:object$reduced.dimension]
    return(new.transformed.data)
  }
  if (task == "c"){
    if(class(object$naivebayes) != "naiveBayes") 
        stop("Classification of newdata can only be done id option 'crule = TRUE' has been chosen at hda() call")
    new.transformed.data <- newdata %*% object$hda.loadings
    new.transformed.data <- new.transformed.data[,1:object$reduced.dimension]
    if (object$reduced.dimension > 1) classification.result <- predict(object$naivebayes, new.transformed.data, ...)
    if (object$reduced.dimension == 1){
      warning("Function predict.naiveBayes {e1071} not implemented for dimension 1 of the reduced space!\n
               Prediction is done without using naiveBayes.\n")
    
      priors <- object$naivebayes$apriori/ sum(object$naivebayes$apriori)  
      classids <- 1:nrow(object$naivebayes$tables[[1]])  ### muss hier anstatt x der variablenname hin?
    
      priodens <- function(i) return(priors[i]*dnorm(new.transformed.data, mean = object$naivebayes$tables[[1]][i,1], sd = object$naivebayes$tables[[1]][i,2]))
      posteriors <- sapply(classids, priodens)
      posteriors <- t(apply(posteriors, 1, function(x) return(x / sum(x))))
      classification.result <- rownames(object$naivebayes$tables[[1]])[apply(posteriors, 1, which.max)]
    }
    
    result <- list(classification.result, new.transformed.data)
    names(result) <- c("classification.result", "new.transformed.data")
    return(result)
    }                     
  }

##############################################################################
### calculation of heteroscedastic discriminant analysis loadings matrix 
compute.loadings <- function(covarray, clsizes, newd, initial = NULL, iters = 7, initers = 10, verb = TRUE){
    
    # covarray: array (of dimension oldd, oldd, classnumber) of the estimated class-specific covariance matrices (of size (oldd, oldd)
    # clsizes : number of corresponding observations per class
    # newd    : dimension of the desired subspace
    # iters   : number of optimization iterations
    # initers : iterations inner loop
    
    if (dim(covarray)[1] != dim(covarray)[2]) stop("Covariance matrices must be quadratic")
    if (newd >= dim(covarray)[2]) stop("Dimension of reduced space should be lower than original")
 
    ncl             <- dim(covarray)[3]
    oldd            <- dim(covarray)[1]
    totN            <- sum(clsizes)  # Total number of observations
    commonwithin    <- matrix(rowSums(sapply(1:ncl, function(i) clsizes[i] * covarray[,,i])), ncol=oldd, nrow=oldd) / totN

    # initialization loadings matrix
    if(!is.null(initial)){
      if (sum(dim(initial)== (dim(covarray)[1:2])) != 2) 
        stop("(Optional) initalization of loading matrix must be quadratic and of size dim(x)[2]")
      Trafo           <- initial  
      }    
    if (is.null(initial)) Trafo   <- diag(oldd)  
    invG            <- array(dim=c(oldd,oldd,oldd))

    for (iter in 1:iters){
      Qu  <- totN * log(det(Trafo)^2) - totN * oldd * (log(2*pi)+1)
      for (i in 1:oldd){
          if (i <= newd){
              G <- matrix(0,oldd,oldd)
              for (m in 1:ncl){
                  sigma_i <- as.numeric(t(Trafo[,i]) %*% covarray[,,m] %*% Trafo[,i])
                  G       <- G + clsizes[m] / sigma_i * covarray[,,m] 
                  Qu       <- Qu - clsizes[m] * log(sigma_i)
                  }
              }
          if (i > newd){
              sigma_i <- as.numeric(t(Trafo[,i]) %*% commonwithin %*% Trafo[,i])
              G       <- totN / sigma_i * commonwithin
              Qu       <- Qu - totN * log(sigma_i)
              }
           
          invG[,,i] <- solve(G)
          }
          
        Qu <- Qu/2
        if (verb) cat(paste("Iteration", iter,"Log Likelihood: ",  Qu, "\n"))
        for (initer in 1:initers){
            for (i in 1:oldd){
                Ce          <- t(solve(t(Trafo)) *  det(t(Trafo)))
                ci_invG     <- Ce[i,] %*% invG[,,i]
                Trafo[,i]  <- t(ci_invG * sqrt(totN / as.numeric((ci_invG %*% t(t(Ce[i,]))))))
                }
            }
        }
     
    Trafo = Trafo / (det(Trafo)^(1/oldd))
    return(Trafo)
    } 

#####################################################################################
### regularize covariance matrices according to input parameters as in Friedman (1989)  
regularize <- function(covarray, clsizes, lambd, gamm){

  if(lambd < 0 || lambd > 1) 
    stop("lambd and gamm must be in [0,1]")
  if(gamm < 0 || gamm > 1) 
    stop("lambd and gamm must be in [0,1]")

  # calculation of common covaiance matrix (over al classes)
  commoncov <- (clsizes[1]-1) * covarray[,,1]
  for(i in 2:dim(covarray)[3]) 
    commoncov <- commoncov + (clsizes[i]-1) * covarray[,,i]
  commoncov <- commoncov / (sum(clsizes)-dim(covarray)[3])
  
  # regularization of covaiances - array
  for(i in 1:dim(covarray)[3]){
    # towards equal covariances
    covarray[,,i] <- (lambd * clsizes[i] * covarray[,,i] + (1-lambd) * sum(clsizes) * commoncov) / 
                     (lambd * clsizes[i] + (1-lambd) * sum(clsizes)) 
    # towards diagonality with average variance
    average.variance.i  <- sum(diag(covarray[,,i]))/dim(covarray)[1]
    shrinkcov.i         <- diag(average.variance.i, dim(covarray)[1])
    covarray[,,i]       <- gamm * (covarray[,,i]) + (1-gamm) * shrinkcov.i     
    }  

  return(covarray)
  }

######################################################  
### main call of heteroscedastic discriminant analysis 
# ... calculation of the input parameters from data set and call of compute.loadings        
hda.default <- function(x, grouping, newdim = 1:(ncol(x)-1), crule = FALSE, 
    reg.lamb = NULL, reg.gamm = NULL, initial.loadings = NULL, 
    sig.levs = c(0.05,0.05), noutit = 7, ninit = 10, verbose = TRUE, ...){
  # x:                  data frame containing the (metric) variables  
  # grouping:           vector of class labels
  # newdim:             dimension of reduced space
  # crule:              TRUE if naiveBayes() model shoulf be built using pkg e1071
  # reg.lamb,reg.gamm:  regularization parameters as in rda() [Friedman, 1989]
  # initial.laodings:   (optional) initial guess of the (quadratic) loadings matrix 
  # sig.levs:           significance levels for tests if no unique newdim is specified
  # noutit:             (optional) number of outer iterations of the algorithm
  # ninit:              (optional) number of inner iterations of the algorithm

  # check for possible errors in function call
  if(length(table(grouping)) < 2) 
    stop("Class vector should contain different levels")
  if(dim(x)[2] < 2) 
    stop("Dimensionality reduction only meaningful if x has at least two coloumns")

  # class parameters
  clsizes <- table(grouping)
  clnb    <- length(clsizes)
  # array of covariance matrices:
  covlist <- by(as.data.frame(x), grouping, cov)
  # dimension of the original data space
  dms     <- dim(covlist[[1]])[1]
  carray  <- array(0,c(dms,dms,clnb))
  for (i in seq_along(covlist)){
    carray[,,i] <- covlist[[i]]
  }
  # regularization ov covariance estimates if specified by input parameters  
  if (sum(c(is.null(reg.lamb),is.null(reg.gamm))) < 2) 
    carray <- regularize(carray, clsizes, reg.lamb, reg.gamm)

  hda.loadings <- NULL
  trace.newdim <- NULL
  # recursive call of hda.default if not a single newdim as specified
  if(length(newdim) > 1){
    pvaleqmeans     <- 0
    pvalhomogcovs   <- 0
    d               <- 1 
    while(((d <= length(newdim)) * ((pvaleqmeans <= sig.levs[1]) + (pvalhomogcovs <= sig.levs[2]))) > 0){
      if (verbose) 
        cat("\nnewdim = ", newdim[d], "\n")
      dummy <- hda(x=x, grouping=grouping, newdim=newdim[d], reg.lamb = reg.lamb, 
                   reg.gamm = reg.gamm, initial.loadings = initial.loadings, 
                   sig.levs = sig.levs, noutit = noutit, ninit = ninit, verbose = verbose)
      if((ncol(x)-newdim[d]) > 1) 
        pvaleqmeans <- dummy$eqmean.test[[3]][6]
      if((ncol(x)-newdim[d]) == 1) 
        pvaleqmeans <- dummy$eqmean.test[[3]][1,5]
      pvalhomogcovs <- dummy$homog.test$pValue
      trace.newdim  <- cbind(trace.newdim, c(pvaleqmeans, pvalhomogcovs)) 
      rownames(trace.newdim) <- c("eqmean.tests","homog.tests")
      colnames(trace.newdim) <- newdim[1:d]
      d <- d+1
    }
    if (verbose) cat("\n")
    hda.loadings  <- dummy$hda.loadings
    newdim    <- newdim[d-1]
  }           


  # calculate transformation matrix if not already done (i.e. if 
  if (is.null(hda.loadings)){
    if(round(newdim) != newdim) 
        stop("newdim must be an integer")
    hda.loadings <- compute.loadings(carray, clsizes, newdim, 
                                     initial = initial.loadings, iters = noutit, 
                                     initers = ninit, verb = verbose)
  }
  rownames(hda.loadings) <- colnames(as.data.frame(x))
  compnames  <- NULL
  for(i in 1:dms) compnames <- c(compnames, paste("comp",i,sep=""))
  colnames(hda.loadings) <- compnames 
  
  # calculate the transformed space
  newspace  <- as.data.frame(as.matrix(x) %*% hda.loadings)
  
  # class distribution parameters in new space
  new.classmeans  <- by(as.data.frame(newspace[,1:newdim]),grouping,colMeans)
  new.classcovs   <- by(as.data.frame(newspace[,1:newdim]),grouping,cov)
  new.classdist   <- list(new.classmeans = new.classmeans, new.classcovs = new.classcovs)
  
  ### several tests of assumptions  
  # test on homogenity of the classwise covariance matrices in remaining dimensions 
  newcovs     <- by(newspace, grouping, cov) # ml-estimates of classes in redundant dimensions
  newcommon   <- matrix(0,dms-newdim,dms-newdim) # initialize common covariance matrix in redundant transformed space
  stat        <- 0 # initialize value of statistic
  for(i in seq_along(clsizes)){ 
    newcovs[[i]]  <- (newcovs[[i]])[(newdim+1):dms,(newdim+1):dms] * (clsizes[i]-1)/clsizes[i]
    stat          <- stat - clsizes[i] * log(det(as.matrix(newcovs[[i]])))
    newcommon     <- newcommon + newcovs[[i]] * clsizes[i] 
  }
  newcommon       <- newcommon / sum(clsizes)
  stat           <- as.numeric(stat + sum(clsizes) * log(det(as.matrix(newcommon))))
  dfs            <- (dms-newdim)*(dms-newdim+1)*(length(clsizes)-1)/2
  pval           <- 1 - pchisq(stat, dfs) # corresponding p-value
  homog.test     <- list(new.common.covariance = newcommon, 
                         new.class.covariances = newcovs, statistic = stat, 
                         dfs = dfs, pValue = pval) 
  
  # test on equal class means in remaining dimension
  new.classmeans2  <- by(as.data.frame(newspace[,(newdim+1):dms]),grouping,colMeans)
  new.classcovs2   <- by(as.data.frame(newspace[,(newdim+1):dms]),grouping,cov)
  if((dms-newdim) > 1) 
    stats <- summary(manova(as.matrix(newspace[,(newdim+1):dms])~grouping), test="Wilks")[[4]][1,]
  if((dms-newdim) == 1) 
    stats <- summary(aov(as.matrix(newspace[,(newdim+1):dms])~grouping))[[1]]
  eqmean.test <- list(new.class.means=new.classmeans2, new.class.covariances=new.classcovs2, stats)
  
  cl <- match.call()
  cl[[1]] <- as.name("hda")
  
  if (crule){
    crule <- naiveBayes(newspace[,1:newdim], grouping)
  } 
  result <- list(hda.loadings = hda.loadings, hda.scores = newspace, 
                 grouping = grouping, class.dist = new.classdist, 
                 reduced.dimension = newdim, naivebayes = crule, 
                 reg.lamb = reg.lamb, reg.gamm = reg.gamm, eqmean.test = eqmean.test, 
                 homog.test = homog.test, hda.call = cl, trace.dimensions = trace.newdim)
                 
  class(result) <- "hda"
  return(result) 
}
