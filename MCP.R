##########################
##### Main Function ######
##########################

# X: Exposure (Vector or Matrix)
# Y: Outcome (Vector or Matrix)
# M: Mediator (Matrix)
# COV.XM & COV.MY: Covariates (Matrix)


hima <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM, 
                 family = c("gaussian", "binomial"), 
                 penalty = c("MCP", "SCAD", "lasso"), 
                 topN = NULL, 
                 parallel = FALSE, 
                 ncore = 1, 
                 verbose = FALSE, 
                 ...) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  
  if (parallel & (ncore == 1)) ncore <- parallel::detectCores()
  
  n <- nrow(M)
  p <- ncol(M)
  
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
  
  d <- min(p, d) # if d > p select all mediators
  
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  #message("Step 1: Sure Independent Screening ...", "     (", Sys.time(), ")")
  
  if(family == "binomial")
  {
    # Screen M using X given the limited information provided by Y (binary)
    # Therefore the family is still gaussian
    if(verbose) message("    Screening M using the association between X and M: ", appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = "gaussian", modelstatement = "Mone ~ X", 
                                   parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else if (family == "gaussian"){
    # Screen M using Y (continuous)
    if(verbose) message("    Screening M using the association between M and Y: ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = family, modelstatement = "Y ~ Mone + X", 
                           parallel = parallel, ncore = ncore, verbose, tag = "Sure Independent Screening")
    SIS_Pvalue <- SIS_Results[2,]
  } else {
    stop(paste0("Family ", family, " is not supported."))
  }
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])  # the index of top mediators
  if(verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ","))
  
  M_SIS <- M[, ID]
  XM <- cbind(M_SIS, X)
  
  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  #message("Step 2: Penalized estimate (", penalty, ") ...", "     (", Sys.time(), ")")
  
  ## Based on the screening results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    fit <- ncvreg(XM, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y, family = family, 
                  penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV.MY))), ...)
  }
  # plot(fit)
  
  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) message("    Tuning parameter lambda selected: ", lam)
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  if(length(ID_1_non) == 0)
  {
    if(verbose) message("    All ", penalty, " beta estimates of the ", length(ID), " mediators are zero. Please check ncvreg package for more options (?ncvreg).")
  } else {
    if(verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
    beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
    ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M
    ## 
    
    if(family == "binomial")
    {
      ## This has been done in step 1 (when Y is binary)
      alpha <- alpha[,ID_test, drop = FALSE]
    } else {
      if(verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
      alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM, glm.family = "gaussian", 
                       modelstatement = "Mone ~ X", parallel = FALSE, ncore = ncore, 
                       verbose, tag = "site-by-site ordinary least squares estimation")
    }
    
    #########################################################################
    ################################ STEP 3 #################################
    #########################################################################
    #message("Step 3: Joint significance test ...", "     (", Sys.time(), ")")
    
    alpha_est_ID_test <- as.numeric(alpha[1, ])  #  the estimator for alpha
    P_adjust_alpha <- length(ID_test) * alpha[2, ]  # the adjusted p-value for alpha (bonferroni)
    P_adjust_alpha[P_adjust_alpha > 1] <- 1
    P_fdr_alpha <- p.adjust(alpha[2, ], "fdr")  # the adjusted p-value for alpha (FDR)
    
    alpha_est <- alpha_est_ID_test
    
    ## Post-test based on the oracle property of the MCP penalty
    if (is.null(COV.MY)) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
    }
    
    res <- summary(glm(Y ~ ., family = family, data = YMX))$coefficients
    est <- res[2:(length(ID_test) + 1), 1]  # the estimator for beta
    P_adjust_beta <- length(ID_test) * res[2:(length(ID_test) + 1), 4]  # the adjused p-value for beta (bonferroni)
    P_adjust_beta[P_adjust_beta > 1] <- 1
    P_fdr_beta <- p.adjust(res[2:(length(ID_test) + 1), 4], "fdr")  # the adjusted p-value for beta (FDR)
    
    ab_est <- alpha_est * beta_est
    
    ## Use the maximum value as p value 
    PA <- rbind(P_adjust_beta, P_adjust_alpha)
    P_value <- apply(PA, 2, max)
    
    FDRA <- rbind(P_fdr_beta, P_fdr_alpha)
    FDR <- apply(FDRA, 2, max)
    
    # Total effect
    if (is.null(COV.MY)) {
      YX <- data.frame(Y = Y, X = X)
    } else {
      YX <- data.frame(Y = Y, X = X, COV.MY)
    }
    
    gamma_est <- coef(glm(Y ~ ., family = family, data = YX))[2]
    
    results <- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                          `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                          `Bonferroni.p` = P_value, `BH.FDR` = FDR, check.names = FALSE)
    
    #message("Done!", "     (", Sys.time(), ")")
    
    doParallel::stopImplicitCluster()
    
    return(results)
  }
}

##############################
##### Utility Functions ######
##############################

# Internal function: parallel computing check

checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel & (ncore > 1)) {
    if (ncore > parallel::detectCores()) {
      message("You requested ", ncore, " cores. There are only ", 
              parallel::detectCores(), " in your machine!")
      ncore <- parallel::detectCores()
    }
    if (verbose) 
      message("    Running ", program.name, " with ", ncore, " cores in parallel...   (", 
              Sys.time(), ")")
    doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              Sys.time(), ")")
    registerDoSEQ()
  }
}



## Internal function: doOne code generater

doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", 
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ", 
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text, 
                   "]};invisible(b)}")
  return(script)
}



## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}



## Internal function: scale data (obsolete function)

scaleto <- function(dat) {
  if (is.null(dat)) 
    return(list(dn = NULL, d = NULL, ds = NULL))
  dat_scale <- scale(dat)
  dat_names <- names(dat)
  if (any(class(dat) %in% c("matrix", "data.frame", "data.table"))) {
    dat_names <- colnames(dat)
    dat <- as.matrix(data.frame(dat_scale))
  } else {
    dat_names <- names(dat)
    dat <- as.numeric(dat_scale)
  }
  dat_scale <- as.numeric(attributes(dat_scale)[["scaled:scale"]])
  return(list(dn = dat_names, d = dat, ds = dat_scale))
}



# Internal function: Sure Independent Screening
# Global variables:
globalVariables("n")
globalVariables("M_chunk")

himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel, ncore, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    conf.names <- colnames(COV)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(paste0("COV.", conf.names), collapse = "+")))))
  }
  
  doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                             glm.family, ", data = datarun))"), "c(1,4)")))
  
  checkParallel(tag, parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}
  
  colnames(results) <- M.names
  return(results)
}



# Internal function: LDPE
# the code of Liu han's JRSSB paper for high-dimensional Cox model
# ID: the index of interested parameter
# X: the covariates matrix with n by p
# OT: the observed time = min(T,C)
# status: the censoring indicator I(T <= C)

LDPE_func <- function(ID, X, OT, status){
  coi <- ID
  x <- X
  d <- dim(x)[2]
  n <- dim(x)[1]
  
  ##Set of tuning parameters
  PF <- matrix(1,1,d)
  PF[ID] <- 1
  fit <- glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  cv.fit <- cv.glmnet(x, survival::Surv(OT, status), family="cox", alpha = 1, standardize = FALSE,penalty.factor=PF)
  betas   <-   coef(fit, s = cv.fit$lambda.min)[1:d]  # the semi-penalized initial estimator  # initial estimator
  
  stime = sort(OT)          # Sorted survival/censored times
  otime = order(OT)         # Order of time
  
  Vs  = matrix(rep(0,d*d),nrow = d)
  Hs  = Vs                                 # Hessian
  ind = 0
  
  la  = rep(0,n)                           # Gradient w.r.t parameter of interest
  lb  = matrix(rep(0,(d-1)*n),nrow = n)    # Gradient w.r.t nuisance parameter (theta)
  i   = 1
  
  while( i<=n)
  {
    if (status[otime[i]]==1)
    {
      ind = which(OT >= stime[i])
      S0  = 0
      S1  = rep(0,d)
      S2  = matrix(rep(0,d*d),nrow = d)
      
      if (length(ind)>0)
      {
        for (j in 1:length(ind))
        {
          tmp = exp(x[ind[j],]%*%betas)
          S0  = S0 + tmp
          
          S1  = S1 + tmp %*%t(x[ind[j],])
          
          tmp = apply(tmp,1,as.numeric)     
          S2  = S2 + tmp*x[ind[j],]%*%t(x[ind[j],])
        }
      }
      S0 = apply(S0,1,as.numeric)
      
      la[i]  = -(x[otime[i],coi] - S1[coi]/S0)
      if (coi == 1)
      {
        lb[i,] = -(x[otime[i],c((coi+1):d)] - S1[c((coi+1):d)]/S0)
      } else if (coi == d){
        lb[i,] = -(x[otime[i],c(1:(coi-1))] - S1[c(1:(coi-1))]/S0)
      } else {
        lb[i,] = -(x[otime[i],c(1:(coi-1), (coi+1):d)] - S1[c(1:(coi-1), (coi+1):d)]/S0)
      }
      V   = S0*S2 - t(S1)%*%(S1)
      Hs  = Hs + V/(n*S0^2)          
    }
    i = i + 1
  }
  
  fit <- glmnet(lb,la,alpha = 1, standardize = FALSE,intercept = FALSE,lambda = sqrt(log(d)/n))
  what <- as.numeric(coef(fit)[2:d])   
  
  if (coi == 1)
  {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c((coi+1):d),coi]
  } else if (coi == d){
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1)),coi]
  } else {
    S = betas[coi] - (mean(la)  - t(what)%*%(colMeans(lb)))/(Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi])
    var   = Hs[coi,coi] - t(what)%*%Hs[c(1:(coi-1),(coi+1):d),coi]
  }
  
  beta_est <- S
  beta_SE  <- sqrt(1/(n*var))
  
  result <- c(beta_est,beta_SE)
  
  return(result)
}
