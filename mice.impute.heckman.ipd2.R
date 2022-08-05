
library(lme4)
library(data.table)
library(mvtnorm)
library(fMultivar)
library(data.table)
library(mgcv)
library(GJRM)
library(Matrix)
library(mixmeta)
library(mvtnorm)
library(pbivnorm)
library(mice)
library(micemd)
library(miceMNAR)
library(data.table)
library(vcd)
library(truncnorm)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(shiny)
library(DT)
library(mice)
library(fMultivar)
library(gridExtra)
library(ggpubr) 
library(reshape)

# Heckman model based imputation method ----

mice.impute.heckman.ipd2 <-function(y,ry,x,wy = NULL, type, full = FALSE, pmm = FALSE, meta_method,...) {
  # 1. Define variables and dataset----
  
  # Rename covariates
  
  colnames(x) <- paste0("x_", 1:length(colnames(x))) #change the covariates name for avoiding conflicts when y is covariate
  bos_name <- colnames(x)[type ==  1] # names of variables present in both outcome and selection model
  sel_name <- colnames(x)[type == -3] # names of variables in selection model alone
  out_name <- colnames(x)[type == -4] # names of variables in outcome model alone
  
  # # Define y type
  if (class(y) == "factor" & nlevels(y) == 2){
    warning("the missing variable is considered binomially distributed")
    family <- "binomial"
  }else{
    warning("the missing variable is considered normally distributed")
    family <- "gaussian"
  }
  
  
  # Check if group variable is defined
  if (length(colnames(x)[type == -2]) == 0) {
    warning("No group variable has been provided, the Heckman imputation model will be applied globally to the dataset.")
    Grp_est <- 0 # Group indicator 0: Heckman on full dataset, 1: Heckman at cluster level
  } else {
    group_name <- colnames(x)[type == -2]
    names.clust <- as.character(unique(x[, group_name]))
    Grp_est <- 1
  }
  
  #Define position of selection and outcome in coefficient vector
  order <- c(sel_name, bos_name, out_name)
  send <- length(sel_name) + length(bos_name) + 1 # Index where selection equation ends in coefficient vector
  oend <- send + length(bos_name) + length(out_name) + 1 # Index where outcome equation ends in coefficient vector
  
  # Define outcome and selection equation
  out <- as.formula(paste0("y", "~", paste(c(bos_name,out_name), collapse = "+")))
  sel <- as.formula(paste0("ry", "~", paste(c(sel_name, bos_name), collapse = "+")))
  
  # Define data & prediction matrix
  data <- data.frame(ry, y, x[, order])
  X <- data.frame(cbind(Int=rep(1, nrow(x)), x[, order]))
  
  # 2. Step 1: Get the theta estimators for each study ----
  Syst_nest <- Heck_est <- 1
  if (Grp_est == 1) { #Heckman at cluster level
    df_list <- split(data, x[, group_name]) 
    res_total <- suppressWarnings( lapply( df_list, copulaIPD,
                                           sel = sel, out = out, send = send, family = family)) #Calculate parameters for each study
    studytype <- sapply(res_total, function(x) x[[2]]) # Specify missing pattern in studies
    fit_list  <- lapply(res_total, function(x) x[[1]]) # Get the parameter estimates
    coef_list <- lapply(fit_list[studytype != 0], `[[`, c('coefficients')) # Get effect vector for each study
    Vb_list   <- lapply(fit_list[studytype != 0], `[[`, c('Vb')) # Get covariance matrix for each study
    coef_mat_s <- do.call(rbind, coef_list)
    
    varnam <- colnames(coef_mat_s)
    selnam <- varnam[grepl("*_s" , varnam)]
    outnam <- varnam[!varnam %in% c(selnam, "sigma.star", "theta.star")]
    
    # 3. Step 2: Get marginal theta and var(theta) ----
    
    if( length(studytype[studytype != 0]) < 2 ){
      Grp_est <- 0
    }else{
      Heck_mod <- get_marginal(coef_mat_s, Vb_list, full = full, selnam, outnam, meta_method)
      Heck_est <- Heck_mod$Mvma_est}
    Sys.nest <-as.numeric(length(studytype[studytype == 0])==0)
    
  } 
  
  if (Grp_est == 0 | Heck_est == 0){ # Heckman on full dataset or Heckman model no estimable
    warning("The Heckman model cannot be estimated marginally, so systematically missing groups will be imputed with the Heckman model based on the full dataset.")
    
    Heck_mod <- copulaIPD( data = data, sel = sel, out = out, family = family, send = send)
    
    if(Heck_mod[[2]] == 0 &(Grp_est==0|Syst_nest==0)){
      stop("There is insufficient information to impute the Heckman model at the marginal or study level.")
    }
  }
  
  
  # 4. Get theta_k and var(theta_k) from full conditional distribution ----
  if (Grp_est == 1) { # Applies imputation at cluster level
    for (i in names.clust) { #Loop across studies
      
      if (studytype[[i]] == 0) { #systematically missing
        star <- star_systematic( Heck_mod, send, oend, family)
        
      } else { # sporadically missing
        star <- star_sporadic( Heck_mod,
                               coef_list_i = coef_list[[i]],
                               Vb_list_i = Vb_list[[i]], selnam, outnam,family)
      }
      
      Xm <- X[!ry & x[, group_name] == as.numeric(i),]
      
      if (nrow(Xm) != 0) { #Cluster with at least one missing value in outcome equation
        y.star <- gen_y_star( Xm= Xm, sel_name = sel_name, bos_name = bos_name,
                              out_name =out_name, beta_s_star = star$beta_s_star,
                              beta_o_star = star$beta_o_star, sigma_star = star$sigma_star,
                              rho_star = star$rho_star,pmm = pmm, y = y, ry = ry)
        
        y[!ry & x[, group_name] == as.numeric(i)] <- y.star
      }
    }
  }else{ # Applies imputation on full dataset
    
    star<-star_systematic(Heck_mod, send, oend, family)
    Xm <- X[!ry,]
    
    if (nrow(Xm) != 0) { 
      y.star <- gen_y_star( Xm= Xm, sel_name = sel_name, bos_name = bos_name,
                            out_name =out_name, beta_s_star = star$beta_s_star,
                            beta_o_star = star$beta_o_star, sigma_star = star$sigma_star,
                            rho_star = star$rho_star, pmm = pmm, y = y, ry = ry)
      y[!ry] <- y.star
    }
  }
  
  return(y[!ry])
}



## F 0.1. CopulaIPD: Apply Binomial or Gaussian model depending on y type

copulaIPD <- function(data, sel, out, family, send) {
  
  fit_ind <- 0 # None model estimable for the cluster
  # A. Estimate Heckman model
  fit <- suppressWarnings(try(GJRM::gjrm( formula = list(sel, out),
                                          data = data,
                                          margins = c("probit", ifelse(family=="binomial","probit","N")),
                                          Model = "BSS",
                                          gamlssfit = TRUE,
                                          extra.regI = "sED",
                                          parscale = TRUE),
                              silent = TRUE))
  
  
  #if (fit_ind!=1) {
  if (!any(inherits(fit, "try-error"))) {
    # model is estimable
    convh <- round(min(eigen(fit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values),6) > 0 # convergence based on hessian positive definiteness
    convg <- round(max(abs(fit$fit$gradient)),1) < 10 # convergence based on abs max gradient
    
    #MAR indication 
    CIcon<-summary(fit)$ CItheta 
    MNAR_ind<-!(abs(CIcon[[1]]-CIcon[[2]])<0.001&CIcon[[1]]<0&CIcon[[2]]>0) # exclusion of cases that blow up variance
    
    if (MNAR_ind){
      fit_ind <- 2
      if (convh & convg) {
        # MNAR estimable
        fit_ind <- 1 # Heckman model estimable for the cluster
      }
    }
  }
  
  if( fit_ind ==2){
    fit<-NULL
    gam1 <- try(mgcv::gam(formula=sel,data = data,family = "binomial", method ="REML"),silent = TRUE)
    gam2 <- try(mgcv::gam(formula=out,data = data,family = family, method ="REML"),silent = TRUE)
    
    if(!any(inherits(gam1, "try-error"))&!any(inherits(gam2, "try-error"))){
      coefficients <- c(gam1$coefficients,gam2$coefficients)
      if(all(!is.na(coefficients))){
        s     <- ifelse(family != "binomial",1,0) 
        ncol1 <- ncol(gam1$Vp)
        ncol2 <- ncol(gam2$Vp)
        Vb    <- matrix(0,ncol = ncol1+ncol2+1+s, nrow = ncol1+ncol2+1+s)
        Vb[1:ncol1,1:ncol1] <- gam1$Vp
        Vb[(ncol1+1):(ncol1+ncol2),(ncol1+1):(ncol1+ncol2)] <- gam2$Vp
        
        if (family != "binomial") {
          coefficients <- c(coefficients, sigma.star = log(sqrt(gam2$scale)))
          Vb[(ncol1+ncol2+1),(ncol1+ncol2+1)] <- gam2$V.sp}
        
        fit$coefficients <- c(coefficients,theta.star=0)
        fit$Vb  <- Vb
        fit_ind <- 2}
    }
    
  }
  
  if (fit_ind != 0) {
    names <- c(paste0(names(fit$coefficients)[1:send], "_s"),
               names(fit$coefficients[(send + 1):length(names(fit$coefficients))]))
    names(fit$coefficients) <- names
    colnames(fit$Vb) <- names
    rownames(fit$Vb) <- names
  }else{
    fit <- NA}  
  
  return(list(fit, fit_ind))
  
}

## F 0.2. cov_mat_vector: Transform covariance matrix in a ordered vector 
cov_mat_vector <- function(cov_mat, vnames) {
  cov_mat[upper.tri(cov_mat)] <- "Up"
  cov_vec <- as.vector(cov_mat[vnames, vnames])
  cov_vec <- as.numeric(cov_vec[cov_vec != "Up"])
  return(cov_vec)
  
}
## F 0.3. draw_theta_psi_mar: Estimate true effect size and draw a marginal theta and psi=var(theta) .

draw_theta_psi_mar <- function(coef_mat_s, Vb_list, meta_method, Mvma_est, vnames = NULL) {
  
  theta_star <- NA
  psi_star <- NA
  
  if (is.null(vnames)) { #use all set of parameters
    vnames <- colnames(coef_mat_s)
  }
  
  # Get covariance matrix
  coef_mat_s <- coef_mat_s[, vnames]
  cov_mat_s <- do.call("rbind", lapply(Vb_list, cov_mat_vector, vnames = vnames))
  
  # Apply multivariate random-effects meta-analysis
  mvma <- suppressWarnings(try(mixmeta::mixmeta(coef_mat_s,cov_mat_s, method = meta_method,
                                                control = list(hessian = TRUE)), silent = TRUE))
  
  if (inherits(mvma,"try-error")) { # Use mm instead
    meta_method = "mm"
    mvma <- suppressWarnings(try(mixmeta::mixmeta(coef_mat_s, cov_mat_s, method = meta_method,
                                                  control = list(hessian = TRUE)),silent = TRUE))
    
    if (inherits(mvma,"try-error")) { # MA can not be estimated
      Mvma_est <- 0}
  }
  
  
  if (Mvma_est == 1) {
    # Draw effects theta_star
    theta_star <- MASS::mvrnorm(n = 1, mu = coef(mvma), Sigma = vcov(mvma))
    
    if (meta_method != "mm") {
      # Draw random effect, psi_star
      if (length(vnames) == 1) {
        qrsh <- 1 / mvma$hessian
      } else {
        Hes <- as.matrix(Matrix::forceSymmetric(mvma$hessian))
        qrsh <- as.matrix(Matrix::nearPD(MASS::ginv(-Hes))$mat)
      }
      
      rpar <- mvtnorm::rmvnorm(1, mean = mvma$par, sigma = qrsh, method = "svd")
      
      if (length(vnames) == 1) {
        psi_star <- rpar ^ 2
      } else {
        psi <- matrix(0, ncol(mvma$Psi), ncol(mvma$Psi))
        psi[lower.tri(psi, diag = TRUE)] <- rpar
        psi_star <- Matrix::tcrossprod(psi)
        
      }
      
    } else # meta_method== reml OR ml
      psi_star <- mvma$Psi
  }
  
  colnames(psi_star) <- names(theta_star)
  rownames(psi_star) <- names(theta_star)
  
  return(list(theta_star, psi_star, Mvma_est))
  
}

## F 0.4 Conditional posterior distribution
draw_cond_theta <- function(theta_mar, theta_k, var_theta_k, vnames) {
  W_m <- MASS::ginv(theta_mar[[2]])
  W_k <- MASS::ginv(var_theta_k[vnames, vnames])
  S <- MASS::ginv(W_m + W_k)
  mu <- S %*% (W_k %*% as.vector(theta_k[vnames]) + W_m %*% as.vector(theta_mar[[1]]))
  theta_star_i <- MASS::mvrnorm(n = 1, mu = as.vector(mu), Sigma = S)
  return(theta_star_i)
}


## F 0.5 Get marginal draws
get_marginal <- function(coef_mat_s, Vb_list, full, selnam, outnam, meta_method ){
  
  beta_s = beta_o = rho_t = sigma_t = NA
  Mvma_est <- 1
  # Draw thetas from the universe of populations
  if (full == TRUE) { # use the entire set of parameters
    
    total_mar <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                     Vb_list = Vb_list,
                                     meta_method = meta_method,
                                     Mvma_est = Mvma_est) # if mvma was estimable on previous parameters (applicable mainly Full=F)
    Mvma_est <- total_mar[[3]]
    
    if (Mvma_est == 1) {# if marginal parameters where estimated proceed with Heckman cluster estimation otherwise opt for the Heckman full model
      
      beta_s <- list(total_mar[[1]][selnam], total_mar[[2]][selnam, selnam])
      beta_o <- list(total_mar[[1]][outnam], total_mar[[2]][outnam, outnam])
      rho_t <- list(total_mar[[1]]["theta.star"],
                    total_mar[[2]]["theta.star", "theta.star"])
      if ("sigma.star"%in%colnames(coef_mat_s)) {
        sigma_t <- list(total_mar[[1]]["sigma.star"],
                        total_mar[[2]]["sigma.star", "sigma.star"])
      }
    }
    
  } else{ # separate the set of parameters in beta_out, beta_s and rho
    
    beta_o <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                  Vb_list = Vb_list,
                                  vnames = outnam,
                                  meta_method = meta_method,
                                  Mvma_est = Mvma_est)
    beta_s <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                  Vb_list = Vb_list,
                                  vnames = selnam,
                                  meta_method = meta_method,
                                  Mvma_est = beta_o[[3]])
    rho_t <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                 Vb_list = Vb_list,
                                 vnames = "theta.star",
                                 meta_method = meta_method,
                                 Mvma_est = beta_s[[3]]) #copula package calls atanh(rho) as theta.star
    
    if ("sigma.star"%in%colnames(coef_mat_s)) {
      sigma_t <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                     Vb_list = Vb_list,
                                     vnames = "sigma.star",
                                     meta_method = meta_method,
                                     Mvma_est = rho_t[[3]]) #copula package calls log(sigma) as sigma.star
      Mvma_est <- sigma_t[[3]]
    }
    
  }
  
  return (list ( beta_s = beta_s,
                 beta_o = beta_o,
                 rho_t = rho_t,
                 sigma_t = sigma_t,
                 Mvma_est = Mvma_est))
}



## F 0.6 Get draw from systematically missing groups
star_systematic <- function(Heck_mod, send, oend, family){
  
  if (is.null(Heck_mod$Mvma_est)){ # From total data model
    star <- mvtnorm::rmvnorm( 1, mean = Heck_mod[[1]]$coefficients,
                              sigma = Heck_mod[[1]]$Vb, method = "svd")
    
    beta_o_star <- star[(send + 1):oend]
    beta_s_star <- star[1:send]
    rho_star <- tanh(star[, "theta.star"])
    sigma_star <- ifelse(family=="gaussian",
                         exp(star[, "sigma.star"]), NA)
    
    
  } else { # From meta model
    
    beta_o_star <- MASS::mvrnorm(n = 1, mu = Heck_mod$beta_o[[1]],Sigma =Heck_mod$beta_o[[2]])
    beta_s_star <- MASS::mvrnorm(n = 1, mu = Heck_mod$beta_s[[1]],Sigma = Heck_mod$beta_s[[2]])
    rho_star <- tanh(MASS::mvrnorm( n = 1, mu = Heck_mod$rho_t[[1]],Sigma = Heck_mod$rho_t[[2]]))
    if (family =="gaussian") {
      sigma_star <- exp(MASS::mvrnorm( n = 1, mu = Heck_mod$sigma_t[[1]],
                                       Sigma = Heck_mod$sigma_t[[2]])) #copula package calls log(sigma) as sigma.star
    } else { #binomial
      sigma_star <- NA}
  }
  
  return (list (beta_s_star = beta_s_star,
                beta_o_star = beta_o_star,
                rho_star = rho_star,
                sigma_star = sigma_star))}


## F 0.7 Get draw from sporadically missing groups
star_sporadic <- function(Heck_mod, coef_list_i, Vb_list_i, selnam, outnam, family){
  
  beta_s_star = beta_o_star = sigma_star = rho_star = NA
  cond <- ifelse(is.null(Heck_mod$Mvma_est),0,Heck_mod$Mvma_est)
  
  if (cond == 1){ # Draw study parameters from conditional distribution
    beta_o_star <- draw_cond_theta( theta_mar = Heck_mod$beta_o, theta_k = coef_list_i,
                                    var_theta_k = Vb_list_i, vnames = outnam)
    beta_s_star <- draw_cond_theta( theta_mar = Heck_mod$beta_s,theta_k = coef_list_i,
                                    var_theta_k = Vb_list_i, vnames = selnam)
    
    rho_star <- tanh( draw_cond_theta( theta_mar = Heck_mod$rho_t,
                                       theta_k = coef_list_i,
                                       var_theta_k = Vb_list_i,
                                       vnames = "theta.star")) #copula package calls atanh(rho) as theta.star
    
    if (family == "gaussian") {
      sigma_star <- exp(draw_cond_theta( theta_mar = Heck_mod$sigma_t,
                                         theta_k = coef_list_i,
                                         var_theta_k = Vb_list_i,
                                         vnames = "sigma.star"))#copula package calls log(sigma) as sigma.star
    }
    
    
  } else { # Draw from study parameters from study distribution
    
    beta_o_star <- MASS::mvrnorm(n = 1,mu = coef_list_i[outnam],Sigma = Vb_list_i[outnam, outnam])
    beta_s_star <- MASS::mvrnorm(n = 1,mu = coef_list_i[selnam],Sigma = Vb_list_i[selnam, selnam])
    rho_star <- tanh(MASS::mvrnorm( n = 1, mu = coef_list_i["theta.star"],
                                    Sigma = Vb_list_i["theta.star", "theta.star"]))
    if (family == "gaussian") {
      sigma_star <- exp(MASS::mvrnorm( n = 1, mu = coef_list_i["sigma.star"],
                                       Sigma = Vb_list_i["sigma.star", "sigma.star"]))#copula package calls log(sigma) as sigma.star
    }
    
  }
  
  return (list ( beta_s_star = beta_s_star,
                 beta_o_star = beta_o_star,
                 rho_star = rho_star,
                 sigma_star = sigma_star))
}


## F 0.8 Generate the imputation values

gen_y_star <- function(Xm, sel_name, bos_name, out_name, beta_s_star, beta_o_star,
                       sigma_star,rho_star, pmm, y, ry) {
  
  XOBO <- data.matrix(Xm[,colnames(Xm) %in% c("Int",bos_name,out_name)]) %*% as.vector(beta_o_star)
  XSBS <- data.matrix(Xm[,colnames(Xm) %in% c("Int",sel_name,bos_name)]) %*% as.vector(beta_s_star)
  
  if (!is.na(sigma_star)) { # normal missing variable
    
    Ratio <- (-stats::dnorm(XSBS) / (stats::pnorm(-XSBS)))
    Ratio[is.na(Ratio) | is.infinite(Ratio)] <- 0.0
    y.star <- XOBO + as.numeric(sigma_star) * as.numeric(rho_star) * Ratio +
      rnorm(nrow(XSBS), 0, sd = sigma_star)
    
    if (pmm == TRUE) {
      idx <- mice::matchindex(y[ry == 1], y.star)
      y.star <- y[ry == 1][idx]
    }
    
  } else { #binomial missing variable
    
    p.star <- pbivnorm::pbivnorm(as.vector(XOBO),-as.vector(XSBS),
                                 -as.numeric(rho_star)) / stats::pnorm(-XSBS)
    p.star[is.na(p.star) | (is.infinite(p.star) & p.star < 0) |p.star < 0.0 |
             p.star == "0"] <- 0.0
    p.star[p.star > 1.0 | p.star == "1" |(is.infinite(p.star) & p.star > 0)] <- 1.0
    
    y.star <-rep(levels(y)[1],nrow(XOBO))
    y.star[runif(nrow(XOBO)) < p.star]<-levels(y)[2]
    
  }
  
  return(y.star)
}


# Rubin rules functions ----

logit <- function(x) {
  log(x/(1-x))
}
inv.logit <- function(x) {
  1/(1+exp(-x))
}
inv.logit.SE<- function(logit.se, c) {
  logit.se * (c*(1-c))
}

f.incidence <- function(variable,n.data) {
  a<-summary(variable)
  b<-prop.test(x = a[[2]], n = n.data, correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
  return(c(prop,ci))
}


f.abs.perstudy<-function(data, outcome_name) {
  m <-length(unique(data$.imp)) # imputed datasets
  dataset.n<-nrow(data[data$.imp==1,])
  #Calculate absolute risk of outcome per study
  for (i in 1:m) {
    d<-data[data$.imp==i,]
    d$outcome <- d[, outcome_name]
    d$outcome<-as.factor(d$outcome)
    inc<-as.data.frame(d %>% 
                         summarise(p=f.incidence(outcome,dataset.n)[1],
                                   ci.l=f.incidence(outcome,dataset.n)[2],
                                   ci.u=f.incidence(outcome,dataset.n)[3]))
    if(i==1) {inc.outcome<-inc}
    if(i>1) {inc.outcome<-rbind(inc.outcome,inc)}
  }
  
  
  #Replace 0's by 0.000001
  colnames(inc.outcome)<-c("incidence","ci.l","ci.u")
  inc.outcome$incidence[inc.outcome$incidence==0]<-0.000001
  inc.outcome$ci.l[inc.outcome$ci.l==0]<-0.000001
  inc.outcome$mse<-(inc.outcome$incidence-mean(inc.outcome$incidence))^2
  
  #Logit transformation
  inc.outcome$logit.incidence<-logit(inc.outcome$incidence)
  inc.outcome$logit.ci.l<-logit(inc.outcome$ci.l)
  inc.outcome$logit.ci.u<-logit(inc.outcome$ci.u)
  inc.outcome$logit.se<-(inc.outcome$logit.ci.u-inc.outcome$logit.ci.l)/(2*1.96)
  
  #Pool rubin
  pool.rubin<- data.frame(matrix(NA, nrow = 1, ncol = 7))
  colnames(pool.rubin)<-c("logit.abs","within","between","logit.var",
                          "logit.se","logit.ci.lb","logit.ci.ub")
  
  pool.rubin$logit.abs <- mean(inc.outcome$logit.incidence)
  pool.rubin$within <- mean(inc.outcome$logit.se^2)
  pool.rubin$between <- (1 + (1/m)) * var(inc.outcome$logit.incidence)
  pool.rubin$logit.var <- pool.rubin$within + pool.rubin$between
  pool.rubin$logit.se <- sqrt(pool.rubin$logit.var)
  pool.rubin$logit.ci.lb <- pool.rubin$logit.abs + qnorm(0.05/2)     * pool.rubin$logit.se
  pool.rubin$logit.ci.ub <- pool.rubin$logit.abs + qnorm(1 - 0.05/2) * pool.rubin$logit.se
  
  pool.rubin$incidence<-inv.logit(pool.rubin$logit.abs)*100
  pool.rubin$ci.lb<-inv.logit(pool.rubin$logit.ci.lb)*100
  pool.rubin$ci.ub<-inv.logit(pool.rubin$logit.ci.ub)*100
  pool.rubin$mcse<- sqrt(sum(inc.outcome$mse)/(m*(m-1))) 
  
  return(pool.rubin)
}

#  Define dataset --

data("hiv")
hiv<-setDT(hiv)[region==5,]
hiv$interviewerID<-as.factor(as.character(hiv$interviewerID))
ID_mixed = glmer(hivconsent ~(1 | interviewerID), data = hiv, family = binomial,
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
reffect<-data.table(ranef(ID_mixed)$interviewerID)
reffect$interviewerID<-levels(hiv$interviewerID)
colnames(reffect)<-c("IDreffect","interviewerID")
hiv<-merge(hiv,reffect,by="interviewerID",all.x=TRUE)

out <- gjrm(list(hivconsent ~ marital +  s(age)+condom+highhiv+IDreffect,
                 hiv ~ marital +highhiv+  condom +s(age)),
            data = hiv, Model = "BSS",
            margins = c("probit", "probit"))
xmatS<-as.matrix(predict(out,1,hiv,type = "lpmatrix"))
xmatO<-as.matrix(predict(out,2,hiv,type = "lpmatrix"))
betaS<-as.matrix(coefficients(out)[1:15])
betaO<-as.matrix(coefficients(out)[16:29])

datasim<-as.data.table(xmatS)
colnames(datasim)[7:15]<-paste0("sage",1:9)
datasim[,XOBO := xmatO%*%betaO]
datasim[,XSBS:= xmatS%*%betaS]
datasim[,marital:=hiv$marital]
datasim[,condom:=hiv$condom]
datasim[,highhiv:=hiv$highhiv]

save(datasim,file="datasim.Rdata")


# Functions for shiny app ----


Bimputationres<-function(ldatasim){
  res1<-Imputationres(datasim=as.data.table(ldatasim[[1]]))
  #res1$error<-"Normal"
  #res2<-Imputationres(datasim=as.data.table(ldatasim[[2]]))
  #res2$error<-"Skewed"
  #res<-rbind(res1)

 
  #colnames(tablefin)<-c("Imputation method","HIV prevalence","LCI","UCI","Missing(%)","Bias","Monte Carlo SE")
  res1
}

Imputationres<-function(datasim){
  
  dataSim0 <- datasim[, c("marital","condom","highhiv","IDreffect","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim0$ybin<-as.factor(datasim$ybin)
  ini <- mice(dataSim0, maxit = 0)
  meth0<-ini$method
  pred0 <- ini$pred
  meth0[c("ybin")]<-"heckman.ipd2"
  pred0["ybin","IDreffect"] <- -3
  
  dataSimNERV <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSimNERV$ybin<-as.factor(datasim$ybin)
  ini <- mice(dataSimNERV, maxit = 0)
  methNERV<-ini$method
  predNERV <- ini$pred
  methNERV[c("ybin")]<-"heckman.ipd2"
  
  
  dataSim1 <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim1$ex1<-ifelse(datasim$incen1!=0,1,0)
  dataSim1$ybin<-as.factor(datasim$ybin1)
  ini <- mice(dataSim1, maxit = 0)
  meth1<-ini$method
  pred1 <- ini$pred
  meth1[c("ybin")]<-"heckman.ipd2"
  pred1["ybin","ex1"] <- -3
  
  dataSim0i <- datasim[, c("marital","condom","highhiv","IDreffect","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim0i$ex1<-ifelse(datasim$incen1!=0,1,0)
  dataSim0i$ybin<-as.factor(datasim$ybin1)
  ini <- mice(dataSim0i, maxit = 0)
  meth0i<-ini$method
  pred0i <- ini$pred
  meth0i[c("ybin")]<-"heckman.ipd2"
  pred0i["ybin","IDreffect"] <- -3
  pred0i["ybin","ex1"] <- -3
  
  #dataSim2 <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  #dataSim2$ex2<-ifelse(datasim$incen2!=0,1,0)
  #dataSim2$ybin<-as.factor(datasim$ybin2)
  #ini <- mice(dataSim2, maxit = 0)
  #meth2<-ini$method
  #pred2 <- ini$pred
  #meth2[c("ybin")]<-"heckman.ipd2"
  #pred2["ybin","ex2"] <- -3
  
  dataSim12 <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim12$incen<-as.factor(datasim$incen)
  dataSim12$ybin<-as.factor(datasim$ybin12)
  ini <- mice(dataSim12, maxit = 0)
  meth12<-ini$method
  pred12 <- ini$pred
  meth12[c("ybin")]<-"heckman.ipd2"
  pred12["ybin","incen"] <- -3
  
  
  imp_mar <- mice( dataSim0, # dataset with missing values
                   m = 15,   # number of imputations
                   seed = 1234, #seed attached to the dataID
                   print = T,
                   maxit=1)
  
  
  imp_heck_NERV <- suppressWarnings(mice( dataSimNERV, # dataset with missing values
                                          m = 15,   # number of imputations
                                          seed = 1234, #seed attached to the dataID
                                          meth = methNERV, #imputation method vector
                                          pred = predNERV, #imputation predictors matrix
                                          print = T,
                                          maxit=1,
                                          full=FALSE,
                                          meta_method="reml",
                                          pmm=FALSE))
  
  imp_heck_iid <- suppressWarnings(mice( dataSim0, # dataset with missing values
                                         m = 15,   # number of imputations
                                         seed = 1234, #seed attached to the dataID
                                         meth = meth0, #imputation method vector
                                         pred = pred0, #imputation predictors matrix
                                         print = T,
                                         maxit=1,
                                         full=FALSE,
                                         meta_method="reml",
                                         pmm=FALSE))
  
  imp_heck_inc1 <- suppressWarnings(mice( dataSim1, # dataset with missing values
                                          m = 15,   # number of imputations
                                          seed = 1234, #seed attached to the dataID
                                          meth = meth1, #imputation method vector
                                          pred = pred1, #imputation predictors matrix
                                          print = T,
                                          maxit=1,
                                          full=FALSE,
                                          meta_method="reml",
                                          pmm=FALSE))
  
  imp_heck_iidi <- suppressWarnings(mice( dataSim0i, # dataset with missing values
                                          m = 15,   # number of imputations
                                          seed = 1234, #seed attached to the dataID
                                          meth = meth0i, #imputation method vector
                                          pred = pred0i, #imputation predictors matrix
                                          print = T,
                                          maxit=1,
                                          full=FALSE,
                                          meta_method="reml",
                                          pmm=FALSE))
  
  
  #imp_heck_inc2 <- suppressWarnings(mice( dataSim2, # dataset with missing values
  #                                        m = 15,   # number of imputations
  #                                        seed = 1234, #seed attached to the dataID
  #                                        meth = meth2, #imputation method vector
  #                                        pred = pred2, #imputation predictors matrix
  #                                        print = T,
  #                                        maxit=1,
  #                                        full=FALSE,
  #                                        meta_method="reml",
  #                                        pmm=FALSE))
  
  imp_heck_inc12 <- suppressWarnings(mice( dataSim12, # dataset with missing values
                                           m = 15,   # number of imputations
                                           seed = 1234, #seed attached to the dataID
                                           meth = meth12, #imputation method vector
                                           pred = pred12, #imputation predictors matrix
                                           print = T,
                                           maxit=1,
                                           full=FALSE,
                                           meta_method="reml",
                                           pmm=FALSE))
  
  
  mar<-f.abs.perstudy(data=complete( imp_mar, "long"),outcome_name<-"ybin")
  heckNERV<-f.abs.perstudy(data=complete( imp_heck_NERV, "long"),outcome_name<-"ybin")
  heckiid<-f.abs.perstudy(data=complete( imp_heck_iid, "long"),outcome_name<-"ybin")
  heckinc1<-f.abs.perstudy(data=complete( imp_heck_inc1, "long"),outcome_name<-"ybin")
  heckiidi<-f.abs.perstudy(data=complete( imp_heck_iidi, "long"),outcome_name<-"ybin")
  #heckinc2<-f.abs.perstudy(data=complete( imp_heck_inc2, "long"),outcome_name<-"ybin")
  heckinc12<-f.abs.perstudy(data=complete( imp_heck_inc12, "long"),outcome_name<-"ybin")
  
  perc<-c(sum(is.na(dataSim0$ybin))/nrow(dataSim0),
          sum(is.na(dataSimNERV$ybin))/nrow(dataSimNERV),
          sum(is.na(dataSim0$ybin))/nrow(dataSim0),
          sum(is.na(dataSim1$ybin))/nrow(dataSim1),
          sum(is.na(dataSim0i$ybin))/nrow(dataSim0i),
          #sum(is.na(dataSim2$ybin))/nrow(dataSim2),
          sum(is.na(dataSim12$ybin))/nrow(dataSim12))
  
  tableimp<-rbind(mar, heckNERV,heckiid,heckinc1,heckiidi,heckinc12)[,c("incidence","ci.lb","ci.ub","mcse")]
  Imputation<-c("MAR", "Heckman no ERV","Heckman interviewer ID","Heckman single incentive low","Heckman interviewer ID+incentive","Heckman scaled incentive")
  tablefin<-as.data.table(cbind(Imputation, tableimp))
  tablefin[,bias:=round(incidence-mean(datasim$ybin0)*100,2)]
  tablefin[,perc:=round(perc*100,2)]
  tablefin[,ci.lb:=round(ci.lb,2)]
  tablefin[,ci.ub:=round(ci.ub,2)]
  tablefin<-tablefin[,c("Imputation","incidence","ci.lb","ci.ub","perc","bias")]
  colnames(tablefin)<-c("Imputation method","HIV prevalence","LCI","UCI","missing(%)","Bias")
  
  return( tablefin)
  }


data_error<-function(datasim, error, prop=NULL, incentv1,incentv2){
  n= nrow(datasim)
  datasim[,eS:=error[,1]]
  datasim[,eO:=error[,2]]
  datasim[,y.star := XOBO+eO]
  datasim[,ry.star:= XSBS+eS]
  datasim[,ybin0:=ifelse(y.star>0,1,0)]
  datasim[,ry:=ifelse(ry.star>0,1,0)]
  datasim[,ybin:=ifelse(ry==1,ybin0,NA)]
  datasim[,observedo:=ifelse(ry.star<=0,"Unobserved","Observed")] 
  datasim[,observedo:=factor(observedo,levels=c("Observed","Unobserved"))]
 
  datasim[, incenc := as.factor(sample(x = c("None", "Small", "Large"), #previous DM treatment efficacy
                                       size = n, replace = TRUE, 
                                       prob = c(prop[1], prop[2]-prop[1], 1-prop[2])))]
  datasim[, incen:= ifelse(incenc=="Small",incentv1,ifelse(incenc=="Large",incentv2,0))]
  datasim[, incen1:= ifelse(incenc%in%c("Small","Large"),incentv1,0)]
  datasim[, incen2:= ifelse(incenc%in%c("Small","Large"),incentv2,0)]
  datasim[,ry1:=ifelse(ry.star+incen1<=0,0,1)]
  datasim[,ybin1:=ifelse(ry1==1,ybin0,NA)]
  datasim[,ry2:=ifelse(ry.star+incen2<=0,0,1)]
  datasim[,ybin2:=ifelse(ry2==1,ybin0,NA)]
  datasim[,ry12:=ifelse(ry.star+incen<=0,0,1)]
  datasim[,ybin12:=ifelse(ry12==1,ybin0,NA)]
  datasim[,yobs:=ifelse(ry12==1,y.star,NA)]
  datasim[,observed:=ifelse(ry.star+incen<=0,"Unobserved","Observed")]
  datasim[,incenc:=factor(incenc,levels=c("None","Small","Large"))]
  datasim[,observed:=factor(observed,levels=c("Observed","Unobserved"))]
  datasim
}


data_mod <- function(rho, prop, incentv1,incentv2) {
  load("datasim.Rdata")
  set.seed(12345)
  rho = ifelse(rho == -1,-0.99,ifelse(rho==1,0.99,rho)) # has to be slightly modified due to simulation of skewed errors
  n= nrow(datasim)
  d <- diag(2)
  d[2,1] <- d[1,2] <- rho
  
  e_norm <- mvtnorm::rmvnorm(n = n, mean = rep(0,2), sigma = d)
  e_skew <- fMultivar::rmvst(n = n, dim = 2, mu = c(0,0), Omega = d, alpha = c(-2,1), df = 8)

  datasim_norm <- data_error(datasim, error=e_norm, prop, incentv1,incentv2)
  datasim_skew <- data_error(datasim, error=e_skew, prop, incentv1,incentv2)
  
  ldatasim<-list(datasim_norm,datasim_skew)
  ldatasim
}


  
## Shiny app1----

ui1 <- fluidPage(
  titlePanel("Variation of rho"),
  
  fluidRow(
    column(6, sliderInput("x", "rho:",
                          min=-1, max=1, value=-1, step=0.1)),
    column(6,tableOutput("table_noinc")) 
  ),
  
  fluidRow(
    column(6, plotOutput("plot_error")),
    column(6, plotOutput("plot_noinc"))
  )
)



server1 <- function(input, output) {
  
  data.sel <- reactive({
    rho <- input$x
    prop<-c(0.85,0.95)
    datasim<-data_mod(rho=rho,prop=prop,incentv1=1,incentv2=2)
    datasim[[1]]
  })
  
  output$plot_error <- renderPlot({ 
    data<-as.data.table(data.sel())
    ggplot(data= as.data.frame(data), aes(eS,eO))+
      geom_point(size=1,colour="lightblue",alpha=0.7)+
      geom_smooth(size=0.5,colour="black",se = F)+
      xlab(expression("Error in the selection equation (" *epsilon[i] **'S'*")"))+
      ylab(expression("Error in the outcome equation (" *epsilon[i] **'O'*")"))+
      theme_minimal()+
      ggtitle("Error terms of both equations")
  })
  
  output$plot_noinc <- renderPlot({ 
    data<-as.data.table(data.sel())
    limry<-max(abs(max(data$ry.star)),abs(min(data$ry.star)))
    limy<-max(abs(max(data$y.star)),abs(min(data$y.star)))
    group.colors <- c("#7FC97F","#386CB0","#F0027F") 
    ggplot(data= as.data.frame(data), aes(ry.star,y.star))+
      labs(shape="Patient")+
      scale_shape_manual(values=c(19,1))+
      scale_y_continuous(breaks = c(0)) +
      scale_x_continuous(breaks = c(0))+ 
      annotate("text", y = -limy/2, x = (-limy-1), label = "HIV-",angle=90,size=2) +
      annotate("text", y =  limy/2, x = (-limy-1), label = "HIV+",angle=90,size=2)+
      coord_cartesian(xlim = c(-limry,limry), ylim = c(-limy,limy),clip = "off")+
      xlab(expression('Latent propensity of observability '(r[i]*'*') %->%''))+
      ylab(expression('Latent test result ' (y[i]*'*') %->%''))+
      theme_minimal()+
      geom_point(size=1,aes(shape=observedo),alpha=0.7,colour="#7FC97F")+
      geom_hline(yintercept=0, linetype="dashed", color = "black")+
      geom_vline(xintercept=0, linetype="dashed", color = "black")+ 
      geom_smooth(size=0.5,colour="black",se = F)+
      annotate("text", y =  limy*2/3, x =0.2, label = "Threshold none",angle=90,size=2,color = "#7FC97F")+
      ggtitle("Latent responses of both equations")+ theme( legend.position = 'right')
  })
  

  
  output$table_noinc <- renderTable({
    data<-as.data.table(data.sel())
    inceno<-as.data.table(table(data$observedo,data$ybin0))
    inceno[,V2:=ifelse(V2==1,"HIVp","HIVtn")]
    inceno<-cast(inceno, V2 ~ V1,value="N")
    inceno$All<-inceno$Observed+inceno$Unobserved
    inceno<-as.data.table(melt(setDT(inceno), id.vars = c("V2")))
    inceno<-as.data.table(cast(inceno, variable ~ V2))
    inceno[,HIV.prev:=round(HIVp/(HIVp+HIVtn)*100,2)]
    colnames(inceno)<-c("Patients","HIV+","HIV-","HIV+(%)")
    inceno
  })
  
}  

imp_shiny1 = function() {
  shinyApp(ui1, server1)
}


## Incentive application

ui2 <- fluidPage(
  titlePanel("Incentive as exclusion restriction variable"),
  
  fluidRow(
    column(6, sliderInput("x", "rho:",min=-1, max=1, value=1, step=0.1)),
    column(6,sliderInput("p","Proportion incentive",min = 0,max = 1,value = c(0.85,0.95),step=0.01)) 
  ),
  fluidRow(
    column(6,tableOutput("table_inc")),
    column(6,tableOutput("table_imp"))       
  ),
  fluidRow(
    column(12, plotOutput("plot_inc"))
  )
)

server2 <- function(input, output) {

  data.sel <- reactive({
    rho  <- input$x
    prop <- input$p
    datasim <- data_mod(rho=rho,prop=prop,incentv1=1,incentv2=2)
    datasim
  })
  
  output$plot_inc <- renderPlot({ 
    datasim<-as.data.table(data.sel()[[1]])
    incentv1=max(datasim$incen1)
    incentv2=max(datasim$incen2)
    limry<-max(abs(max(datasim$ry.star)),abs(min(datasim$ry.star)))
    limy<-max(abs(max(datasim$y.star)),abs(min(datasim$y.star)))
    group.colors <- c("#7FC97F","#386CB0","#F0027F") 
    ggplot(data= as.data.frame(datasim), aes(ry.star,y.star))+
      labs(x ="latent propensity of observability(ry*) -->", y = "latent outcome response (y*) --->")+
      labs(shape="Measurement")+
      scale_shape_manual(values=c(19,1))+
      scale_y_continuous(breaks = c(0)) +
      scale_x_continuous(breaks = c(0))+ 
      annotate("text", y = -limy/2, x = (-limy-1), label = "HIV-",angle=90,size=2) +
      annotate("text", y =  limy/2, x = (-limy-1), label = "HIV+",angle=90,size=2)+
      coord_cartesian(xlim = c(-limry,limry), ylim = c(-limy,limy),clip = "off")+
      xlab(expression('Latent propensity of observability '(R[i]*'*') %->%''))+
      ylab(expression('Latent test result ' (Y[i]*'*') %->%''))+
      theme_minimal()+
      geom_point(size=1,aes(colour=incenc,shape=observed),alpha=0.7)+
      geom_smooth(size=0.5,colour="black",se = F)+
      geom_hline(yintercept=0, linetype="dashed", color = "black")+
      geom_vline(xintercept=0, linetype="dashed", color = "black")+ 
      geom_vline(xintercept=-incentv1, linetype="dashed", color = "black")+ 
      geom_vline(xintercept=-incentv2, linetype="dashed", color = "black")+
      annotate("text", y =  limy*2/3, x =0.2, label = "Threshold none",angle=90,size=2,color = "#7FC97F")+
      annotate("text", y =  limy*2/3, x =(0.2-incentv1), label = "Threshold small",angle=90,size=2, color = "#386CB0")+
      annotate("text", y =  limy*2/3, x =(0.2-incentv2), label = "Threshold large",angle=90,size=2, color = "#F0027F")+
      labs(shape="Measurement",color="Incentive")+
      scale_colour_manual(values=group.colors)+
      ggtitle("Survey with incentive")+theme_minimal()+
      theme( legend.position = 'bottom')})

  output$table_inc <- renderTable({
    datasim<-as.data.table(data.sel()[[1]])
    incen.table<-as.data.table(table(datasim$observed,datasim$ybin0,datasim$incenc))
    colnames(incen.table)<-c("Measurement","HIVtest","Incentive","N")
    incen.table[,HIVtest:=ifelse(HIVtest==1,"HIVp","HIVtn")]
    incen.table<-cast(incen.table, Incentive+HIVtest ~ Measurement)
    incen.table$Total<-incen.table$Observed+incen.table$Unobserved
    incen.table<-as.data.table(melt(setDT(incen.table), id.vars = c("Incentive","HIVtest")))
    incen.table[,Incentive:=factor(Incentive,levels=c("None","Small","Large"))]
    incen.table<-data.table(cast(incen.table, Incentive+variable ~ HIVtest))
    rowall<-data.table("-","All",sum(incen.table[variable=="Total",HIVp],na.rm=TRUE),sum(incen.table[variable=="Total",HIVtn],na.rm=TRUE))
    names(rowall) <- names(incen.table) 
    incen.table<-rbind(incen.table,rowall)
    incen.table[,HIV.prev:=round(HIVp/(HIVp+HIVtn)*100,2)]
    colnames(incen.table)<-c("Incentive","Measurement","HIV+","HIV-","HIV+(%)")
    incen.table})

  output$table_imp <- renderTable(Bimputationres(ldatasim=data.sel()))

}

imp_shiny2 = function() {
  shinyApp(ui2, server2)
}
