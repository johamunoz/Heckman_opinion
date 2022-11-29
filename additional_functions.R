
library(lme4)
library(data.table)
library(fMultivar)
library(mice)

library(ggplot2)
library(dplyr)
library(shiny)
library(reshape)

# Prevalence functions ----
logit <- function(x){log(x/(1-x))}
inv_logit<- function(x){exp(x)/(1+exp(x))}

#Function for calculating prevalence CI
f.prevalence <- function(data,outcome_name) {
  a<-summary(data[,get(outcome_name)])
  b<-prop.test(x = a[[2]], n = sum(a), correct = TRUE)
  prop<-b$estimate[[1]]
  ci<-b$conf
  se<-(ci[2]-ci[1])/(2*1.96)
  #Replace 0's by 0.000001
  prevalence<-ifelse(prop==0,1e-10,ifelse(prop==1,1-1e-10,prop))
  return(list(prev= prevalence,prev_se=se))
}

#Function to calculate pooled prevalence across imputed datasets 
prevalence.pool<-function(data, outcome_name) {
  m <-length(unique(data$.imp))
  #Calculate absolute risk of outcome per study
  prev_group<-setDT(data)[, f.prevalence(data=.SD, outcome_name=outcome_name), by = list(.imp)]
  
  #Logit transformation
  prev_group[,prev_logit:=logit(prev)] 
  prev_group[,prev_se_logit:=prev_se/(prev*(1-prev))]
  
  # Pool cluster estimates
  pool_est_logit <- mean(prev_group$prev_logit)# pool est
  w_var_logit <- mean(prev_group$prev_se_logit^2) #within var
  b_var_logit <- var(prev_group$prev_logit) #between var
  pool_se_logit <- sqrt(w_var_logit + (1 + (1/m)) * b_var_logit) #pool se
  r <- (1 + 1 / m) * (b_var_logit / w_var_logit)
  v <- (m - 1) * (1 + (1/r))^2
  t <- qt(0.975, v) # t critical value for 95% CI
  pool_LC_logit=pool_est_logit-pool_se_logit*t
  pool_UC_logit=pool_est_logit+pool_se_logit*t
  
  prevalence<-inv_logit(pool_est_logit)*100
  ci.lb<-inv_logit(pool_LC_logit)*100
  ci.ub<-inv_logit(pool_UC_logit)*100
  return(data.frame(cbind(prevalence=prevalence,ci.lb=ci.lb,ci.ub=ci.ub)))
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

Imputationres<-function(datasim){
  
  dataSim0 <- datasim[, c("marital","condom","highhiv","IDreffect","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim0$ybin<-as.factor(datasim$ybin)
  ini <- mice(dataSim0, maxit = 0)
  meth0<-ini$method
  pred0 <- ini$pred
  meth0[c("ybin")]<-"2l.heckman"
  pred0["ybin","IDreffect"] <- -3
  
  dataSimNERV <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSimNERV$ybin<-as.factor(datasim$ybin)
  ini <- mice(dataSimNERV, maxit = 0)
  methNERV<-ini$method
  predNERV <- ini$pred
  methNERV[c("ybin")]<-"2l.heckman"
  
  
  dataSim1 <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim1$ex1<-ifelse(datasim$incen1!=0,1,0)
  dataSim1$ybin<-as.factor(datasim$ybin1)
  ini <- mice(dataSim1, maxit = 0)
  meth1<-ini$method
  pred1 <- ini$pred
  meth1[c("ybin")]<-"2l.heckman"
  pred1["ybin","ex1"] <- -3
  
  dataSim0i <- datasim[, c("marital","condom","highhiv","IDreffect","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim0i$ex1<-ifelse(datasim$incen1!=0,1,0)
  dataSim0i$ybin<-as.factor(datasim$ybin1)
  ini <- mice(dataSim0i, maxit = 0)
  meth0i<-ini$method
  pred0i <- ini$pred
  meth0i[c("ybin")]<-"2l.heckman"
  pred0i["ybin","IDreffect"] <- -3
  pred0i["ybin","ex1"] <- -3
  
  
  dataSim12 <- datasim[, c("marital","condom","highhiv","sage1","sage2","sage3","sage4","sage5","sage6","sage7","sage8","sage9")]
  dataSim12$incen<-as.factor(datasim$incen)
  dataSim12$ybin<-as.factor(datasim$ybin12)
  ini <- mice(dataSim12, maxit = 0)
  meth12<-ini$method
  pred12 <- ini$pred
  meth12[c("ybin")]<-"2l.heckman"
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
                                          maxit=1,
                                          meta_method="reml"))
  
  imp_heck_iid <- suppressWarnings(mice( dataSim0, # dataset with missing values
                                         m = 15,   # number of imputations
                                         seed = 1234, #seed attached to the dataID
                                         meth = meth0, #imputation method vector
                                         pred = pred0, #imputation predictors matrix
                                         maxit=1,
                                         meta_method="reml"))
  
  imp_heck_inc1 <- suppressWarnings(mice( dataSim1, # dataset with missing values
                                          m = 15,   # number of imputations
                                          seed = 1234, #seed attached to the dataID
                                          meth = meth1, #imputation method vector
                                          pred = pred1, #imputation predictors matrix
                                          maxit=1,
                                          meta_method="reml"))
  
  imp_heck_iidi <- suppressWarnings(mice( dataSim0i, # dataset with missing values
                                          m = 15,   # number of imputations
                                          seed = 1234, #seed attached to the dataID
                                          meth = meth0i, #imputation method vector
                                          pred = pred0i, #imputation predictors matrix
                                          maxit=1,
                                          meta_method="reml"))
  
  imp_heck_inc12 <- suppressWarnings(mice( dataSim12, # dataset with missing values
                                           m = 15,   # number of imputations
                                           seed = 1234, #seed attached to the dataID
                                           meth = meth12, #imputation method vector
                                           pred = pred12, #imputation predictors matrix
                                           maxit=1,
                                           meta_method="reml"))
  
  
  mar<-prevalence.pool(data=complete( imp_mar, "long"),outcome_name<-"ybin")
  heckNERV<-prevalence.pool(data=complete( imp_heck_NERV, "long"),outcome_name<-"ybin")
  heckiid<-prevalence.pool(data=complete( imp_heck_iid, "long"),outcome_name<-"ybin")
  heckinc1<-prevalence.pool(data=complete( imp_heck_inc1, "long"),outcome_name<-"ybin")
  heckiidi<-prevalence.pool(data=complete( imp_heck_iidi, "long"),outcome_name<-"ybin")
  heckinc12<-prevalence.pool(data=complete( imp_heck_inc12, "long"),outcome_name<-"ybin")
  
  perc<-c(sum(is.na(dataSim0$ybin))/nrow(dataSim0),
          sum(is.na(dataSimNERV$ybin))/nrow(dataSimNERV),
          sum(is.na(dataSim0$ybin))/nrow(dataSim0),
          sum(is.na(dataSim1$ybin))/nrow(dataSim1),
          sum(is.na(dataSim0i$ybin))/nrow(dataSim0i),
          sum(is.na(dataSim12$ybin))/nrow(dataSim12))
  
  tableimp<-rbind(mar, heckNERV,heckiid,heckinc1,heckiidi,heckinc12)[,c("prevalence","ci.lb","ci.ub")]
  Imputation<-c("MAR", "Heckman no ERV","Heckman interviewer ID","Heckman small incentive","Heckman interviewer ID+incentive","Heckman staggered incentive")
  tablefin<-as.data.table(cbind(Imputation, tableimp))
  tablefin[,bias:=round(prevalence-mean(datasim$ybin0)*100,2)]
  tablefin[,perc:=round(perc*100,2)]
  tablefin[,ci.lb:=round(ci.lb,2)]
  tablefin[,ci.ub:=round(ci.ub,2)]
  tablefin<-tablefin[,c("Imputation","prevalence","ci.lb","ci.ub","perc","bias")]
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
  rho = ifelse(rho == -1,-0.99,ifelse(rho==1,0.99,rho)) # it has to be slightly modified due to simulation of skewed errors
  n= nrow(datasim)
  d <- diag(2)
  d[2,1] <- d[1,2] <- rho
  
  e_norm <- mvtnorm::rmvnorm(n = n, mean = rep(0,2), sigma = d)
  e_skew <- fMultivar::rmvst(n = n, dim = 2, mu = c(0,0), Omega = d, alpha = c(-2,1), df = 8)
  
  datasim_norm <- data_error(datasim, error=e_norm, prop, incentv1,incentv2)
  datasim_skew <- data_error(datasim, error=e_skew, prop, incentv1,incentv2)
  
  ldatasim<-list(datasim_norm,datasim_skew)
  ldatasim[[1]]
}

