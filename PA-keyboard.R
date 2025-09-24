######################################################################
#------------------Part 1: Function for PA-keyboard------------------#
######################################################################
#Load the dependent R packages.
library(purrr)
library(R2jags)
library(tidyverse)
library(parallel)
library(doParallel)

#The sub-function called within the PA_keyboard_function.------------------------
{## Calculate the area under curve (AUC) of plasma concentration.
cal_p <- function(CL,dose){
  AUC.data<- dose/CL
  return(AUC.data)
}
## Performing posterior sampling of PD parameters in the PD model.
generate <- function(toxicity,cohortsize,AUC_data,bug.model.path){
  data <- list (toxicity=toxicity,
                cohortsize=cohortsize,
                AUC_data=AUC_data
  )
  
  pd.sim <- jags(data,model.file = bug.model.path,
                 parameters.to.save = c("beta"),
                 n.chains = 3, n.iter = 10000,DIC=TRUE,progress.bar="none",quiet =TRUE
  )
  
  beta_mean <- pd.sim[["BUGSoutput"]][["mean"]]$beta
  beta <-pd.sim[["BUGSoutput"]][["sims.matrix"]][,1]
  result <- list(
    result=beta,
    mean=beta_mean
  )
  return(result)
}
## Calculate the overlap coefficient (OVL) between two distributions.
similarity_ov <- function(x, y, from = NULL, to = NULL, n = 4096, bw = c("max","mean","x","y")) {
  stopifnot(is.numeric(x), is.numeric(y))
  
  if (is.null(from)) from <- min(min(x), min(y))
  if (is.null(to))   to   <- max(max(x), max(y))
  
  d1 <- density(x, from = from, to = to, n = n)
  d2 <- density(y, from = from, to = to, n = n)
  
  bw <- match.arg(bw)
  if (bw %in% c("max","mean")) {
    bw_use <- if (bw == "max") max(d1$bw, d2$bw) else mean(c(d1$bw, d2$bw))
    d1 <- density(x, from = from, to = to, n = n, bw = bw_use)
    d2 <- density(y, from = from, to = to, n = n, bw = bw_use)
  } else if (bw == "x") {
    d2 <- density(y, from = from, to = to, n = n, bw = d1$bw)
  } else if (bw == "y") {
    d1 <- density(x, from = from, to = to, n = n, bw = d2$bw)
  }
  
  dx  <- d1$x[2] - d1$x[1]
  ovl <- sum(pmin(d1$y, d2$y)) * dx   
  ovl <- max(min(ovl, 1), 0)         
  return(ovl)
}
## Splitting the keys in the keyboard method according to the target toxicity rate.
getkey <- function(target,epi1,epi2){
  marginL <- target-epi1  
  marginU <- target+epi2
  delta=marginU-marginL
  lkey=NULL; rkey=NULL
  i=0; cutoff=0.3
  while(cutoff>0)
  {
    i=i+1
    cutoff = marginL-i*delta
    lkey = c(cutoff, lkey)
  }
  lkey[lkey<0]=0
  
  i=0; cutoff=0.3
  while(cutoff<1)
  {
    i=i+1
    cutoff = marginU+i*delta
    rkey = c(rkey, cutoff)
  }
  rkey[rkey>1]=1
  key=c(lkey, marginL, marginU, rkey)
  
  return(key)
}
## Find the key with the maximum posterior probability in the Keyboard method.
find_d <- function(keys,n,y,n_borrow=0,y_borrow=0,delta=0,target=0.30,epi1=0.05,alpha=0.05,beta=0.05){
  nkeys <- length(keys)-1
  marginL <- target-epi1
  p <- vector("numeric",length=nkeys)
  descision <- vector("integer",length=1)
  if (delta == 0) {
    for (i in 1:length(keys)-1){
      p[i] <- pbeta(keys[i+1],y+alpha,n-y+beta)-pbeta(keys[i],y+alpha,n-y+beta)
    }
  }else{
    for (i in 1:length(keys)-1){
      p[i] <- pbeta(keys[i+1],y+delta*y_borrow+alpha,n-y+delta*(n_borrow-y_borrow)+beta)-pbeta(keys[i],y+delta*y_borrow+alpha,n-y+delta*(n_borrow-y_borrow)+beta)
    }
  }
  
  N <- which(keys==marginL)
  if (which.max(p)> N){
    descision <- -1
  }else if (which.max(p)== N){
    descision <- 0
  }else{
    descision <- 1
  }
  return(descision)
}
## PAVA Algorithm: Adjustment for monotonicity of toxicity-dose relationship during posterior estimation.
pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1)
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol)))
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}
}
#Notion: seed--Set the number of seeds.
#        p.true.a--The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
#        p.true--The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
#        pk--PK parameters for adult virtual populations in simulation studies.
#        weightbsa--The drug clearance rate in pediatric trial.
#        true_dose_a--Dose setting for adult virtual populations in simulation studies.
#        target--Target toxicity rate.
#        ncohort--Maximum enrolment cohort size.
#        cohortsize_a--Number of patients enrolled in each cohort for the adult trial.
#        cohortsize_p--Number of patients enrolled in each cohort for the Pediatrics trial.
#        n.earlystop--Early termination occurs when the maximum number of patients accumulating at a given dose is reached.
#        startdose--Starting dose for dose escalation.
#        epi1,epi2--Together they indicate the width of the target key.
#        overtoxicity--The threshold for acceptable toxicity levels, exceeding which indicates overdose.
#        AUC_LOW,AUC_HIGH--The range for dose matching between adult and paediatric trials based on AUC.
#        c_A_adv--Number of cohorts enrolled in the adult trial ahead of the paediatric trial.
#        bug.model.path--The storage path for the bug model to be invoked during method execution.
PA_keyboard_function <- function(seed,p.true.a,p.true,pk,weightbsa,true_dose_a,target=0.30,ncohort = 10,
                                 cohortsize_a = 3,cohortsize_p = 3,n.earlystop = 15,startdose = 1,
                                 epi1=0.05,epi2=0.05,overtoxicity=0.95,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2,
                                 bug.model.path)
{
  set.seed(123+seed)
  keys <- getkey(target=target,epi1=epi1,epi2=epi2)
  nkeys <-  length(keys)-1
  d  <-  startdose
  d_a <- startdose
  marginL <- target-epi1
  marginU <- target+epi2
  ndose_p  <-  length(true_dose_a)
  ndose_a  <-  length(true_dose_a)
  dose_p <- 1:ndose_p
  dose_p_2 <- 1:ndose_p   
  dose_a <- 1:ndose_a
  dose_a_2 <- 1:ndose_a
  x <- matrix(rep(0, (1+cohortsize_p) * ncohort), ncol = cohortsize_p+1)
  x_a <- matrix(rep(0,(1+cohortsize_a) * ncohort), ncol = cohortsize_a+1)
  data <- matrix(rep(0, (1*cohortsize_p+1) * ncohort), ncol = 1*cohortsize_p+1)
  data_a <- matrix(rep(0, (1*cohortsize_a+1) * ncohort), ncol =1*cohortsize_a+1)
  pi_a <- matrix(rep(0, (2*cohortsize_a+1) * ncohort), ncol = 2*cohortsize_p+1)
  pi_p <- matrix(rep(0, (2*cohortsize_p+1) * ncohort), ncol = 2*cohortsize_p+1)
  CL_p <-  vector("numeric",length = cohortsize_p*ncohort)
  acutal_dose_p <- vector("numeric",length = cohortsize_p*ncohort)
  age_p <- vector("integer",length = cohortsize_p*ncohort)
  y <- rep(0, ndose_p)
  n <- rep(0, ndose_p)
  y_a <- rep(0, ndose_a)
  n_a <- rep(0, ndose_a)
  AUC_a <- vector("list",length = ndose_a)
  AUC_p <- vector("list",length = ndose_p)
  d_s_store <- vector("list",length = ndose_p)
  
  n_earlystop_p <- vector("integer",length = 1)
  n_stop <- vector("integer",length = 1)
  n_earlystop_a <- vector("integer",length = 1)
  data_delta<- list(0,0,0,0,0)
  data_similarity <- list(0,0,0,0,0)
  data_beta_p <- list(0,0,0,0,0)
  data_beta_a <- list(0,0,0,0,0)
  alpha <- 0.05
  beta <- 0.05
  delta <- 0
  
  for (j in 1:c_A_adv) {
    data_a[j,] <- c(d_a,rlnorm(cohortsize_a, meanlog = pk$a_mu_CL, sdlog =pk$a_sigma_CL))
    pi_a[j,] <- c(d_a,p.true.a[d_a],p.true.a[d_a],p.true.a[d_a],cal_p(data_a[j,2],true_dose_a[d_a]),cal_p(data_a[j,3],true_dose_a[d_a]),cal_p(data_a[j,4],true_dose_a[d_a]))
    x_a[j,1] <- d_a
    x_a[j,2:4] <- rbinom(3,1,p.true.a[d_a])
    
    y_a[d_a] <- y_a[d_a]+sum(x_a[j,2:4])
    n_a[d_a]  <-  n_a[d_a] + cohortsize_a
    if (n_a[d_a]>=2*cohortsize_a){
      p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
      if (p_overtoxicity>overtoxicity) {
        dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
        if(is_empty(dose_a_2)){
          n_earlystop_a <- 1
          break
        }
        if(length(dose_a_2)==0){
          n_earlystop_a <- 1
          break
        }
      }
    }
    d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
    d_a <- d_a+d_a_current
    if (d_a == 0) d_a <- 1
    if (!identical(dose_a_2,dose_a)) {
      if (d_a > which.max(dose_a_2)){
        d_a <- which.max(dose_a_2)
      }
    }else{
      if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
    }
  }
  for (i in 1:ncohort) {
    a_cohort <- ncohort
    if (n_earlystop_a == 1) {
      a_cohort <- c_A_adv
    }
    if (n_earlystop_a ==0 & n_stop == 0){
      if (i<= a_cohort-c_A_adv){
        data_a[i+c_A_adv,] <- c(d_a,rlnorm(cohortsize_a, meanlog = pk$a_mu_CL, sdlog =pk$a_sigma_CL))
        pi_a[i+c_A_adv,] <- c(d_a,p.true.a[d_a],p.true.a[d_a],p.true.a[d_a],cal_p(data_a[i+c_A_adv,2],true_dose_a[d_a]),cal_p(data_a[i+c_A_adv,3],true_dose_a[d_a]),cal_p(data_a[i+c_A_adv,4],true_dose_a[d_a]))
        x_a[i+c_A_adv,1] <- d_a
        x_a[i+c_A_adv,2:4] <- rbinom(3,1,p.true.a[d_a])
        
        y_a[d_a] <- y_a[d_a]+sum(x_a[i+c_A_adv,2:4])
        n_a[d_a]  <-  n_a[d_a] + cohortsize_a
        
        if (n_a[d_a] >= n.earlystop) {n_stop <- 1}
        else{
          if (n_a[d_a]>=2*cohortsize_a){
            p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
            if (p_overtoxicity > overtoxicity) {
              dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
              if(is_empty(dose_a_2)){
                n_earlystop_a <- 1
              }
              if(length(dose_a_2)==0){
                n_earlystop_a <- 1
              }
            }
          }
          d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
          d_a <- d_a+d_a_current
          if (n_earlystop_a == 0){
            if (d_a == 0) {d_a <- 1}
            if (!identical(dose_a_2,dose_a)) {
              if (d_a > which.max(dose_a_2)){
                d_a <- which.max(dose_a_2)
              }
            }else{
              if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
            }
          }
          
        }
      }}
    age_p[(i-1)*3+1]<-floor(runif(1,min=10,max=15))
    age_p[(i-1)*3+2]<-floor(runif(1,min=10,max=15))
    age_p[(i-1)*3+3]<-floor(runif(1,min=10,max=15))
    CL_p[(i-1)*3+1] <- weightbsa %>% filter(Age==age_p[(i-1)*3+1]) %>% select(pk_cl_p) %>% c() %>% unname() %>% unlist()
    CL_p[(i-1)*3+2] <- weightbsa %>% filter(Age==age_p[(i-1)*3+2]) %>% select(pk_cl_p) %>% c() %>% unname() %>% unlist()
    CL_p[(i-1)*3+3] <- weightbsa %>% filter(Age==age_p[(i-1)*3+3]) %>% select(pk_cl_p) %>% c() %>% unname() %>% unlist()
    data[i,] <- c(d,rlnorm(1, meanlog = log(CL_p[(i-1)*3+1]), sdlog =pk$p_sigma_CL),rlnorm(1, meanlog = log(CL_p[(i-1)*3+2]), sdlog =pk$p_sigma_CL),rlnorm(1, meanlog = log(CL_p[(i-1)*3+3]), sdlog =pk$p_sigma_CL))
    acutal_dose_p[(i-1)*3+1] <- weightbsa %>% filter(Age==age_p[(i-1)*3+1]) %>% select(all_of(d)) %>% c() %>% unname() %>% unlist()
    acutal_dose_p[(i-1)*3+2] <- weightbsa %>% filter(Age==age_p[(i-1)*3+2]) %>% select(all_of(d)) %>% c() %>% unname() %>% unlist()
    acutal_dose_p[(i-1)*3+3] <- weightbsa %>% filter(Age==age_p[(i-1)*3+3]) %>% select(all_of(d)) %>% c() %>% unname() %>% unlist()
    pi_p[i,] <- c(d,p.true[d],p.true[d],p.true[d],cal_p(data[i,2],acutal_dose_p[(i-1)*3+1]),cal_p(data[i,3],acutal_dose_p[(i-1)*3+2]),cal_p(data[i,4],acutal_dose_p[(i-1)*3+3]))
    x[i,1] <- d
    x[i,2:4] <- rbinom(3,1,p.true[d])
    
    y[d]  <-  y[d] + sum(x[i,2:4])
    n[d]  <-  n[d] + cohortsize_p
    
    AUC_p[[1]] <- as.data.frame(pi_p) %>% dplyr::filter(V1==1)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric() 
    AUC_p[[2]] <- as.data.frame(pi_p) %>% dplyr::filter(V1==2)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric() 
    AUC_p[[3]] <- as.data.frame(pi_p) %>% dplyr::filter(V1==3)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric() 
    AUC_p[[4]] <- as.data.frame(pi_p) %>% dplyr::filter(V1==4)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric() 
    AUC_p[[5]] <- as.data.frame(pi_p) %>% dplyr::filter(V1==5)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric() 
    
    AUC_a[[1]] <- as.data.frame(pi_a) %>% dplyr::filter(V1==1)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric()
    AUC_a[[2]] <- as.data.frame(pi_a) %>% dplyr::filter(V1==2)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric()
    AUC_a[[3]] <- as.data.frame(pi_a) %>% dplyr::filter(V1==3)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric()
    AUC_a[[4]] <- as.data.frame(pi_a) %>% dplyr::filter(V1==4)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric()
    AUC_a[[5]] <- as.data.frame(pi_a) %>% dplyr::filter(V1==5)%>% dplyr::select(V5,V6,V7) %>% unlist() %>% as.numeric()
    
    if (n[d] >= n.earlystop){break}

    if (n[d]> 3 ) {
      mean_auc_p <- mean(AUC_p[[d]])
      apply_sort <- c(mean(AUC_a[[1]]),mean(AUC_a[[2]]),mean(AUC_a[[3]]),mean(AUC_a[[4]]),mean(AUC_a[[5]]))
      d_s <-  sort(abs(apply_sort - mean_auc_p), index.return = T)[['ix']][1]
      d_s_store[[d]] <- c(d_s_store[[d]],d_s)
    if (apply_sort[d_s] >AUC_LOW*mean_auc_p & apply_sort[d_s]< AUC_HIGH*mean_auc_p & n_a[d_s] > 0){
        toxicity <- as.data.frame(x_a) %>% dplyr::filter(V1==d_s)%>% dplyr::select(V2,V3,V4) %>% unlist() %>% as.numeric()
        toxicity.p <- as.data.frame(x) %>% dplyr::filter(V1==d) %>% dplyr::select(V2,V3,V4) %>% unlist() %>% as.numeric()
        cycle_a <- n_a[d_s]
        cycle_p <- n[d]
        beta1_a <- generate(toxicity = toxicity,cohortsize=cycle_a,AUC_data=AUC_a[[d_s]],bug.model.path=bug.model.path)
        beta1_p <- generate(toxicity = toxicity.p,cohortsize=cycle_p,AUC_data=AUC_p[[d]],bug.model.path=bug.model.path)
        generate_p <- beta1_p[["result"]]
        generate_a <- beta1_a[["result"]]
        data_beta_p[[d]] <- c(data_beta_p[[d]],beta1_p[["mean"]])
        data_beta_a[[d_s]]<- c(data_beta_a[[d_s]],beta1_a[["mean"]])
        delta <- round(similarity_ov(generate_p,generate_a),2)
        # delta <- round(elastic(T=T),2)
        data_delta[[d]] <- c(data_delta[[d]],delta)
        data_similarity[[d]] <- c(data_similarity[[d]],T)
        if (delta > 0.01){
          d_current <- find_d(keys=keys,n=n[d],y=y[d],n_borrow=n_a[d_s],y_borrow=y_a[d_s],delta=delta,target=target,epi1=0.05,beta = beta,alpha = alpha)
        }else{
          d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
          delta  <-  0
        }
      }else{
        d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
        delta  <-  0 
      }
      
    } else{
      d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
    }
    
    data_delta[[d]] <- c(data_delta[[d]],delta)
    if (n[d]>=2*cohortsize_p){
      if (delta==0){
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)
        if (p_overtoxicity > overtoxicity) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
      }else{
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)
        p_overtoxicity_b <- 1-pbeta(target,y[d]+delta*y_a[d_s]+alpha,n[d]-y[d]+delta*(n_a[d_s]-y_a[d_s])+beta)
        if (p_overtoxicity > overtoxicity | p_overtoxicity_b > overtoxicity) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
        
      }
    }
    d <- d+d_current
    if (d == 0 ) {d <- 1}
    if (!identical(dose_p_2,dose_p)) {
      if (d > which.max(dose_p_2)){
        d <- which.max(dose_p_2)
      }
    }else{
      if (d > which.max(dose_p)) {d <- which.max(dose_p)}
    }
    delta <- 0
  }
  ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
  phat = (y + alpha)/(n+alpha+beta)
  phat.var = (y +alpha) * (n - y +beta)/((n +alpha+beta)^2 * (n +alpha+beta + 1))
  
  ## perform the isotonic transformation using PAVA
  phat1 = pava(phat, wt = 1/phat.var)
  phat = phat1 + (1:length(phat1)) * 1e-10  ## break ties by adding an increasingly small number
  
  selectd = sort(abs(phat - target), index.return = T)[['ix']][1] ## select dose closest to the target as the MTD
  dose_select <- 1:ndose_p
  selectdose = dose_select[selectd]
  
  phat_b <- vector("numeric",length = ndose_p)
  phat.var_b <- vector("numeric",length = ndose_p)
  final_delta <- vector("numeric",length = ndose_p)
  d_s_final <- vector("numeric",length = ndose_p)
  for(i in 1:5){
    mean_auc_p <- mean(AUC_p[[i]])
    apply_sort_a <- c(mean(AUC_a[[1]]),mean(AUC_a[[2]]),mean(AUC_a[[3]]),mean(AUC_a[[4]]),mean(AUC_a[[5]]))
    d_s_final[i]  <-  sort(abs(apply_sort_a - mean_auc_p), index.return = T)[['ix']][1]
    if(!is.na(mean_auc_p)){
      if (apply_sort_a[d_s_final[i]] >AUC_LOW*mean_auc_p & apply_sort_a[d_s_final[i]]< AUC_HIGH*mean_auc_p){
        d_s_final[i] <- d_s_final[i]
      }else{
        d_s_final[i] <- NA
      }
      if (!is.na(d_s_final[i])){
        if ( apply_sort_a[d_s_final[i]] >AUC_LOW*mean_auc_p & apply_sort_a[d_s_final[i]]< AUC_HIGH*mean_auc_p &n[i]> 3 ){
          toxicity <- as.data.frame(x_a) %>% dplyr::filter(V1==d_s_final[i])%>% dplyr::select(V2,V3,V4) %>% unlist() %>% as.numeric()
          toxicity.p <- as.data.frame(x) %>% dplyr::filter(V1==i) %>% dplyr::select(V2,V3,V4) %>% unlist() %>% as.numeric()
          cycle_a <- n_a[d_s_final[i]]
          cycle_p <- n[i]
          beta1_a <- generate(toxicity = toxicity,cohortsize=cycle_a,AUC_data=AUC_a[[d_s_final[i]]],bug.model.path=bug.model.path)
          beta1_p <- generate(toxicity = toxicity.p,cohortsize=cycle_p,AUC_data=AUC_p[[i]],bug.model.path=bug.model.path)
          generate_p <- beta1_p[["result"]]
          generate_a <- beta1_a[["result"]]
          final_delta[i] <- round(similarity_ov(generate_p,generate_a),2)
          
        }else{
          final_delta[i] <- 0
        }}else{
          final_delta[i] <- 0
        }
    }
  }

  d_s_final[is.na(d_s_final)] <- 0
  for (i in 1:5){
    if (d_s_final[i]!=0){
      phat_b[i] = (y[i] +y_a[d_s_final[i]]*final_delta[i]+alpha)/(n[i]+n_a[d_s_final[i]]*final_delta[i]+alpha+beta)
      phat.var_b[i] = (y[i] +y_a[d_s_final[i]]*final_delta[i]+alpha) * (n[i]+n_a[d_s_final[i]]*final_delta[i] - y[i] -y_a[d_s_final[i]]*final_delta[i]+beta)/((n[i]+n_a[d_s_final[i]]*final_delta[i]+alpha+beta)^2 * (n[i]+n_a[d_s_final[i]]*final_delta[i]+alpha+beta + 1))
    }else{
      phat_b[i] = (y[i] +alpha)/(n[i]+alpha+beta)
      phat.var_b[i] = (y[i] +alpha) * (n[i]-y[i]+beta)/((n[i]+alpha+beta)^2 * (n[i]+alpha+beta+ 1))
    }
  }
  
  phat_b_1 = pava(phat_b, wt = 1/phat.var_b)
  phat_b = phat_b_1 + (1:length(phat_b_1)) * 1e-10  ## break ties by adding an increasingly small number
  
  selectd_b = sort(abs(phat_b - target), index.return = T)[['ix']][1] ## select dose closest to the target as the MTD
  dose_select_b <- 1:ndose_p
  selectdose_b = dose_select_b[selectd_b]
  
  N <- sum(n)
  Num_MTD <- which(p.true==target)
  if (n_earlystop_p==1){
    selectdose  <-  0  
    selectdose_b  <- 0 
  }
  if (Num_MTD==5){
    overdose <- 0
  }else{
    if(sum(n[(Num_MTD+1):length(dose_select)])> 0.6*sum(n)){
      overdose <- 1 
    }else{
      overdose <- 0
    }
  }
  
  data_final <- data.frame(
    n_earlystop_a = n_earlystop_a,
    n_earlystop_p = n_earlystop_p,
    selectdose=selectdose,
    selectdose_b=selectdose_b,
    N=N,
    overdose=overdose
  )
  sim <- list(
    data=data_final,
    data_delta=data_delta,
    data_similarity=data_similarity,
    data_beta_a=data_beta_a,
    data_beta_p=data_beta_p,
    N=n,
    d_s_final=d_s_final,
    phat=phat,
    phat_b=phat_b
  )
  return(sim)
}
######################################################################
#--------------Part 2: Function for Comparative method---------------#
######################################################################
#The original keyboard method without borrowing information.
Original_keyboard_function <- function(seed,p.true,p.true.a,target=0.30,ncohort = 10,
                                       cohortsize = 3,n.earlystop = 15,startdose = 1,epi1=0.05,epi2=0.05,
                                       overtoxicity=0.95)
{
  set.seed(123+seed)
  keys <- getkey(target=target,epi1=epi1,epi2=epi2)
  d  <-  startdose
  d_a <- startdose
  marginL <- target-epi1
  marginU <- target+epi2
  ndose_p = length(p.true)
  ndose_a = length(p.true.a)
  dose_p <- 1:ndose_p
  dose_p_2 <- 1:ndose_p   
  dose_a <- 1:ndose_a
  dose_a_2 <- 1:ndose_a
  x <- matrix(rep(0, cohortsize * ncohort), ncol = cohortsize)
  x_a <- matrix(rep(0, cohortsize * ncohort), ncol = cohortsize)
  data <- matrix(rep(0, (1*cohortsize+1) * ncohort), ncol = 1*cohortsize+1)
  data_a <- matrix(rep(0, (1*cohortsize+1) * ncohort), ncol =1*cohortsize+1)
  y <- rep(0, ndose_p)
  n <- rep(0, ndose_p)
  y_a <- rep(0, ndose_a)
  n_a <- rep(0, ndose_a) 
  n_earlystop_p <- vector("integer",length = 1)
  n_stop <- vector("integer",length = 1)
  n_earlystop_a <- vector("integer",length = 1)
  data_delta<- list(0,0,0,0,0)
  alpha <- 0.05
  beta <- 0.05
  delta <- 0
  
  for (j in 1:2) {
    x_a[j,] <- rbinom(cohortsize,1,p.true.a[d_a])
    data_a[j,] <- c(x_a[j,],d_a)
    
    y_a[d_a] <- y_a[d_a]+sum(x_a[j,])
    n_a[d_a]  <-  n_a[d_a] + cohortsize
    if (n_a[d_a]>= 2*cohortsize ){
      p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
      if (p_overtoxicity>overtoxicity) {
        dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
        if(is_empty(dose_a_2)){
          n_earlystop_a <- 1
          break
        }
        if(length(dose_a_2)==0){
          n_earlystop_a <- 1
          break
        }
      }
    }
    d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
    d_a <- d_a+d_a_current  
    if (d_a == 0) d_a <- 1
    if (!identical(dose_a_2,dose_a)) {
      if (d_a > which.max(dose_a_2)){
        d_a <- which.max(dose_a_2)
      }
    }else{
      if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
    }
  }
  for (i in 1:ncohort) {
    a_cohort <- ncohort
    if (n_earlystop_a == 1) {
      a_cohort <- 2
    }
    if (n_earlystop_a ==0 & n_stop == 0){
      if (i<= a_cohort-2){    
        x_a[i+2,] <- rbinom(cohortsize,1,p.true.a[d_a])
        data_a[i+2,] <- c(x_a[i+2,],d_a)
        
        y_a[d_a] <- y_a[d_a]+sum(x_a[i+2,])
        n_a[d_a]  <-  n_a[d_a] + cohortsize
        if (n_a[d_a] >= n.earlystop) {n_stop <- 1}
        else{
          if (n_a[d_a]>=2*cohortsize){
            p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
            if (p_overtoxicity>overtoxicity) {
              dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
              if(is_empty(dose_a_2)){
                n_earlystop_a <- 1
              }
              if(length(dose_a_2)==0){
                n_earlystop_a <- 1
              }
            }
          }
          d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
          d_a <- d_a+d_a_current  
          if (n_earlystop_a == 0){
            if (d_a == 0) {d_a <- 1}
            if (!identical(dose_a_2,dose_a)) {
              if (d_a > which.max(dose_a_2)){
                d_a <- which.max(dose_a_2)
              }
            }else{
              if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
            }
          }
        }
      }}
    x[i,] <- rbinom(cohortsize,1,p.true[d])
    data[i,] <- c(x[i,],d)
    
    y[d]  <-  y[d] + sum(x[i,])
    n[d]  <-  n[d] + cohortsize
    if (n[d] >= n.earlystop){break}
    d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
    if (n[d]>=2*cohortsize){
      if (delta==0){
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)
        if (p_overtoxicity > overtoxicity) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
      }else{
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)    
        if (p_overtoxicity > overtoxicity ) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
        
      }
    }
    d <- d+d_current
    if (d == 0 ) {d <- 1}
    if (!identical(dose_p_2,dose_p)) {
      if (d > which.max(dose_p_2)){
        d <- which.max(dose_p_2)
      }
    }else{
      if (d > which.max(dose_p)) {d <- which.max(dose_p)}
    }
    
  }
  ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
  phat = (y + alpha)/(n+alpha+beta)
  phat.var = (y +alpha+beta) * (n - y +alpha+beta)/((n +alpha+beta)^2 * (n +alpha+beta + 1))
  
  ## perform the isotonic transformation using PAVA
  phat1 = pava(phat, wt = 1/phat.var)
  phat = phat1 + (1:length(phat1)) * 1e-10  ## break ties by adding an increasingly small number
  
  selectd = sort(abs(phat - target), index.return = T)[['ix']][1] ## select dose closest to the target as the MTD
  dose_select <- 1:ndose_p
  selectdose = dose_select[selectd]
  
  
  
  N <- sum(n)
  if (is_empty(n[which(p.true==target)])){
    selectd_nomatch = sort(abs(p.true - target), index.return = T)[['ix']][1]
    N_MTD <- n[selectd_nomatch] 
  }
  
  if (n_earlystop_p==1){
    selectdose  <-  0  
  }
  
  Num_MTD <- which(p.true==target)
  if (Num_MTD==5){
    overdose <- 0
  }else{
    if(sum(n[(Num_MTD+1):length(dose_select)])>0.6*sum(n)){
      overdose <- 1 
    }else{
      overdose <- 0
    }
  }
  
  data <- data.frame(
    n_earlystop_a = n_earlystop_a,
    n_earlystop_p = n_earlystop_p,
    selectdose=selectdose,
    overdose=overdose,
    N=N
  )
  sim <- list(
    data=data,
    N=n,
    phat=phat
  )
  return(sim)
}
#p-Keyboard method, which borrows information only based on fixed doses after dose BSA conversion.
p_keyboard_function <- function(seed,p.true,p.true.a,true_dose_a,true_dose_p,target=0.30,ncohort = 10,
                                cohortsize_p = 3,cohortsize_a = 3,n.earlystop = 15,startdose = 1,
                                epi1=0.05,epi2=0.05,overtoxicity=0.95)
{
  set.seed(123+seed)
  keys <- getkey(target=target,epi1=epi1,epi2=epi2)
  nkeys <-  length(keys)-1
  d  <-  startdose
  d_a <- startdose
  marginL <- target-epi1
  marginU <- target+epi2
  ndose_p = length(p.true)
  ndose_a = length(p.true.a)
  dose_p <- 1:ndose_p
  dose_p_2 <- 1:ndose_p   
  dose_a <- 1:ndose_a
  dose_a_2 <- 1:ndose_a
  x <- matrix(rep(0, cohortsize_p * ncohort), ncol = cohortsize_p)
  x_a <- matrix(rep(0, cohortsize_a * ncohort), ncol = cohortsize_a)
  data <- matrix(rep(0, (cohortsize_p+1) * ncohort), ncol = cohortsize_p+1)
  data_a <- matrix(rep(0, (cohortsize_a+1) * ncohort), ncol =cohortsize_a+1)
  y <- rep(0, ndose_p)
  n <- rep(0, ndose_p)
  y_a <- rep(0, ndose_a)
  n_a <- rep(0, ndose_a)
  n_earlystop_p <- vector("integer",length = 1)
  n_stop <- vector("integer",length = 1)
  n_earlystop_a <- vector("integer",length = 1)
  data_delta<- list(0,0,0,0,0)
  data_similarity <- list(0,0,0,0,0)
  alpha <- 0.05
  beta <- 0.05
  delta <- 0
  similarity <- function (x1,x2){
    kldiv <- function(x1, x2, nbreaks = 100, minx = min(c(x1, x2)),
                      maxx = max(c(x1, x2)), small = 0.001){
      dy <- (maxx-minx)/nbreaks
      breaks <- seq(minx, maxx, by = dy)
      x1d <- hist(x1, plot = F, breaks = breaks)
      x2d <- hist(x2, plot = F, breaks = breaks)
      x1p <- (x1d$counts + small) / sum(x1d$counts + small)
      x2p <- (x2d$counts + small) / sum(x2d$counts + small)
      vals <- x1p * log2(x1p/x2p)
      kd <- sum(na.omit(vals))
      
      hout <- list(kd = kd,
                   histograms = list(h1 = x1d, h2 = x2d))
      class(hout) <- "kldiverg"
      return(hout)
    }
    distance1 <- kldiv(x1, x2, nbreaks = 100, minx = min(c(x1, x2)),maxx = max(c(x1, x2)), small = 0.01)[["kd"]]
    distance2 <-kldiv(x1, x2, nbreaks = 100, minx = min(c(x1, x2)),maxx = max(c(x1, x2)), small = 0.01)[["kd"]]
    distance <-0.5*(distance1+distance2)
    return(distance)
  }
  elastic <- function(a=1,b=1,T){
    delta <- 1/(1+exp(1)^(a+b*log(T)))
  }
  for (j in 1:2) {
    x_a[j,] <- rbinom(cohortsize_a,1,p.true.a[d_a])
    data_a[j,] <- c(x_a[j,],d_a)
    y_a[d_a] <- y_a[d_a]+sum(x_a[j,])
    n_a[d_a]  <-  n_a[d_a] + cohortsize_a
    if (n_a[d_a]>=2*cohortsize_a){
      p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
      if (p_overtoxicity>overtoxicity) {
        dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
        if(is_empty(dose_a_2)){
          n_earlystop_a <- 1
          break
        }
        if(length(dose_a_2)==0){
          n_earlystop_a <- 1
          break
        }
      }
    }
    d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
    d_a <- d_a+d_a_current
    if (d_a == 0) d_a <- 1
    if (!identical(dose_a_2,dose_a)) {
      if (d_a > which.max(dose_a_2)){
        d_a <- which.max(dose_a_2)
      }
    }else{
      if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
    }
  }
  for (i in 1:ncohort) {
    a_cohort <- ncohort
    if (n_earlystop_a == 1) {
      a_cohort <- 2
    }
    if (n_earlystop_a ==0 & n_stop == 0){
      if (i<= a_cohort-2){
        x_a[i+2,] <- rbinom(cohortsize_a,1,p.true.a[d_a])
        data_a[i+2,] <- c(x_a[i+2,],d_a)
        
        y_a[d_a] <- y_a[d_a]+sum(x_a[i+2,])
        n_a[d_a]  <-  n_a[d_a] + cohortsize_a
        if (n_a[d_a] >= n.earlystop) {n_stop <- 1}
        else{
          if (n_a[d_a]>=2*cohortsize_a){
            p_overtoxicity <- 1-pbeta(target,y_a[d_a]+alpha,n_a[d_a]-y_a[d_a]+beta)
            if (p_overtoxicity>overtoxicity) {
              dose_a_2 <- dose_a_2[-d_a:-length(dose_a_2)]
              if(is_empty(dose_a_2)){
                n_earlystop_a <- 1
              }
              if(length(dose_a_2)==0){
                n_earlystop_a <- 1
              }
            }
          }
          d_a_current <- find_d(keys=keys,n=n_a[d_a],y=y_a[d_a],target=target,epi1=0.05,beta = beta,alpha = alpha)
          d_a <- d_a+d_a_current
          if (n_earlystop_a == 0){
            if (d_a == 0) {d_a <- 1}
            if (!identical(dose_a_2,dose_a)) {
              if (d_a > which.max(dose_a_2)){
                d_a <- which.max(dose_a_2)
              }
            }else{
              if (d_a > which.max(dose_a)) {d_a <- which.max(dose_a)}
            }
          }
          
        }
      }}
    x[i,] <- rbinom(cohortsize_p,1,p.true[d])
    data[i,] <- c(x[i,],d)
    
    y[d]  <-  y[d] + sum(x[i,])
    n[d]  <-  n[d] + cohortsize_p
    if (n[d] >= n.earlystop){break}
    if (n[d]> 1*cohortsize_p) {
      if (n_a[d] > 0){
        generate_a <- rbeta(100000,y_a[d]+alpha,n_a[d]-y_a[d]+beta)
        generate_p <- rbeta(100000,y[d]+alpha,n[d]-y[d]+beta)
        T <- similarity(generate_p,generate_a)
        delta <- round(elastic(T=T),2)
        data_delta[[d]] <- c(data_delta[[d]],delta)
        data_similarity[[d]] <- c(data_similarity[[d]],T)
        if (delta > 0.01){
          d_current <- find_d(keys=keys,n=n[d],y=y[d],n_borrow=n_a[d],y_borrow=y_a[d],delta=delta,target=target,epi1=0.05,beta = beta,alpha = alpha)
        }else{
          d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
          delta  <-  0
        }
      }
      
    } else{
      d_current <- find_d(keys=keys,n=n[d],y=y[d],target=target,epi1=0.05,beta = beta,alpha = alpha)
    }
    if (n[d]>=2*cohortsize_p){
      if (delta==0){
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)
        if (p_overtoxicity > overtoxicity) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
      }else{
        p_overtoxicity <- 1-pbeta(target,y[d]+alpha,n[d]-y[d]+beta)
        p_overtoxicity_b <- 1-pbeta(target,y[d]+delta*y_a[d]+alpha,n[d]-y[d]+delta*(n_a[d]-y_a[d])+beta)
        if (p_overtoxicity > overtoxicity | p_overtoxicity_b > overtoxicity) {
          dose_p_2 <- dose_p_2[-d:-length(dose_p_2)]
          if(is_empty(dose_p_2)){
            n_earlystop_p <- 1
            break
          }
          if(length(dose_p_2)==0 ){
            n_earlystop_p <- 1
            break
          }
        }
        
      }
    }
    d <- d+d_current
    if (d == 0 ) {d <- 1}
    if (!identical(dose_p_2,dose_p)) {
      if (d > which.max(dose_p_2)){
        d <- which.max(dose_p_2)
      }
    }else{
      if (d > which.max(dose_p)) {d <- which.max(dose_p)}
    }
    delta <- 0
  }
  phat = (y + alpha)/(n+alpha+beta)
  phat.var = (y +alpha) * (n - y +beta)/((n +alpha+beta)^2 * (n +alpha+beta + 1))
  
  phat1 = pava(phat, wt = 1/phat.var)
  phat = phat1 + (1:length(phat1)) * 1e-10  
  
  selectd = sort(abs(phat - target), index.return = T)[['ix']][1] 
  dose_select <- 1:ndose_p
  selectdose = dose_select[selectd]

  phat_b <- vector("numeric",length = ndose_p)
  phat.var_b <- vector("numeric",length = ndose_p)
  final_delta <- vector("numeric",length = ndose_p)
  for(i in 1:5){
    if (n[i]> 3 & n_a[i]> 0){
      generate_a <- rbeta(100000,y_a[d]+alpha,n_a[d]-y_a[d]+beta)
      generate_p <- rbeta(100000,y[d]+alpha,n[d]-y[d]+beta)
      T <- similarity(generate_p,generate_a)
      final_delta[i] <- round(elastic(T=T),2)
    }else{
      final_delta[i] <- 0
    }
  }

  for (i in 1:5){
    phat_b[i] = (y[i] +y_a[i]*final_delta[i]+alpha)/(n[i]+n_a[i]*final_delta[i]+alpha+beta)
    phat.var_b[i] = (y[i] +y_a[i]*final_delta[i]+alpha) * (n[i]+n_a[i]*final_delta[i] - y[i] -y_a[i]*final_delta[i]+beta)/((n[i]+n_a[i]*final_delta[i]+alpha+beta)^2 * (n[i]+n_a[i]*final_delta[i]+alpha+beta + 1))
  }
  
  ## perform the isotonic transformation using PAVA
  phat_b_1 = pava(phat_b, wt = 1/phat.var_b)
  phat_b = phat_b_1 + (1:length(phat_b_1)) * 1e-10  ## break ties by adding an increasingly small number
  
  selectd_b = sort(abs(phat_b - target), index.return = T)[['ix']][1] ## select dose closest to the target as the MTD
  dose_select_b <- 1:ndose_p
  selectdose_b = dose_select_b[selectd_b]
  
  N <- sum(n)
  Num_MTD <- which(p.true==target)
  if (n_earlystop_p==1){
    selectdose  <-  0  
    selectdose_b  <- 0 
  }
  if (Num_MTD==5){
    overdose <- 0
  }else{
    if(sum(n[(Num_MTD+1):length(dose_select)])> 0.6*sum(n)){
      overdose <- 1 
    }else{
      overdose <- 0
    }
  }
  
  data <- data.frame(
    n_earlystop_a = n_earlystop_a,
    n_earlystop_p = n_earlystop_p,
    selectdose=selectdose,
    selectdose_b=selectdose_b,
    N=N,
    overdose=overdose
  )
  sim <- list(
    data=data,
    data_delta=data_delta,
    data_similarity=data_similarity,
    N=n,
    N_a=n_a,
    phat=phat,
    phat_b=phat_b
  )
  return(sim)
}

######################################################################
#-----------------Part 3: Simulation study code----------------------#
######################################################################
#Load the dependent R packages.
library(Keyboard)
library(beepr)
library(readxl)

#The drug clearance rate in pediatric trial.
weightbsa <- read_excel("~/PA-keyboard/bodyweightbsa.xlsx")
weightbsa2 <- weightbsa %>% mutate(pk_cl_a=pk$a_mu_CL,p_sigma_CL=pk$p_sigma_CL) %>% mutate(pk_cl_p=round(exp(pk_cl_a)*(Weight/70)^0.75,2))
# PK parameters for adult virtual populations in simulation studies.
pk <- data.frame(
  a_mu_CL=2.5,
  a_sigma_CL=0.01,
  p_sigma_CL=0.01
)
# Dose setting for adult virtual populations in simulation studies.
true_dose_a <- c(0.05,0.085,0.11,0.14,0.17)
# Number of computer simulations
times <- 5000

#Scenarios 1
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.2,0.3,0.45,0.6)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.45,0.6)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_1 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                                  weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                                  n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_1 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                    target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                    n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_1 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                    target=0.3,ncohort = 10,cohortsize = 3,
                                    n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})


#Scenarios 2
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.2,0.25,0.3)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.05,0.1,0.3,0.4,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_2 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_2 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_2 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})


#Scenarios 3
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.3,0.4,0.5)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.3,0.4,0.5,0.6,0.7)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_3 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_3 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})



#Scenarios 4
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.05,0.1,0.3,0.35)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.35,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_4 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_4 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 5
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.3,0.4,0.5)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.15,0.3,0.5,0.6,0.65)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_5 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_5 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_5 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 6
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.05,0.1,0.2,0.3)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.05,0.1,0.2,0.3,0.45)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_6 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_6 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_6 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})


#Scenarios 7
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.05,0.08,0.2,0.3,0.4)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.4,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_7 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_7 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_7 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})


#Scenarios 8
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.15,0.3,0.45,0.55)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.15,0.3,0.45,0.55,0.6)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  results_PA_keyboard_8 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=2)
  results_p_keyboard_8 <- parallel::parLapply(cl,seed,p_keyboard_function,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                              target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                              n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  results_Original_keyboard_8 <- parallel::parLapply(cl,seed,Original_keyboard_function,p.true.a=p.true.a,p.true=p.true,
                                                     target=0.3,ncohort = 10,cohortsize = 3,
                                                     n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})


######################################################################
#--Part 4:Calculation of Evaluation Indicators of Simulation study---#
######################################################################
output_fun <- function(results_PA_keyboard,results_p_keyboard,results_Original_keyboard){
  PofMTD_1_nob <- (results_Original_keyboard %>% map_dbl(~.x$data[[3]]) %>% table()/times)[c("1","2","3","4","5")]*100 %>% round(digits = 2)
  N_1_nob <- apply( do.call(rbind, lapply(results_Original_keyboard, `[[`, "N")),2,mean)%>% round(digits = 2)
  PofMTD_1_tox <- (results_p_keyboard %>% map_dbl(~.x$data[[4]]) %>% table()/times)[c("1","2","3","4","5")]*100 %>% round(digits = 2)
  N_1_tox <- apply( do.call(rbind, lapply(results_p_keyboard, `[[`, "N")),2,mean)%>% round(digits = 2)
  PofMTD_1_PA <- (results_PA_keyboard %>% map_dbl(~.x$data[[4]]) %>% table()/times)[c("1","2","3","4","5")]*100 %>% round(digits = 2)
  N_1_PA <- apply( do.call(rbind, lapply(results_PA_keyboard, `[[`, "N")),2,mean)%>% round(digits = 2)

  output1_PofMTD_N <- rbind(PofMTD_1_nob,N_1_nob,
                            PofMTD_1_tox,N_1_tox,
                            PofMTD_1_PA,N_1_PA)
  
  Pofearlystop_nob <- results_Original_keyboard %>% map_dbl(~.x$data[[2]]) %>% mean*100 %>% round(digits = 2)
  Pofoverdose_nob <- results_Original_keyboard %>% map_dbl(~.x$data[[4]]) %>% mean*100 %>% round(digits = 2)
  Pofearlystop_tox <- results_p_keyboard %>% map_dbl(~.x$data[[2]]) %>% mean*100 %>% round(digits = 2)
  Pofoverdose_tox <- results_p_keyboard %>% map_dbl(~.x$data[[6]]) %>% mean*100 %>% round(digits = 2)
  Pofearlystop_PA <- results_PA_keyboard %>% map_dbl(~.x$data[[2]]) %>% mean*100 %>% round(digits = 2)
  Pofoverdose_PA <- results_PA_keyboard %>% map_dbl(~.x$data[[6]]) %>% mean*100 %>% round(digits = 2)

  ALLN_nob <- results_Original_keyboard %>% map_dbl(~.x$data[["N"]]) %>% mean
  ALLN_tox <- results_p_keyboard %>% map_dbl(~.x$data[["N"]]) %>% mean
  ALLN_PA <- results_PA_keyboard %>% map_dbl(~.x$data[["N"]]) %>% mean
  
  output2_early_over <- rbind(c(ALLN_nob,Pofearlystop_nob,Pofoverdose_nob),
                              c(ALLN_tox,Pofearlystop_tox,Pofoverdose_tox),
                              c(ALLN_PA,Pofearlystop_PA,Pofoverdose_PA))
  colnames(output2_early_over) <- c("ALL_N","earlystop","overdose")
  rownames(output2_early_over) <- c("Original-keyboard","p-keyboard","PA-keyboard")
  
  result_1 <- list(MTD_N=output1_PofMTD_N,early_over=output2_early_over)
  return(result_1)
}
result_1<- output_fun(results_Original_keyboard_1,results_p_keyboard_1,results_PA_keyboard_1)
result_2<- output_fun(results_Original_keyboard_2,results_p_keyboard_2,results_PA_keyboard_2)
result_3<- output_fun(results_Original_keyboard_3,results_p_keyboard_3,results_PA_keyboard_3)
result_4<- output_fun(results_Original_keyboard_4,results_p_keyboard_4,results_PA_keyboard_4)
result_5<- output_fun(results_Original_keyboard_5,results_p_keyboard_5,results_PA_keyboard_5)
result_6<- output_fun(results_Original_keyboard_6,results_p_keyboard_6,results_PA_keyboard_6)
result_7<- output_fun(results_Original_keyboard_7,results_p_keyboard_7,results_PA_keyboard_7)
result_8<- output_fun(results_Original_keyboard_8,results_p_keyboard_8,results_PA_keyboard_8)
MTD_N_ALL <- rbind(result_1$MTD_N,result_2$MTD_N,result_3$MTD_N,result_4$MTD_N,
                   result_5$MTD_N,result_6$MTD_N,result_7$MTD_N,result_8$MTD_N
) %>% as.data.frame()
early_over_ALL <- rbind(result_1$early_over,result_2$early_over,result_3$early_over,result_4$early_over,
                        result_5$early_over,result_6$early_over,result_7$early_over,result_8$early_over
)%>% as.data.frame()


######################################################################
#-------------------Part 5: Sensitivity Analysis---------------------#
######################################################################
#Load the dependent R packages.
library(Keyboard)
library(beepr)
library(readxl)

#The drug clearance rate in pediatric trial.
weightbsa <- read_excel("~/PA-keyboard/bodyweightbsa.xlsx")
weightbsa2 <- weightbsa %>% mutate(pk_cl_a=pk$a_mu_CL,p_sigma_CL=pk$p_sigma_CL) %>% mutate(pk_cl_p=round(exp(pk_cl_a)*(Weight/70)^0.75,2))
# PK parameters for adult virtual populations in simulation studies.
pk <- data.frame(
  a_mu_CL=2.5,
  a_sigma_CL=0.01,
  p_sigma_CL=0.01
)
# Dose setting for adult virtual populations in simulation studies.
true_dose_a <- c(0.05,0.085,0.11,0.14,0.17)
# Number of computer simulations
times <- 5000


#Scenarios 1
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.2,0.3,0.45,0.6)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.45,0.6)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_1_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                                 weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                                 n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_1_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                           weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                           n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_1_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                           weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                           n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 2
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.2,0.25,0.3)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.05,0.1,0.3,0.4,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_2_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_2_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_2_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 3
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.3,0.4,0.5)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.3,0.4,0.5,0.6,0.7)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_3_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_3_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_3_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 4
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.05,0.1,0.3,0.35)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.35,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_4_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_4_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_4_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 5
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.1,0.15,0.3,0.4,0.5)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.15,0.3,0.5,0.6,0.65)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_5_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_5_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_5_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 6
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.05,0.1,0.2,0.3)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.05,0.1,0.2,0.3,0.45)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_6_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_6_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_6_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 7
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.05,0.08,0.2,0.3,0.4)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.1,0.2,0.3,0.4,0.5)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_7_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_7_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_7_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

#Scenarios 8
#The actual toxicity rate of the adult virtual population at each dose level in the simulation study.
p.true.a <- c(0.01,0.15,0.3,0.45,0.55)
#The actual toxicity rate of the Pediatrics virtual population at each dose level in the simulation study.
p.true <- c(0.15,0.3,0.45,0.55,0.6)

system.time({
  cl<- makeCluster(30)
  clusterSetRNGStream(cl, iseed =5)
  clusterEvalQ(cl,library(purrr))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(R2jags))
  clusterExport(cl,c('p.true.a','pk','p.true',"true_dose_a","weightbsa2"))
  clusterExport(cl,c("pava","generate","find_d","similarity_ov","cal_p","getkey"))
  seed <- 1:times
  data_sensi_8_0.8_1.25 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                               weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                               n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.8,AUC_HIGH=1.25,c_A_adv=2)
  data_sensi_8_c3 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=3)
  data_sensi_8_c4 <- parallel::parLapply(cl,seed,PA_keyboard_function,pk=pk,p.true.a=p.true.a,p.true=p.true,true_dose_a=true_dose_a,
                                         weightbsa=weightbsa2,target=0.3,ncohort = 10,cohortsize_a = 3,cohortsize_p = 3,
                                         n.earlystop = 100,startdose = 1,epi1=0.05,epi2=0.05,AUC_LOW=0.9,AUC_HIGH=1.1,c_A_adv=4)
  stopCluster(cl)
  beep(sound = 8, expr = NULL)
})

output_fun_sensi <- function(M1_data,data_S1,data_S2.1,data_S2.2){
  calulation_sensi<- function(M1_data){
    PofMTD_M<- (M1_data %>% map_dbl(~.x$data[[4]]) %>% table()/times)[c("1","2","3","4","5")]*100 %>% round(digits = 2)
    N_M <- apply( do.call(rbind, lapply(M1_data, `[[`, "N")),2,mean)%>% round(digits = 2)
    Pofearlystop_M <- M1_data %>% map_dbl(~.x$data[[2]]) %>% mean*100 %>% round(digits = 2)
    Pofoverdose_M <- M1_data %>% map_dbl(~.x$data[[6]]) %>% mean*100 %>% round(digits = 2)
    ALLN_M <- M1_data %>% map_dbl(~.x$data[["N"]]) %>% mean
    return(list(PofMTD=PofMTD_M,
                N=N_M,
                Pofearlystop=Pofearlystop_M,
                Pofoverdose=Pofoverdose_M,
                ALLN=ALLN_M))
  }
  outp_M <- calulation_sensi(M1_data)
  outp_S1 <- calulation_sensi(data_S1)
  outp_S2.1 <- calulation_sensi(data_S2.1)
  outp_S2.2 <- calulation_sensi(data_S2.2)
  out_name <- c("Main analysis","Sensitive analysis 1","Sensitive analysis 2.1","Sensitive analysis 2.2")
  PofMTD <- cbind(out_name, rbind(outp_M$PofMTD,outp_S1$PofMTD,outp_S2.1$PofMTD,outp_S2.2$PofMTD))
  N <- cbind(out_name, rbind(outp_M$N,outp_S1$N,outp_S2.1$N,outp_S2.2$N))
  Pofearlystop <- cbind(out_name, rbind(outp_M$Pofearlystop,outp_S1$Pofearlystop,outp_S2.1$Pofearlystop,outp_S2.2$Pofearlystop))
  Pofoverdose <- cbind(out_name, rbind(outp_M$Pofoverdose,outp_S1$Pofoverdose,outp_S2.1$Pofoverdose,outp_S2.2$Pofoverdose))
  ALLN <- cbind(out_name, rbind(outp_M$ALLN,outp_S1$ALLN,outp_S2.1$ALLN,outp_S2.2$ALLN))
  return(list(PofMTD=PofMTD,
              N=N,
              Pofearlystop=Pofearlystop,
              Pofoverdose=Pofoverdose,
              ALLN=ALLN))
}
result_1_sensi <- output_fun_sensi(results_PA_keyboard_1,data_sensi_1_0.8_1.25,data_sensi_1_c3,data_sensi_1_c4)
result_2_sensi <- output_fun_sensi(results_PA_keyboard_2,data_sensi_2_0.8_1.25,data_sensi_2_c3,data_sensi_2_c4)
result_3_sensi <- output_fun_sensi(results_PA_keyboard_3,data_sensi_3_0.8_1.25,data_sensi_3_c3,data_sensi_3_c4)
result_4_sensi <- output_fun_sensi(results_PA_keyboard_4,data_sensi_4_0.8_1.25,data_sensi_4_c3,data_sensi_4_c4)
result_5_sensi <- output_fun_sensi(results_PA_keyboard_5,data_sensi_5_0.8_1.25,data_sensi_5_c3,data_sensi_5_c4)
result_6_sensi <- output_fun_sensi(results_PA_keyboard_6,data_sensi_6_0.8_1.25,data_sensi_6_c3,data_sensi_6_c4)
result_7_sensi <- output_fun_sensi(results_PA_keyboard_7,data_sensi_7_0.8_1.25,data_sensi_7_c3,data_sensi_7_c4)
result_8_sensi <- output_fun_sensi(results_PA_keyboard_8,data_sensi_8_0.8_1.25,data_sensi_8_c3,data_sensi_8_c4)
Scenarios_name <- paste("Scenarios",1:8)
result_sensi_PofMTD <- cbind(rep(Scenarios_name,each=4),rbind(result_1_sensi$PofMTD,result_2_sensi$PofMTD,
                                                              result_3_sensi$PofMTD,result_4_sensi$PofMTD,
                                                              result_5_sensi$PofMTD,result_6_sensi$PofMTD,
                                                              result_7_sensi$PofMTD,result_8_sensi$PofMTD))
result_sensi_N <- cbind(rep(Scenarios_name,each=4),rbind(result_1_sensi$N,result_2_sensi$N,
                                                         result_3_sensi$N,result_4_sensi$N,
                                                         result_5_sensi$N,result_6_sensi$N,
                                                         result_7_sensi$N,result_8_sensi$N))
result_sensi_Pofearlystop <- cbind(rep(Scenarios_name,each=4),rbind(result_1_sensi$Pofearlystop,result_2_sensi$Pofearlystop,
                                                                    result_3_sensi$Pofearlystop,result_4_sensi$Pofearlystop,
                                                                    result_5_sensi$Pofearlystop,result_6_sensi$Pofearlystop,
                                                                    result_7_sensi$Pofearlystop,result_8_sensi$Pofearlystop))
result_sensi_Pofoverdose <- cbind(rep(Scenarios_name,each=4),rbind(result_1_sensi$Pofoverdose,result_2_sensi$Pofoverdose,
                                                                   result_3_sensi$Pofoverdose,result_4_sensi$Pofoverdose,
                                                                   result_5_sensi$Pofoverdose,result_6_sensi$Pofoverdose,
                                                                   result_7_sensi$Pofoverdose,result_8_sensi$Pofoverdose))
result_sensi_ALLN <- cbind(rep(Scenarios_name,each=4),rbind(result_1_sensi$ALLN,result_2_sensi$ALLN,
                                                            result_3_sensi$ALLN,result_4_sensi$ALLN,
                                                            result_5_sensi$ALLN,result_6_sensi$ALLN,
                                                            result_7_sensi$ALLN,result_8_sensi$ALLN))



