
#------------------------------------------------------------
# function to estimate stage 1 alpha for ar1 structure
#------------------------------------------------------------
# `mdat`: a data frame containing the longitudinal data, with columns for the pair ID, 
#     cluster variable, time variable, and outcome variable.
# `id`: a vector of pair IDs.
# `Z`: a matrix (Y-mu)/V^1/2
# `Qinv`: the inverse of the correlation(between cluster) matrix.

estalpha1_ar1 <- function(mdat, Z, Qinv, nvisits) {
  Fa <- Fb <- 0
  for (i in mdat$cluster_id) { # for each pair
    #cveci <- unique(mdat[mdat$subject_id == i,]$cluster.var) #(1,2) 2 subejcts in each pair
    S1_j <- S2_j <- S1_ja <- S1_jb <- 0
    Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
    matZ_i1 <- matrix(Z_i1, nrow = nvisits) 
    Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
    matZ_i2 <- matrix(Z_i2, nrow = nvisits) 
    for (j in 1:2) { # for each subject in ith pair
      t_ij <- nvisits
      if (t_ij > 1) {
        for (k in 1:(t_ij-1)) {
          matZ1 <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
          matZ2 <- matrix(c(matZ_i1[k+1],matZ_i2[k+1]),nrow=2)
          S2_j <- S2_j + t(matZ1) %*% Qinv %*% matZ2
        }
        if (t_ij ==2){
          for(k in 1:t_ij){
            matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
            S1_j <- S1_j + t(matZ) %*% Qinv %*% matZ
          }
        }else{
          for(k in 1:t_ij){
            matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
            S1_ja <- S1_ja + t(matZ) %*% Qinv %*% matZ
          }
          for(k in 2:(t_ij-1)){
            matZ <- matrix(c(matZ_i1[k],matZ_i2[k]),nrow=2)
            S1_jb <- S1_jb + t(matZ) %*% Qinv %*% matZ
          }
          S1_j <- S1_ja + S1_jb
        }
      }
    }
    Fa <- Fa + S1_j
    Fb <- Fb + S2_j
  }
  ### stage 1 estimate of alpha
  alpha0 <- ( Fa - sqrt( ( Fa - 2 * Fb ) * ( Fa + 2 * Fb ) ) ) / ( 2 * Fb )
  return(alpha0)
}



#------------------------------------------------------------
# function to estimate stage 2 alpha for ar1 structure
#------------------------------------------------------------

estalpha2_ar1 <- function(alpha0){
  alpha <- as.numeric( 2 * alpha0 / ( 1 + alpha0 ^ 2 ) )
  return(alpha)
}




#------------------------------------------------------------
# function to estimate stage 1 alpha for exchangeable structure
#------------------------------------------------------------
estalpha1_exch <- function(mdat, Z, Qinv, nvisits){
  match.call()
  alphafun <- function(alpha){
    GG1 <- GG2 <- 0
    for (i in mdat$cluster_id){
      #cveci <- unique(mdat[mdat$cluster_id == i,]$cluster.var) #(1,2) 2 subejcts in each pair
      GG1j <- GG2j <- 0
      Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
      matZ_i1 <- matrix(Z_i1, nrow = nvisits) 
      Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
      matZ_i2 <- matrix(Z_i2, nrow = nvisits) 
      t_ij <- nvisits
      if(t_ij > 1){
        g1 <- vector()
        for(t in 1:t_ij){
          matZ <- matrix(c(matZ_i1[t],matZ_i2[t]),nrow=2)
          g1[t] <- t(matZ) %*% Qinv %*% matZ
        } 
        G1 <- sum(g1)
        
        g2 <- vector()
        for(t in 1:(t_ij - 1)){
          for(tt in (t+1):t_ij){
            matZ1 <- matrix(c(matZ_i1[t],matZ_i2[t]),nrow=2)
            matZ2 <- matrix(c(matZ_i1[tt],matZ_i2[tt]),nrow=2)
            g2 <- c(g2, t(matZ1) %*% Qinv %*% matZ2) 
          }
        }
        G2 <- sum(g2)
        
        denom <- ( 1 + ( t_ij - 1 ) * alpha ) ^ 2
        num1 <- alpha ^ 2 * ( t_ij - 1 ) * ( t_ij - 2 ) + 2 * alpha * ( t_ij - 1 )
        num2 <- ( 1 + alpha ^ 2 * ( t_ij - 1 ) )
      }
      GG1j <- GG1j + ( G1 * num1 ) / denom
      GG2j <- GG2j + ( G2 * num2 ) / denom
      
      GG1 <- GG1 + GG1j
      GG2 <- GG2 + GG2j
    }
    GG1 - 2 * GG2
  }
  ### stage 1 estimate of alpha
  alpha0 <- uniroot(alphafun, c(0,1), tol = 1e-10, extendInt = "yes")$root
  return(alpha0)
}


#------------------------------------------------------------
# function to estimate stage 2 alpha for exch structure
#------------------------------------------------------------

estalpha2_exch <- function(alpha0, mdat){
  
  match.call()
  
  alphapart1 <- alphapart2 <- 0
  for (i in mdat$cluster_id){
    #cveci <- unique(mdat[mdat$cluster_id == i,]$cluster.var)
    alphapart1j <- alphapart2j <- 0
    for (j in 1:2){
      t_ij <- 5#nlevels(as.factor(mdat[mdat$cluster_id == i & mdat$cluster.var == j,]$time.var))
      if(t_ij > 1){
        alphapart1num <- alpha0 * ( t_ij - 1 )* ( alpha0 * (t_ij - 2) + 2 )
        alphapart2num <- ( t_ij - 1 ) * ( 1 + alpha0 ^ 2 * (t_ij - 1) )
        alphaden <- ( 1 + alpha0 * ( t_ij - 1 ) ) ^ 2
        
        alphapart1j <- alphapart1j + alphapart1num / alphaden
        alphapart2j <- alphapart2j + alphapart2num / alphaden
      }
    }
    alphapart1 <- alphapart1 + alphapart1j
    alphapart2 <- alphapart2 + alphapart2j
  }
  alpha <- alphapart1 / alphapart2
  return(alpha)
}



#------------------------------------------------------------
# function to estimate stage 1 tau for exchangeable structure
#------------------------------------------------------------
esttau1_exch <- function(mdat, Z, Rinv,nvisits){
  Fa <- Fb <- 0
  for (i in mdat$cluster_id){
    a_1 <- a_2 <- 0
    Z_i1 <- Z[mdat$cluster_id == i & mdat$cluster.var == 1]
    matZ_i1 <- matrix(Z_i1, nrow = nvisits)
    Z_i2 <- Z[mdat$cluster_id == i & mdat$cluster.var == 2]
    matZ_i2 <- matrix(Z_i2, nrow = nvisits)
    a_1 <- a_1 + t(matZ_i1) %*% Rinv %*% matZ_i1 + t(matZ_i2) %*% Rinv %*% matZ_i2
    a_2 <- a_2 + t(matZ_i1) %*% Rinv %*% matZ_i2
    
    Fa <- Fa + a_1
    Fb <- Fb + a_2
  }
  ### stage 1 estimate of tau
  tau0 <- ( Fa - sqrt( ( Fa - 2 * Fb ) * ( Fa + 2 * Fb ) ) ) / ( 2 * Fb )
  return(tau0)
}

#------------------------------------------------------------
# function to estimate stage 2 tau for exchangeable structure
#------------------------------------------------------------
esttau2_exch <- function(tau0){
  tau <- as.numeric( 2 * tau0 / ( 1 + tau0 ^ 2 ) )
  return(tau)
}



#------------------------------------------------------------
# function for exchangeable correlaiton structure
#------------------------------------------------------------
exch_cor <- function(rho, n) {
  # Create an nxn matrix filled with tau
  cor_matrix <- matrix(rho, n, n)
  # Set the diagonal to 1
  diag(cor_matrix) <- 1
  return(cor_matrix)
}

#------------------------------------------------------------
# function for ar1 correlaiton structure
#------------------------------------------------------------

ar1_cor <- function(rho,n) {
  rho <- as.numeric(rho)
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}



#------------------------------------------------------------
# function to fit QLS
#------------------------------------------------------------
qls <- function(formula, data, subject_id, cluster_id , cluster.var, time.var, time.str){
  iter <- 0
  bdiff <- 1
  alpha0  <- 0 # initial alpha estimate
  n <- length(unique(time.var)) # assume balanced time points
  
  # use independent GEE to get initial beta estimates 
  init_mod <- geeglm(formula = formula, data = data, family = gaussian,
                     id = subject_id, waves = visit,corstr = "independence") 
  #summary(init_mod)
  beta0 <- as.vector(coef(init_mod))
  Z0 <- residuals(init_mod,"pearson") #init_mod$residuals
  
  # compute initial tau estimate
  if (time.str == "ind") {Rinv <- solve(diag(n))} # n=5 
  if (time.str == "ar1") {Rinv <- solve(ar1_cor(alpha0, n))} 
  if (time.str == "exch") {Rinv <- solve(exch_cor(alpha0, n))}
  tau0 <- esttau1_exch(mdat=data, Z = Z0, Rinv = Rinv,nvisits = n) 
  
  while(max(abs(bdiff)) > .00000001){
    
    # update beta estimates with alpha0, tau0
    Qi <- exch_cor(tau0,2) 
    if (time.str == "ind") {Ri <- diag(n)}
    if (time.str == "ar1") {Ri <- ar1_cor(alpha0, n)} 
    if (time.str == "exch") {Ri <- exch_cor(alpha0, n)}
    Fi <- kronecker(Qi,Ri)
    zcor <- fixed2Zcor(Fi, id=data$cluster_id, waves=data$order) # ？ 
    mod1 <- geeglm(formula = formula, data = matched_pair, family = gaussian,
                   id = cluster_id,corstr = "userdefined",zcor = zcor)
    #summary(mod1)
    beta1 <- as.vector(coef(mod1))
    bdiff <- beta1 - beta0
    
    # update tau0
    Z1 <- residuals(mod1,"pearson") #mod1$residuals
    if (time.str == "ind") {Rinv <- solve(diag(n))}
    if (time.str == "ar1") {Rinv <- solve(ar1_cor(alpha0, n))} 
    if (time.str == "exch") {Rinv <- solve(exch_cor(alpha0, n))}
    tau0 <- esttau1_exch(mdat=data, Z = Z1, Rinv = Rinv,nvisits = n) 
    # update alpha0 (initial alpha0 for the next iteration)
    Qinv <- solve(exch_cor(tau0, 2))
    if (time.str == "ind") {alpha0 <- alpha0}
    if (time.str == "ar1") {alpha0 <- estalpha1_ar1(mdat=data, Z=Z1, Qinv=Qinv, nvisits = n)}
    if (time.str == "exch") {alpha0 <- estalpha1_exch(mdat=data, Z=Z1, Qinv=Qinv, nvisits = n)}
    
    iter <- iter + 1
    beta0 <- beta1
    print(paste("iter:", iter, sep = " "))
    print(paste("alpha0:",alpha0, sep = " "))
    print(paste("tau0:",as.numeric(tau0), sep = " "))
    print(paste("bdiff:",max(abs(bdiff)), sep = " "))
  }
  
  # after converge, get stage 2 estimates
  tau2 <- esttau2_exch(tau0)
  if (time.str == "ind") {alpha2 <- alpha0}
  if (time.str == "ar1") {alpha2 <- estalpha2_ar1(alpha0)}
  if (time.str == "exch") {alpha2 <- estalpha2_exch(alpha0, mdat = data)}
  
  # obtain final beta estimates
  Qi <- exch_cor(tau2,2)
  if (time.str == "ind") {Ri <- diag(n)}
  if (time.str == "ar1") {Ri <- ar1_cor(alpha2, n)}
  if (time.str == "exch") {Ri <- exch_cor(alpha2, n)}
  Fi <- kronecker(Qi,Ri)
  zcor <- fixed2Zcor(Fi, id=data$cluster_id, waves=data$order) 
  mod <- geeglm(formula = formula, data = data, family = gaussian,
                id = cluster_id,corstr = "userdefined", zcor = zcor)
  beta <- as.vector(coef(mod))
  se <- summary(mod)$coefficient[,2] 
  
  fit <- list()
  fit$call <- match.call()
  fit$coefficients <- beta
  fit$se <- se
  fit$alpha <- alpha2
  fit$tau <- tau2
  fit$niter <- iter
  fit
}




#------------------------------------------------------------
### generate big cohort 
set.seed(2)
simnum <- 1000
# total sample size
N <- 250
# total number of visits
maxT <- 5
# vector of visits for all patients
visit <- rep(1:maxT, N)
# vector of patient ids
id <- rep(1:N, each = maxT)

# proportion of bav = 0 in hospital cohort
pbav0 = 0.1

# coef for baseline covariates for bav = 0 (Moderate)
b01 <- 6 # 4-->6
b02 <- -0.1
b03 <- -1.3
b04 <- 1.5
# coef for baseline covariates for bav = 1
b11 <- 6 # 4-->6
b12 <- -0.1
b13 <- -1.3
b14 <- 1.5

# coef for generating outcome values
b1_r <- 20 #intercept
b2_r <- -2 #bav
b3_r <- -0.5 #visit
b4_r <- -0.1 #age
b5_r <- -2 #sex_female
b6_r <- 2 #bsa baseline
b7_r <- 0.5 #visit:bav
beta_true <- cbind(b1_r, b2_r, b3_r, b4_r, b5_r, b6_r, b7_r)
beta_true_s <- cbind(b1_r,b3_r,b7_r)
p <- length(beta_true_s) # number of regression parameters

#Simulation----
# baseline covariates
age0i <- rnorm(pbav0 * N, mean = 60, sd = 10)
female0i <- rbinom(pbav0 * N, size = 1, prob = 0.3)
bsa_bl0i <- rnorm(pbav0 * N, mean = 2, sd = 0.2)

age1i <- rnorm((1 - pbav0) * N, mean = 60, sd = 10)
female1i <- rbinom((1 - pbav0) * N, size = 1, prob = 0.3)
bsa_bl1i <- rnorm((1 - pbav0) * N, mean = 2, sd = 0.2)

ps0_xbeta <- b01 + b02 * age0i + b03 * female0i + b04 * bsa_bl0i
pscore0 <- exp(ps0_xbeta)/ (1 + exp(ps0_xbeta))
ps1_xbeta <- b11 + b12 * age1i + b13 * female1i + b14 * bsa_bl1i
pscore1 <- exp(ps1_xbeta)/ (1 + exp(ps1_xbeta))
ps_xbeta <- c(ps0_xbeta,ps1_xbeta)
pscore <- c(pscore0, pscore1)
bavi <- rbinom(N, size = 1, prob = pscore)

# repeat baseline covariates maxT times
agei <- c(age0i, age1i)
femalei <- c(female0i, female1i)
bsa_bli <- c(bsa_bl0i, bsa_bl1i)
age <- rep(agei, each = maxT)
female <- rep(femalei, each = maxT)
bsa_bl <- rep(bsa_bli, each = maxT)
bav <- rep(bavi, each = maxT)

# generate outcome values
bi <- rnorm(N, mean = 0, sd = 5)
b <- rep(bi, each = maxT)
e <- rnorm(N*maxT, mean = 0, sd = 1)
# add the confounders (age, female, bsa_bl)
root <- b1_r + b2_r*bav +b3_r * visit + b4_r * age + b5_r * female + b6_r * bsa_bl + b7_r * bav * visit + b + e

simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))

# * Matching----
# create sample data by matching patients based on ps
simdat_base <- simdat %>% group_by(id) %>% slice(1)
ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
#summary(ps)
ps.pm <- pairmatch(ps, data = simdat_base)
#summary(ps.pm) 
matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
  filter(!is.na(matches))

K <- nrow(matched_base)

#bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
matched_long <- simdat[simdat$id %in% matched_base$id,] # sample data

# join the matches column from matched_base to matched_long
matched_long <- matched_long %>% 
  left_join(matched_base %>% select(id, matches), by = "id") %>%
  mutate(matchid = matches)

# Re-arrange the dataset with cluster id and subject id
matched_pair <- matched_long[order(matched_long$matches),] %>%
  mutate(cluster_id = rep(1:(K/2), each=(maxT*2)),
         subject_id = rep(1:K,each=maxT),
         cluster.var = rep(rep(1:2,each = maxT), (K/2)),
         order = rep(seq(1:10),(K/2))) %>%
  select(-matches)



qls(root ~  visit + bav:visit, data = matched_pair, subject_id = subject_id, 
    cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, time.str="ind")

qls(root ~  visit + bav:visit, data = matched_pair, subject_id = subject_id, 
    cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, time.str="ar1")

qls(root ~  visit + bav:visit, data = matched_pair, subject_id = subject_id, 
    cluster_id = cluster_id, cluster.var = cluster.var, time.var=visit, time.str="exch")

