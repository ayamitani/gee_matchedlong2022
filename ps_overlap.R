library(overlapping)
set.seed(2)
simnum <- 10
# total sample size
N <- 250
p <- 3 # number of regression parameters
# total number of visits
maxT <- 5
# vector of visits for all patients
visit <- rep(1:maxT, N)
# vector of patient ids
id <- rep(1:N, each = maxT)

# proportion of bav = 0 in hospital cohort
pbav0 = 0.1

overlap_ch <-rep(NA, simnum)
overlap_cm <-rep(NA, simnum)
overlap_cl <-rep(NA, simnum)
overlap_sh <-rep(NA, simnum)
overlap_sm <-rep(NA, simnum)
overlap_sl <-rep(NA, simnum)

# coef for baseline covariates (high overlap)
b1_h <- 6 # 4-->6
b2_h <- -0.1
b3_h <- -1.3
b4_h <- 1.5

# coef for baseline covariates for bav = 0 (moderate overlap)
b01_m <- 3 
b02_m <- -0.1
b03_m <- -1.5
b04_m <- 1.5
# coef for baseline covariates for bav = 1 (moderate overlap)
b11_m <- 6 
b12_m <- -0.1
b13_m <- -1.3
b14_m <- 2

# coef for baseline covariates for bav = 0 (low overlap)
b01_l <- 2
b02_l <- -0.2
b03_l <- -1.5
b04_l <- 1.5
# coef for baseline covariates for bav = 1 (low overlap)
b11_l <- 10
b12_l <- -0.1
b13_l <- -1.3
b14_l <- 1.5

# coef for generating outcome values
b1_r <- 18 #intercept
b2_r <- -0.1 #visit
b3_r <- 0.01 #age
b4_r <- -2 #sex_female
b5_r <- 8 #bsa baseline
b6_r <- 0.3 #visit:bav

for (i in 1:simnum) {
  ## High overlap:-------------------------------------------------------------
  # baseline covariates
  agei <- rnorm(N, mean = 60, sd = 10)
  femalei <- rbinom(N, size = 1, prob = 0.3)
  bsa_bli <- rnorm(N, mean = 2, sd = 0.2)
  ps_xbeta <- b1_h + b2_h * agei + b3_h * femalei + b4_h * bsa_bli
  pscore <- exp(ps_xbeta)/ (1 + exp(ps_xbeta))
  bavi <- rbinom(N, size = 1, prob = pscore)
  # repeat baseline covariates maxT times
  age <- rep(agei, each = maxT)
  female <- rep(femalei, each = maxT)
  bsa_bl <- rep(bsa_bli, each = maxT)
  bav <- rep(bavi, each = maxT)
  prop_score <- rep(ps_xbeta, each = maxT)
  
  # generate outcome values
  bi <- rnorm(N, mean = 0, sd = 5)
  b <- rep(bi, each = maxT)
  e <- rnorm(N*maxT, mean = 0, sd = 1)
  root <- b1_r + b2_r * visit + b3_r * age + b4_r * female + b5_r * bsa_bl + b6_r * bav * visit + b + e
  
  simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root, prop_score))
  simdat_base <- simdat %>% group_by(id) %>% slice(1)
  
  # calculate proportion of overlapping btw cohort p.score densities 
  ps_tav_ch <- subset(simdat_base, bav == 0, select = prop_score) # Subset the data for tav
  ps_bav_ch <- subset(simdat_base, bav == 1, select = prop_score) # Subset the data for bav
  xch <- list(ps_tav_ch = ps_tav_ch$prop_score,ps_bav_ch=ps_bav_ch$prop_score)
  overlap_ch[i] <- overlap(xch,type="2",plot=FALSE)
  
  #Matching
  ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
  ps.pm <- pairmatch(ps, data = simdat_base)
  matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
    filter(!is.na(matches))
  #bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
  
  # calculate proportion of overlapping btw sample p.score densities 
  ps_tav <- subset(matched_base, bav == 0, select = prop_score) 
  ps_bav <- subset(matched_base, bav == 1, select = prop_score) 
  xh <- list(ps_tav = ps_tav$prop_score,ps_bav=ps_bav$prop_score)
  overlap_sh[i] <- overlap(xh,type="2",plot=FALSE)
  
  # Moderate overlap:--------------------------------------------------
  # baseline covariates
  age0i <- rnorm(pbav0 * N, mean = 60, sd = 10)
  female0i <- rbinom(pbav0 * N, size = 1, prob = 0.3)
  bsa_bl0i <- rnorm(pbav0 * N, mean = 2, sd = 0.2)
  age1i <- rnorm((1 - pbav0) * N, mean = 60, sd = 10)
  female1i <- rbinom((1 - pbav0) * N, size = 1, prob = 0.3)
  bsa_bl1i <- rnorm((1 - pbav0) * N, mean = 2, sd = 0.2)
  
  ps0_xbeta <- b01_m + b02_m * age0i + b03_m * female0i + b04_m * bsa_bl0i
  pscore0 <- exp(ps0_xbeta)/ (1 + exp(ps0_xbeta))
  ps1_xbeta <- b11_m + b12_m * age1i + b13_m * female1i + b14_m * bsa_bl1i
  pscore1 <- exp(ps1_xbeta)/ (1 + exp(ps1_xbeta))
  ps_xbeta <- c(ps0_xbeta,ps1_xbeta)
  pscore <- c(pscore0, pscore1)
  bavi_m <- rbinom(N, size = 1, prob = pscore)
  
  # repeat baseline covariates maxT times
  agei <- c(age0i, age1i)
  femalei <- c(female0i, female1i)
  bsa_bli <- c(bsa_bl0i, bsa_bl1i)
  age <- rep(agei, each = maxT)
  female <- rep(femalei, each = maxT)
  bsa_bl <- rep(bsa_bli, each = maxT)
  bav <- rep(bavi_m, each = maxT)
  prop_score <- rep(ps_xbeta, each = maxT)
  
  # generate outcome values
  bi <- rnorm(N, mean = 0, sd = 5)
  b <- rep(bi, each = maxT)
  e <- rnorm(N*maxT, mean = 0, sd = 1)
  root <- b1_r + b2_r * visit + b3_r * age + b4_r * female + b5_r * bsa_bl + b6_r * bav * visit + b + e
  
  simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root,prop_score))
  simdat_base <- simdat %>% group_by(id) %>% slice(1)
  
  # calculate proportion of overlapping btw cohort p.score densities 
  ps_tav_cm <- subset(simdat_base, bav == 0, select = prop_score) # Subset the data for tav
  ps_bav_cm <- subset(simdat_base, bav == 1, select = prop_score) # Subset the data for bav
  xcm <- list(ps_tav_cm = ps_tav_cm$prop_score,ps_bav_cm=ps_bav_cm$prop_score)
  overlap_cm[i] <- overlap(xcm,type="2",plot=FALSE)
  
  # Matching
  ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
  ps.pm <- pairmatch(ps, data = simdat_base)
  matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
    filter(!is.na(matches))
  #bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
  
  # calculate proportion of overlapping btw sample p.score densities   
  ps_tav <- subset(matched_base, bav == 0, select = prop_score) 
  ps_bav <- subset(matched_base, bav == 1, select = prop_score) 
  xm <- list(ps_tav = ps_tav$prop_score,ps_bav=ps_bav$prop_score)
  overlap_sm[i] <- overlap(xm,type="2",plot=FALSE)
  
  # Low overlap:-------------------------------------------------------------
  # baseline covariates
  age0i <- rnorm(pbav0 * N, mean = 60, sd = 10)
  female0i <- rbinom(pbav0 * N, size = 1, prob = 0.3)
  bsa_bl0i <- rnorm(pbav0 * N, mean = 2, sd = 0.2)
  age1i <- rnorm((1 - pbav0) * N, mean = 60, sd = 10)
  female1i <- rbinom((1 - pbav0) * N, size = 1, prob = 0.3)
  bsa_bl1i <- rnorm((1 - pbav0) * N, mean = 2, sd = 0.2)
  
  ps0_xbeta <- b01_l + b02_l * age0i + b03_l * female0i + b04_l * bsa_bl0i
  pscore0 <- exp(ps0_xbeta)/ (1 + exp(ps0_xbeta))
  ps1_xbeta <- b11_l + b12_l * age1i + b13_l * female1i + b14_l * bsa_bl1i
  pscore1 <- exp(ps1_xbeta)/ (1 + exp(ps1_xbeta))
  ps_xbeta <- c(ps0_xbeta,ps1_xbeta)
  pscore <- c(pscore0, pscore1)
  bavi_l <- rbinom(N, size = 1, prob = pscore)
  
  # repeat baseline covariates maxT times
  agei <- c(age0i, age1i)
  femalei <- c(female0i, female1i)
  bsa_bli <- c(bsa_bl0i, bsa_bl1i)
  age <- rep(agei, each = maxT)
  female <- rep(femalei, each = maxT)
  bsa_bl <- rep(bsa_bli, each = maxT)
  bav <- rep(bavi_l, each = maxT)
  prop_score <- rep(ps_xbeta, each = maxT)
  
  # generate outcome values
  bi <- rnorm(N, mean = 0, sd = 5)
  b <- rep(bi, each = maxT)
  e <- rnorm(N*maxT, mean = 0, sd = 1)
  root <- b1_r + b2_r * visit + b3_r * age + b4_r * female + b5_r * bsa_bl + b6_r * bav * visit + b + e
  
  simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root,prop_score))
  simdat_base <- simdat %>% group_by(id) %>% slice(1)
  
  # calculate proportion of overlapping btw cohort p.score densities 
  ps_tav_cl <- subset(simdat_base, bav == 0, select = prop_score) # Subset the data for tav
  ps_bav_cl <- subset(simdat_base, bav == 1, select = prop_score) # Subset the data for bav
  xcl <- list(ps_tav_cl = ps_tav_cl$prop_score,ps_bav_cl=ps_bav_cl$prop_score)
  overlap_cl[i] <- overlap(xcl,type="2",plot=FALSE)
  
  # Matching
  ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
  ps.pm <- pairmatch(ps, data = simdat_base)
  matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
    filter(!is.na(matches))
  #bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
  
  # calculate proportion of overlapping btw sample p.score densities 
  ps_tav <- subset(matched_base, bav == 0, select = prop_score) # Subset the data for tav
  ps_bav <- subset(matched_base, bav == 1, select = prop_score) # Subset the data for bav
  xl <- list(ps_tav = ps_tav$prop_score,ps_bav=ps_bav$prop_score)
  overlap_sl[i] <- overlap(xl,type="2",plot=FALSE)
  
  print(i)
}
# type=2 returns the proportion of overlapping area between two densities
prop <- c("High","Moderate","Minimal")
cohort <- round(c(mean(unlist(overlap_ch)),mean(unlist(overlap_cm)),mean(unlist(overlap_cl))),3)
sample <- round(c(mean(unlist(overlap_sh)),mean(unlist(overlap_sm)),mean(unlist(overlap_sl))),3)
(overlap <- as.data.frame(cbind(prop, cohort, sample)))

