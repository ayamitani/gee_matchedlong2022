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

# coef for baseline covariates
b1 <- 6 # 4-->6
b2 <- -0.1
b3 <- -1.3
b4 <- 1.5

# coef for generating outcome values
b1_r <- 30
b2_r <- 0.5
b3_r <- -0.2
b4_r <- 0.5
beta_true <- cbind(b1_r, b3_r, b4_r)

# create lists for output(cohort)
estbeta_ind <- vector("list", length = simnum)
estbeta_exch <- vector("list", length = simnum)
estbeta_ar1 <- vector("list", length = simnum)
estbeta_unstr <- vector("list", length = simnum)
se_ind <- vector("list", length = simnum)
se_exch <- vector("list", length = simnum)
se_ar1 <- vector("list", length = simnum)
se_unstr <- vector("list", length = simnum)

relbias_ind <- vector("list", length = simnum) # Relative bias
relbias_exch <- vector("list", length = simnum)
relbias_ar1 <- vector("list", length = simnum)
relbias_unstr <- vector("list", length = simnum)
mse_ind <- vector("list", length = simnum) # MSE
mse_exch <- vector("list", length = simnum)
mse_ar1 <- vector("list", length = simnum)
mse_unstr <- vector("list", length = simnum)
cp_ind <- vector("list", length = simnum) # list for for calculating coverage probability
cp_exch <- vector("list", length = simnum)
cp_ar1 <- vector("list", length = simnum)
cp_unstr <- vector("list", length = simnum)
outvec_ind <- vector("list", length = simnum)
outvec_exch <- vector("list", length = simnum)
outvec_ar1 <- vector("list", length = simnum)
outvec_unstr <- vector("list", length = simnum)

# create lists for output(sample)
estbeta_inds <- vector("list", length = simnum)
estbeta_exchs <- vector("list", length = simnum)
estbeta_ar1s <- vector("list", length = simnum)
estbeta_unstrs <- vector("list", length = simnum)
se_inds <- vector("list", length = simnum)
se_exchs <- vector("list", length = simnum)
se_ar1s <- vector("list", length = simnum)
se_unstrs <- vector("list", length = simnum)

relbias_inds <- vector("list", length = simnum) # Relative bias
relbias_exchs <- vector("list", length = simnum)
relbias_ar1s <- vector("list", length = simnum)
relbias_unstrs <- vector("list", length = simnum)
mse_inds <- vector("list", length = simnum) # MSE
mse_exchs <- vector("list", length = simnum)
mse_ar1s <- vector("list", length = simnum)
mse_unstrs <- vector("list", length = simnum)
cp_inds <- vector("list", length = simnum) # list for for calculating coverage probability
cp_exchs <- vector("list", length = simnum)
cp_ar1s <- vector("list", length = simnum)
cp_unstrs <- vector("list", length = simnum)
outvec_inds <- vector("list", length = simnum)
outvec_exchs <- vector("list", length = simnum)
outvec_ar1s <- vector("list", length = simnum)
outvec_unstrs <- vector("list", length = simnum)

geeglm_ci <- function(model, level = 0.95) {
  # Calculate the lower and upper bounds of the confidence interval for each parameter
  lower <- coef(model) - qnorm(1-(1-0.95)/2) * summary(model)$coefficients[, "Std.err"]
  upper <- coef(model) + qnorm(1-(1-0.95)/2) * summary(model)$coefficients[, "Std.err"]
  results <- data.frame(lower = lower, upper = upper)
  return(results)
}

#Simulation----
for (s in 1:simnum) {
# baseline covariates
agei <- rnorm(N, mean = 60, sd = 10)
femalei <- rbinom(N, size = 1, prob = 0.3)
bsa_bli <- rnorm(N, mean = 2, sd = 0.2)
ps_xbeta <- b1 + b2 * agei + b3 * femalei + b4 * bsa_bli
pscore <- exp(ps_xbeta)/ (1 + exp(ps_xbeta))
bavi <- rbinom(N, size = 1, prob = pscore)

# repeat baseline covariates maxT times
age <- rep(agei, each = maxT)
female <- rep(femalei, each = maxT)
bsa_bl <- rep(bsa_bli, each = maxT)
bav <- rep(bavi, each = maxT)

# generate outcome values
bi <- rnorm(N, mean = 0, sd = 5)
b <- rep(bi, each = maxT)
e <- rnorm(N*maxT, mean = 0, sd = 1)
root <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e

simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))
# * Matching----
# create sample data by matching patients based on ps
simdat_base <- simdat %>% group_by(id) %>% slice(1)
ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
#summary(ps)
ps.pm <- pairmatch(ps, data = simdat_base)
#summary(ps.pm) # 190 matched pairs
matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
  filter(!is.na(matches))
#bal.tab(ps.pm, covs = subset(simdat_base, select = c(age, female, bsa_bl)),distance = ps$fitted.values) 
matched_long <- simdat[simdat$id %in% matched_base$id,] # sample data


# fit gee with indep, exch, ar1, unstr to entire cohort and sample data and save estimates and ses 
# repeat 1000 times and get mean coefs and ses

# ** a) entire cohort----
  #1) independece
  gee_ind <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "independence")
  estbeta_ind[[s]] <- coef(gee_ind) # beta estimates
  se_ind[[s]] <- summary(gee_ind)$coefficient[,2] # standard errors
  ci_ind <- geeglm_ci(gee_ind) # 95% CI
  for (i in 1:length(beta_true)) {
    relbias_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ind[[s]][i] <- (estbeta_ind[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ind[i,1]&beta_true[i]<=ci_ind[i,2], cp_ind[[s]][i]<-1, cp_ind[[s]][i]<-0)
  }
  outvec_ind[[s]] <- c(s, N, estbeta_ind[[s]][1],estbeta_ind[[s]][2],estbeta_ind[[s]][3],
                       se_ind[[s]][1],se_ind[[s]][2],se_ind[[s]][3],relbias_ind[[s]][1],
                       relbias_ind[[s]][2],relbias_ind[[s]][3],mse_ind[[s]][1],mse_ind[[s]][2],
                       mse_ind[[s]][3],cp_ind[[s]][1],cp_ind[[s]][2],cp_ind[[s]][3])
  
  #2) exchangeable
  gee_exch <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "exchangeable")
  estbeta_exch[[s]] <- coef(gee_exch)
  se_exch[[s]] <- summary(gee_exch)$coefficient[,2]
  ci_exch <- geeglm_ci(gee_exch)
  for (i in 1:length(beta_true)) {
    relbias_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i]) / beta_true[i]
    mse_exch[[s]][i] <- (estbeta_exch[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_exch[i,1]&beta_true[i]<=ci_exch[i,2], cp_exch[[s]][i]<-1, cp_exch[[s]][i]<-0)
  }
  outvec_exch[[s]] <- c(s, N, estbeta_exch[[s]][1],estbeta_exch[[s]][2],estbeta_exch[[s]][3],
                        se_exch[[s]][1],se_exch[[s]][2],se_exch[[s]][3],relbias_exch[[s]][1],
                        relbias_exch[[s]][2],relbias_exch[[s]][3],mse_exch[[s]][1],mse_exch[[s]][2],
                        mse_exch[[s]][3],cp_exch[[s]][1],cp_exch[[s]][2],cp_exch[[s]][3])
  
  #3) AR1
  gee_ar1 <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                     id = id, waves = visit, corstr = "ar1")
  estbeta_ar1[[s]] <- coef(gee_ar1)
  se_ar1[[s]] <- summary(gee_ar1)$coefficient[,2] 
  ci_ar1 <- geeglm_ci(gee_ar1)
  for (i in 1:length(beta_true)) {
    relbias_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ar1[[s]][i] <- (estbeta_ar1[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ar1[i,1]&beta_true[i]<=ci_ar1[i,2], cp_ar1[[s]][i]<-1, cp_ar1[[s]][i]<-0)
  }
  outvec_ar1[[s]] <- c(s, N, estbeta_ar1[[s]][1],estbeta_ar1[[s]][2],estbeta_ar1[[s]][3],
                       se_ar1[[s]][1],se_ar1[[s]][2],se_ar1[[s]][3],relbias_ar1[[s]][1],
                       relbias_ar1[[s]][2],relbias_ar1[[s]][3],mse_ar1[[s]][1],mse_ar1[[s]][2],
                       mse_ar1[[s]][3],cp_ar1[[s]][1],cp_ar1[[s]][2],cp_ar1[[s]][3])
  
  #4) unstructured
  gee_unstr <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "unstructured")
  estbeta_unstr[[s]] <- coef(gee_unstr)
  se_unstr[[s]] <- summary(gee_unstr)$coefficient[,2] 
  ci_unstr <- geeglm_ci(gee_unstr)
  for (i in 1:length(beta_true)) {
    relbias_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i]) / beta_true[i]
    mse_unstr[[s]][i] <- (estbeta_unstr[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_unstr[i,1]&beta_true[i]<=ci_unstr[i,2], cp_unstr[[s]][i]<-1, cp_unstr[[s]][i]<-0)
  }
  outvec_unstr[[s]] <- c(s, N, estbeta_unstr[[s]][1],estbeta_unstr[[s]][2],estbeta_unstr[[s]][3],
                         se_unstr[[s]][1],se_unstr[[s]][2],se_unstr[[s]][3],relbias_unstr[[s]][1],
                         relbias_unstr[[s]][2],relbias_unstr[[s]][3],mse_unstr[[s]][1],mse_unstr[[s]][2],
                         mse_unstr[[s]][3],cp_unstr[[s]][1],cp_unstr[[s]][2],cp_unstr[[s]][3])
  
  # ** b) sample data----
  #1) independece
  gee_inds <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                     id = id, waves = visit, corstr = "independence")
  estbeta_inds[[s]] <- coef(gee_inds) # parameter estimates
  se_inds[[s]] <- summary(gee_inds)$coefficient[,2] # standard errors
  ci_inds <- geeglm_ci(gee_inds) # 95% CI
  for (i in 1:length(beta_true)) {
    relbias_inds[[s]][i] <- (estbeta_inds[[s]][i] - beta_true[i]) / beta_true[i]
    mse_inds[[s]][i] <- (estbeta_inds[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_inds[i,1]&beta_true[i]<=ci_inds[i,2], cp_inds[[s]][i]<-1, cp_inds[[s]][i]<-0)
  }
  outvec_inds[[s]] <- c(s, nrow(matched_base), estbeta_inds[[s]][1],estbeta_inds[[s]][2],estbeta_inds[[s]][3],
                       se_inds[[s]][1],se_inds[[s]][2],se_inds[[s]][3],relbias_inds[[s]][1],
                       relbias_inds[[s]][2],relbias_inds[[s]][3],mse_inds[[s]][1],mse_inds[[s]][2],
                       mse_inds[[s]][3],cp_inds[[s]][1],cp_inds[[s]][2],cp_inds[[s]][3])
  #2) exchangeable
  gee_exchs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                      id = id, waves = visit, corstr = "exchangeable")
  estbeta_exchs[[s]] <- coef(gee_exchs)
  se_exchs[[s]] <- summary(gee_exchs)$coefficient[,2]
  ci_exchs <- geeglm_ci(gee_exchs)
  for (i in 1:length(beta_true)) {
    relbias_exchs[[s]][i] <- (estbeta_exchs[[s]][i] - beta_true[i]) / beta_true[i]
    mse_exchs[[s]][i] <- (estbeta_exchs[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_exchs[i,1]&beta_true[i]<=ci_exchs[i,2], cp_exchs[[s]][i]<-1, cp_exchs[[s]][i]<-0)
  }
  outvec_exchs[[s]] <- c(s, nrow(matched_base), estbeta_exchs[[s]][1],estbeta_exchs[[s]][2],estbeta_exchs[[s]][3],
                        se_exchs[[s]][1],se_exchs[[s]][2],se_exchs[[s]][3],relbias_exchs[[s]][1],
                        relbias_exchs[[s]][2],relbias_exchs[[s]][3],mse_exchs[[s]][1],mse_exchs[[s]][2],
                        mse_exchs[[s]][3],cp_exchs[[s]][1],cp_exchs[[s]][2],cp_exchs[[s]][3])
  #3) AR1
  gee_ar1s <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                     id = id, waves = visit, corstr = "ar1")
  estbeta_ar1s[[s]] <- coef(gee_ar1s)
  se_ar1s[[s]] <- summary(gee_ar1s)$coefficient[,2]  
  ci_ar1s <- geeglm_ci(gee_ar1s)
  for (i in 1:length(beta_true)) {
    relbias_ar1s[[s]][i] <- (estbeta_ar1s[[s]][i] - beta_true[i]) / beta_true[i]
    mse_ar1s[[s]][i] <- (estbeta_ar1s[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_ar1s[i,1]&beta_true[i]<=ci_ar1s[i,2], cp_ar1s[[s]][i]<-1, cp_ar1s[[s]][i]<-0)
  }
  outvec_ar1s[[s]] <- c(s, nrow(matched_base), estbeta_ar1s[[s]][1],estbeta_ar1s[[s]][2],estbeta_ar1s[[s]][3],
                       se_ar1s[[s]][1],se_ar1s[[s]][2],se_ar1s[[s]][3],relbias_ar1s[[s]][1],
                       relbias_ar1s[[s]][2],relbias_ar1s[[s]][3],mse_ar1s[[s]][1],mse_ar1s[[s]][2],
                       mse_ar1s[[s]][3],cp_ar1s[[s]][1],cp_ar1s[[s]][2],cp_ar1s[[s]][3])
  
  #4) Unstructur
  gee_unstrs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                       id = id, waves = visit, corstr = "unstructured")
  estbeta_unstrs[[s]] <- coef(gee_unstrs)
  se_unstrs[[s]] <- summary(gee_unstrs)$coefficient[,2] 
  ci_unstrs <- geeglm_ci(gee_unstrs)
  for (i in 1:length(beta_true)) {
    relbias_unstrs[[s]][i] <- (estbeta_unstrs[[s]][i] - beta_true[i]) / beta_true[i]
    mse_unstrs[[s]][i] <- (estbeta_unstrs[[s]][i] - beta_true[i])^2
    ifelse(beta_true[i]>=ci_unstrs[i,1]&beta_true[i]<=ci_unstrs[i,2], cp_unstrs[[s]][i]<-1, cp_unstrs[[s]][i]<-0)
  }
  outvec_unstrs[[s]] <- c(s, nrow(matched_base), estbeta_unstrs[[s]][1],estbeta_unstrs[[s]][2],estbeta_unstrs[[s]][3],
                         se_unstrs[[s]][1],se_unstrs[[s]][2],se_unstrs[[s]][3],relbias_unstrs[[s]][1],
                         relbias_unstrs[[s]][2],relbias_unstrs[[s]][3],mse_unstrs[[s]][1],mse_unstrs[[s]][2],
                         mse_unstrs[[s]][3],cp_unstrs[[s]][1],cp_unstrs[[s]][2],cp_unstrs[[s]][3])
  print(s)
}

# output ----
# 1)  entire cohort outcome
names <- c("simnum", "sample_size", "estbeta0", "estbeta1", "estbeta2",
           "SE(beta0)", "SE(beta1)", "SE(beta2)","RelBias(b0)","RelBias(b1)",
           "RelBias(b2)", "MSE(b0)","MSE(b1)","MSE(b2)","CovProb(b0)",
           "CovProb(b1)","CovProb(b2)")

outvec_ind <- do.call("rbind", outvec_ind)
colnames(outvec_ind) <- names
outvec_exch <- do.call("rbind",outvec_exch)
colnames(outvec_exch) <- names
outvec_ar1 <- do.call("rbind", outvec_ar1)
colnames(outvec_ar1) <- names
outvec_unstr <- do.call("rbind", outvec_unstr)
colnames(outvec_unstr) <- names

outmean_ind <- c(simnum,colMeans(outvec_ind[,-1]), sd(outvec_ind[,"estbeta0"]),
                 sd(outvec_ind[,"estbeta1"]),sd(outvec_ind[,"estbeta2"]))
outmean_exch <- c(simnum,colMeans(outvec_exch[,-1]),sd(outvec_exch[,"estbeta0"]),
                  sd(outvec_exch[,"estbeta1"]),sd(outvec_exch[,"estbeta2"]))
outmean_ar1 <- c(simnum,colMeans(outvec_ar1[,-1]),sd(outvec_ar1[,"estbeta0"]),
                 sd(outvec_ar1[,"estbeta1"]),sd(outvec_ar1[,"estbeta2"]))
outmean_unstr <- c(simnum,colMeans(outvec_unstr[,-1]),sd(outvec_unstr[,"estbeta0"]),
                   sd(outvec_unstr[,"estbeta1"]),sd(outvec_unstr[,"estbeta2"]))
out_cohort <- rbind(outmean_ind,outmean_exch,outmean_ar1,outmean_unstr)
colnames(out_cohort)[c(1,(ncol(out_cohort)-2):ncol(out_cohort))] <- c("simnum","SD(b0)","SD(b1)","SD(b2)")
rownames(out_cohort) <- c("Independence", "Exchangeable","AR(1)","Unstructured")

table1 <- kable(out_cohort, digits = 3, escape = FALSE,
    caption = "GEE models comparison for the entire cohort",
    col.names = c("simnum","sample size", "$\\hat \\beta_0$", "$\\hat \\beta_1$", "$\\hat \\beta_2$",
                  "SE(b0)", "SE(b1)", "SE(b2)","RelBias(b0)","RelBias(b1)",
                  "RelBias(b2)", "MSE(b0)","MSE(b1)","MSE(b2)","CovProb(b0)",
                  "CovProb(b1)","CovProb(b2)","SD(b0)","SD(b1)","SD(b2)")) %>%
  kable_styling(full_width = F, position = "center")


# 2) sample data outcome
outvec_inds <- do.call("rbind", outvec_inds)
colnames(outvec_inds) <- names
outvec_exchs <- do.call("rbind",outvec_exchs)
colnames(outvec_exchs) <- names
outvec_ar1s <- do.call("rbind", outvec_ar1s)
colnames(outvec_ar1s) <- names
outvec_unstrs <- do.call("rbind", outvec_unstrs)
colnames(outvec_unstrs) <- names

outmean_inds <- c(simnum,colMeans(outvec_inds[,-1]), sd(outvec_inds[,"estbeta0"]),
                 sd(outvec_inds[,"estbeta1"]),sd(outvec_inds[,"estbeta2"]))
outmean_exchs <- c(simnum,colMeans(outvec_exchs[,-1]),sd(outvec_exchs[,"estbeta0"]),
                  sd(outvec_exchs[,"estbeta1"]),sd(outvec_exchs[,"estbeta2"]))
outmean_ar1s <- c(simnum,colMeans(outvec_ar1s[,-1]),sd(outvec_ar1s[,"estbeta0"]),
                 sd(outvec_ar1s[,"estbeta1"]),sd(outvec_ar1s[,"estbeta2"]))
outmean_unstrs <- c(simnum,colMeans(outvec_unstrs[,-1]),sd(outvec_unstrs[,"estbeta0"]),
                   sd(outvec_unstrs[,"estbeta1"]),sd(outvec_unstrs[,"estbeta2"]))
out_sample <- rbind(outmean_inds,outmean_exchs,outmean_ar1s,outmean_unstrs)
colnames(out_sample)[c(1,(ncol(out_sample)-2):ncol(out_sample))] <- c("simnum","SD(b0)","SD(b1)","SD(b2)")
rownames(out_sample) <- c("Independence", "Exchangeable","AR(1)","Unstructured")

table2 <- kable(out_sample, digits = 3, escape = FALSE,
              caption = "GEE models comparison for the sample data(matched pairs)",
              col.names = c("simnum","sample size", "$\\hat \\beta_0$", "$\\hat \\beta_1$", "$\\hat \\beta_2$",
                            "SE(b0)", "SE(b1)", "SE(b2)","RelBias(b0)","RelBias(b1)",
                            "RelBias(b2)", "MSE(b0)","MSE(b1)","MSE(b2)","CovProb(b0)",
                            "CovProb(b1)","CovProb(b2)","SD(b0)","SD(b1)","SD(b2)")) %>%
  kable_styling(full_width = F, position = "center")


comb <- rbind(c(rep("",ncol(out_cohort))), round(out_cohort,3),
      c(rep("",ncol(out_sample))) ,round(out_sample,3))
rownames(comb)[c(1,6)] <- c("Entire cohort","Sample data(matched pairs)")

knitr::kable(comb, digits = 3, escape = FALSE,
      caption = "GEE models comparison for difference correlation structures",
      col.names = c("simulation number","sample size", "$\\hat \\beta_0$", "$\\hat \\beta_1$", "$\\hat \\beta_2$",
                    "SE(b0)", "SE(b1)", "SE(b2)","Relative Bias(b0)","Relative Bias(b1)",
                    "Relative Bias(b2)", "MSE(b0)","MSE(b1)","MSE(b2)","Coverage Prob(b0)",
                    "Coverage Prob(b1)","Coverage Prob(b2)","SD(b0)","SD(b1)","SD(b2)")) %>%
  kable_styling(full_width = F, position = "center")

test <- simdat %>% group_by(id) %>% slice(1)
table(test$bav)
mean(test$bav)

