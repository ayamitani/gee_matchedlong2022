### generate big cohort 

# total sample size
N <- 500
# total number of visits
maxT <- 5
# vector of visits for all patients
visit <- rep(1:maxT, 500)
# vector of patient ids
id <- rep(1:N, each = maxT)

# baseline covariates
b1 <- 4
b2 <- -0.1
b3 <- -1.3
b4 <- 1.5
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
b1_r <- 30
b2_r <- 0.5
b3_r <- -0.2
b4_r <- 0.5
bi <- rnorm(N, mean = 0, sd = 5)
b <- rep(bi, each = maxT)
e <- rnorm(N*maxT, mean = 0, sd = 1)
root <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e

simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))


# create sample data by matching patients based on ps
simdat_base <- simdat %>% group_by(id) %>% slice(1)
ps <- glm(bav ~ age + female + bsa_bl, family = binomial, data = simdat_base)
summary(ps)

ps.pm <- pairmatch(ps, data = simdat_base)
summary(ps.pm) # 190 matched pairs
matched_base <- data.frame(simdat_base, matches = ps.pm, check.rows = TRUE) %>%
  filter(!is.na(matches))

# check balance
covs <- subset(simdat_base, select = c(age, female, bsa_bl))
p.score <- ps$fitted.values
bal.tab(ps.pm, covs = covs, distance = ~ p.score)

matched_long <- simdat[simdat$id %in% matched_base$id,] # sample data


# fit gee with indep, exch, ar1, unstr to entire cohort and sample data and save estimates and ses 
# repeat 1000 times and get mean coefs and ses

# 1) entire cohort
simnum <- 1000
estbeta_ind <- vector("list", length = simnum)
estbeta_exch <- vector("list", length = simnum)
estbeta_ar1 <- vector("list", length = simnum)
estbeta_unstr <- vector("list", length = simnum)
se_ind <- vector("list", length = simnum)
se_exch <- vector("list", length = simnum)
se_ar1 <- vector("list", length = simnum)
se_unstr <- vector("list", length = simnum)

for (s in 1:simnum) {
  gee_ind <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "independence")
  estbeta_ind[[s]] <- coef(gee_ind) # parameter estimates
  se_ind[[s]] <- summary(gee_ind)$coefficient[,2] # standard errors
  
  gee_exch <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "exchangeable")
  estbeta_exch[[s]] <- coef(gee_exch)
  se_exch[[s]] <- summary(gee_exch)$coefficient[,2]

  gee_ar1 <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                     id = id, waves = visit, corstr = "ar1")
  estbeta_ar1[[s]] <- coef(gee_ar1)
  se_ar1[[s]] <- summary(gee_ar1)$coefficient[,2]  
  
  gee_unstr <- geeglm(root ~ visit + bav:visit, data = simdat, family = gaussian,
                    id = id, waves = visit, corstr = "unstructured")
  estbeta_unstr[[s]] <- coef(gee_unstr)
  se_unstr[[s]] <- summary(gee_unstr)$coefficient[,2] 
}

beta_mean_ind <- colMeans(do.call("rbind", estbeta_ind)) 
se_mean_ind <- colMeans(do.call("rbind",se_ind))
beta_mean_exch <- colMeans(do.call("rbind", estbeta_exch)) 
se_mean_exch <- colMeans(do.call("rbind",se_exch))
beta_mean_ar1 <- colMeans(do.call("rbind", estbeta_ar1)) 
se_mean_ar1 <- colMeans(do.call("rbind",se_ar1))
beta_mean_unstr <- colMeans(do.call("rbind", estbeta_unstr)) 
se_mean_unstr <- colMeans(do.call("rbind",se_unstr))

out_cohort <- cbind(beta_mean_ind,se_mean_ind,beta_mean_exch, se_mean_exch, 
                    beta_mean_ar1,se_mean_ar1, beta_mean_unstr,se_mean_unstr)
kbl(out_cohort,digit = 2, escape = FALSE, caption = "GEE estimates and SEs for the entire cohort",
    col.names = rep(c("$\\hat \\beta$","SE"), 4)) %>%
  kable_styling(full_width = F, position = "center") %>%
  add_header_above(c(" ", "Ind" = 2, "Exch" = 2, "AR1" = 2, "Unstr"=2))


# 2) sample data
estbeta_inds <- vector("list", length = simnum)
estbeta_exchs <- vector("list", length = simnum)
estbeta_ar1s <- vector("list", length = simnum)
estbeta_unstrs <- vector("list", length = simnum)
se_inds <- vector("list", length = simnum)
se_exchs <- vector("list", length = simnum)
se_ar1s <- vector("list", length = simnum)
se_unstrs <- vector("list", length = simnum)

for (s in 1:simnum) {
  gee_inds <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                    id = id, waves = visit, corstr = "independence")
  estbeta_inds[[s]] <- coef(gee_inds) # parameter estimates
  se_inds[[s]] <- summary(gee_inds)$coefficient[,2] # standard errors
  
  gee_exchs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                     id = id, waves = visit, corstr = "exchangeable")
  estbeta_exchs[[s]] <- coef(gee_exchs)
  se_exchs[[s]] <- summary(gee_exchs)$coefficient[,2]
  
  gee_ar1s <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                    id = id, waves = visit, corstr = "ar1")
  estbeta_ar1s[[s]] <- coef(gee_ar1s)
  se_ar1s[[s]] <- summary(gee_ar1s)$coefficient[,2]  
  
  gee_unstrs <- geeglm(root ~  visit + bav:visit, data = matched_long, family = gaussian,
                      id = id, waves = visit, corstr = "unstructured")
  estbeta_unstrs[[s]] <- coef(gee_unstrs)
  se_unstrs[[s]] <- summary(gee_unstrs)$coefficient[,2] 
}

beta_mean_inds <- colMeans(do.call("rbind", estbeta_inds)) 
se_mean_inds <- colMeans(do.call("rbind",se_inds))
beta_mean_exchs <- colMeans(do.call("rbind", estbeta_exchs)) 
se_mean_exchs <- colMeans(do.call("rbind",se_exchs))
beta_mean_ar1s <- colMeans(do.call("rbind", estbeta_ar1s)) 
se_mean_ar1s <- colMeans(do.call("rbind",se_ar1s))
beta_mean_unstrs <- colMeans(do.call("rbind", estbeta_unstrs)) 
se_mean_unstrs <- colMeans(do.call("rbind",se_unstrs))

out_sample <- cbind(beta_mean_inds,se_mean_inds,beta_mean_exchs, se_mean_exchs, 
                    beta_mean_ar1s,se_mean_ar1s, beta_mean_unstrs,se_mean_unstrs)
kbl(out_sample,digit = 2, escape = FALSE, caption = "GEE estimates and SEs for sample data",
    col.names = rep(c("$\\hat \\beta$","SE"), 4)) %>%
  kable_styling(full_width = F, position = "center") %>%
  add_header_above(c(" ", "Ind" = 2, "Exch" = 2, "AR1" = 2, "Unstr"=2))

