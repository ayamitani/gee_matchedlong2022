### generate matched sample

# matched sample 
N <- 25
# total sample
M <- 2*N
# total number of visits
maxT <- 5
# vector of visits for all matches
visit <- rep(1:maxT, N)
# vector of matched ids
matchid <- rep(1:N, each = maxT)

# true baseline coefficients
b1_r <- 30
b2_r <- 0.5
b3_r <- -0.2
b4_r <- 0.5

bi <- rnorm(N, mean = 0, sd = 5)
b <- rep(bi, each = maxT)

# generate outcome for exposed group (bav = 1)
bav <- 1
e <- rnorm(N*maxT, mean = 0, sd = 3)
root1 <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e

# generate outcome for exposed group (bav = 0)
bav <- 0
e <- rnorm(N*maxT, mean = 0, sd = 3)
root0 <- b1_r + b2_r * bav + b3_r * visit + b4_r * bav * visit + b + e

simdat <- as.data.frame(cbind(matchid, visit, root0, root1)) %>%
  pivot_longer(cols = starts_with("root"),
               names_to = "bav",
               names_prefix = "root",
               values_to = "root") %>%
  mutate(patientid = ) # ADD PATIENT ID

# FIT THE FOUR GEE MODELS 
# MEAN BETA
# SD of BETA
# MEAN SE(BETA)
# MEAN BIAS = mean(betahat - truebeta)
# MSE = mean((betahat - truebeta)^2)
# COVERAGE PROBABILITY 
# compute 95% confidence interval for each betahat
# see if truebeta lies in the 95% CI
# how many times does truebeta lies in the 95% CI out of the 1000 simulations

# for each simulation
# outputvec <- c(s, sample_sample, estbeta0ind, estbeta1ind, estbeta2ind, sebeta0ind, sebeta1ind, sebeta2ind, biasbeta0ind, sqerrbeta0ind, covbeta0ind,
# give names names(outputvec) <- c("simnum", "sample_size", "estbeta0id",...)
# repeat for ar1, exch, unstr)
# simresult <- outside of forloop create a list of length 1000
# simresult[[s]] <- outputvec
# output simresult
# outside of the forloop, do.call rbind
