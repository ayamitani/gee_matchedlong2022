# Reproduce Simulation

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