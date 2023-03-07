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
agei <- rnorm(N, mean = 60, sd = 10)
femalei <- rbinom(N, size = 1, prob = 0.3)
bsa_bli <- rnorm(N, mean = 2, sd = 0.2)
ps_xbeta <- 4 + -0.1 * agei + -1.3 * femalei + 1.5 * bsa_bli
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
root <- 30 + 0.5 * bav + -0.2 * visit + 0.5 * bav * visit + b + e

simdat <- as.data.frame(cbind(id, visit, age, female, bsa_bl, bav, root))

# create sample data by matching patients based on ps


# fit gee with indep, exch, ar1, unstr to entire cohort and sample data and save estimates and ses 

# repeat 1000 times and get mean coefs and ses
