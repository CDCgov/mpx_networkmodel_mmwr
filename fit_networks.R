

## Libraries ----------------------------------------------------------------
library(EpiModel)
library(EpiModelHIV)

## Network characteristics ---------------------------------------------------

# Original model was built to have ability to have different behavior among White and Black MSM in Atlanta
# While there are parameter names that for W and B characteristics, they are the same
# (race is not actually in this model)
# Differences in behavior among come from sexual activity groups not defined by race 
# Not all targets outlined in this section are used final networks

pop.size <- 10000
num.B <- num.W <- pop.size/2

# mean/pers degree distributions matrices.
deg.mp.B <- deg.mp.W <-
  (matrix(c(0.471, 0.167, 0.074,
            0.22, 0.047, 0.021), byrow = TRUE, nrow = 2))

# Instant rates
mdeg.inst.B <- mdeg.inst.W <-
  (matrix(c(0.065/7, 0.087/7, 0.086/7,
            0.056/7, 0.055/7, 0.055/7), byrow = TRUE, nrow = 2))

top5 <- TRUE 
# Instant rates (with additional partners based on highest risk group)
# We add to each category weighted on original dist of instant rates by other partnerships 
if (top5==TRUE) {
  mdeg.inst.B <- mdeg.inst.W <-
    (matrix(c(0.022542267, 0.030171957,0.029825153,
              0.01942103, 0.019074226, 0.019074226), byrow = TRUE, nrow = 2))
}

# Quintile distribution of overall AI rates
# adding risk group for top 5% based on Weiss et al. 2020
qnts.W <- qnts.B <- c(0.0000,
                      0.007/7,
                      0.038/7,
                      0.071/7,
                      0.221/7)

if (top5==TRUE) {
  qnts.W <- qnts.B <- c(0.0000,
                        0.007/7,
                        0.038/7,
                        0.071/7,
                        0.221/7,
                        2/7)
}

# Proportion in same-race partnerships (main, casl, inst)
prop.hom.mpi.B <- prop.hom.mpi.W <- (c(0.9484, 0.9019, 0.9085) +
                                       c(0.9154, 0.8509, 0.8944))/2

# Mean age diffs (main, casl, inst)
sqrt.adiff.BB <- c(0.464, 0.586, 0.544)
sqrt.adiff.BW <- c(0.464, 0.586, 0.544)
sqrt.adiff.WW <- c(0.464, 0.586, 0.544)

# Mean durations
rates.main <- mean(c(1/407,
                     1/407,
                     1/407))
rates.pers <- mean(c(1/166,
                     1/166,
                     1/166))

durs.main <- 1/rates.main 
durs.pers <- 1/rates.pers 

# Age-sex-specific mortality rates
ages <- 18:39
asmr.B <- c(rep(0, 17),
            1-(1-c(rep(0.00159, 7),
                   rep(0.00225, 10),
                   rep(0.00348, 5)))^(1/(365/1)), 1)

asmr.W <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 5)))^(1/(365/1)), 1)

# I, R, V role frequencies
role.B.prob <- role.W.prob <- (c(0.242, 0.321, 0.437))

time.unit = 1 # days

mdeg.hiv.pers = c(0.38631, 0.54083)
mdeg.hiv.inst = c(0.022021, 0.034353)


## Create Target Statistics -------------------------------------------
source("setup_functions.R")

st <- calc_nwstats_msm(
  method = 1, #1 for 1 race models, 2 for two race models
  top5 = top5,
  time.unit = time.unit,
  num.B = num.B,
  num.W = num.W,
  deg.mp.B = deg.mp.B,
  deg.mp.W = deg.mp.W,
  mdeg.inst.B = mdeg.inst.B,
  mdeg.inst.W = mdeg.inst.W,
  qnts.B = qnts.B,
  qnts.W = qnts.W,
  prop.hom.mpi.B = prop.hom.mpi.B,
  prop.hom.mpi.W = prop.hom.mpi.W,
  balance = "mean",
  sqrt.adiff.BB = sqrt.adiff.BB,
  sqrt.adiff.WW = sqrt.adiff.WW,
  sqrt.adiff.BW = sqrt.adiff.BW,
  diss.main = ~offset(edges),
  diss.pers = ~offset(edges),
  durs.main = durs.main,
  durs.pers = durs.pers,
  ages = ages,
  asmr.B = asmr.B,
  asmr.W = asmr.W,
  role.B.prob = role.B.prob,
  role.W.prob = role.W.prob,
  mdeg.hiv.pers = mdeg.hiv.pers,
  mdeg.hiv.inst = mdeg.hiv.inst)


# 1. Main Model -----------------------------------------------------------

# Initialize network
nw.main <- base_nw_msm(st)

# Assign degree
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)

# Formulas
formation.m <- ~edges +
  nodefactor("deg.pers") +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.m,
                coef.diss = st$coef.diss.m,
                constraints = ~bd(maxout = 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))


# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)

# Formulas
formation.p <- ~edges +
  nodefactor("deg.main") +
  concurrent +
  absdiff("sqrt.age") +
  nodefactor("hiv", levels=2) +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

# Fit model
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.p,
                coef.diss = st$coef.diss.p,
                constraints = ~bd(maxout = 2),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))



# Fit inst model ----------------------------------------------------------

# Initialize network
nw.inst <- nw.main

# Assign degree
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main")
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers")
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers")

# Formulas
formation.i <- ~edges +
  nodefactor(c("deg.main", "deg.pers")) +
  nodefactor("riskg", levels = -2) +
  absdiff("sqrt.age") +
  nodefactor("hiv", levels=2) +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = st$stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e7,
                                                MCMLE.maxit = 250))

est <- list(fit.m, fit.p, fit.i) 


save(est, file = "mpx_network_object.rda")

#m.dx <- netdx(fit.m, dynamic=TRUE, nsims=5, nsteps=500)
#p.dx <- netdx(fit.p, dynamic=TRUE, nsims=5, nsteps=500)
#i.dx <- netdx(fit.i, dynamic=FALSE, nsims=500)
