#############################################################################
# Example: Weightloss RE study level
# A numerical solution to the threshold problem, using INLA to provide fast
# posterior approximations.
#############################################################################

# Install INLA package (not available on CRAN)
if (!require(INLA)) {
  install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")
}


# Data setup --------------------------------------------------------------

# Read in the raw data
dat.raw <- read.delim("../WinBUGS/WeightlossRE_dataset.txt")
dat.raw$b <- dat.raw$t.1

# Reshape the data to have one row per arm per study
dat.raw2 <- reshape(dat.raw, varying=c("t.1","y.1","se.1","t.2","y.2","se.2",
                                       "t.3","y.3","se.3","t.4","y.4","se.4"), 
                    timevar="arm", idvar="studyID", direction="long")

# Sort data by study and contrast, removing NA rows
dat.raw2 <- dat.raw2[order(dat.raw2$studyID, dat.raw2$arm, dat.raw2$y, na.last=NA),]

# Likelihood covariance matrix (note arms are independent)
lik.cov <- diag(dat.raw2$se^2)

K <- length(unique(dat.raw2$t)) # number of treatments
n <- nrow(dat.raw)  # number of studies
N <- nrow(dat.raw2)  # number of data points

# Flag for baseline rows
dat.raw2$bl <- rep(FALSE,nrow(dat.raw2))
dat.raw2$bl[1] <- TRUE
dat.raw2$bl[cumsum(table(dat.raw2$studyID))[-n]+1] <- TRUE

# Design matrix for mu parameters
M <- matrix(0, ncol=n, nrow=N)
for (i in 1:N) M[i,dat.raw2$studyID[i]] <- 1

# # Design matrix for delta parameters
# L <- matrix(0, ncol=sum(dat.raw$n_arms)-nrow(dat.raw), nrow=N)
# j <- 1
# for (i in 1:N) {
#   if (dat.raw2$arm[i] != 1) {
#     L[i,j] <- 1
#     j<-j+1
#   }
# }
# 
# # Design matrix for d parameters
# X <- matrix(0, ncol=K-1, nrow=sum(L))
# temp <- dat.raw2[rowSums(L) > 0,c("b","t")]
# X[cbind(which(temp$t>1),temp[temp$t>1,"t"]-1)] <- 1
# X[cbind(which(temp$b>1),temp[temp$b>1,"b"]-1)] <- -1
# 
# # Write matrix L%*%X with NA rows for the baseline arms
# LX <- L %*% X
# LX[dat.raw2$bl,] <- NA

# Design matrix for d parameters, equivalent to LX
#   - Baseline rows are all zero
#   - All other rows just look like rows of X
#   - We are careful that rows for arms comparing
#     baseline vs. baseline are also all zero

LX <- matrix(0, ncol=K-1, nrow=N)
attach(dat.raw2)
for(i in 1:N){
  if(t[i] != b[i] && t[i] > 1) LX[i,t[i] - 1] <- 1
  if(t[i] != b[i] && b[i] > 1) LX[i,b[i] - 1] <- -1
}
detach(dat.raw2)

# "Design" matrix A for Sigma_n = A*sigma^2
require(Matrix)
A.diag <- vector("list",n)
attach(dat.raw)
for (i in 1:n){
  if(n_arms[i]==2){
    A.diag[[i]] <- 1
  }
  else if(n_arms[i]>2){
    A.diag[[i]] <- matrix(.5,nrow=n_arms[i]-1,ncol=n_arms[i]-1)
    diag(A.diag[[i]]) <- 1
  }
}
detach(dat.raw)
A <- bdiag(A.diag)

# # Create random effects index - modified from Sauter and Held (2015)
# dat.raw2$re <- dat.raw2$studyID
# dat.raw2[dat.raw2$bl,"re"] <- NA   # baseline arms are NA
# 
# # Create random effects grouping - modified from Sauter and Held (2015)
# dat.raw2$grp <- NA
# dat.raw2[!dat.raw2$bl,"grp"] <- unlist(sapply(dat.raw$n_arms-1,seq))

# Create random effects index
# Note that all non-baseline arms have a random effect, even if the arm
# comparator is the same as the baseline - then the RE is zero mean, and
# the study just adds to the estimation of sigma^2.
dat.raw2$re2 <- NA
dat.raw2$re2[!dat.raw2$bl] <- 1:sum(!dat.raw2$bl)

# # Create random effects treatment grouping
# dat.raw2$trtgrp <- NA
# dat.raw2$trtgrp[!dat.raw2$bl] <- dat.raw2$t[!dat.raw2$bl]


# INLA Model --------------------------------------------------------------

# Define uniform prior for between studies sd - code from Sauter and Held (2015)

# Upper limit for uniform distribution:
ul <- 5
# Function for Uniform distribution:
hyperunif.function <- function(x) {
  if (exp(x) ^ -0.5 < ul & exp(x) ^ -0.5 > 0) {
    logdens <- log(1 / ul)
  }else{
    logdens <- log(0.1e-320)
  }
  logdenst <- logdens + log(0.5 * exp(-x / 2))
  return(logdenst)
}

# Define grid:
lprec <- seq(from = -40, to = 40, len = 20000)
# Create table:
prior.table <-
  paste(c("table:", cbind(
    lprec,  sapply(lprec, FUN = hyperunif.function)
  )),
  sep = "", collapse = " ")


# # Compute sig2 correlations for each study on INLA internal scale
# # Sauter and Held (2015) cite Riebler et al. (2012) for the transform
# rho <- 0.5 # Intra-arm correlation is 0.5
# max_narms <- max(dat.raw$n_arms) # Maximum number of arms
# cor_sig2 <- log((1 + rho * (max_narms - 1 )) / (1 - rho))


# For pairwise MA:
# - Replace M with a studyID variable(?). Or diag(Nstudy). Not needed - contrast level data.
# - Replace LX with a column of 1s, or equivalently add in the global intercept which will
#   be interpreted as the treatment contrast
# - A will be a diagonal matrix of 1s, dimension NstudyxNstudy

# prior=prior.table specifies the uniform prior for the heterogeneity standard deviation
#  For simplicity, try (with a trunc. normal prior on tau)
#    f(studyID, model = "iid",
#      hyper = list(prec = list(prior = "logtnormal", param = c(mean0, prec0))))

#### SUGGESTED MODEL FOR DAVID ####
# Model for a pairwise meta-analysis with 2-arm trials

# RE_mod <- y ~ -1 +   # no global intercept
#               d +    # treatment effect parameter, vector of ones in data. Or leave out, and use global intercept instead.
#                      # random effects, with truncated normal prior on standard deviation
#               f(studyID, model = "iid",
#                 hyper = list(prec = list(
#                   prior = "logtnormal", param = c(mean0, prec0))
#                   )
#                 )

###################################

# Model Formula
RE_mod <- y ~ -1 +    # no global intercept
              M +     # study-level baseline mu
              LX +    # treatment contrasts delta
                      # random effects
              f(re2, model="generic0", Cmatrix=solve(A),
                hyper=list(prec=list(prior=prior.table)))

####################################################################
## This RE code from Sauter and Held (2015) yields
## identical results, but is harder to understand!!
#               f(re, model = "iid",
#                   hyper = list(prec = list(prior = prior.table)),
#                   group = grp,
#                   control.group = list(model = "exchangeable",
#                                        hyper = list(rho = list(
#                                          fixed = TRUE,
#                                          initial = cor_sig2
#                                          ))
#                                        )
#                 )
####################################################################

# Run INLA
RE_INLA <- inla(RE_mod, data = dat.raw2, 
            # Likelihood distribution
                family = "gaussian",
            # Fix likelihood hyperpars as the data is fixed with known precision.
                control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
            # Likelihood precisions
                scale = 1/diag(lik.cov),
            # Prior distribution for "fixed" effects - really for mu and d
                control.fixed = list(mean = 0, prec = 0.0001, 
                                   expand.factor.strategy = "inla",
            # Optionally compute correlation matrix of mu and d
                                   correlation.matrix = TRUE),
            # Optionally compute DIC
                control.compute = list(dic=TRUE))

# OPTIONAL: Increase accuracy of hyperparameter (tau) estimation. Doesn't do much for
# us here.
#RE_INLA <- inla.hyperpar(RE_INLA)

summary(RE_INLA)


# Compare to WinBUGS ------------------------------------------------------

# Get WinBUGS CODA
require(coda)
# Parallelise to reduce computation time
require(parallelsugar) # Fix mclapply on Windows, from GitHub nathanvan/parallelsugar

system.time(
  dat.CODA <- mcmc.list(mclapply(c("../WinBUGS/WeightlossRE/Coda1.txt",
                                   "../WinBUGS/WeightlossRE/Coda2.txt",
                                   "../WinBUGS/WeightlossRE/Coda3.txt"), 
                                 read.coda, 
                                 index.file = "../WinBUGS/WeightlossRE/CodaIndex.txt"))
)

print(codastat <- summary(dat.CODA))

# Get indices of variables
vnames <- sub("(.*)\\[.*","\\1",varnames(dat.CODA))
ind.d <- which(vnames=="d")
ind.delta <- which(vnames=="delta")
ind.mu <- which(vnames=="mu")
ind.theta <- which(vnames=="theta")
ind.sd <- which(vnames=="sd")

vnames.inla <- RE_INLA$names.fixed
ind.d.inla <- which(stringr::str_detect(vnames.inla,"^LX"))
ind.mu.inla <- which(stringr::str_detect(vnames.inla,"^M"))

## Plot and compare densities
gp <- par(mfrow=c(2,2), mar=c(3,4,3,1))

marginals_dk <- RE_INLA$marginals.fixed[ind.d.inla]

for (i in 1:(K-1)){
  # Histogram
  temphist <- hist(as.matrix(dat.CODA[,ind.d[i]]), breaks=50, plot=F)
  plot(temphist, ylim=c(0,max(temphist$density)*1.1), freq=F, 
       main=paste0("d[",i+1,"]"), xlab="", col="lightblue", border="lightgrey")

  # INLA
  lines(inla.smarginal(marginals_dk[[i]]), col="red", lwd=1)
}

par(gp)

## Study baseline parameters mu
gp <- par(mfrow=c(5,5), mar=c(0,0.5,2,0.5))

marginals_mu <- RE_INLA$marginals.fixed[ind.mu.inla]

for (i in 1:ncol(M)){
  # Histogram
  temphist <- hist(as.matrix(dat.CODA[,ind.mu[i]]), breaks=50, plot=F)
  plot(temphist, ylim=c(0,max(temphist$density)*1.2), freq=F, 
       main=paste0("mu[",i,"]"), xlab="", ylab="", 
       xaxt="n", yaxt="n",
       col="lightblue", border="lightgrey")
  
  # INLA
  lines(inla.smarginal(marginals_mu[[i]]), col="red", lwd=1)
}

par(gp)

gp <- par(mfrow=c(1,2))
## Between studies variance
# Histogram
temphist <- hist(as.matrix(dat.CODA[,ind.sd])^2, breaks=50, plot=F)
plot(temphist, xlim=c(0,15), ylim=c(0,max(temphist$density)*1.1), freq=F, 
     main=expression(sigma^2), xlab="", col="lightblue", border="lightgrey")

# INLA
lines(inla.tmarginal(function(x) 1/x, RE_INLA$marginals.hyperpar[[1]]), 
      col="red", lwd=1, type="l")

## Between studies precision
# Histogram
temphist <- hist(as.matrix(dat.CODA[,ind.sd])^-2, breaks=50, plot=F)
plot(temphist, xlim=c(0,1), ylim=c(0,max(temphist$density)*1.1), freq=F, 
     main=expression(1/sigma^2), xlab="", col="lightblue", border="lightgrey")

# INLA
lines(inla.smarginal(RE_INLA$marginals.hyperpar[[1]]), 
      col="red", lwd=1, type="l")
par(gp)

## Check differences in correlation matrices
cov.diff <- RE_INLA$misc$lincomb.derived.covariance.matrix - cov(as.matrix(dat.CODA[,c(ind.mu,ind.d)]))
# log10(abs(cov.diff))
source("../../../../utils/matrixplot.R")
matrixplot(cov.diff)
matrixplot(ceiling(log10(abs(cov.diff))))

## Check differences in means
mean.diff <- c(RE_INLA$summary.fixed$mean, 
               inla.emarginal(function(x) 1/x, RE_INLA$marginals.hyperpar[[1]])) - 
              c(codastat$statistics[c(ind.mu, ind.d),"Mean"], mean(as.matrix(dat.CODA[,ind.sd])^2))
