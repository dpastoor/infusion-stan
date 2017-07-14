# dev: edit r script which generates the data
rm(list = ls())
gc()

modelName <- "oneCptPop"
# modelName <- "TwoCptModelPopulation"

## directory paths, assuming we are in the lab-notebook
## directory.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
modelDir <- file.path(projectDir, "stanModel")
toolsDir <- file.path("tools")
stanDir <- file.path(projectDir, "cmdstan")
tempDir <- file.path(modelDir, modelName, "temp")

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "cmdStanTools.R"))

## load helpful libraries
# .libPaths("~/svn-StanPmetrics/script/lib")
library(rstan)
library(parallel)

set.seed(11191989)

parametersToPlot <- c("CLHat", "VHat", "sigma_add", "sigma_prop", "lp__")

## Randomly generate initial estimates
nIIV <- 2
nSubjects <- 50
init <- function()
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       VHat = exp(rnorm(1, log(70), 0.2)),
       sigma = runif(1, 0.5, 2),
       sigma_prop = runif(1, 0.5, 2),
       L = diag(nIIV),
       etaStd = matrix(rep(0, nIIV * nSubjects), nrow = nIIV),
       omega = runif(nIIV, 0.5, 2))

# nIIV <- 5
# nSubjects <- 50
# init <- function()
#   list(CLHat = exp(rnorm(1, log(10), 0.2)),
#        QHat = exp(rnorm(1, log(20), 0.2)),
#        V1Hat = exp(rnorm(1, log(70), 0.2)),
#        V2Hat = exp(rnorm(1, log(70), 0.2)),
#        kaHat = exp(rnorm(1, log(1), 0.2)),
#        sigma = runif(1, 0.5, 2),
#        L = diag(nIIV),
#        etaStd = matrix(rep(0, nIIV * nSubjects), nrow = nIIV),
#        omega = runif(nIIV, 0.5, 2),
#        logtheta = matrix(rep(log(c(exp(rnorm(1, log(10), 0.2)),
#                                    exp(rnorm(1, log(20), 0.2)),
#                                    exp(rnorm(1, log(70), 0.2)),
#                                    exp(rnorm(1, log(70), 0.2)),
#                                    exp(rnorm(1, log(1), 0.2)))),
#                              ea = nSubjects),
#                          nrow = nSubjects))

###############################################################################
## run stan
nChains <- 1  # 4
chains <- 1:nChains
nPost <- 1  # 1000
nBurn <- 0  # 1000
nThin <- 1

nIter <- nPost * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init) {
           tempDir <- file.path(tempDir, chain)
           dir.create(tempDir)
           inits <- init()
           with(inits, stan_rdump(ls(inits), file = file.path(tempDir,
                                                              "init.R")))
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = file.path(tempDir, "init.R"), 
                    seed = sample(1:999999, 1),
                    chain = chain, refresh = 100,
                    adapt_delta = 0.95, stepsize = 0.01)
         },
         model = file.path(modelDir, modelName),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = init,
         iter = nIter, warmup = nBurnin, thin = nThin,
         mc.cores = min(nChains, detectCores()))
