# ------------------------------------ #
# Script de test des différents priors #
# Auteur : G. Mulier                   #
# Créé le 02/04/2024, modifié 02/04/24 #
# ------------------------------------ #

# Packages et helpers ----

library(rstan)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)



# Définition des modèles MCMC ----

## Power prior ----

## Un seul essai ----

Modele1 <- "
data{
  int<lower = 1> N;
  int<lower = 0> x;
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = 0> gamm;
}
parameters{
  real<lower = 0, upper = 1> theta;
}
transformed parameters{
  real logVrais = x * log(theta) + (N - x) * log1m(theta);
}
model{
  target += beta_lpdf(theta | a, b); // Prior
  target += gamm * logVrais; // Vraissemblance avec puissance
}
"
Compil1 <- stan_model(model_code = Modele1)

## Essai multiples ----
ModeleM <- "
data{
  int<lower = 1> Nt;
  int<lower = 1> N[Nt];
  int<lower = 0> x[Nt];
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = 0> gamm;
}
parameters{
  real<lower = 0, upper = 1> theta[Nt];
}
transformed parameters{
}
model{
  target += beta_lpdf(theta | a, b); // Prior
  for (i in 1:Nt) {
    target += gamm * (x[i] * log(theta[i]) + (N[i] - x[i]) * log1m(theta[i])); // Vraissemblance avec puissance
  } 
}
"
CompilM <- stan_model(model_code = ModeleM)

## Hiérarchique ----
ModeleH <- "
data{
  int<lower = 1> Nt;
  int<lower = 1> N[Nt];
  int<lower = 0> x[Nt];
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = 0> gamm;
}
parameters{
  real<lower = 0, upper = 1> theta[Nt];
  real<lower = 0, upper = 1> mu;
  real<lower = 0> sigma;
}
transformed parameters{
  real alpha = - mu * (sigma * sigma + mu * mu - mu) / (sigma * sigma);
  real beta = (sigma * sigma + mu * mu - mu) * (mu - 1) / (sigma * sigma);
}
model{
  // Hyperpriors
  mu  ~ normal(0, 1000);
  sigma ~ normal(0, 1000);
  for (i in 1:Nt) {
    target += beta_lpdf(theta[i] | alpha, beta); // Prior
    target += gamm * (x[i] * log(theta[i]) + (N[i] - x[i]) * log1m(theta[i])); // Vraissemblance avec puissance
  } 
}
"
CompilH <- stan_model(model_code = ModeleH)

# MAP Prior ----
ModeleMAP <- "
data {
  int <lower = 2> Nt; // Nombre d'essais
  int<lower = 1> n[Nt];
  int<lower = 0> y[Nt];
}

parameters{
  real p_raw[Nt];
  real mu;
  real<lower = 0> sigma;
}

transformed parameters{
  real p[Nt];
  real logit_p[Nt];
  for(j in 1:Nt){
    logit_p[j] = mu + sigma * p_raw[j]; // Non central parametrization
  }
  p = inv_logit(logit_p);
}

model{

  // prior distributions
  sigma ~ normal(0, 10);
  mu  ~ normal(0, 10);

  for (i in 1:Nt) {
    // binomial likelihood
    p_raw[i] ~ normal(0, 1);
    y[i] ~ binomial(n[i], p[i]);
  }
 
}

generated quantities {
  real pnew;
  pnew = normal_rng(mu, sigma);
}
"
CompilMAP <- stan_model(model_code = ModeleMAP)

mixBeta <- function(x, ab, w=NULL){
  # ab list of vectors with elements a, b (shape parameters of each beta component)
  n.comp <- length(ab)
  if (n.comp == 1) {out <- dbeta(x, shape1 = ab[[1]][[1]], shape2 = ab[[1]][[2]])}
  else {
    if (length(w) == (n.comp - 1)) { w <- c(w, 1-sum(w)) }
    
    out <- 0
    for (i in 1:n.comp){ out <- out + w[i]*dbeta(x, shape1 = ab[[i]][[1]], shape2 = ab[[i]][[2]])}
  }
  return(out)
}


MAP.Binary2PriorML2C <- function(fit, ...){
  
  mix.2c <- fitdistr(x = round(extract(fit, pars = "pnew")$pnew, 4),
                     densfun = function(x, a1, b1, a2, b2, w1){
                       mixBeta(x, ab = list(c(a1, b1), c(a2, b2)), w = w1)
                     },
                     start = list(
                       a1 = 5,
                       b1 = 45,
                       a2 = 1,
                       b2 = 5,
                       w1 = 0.8
                     ),
                     lower = c(0.001, 0.001, 0.001, 0.001, 0.001),
                     upper = c(10000, 10000, 10000, 10000, 0.999))
  
  ab <- list(c(mix.2c$estimate[['a1']], mix.2c$estimate[['b1']]),
             c(mix.2c$estimate[['a2']], mix.2c$estimate[['b2']]))
  w <- mix.2c$estimate[['w1']]
  
  prior <- function(x,...){
    y <- mixBeta(x, ab = ab, w = w)
    return(y)
  }
  
  attributes(prior) <- list(params = list(ab = ab, w = w))
  
  return(prior)
}

MAP.Binary2PriorML <- function(fit, ...){
  mle.fit <- fitdistr(x = 1 / (1 + exp(-rstan::extract(fit, pars = "pnew")$pnew)),
                      densfun = 'beta',
                      ...)
  
  prior <- function(x,...){
    y <- dbeta(x, shape1 = mle.fit$estimate[['shape1']], shape2 = mle.fit$estimate[['shape2']])
    return(y)
  }
  
  attributes(prior) <- list(params = mle.fit$estimate)
  
  return(prior)
}

# Générer des données et appliquer les différents priors ----

NEssais <- 3
NSimul <- 10000
NPtsBras <- c(40, 60, 20)
ParamsBeta <- list(c(10, 10), c(10, 30), c(60, 20))
NPtsEssais <- matrix(rep(NPtsBras, each = NSimul), ncol = NEssais, nrow = NSimul)
Resultats <- list()

cl <- makeCluster(50)
registerDoParallel(cl)
set.seed(121221, kind = "L'Ecuyer-CMRG")


for (loi in seq_along(ParamsBeta)) {
  # En quelque sorte les hyperparamètres, et la simulations des paramètres p de la loi binomiale
  cat(paste0("\n\nScénario n°", loi, "\n\n"), file = "~/R/log.txt", append = TRUE)
  AlphaB <- ParamsBeta[[loi]][1]
  BetaB <- ParamsBeta[[loi]][2]
  ProbasEssais <- matrix(rbeta(NSimul * NEssais, AlphaB, BetaB), ncol = NEssais, nrow = NSimul)
  ResultatsEssais <- matrix(rbinom(NSimul * NEssais, NPtsEssais, ProbasEssais), ncol = NEssais, nrow = NSimul)
  Resultats[[loi]] <- foreach(essai = seq_len(NSimul),
                              .combine = "rbind",
                              .packages = c("MASS", "rstan"),
                              .export = c("Compil1", "CompilM", "CompilH", "CompilMAP", "MAP.Binary2PriorML", "AlphaB", "BetaB")) %dopar% {
                                if (essai %% 50 == 0) cat(paste0("simu ", essai, "\n"), file = "~/R/log.txt", append = TRUE)
                                # Power prior données poolées
                                DonneesRegroup <- list(N = sum(NPtsBras), x = sum(ResultatsEssais[essai, ]), a = .5, b = .5, gamm = .75)
                                FitPPrior <- sampling(Compil1, 
                                                      data = DonneesRegroup,         
                                                      chains = 4,             
                                                      warmup = 5000,          
                                                      iter = 25000,
                                                      thin = 10,
                                                      cores = 1)
                                mcmcpp1 <- as.data.frame(cbind(simu = essai, var = "theta", algo = "MCMC", prior = "PP1",
                                                               t(summary(FitPPrior, pars = "theta")$summary[, c(1, 4, 8)])))
                                names(mcmcpp1)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                APP1 <- DonneesRegroup$a + DonneesRegroup$gamm * DonneesRegroup$x
                                BPP1 <- DonneesRegroup$b + DonneesRegroup$gamm * (DonneesRegroup$N - DonneesRegroup$x)
                                theopp1 <- data.frame("simu" = essai, "var" = "theta", "algo" = "analytique", "prior" = "PP1", 
                                                      "moy" = APP1 / (APP1 + BPP1), "icinf" = qbeta(.025, APP1, BPP1), 
                                                      "icsup" = qbeta(.975, APP1, BPP1))
                                # Power prior multiples
                                Donnees <- list(Nt = NEssais, N = NPtsBras, x = ResultatsEssais[essai, ], a = .5, b = .5, gamm = .75)
                                FitPPrior2 <- sampling(CompilM, 
                                                       data = Donnees,         
                                                       chains = 4,             
                                                       warmup = 5000,          
                                                       iter = 25000,
                                                       thin = 10,
                                                       cores = 1)
                                mcmcpp2 <- as.data.frame(cbind("simu" = essai, 
                                                               "var" = paste0("theta", seq_len(NEssais)),
                                                               "algo" = "MCMC", "prior" = "PP2",
                                                               data.frame(summary(FitPPrior2, pars = "theta")$summary[, c(1, 4, 8)])))
                                names(mcmcpp2)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                APP2 <- Donnees$a + Donnees$gamm * Donnees$x
                                BPP2 <- Donnees$b + Donnees$gamm * (Donnees$N - Donnees$x)
                                theopp2 <- data.frame("simu" = essai, "var" = paste0("theta", seq(1, NEssais)), "algo" = "analytique", "prior" = "PP2", 
                                                      "moy" = APP2 / (APP2 + BPP2), "icinf" = qbeta(.025, APP2, BPP2), 
                                                      "icsup" = qbeta(.975, APP2, BPP2))
                                # Power prior hiérarchique
                                FitPPriorH <- sampling(CompilH, 
                                                       data = Donnees,         
                                                       chains = 4,             
                                                       warmup = 5000,          
                                                       iter = 25000,
                                                       thin = 10,
                                                       cores = 1)
                                mcmcpp3 <- as.data.frame(cbind("simu" = essai, 
                                                               "var" = c(paste0("theta", seq_len(NEssais)), "mu", "sigma"),
                                                               "algo" = "MCMC", "prior" = "PP3",
                                                               data.frame(summary(FitPPriorH, pars = c("theta", "mu", "sigma"))$summary[, c(1, 4, 8)])))
                                names(mcmcpp3)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                theopp3 <- data.frame("simu" = essai, "var" = c(paste0("theta", seq(1, NEssais)), "mu"), 
                                                      "algo" = "analytique", "prior" = "PP3", 
                                                      "moy" = c(APP2, APP1) / (c(APP2, APP1) + c(BPP2, BPP1)), "icinf" = qbeta(.025, c(APP2, APP1), c(BPP2, BPP1)), 
                                                      "icsup" = qbeta(.975, c(APP2, APP1), c(BPP2, BPP1)))
                                # MAP prior
                                DonneesMAP <- list(Nt = NEssais, n = NPtsBras, y = ResultatsEssais[essai, ])
                                FitPPriorMAP <- sampling(CompilMAP, 
                                                         data = DonneesMAP,         
                                                         chains = 4,             
                                                         warmup = 50000,          
                                                         iter = 100000,
                                                         thin = 20,
                                                         cores = 1)
                                Parametres <- tryCatch(attr(MAP.Binary2PriorML(FitPPriorMAP, start = list(shape1 = AlphaB, shape2 = BetaB), lower = c(0.001, 0.001)), "params"),
                                                       error = function(e) {
                                                         return(c(NA, NA))
                                                       })
                                MAPPrior <- data.frame("simu" = essai, "var" = "theta", 
                                                       "algo" = "MCMC", "prior" = "MAP", 
                                                       "moy" = Parametres[[1]] / (Parametres[[1]] + Parametres[[2]]), "icinf" = qbeta(.025, Parametres[[1]], Parametres[[2]]), 
                                                       "icsup" = qbeta(.975, Parametres[[1]], Parametres[[2]]))
                                return(rbind(mcmcpp1, theopp1, mcmcpp2, theopp2, mcmcpp3, theopp3, MAPPrior))
                                
                              }
}

NEssais <- 10
NSimul <- 10000
NPtsBras <- rep(40, NEssais)
ParamsBeta <- list(c(10, 10), c(10, 30), c(60, 20))
NPtsEssais <- matrix(rep(NPtsBras, each = NSimul), ncol = NEssais, nrow = NSimul)
Resultats2 <- list()

for (loi in seq_along(ParamsBeta)) {
  # En quelque sorte les hyperparamètres, et la simulations des paramètres p de la loi binomiale
  cat(paste0("\n\nScénario n°", loi, "\n\n"), file = "~/R/log.txt", append = TRUE)
  AlphaB <- ParamsBeta[[loi]][1]
  BetaB <- ParamsBeta[[loi]][2]
  ProbasEssais <- matrix(rbeta(NSimul * NEssais, AlphaB, BetaB), ncol = NEssais, nrow = NSimul)
  ResultatsEssais <- matrix(rbinom(NSimul * NEssais, NPtsEssais, ProbasEssais), ncol = NEssais, nrow = NSimul)
  Resultats2[[loi]] <- foreach(essai = seq_len(NSimul),
                               .combine = "rbind",
                               .packages = c("MASS", "rstan"),
                               .export = c("Compil1", "CompilM", "CompilH", "CompilMAP", "MAP.Binary2PriorML", "AlphaB", "BetaB")) %dopar% {
                                 if (essai %% 50 == 0) cat(paste0("simu ", essai, "\n"), file = "~/R/log.txt", append = TRUE)
                                 # Power prior données poolées
                                 DonneesRegroup <- list(N = sum(NPtsBras), x = sum(ResultatsEssais[essai, ]), a = .5, b = .5, gamm = .75)
                                 FitPPrior <- sampling(Compil1, 
                                                       data = DonneesRegroup,         
                                                       chains = 4,             
                                                       warmup = 5000,          
                                                       iter = 25000,
                                                       thin = 10,
                                                       cores = 1)
                                 mcmcpp1 <- as.data.frame(cbind(simu = essai, var = "theta", algo = "MCMC", prior = "PP1",
                                                                t(summary(FitPPrior, pars = "theta")$summary[, c(1, 4, 8)])))
                                 names(mcmcpp1)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                 APP1 <- DonneesRegroup$a + DonneesRegroup$gamm * DonneesRegroup$x
                                 BPP1 <- DonneesRegroup$b + DonneesRegroup$gamm * (DonneesRegroup$N - DonneesRegroup$x)
                                 theopp1 <- data.frame("simu" = essai, "var" = "theta", "algo" = "analytique", "prior" = "PP1", 
                                                       "moy" = APP1 / (APP1 + BPP1), "icinf" = qbeta(.025, APP1, BPP1), 
                                                       "icsup" = qbeta(.975, APP1, BPP1))
                                 # Power prior multiples
                                 Donnees <- list(Nt = NEssais, N = NPtsBras, x = ResultatsEssais[essai, ], a = .5, b = .5, gamm = .75)
                                 FitPPrior2 <- sampling(CompilM, 
                                                        data = Donnees,         
                                                        chains = 4,             
                                                        warmup = 5000,          
                                                        iter = 25000,
                                                        thin = 10,
                                                        cores = 1)
                                 mcmcpp2 <- as.data.frame(cbind("simu" = essai, 
                                                                "var" = paste0("theta", seq_len(NEssais)),
                                                                "algo" = "MCMC", "prior" = "PP2",
                                                                data.frame(summary(FitPPrior2, pars = "theta")$summary[, c(1, 4, 8)])))
                                 names(mcmcpp2)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                 APP2 <- Donnees$a + Donnees$gamm * Donnees$x
                                 BPP2 <- Donnees$b + Donnees$gamm * (Donnees$N - Donnees$x)
                                 theopp2 <- data.frame("simu" = essai, "var" = paste0("theta", seq(1, NEssais)), "algo" = "analytique", "prior" = "PP2", 
                                                       "moy" = APP2 / (APP2 + BPP2), "icinf" = qbeta(.025, APP2, BPP2), 
                                                       "icsup" = qbeta(.975, APP2, BPP2))
                                 # Power prior hiérarchique
                                 FitPPriorH <- sampling(CompilH, 
                                                        data = Donnees,         
                                                        chains = 4,             
                                                        warmup = 5000,          
                                                        iter = 25000,
                                                        thin = 10,
                                                        cores = 1)
                                 mcmcpp3 <- as.data.frame(cbind("simu" = essai, 
                                                                "var" = c(paste0("theta", seq_len(NEssais)), "mu", "sigma"),
                                                                "algo" = "MCMC", "prior" = "PP3",
                                                                data.frame(summary(FitPPriorH, pars = c("theta", "mu", "sigma"))$summary[, c(1, 4, 8)])))
                                 names(mcmcpp3)[seq(5, 7)] <- c("moy", "icinf", "icsup")
                                 theopp3 <- data.frame("simu" = essai, "var" = c(paste0("theta", seq(1, NEssais)), "mu"), 
                                                       "algo" = "analytique", "prior" = "PP3", 
                                                       "moy" = c(APP2, APP1) / (c(APP2, APP1) + c(BPP2, BPP1)), "icinf" = qbeta(.025, c(APP2, APP1), c(BPP2, BPP1)), 
                                                       "icsup" = qbeta(.975, c(APP2, APP1), c(BPP2, BPP1)))
                                 # MAP prior
                                 DonneesMAP <- list(Nt = NEssais, n = NPtsBras, y = ResultatsEssais[essai, ])
                                 FitPPriorMAP <- sampling(CompilMAP, 
                                                          data = DonneesMAP,         
                                                          chains = 4,             
                                                          warmup = 50000,          
                                                          iter = 100000,
                                                          thin = 20,
                                                          cores = 1)
                                 Parametres <- tryCatch(attr(MAP.Binary2PriorML(FitPPriorMAP, start = list(shape1 = AlphaB, shape2 = BetaB), lower = c(0.001, 0.001)), "params"),
                                                        error = function(e) {
                                                          return(c(NA, NA))
                                                        })
                                 MAPPrior <- data.frame("simu" = essai, "var" = "theta", 
                                                        "algo" = "MCMC", "prior" = "MAP", 
                                                        "moy" = Parametres[[1]] / (Parametres[[1]] + Parametres[[2]]), "icinf" = qbeta(.025, Parametres[[1]], Parametres[[2]]), 
                                                        "icsup" = qbeta(.975, Parametres[[1]], Parametres[[2]]))
                                 return(rbind(mcmcpp1, theopp1, mcmcpp2, theopp2, mcmcpp3, theopp3, MAPPrior))
                                 
                               }
}

stopCluster(cl)

save(Resultats, Resultats2, file = "~/R/res_tests_priors.RData")


