# -------------------------------------------------------------------- #
# Exploration of the influence of the random seed for CBHM FWER        #
# Auteur : G. Mulier                                                   #
# Créé le 05/11/2025, modifié le 05/11/2025                            #
# -------------------------------------------------------------------- #

# Packages  and helpers ----

# github.com/GuillaumeMulier/multibrasBOP2
# devtools::load_all("~/pkgs/multibrasBOP2/")
cat("Charging packages and models...\n\n", file = "~/simu_priors/log.txt")
library(dplyr)
library(rlang)
library(stringr)
library(purrr)
library(parallel)
library(foreach)
library(doParallel)
library(rstan)
source("~/R/GetOC.R")
# source("R/GetOC.R")

# Stan models ----

# Precompile stan models 
cat("----\nCompiling efficacy/toxicity STAN\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
CompiledModelsEffTox <- list(
  "cbhmBOP_efftox" = CompilCBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, 
                                       moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5)
)


# 3 arms simulation ----

## Parameters for the simulation ----

Alpha <- .1
NSimu <- 5000
AnaEff <- rep(29, 2)
AnaTox <- rep(29, 2)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)
NBras <- 3

# Scénarios
Scenarios <- list(
  "Sc1"  = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45))
)
NomScenars <- names(Scenarios)

# Parameters for CBHM
CBHM_eff <- CalibrateCBHM(NPts = rep(sum(AnaEff), NBras), Q0 = sum(PN[c(1, 2)]), Q1 = sum(PA[c(1, 2)]))
CBHM_tox <- CalibrateCBHM(NPts = rep(sum(AnaTox), NBras), Q0 = sum(PN[c(1, 3)]), Q1 = sum(PA[c(1, 3)]))

# Seuil
if (FALSE) {
  SeuilBOP <- deter_cutoff(alpha = Alpha, 
                           n_bras = NBras,
                           nsim_oc = 10000,
                           ana_inter = AnaEff, ana_inter_tox = AnaTox, 
                           p_n = PN, p_a = PA,
                           power_seq = seq(0, 1.5, .01), 
                           seed = 121221, affich_mat = "No", methode = 4L)
} else {
  SeuilBOP <- list(c("C_" = .76, "gamma" = 1.31, "alpha_calc" = .0913, "puissance_calc" = .8064))
}
# The explored methods
Methodes <- list(
  "cbhmBOP_efftox" = list(methode = "cbhm_efftox", A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)


##  Simulations ----

cl <- makeCluster(10)
registerDoParallel(cl)

if (TRUE) {
  cat("----\nMain simulation: 3 arms\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  # Parameters of the method
  Params <- Methodes[[1]]
  
  # Simulate the scenarios 
  set.seed(121221)
  VecSeeds <- c(121221, sample.int(1e+9, 9))
  # VecScenarios <- seq_along(Scenarios)
  ResT <- foreach(i = VecSeeds,
                  .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                  .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "CBHM_eff",
                              "CBHM_tox", "gen_patients_multinom",
                              "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                              "opcharac", "PA", "Params", "PN", "real_essai_bayeslog_efftox",
                              "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                              "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                              "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                              "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                              "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                              "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                              "real_essai_modhier_efftox", "real_essai_modhier_tox", "simu_simon",
                              "Scenarios", "SeuilBOP", "summarise_decision", "summarise_detect",
                              "summarise_ttt")) %dopar% {
                                
                                NScenar <- 1
                                NumMethode <- 1
                                cat(paste0("Seed n°", i, "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                
                                # Generate trials
                                tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                        multinom_ttt = Scenarios[[NScenar]],
                                                                        rand_ratio = rep(1, NBras), seed = i)
                                
                                # Simulate the result with the different methods
                                Resultat <- tryCatch(opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                              p_n = PN, p_a = PA,
                                                              CPar = SeuilBOP[[1]][["C_"]], PPar = SeuilBOP[[1]][["gamma"]],
                                                              methode = Params$methode,
                                                              A0_tox = Params$A0, SeuilP_tox = Params$SeuilP,
                                                              A0_eff = Params$A0, SeuilP_eff = Params$SeuilP,
                                                              a_tox = CBHM_tox$a, b_tox = CBHM_tox$b,
                                                              a_eff = CBHM_eff$a, b_eff = CBHM_eff$b,
                                                              p_mix_tox = rep(.5, NBras), p_mix_eff = rep(.5, NBras),
                                                              tableau_essais = tableau_essais),
                                                     error = function(e) list(paste0(e)))
                                
                                Resultat <- append(Resultat,
                                                   values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                   after = 2)
                                return(Resultat)
                                
                              }
  
  # Save results
  save(ResT, file = "~/simu_priors/resultats_cbhmseed_20251105.RData")
  
} else {
  cat("Decided not to do it (main analysis). Continuing with next simulation scenarios...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
}

stopCluster(cl)

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

