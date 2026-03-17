# -------------------------------------------------------------------- #
# Simulations to answer reviewers' questions                           #
# Author : G. Mulier                                                   #
# Created the 16/03/2026, modified the 13/03/2026                      #
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
  # "hBOP_efftox" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = 1,
  #                                  moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = 1),
  # "cbhmBOP_efftox" = CompilCBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, 
  #                                      moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5),
  # "log1_efftox" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 2.5,
  #                                     mu_inter_tox = log(.4 / .6), sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 2.5),
  # "log2_efftox" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 2.5,
  #                                     mu_inter_tox = log(.4 / .6), sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 2.5,
  #                                     PentePos_eff = TRUE, PentePos_tox = TRUE),
  "normpowBOP_efftox" = CompilNPower_efftox(phi_eff = .3, phi_tox = .4)
)
cat("----\nCompiling toxicity STAN\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
CompiledModelsTox <- list(
  # "hBOP_tox" = CompilBHM_tox(moy_mu = log(.4 / .6), sigma_mu = 2.5, moy_sig = 0, sigma_sig = 1),
  # "cbhmBOP_tox" = CompilCBHM_tox(moy_mu = log(.4 / .6), sigma_mu = 2.5),
  # "log1_tox" = CompilModLog_tox(mu_inter = log(.4 / .6), sigma_inter = 2.5, mu_coef = 0.22, sigma_coef = 2.5),
  # "log2_tox" = CompilModLog_tox(mu_inter = log(.4 / .6), sigma_inter = 2.5, mu_coef = 0.22, sigma_coef = 2.5, PentePos = TRUE),
  "normpowBOP_tox" = CompilNPower_tox(phi_tox = .4, a = 1, b = 1)
)

# Normalized power prior ----

## 3 arms simulation ----

### Parameters for the simulation ----

Alpha <- .1
NSimu <- 5000
AnaEff <- rep(29, 2)
AnaTox <- rep(29, 2)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)
NBras <- 3

# Scénarios
Scenarios <- list(
  "Sc1"  = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
  "Sc2"  = list(ttt1 = c(0.13, 0.12, 0.27, 0.48), ttt2 = c(0.15, 0.13, 0.27, 0.45), ttt3 = c(0.16, 0.14, 0.29, 0.41)),
  "Sc3"  = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
  "Sc4"  = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.17, 0.35, 0.11, 0.37), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI1" = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.19, 0.36, 0.11, 0.34), ttt3 = c(0.18, 0.34, 0.12, 0.36)),
  "ScI3" = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.18, 0.34, 0.12, 0.36), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  "Sc5"  = list(ttt1 = c(0.11, 0.19, 0.17, 0.53), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  "Sc6"  = list(ttt1 = c(0.11, 0.19, 0.17, 0.53), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.21, 0.21, 0.39)),
  "Sc7"  = list(ttt1 = c(0.14, 0.26, 0.14, 0.46), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  "Sc8"  = list(ttt1 = c(0.18, 0.32, 0.12, 0.38), ttt2 = c(0.22, 0.28, 0.15, 0.35), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
  "Sc9"  = list(ttt1 = c(0.12, 0.18, 0.18, 0.52), ttt2 = c(0.17, 0.23, 0.18, 0.42), ttt3 = c(0.23, 0.27, 0.17, 0.33))
)
NomScenars <- names(Scenarios)

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
  "normpowBOP_efftox" = list(methode = "bop_normpower_efftox", A0 = NA, SeuilP = NA),
  "normpowBOP_tox" = list(methode = "bop_normpower_tox", A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)


###  Simulations ----

cl <- makeCluster(20)
registerDoParallel(cl)

if (TRUE) {
  cat("----\nMain simulation: 3 arms (normalized power prior)\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  for (m in seq_len(ceiling(length(Methodes) / 2))) {
    
    # Load 
    if (2 * m <= length(Methodes)) {
      cat(paste0("----\n", NomMethodes[2 * m - 1], " and ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    } else {
      cat(paste0("----\n", NomMethodes[2 * m - 1], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    }
    
    # Parameters of the method
    Params1 <- Methodes[[2 * m - 1]]
    if (2 * m <= length(Methodes)) Params2 <- Methodes[[2 * m]]
    
    # Simulate the scenarios 
    VecScenarios <- if (2 * m <= length(Methodes)) seq_len(2 * length(Scenarios)) else seq_along(Scenarios)
    ResT <- foreach(i = VecScenarios,
                    .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                    .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "gen_patients_multinom",
                                "m", "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                                "opcharac", "PA", "Params1", "PN", "real_essai_bayeslog_efftox",
                                "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                                "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                                "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                                "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                                "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                                "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                                "real_essai_modhier_efftox", "real_essai_modhier_tox", "simu_simon",
                                "Scenarios", "SeuilBOP", "summarise_decision", "summarise_detect",
                                "summarise_ttt")) %dopar% {
                                  
                                  NScenar <- i %% length(Scenarios)
                                  if (NScenar == 0) NScenar <- length(Scenarios)
                                  NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                  Params <- if (i <= length(Scenarios)) Params1 else Params2
                                  cat(paste0("Scenario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                  
                                  # Generate trials
                                  tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                          multinom_ttt = Scenarios[[NScenar]],
                                                                          rand_ratio = rep(1, NBras), seed = 121221)
                                  
                                  # Simulate the result with the different methods
                                  if (Params$methode == "simoniva") {
                                    Resultat <- simu_simon(p_n = PN, p_a = PA,
                                                           tableau_essais = tableau_essais,
                                                           CaracSeuilSimon = SeuilSimon,
                                                           CaracSeuilIva = SeuilIva)
                                  } else {
                                    Resultat <- tryCatch(opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                                  p_n = PN, p_a = PA,
                                                                  CPar = SeuilBOP[[1]][["C_"]], PPar = SeuilBOP[[1]][["gamma"]],
                                                                  methode = Params$methode,
                                                                  A0_tox = Params$A0, SeuilP_tox = Params$SeuilP,
                                                                  A0_eff = Params$A0, SeuilP_eff = Params$SeuilP,
                                                                  a_tox = NA, b_tox = NA,
                                                                  a_eff = NA, b_eff = NA,
                                                                  p_mix_tox = rep(.5, NBras), p_mix_eff = rep(.5, NBras),
                                                                  tableau_essais = tableau_essais),
                                                         error = function(e) list(paste0(e)))
                                  }
                                  
                                  Resultat <- append(Resultat,
                                                     values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                     after = 2)
                                  return(Resultat)
                                  
                                }
    
    # Save results
    save(ResT, file = paste0("~/simu_priors/resultats_priorsppalnorm_20260316_", m, ".RData"))
    
  }
  
} else {
  cat("Decided not to do it (main analysis for normalized power prior). Continuing with next simulation scenarios...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
}

stopCluster(cl)


## 4 arms (sensitivity analysis) ----

cl <- makeCluster(8)
registerDoParallel(cl)

### Parameters for the simulation ----

Alpha <- .1
NSimu <- 5000
AnaEff <- rep(35, 2)
AnaTox <- rep(35, 2)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)
NBras <- 4

# Scénarios
Scenarios <- list(
  "Sc2"  = list(ttt1 = c(0.12, 0.12, 0.28, 0.48), ttt2 = c(0.13, 0.13, 0.29, 0.45), ttt3 = c(0.15, 0.13, 0.29, 0.43), ttt4 = c(0.17, 0.13, 0.33, 0.37)),
  "Sc4"  = list(ttt1 = c(0.15, 0.35, 0.09, 0.41), ttt2 = c(0.16, 0.36, 0.10, 0.38), ttt3 = c(0.18, 0.36, 0.10, 0.36), ttt4 = c(0.20, 0.36, 0.10, 0.34)),
  "ScI1" = list(ttt1 = c(0.07, 0.18, 0.13, 0.62), ttt2 = c(0.10, 0.20, 0.15, 0.55), ttt3 = c(0.18, 0.32, 0.12, 0.38), ttt4 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.12, 0.38, 0.08, 0.42), ttt2 = c(0.16, 0.36, 0.09, 0.39), ttt3 = c(0.19, 0.35, 0.11, 0.35), ttt4 = c(0.25, 0.31, 0.15, 0.29))
)
NomScenars <- names(Scenarios)

# Seuil
if (FALSE) {
  SeuilBOP4 <- deter_cutoff(alpha = Alpha, 
                            n_bras = NBras,
                            nsim_oc = 10000,
                            ana_inter = AnaEff, ana_inter_tox = AnaTox, 
                            p_n = PN, p_a = PA,
                            power_seq = seq(0, 1.5, .01), 
                            seed = 121221, affich_mat = "No", methode = 4L)
} else {
  SeuilBOP4 <- list(c("C_" = .775, "gamma" = 1.05, "alpha_calc" = .099, "puissance_calc" = .8075))
}


### Simulations ----

if (TRUE) {
  cat("----\nSensitivity analysis (normalized power prior) : 4 arms\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  for (m in seq_len(ceiling(length(Methodes) / 2))) {
    
    # Load 
    if (2 * m <= length(Methodes)) {
      cat(paste0("----\n", NomMethodes[2 * m - 1], " and ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    } else {
      cat(paste0("----\n", NomMethodes[2 * m - 1], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    }
    
    # Parameters of the method
    Params1 <- Methodes[[2 * m - 1]]
    if (2 * m <= length(Methodes)) Params2 <- Methodes[[2 * m]]
    
    # Simulate the scenarios 
    VecScenarios <- if (2 * m <= length(Methodes)) seq_len(2 * length(Scenarios)) else seq_along(Scenarios)
    ResT <- foreach(i = VecScenarios,
                    .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                    .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "gen_patients_multinom",
                                "m", "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                                "opcharac", "PA", "Params1", "PN", "real_essai_bayeslog_efftox",
                                "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                                "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                                "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                                "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                                "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                                "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                                "real_essai_modhier_efftox", "real_essai_modhier_tox", "simu_simon",
                                "Scenarios", "SeuilBOP4", "summarise_decision", "summarise_detect",
                                "summarise_ttt")) %dopar% {
                                  
                                  NScenar <- i %% length(Scenarios)
                                  if (NScenar == 0) NScenar <- length(Scenarios)
                                  NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                  Params <- if (i <= length(Scenarios)) Params1 else Params2
                                  cat(paste0("Scenario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                  
                                  # Générer essais
                                  tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                          multinom_ttt = Scenarios[[NScenar]],
                                                                          rand_ratio = rep(1, NBras), seed = 121221)
                                  
                                  # Simuler le résultat
                                  if (Params$methode == "simoniva") {
                                    Resultat <- simu_simon(p_n = PN, p_a = PA,
                                                           tableau_essais = tableau_essais,
                                                           CaracSeuilSimon = SeuilSimon4,
                                                           CaracSeuilIva = SeuilIva4)
                                  } else {
                                    Resultat <- tryCatch(opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                                  p_n = PN, p_a = PA,
                                                                  CPar = SeuilBOP4[[1]][["C_"]], PPar = SeuilBOP4[[1]][["gamma"]],
                                                                  methode = Params$methode,
                                                                  A0_tox = Params$A0, SeuilP_tox = Params$SeuilP,
                                                                  A0_eff = Params$A0, SeuilP_eff = Params$SeuilP,
                                                                  a_tox = CBHM_tox4$a, b_tox = CBHM_tox4$b,
                                                                  a_eff = CBHM_eff4$a, b_eff = CBHM_eff4$b,
                                                                  p_mix_tox = rep(.5, NBras), p_mix_eff = rep(.5, NBras),
                                                                  tableau_essais = tableau_essais),
                                                         error = function(e) list(paste0(e)))
                                  }
                                  Resultat <- append(Resultat,
                                                     values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                     after = 2)
                                  return(Resultat)
                                  
                                }
    
    # Save results
    save(ResT, file = paste0("~/simu_priors/resultats_priorssens4norm_20260316_", m, ".RData"))
    
  }
  
} else {
  cat("Decided not to do it (4 arms sensitivity for normalized power prior). Continuing with next simulation scenarios...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
}


## 5 arms (sensitivity analysis, normalized power prior) ----

### Parameters for the simulation ----

Alpha <- .1
NSimu <- 5000
AnaEff <- rep(35, 2)
AnaTox <- rep(35, 2)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)
NBras <- 5

# Scénarios
Scenarios <- list(
  "Sc2"  = list(ttt1 = c(0.11, 0.11, 0.29, 0.49), ttt2 = c(0.13, 0.11, 0.29, 0.47), ttt3 = c(0.14, 0.12, 0.30, 0.44), ttt4 = c(0.15, 0.13, 0.31, 0.41), ttt5 = c(0.17, 0.13, 0.31, 0.39)),
  "Sc4"  = list(ttt1 = c(0.13, 0.37, 0.09, 0.41), ttt2 = c(0.15, 0.36, 0.09, 0.40), ttt3 = c(0.16, 0.36, 0.10, 0.38), ttt4 = c(0.18, 0.35, 0.10, 0.37), ttt5 = c(0.19, 0.35, 0.11, 0.35)),
  "ScI1" = list(ttt1 = c(0.08, 0.17, 0.16, 0.59), ttt2 = c(0.10, 0.20, 0.16, 0.54), ttt3 = c(0.17, 0.33, 0.11, 0.39), ttt4 = c(0.19, 0.34, 0.11, 0.36), ttt5 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.12, 0.38, 0.08, 0.42), ttt2 = c(0.14, 0.37, 0.08, 0.41), ttt3 = c(0.16, 0.36, 0.09, 0.39), ttt4 = c(0.19, 0.34, 0.11, 0.36), ttt5 = c(0.25, 0.29, 0.15, 0.31))
)
NomScenars <- names(Scenarios)

# Seuil
if (FALSE) {
  SeuilBOP5 <- deter_cutoff(alpha = Alpha, 
                            n_bras = NBras,
                            nsim_oc = 10000,
                            ana_inter = AnaEff, ana_inter_tox = AnaTox, 
                            p_n = PN, p_a = PA,
                            power_seq = seq(0, 1.5, .01), 
                            seed = 121221, affich_mat = "No", methode = 4L)
} else {
  SeuilBOP5 <- list(c("C_" = .79, "gamma" = 1.05, "alpha_calc" = .0984, "puissance_calc" = .7997))
}

### Simulations ----

if (TRUE) {
  cat("----\nSensitivity analysis : 5 arms (normalized power prior)\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  for (m in seq_len(ceiling(length(Methodes) / 2))) {
    
    # Load 
    if (2 * m <= length(Methodes)) {
      cat(paste0("----\n", NomMethodes[2 * m - 1], " and ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    } else {
      cat(paste0("----\n", NomMethodes[2 * m - 1], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    }
    
    # Parameters of the method
    Params1 <- Methodes[[2 * m - 1]]
    if (2 * m <= length(Methodes)) Params2 <- Methodes[[2 * m]]
    
    # Simulate the scenarios 
    VecScenarios <- if (2 * m <= length(Methodes)) seq_len(2 * length(Scenarios)) else seq_along(Scenarios)
    # VecScenarios <- seq_along(Scenarios)
    ResT <- foreach(i = VecScenarios,
                    .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                    .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "gen_patients_multinom",
                                "m", "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                                "opcharac", "PA", "Params1", "PN", "real_essai_bayeslog_efftox",
                                "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                                "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                                "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                                "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                                "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                                "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                                "real_essai_modhier_efftox", "real_essai_modhier_tox", "simu_simon",
                                "Scenarios", "SeuilBOP5", "summarise_decision", "summarise_detect",
                                "summarise_ttt")) %dopar% {
                                  
                                  NScenar <- i %% length(Scenarios)
                                  if (NScenar == 0) NScenar <- length(Scenarios)
                                  NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                  Params <- if (i <= length(Scenarios)) Params1 else Params2
                                  cat(paste0("Scenario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                  
                                  # Générer essais
                                  tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                          multinom_ttt = Scenarios[[NScenar]],
                                                                          rand_ratio = rep(1, NBras), seed = 121221)
                                  
                                  # Simuler le résultat
                                  if (Params$methode == "simoniva") {
                                    Resultat <- simu_simon(p_n = PN, p_a = PA,
                                                           tableau_essais = tableau_essais,
                                                           CaracSeuilSimon = SeuilSimon5,
                                                           CaracSeuilIva = SeuilIva5)
                                  } else {
                                    Resultat <- tryCatch(opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                                  p_n = PN, p_a = PA,
                                                                  CPar = SeuilBOP5[[1]][["C_"]], PPar = SeuilBOP5[[1]][["gamma"]],
                                                                  methode = Params$methode,
                                                                  A0_tox = Params$A0, SeuilP_tox = Params$SeuilP,
                                                                  A0_eff = Params$A0, SeuilP_eff = Params$SeuilP,
                                                                  a_tox = CBHM_tox5$a, b_tox = CBHM_tox5$b,
                                                                  a_eff = CBHM_eff5$a, b_eff = CBHM_eff5$b,
                                                                  p_mix_tox = rep(.5, NBras), p_mix_eff = rep(.5, NBras),
                                                                  tableau_essais = tableau_essais),
                                                         error = function(e) list(paste0(e)))
                                  }
                                  Resultat <- append(Resultat,
                                                     values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                     after = 2)
                                  return(Resultat)
                                  
                                }
    
    # Save results
    save(ResT, file = paste0("~/simu_priors/resultats_priorssens5norm_2026032026_", m, ".RData"))
    
  }
  
} else {
  cat("Decided not to do it (5 arms sensitivity, normalized power prior). Continuing with next simulation scenarios...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
}

stopCluster(cl)

# Inverted U-shaped scenarios ----

## 3 arms simulation ----

### Parameters for the simulation ----

Alpha <- .1
NSimu <- 5000
AnaEff <- rep(29, 2)
AnaTox <- rep(29, 2)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)
NBras <- 3

# Scénarios
Scenarios <- list(
  "ScI2" = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.19, 0.36, 0.11, 0.34), ttt3 = c(0.18, 0.34, 0.12, 0.36)),
  "Sc6"  = list(ttt1 = c(0.11, 0.19, 0.17, 0.53), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.21, 0.21, 0.39))
)
NomScenars <- names(Scenarios)

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
  "powBOP_efftox" = list(methode = "bop_power_efftox", A0 = .5, SeuilP = NA),
  "powBOP_tox" = list(methode = "bop_power_tox", A0 = .5, SeuilP = NA),
  "hBOP_efftox" = list(methode = "hier_efftox", A0 = NA, SeuilP = NA),
  "hBOP_tox" = list(methode = "hier_tox", A0 = NA, SeuilP = NA),
  "cbhmBOP_efftox" = list(methode = "cbhm_efftox", A0 = NA, SeuilP = NA),
  "cbhmBOP_tox" = list(methode = "cbhm_tox", A0 = NA, SeuilP = NA),
  "log1BOP_efftox" = list(methode = "bop_log1_efftox", A0 = NA, SeuilP = NA),
  "log1BOP_tox" = list(methode = "bop_log1_tox", A0 = NA, SeuilP = NA),
  "log2BOP_efftox" = list(methode = "bop_log2_efftox", A0 = NA, SeuilP = NA),
  "log2BOP_tox" = list(methode = "bop_log2_tox", A0 = NA, SeuilP = NA),
  "mBOP" = list(methode = "bop", efftox = 1, tox = 1, A0 = NA, SeuilP = NA),
  "Simon+Iva" = list(methode = "simoniva", A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)


###  Simulations ----

cl <- makeCluster(20)
registerDoParallel(cl)

if (TRUE) {
  cat("----\nMain simulation: 3 arms (non-monotonic scenarios)\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  for (m in seq_len(ceiling(length(Methodes) / 2))) {
    
    # Load 
    if (2 * m <= length(Methodes)) {
      cat(paste0("----\n", NomMethodes[2 * m - 1], " and ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    } else {
      cat(paste0("----\n", NomMethodes[2 * m - 1], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
    }
    
    # Parameters of the method
    Params1 <- Methodes[[2 * m - 1]]
    if (2 * m <= length(Methodes)) Params2 <- Methodes[[2 * m]]
    
    # Simulate the scenarios 
    VecScenarios <- if (2 * m <= length(Methodes)) seq_len(2 * length(Scenarios)) else seq_along(Scenarios)
    ResT <- foreach(i = VecScenarios,
                    .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                    .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "CBHM_eff",
                                "CBHM_tox", "gen_patients_multinom",
                                "m", "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                                "opcharac", "PA", "Params1", "PN", "real_essai_bayeslog_efftox",
                                "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                                "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                                "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                                "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                                "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                                "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                                "real_essai_modhier_efftox", "real_essai_modhier_tox", "simu_simon",
                                "Scenarios", "SeuilBOP", "SeuilSimon", "SeuilIva", "summarise_decision", "summarise_detect",
                                "summarise_ttt")) %dopar% {
                                  
                                  NScenar <- i %% length(Scenarios)
                                  if (NScenar == 0) NScenar <- length(Scenarios)
                                  NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                  Params <- if (i <= length(Scenarios)) Params1 else Params2
                                  cat(paste0("Scenario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                  
                                  # Generate trials
                                  tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                          multinom_ttt = Scenarios[[NScenar]],
                                                                          rand_ratio = rep(1, NBras), seed = 121221)
                                  
                                  # Simulate the result with the different methods
                                  if (Params$methode == "simoniva") {
                                    Resultat <- simu_simon(p_n = PN, p_a = PA,
                                                           tableau_essais = tableau_essais,
                                                           CaracSeuilSimon = SeuilSimon,
                                                           CaracSeuilIva = SeuilIva)
                                  } else {
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
                                  }
                                  
                                  Resultat <- append(Resultat,
                                                     values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                     after = 2)
                                  return(Resultat)
                                  
                                }
    
    # Save results
    save(ResT, file = paste0("~/simu_priors/resultats_priorsppalmonoton_20260316_", m, ".RData"))
    
  }
  
} else {
  cat("Decided not to do it (main analysis for non-monotonic scenarios). Continuing with next simulation scenarios...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
}

stopCluster(cl)

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

