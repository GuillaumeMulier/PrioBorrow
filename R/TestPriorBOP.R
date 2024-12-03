# -------------------------------------------------------------------- #
# Script de simulation test pour l'inclusion de borrowing dans le BOP2 #
# Auteur : G. Mulier                                                   #
# Créé le 21/05/2024, modifié le 30/08/2024                            #
# -------------------------------------------------------------------- #

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

# Paramètres de simulation ----

Alpha <- .1
NSimu <- 5000
# NSimu <- 2
AnaEff <- rep(15, 3)
AnaTox <- rep(15, 3)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)


# 3 bras comme dans la figure ----

NBras <- 3

# Scénarios
Scenarios <- list(
  Sc1 = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
  Sc2 = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
  Sc3 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.19, 0.23, 0.16, 0.42), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
  Sc4 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.22, 0.18, 0.40), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
  Sc5 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.16, 0.26, 0.15, 0.43), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
  Sc6 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.24, 0.15, 0.43), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
  Sc7 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.22, 0.18, 0.40), ttt3 = c(0.22, 0.28, 0.16, 0.34)),
  Sc8 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.15, 0.40), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
  Sc9 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.18, 0.37), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
  Sc10 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.27, 0.13, 0.42), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
  Sc11 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.27, 0.15, 0.40), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
  Sc12 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.18, 0.37), ttt3 = c(0.22, 0.28, 0.16, 0.34))
)
NomScenars <- names(Scenarios)

# Precompile stan models 
cat("----\nCompiling toxicity STAN\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
CompiledModelsTox <- list(
  "hBOP_tox" = CompilBHM_tox(),
  "cbhmBOP_tox" = CompilCBHM_tox(),
  "exnexBOP_tox" = CompilEXNEX_tox(moy_muex = -.62, sigma_muex = 5, sigma_sigex = 5, mu_nex = -.62, sigma_nex = 3),
  "log1_tox" = CompilModLog_tox(),
  "log2_tox" = CompilModLog_tox(PentePos = TRUE),
  "log3_tox" = CompilModLog_tox(mu_coef = .5),
  "log4_tox" = CompilModLog_tox(SecondCov = TRUE),
  "log5_tox" = CompilModLog_tox(PentePos = TRUE, SecondCov = TRUE),
  "log6_tox" = CompilModLog_tox(SecondCov = TRUE)
)
cat("----\nCompiling efficacy/toxicity STAN\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
CompiledModelsEffTox <- list(
  "hBOP_efftox" = CompilBHM_efftox(),
  "cbhmBOP_efftox" = CompilCBHM_efftox(),
  "exnexBOP_efftox" = CompilEXNEX_efftox(moy_muex_tox = -.62, sigma_muex_tox = 5, sigma_sigex_tox = 5, mu_nex_tox = -.62, sigma_nex_tox = 3, moy_muex_eff = -.41, sigma_muex_eff = 5, sigma_sigex_eff = 5, mu_nex_eff = -.41, sigma_nex_eff = 3),
  "log1_efftox" = CompilModLog_efftox(),
  "log2_efftox" = CompilModLog_efftox(PentePos_eff = TRUE, PentePos_tox = TRUE),
  "log3_efftox" = CompilModLog_efftox(mu_coef_eff = .5, mu_coef_tox = .5),
  "log4_efftox" = CompilModLog_efftox(SecondCov_eff = TRUE, SecondCov_tox = TRUE),
  "log5_efftox" = CompilModLog_efftox(PentePos_eff = TRUE, SecondCov_eff = TRUE, PentePos_tox = TRUE, SecondCov_tox = TRUE),
  "log6_efftox" = CompilModLog_efftox(SecondCov_eff = TRUE, SecondCov_tox = TRUE)
)


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
                           power_seq = seq(.25, 1.5, .005), 
                           seed = 121221, affich_mat = "No", methode = 4L)
} else {
  SeuilBOP <- list(c("C_" = .78, "gamma" = 1.455, "alpha_calc" = .0695, "puissance_calc" = .7093))
}
# The explored methods
Methodes <- list(
  "seqBOP_tox" = list(methode = "bop_seq_tox", efftox = 0, tox = 1, A0 = NA, SeuilP = NA),
  "seqBOP_efftox" = list(methode = "bop_seq_efftox", efftox = 1, tox = 0, A0 = NA, SeuilP = NA),
  "tBOP a=.5 s=.1_tox" = list(methode = "bop_power_test_tox", A0 = .5, SeuilP = .1),
  "tBOP a=.5 s=.5_tox" = list(methode = "bop_power_test_tox", A0 = .5, SeuilP = .5),
  "tBOP a=.5 s=.1_efftox" = list(methode = "bop_power_test_efftox", A0 = .5, SeuilP = .1),
  "tBOP a=.5 s=.5_efftox" = list(methode = "bop_power_test_efftox", A0 = .5, SeuilP = .5),
  "hBOP_tox" = list(methode = "hier_tox", A0 = NA, SeuilP = NA),
  "hBOP_efftox" = list(methode = "hier_efftox", A0 = NA, SeuilP = NA),
  "cbhmBOP_tox" = list(methode = "cbhm_tox", A0 = NA, SeuilP = NA),
  "cbhmBOP_efftox" = list(methode = "cbhm_efftox", A0 = NA, SeuilP = NA),
  "exnexBOP_tox" = list(methode = "exnex_tox", A0 = NA, SeuilP = NA),
  "exnexBOP_efftox" = list(methode = "exnex_efftox", A0 = NA, SeuilP = NA),
  "bop_log1_tox" = list(methode = "bop_log1_tox", A0 = NA, SeuilP = NA),
  "bop_log2_tox" = list(methode = "bop_log2_tox", A0 = NA, SeuilP = NA),
  "bop_log5_tox" = list(methode = "bop_log5_tox", A0 = NA, SeuilP = NA),
  "bop_log6_tox" = list(methode = "bop_log6_tox", A0 = NA, SeuilP = NA),
  "bop_log1_efftox" = list(methode = "bop_log1_efftox", A0 = NA, SeuilP = NA),
  "bop_log2_efftox" = list(methode = "bop_log2_efftox", A0 = NA, SeuilP = NA),
  "bop_log5_efftox" = list(methode = "bop_log5_efftox", A0 = NA, SeuilP = NA),
  "bop_log6_efftox" = list(methode = "bop_log6_efftox", A0 = NA, SeuilP = NA),
  "mBOP" = list(methode = "bop", efftox = 1, tox = 1, A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)

# Simulations
cl <- makeCluster(24)
registerDoParallel(cl)

cat("----\nSimulations à 3 bras\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)

# for (m in seq_len(ceiling(length(Methodes) / 2))) {
for (m in c(1)) {
  
  # Load 
  if (2 * m <= length(Methodes)) {
    cat(paste0("----\n", NomMethodes[2 * m - 1], " et ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
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
                              "opcharac", "PA", "Params1", "Params2", "PN", "real_essai_bayeslog_efftox",
                              "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                              "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                              "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                              "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                              "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                              "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                              "real_essai_modhier_efftox", "real_essai_modhier_tox",
                              "Scenarios", "SeuilBOP", "summarise_decision", "summarise_detect",
                              "summarise_ttt")) %dopar% {
                                
                                NScenar <- i %% length(Scenarios)
                                if (NScenar == 0) NScenar <- length(Scenarios)
                                NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                Params <- if (i <= length(Scenarios)) Params1 else Params2
                                cat(paste0("Scénario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                
                                # Générer essais
                                tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                        multinom_ttt = Scenarios[[NScenar]],
                                                                        rand_ratio = rep(1, NBras), seed = 121221)
                                
                                # Simuler le résultat
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
  save(ResT, file = paste0("~/simu_priors/resultats_priors_20241007_", m, ".RData"))
       
}

stopCluster(cl)

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

