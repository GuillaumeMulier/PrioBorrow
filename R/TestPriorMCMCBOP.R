# -------------------------------------------------------------------- #
# Script de simulation test pour l'inclusion de borrowing dans le BOP2 #
# Auteur : G. Mulier                                                   #
# Créé le 31/03/2025, modifié le 31/03/2025                            #
# -------------------------------------------------------------------- #

# devtools::load_all("~/pkgs/multibrasBOP2/")
cat("Charging packages and models...\n\n", file = "~/simu_priors/log.txt", append = TRUE)
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
  Sc1  = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
  Sc2  = list(ttt1 = c(0.13, 0.12, 0.27, 0.48), ttt2 = c(0.15, 0.13, 0.27, 0.45), ttt3 = c(0.16, 0.14, 0.29, 0.41)),
  Sc3  = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
  Sc4  = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.17, 0.35, 0.11, 0.37), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  Sc5  = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  Sc6  = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.18, 0.34, 0.12, 0.36), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  Sc7  = list(ttt1 = c(0.11, 0.19, 0.17, 0.53), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  Sc8  = list(ttt1 = c(0.14, 0.26, 0.14, 0.46), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  Sc9  = list(ttt1 = c(0.18, 0.32, 0.12, 0.38), ttt2 = c(0.22, 0.28, 0.15, 0.35), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
  Sc10 = list(ttt1 = c(0.12, 0.18, 0.18, 0.52), ttt2 = c(0.17, 0.23, 0.18, 0.42), ttt3 = c(0.23, 0.27, 0.17, 0.33))
)
NomScenars <- names(Scenarios)

# Precompile stan models 
cat("----\nCompiling efficacy/toxicity STAN\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
CompiledModelsEffTox <- list(
  "hBOP_h0_peuinf_var1_NA" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = 1,
                                              moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = 1),
  "hBOP_h1_peuinf_var1_NA" = CompilBHM_efftox(moy_mu_eff = log(.5 / .5), sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = 1,
                                              moy_mu_tox = log(.3 / .7), sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = 1),
  "hBOP_50p_peuinf_var1_NA" = CompilBHM_efftox(moy_mu_eff = 0, sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = 1,
                                               moy_mu_tox = 0, sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = 1),
  "hBOP_h0_inf_var1_NA" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 1.94, moy_sig_eff = 0, sigma_sig_eff = 1,
                                           moy_mu_tox = log(.4 / .6), sigma_mu_tox = 1.78, moy_sig_tox = 0, sigma_sig_tox = 1),
  "hBOP_h0_noninf_var1_NA" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 10, moy_sig_eff = 0, sigma_sig_eff = 1,
                                              moy_mu_tox = log(.4 / .6), sigma_mu_tox = 10, moy_sig_tox = 0, sigma_sig_tox = 1),
  "hBOP_h0_peuinf_var5_NA" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = 5,
                                              moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = 5),
  "hBOP_h0_peuinf_var0p5_NA" = CompilBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5, moy_sig_eff = 0, sigma_sig_eff = .5,
                                                moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5, moy_sig_tox = 0, sigma_sig_tox = .5),
  "cbhmBOP_h0_peuinf_NA_NA" = CompilCBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 2.5,
                                                moy_mu_tox = log(.4 / .6), sigma_mu_tox = 2.5),
  "cbhmBOP_h1_peuinf_NA_NA" = CompilCBHM_efftox(moy_mu_eff = log(.5 / .5), sigma_mu_eff = 2.5,
                                                moy_mu_tox = log(.3 / .7), sigma_mu_tox = 2.5),
  "cbhmBOP_50p_peuinf_NA_NA" = CompilCBHM_efftox(moy_mu_eff = 0, sigma_mu_eff = 2.5,
                                                 moy_mu_tox = 0, sigma_mu_tox = 2.5),
  "cbhmBOP_h0_inf_NA_NA" = CompilCBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 1.94,
                                             moy_mu_tox = log(.4 / .6), sigma_mu_tox = 1.78),
  "cbhmBOP_h0_noninf_NA_NA" = CompilCBHM_efftox(moy_mu_eff = log(.3 / .7), sigma_mu_eff = 10,
                                                moy_mu_tox = log(.4 / .6), sigma_mu_tox = 10),
  "logBOP_h0_inf_crois_inf" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 1, mu_coef_eff = 0.42, sigma_coef_eff = 1,
                                                  mu_inter_tox = log(.4 / .6), sigma_inter_tox = 1, mu_coef_tox = 0.22, sigma_coef_tox = 1),
  "logBOP_h0_peuinf_crois_inf" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 1,
                                                     mu_inter_tox = log(.4 / .6), sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 1),
  "logBOP_h1_inf_crois_inf" = CompilModLog_efftox(mu_inter_eff = log(.5 / .5), sigma_inter_eff = 1, mu_coef_eff = 0.42, sigma_coef_eff = 1,
                                                  mu_inter_tox = log(.3 / .7), sigma_inter_tox = 1, mu_coef_tox = 0.22, sigma_coef_tox = 1),
  "logBOP_h1_peuinf_crois_inf" = CompilModLog_efftox(mu_inter_eff = log(.5 / .5), sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 1,
                                                     mu_inter_tox = log(.3 / .7), sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 1),
  "logBOP_jef_jef_crois_inf" = CompilModLog_efftox(mu_inter_eff = 0, sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 1,
                                                   mu_inter_tox = 0, sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 1),
  "logBOP_jef_jef_trescrois_inf" = CompilModLog_efftox(mu_inter_eff = 0, sigma_inter_eff = 2.5, mu_coef_eff = 2, sigma_coef_eff = 1,
                                                       mu_inter_tox = 0, sigma_inter_tox = 2.5, mu_coef_tox = 2, sigma_coef_tox = 1),
  "logBOP_jef_jef_crois_peuinf" = CompilModLog_efftox(mu_inter_eff = 0, sigma_inter_eff = 2.5, mu_coef_eff = 0.42, sigma_coef_eff = 2.5,
                                                      mu_inter_tox = 0, sigma_inter_tox = 2.5, mu_coef_tox = 0.22, sigma_coef_tox = 2.5),
  "logBOP_h0_inf_jef_jef" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 1, mu_coef_eff = 0, sigma_coef_eff = 2.5,
                                                mu_inter_tox = log(.4 / .6), sigma_inter_tox = 1, mu_coef_tox = 0, sigma_coef_tox = 2.5),
  "logBOP_h0_inf_jef_inf" = CompilModLog_efftox(mu_inter_eff = log(.3 / .7), sigma_inter_eff = 1, mu_coef_eff = 0, sigma_coef_eff = 1,
                                                mu_inter_tox = log(.4 / .6), sigma_inter_tox = 1, mu_coef_tox = 0, sigma_coef_tox = 1)
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
NomMethodes <- names(CompiledModelsEffTox)

# Simulations
cl <- makeCluster(20)
registerDoParallel(cl)

cat("----\nSimulations des essais\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)

# for (m in c(1)) {
for (m in seq_len(ceiling(length(NomMethodes) / 2))) {
  
  # Load 
  if (2 * m <= length(NomMethodes)) {
    cat(paste0("----\n", NomMethodes[2 * m - 1], " et ", NomMethodes[2 * m], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
  } else {
    cat(paste0("----\n", NomMethodes[2 * m - 1], "\n----\n\n"), file = "~/simu_priors/log.txt", append = TRUE)
  }
  
  # Parameters of the method
  Modele1 <- NomMethodes[[2 * m - 1]]
  if (2 * m <= length(NomMethodes)) Modele2 <- NomMethodes[[2 * m]]
  
  # Simulate the scenarios 
  VecScenarios <- if (2 * m <= length(NomMethodes)) seq_len(2 * length(Scenarios)) else seq_along(Scenarios)
  ResT <- foreach(i = VecScenarios,
                  .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                  .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "CBHM_eff",
                              "CBHM_tox", "gen_patients_multinom",
                              "m", "NBras", "NomMethodes", "NomScenars", "NSimu",
                              "opcharac", "PA", "Modele1", "Modele2", "PN", "real_essai_bayeslog_efftox",
                              "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                              "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                              "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                              "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                              "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                              "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                              "real_essai_modhier_efftox", "real_essai_modhier_tox",
                              "Scenarios", "SeuilBOP", "summarise_decision", "summarise_detect",
                              "summarise_ttt", "CompiledModelsEffTox")) %dopar% {
                                
                                NScenar <- i %% length(Scenarios)
                                if (NScenar == 0) NScenar <- length(Scenarios)
                                NumMethode <- if (i <= length(Scenarios)) 2 * m - 1 else 2 * m
                                Modele <- if (i <= length(Scenarios)) Modele1 else Modele2
                                cat(paste0("Scénario n°", NScenar, "/", length(Scenarios), " : ", NomScenars[NScenar], " with ", NomMethodes[NumMethode], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                
                                # Générer essais
                                tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox,
                                                                        multinom_ttt = Scenarios[[NScenar]],
                                                                        rand_ratio = rep(1, NBras), seed = 121221)
                                
                                # Simuler le résultat
                                if (grepl("^hBOP", Modele)) { # Pas trop envie de changer mes fonction donc je triche en attribuant le bon élément de la liste de modèle
                                  Methode <- "hier_efftox"
                                  CompiledModelsEffTox[["hBOP_efftox"]] <- CompiledModelsEffTox[[Modele]]
                                } else if (grepl("^cbhmBOP", Modele)) {
                                  Methode <- "cbhm_efftox"
                                  CompiledModelsEffTox[["cbhmBOP_efftox"]] <- CompiledModelsEffTox[[Modele]]
                                } else if (grepl("^logBOP", Modele)) {
                                  Methode <- "bop_log1_efftox"
                                  CompiledModelsEffTox[["log1_efftox"]] <- CompiledModelsEffTox[[Modele]]
                                }
                                Resultat <- tryCatch(opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                              p_n = PN, p_a = PA,
                                                              CPar = SeuilBOP[[1]][["C_"]], PPar = SeuilBOP[[1]][["gamma"]],
                                                              methode = Methode,
                                                              A0_tox = NA, SeuilP_tox = NA,
                                                              A0_eff = NA, SeuilP_eff = NA,
                                                              a_tox = CBHM_tox$a, b_tox = CBHM_tox$b,
                                                              a_eff = CBHM_eff$a, b_eff = CBHM_eff$b,
                                                              p_mix_tox = rep(.5, NBras), p_mix_eff = rep(.5, NBras),
                                                              tableau_essais = tableau_essais),
                                                     error = function(e) list(paste0(e)))
                                Resultat <- append(Resultat,
                                                   values = list(data.frame(methode = NomMethodes[NumMethode], scenar = NomScenars[NScenar])),
                                                   after = 2)
                                if (grepl("^hBOP", Modele)) {
                                  CompiledModelsEffTox[["hBOP_efftox"]] <- NULL
                                } else if (grepl("cbhmBOP", Modele)) {
                                  CompiledModelsEffTox[["cbhmBOP_efftox"]] <- NULL
                                } else if (grepl("^logBOP", Modele)) {
                                  CompiledModelsEffTox[["log1_efftox"]] <- NULL
                                }
                                return(Resultat)
                                
                              }
  
  # Save results
  save(ResT, file = paste0("~/simu_priors/resultats_priors_20250331_", m, ".RData"))
  
}

stopCluster(cl)

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

