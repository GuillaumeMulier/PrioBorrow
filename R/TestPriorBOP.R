# -------------------------------------------------------------------- #
# Script de simulation test pour l'inclusion de borrowing dans le BOP2 #
# Auteur : G. Mulier                                                   #
# Créé le 21/05/2024, modifié le 05/12/2024                            #
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
  "log6_tox" = CompilModLog_tox(SecondCov = TRUE),
  "crm_fixed_tox" = CompilModLogCrm_tox(fixed_intercept = TRUE),
  "crm_unfixed_tox" = CompilModLogCrm_tox(fixed_intercept = FALSE)
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
  "log6_efftox" = CompilModLog_efftox(SecondCov_eff = TRUE, SecondCov_tox = TRUE),
  "crm_fixed_efftox" = CompilModLogCrm_efftox(fixed_intercept = TRUE),
  "crm_unfixed_efftox" = CompilModLogCrm_efftox(fixed_intercept = FALSE)
)

# Parameters for CBHM
CBHM_eff <- CalibrateCBHM(NPts = rep(sum(AnaEff), NBras), Q0 = sum(PN[c(1, 2)]), Q1 = sum(PA[c(1, 2)]))
CBHM_tox <- CalibrateCBHM(NPts = rep(sum(AnaTox), NBras), Q0 = sum(PN[c(1, 3)]), Q1 = sum(PA[c(1, 3)]))

# Skeletons for CRM
SkelToxPPal <- MakeSkeleton(c(.3, .35, .4), A0 = 3, B0 = 1, MTD = NA, Delta = NULL)
SkelToxFaible <- MakeSkeleton(c(.25, .27, .3), A0 = 3, B0 = 1, MTD = NA, Delta = NULL)
SkelToxForte <- MakeSkeleton(c(.3, .35, .4), A0 = 3, B0 = 1, MTD = NA, Delta = NULL)
SkelToxDel <- MakeSkeleton(c(.3, .35, .4), A0 = 3, B0 = 1, MTD = 3, Delta = .1) # http://www.columbia.edu/~yc632/pub/crmcal.pdf
SkelEffPpal <- MakeSkeleton(c(.3, .4, .5), A0 = 3, B0 = 1, MTD = 3, Delta = NULL)
SkelEffFaible <- MakeSkeleton(c(.25, .27, .3), A0 = 3, B0 = 1, MTD = NA, Delta = NULL)
SkelEffForte <- MakeSkeleton(c(.5, .55, .6), A0 = 3, B0 = 1, MTD = NA, Delta = NULL)
SkelEffDel <- MakeSkeleton(c(.3, .4, .5), A0 = 3, B0 = 1, MTD = 3, Delta = .15)

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
  "powBOP_tox" = list(methode = "bop_power_tox", A0 = .5, SeuilP = NA),
  "powBOP_efftox" = list(methode = "bop_power_efftox", A0 = .5, SeuilP = NA),
  # "tBOP a=.5 s=.1_tox" = list(methode = "bop_power_test_tox", A0 = .5, SeuilP = .1),
  # "tBOP a=.5 s=.5_tox" = list(methode = "bop_power_test_tox", A0 = .5, SeuilP = .5),
  # "tBOP a=.5 s=.1_efftox" = list(methode = "bop_power_test_efftox", A0 = .5, SeuilP = .1),
  # "tBOP a=.5 s=.5_efftox" = list(methode = "bop_power_test_efftox", A0 = .5, SeuilP = .5),
  "hBOP_tox" = list(methode = "hier_tox", A0 = NA, SeuilP = NA),
  "hBOP_efftox" = list(methode = "hier_efftox", A0 = NA, SeuilP = NA),
  "cbhmBOP_tox" = list(methode = "cbhm_tox", A0 = NA, SeuilP = NA),
  "cbhmBOP_efftox" = list(methode = "cbhm_efftox", A0 = NA, SeuilP = NA),
  # "exnexBOP_tox" = list(methode = "exnex_tox", A0 = NA, SeuilP = NA),
  # "exnexBOP_efftox" = list(methode = "exnex_efftox", A0 = NA, SeuilP = NA),
  "bop_log1_tox" = list(methode = "bop_log1_tox", A0 = NA, SeuilP = NA),
  "bop_log2_tox" = list(methode = "bop_log2_tox", A0 = NA, SeuilP = NA),
  # "bop_log5_tox" = list(methode = "bop_log5_tox", A0 = NA, SeuilP = NA),
  # "bop_log6_tox" = list(methode = "bop_log6_tox", A0 = NA, SeuilP = NA),
  "bop_log1_efftox" = list(methode = "bop_log1_efftox", A0 = NA, SeuilP = NA),
  "bop_log2_efftox" = list(methode = "bop_log2_efftox", A0 = NA, SeuilP = NA),
  # "bop_log5_efftox" = list(methode = "bop_log5_efftox", A0 = NA, SeuilP = NA),
  # "bop_log6_efftox" = list(methode = "bop_log6_efftox", A0 = NA, SeuilP = NA),
  "mBOP" = list(methode = "bop", efftox = 1, tox = 1, A0 = NA, SeuilP = NA),
  # After tests on 3 scenarios, crm_unfixed is the most promising and thus we decided to explore it on different skeletons on all scenarios
  # "crmfixeda3_tox" = list(methode = "crm_tox", efftox = 1, tox = 1, A0 = 3, SeuilP = SkelTox),
  # "crmunfixed_tox" = list(methode = "crm_tox", efftox = 1, tox = 1, A0 = NA, SeuilP = SkelTox),
  # "crmfixeddela3_tox" = list(methode = "crm_tox", efftox = 0, tox = 1, A0 = 3, SeuilP = SkelToxDel),
  # "crmunfixeddel_tox" = list(methode = "crm_tox", efftox = 0, tox = 1, A0 = NA, SeuilP = SkelToxDel),
  # "crmfixeda3_efftox" = list(methode = "crm_efftox", efftox = 1, tox = 1, A0 = 3, SeuilP = list(SkelTox, SkelEff)),
  # "crmunfixed_efftox" = list(methode = "crm_efftox", efftox = 1, tox = 1, A0 = NA, SeuilP = list(SkelTox, SkelEff)),
  # "crmfixeddela3_efftox" = list(methode = "crm_efftox", efftox = 0, tox = 1, A0 = 3, SeuilP = list(SkelToxDel, SkelEffDel)),
  # "crmunfixeddel_efftox" = list(methode = "crm_efftox", efftox = 0, tox = 1, A0 = NA, SeuilP = list(SkelToxDel, SkelEffDel))
  "crmbopppal_tox" = list(methode = "crm_tox", efftox = 0, tox = 1, A0 = NA, SeuilP = SkelToxPPal),
  "crmbopcons_tox" = list(methode = "crm_tox", efftox = 0, tox = 1, A0 = NA, SeuilP = SkelToxForte),
  "crmbopopt_tox" = list(methode = "crm_tox", efftox = 0, tox = 1, A0 = NA, SeuilP = SkelToxFaible),
  "crmbopppal_efftox" = list(methode = "crm_efftox", efftox = 1, tox = 1, A0 = NA, SeuilP = list(SkelToxPPal, SkelEffPpal)),
  "crmbopcons_efftox" = list(methode = "crm_efftox", efftox = 1, tox = 1, A0 = NA, SeuilP = list(SkelToxForte, SkelEffFaible)),
  "crmbopopt_efftox" = list(methode = "crm_efftox", efftox = 1, tox = 1, A0 = NA, SeuilP = list(SkelToxFaible, SkelEffForte))
)
NomMethodes <- names(Methodes)

# Simulations
cl <- makeCluster(20)
registerDoParallel(cl)

if (FALSE) {
  cat("----\nSimulations à 3 bras\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  
  # for (m in c(1)) {
  for (m in seq_len(ceiling(length(Methodes) / 2))) {
    
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
                                "summarise_ttt", "SkelToxPPal", "SkelToxFaible", "SkelToxForte", 
                                "SkelEffPpal", "SkelEffFaible", "SkelEffForte")) %dopar% {
                                  
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
    save(ResT, file = paste0("~/simu_priors/resultats_priors_20241205_", m, ".RData"))
    
  }
  
}

if (TRUE) {
  # Additional simulations for logistic regression with crm skeleton
  cat("----\nSimulations à 3 bras avec squelette CRM\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
  VecScenars <- seq_along(Scenarios)
  VecMethodes <- c("crmppal_tox", "crmcons_tox", "crmopt_tox", "crmppal_efftox", "crmcons_efftox", "crmopt_efftox")
  GrilleSimu <- expand.grid(sc = VecScenars, meth = VecMethodes)
  GrilleSimu$meth <- as.character(GrilleSimu$meth)
  ResT <- foreach(i = seq_len(nrow(GrilleSimu)),
                  .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                  .export = c("%nin%", "Alpha", "AnaEff", "AnaTox", "CBHM_eff",
                              "CBHM_tox", "gen_patients_multinom",
                              "Methodes", "NBras", "NomMethodes", "NomScenars", "NSimu",
                              "opcharac", "PA", "PN", "real_essai_bayeslog_efftox",
                              "real_essai_bayeslog_tox", "real_essai_bop", "real_essai_bop_borrow",
                              "real_essai_bop_borrow_test_efftox", "real_essai_bop_borrow_test_tox",
                              "real_essai_bop_power_efftox", "real_essai_bop_power_test_efftox",
                              "real_essai_bop_power_test_tox", "real_essai_bop_power_tox",
                              "real_essai_bop_seq_efftox", "real_essai_bop_seq_tox", "real_essai_modcbhm_efftox",
                              "real_essai_modcbhm_tox", "real_essai_modexnex_efftox", "real_essai_modexnex_tox",
                              "real_essai_modhier_efftox", "real_essai_modhier_tox",
                              "Scenarios", "SeuilBOP", "summarise_decision", "summarise_detect",
                              "summarise_ttt", "SkelToxPpal", "SkelToxF", "SkelEff", "SkelEffDel")) %dopar% {
                                
                                NScenar <- GrilleSimu$sc[i]
                                Params <- Methodes[[GrilleSimu$meth[i]]]
                                cat(paste0("Scénario n°", NScenar, "/", nrow(GrilleSimu), " : ", NomScenars[NScenar], " with ", GrilleSimu$meth[i], "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                                
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
                                                   values = list(data.frame(methode = GrilleSimu$meth[i], scenar = NomScenars[NScenar])),
                                                   after = 2)
                                return(Resultat)
                                
                              }
  
  # Save results
  save(ResT, file = paste0("~/simu_priors/resultats_priors_crm_20250110.RData"))
}

stopCluster(cl)

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

