# -------------------------------------------------------------------- #
# Script de simulation test pour l'inclusion de borrowing dans le BOP2 #
# Auteur : G. Mulier                                                   #
# Créé le 21/05/2024, modifié le 08/07/2024                            #
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
AnaEff <- rep(15, 3)
AnaTox <- rep(15, 3)
PN <- c(.15, .15, .25, .45)
PA <- c(.2, .3, .1, .4)

Methodes <- list(
  "mBOP" = list(methode = "bop", A0 = NA, SeuilP = NA),
  "seqBOP" = list(methode = "bop_seq", A0 = NA, SeuilP = NA),
  "tBOP a=.5 s=.05" = list(methode = "bop_power_test", A0 = .5, SeuilP = .05),
  "tBOP a=.5 s=.1" = list(methode = "bop_power_test", A0 = .5, SeuilP = .1),
  "tBOP a=.5 s=.25" = list(methode = "bop_power_test", A0 = .5, SeuilP = .25),
  "tBOP a=.5 s=.5" = list(methode = "bop_power_test", A0 = .5, SeuilP = .5),
  "tBOP a=.25 s=.05" = list(methode = "bop_power_test", A0 = .25, SeuilP = .05),
  "tBOP a=.25 s=.1" = list(methode = "bop_power_test", A0 = .25, SeuilP = .1),
  "tBOP a=.25 s=.25" = list(methode = "bop_power_test", A0 = .25, SeuilP = .25),
  "tBOP a=.25 s=.5" = list(methode = "bop_power_test", A0 = .25, SeuilP = .5),
  "hBOP" = list(methode = "hier_tox", A0 = NA, SeuilP = NA),
  "log4" = list(methode = "bop_log4", A0 = NA, SeuilP = NA),
  "log5" = list(methode = "bop_log5", A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)

# 3 bras comme dans la figure ----

NBras <- 3

# Scénarios
Scenarios <- list(
  Sc1 = list(c(0.15, 0.15, 0.25, 0.45), c(0.15, 0.15, 0.25, 0.45), c(0.15, 0.15, 0.25, 0.45)),
  Sc2 = list(c(0.20, 0.30, 0.10, 0.40), c(0.20, 0.30, 0.10, 0.40), c(0.20, 0.30, 0.10, 0.40)),
  Sc3 = list(c(0.10, 0.20, 0.15, 0.55), c(0.20, 0.30, 0.10, 0.40), c(0.25, 0.35, 0.10, 0.30)),
  Sc4 = list(c(0.15, 0.25, 0.10, 0.50), c(0.22, 0.38, 0.08, 0.32), c(0.25, 0.35, 0.10, 0.30)),
  Sc5 = list(c(0.18, 0.22, 0.17, 0.43), c(0.25, 0.25, 0.15, 0.35), c(0.27, 0.23, 0.18, 0.32)),
  Sc6 = list(c(0.20, 0.30, 0.15, 0.35), c(0.25, 0.30, 0.15, 0.30), c(0.30, 0.30, 0.15, 0.25)),
  Sc7 = list(c(0.20, 0.30, 0.12, 0.38), c(0.22, 0.28, 0.13, 0.37), c(0.24, 0.26, 0.12, 0.38))
)
if (FALSE) {
  library(tidyr)
  library(ggplot2)
  Graphe <- do.call("rbind", lapply(names(Scenarios), \(x) {
    do.call("rbind", lapply(seq_along(Scenarios[[x]]), \(y) {
      data.frame(
        scenar = x,
        dose = y,
        'Efficacité' = Scenarios[[x]][[y]][1] + Scenarios[[x]][[y]][2],
        'Toxicité' = Scenarios[[x]][[y]][1] + Scenarios[[x]][[y]][3]
      )
    }))
  })) %>% 
    pivot_longer(-c(1, 2)) %>% 
    mutate(scenar = factor(scenar, levels = paste0("Sc", 1:10))) %>% 
    ggplot(aes(dose, value, color = name, shape = name)) +
    geom_point() +
    geom_line() +
    facet_wrap(vars(scenar), nrow = 2) +
    theme_light(base_size = 14) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
          strip.text = element_text(color = "black")) +
    scale_x_continuous(breaks = 1:3) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_discrete(type = c("steelblue", "darkred")) +
    labs(x = "Dose", y = "Proportion", color = NULL, shape = NULL)
  ggsave(Graphe, filename = "Figures/scenar_simul.png", device = "png", height = 6, width = 8)
}
NomScenars <- names(Scenarios)
GrilleSimu <- expand.grid(
  scenar = seq_along(Scenarios),
  methode = seq_along(Methodes)
)

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


# Simulations
cl <- makeCluster(21)
registerDoParallel(cl)

cat("----\nSimulations à 3 bras\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
Res1 <- foreach(i = seq_len(nrow(GrilleSimu)),
                .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                .export = c("%nin%", "AnaEff", "AnaTox", "GrilleSimu", "Methodes", "NBras", 
                            "NomMethodes", "NomScenars", "NSimu", "PA", "PN", 
                            "real_essai_bop", "real_essai_bop_borrow", 
                            "real_essai_bop_borrow_test", "real_essai_bop_power", "real_essai_bop_power_test",
                            "real_essai_bayeslog", "real_essai_bop_seq", "real_essai_modhier",
                            "CompilHier", "CompilLog1", "CompilLog2", "CompilLog3", "CompilLog4", "CompilLog5",
                            "Scenarios", "SeuilBOP", "summarise_decision", 
                            "summarise_detect", "summarise_ttt", "gen_patients_multinom")) %dopar% {
                              
                              cat(paste0("Scénario n°", i, "/", nrow(GrilleSimu), " : ", NomScenars[GrilleSimu$scenar[i]], " with ", Methodes[[GrilleSimu$methode[i]]]$methode, "\n"), file = "~/simu_priors/log.txt", append = TRUE)
                              
                              # Générer essais
                              tableau_essais <- gen_patients_multinom(NSimu, AnaEff, AnaTox, 
                                                                      multinom_ttt = Scenarios[[GrilleSimu$scenar[i]]], 
                                                                      rand_ratio = rep(1, NBras), seed = 121221)
                              
                              # Simuler le résultat
                              Params <- Methodes[[GrilleSimu$methode[i]]]
                              Resultat <- opcharac(ana_inter = AnaEff, ana_inter_tox = AnaTox,
                                                   p_n = PN, p_a = PA,
                                                   CPar = SeuilBOP[[1]][["C_"]], PPar = SeuilBOP[[1]][["gamma"]],
                                                   methode = Params$methode, A0 = Params$A0, SeuilP = Params$SeuilP,
                                                   tableau_essais = tableau_essais)
                              Resultat <- append(Resultat, 
                                                 values = list(data.frame(methode = NomMethodes[GrilleSimu$methode[i]], scenar = NomScenars[GrilleSimu$scenar[i]])),
                                                 after = 2)
                              return(Resultat)
                              
                            }

stopCluster(cl)

# Sauvegarder les résultats ----

save(Res1, file = "~/simu_priors/res_test_priors_v3.RData")

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

