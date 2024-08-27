# -------------------------------------------------------------------- #
# Script de simulation test pour l'inclusion de borrowing dans le BOP2 #
# Auteur : G. Mulier                                                   #
# Créé le 21/05/2024, modifié le 05/07/2024                            #
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
AnaEff <- rep(20, 4)
AnaTox <- rep(20, 4)
PN <- c(.08, .07, .22, .63)
PA <- c(.15, .2, .05, .6)

Methodes <- list(
  "mBOP" = list(methode = "bop", A0 = NA, SeuilP = NA),
  "mBOP borrow" = list(methode = "bop_borrow", A0 = NA, SeuilP = NA),
  "pBOP a=.2" = list(methode = "bop_power", A0 = .2, SeuilP = NA),
  "pBOP a=.4" = list(methode = "bop_power", A0 = .4, SeuilP = NA),
  "pBOP a=.6" = list(methode = "bop_power", A0 = .6, SeuilP = NA),
  "pBOP a=.8" = list(methode = "bop_power", A0 = .8, SeuilP = NA),
  "pBOP a=1" = list(methode = "bop_power", A0 = 1, SeuilP = NA),
  "tBOP a=1 s=.05" = list(methode = "bop_power_test", A0 = 1, SeuilP = .05),
  "tBOP a=1 s=.1" = list(methode = "bop_power_test", A0 = 1, SeuilP = .1),
  "tBOP a=1 s=.25" = list(methode = "bop_power_test", A0 = 1, SeuilP = .25),
  "tBOP a=1 s=.5" = list(methode = "bop_power_test", A0 = 1, SeuilP = .5),
  "tBOP a=.5 s=.05" = list(methode = "bop_power_test", A0 = .5, SeuilP = .05),
  "tBOP a=.5 s=.1" = list(methode = "bop_power_test", A0 = .5, SeuilP = .1),
  "tBOP a=.5 s=.25" = list(methode = "bop_power_test", A0 = .5, SeuilP = .25),
  "tBOP a=.5 s=.5" = list(methode = "bop_power_test", A0 = .5, SeuilP = .5),
  "tBOP borrow" = list(methode = "bop_borrow_test", A0 = 1, SeuilP = NA),
  "hBOP" = list(methode = "hier_tox", A0 = NA, SeuilP = NA),
  "log1" = list(methode = "bop_log1", A0 = NA, SeuilP = NA),
  "log2" = list(methode = "bop_log2", A0 = NA, SeuilP = NA),
  "log3" = list(methode = "bop_log3", A0 = NA, SeuilP = NA),
  "seqBOP" = list(methode = "bop_seq", A0 = NA, SeuilP = NA)
)
NomMethodes <- names(Methodes)

# 3 bras comme dans la figure ----

NBras <- 3

# Scénarios
Scenarios <- list(
  "Sc1" = list(PN, PN, PN),
  "Sc2" = list(PA, PA, PA),
  "Sc3" = list(c(.1, .15, .05, .7), c(.15, .2, .05, .6), c(.18, .27, .07, .48)),
  "Sc4" = list(c(.1, .15, .05, .7), c(.15, .2, .1, .55), c(.18, .27, .07, .48)),
  "Sc5" = list(c(.1, .1, .05, .75), c(.15, .2, .1, .55), c(.15, .2, .1, .55)),
  "Sc6" = list(c(.15, .1, .1, .65), c(.2, .15, .1, .55), c(.25, .2, .1, .45)),
  "Sc7" = list(c(.15, .1, .1, .65), c(.22, .18, .08, .52), c(.22, .18, .08, .52)),
  "Sc8" = list(c(.15, .2, .1, .55), c(.25, .2, .1, .45), c(.25, .3, .1, .35)),
  "Sc9" = list(c(.07, .03, .08, .82), c(.08, .07, .12, .73), c(.12, .08, .13, .67)),
  "Sc10" = list(c(.07, .03, .03, .87), c(.1, .05, .05, .8), c(.13, .07, .07, .73))
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
                           nsim_oc = NSimu,
                           ana_inter = AnaEff, ana_inter_tox = AnaTox, 
                           p_n = PN, p_a = PA,
                           power_seq = seq(.25, 1.5, .005), 
                           seed = 121221, affich_mat = "No", methode = 4L)[[1]]
} else {
  SeuilBOP <- list(c("C_" = .71, "gamma" = 1.5, "alpha_calc" = .0963, "puissance_calc" = .9258))
}


# Simulations
cl <- makeCluster(20)
registerDoParallel(cl)

cat("----\nSimulations à 3 bras\n----\n\n", file = "~/simu_priors/log.txt", append = TRUE)
Res1 <- foreach(i = seq_len(nrow(GrilleSimu)),
                .packages = c("dplyr", "purrr", "rlang", "stringr", "rstan"),
                .export = c("%nin%", "AnaEff", "AnaTox", "GrilleSimu", "Methodes", "NBras", 
                            "NomMethodes", "NomScenars", "NSimu", "PA", "PN", 
                            "real_essai_bop", "real_essai_bop_borrow", 
                            "real_essai_bop_borrow_test", "real_essai_bop_power", "real_essai_bop_power_test",
                            "real_essai_bayeslog", "real_essai_bop_seq", "real_essai_modhier",
                            "CompilHier", "CompilLog1", "CompilLog2", "CompilLog3",
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

save(Res1, file = "~/simu_priors/res_test_priors_v2.RData")

cat("\nFini !\n", file = "~/simu_priors/log.txt", append = TRUE)

