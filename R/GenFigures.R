# --------------------------- #
# Create illustrative figures #
# Author: G. Mulier           #
# Created 2024/07/05          #
# Modified 2025/11/06         #
# --------------------------- #

# Packages and helpers ----

library(tidyverse)
library(patchwork)
library(ggtext)
library(rlang)

theme_set(theme_light(base_size = 26) +
            theme(strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                  strip.text = element_text(face = "bold", color = "black", size = 18),
                  axis.text.x = element_markdown()))


# Data ----

CaracGlobales <- list()
CaracBras <- list()
CaracEssais <- list()
walk(1:6, \(fich_num) {
  Fichiers <- paste0("Data/SimuPpal20250611/resultats_priorsppal_20250611_", 1:6, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobales", append(CaracGlobales, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBras", append(CaracBras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssais", append(CaracEssais, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobales <- do.call("rbind", CaracGlobales)
CaracGlobales$methode[CaracGlobales$methode == "mBOP"] <- "mBOP_both"
CaracGlobales$methode[CaracGlobales$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracGlobales <- separate(CaracGlobales, methode, c("methode", "cible"), "_")
CaracBras <- do.call("rbind", CaracBras)
CaracBras$methode[CaracBras$methode == "mBOP"] <- "mBOP_both"
CaracBras$methode[CaracBras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracBras <- separate(CaracBras, methode, c("methode", "cible"), "_")
CaracEssais <- do.call("rbind", CaracEssais)
CaracEssais$methode[CaracEssais$methode == "mBOP"] <- "mBOP_both"
CaracEssais$methode[CaracEssais$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracEssais <- separate(CaracEssais, methode, c("methode", "cible"), "_")
CaracEssais$larg_ic_eff <- CaracEssais$icsup_eff - CaracEssais$icinf_eff
CaracEssais$larg_ic_tox <- CaracEssais$icsup_tox - CaracEssais$icinf_tox
CaracGlobales$n_bras <- "3 arms"
CaracBras$n_bras <- "3 arms"
CaracEssais$n_bras <- "3 arms"

CaracGlobales4Bras <- list()
CaracBras4Bras <- list()
CaracEssais4Bras <- list()
walk(1:6, \(fich_num) {
  Fichiers <- paste0("Data/SimuSensi20250611/resultats_priorssens4_20250611_", 1:6, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobales4Bras", append(CaracGlobales4Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBras4Bras", append(CaracBras4Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssais4Bras", append(CaracEssais4Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobales4Bras <- do.call("rbind", CaracGlobales4Bras)
CaracGlobales4Bras$methode[CaracGlobales4Bras$methode == "mBOP"] <- "mBOP_both"
CaracGlobales4Bras$methode[CaracGlobales4Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracGlobales4Bras <- separate(CaracGlobales4Bras, methode, c("methode", "cible"), "_")
CaracBras4Bras <- do.call("rbind", CaracBras4Bras)
CaracBras4Bras$methode[CaracBras4Bras$methode == "mBOP"] <- "mBOP_both"
CaracBras4Bras$methode[CaracBras4Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracBras4Bras <- separate(CaracBras4Bras, methode, c("methode", "cible"), "_")
CaracEssais4Bras <- do.call("rbind", CaracEssais4Bras)
CaracEssais4Bras$methode[CaracEssais4Bras$methode == "mBOP"] <- "mBOP_both"
CaracEssais4Bras$methode[CaracEssais4Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracEssais4Bras <- separate(CaracEssais4Bras, methode, c("methode", "cible"), "_")
CaracEssais4Bras$larg_ic_eff <- CaracEssais4Bras$icsup_eff - CaracEssais4Bras$icinf_eff
CaracEssais4Bras$larg_ic_tox <- CaracEssais4Bras$icsup_tox - CaracEssais4Bras$icinf_tox
CaracGlobales4Bras$n_bras <- "4 arms"
CaracBras4Bras$n_bras <- "4 arms"
CaracEssais4Bras$n_bras <- "4 arms"

CaracGlobales5Bras <- list()
CaracBras5Bras <- list()
CaracEssais5Bras <- list()
walk(1:6, \(fich_num) {
  Fichiers <- paste0("Data/SimuSensi20250611/resultats_priorssens5_20250611_", 1:6, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobales5Bras", append(CaracGlobales5Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBras5Bras", append(CaracBras5Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssais5Bras", append(CaracEssais5Bras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobales5Bras <- do.call("rbind", CaracGlobales5Bras)
CaracGlobales5Bras$methode[CaracGlobales5Bras$methode == "mBOP"] <- "mBOP_both"
CaracGlobales5Bras$methode[CaracGlobales5Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracGlobales5Bras <- separate(CaracGlobales5Bras, methode, c("methode", "cible"), "_")
CaracBras5Bras <- do.call("rbind", CaracBras5Bras)
CaracBras5Bras$methode[CaracBras5Bras$methode == "mBOP"] <- "mBOP_both"
CaracBras5Bras$methode[CaracBras5Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracBras5Bras <- separate(CaracBras5Bras, methode, c("methode", "cible"), "_")
CaracEssais5Bras <- do.call("rbind", CaracEssais5Bras)
CaracEssais5Bras$methode[CaracEssais5Bras$methode == "mBOP"] <- "mBOP_both"
CaracEssais5Bras$methode[CaracEssais5Bras$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracEssais5Bras <- separate(CaracEssais5Bras, methode, c("methode", "cible"), "_")
CaracEssais5Bras$larg_ic_eff <- CaracEssais5Bras$icsup_eff - CaracEssais5Bras$icinf_eff
CaracEssais5Bras$larg_ic_tox <- CaracEssais5Bras$icsup_tox - CaracEssais5Bras$icinf_tox
CaracGlobales5Bras$n_bras <- "5 arms"
CaracBras5Bras$n_bras <- "5 arms"
CaracEssais5Bras$n_bras <- "5 arms"

CaracGlobalesPriors <- list()
CaracBrasPriors <- list()
CaracEssaisPriors <- list()
walk(1:6, \(fich_num) {
  Fichiers <- paste0("Data/SimuSensi20250611/resultats_priorssensprio_20250611_", 1:6, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobalesPriors", append(CaracGlobalesPriors, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBrasPriors", append(CaracBrasPriors, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssaisPriors", append(CaracEssaisPriors, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobalesPriors <- do.call("rbind", CaracGlobalesPriors)
CaracBrasPriors <- do.call("rbind", CaracBrasPriors)
CaracEssaisPriors <- do.call("rbind", CaracEssaisPriors)
CaracEssaisPriors$larg_ic_eff <- CaracEssaisPriors$icsup_eff - CaracEssaisPriors$icinf_eff
CaracEssaisPriors$larg_ic_tox <- CaracEssaisPriors$icsup_tox - CaracEssaisPriors$icinf_tox

CaracGlobalesCrm <- list()
CaracBrasCrm <- list()
CaracEssaisCrm <- list()
walk(1:3, \(fich_num) {
  Fichiers <- paste0("Data/SimuSensi20250611/resultats_priorssenscrm_20250611_", 1:3, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobalesCrm", append(CaracGlobalesCrm, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBrasCrm", append(CaracBrasCrm, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssaisCrm", append(CaracEssaisCrm, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobalesCrm <- do.call("rbind", CaracGlobalesCrm)
CaracBrasCrm <- do.call("rbind", CaracBrasCrm)
CaracEssaisCrm <- do.call("rbind", CaracEssaisCrm)
CaracEssaisCrm$larg_ic_eff <- CaracEssaisCrm$icsup_eff - CaracEssaisCrm$icinf_eff
CaracEssaisCrm$larg_ic_tox <- CaracEssaisCrm$icsup_tox - CaracEssaisCrm$icinf_tox

CaracGlobalesLogit <- list()
CaracBrasLogit <- list()
CaracEssaisLogit <- list()
walk(1:2, \(fich_num) {
  Fichiers <- paste0("Data/SimuLogistique20250611/resultats_priorslog_20250611_", 1:2, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobalesLogit", append(CaracGlobalesLogit, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBrasLogit", append(CaracBrasLogit, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssaisLogit", append(CaracEssaisLogit, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobalesLogit <- do.call("rbind", CaracGlobalesLogit)
CaracGlobalesLogit$methode[CaracGlobalesLogit$methode == "mBOP"] <- "mBOP_both"
CaracGlobalesLogit$methode[CaracGlobalesLogit$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracGlobalesLogit <- separate(CaracGlobalesLogit, methode, c("methode", "cible"), "_")
CaracBrasLogit <- do.call("rbind", CaracBrasLogit)
CaracBrasLogit$methode[CaracBrasLogit$methode == "mBOP"] <- "mBOP_both"
CaracBrasLogit$methode[CaracBrasLogit$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracBrasLogit <- separate(CaracBrasLogit, methode, c("methode", "cible"), "_")
CaracEssaisLogit <- do.call("rbind", CaracEssaisLogit)
CaracEssaisLogit$methode[CaracEssaisLogit$methode == "mBOP"] <- "mBOP_both"
CaracEssaisLogit$methode[CaracEssaisLogit$methode == "Simon+Iva"] <- "Simon+TM_both"
CaracEssaisLogit <- separate(CaracEssaisLogit, methode, c("methode", "cible"), "_")
CaracEssaisLogit$larg_ic_eff <- CaracEssaisLogit$icsup_eff - CaracEssaisLogit$icinf_eff
CaracEssaisLogit$larg_ic_tox <- CaracEssaisLogit$icsup_tox - CaracEssaisLogit$icinf_tox
CaracGlobalesLogit$n_bras <- "3 arms"
CaracBrasLogit$n_bras <- "3 arms"
CaracEssaisLogit$n_bras <- "3 arms"
CaracGlobalesLogit$methode <- paste0("ver", CaracGlobalesLogit$methode)
CaracBrasLogit$methode <- paste0("ver", CaracBrasLogit$methode)
CaracEssaisLogit$methode <- paste0("ver", CaracEssaisLogit$methode)


# Scenarios ----

ListeScenars <- list(
  Sc1  = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
  Sc2  = list(ttt1 = c(0.13, 0.12, 0.27, 0.48), ttt2 = c(0.15, 0.13, 0.27, 0.45), ttt3 = c(0.16, 0.14, 0.29, 0.41)),
  Sc3  = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
  Sc4  = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.17, 0.35, 0.11, 0.37), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  Sc5  = list(ttt1 = c(0.11, 0.19, 0.17, 0.53), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  Sc6  = list(ttt1 = c(0.14, 0.26, 0.14, 0.46), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
  Sc7  = list(ttt1 = c(0.18, 0.32, 0.12, 0.38), ttt2 = c(0.22, 0.28, 0.15, 0.35), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
  Sc8  = list(ttt1 = c(0.12, 0.18, 0.18, 0.52), ttt2 = c(0.17, 0.23, 0.18, 0.42), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
  ScI1 = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  ScI2 = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.18, 0.34, 0.12, 0.36), ttt3 = c(0.25, 0.30, 0.15, 0.30))
)
NomsScenars <- c(
  "Sc1"  = "Sc1 : Global Null",
  "Sc2"  = "Sc2 : Global Null increasing",
  "Sc3"  = "Sc3 : Global Alternative",
  "Sc4"  = "Sc4 : Global Alternative increasing",
  "ScI1" = "ScI1 : 1st dosis futile",
  "ScI2" = "ScI2 : 3rd dosis toxic",
  "Sc5"  = "Sc5 : 1st dosis futile and 3rd dosis toxic",
  "Sc6"  = "Sc6 : 3rd dosis toxic and 1st dosis intermediate",
  "Sc7"  = "Sc7 : 2nd dosis intermediate and 3rd dosis toxic",
  "Sc8"  = "Sc8 : No promising dosis"
)

# Figure 1 ----

# Couleurs <- c("prometteur" = "green", "futile/toxique" = "red", "futile" = "coral", "toxique" = "darkred", "intermediaire" = "steelblue")
Graphe <- do.call("rbind", lapply(names(ListeScenars), \(l) {
  do.call("rbind", lapply(names(ListeScenars[[l]]), \(x) {
    data.frame(scenar = l, 
               bras = x,
               eff = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][2],
               tox = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][3])
  }))
})) %>% 
  # mutate(scenar = factor(scenar, levels = c(paste0("Sc", 1:8), "ScI1", "ScI2"), labels = NomsScenars[c(paste0("Sc", 1:8), "ScI1", "ScI2")]),
  mutate(scenar = factor(scenar, levels = c(paste0("Sc", 1:8), "ScI1", "ScI2")),
         efficacite = case_when(eff < .31 ~ "futile",
                                eff < .5 ~ "intermediaire",
                                TRUE ~ "efficace"),
         toxicite = case_when(tox < .31 ~ "non toxique",
                              tox < .4 ~ "intermediaire",
                              TRUE ~ "toxique"),
         couleur = case_when(efficacite == "efficace" & toxicite == "non toxique" ~ "Promising",
                             efficacite == "futile" & toxicite == "toxique" ~ "Stopping",
                             efficacite == "futile" ~ "Stopping",
                             toxicite == "toxique" ~ "Stopping",
                             TRUE ~ "Intermediate"),
         bras = gsub("ttt", "D", bras)) %>% 
  # mutate(couleur2 = Couleurs[couleur],
  #        bras = paste0("<span style='color:", couleur2, "'>", bras, "</span>"),
  #        bras = factor(bras, levels = c(unique(bras)[grepl("ttt1", unique(bras))], unique(bras)[grepl("ttt2", unique(bras))], unique(bras)[grepl("ttt3", unique(bras))]))) %>% 
  ggplot(aes(bras)) +
  geom_line(aes(y = eff, color = "Efficacy", group = scenar)) +
  geom_point(aes(y = eff, color = "Efficacy", shape = couleur), size = 3) +
  geom_line(aes(y = tox, color = "Toxicity", group = scenar)) +
  geom_point(aes(y = tox, color = "Toxicity", shape = couleur), size = 3) +
  facet_wrap(vars(scenar), ncol = 2, scales = "free_x") +
  expand_limits(y = 0) +
  theme_light(base_size = 18) +
  theme(strip.background = element_rect(fill = "white", color = "black", size = 1.2),
        strip.text = element_text(face = "bold", color = "black", size = 15),
        axis.text.x = element_markdown(),
        legend.position = "right") +
  scale_color_discrete(type = c("darkblue", "darkred")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 2))) +
  labs(x = "Treatment dose", y = "Probability", color = "Endpoint", shape = "Decision")
# caption = "Vert=prometteur / rouge=futile et toxique / beige=futile / rouge foncé=toxique / bleu=intermédiaire")
g <- ggplot_gtable(ggplot_build(Graphe))
stripes <- which(grepl("strip-t", g$layout$name))
for (i in stripes[1:2]) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <- "steelblue"
}
grid::grid.draw(g)
ggsave(g, filename = "Figures/scenar_simul_v5.png", height = 300, width = 230, units = "mm", dpi = 300)

# Figure 2 ----

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4"), methode %in% c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("Sc1", "Sc2", "Sc3", "Sc4")),
         methode = factor(methode, levels = rev(c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_sc1234.png", device = "png", 
       height = 27.33 / 2, width = 40.99 / 2, units = "cm", dpi = 300)


# Figure 3 ----

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8"), methode %in% c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("Sc5", "Sc6", "Sc7", "Sc8")),
         methode = factor(methode, levels = rev(c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_sc5678.png", device = "png", 
       height = 130, width = 170, units = "mm", dpi = 300)

# Figure 4 and table 2 ----

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("ScI1", "ScI2"), methode %in% c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("ScI1", "ScI2")),
         methode = factor(methode, levels = rev(c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_scibru.png", device = "png", 
       height = 100, width = 170, units = "mm", dpi = 300)

CaracEssais %>% 
  filter(cible != "tox", scenar == "ScI1") %>% 
  group_by(methode, n_simu) %>% 
  summarise(choix = ifelse(all(decision != "Accept the treatment"), "", paste(ttt[decision == "Accept the treatment"], collapse = "-"))) %>% 
  ungroup() %>% 
  mutate(prem_dose = gsub("^(ttt\\d)-.*$", "\\1", choix)) %>% 
  count(methode, prem_dose) %>% 
  mutate(pct = sprintf("%.1f%%", 100 * n / sum(n)), .by = methode) %>% 
  print(n = Inf)
CaracEssais %>% 
  filter(cible != "tox", scenar == "ScI2") %>% 
  group_by(methode, n_simu) %>% 
  summarise(choix = ifelse(all(decision != "Accept the treatment"), "", paste(ttt[decision == "Accept the treatment"], collapse = "-"))) %>% 
  ungroup() %>% 
  mutate(prem_dose = gsub("^(ttt\\d)-.*$", "\\1", choix)) %>% 
  count(methode, prem_dose) %>% 
  mutate(pct = sprintf("%.1f%%", 100 * n / sum(n)), .by = methode) %>% 
  print(n = Inf)
CaracEssais %>% 
  filter(cible != "tox", scenar %in% c("ScI1", "ScI2")) %>% 
  group_by(scenar, methode, n_simu) %>% 
  summarise(choix = ifelse(all(decision != "Accept the treatment"), "", paste(ttt[decision == "Accept the treatment"], collapse = "-"))) %>% 
  ungroup() %>% 
  mutate(verite = ifelse(scenar == "ScI1", "ttt2-ttt3", "ttt1-ttt2"),
         final = as.numeric(verite == choix)) %>% 
  count(scenar, methode, final) %>% 
  mutate(pct = sprintf("%.1f", 100 * n / sum(n)), .by = c(scenar, methode)) %>% 
  filter(final == 1) %>% 
  complete(scenar, methode, fill = list(pct = "0.0")) %>% 
  select(methode, scenar, pct) %>% 
  pivot_wider(names_from = "scenar", values_from = "pct") %>% 
  slice(c(5, 7, 6, 2, 1, 3, 4)) %>% 
  apply(1, \(x) paste(x, collapse = " & ")) %>% 
  paste(collapse = " \\\\ \n") %>% cat()

CaracBras %>% 
  filter(cible != "tox") %>% 
  mutate(Scenario = scenar, 
         Dose = gsub("^ttt", "D", ttt),
         arret_precoce = sprintf("%.1f%%", 100 * arret_precoce)) %>% 
  select(methode, Scenario, Dose, arret_precoce) %>% 
  pivot_wider(names_from = "methode", values_from = "arret_precoce") %>% 
  filter(Scenario %in% c("ScI1", "ScI2")) %>% 
  print(n = Inf)

# Table S1 ----

CaracBras %>% 
  filter(cible != "tox") %>% 
  mutate(Scenario = scenar, 
         Dose = gsub("^ttt", "Dose ", ttt),
         AP = sprintf("%.1f%%", 100 * arret_precoce),
         MPts = sprintf("%.1f", tot_pat)) %>% 
  select(methode, Scenario, Dose, AP, MPts) %>% 
  pivot_wider(names_from = "methode", values_from = c("AP", "MPts"), names_glue = "{methode}_{.value}") %>% 
  select(Scenario, Dose, mBOP_AP, mBOP_MPts, `Simon+TM_AP`, `Simon+TM_MPts`, powBOP_AP, powBOP_MPts, 
         hBOP_AP, hBOP_MPts, cbhmBOP_AP, cbhmBOP_MPts, log1BOP_AP, log1BOP_MPts, log2BOP_AP, log2BOP_MPts) %>% 
  xtable::xtable()

# Efficacy and toxicity estimations ----

TabScenars <- do.call("rbind", lapply(names(ListeScenars), \(nom_scenar) {
  do.call("rbind", lapply(names(ListeScenars[[nom_scenar]]), \(nom_bras) {
    VecProba <- ListeScenars[[nom_scenar]][[nom_bras]]
    data.frame("scenar" = nom_scenar, "ttt" = nom_bras, "eff_true" = VecProba[1] + VecProba[2], "tox_true" = VecProba[1] + VecProba[3]) %>% 
      mutate(ttt = gsub("^ttt", "D", ttt))
  }))
}))
Graphe <- ((CaracEssais %>% 
              filter(cible != "tox", scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4")) %>% 
              mutate(ttt = gsub("^ttt", "D", ttt),
                     methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
              ggplot(aes(y = est_eff, x = methode, fill = methode)) +
              geom_boxplot(alpha = .6) +
              geom_hline(data = TabScenars %>% filter(scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4")), aes(yintercept = eff_true), 
                         color = "darkred", linetype = "solid", linewidth = .7) +
              facet_grid(scenar ~ ttt) +
              scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
              theme_light(base_size = 13) +
              scale_y_continuous(labels = scales::percent_format()) +
              theme(legend.position = "none",
                    strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                    strip.text = element_text(face = "bold", color = "black", size = 12),
                    axis.text.x = element_markdown(angle = 45, hjust = 1)) +
              labs(x = NULL, y = "Efficacy")) /
             (CaracEssais %>% 
                filter(cible != "tox", scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4")) %>% 
                mutate(ttt = gsub("^ttt", "D", ttt),
                       methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
                ggplot(aes(y = est_tox, x = methode, fill = methode)) +
                geom_boxplot(alpha = .6) +
                geom_hline(data = TabScenars %>% filter(scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4")), aes(yintercept = tox_true), 
                           color = "darkred", linetype = "solid", linewidth = .7) +
                facet_grid(scenar ~ ttt) +
                scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
                theme_light(base_size = 13) +
                scale_y_continuous(labels = scales::percent_format()) +
                theme(legend.position = "none",
                      strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                      strip.text = element_text(face = "bold", color = "black", size = 12),
                      axis.text.x = element_markdown(angle = 45, hjust = 1)) +
                labs(x = NULL, y = "Toxicity"))) +
  plot_annotation(tag_levels = "A")
ggsave(Graphe, filename = "Figures/estimations_sc1234.png", device = "png", 
       height = 250, width = 150, units = "mm", dpi = 300)
Graphe <- ((CaracEssais %>% 
              filter(cible != "tox", scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8")) %>% 
              mutate(ttt = gsub("^ttt", "D", ttt),
                     methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
              ggplot(aes(y = est_eff, x = methode, fill = methode)) +
              geom_boxplot(alpha = .6) +
              geom_hline(data = TabScenars %>% filter(scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8")), aes(yintercept = eff_true), 
                         color = "darkred", linetype = "solid", linewidth = .7) +
              facet_grid(scenar ~ ttt) +
              scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
              theme_light(base_size = 13) +
              scale_y_continuous(labels = scales::percent_format()) +
              theme(legend.position = "none",
                    strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                    strip.text = element_text(face = "bold", color = "black", size = 12),
                    axis.text.x = element_markdown(angle = 45, hjust = 1)) +
              labs(x = NULL, y = "Efficacy")) /
             (CaracEssais %>% 
                filter(cible != "tox", scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8")) %>% 
                mutate(ttt = gsub("^ttt", "D", ttt),
                       methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
                ggplot(aes(y = est_tox, x = methode, fill = methode)) +
                geom_boxplot(alpha = .6) +
                geom_hline(data = TabScenars %>% filter(scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8")), aes(yintercept = tox_true), 
                           color = "darkred", linetype = "solid", linewidth = .7) +
                facet_grid(scenar ~ ttt) +
                scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
                theme_light(base_size = 13) +
                scale_y_continuous(labels = scales::percent_format()) +
                theme(legend.position = "none",
                      strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                      strip.text = element_text(face = "bold", color = "black", size = 12),
                      axis.text.x = element_markdown(angle = 45, hjust = 1)) +
                labs(x = NULL, y = "Toxicity"))) +
  plot_annotation(tag_levels = "A")
ggsave(Graphe, filename = "Figures/estimations_sc5678.png", device = "png", 
       height = 250, width = 150, units = "mm", dpi = 300)
Graphe <- ((CaracEssais %>% 
              filter(cible != "tox", scenar %in% c("ScI1", "ScI2")) %>% 
              mutate(ttt = gsub("^ttt", "D", ttt),
                     methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
              ggplot(aes(y = est_eff, x = methode, fill = methode)) +
              geom_boxplot(alpha = .6) +
              geom_hline(data = TabScenars %>% filter(scenar %in% c("ScI1", "ScI2")), aes(yintercept = eff_true), 
                         color = "darkred", linetype = "solid", linewidth = .7) +
              facet_grid(scenar ~ ttt) +
              scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
              theme_light(base_size = 13) +
              scale_y_continuous(labels = scales::percent_format()) +
              theme(legend.position = "none",
                    strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                    strip.text = element_text(face = "bold", color = "black", size = 12),
                    axis.text.x = element_markdown(angle = 45, hjust = 1)) +
              labs(x = NULL, y = "Efficacy")) /
             (CaracEssais %>% 
                filter(cible != "tox", scenar %in% c("ScI1", "ScI2")) %>% 
                mutate(ttt = gsub("^ttt", "D", ttt),
                       methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
                ggplot(aes(y = est_tox, x = methode, fill = methode)) +
                geom_boxplot(alpha = .6) +
                geom_hline(data = TabScenars %>% filter(scenar %in% c("ScI1", "ScI2")), aes(yintercept = tox_true), 
                           color = "darkred", linetype = "solid", linewidth = .7) +
                facet_grid(scenar ~ ttt) +
                scale_fill_discrete(type = c("#0f3d9b", "#08a448", "#7846ab", "#c33215", "#dd790e", "#11caa6", "#07772e")) + 
                theme_light(base_size = 13) +
                scale_y_continuous(labels = scales::percent_format()) +
                theme(legend.position = "none",
                      strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                      strip.text = element_text(face = "bold", color = "black", size = 12),
                      axis.text.x = element_markdown(angle = 45, hjust = 1)) +
                labs(x = NULL, y = "Toxicity"))) +
  plot_annotation(tag_levels = "A")
ggsave(Graphe, filename = "Figures/estimations_scibru.png", device = "png", 
       height = 250, width = 150, units = "mm", dpi = 300)

# Sensitivity analysis: 4 and 5 arms ----

Scenarios <- list(
  "Sc2"  = list(ttt1 = c(0.13, 0.12, 0.27, 0.48), ttt2 = c(0.15, 0.13, 0.27, 0.45), ttt3 = c(0.16, 0.14, 0.29, 0.41)),
  "Sc4"  = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.17, 0.35, 0.11, 0.37), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI1" = list(ttt1 = c(0.10, 0.20, 0.15, 0.55), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.15, 0.35, 0.10, 0.40), ttt2 = c(0.18, 0.34, 0.12, 0.36), ttt3 = c(0.25, 0.30, 0.15, 0.30))
)
TabScenars3 <- do.call("rbind", lapply(names(Scenarios), \(nom_scenar) {
  do.call("rbind", lapply(names(Scenarios[[nom_scenar]]), \(nom_bras) {
    VecProba <- Scenarios[[nom_scenar]][[nom_bras]]
    data.frame("nb_bras" = "3 arms", "scenar" = nom_scenar, "ttt" = nom_bras, "eff_true" = VecProba[1] + VecProba[2], "tox_true" = VecProba[1] + VecProba[3])
  }))
}))
Scenarios <- list(
  "Sc2"  = list(ttt1 = c(0.12, 0.12, 0.28, 0.48), ttt2 = c(0.13, 0.13, 0.29, 0.45), ttt3 = c(0.15, 0.13, 0.29, 0.43), ttt4 = c(0.17, 0.13, 0.33, 0.37)),
  "Sc4"  = list(ttt1 = c(0.15, 0.35, 0.09, 0.41), ttt2 = c(0.16, 0.36, 0.10, 0.38), ttt3 = c(0.18, 0.36, 0.10, 0.36), ttt4 = c(0.20, 0.36, 0.10, 0.34)),
  "ScI1" = list(ttt1 = c(0.07, 0.18, 0.13, 0.62), ttt2 = c(0.10, 0.20, 0.15, 0.55), ttt3 = c(0.18, 0.32, 0.12, 0.38), ttt4 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.12, 0.38, 0.08, 0.42), ttt2 = c(0.16, 0.36, 0.09, 0.39), ttt3 = c(0.19, 0.35, 0.11, 0.35), ttt4 = c(0.25, 0.31, 0.15, 0.29))
)
TabScenars4 <- do.call("rbind", lapply(names(Scenarios), \(nom_scenar) {
  do.call("rbind", lapply(names(Scenarios[[nom_scenar]]), \(nom_bras) {
    VecProba <- Scenarios[[nom_scenar]][[nom_bras]]
    data.frame("nb_bras" = "4 arms", "scenar" = nom_scenar, "ttt" = nom_bras, "eff_true" = VecProba[1] + VecProba[2], "tox_true" = VecProba[1] + VecProba[3])
  }))
}))
Scenarios <- list(
  "Sc2"  = list(ttt1 = c(0.11, 0.11, 0.29, 0.49), ttt2 = c(0.13, 0.11, 0.29, 0.47), ttt3 = c(0.14, 0.12, 0.30, 0.44), ttt4 = c(0.15, 0.13, 0.31, 0.41), ttt5 = c(0.17, 0.13, 0.31, 0.39)),
  "Sc4"  = list(ttt1 = c(0.13, 0.37, 0.09, 0.41), ttt2 = c(0.15, 0.36, 0.09, 0.40), ttt3 = c(0.16, 0.36, 0.10, 0.38), ttt4 = c(0.18, 0.35, 0.10, 0.37), ttt5 = c(0.19, 0.35, 0.11, 0.35)),
  "ScI1" = list(ttt1 = c(0.08, 0.17, 0.16, 0.59), ttt2 = c(0.10, 0.20, 0.16, 0.54), ttt3 = c(0.17, 0.33, 0.11, 0.39), ttt4 = c(0.19, 0.34, 0.11, 0.36), ttt5 = c(0.19, 0.36, 0.11, 0.34)),
  "ScI2" = list(ttt1 = c(0.12, 0.38, 0.08, 0.42), ttt2 = c(0.14, 0.37, 0.08, 0.41), ttt3 = c(0.16, 0.36, 0.09, 0.39), ttt4 = c(0.19, 0.34, 0.11, 0.36), ttt5 = c(0.25, 0.29, 0.15, 0.31))
)
TabScenars5 <- do.call("rbind", lapply(names(Scenarios), \(nom_scenar) {
  do.call("rbind", lapply(names(Scenarios[[nom_scenar]]), \(nom_bras) {
    VecProba <- Scenarios[[nom_scenar]][[nom_bras]]
    data.frame("nb_bras" = "5 arms", "scenar" = nom_scenar, "ttt" = nom_bras, "eff_true" = VecProba[1] + VecProba[2], "tox_true" = VecProba[1] + VecProba[3])
  }))
}))
TabScenarRegroup <- bind_rows(TabScenars3, TabScenars4, TabScenars5)
Graphe <- TabScenarRegroup %>% 
  mutate(scenar = factor(scenar, 
                         levels = c(paste0("Sc", c(2, 4)), "ScI1", "ScI2"))) %>%
  mutate(efficacite = case_when(eff_true < .31 ~ "futile",
                                eff_true < .5 ~ "intermediaire",
                                TRUE ~ "efficace"),
         toxicite = case_when(tox_true < .31 ~ "non toxique",
                              tox_true < .4 ~ "intermediaire",
                              TRUE ~ "toxique"),
         couleur = case_when(efficacite == "efficace" & toxicite == "non toxique" ~ "Promising",
                             efficacite == "futile" & toxicite == "toxique" ~ "Stopping",
                             efficacite == "futile" ~ "Stopping",
                             toxicite == "toxique" ~ "Stopping",
                             TRUE ~ "Intermediate"),
         bras = gsub("ttt", "D", ttt)) %>% 
  ggplot(aes(bras, group = nb_bras)) +
  geom_line(aes(y = eff_true, color = "Efficacy")) +
  geom_point(aes(y = eff_true, color = "Efficacy", shape = couleur), size = 3) +
  geom_line(aes(y = tox_true, color = "Toxicity")) +
  geom_point(aes(y = tox_true, color = "Toxicity", shape = couleur), size = 3) +
  facet_grid(nb_bras ~ scenar, scales = "free_x") +
  expand_limits(y = 0) +
  scale_color_discrete(type = c("darkblue", "darkred")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "Treatment dose", y = "Probability", color = "Endpoint", shape = "Decision")
ggsave(Graphe, filename = "Figures/scenar_sensi_arms.png", device = "png", height = 7, width = 14)

CaracBrasSensiBras <- rbind(CaracBras %>% filter(scenar %in% c("Sc2", "Sc4", "ScI1", "ScI2")), CaracBras4Bras, CaracBras5Bras)
Graphe <- CaracBrasSensiBras %>%
  filter(cible %in% c("efftox", "both"), scenar == "Sc2") %>%
  mutate(ttt = gsub("^ttt", "D", ttt),
         methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
  ggplot(aes(ttt, rejet_h0, group = n_bras, color = n_bras)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(vars(methode)) +
  labs(x = "Dose", y = "Proportion of conclusion to promising dose", color = "# arms") +
  scale_color_discrete(type = c("black", "darkblue", "darkred")) +
  expand_limits(x = 0) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(Graphe, filename = "Figures/resultbras_sensi_arms_2.png", device = "png", height = 8, width = 11)
Graphe <- CaracBrasSensiBras %>%
  filter(cible %in% c("efftox", "both"), scenar == "Sc4") %>%
  mutate(ttt = gsub("^ttt", "D", ttt),
         methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
  ggplot(aes(ttt, rejet_h0, group = n_bras, color = n_bras)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(vars(methode)) +
  labs(x = "Dose", y = "Proportion of conclusion to promising dose", color = "# arms") +
  scale_color_discrete(type = c("black", "darkblue", "darkred")) +
  expand_limits(x = 0) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(Graphe, filename = "Figures/resultbras_sensi_arms_4.png", device = "png", height = 8, width = 11)
Graphe <- CaracBrasSensiBras %>%
  filter(cible %in% c("efftox", "both"), scenar == "ScI1") %>%
  mutate(ttt = gsub("^ttt", "D", ttt),
         methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
  ggplot(aes(ttt, rejet_h0, group = n_bras, color = n_bras)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(vars(methode)) +
  labs(x = "Dose", y = "Proportion of conclusion to promising dose", color = "# arms") +
  scale_color_discrete(type = c("black", "darkblue", "darkred")) +
  expand_limits(x = 0) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(Graphe, filename = "Figures/resultbras_sensi_arms_I1.png", device = "png", height = 8, width = 11)
Graphe <- CaracBrasSensiBras %>%
  filter(cible %in% c("efftox", "both"), scenar == "ScI2") %>%
  mutate(ttt = gsub("^ttt", "D", ttt),
         methode = factor(methode, levels = c("mBOP", "Simon+TM", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))) %>% 
  ggplot(aes(ttt, rejet_h0, group = n_bras, color = n_bras)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(vars(methode)) +
  labs(x = "Dose", y = "Proportion of conclusion to promising dose", color = "# arms") +
  scale_color_discrete(type = c("black", "darkblue", "darkred")) +
  expand_limits(x = 0) +
  scale_y_continuous(labels = scales::percent_format())
ggsave(Graphe, filename = "Figures/resultbras_sensi_arms_I2.png", device = "png", height = 8, width = 11)

# Sensitivity analysis: priors ----

Graphe <- CaracBrasPriors %>% 
  filter(grepl("^h|^H", methode)) %>% 
  mutate(methode = ifelse(methode == "hBOP_1", "hBOP", methode),
         methode = factor(methode, levels = c("H3_2", "H3_1", "H2_2", "H2_1", "H1_2", "H1_1", "hBOP")),
         ttt = gsub("^ttt", "D", ttt)) %>% 
  ggplot(aes(x = rejet_h0, y = methode, color = ttt, shape = ttt)) +
  geom_point(size = 3) +
  facet_wrap(vars(scenar), scales = "free_x") +
  expand_limits(x = 0) +
  scale_x_continuous(labels = scales:: percent_format()) +
  scale_color_discrete(type = c("orange", "darkred", "darkblue")) + 
  labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose")
ggsave(Graphe, filename = "Figures/resultbras_sensi_priors.png", device = "png", height = 8, width = 11)

Graphe <- CaracBrasPriors %>% 
  filter(grepl("^c|^C", methode)) %>% 
  mutate(methode = ifelse(methode == "cbhmBOP_1", "cbhmBOP", methode),
         methode = factor(methode, levels = c("C2_2", "C2_1", "C1_2", "C1_1", "cbhmBOP")),
         ttt = gsub("^ttt", "D", ttt)) %>% 
  ggplot(aes(x = rejet_h0, y = methode, color = ttt, shape = ttt)) +
  geom_point(size = 3) +
  facet_wrap(vars(scenar), scales = "free_x") +
  expand_limits(x = 0) +
  scale_x_continuous(labels = scales:: percent_format()) +
  scale_color_discrete(type = c("orange", "darkred", "darkblue")) + 
  labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose")
ggsave(Graphe, filename = "Figures/resultbras_sensi_priors2.png", device = "png", height = 8, width = 11)

Graphe <- CaracBrasPriors %>% 
  filter(grepl("^l|^L", methode)) %>% 
  mutate(methode = ifelse(methode == "log1BOP_1", "log1BOP", methode),
         methode = factor(methode, levels = c("L5_2", "L5_1", "L4_2", "L4_1", "L3_2", "L3_1", "L2_2", "L2_1", "L1_2", "L1_1", "log1BOP")),
         ttt = gsub("^ttt", "D", ttt)) %>% 
  ggplot(aes(x = rejet_h0, y = methode, color = ttt, shape = ttt)) +
  geom_point(size = 3) +
  facet_wrap(vars(scenar), scales = "free_x") +
  expand_limits(x = 0) +
  scale_x_continuous(labels = scales:: percent_format()) +
  scale_color_discrete(type = c("orange", "darkred", "darkblue")) + 
  labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose")
ggsave(Graphe, filename = "Figures/resultbras_sensi_priors3.png", device = "png", height = 8, width = 11)

# Sensitivity analysis: CRM models ----

Couleurs <- c("class" = "black", "opti" = "darkgreen", "cons" = "darkred", "ppal" = "darkblue")
TabBrasCrm <- CaracBrasCrm %>% 
  bind_rows(CaracBras %>% filter(methode == "mBOP", scenar %in% c("Sc2", "Sc4", "ScI1", "ScI2"))) %>% 
  mutate(ttt = gsub("^ttt", "D", ttt), 
         label = methode,
         label = gsub("_efftox$", "", label),
         label = ifelse(label == "mBOP", "mBOP_class", label)) %>%
  separate(label, into = c("lablab", "skel")) %>% 
  mutate(skel = factor(skel, levels = c("class", "opti", "cons", "ppal"), labels = c("class", "opti", "pess", "main")),
         lablab = ifelse(skel == "class",
                         paste0("<span style='color:", Couleurs[skel], ";'>", lablab, "</span>"),
                         paste0(lablab, " (<span style='color:", Couleurs[skel], ";'>", skel, "</span>)")))
LabelsVec <- TabBrasCrm %>% 
  count(methode, lablab) %>% 
  mutate(methode = gsub("unfixed", "gunfixed", methode)) %>% 
  arrange(methode) %>% 
  pull(lablab)
Graphe <- TabBrasCrm %>% 
  mutate(lablab = factor(lablab, levels = LabelsVec)) %>% 
  ggplot(aes(x = rejet_h0, y = lablab, color = ttt, shape = ttt)) +
  geom_point(size = 3) +
  facet_wrap(vars(scenar), scales = "free_x") +
  expand_limits(x = 0) +
  theme(axis.text.y = element_markdown()) +
  scale_x_continuous(labels = scales:: percent_format()) +
  scale_color_discrete(type = c("orange", "darkred", "darkblue")) + 
  labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose")
ggsave(Graphe, filename = "Figures/resultbras_sensi_crm_v2.png", device = "png", height = 10, width = 12)


# Sensitivity analysis : Logistic regression with log(d/d*) instead of just d/d* ----

CaracGlobalesLogit <- bind_rows(CaracGlobalesLogit, CaracGlobales %>% filter(methode %in% c("log1BOP", "log2BOP")))
CaracBrasLogit <- bind_rows(CaracBrasLogit, CaracBras %>% filter(methode %in% c("log1BOP", "log2BOP")))
CaracEssaisLogit <- bind_rows(CaracEssaisLogit, CaracEssais %>% filter(methode %in% c("log1BOP", "log2BOP")))

Graphes <- CaracBrasLogit %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("Sc1", "Sc2", "Sc3", "Sc4")),
         methode = factor(methode, levels = rev(c("log1BOP", "verlog1BOP", "log2BOP", "verlog2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_sc1234_veriflog.png", device = "png", 
       height = 27.33 / 2, width = 40.99 / 2, units = "cm", dpi = 300)

Graphes <- CaracBrasLogit %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc5", "Sc6", "Sc7", "Sc8")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("Sc5", "Sc6", "Sc7", "Sc8")),
         methode = factor(methode, levels = rev(c("log1BOP", "verlog1BOP", "log2BOP", "verlog2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_sc5678_veriflog.png", device = "png", 
       height = 27.33 / 2, width = 40.99 / 2, units = "cm", dpi = 300)

Graphes <- CaracBrasLogit %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("ScI1", "ScI2")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         scenar = factor(scenar, levels = c("ScI1", "ScI2")),
         methode = factor(methode, levels = rev(c("log1BOP", "verlog1BOP", "log2BOP", "verlog2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt, shape = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 4) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "orange")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose", shape = "Dose") +
      theme_light(base_size = 15) +
      theme(strip.background = element_rect(fill = "white", color = "black", size = 1.6),
            strip.text = element_text(face = "bold", color = "black"))
  })
wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect")
ggsave(wrap_plots(Graphes, guides = "collect", axes = "collect_y", axis_titles = "collect"), filename = "Figures/rejet_bras_scibru_veriflog.png", device = "png", 
       height = 27.33 / 2, width = 40.99 / 2, units = "cm", dpi = 300)


