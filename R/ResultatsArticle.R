library(tidyverse)
library(rlang)
library(flextable)
library(readxl)
library(patchwork)

theme_set(theme_light(base_size = 14) +
            theme(strip.background = element_rect(fill = "white", color = "black", size = 1.2),
                  strip.text = element_text(face = "bold", color = "black")))
CaracGlobales <- list()
CaracBras <- list()
CaracEssais <- list()

walk(1:6, \(fich_num) {
  Fichiers <- paste0("Data/Simu20250416/resultats_priors_20250416_", 1:6, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobales", append(CaracGlobales, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]]))))), global_env())
  assign("CaracBras", append(CaracBras, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]]))))), global_env())
  assign("CaracEssais", append(CaracEssais, list(do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]]))))), global_env())
})
CaracGlobales <- do.call("rbind", CaracGlobales)
CaracGlobales$methode[CaracGlobales$methode == "mBOP"] <- "mBOP_both"
CaracGlobales$methode <- gsub("bop_log", "boplog", CaracGlobales$methode)
CaracGlobales <- separate(CaracGlobales, methode, c("methode", "cible"), "_")
CaracBras <- do.call("rbind", CaracBras)
CaracBras$methode <- gsub("bop_log", "boplog", CaracBras$methode)
CaracBras$methode[CaracBras$methode == "mBOP"] <- "mBOP_both"
CaracBras <- separate(CaracBras, methode, c("methode", "cible"), "_")
CaracEssais <- do.call("rbind", CaracEssais)
CaracEssais$methode[CaracEssais$methode == "mBOP"] <- "mBOP_both"
CaracEssais$methode <- gsub("bop_log", "boplog", CaracEssais$methode)
CaracEssais <- separate(CaracEssais, methode, c("methode", "cible"), "_")
CaracEssais$larg_ic_eff <- CaracEssais$icsup_eff - CaracEssais$icinf_eff
CaracEssais$larg_ic_tox <- CaracEssais$icsup_tox - CaracEssais$icinf_tox

load("Data/Simu20241205/resultats_priors_crm_20250110.RData")
CaracGlobalesCrm <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]])))
CaracGlobalesCrm <- separate(CaracGlobalesCrm, methode, c("methode", "cible"), "_")
CaracBrasCrm <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]])))
CaracBrasCrm <- separate(CaracBrasCrm, methode, c("methode", "cible"), "_")
CaracEssaisCrm <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]])))
CaracEssaisCrm <- separate(CaracEssaisCrm, methode, c("methode", "cible"), "_")
CaracGlobalesCrm <- bind_rows(CaracGlobalesCrm, CaracGlobales %>% filter(methode %in% c("mBOP", "boplog1")))
CaracBrasCrm <- bind_rows(CaracBrasCrm, CaracBras %>% filter(methode %in% c("mBOP", "boplog1")))
CaracEssaisCrm <- bind_rows(CaracEssaisCrm, CaracEssais %>% filter(methode %in% c("mBOP", "boplog1")))

load("Data/Simu20241205/resultats_priors_crm_20250127.RData")
CaracGlobalesCrmPow <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[1]])))
CaracGlobalesCrmPow <- separate(CaracGlobalesCrmPow, methode, c("methode", "cible"), "_")
CaracBrasCrmPow <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[2]])))
CaracBrasCrmPow <- separate(CaracBrasCrmPow, methode, c("methode", "cible"), "_")
CaracEssaisCrmPow <- do.call("rbind", lapply(ResT, \(x) cbind(x[[3]], x[[4]])))
CaracEssaisCrmPow <- separate(CaracEssaisCrmPow, methode, c("methode", "cible"), "_")
CaracGlobalesCrmPow <- bind_rows(CaracGlobalesCrmPow, CaracGlobales %>% filter(methode %in% c("mBOP", "boplog1")))
CaracBrasCrmPow <- bind_rows(CaracBrasCrmPow, CaracBras %>% filter(methode %in% c("mBOP", "boplog1")))
CaracEssaisCrmPow <- bind_rows(CaracEssaisCrmPow, CaracEssais %>% filter(methode %in% c("mBOP", "boplog1")))

CaracGlobalesPriors <- list()
CaracBrasPriors <- list()
CaracEssaisPriors <- list()
walk(1:11, \(fich_num) {
  Fichiers <- paste0("Data/SimuPrior20250403/resultats_priors_20250331_", 1:11, ".RData")
  fichier_temp <- Fichiers[fich_num]
  load(fichier_temp, envir = current_env())
  assign("CaracGlobalesPriors", append(CaracGlobalesPriors, list(do.call("rbind", lapply(ResT, \(x) if (length(x) == 4) {cbind(x[[3]], x[[1]])})))), global_env())
  assign("CaracBrasPriors", append(CaracBrasPriors, list(do.call("rbind", lapply(ResT, \(x) if (length(x) == 4) {cbind(x[[3]], x[[2]])})))), global_env())
  assign("CaracEssaisPriors", append(CaracEssaisPriors, list(do.call("rbind", lapply(ResT, \(x) if (length(x) == 4) {cbind(x[[3]], x[[4]])})))), global_env())
})
CaracGlobalesPriors <- do.call("rbind", CaracGlobalesPriors)
CaracGlobalesPriors <- separate(CaracGlobalesPriors, methode, c("methode", "centre", "inform", "hypervar_coef", "inform_coef"), "_")
CaracBrasPriors <- do.call("rbind", CaracBrasPriors)
CaracBrasPriors <- separate(CaracBrasPriors, methode, c("methode", "centre", "inform", "hypervar_coef", "inform_coef"), "_")
# CaracEssaisPriors <- do.call("rbind", CaracEssaisPriors)

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc1", "Sc2", "Sc3", "Sc4"), methode %in% c("mBOP", "powBOP", "hBOP", "cbhmBOP", "boplog1", "boplog2")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         methode = gsub("boplog1", "log1BOP", methode),
         methode = gsub("boplog2", "log2BOP", methode),
         methode = factor(methode, levels = rev(c("mBOP", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP")))) %>% 
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] %in% c("Sc1", "Sc2")) {
      LimiteSup <- .05
    } else {
      LimiteSup <- 1
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 2) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_x_continuous(labels = scales::percent_format()) +
      scale_color_discrete(type = c("darkred", "steelblue", "darkgreen")) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose")
  })
wrap_plots(Graphes, guides = "collect")
ggsave(wrap_plots(Graphes, guides = "collect"), filename = "Figures/rejet_bras_sc1234.png", device = "png", height = 7, width = 10)

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc7", "Sc8", "Sc9", "Sc10"), methode %in% c("mBOP", "powBOP", "hBOP", "cbhmBOP", "boplog1", "boplog2")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         methode = gsub("boplog1", "log1BOP", methode),
         methode = gsub("boplog2", "log2BOP", methode),
         methode = factor(methode, levels = rev(c("mBOP", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))), 
         scenar = factor(scenar, levels = c("Sc7", "Sc8", "Sc9", "Sc10"), labels = c("Sc5", "Sc6", "Sc7", "Sc8"))) %>%
  split(.$scenar) %>% 
  map(\(Sc) {
    if (Sc$scenar[1] == "Sc8") {
      LimiteSup <- .9
    } else if (Sc$scenar[1] == "Sc6") {
      LimiteSup <- .9
    } else {
      LimiteSup <- .9
    }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 2) +
      facet_wrap(vars(scenar), scales = "free_x") +
      coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_color_discrete(type = c("darkred", "steelblue", "darkgreen")) +
      scale_x_continuous(breaks = seq(0, .8, .2), labels = scales::percent_format()) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose")
  })
wrap_plots(Graphes, guides = "collect")
ggsave(wrap_plots(Graphes, guides = "collect"), filename = "Figures/rejet_bras_sc5678.png", device = "png", height = 7, width = 10.5)

Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc5", "Sc6"), methode %in% c("mBOP", "powBOP", "hBOP", "cbhmBOP", "boplog1", "boplog2")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         methode = gsub("boplog1", "log1BOP", methode),
         methode = gsub("boplog2", "log2BOP", methode),
         methode = factor(methode, levels = rev(c("mBOP", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))), 
         scenar = factor(scenar, levels = c("Sc5", "Sc6"), labels = c("ScI1", "ScI2"))) %>%
  split(.$scenar) %>% 
  map(\(Sc) {
    # if (Sc$scenar[1] == "Sc8") {
    #   LimiteSup <- .6
    # } else if (Sc$scenar[1] == "Sc6") {
    #   LimiteSup <- .9
    # } else {
    #   LimiteSup <- .8
    # }
    ggplot(Sc, aes(rejet_h0, methode, color = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 2) +
      facet_wrap(vars(scenar), scales = "free_x") +
      expand_limits(x = 0) +
      # coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_color_discrete(type = c("darkred", "steelblue", "darkgreen")) +
      scale_x_continuous(breaks = seq(0, .8, .2), labels = scales::percent_format()) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose")
  })
wrap_plots(Graphes, guides = "collect")
ggsave(wrap_plots(Graphes, guides = "collect"), filename = "Figures/rejet_bras_scibru.png", device = "png", height = 6, width = 10)

# PCS
CaracEssais %>% 
  filter(scenar == "Sc5") %>% 
  group_by(methode, n_simu) %>% 
  summarise(decision = ifelse(all(decision %in% c("Stopping", "Early stopping")), "Aucun", paste(ttt[decision == "Accept the treatment"], collapse = "-"))) %>% 
  ungroup() %>% 
  summarise(pct = sprintf("%.2f%%", 100 * mean(decision == "ttt2-ttt3")), .by = methode)
CaracEssais %>% 
  filter(scenar == "Sc6") %>% 
  group_by(methode, n_simu) %>% 
  summarise(decision = ifelse(all(decision %in% c("Stopping", "Early stopping")), "Aucun", paste(ttt[decision == "Accept the treatment"], collapse = "-"))) %>% 
  ungroup() %>% 
  summarise(pct = sprintf("%.2f%%", 100 * mean(decision == "ttt1-ttt2")), .by = methode)

# Sensitivity analysis
## FWER and power
Donnees <- tribble(
  ~methode, ~centre, ~inform, ~hypervar_coef, ~inform_coef, ~denomination,
  "hBOP", "h0", "peuinf", "var1", "NA", "hBOP",
  "hBOP", "h1", "peuinf", "var1", "NA", "H1",
  "hBOP", "50p", "peuinf", "var1", "NA", "H2",
  "hBOP", "h0", "inf", "var1", "NA", "H3",
  "hBOP", "h0", "noninf", "var1", "NA", "H4",
  "hBOP", "h0", "peuinf", "var5", "NA", "H5",
  "hBOP", "h0", "peuinf", "var0p5", "NA", "H6",
  "cbhmBOP", "h0", "peuinf", "NA", "NA", "cbhmBOP",
  "cbhmBOP", "h1", "peuinf", "NA", "NA", "C1", 
  "cbhmBOP", "50p", "peuinf", "NA", "NA", "C2",
  "cbhmBOP", "h0", "inf", "NA", "NA", "C3",
  "cbhmBOP", "h0", "noninf", "NA", "NA", "C4",
  "logBOP", "h0", "inf", "crois", "inf", "log1BOP",
  "logBOP", "h0", "peuinf", "crois", "inf", "L1",
  "logBOP", "h1", "inf", "crois", "inf", "L2",
  "logBOP", "h1", "peuinf", "crois", "inf", "L3",
  "logBOP", "jef", "jef", "crois", "inf", "L4",
  "logBOP", "jef", "jef", "crois", "peuinf", "L5",
  "logBOP", "jef", "jef", "trescrois", "inf", "L6"
)
map_dfr(seq_len(nrow(Donnees)), \(x) {
  StringRes <- CaracGlobalesPriors %>% 
    filter(methode == Donnees$methode[[x]],
           centre == Donnees$centre[[x]],
           inform == Donnees$inform[[x]],
           hypervar_coef == Donnees$hypervar_coef[[x]],
           inform_coef == Donnees$inform_coef[[x]],
           scenar %in% paste0("Sc", 1:4)) %>% 
    arrange(scenar) %>% 
    mutate(rejet_glob = sprintf("%.2f", 100 * rejet_glob)) %>% 
    pull(rejet_glob) %>% 
    paste(collapse = " & ")
  return(data.frame(
    methode = Donnees$methode[[x]],
    centre = Donnees$centre[[x]],
    inform = Donnees$inform[[x]],
    hypervar_coef = Donnees$hypervar_coef[[x]],
    inform_coef = Donnees$inform_coef[[x]],
    txt = StringRes
  ))
})

ListeGraphes <- map_dfr(seq_len(nrow(Donnees)), \(x) {
  CaracBrasPriors %>% 
    filter(methode == Donnees$methode[[x]],
           centre == Donnees$centre[[x]],
           inform == Donnees$inform[[x]],
           hypervar_coef == Donnees$hypervar_coef[[x]],
           inform_coef == Donnees$inform_coef[[x]],
           scenar %in% paste0("Sc", 5:6)) %>% 
    mutate(scenar = factor(scenar, levels = paste0("Sc", 5:6), labels = paste0("ScI", 1:2)),
           denomination = Donnees$denomination[[x]])
}) %>% 
  split(.$methode) %>% 
  map(\(.Data) {
    Ordre <- unique(.Data$denomination)
    Chiffres <- grepl("[0-9]$", Ordre)
    Ordre <- c(rev(sort(Ordre[Chiffres])), Ordre[!Chiffres])
    ggplot(.Data %>% 
             mutate(ttt = factor(ttt, levels = paste0("ttt", 1:3), labels = paste0("D", 1:3)),
                    denomination = factor(denomination, levels = Ordre)), 
           aes(x = rejet_h0, y = denomination, color = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 2) +
      facet_wrap(vars(scenar), scales = "free_x") +
      expand_limits(x = 0) +
      scale_color_discrete(type = c("darkred", "steelblue", "darkgreen")) +
      scale_x_continuous(breaks = seq(0, .8, .2), labels = scales::percent_format()) +
      labs(y = NULL, x = "Proportion of conclusion to promising dose", color = "Dose")
  })
DesignPlot <- "
AABB
#CC#
"
Graphe <- wrap_plots(ListeGraphes) +
  plot_layout(design = DesignPlot, guides = "collect", widths = c(1, 1, 1.5, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
ggsave(Graphe, filename = "Figures/sensi_prior.png", device = "png", height = 10, width = 12)  


Graphes <- CaracBras %>% 
  filter(cible %in% c("efftox", "both"), scenar %in% c("Sc2", "Sc4", "Sc5", "Sc6"), methode %in% c("cbhmBOP")) %>% 
  mutate(ttt = gsub("ttt", "D", ttt),
         methode = gsub("boplog1", "log1BOP", methode),
         methode = gsub("boplog2", "log2BOP", methode),
         methode = factor(methode, levels = rev(c("mBOP", "powBOP", "hBOP", "cbhmBOP", "log1BOP", "log2BOP"))), 
         scenar = factor(scenar, levels = c("Sc2", "Sc4", "Sc5", "Sc6"), labels = c(paste0("Sc", 1:4)))) %>% 
    ggplot(aes(rejet_h0, scenar, color = ttt)) +
      geom_point(position = position_dodge2(width = .2), size = 2) +
      expand_limits(x = 0) +
      # coord_cartesian(xlim = c(0, LimiteSup)) +
      scale_color_discrete(type = c("darkred", "steelblue", "darkgreen")) +
      scale_x_continuous(breaks = seq(0, .8, .2), labels = scales::percent_format()) +
      labs(y = NULL, x = "Proportion de conclusion\n à une dose prometteuse", color = "Dose")
ggsave(Graphes, filename = "Figures/sensi_cbhm.png", device = "png", height = 6, width = 10)




