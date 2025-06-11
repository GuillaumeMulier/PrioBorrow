# --------------------------- #
# Create illustrative figures #
# Author: G. Mulier           #
# Created 2024/07/05          #
# Modified 2024/12/04         #
# --------------------------- #

# Packages and helpers ----

library(tidyverse)
library(patchwork)
library(ggtext)

theme_set(theme_light(base_size = 22) +
            theme(axis.text.x = element_markdown()))


# Scenarios ----

# ListeScenars <- list(Sc1 = list(Bras1 = c(0.15, 0.15, 0.25, 0.45), Bras2 = c(0.15, 0.15, 0.25, 0.45), Bras3 = c(0.15, 0.15, 0.25, 0.45)),
#                      Sc2 = list(Bras1 = c(0.20, 0.30, 0.10, 0.40), Bras2 = c(0.20, 0.30, 0.10, 0.40), Bras3 = c(0.20, 0.30, 0.10, 0.40)),
#                      Sc3 = list(Bras1 = c(0.10, 0.20, 0.15, 0.55), Bras2 = c(0.20, 0.30, 0.10, 0.40), Bras3 = c(0.25, 0.35, 0.10, 0.30)),
#                      Sc4 = list(Bras1 = c(0.15, 0.25, 0.10, 0.50), Bras2 = c(0.22, 0.38, 0.08, 0.32), Bras3 = c(0.25, 0.35, 0.10, 0.30)),
#                      Sc5 = list(Bras1 = c(0.18, 0.22, 0.17, 0.43), Bras2 = c(0.25, 0.25, 0.15, 0.35), Bras3 = c(0.27, 0.23, 0.18, 0.32)),
#                      Sc6 = list(Bras1 = c(0.20, 0.30, 0.15, 0.35), Bras2 = c(0.25, 0.30, 0.15, 0.30), Bras3 = c(0.30, 0.30, 0.15, 0.25)),
#                      Sc7 = list(Bras1 = c(0.20, 0.30, 0.12, 0.38), Bras2 = c(0.22, 0.28, 0.13, 0.37), Bras3 = c(0.24, 0.26, 0.12, 0.38)))
# ListeScenars <- list(
#   Sc1 = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
#   Sc2 = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
#   Sc3 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.19, 0.23, 0.16, 0.42), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
#   Sc4 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.22, 0.18, 0.40), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
#   Sc5 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.16, 0.26, 0.15, 0.43), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
#   Sc6 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.24, 0.15, 0.43), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
#   Sc7 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.22, 0.18, 0.40), ttt3 = c(0.22, 0.28, 0.16, 0.34)),
#   Sc8 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.15, 0.40), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
#   Sc9 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.18, 0.37), ttt3 = c(0.25, 0.25, 0.15, 0.35)),
#   Sc10 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.27, 0.13, 0.42), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
#   Sc11 = list(ttt1 = c(0.13, 0.22, 0.15, 0.50), ttt2 = c(0.18, 0.27, 0.15, 0.40), ttt3 = c(0.20, 0.30, 0.15, 0.35)),
#   Sc12 = list(ttt1 = c(0.15, 0.20, 0.15, 0.50), ttt2 = c(0.20, 0.25, 0.18, 0.37), ttt3 = c(0.22, 0.28, 0.16, 0.34))
# )
# ListeScenars <- list(
#   Sc1 = list(ttt1 = c(0.15, 0.15, 0.25, 0.45), ttt2 = c(0.15, 0.15, 0.25, 0.45), ttt3 = c(0.15, 0.15, 0.25, 0.45)),
#   Sc2 = list(ttt1 = c(0.20, 0.30, 0.10, 0.40), ttt2 = c(0.20, 0.30, 0.10, 0.40), ttt3 = c(0.20, 0.30, 0.10, 0.40)),
#   Sc3 = list(ttt1 = c(0.12, 0.18, 0.18, 0.52), ttt2 = c(0.17, 0.23, 0.18, 0.42), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
#   Sc4 = list(ttt1 = c(0.18, 0.32, 0.12, 0.38), ttt2 = c(0.22, 0.28, 0.15, 0.35), ttt3 = c(0.23, 0.27, 0.17, 0.33)),
#   Sc5 = list(ttt1 = c(0.15, 0.25, 0.15, 0.45), ttt2 = c(0.20, 0.30, 0.13, 0.37), ttt3 = c(0.25, 0.30, 0.15, 0.30)),
#   Sc6 = list(ttt1 = c(0.16, 0.34, 0.11, 0.39), ttt2 = c(0.18, 0.32, 0.12, 0.38), ttt3 = c(0.20, 0.30, 0.13, 0.37)),
#   Sc7 = list(ttt1 = c(0.12, 0.18, 0.18, 0.52), ttt2 = c(0.13, 0.17, 0.22, 0.48), ttt3 = c(0.19, 0.21, 0.21, 0.39)),
#   Sc8 = list(ttt1 = c(0.17, 0.23, 0.18, 0.42), ttt2 = c(0.23, 0.27, 0.17, 0.33), ttt3 = c(0.23, 0.27, 0.17, 0.33))
# )
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
# Couleurs <- c("prometteur" = "green", "futile/toxique" = "red", "futile" = "coral", "toxique" = "darkred", "intermediaire" = "steelblue")
Graphe <- do.call("rbind", lapply(names(ListeScenars), \(l) {
  do.call("rbind", lapply(names(ListeScenars[[l]]), \(x) {
    data.frame(scenar = l, 
               bras = x,
               eff = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][2],
               tox = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][3])
  }))
})) %>% 
  mutate(scenar = factor(scenar, 
                         levels = c(paste0("Sc", 1:8), "ScI1", "ScI2"),
                         labels = NomsScenars[c(paste0("Sc", 1:8), "ScI1", "ScI2")])) %>%
  mutate(efficacite = case_when(eff < .31 ~ "futile",
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
  scale_color_discrete(type = c("darkblue", "darkred")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "Treatment dosis", y = "Probability", color = "Endpoint", shape = "Decision")
       # caption = "Vert=prometteur / rouge=futile et toxique / beige=futile / rouge foncé=toxique / bleu=intermédiaire")
g <- ggplot_gtable(ggplot_build(Graphe))
stripes <- which(grepl("strip-t", g$layout$name))
for (i in stripes[1:2]) {
  g$grobs[[i]]$grobs[[1]]$children[[1]]$gp$fill <- "steelblue"
}
grid::grid.draw(g)
ggsave(g, filename = "Figures/scenar_simul.png", height = 8, width = 10)
