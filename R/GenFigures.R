# --------------------------- #
# Create illustrative figures #
# Author: G. Mulier           #
# Created 2024/07/05          #
# Modified 2024/07/05         #
# --------------------------- #

# Packages and helpers ----

library(tidyverse)
library(patchwork)

theme_set(theme_light(base_size = 16))


# Scenarios ----

ListeScenars <- list(Sc1 = list(Bras1 = c(0.15, 0.15, 0.25, 0.45), Bras2 = c(0.15, 0.15, 0.25, 0.45), Bras3 = c(0.15, 0.15, 0.25, 0.45)),
                     Sc2 = list(Bras1 = c(0.20, 0.30, 0.10, 0.40), Bras2 = c(0.20, 0.30, 0.10, 0.40), Bras3 = c(0.20, 0.30, 0.10, 0.40)),
                     Sc3 = list(Bras1 = c(0.10, 0.20, 0.15, 0.55), Bras2 = c(0.20, 0.30, 0.10, 0.40), Bras3 = c(0.25, 0.35, 0.10, 0.30)),
                     Sc4 = list(Bras1 = c(0.15, 0.25, 0.10, 0.50), Bras2 = c(0.22, 0.38, 0.08, 0.32), Bras3 = c(0.25, 0.35, 0.10, 0.30)),
                     Sc5 = list(Bras1 = c(0.18, 0.22, 0.17, 0.43), Bras2 = c(0.25, 0.25, 0.15, 0.35), Bras3 = c(0.27, 0.23, 0.18, 0.32)),
                     Sc6 = list(Bras1 = c(0.20, 0.30, 0.15, 0.35), Bras2 = c(0.25, 0.30, 0.15, 0.30), Bras3 = c(0.30, 0.30, 0.15, 0.25)),
                     Sc7 = list(Bras1 = c(0.20, 0.30, 0.12, 0.38), Bras2 = c(0.22, 0.28, 0.13, 0.37), Bras3 = c(0.24, 0.26, 0.12, 0.38)))
Graphe <- do.call("rbind", lapply(names(ListeScenars), \(l) {
  do.call("rbind", lapply(names(ListeScenars[[l]]), \(x) {
    data.frame(scenar = l, 
               bras = x,
               eff = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][2],
               tox = ListeScenars[[l]][[x]][1] + ListeScenars[[l]][[x]][3])
  }))
})) %>% 
  ggplot(aes(bras)) +
  geom_line(aes(y = eff, color = "Efficacité", group = scenar)) +
  geom_point(aes(y = eff, color = "Efficacité")) +
  geom_line(aes(y = tox, color = "Toxicité", group = scenar)) +
  geom_point(aes(y = tox, color = "Toxicité")) +
  facet_wrap(vars(scenar)) +
  expand_limits(y = 0) +
  labs(x = "Bras de traitement", y = "Probabilité", color = "Critère")
ggsave(Graphe, filename = "Figures/scenar_simul_v2.png", height = 8, width = 10)
