# ----------------------------------------------------------- #
# Fonctions utilisées pour analyser les essais et les simuler #
# Auteur : G. Mulier                                          #
# Créé le 20/05/2024, modifié le 05/07/2024                   #
# ----------------------------------------------------------- #

"%nin%" <- function(x, y) match(x, y, nomatch = 0L) == 0L

summarise_ttt <- function(data, groupe, colonne, intitule = {{colonne}}) {
  data %>%
    group_by({{groupe}}) %>%
    summarise("{{intitule}}" := mean({{colonne}}))
}

summarise_decision <- function(data, groupe, colonne, char_decision, intitule) {
  data %>%
    group_by({{groupe}}) %>%
    summarise("{{intitule}}" := mean({{colonne}} == char_decision))
}

summarise_detect <- function(data, groupe, col_decision, char_decision, intitule) {
  data %>%
    group_by({{groupe}}) %>%
    summarise("{{intitule}}" := mean(str_detect({{col_decision}}, char_decision)))
}


gen_patients_multinom <- function(n_sim,
                                  ana_inter,
                                  ana_inter_tox = NULL,
                                  rand_ratio = NULL,
                                  multinom_cont = NULL,
                                  multinom_ttt = list(),
                                  mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                                  seed = 1024,
                                  methode = 2L) {
  # Checks arguments
  if (length(n_sim) != 1 || !is.numeric(n_sim) || n_sim < 0 || n_sim %% 1 != 0)
    stop("\"n_sim\" should be a positive integer.", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (methode %nin% c(1L, 2L, 3L, 4L))
    stop(
      "Choice between 3 methods with an integer:
              * 1L = trial by trial;
              * 2L = trial by trial, patient by patient;
              * 3L = whole in 1;
              * 4L = like method 2, but in C++ to be faster.
      Speed is following: 4L > 3L > 1L > 2L.",
      call. = FALSE
    )
  
  # Checking the specified probabilities
  if (length(multinom_ttt) == 0)
    stop("No specified treatment.", call. = FALSE)
  if (!is.list(multinom_ttt)) {
    multinom_ttt <- list(ttt1 = multinom_ttt)
    warning("\"multinom_ttt\" should be a list. Converted with list(multinom_ttt) assuming that there is only 1 treatment group.", call. = FALSE, immediate. = TRUE)
  }
  if (any(lapply(multinom_ttt, length) != 4))
    stop("Each law in \"multinom_ttt\" should be a vector of length 4.", call. = FALSE)
  if (any(unlist(lapply(multinom_ttt, FUN = function(x) !dplyr::near(sum(x), 1)))))
    stop("For each treatment the law shoudl sum to 1.", call. = FALSE)
  if (is.null(multinom_cont)) {
    message("No control group.")
    if (is.null(names(multinom_ttt))) {
      noms_ttt <- c(paste0("ttt", seq_len(length(multinom_ttt))))
    } else {
      noms_ttt <- c(names(multinom_ttt))
    }
    proba <- multinom_ttt
    names(proba) <- noms_ttt
  } else {
    if (length(multinom_cont) != 4)
      stop("A vector of length 4 must be provided if you want to specify the control group.", call. = FALSE)
    if (!dplyr::near(sum(multinom_cont), 1))
      stop("Control law should sum to 1.", call. = FALSE)
    if (is.null(names(multinom_ttt))) {
      noms_ttt <- c("ttt0", paste0("ttt", seq_len(length(multinom_ttt))))
    } else {
      noms_ttt <- c("ttt0", names(multinom_ttt))
    }
    proba <- append(list(multinom_cont), multinom_ttt)
    names(proba) <- noms_ttt
  }
  # ttt0 is always for control group, and ttt1, 2, ... will be evaluated treatments.
  # If no control group, start at ttt1.
  
  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.")
    anas_inters_cum <- sort(union(cumsum(ana_inter), cumsum(ana_inter_tox)))
    ana_inter     <- c(anas_inters_cum[1], diff(anas_inters_cum))
  }
  
  # Get the sample sizes at each interim analysis
  if (is.null(rand_ratio)) {
    message("Without specifying \"rand_ratio\", all groups are assumed to be of same size.")
    rand_ratio <- rep(1, length(multinom_ttt) + !is.null(multinom_cont))
  } else if (length(rand_ratio) == 1 & length(multinom_ttt) != 1) {
    message("Each arm will have \"rand_ratio\" x \"ana_inter\" patients at each interim analysis.")
    if (!is.numeric(rand_ratio) || rand_ratio < 0)
      stop("\"rand_ratio\" should be a positive real number indicating the ratio between the number of patients in treatment arms and in control arm.", call. = FALSE)
    if (is.null(multinom_cont)) {
      rand_ratio <- rep(rand_ratio, length(multinom_ttt))
    } else {
      rand_ratio <- c(1, rep(rand_ratio, length(multinom_ttt)))
    }
  } else {
    if (!is.numeric(rand_ratio) || any(rand_ratio < 0))
      stop("\"rand_ratio\" should be a vector of positive real numbers indicating the ratio between the number of patients in treatment arms and in control arm.", call. = FALSE)
    if (length(rand_ratio) != length(multinom_ttt))
      stop("\"rand_ratio\" and \"multinom_ttt\"' must have the same length.", call. = FALSE)
    if (!is.null(multinom_cont)) {
      rand_ratio <- c(1, rand_ratio)
    }
  }
  
  anas_inters <- lapply(
    rand_ratio,
    FUN = function(x, y) {x * y},
    y = ana_inter
  ) # Number of additional patients at each interim analysis for each group
  if (lapply(anas_inters, function(x) any(!dplyr::near(x %% 1, 0))) %>% unlist() %>% any()) {
    warning("\"rand_ratio\" x \"ana_inter\" must only return integers.
         The numbers were converted to integers using the ceiling function.",
         call. = FALSE, immediate. = TRUE)
    anas_inters <- lapply(anas_inters,  ceiling)
  }
  names(anas_inters) <- noms_ttt
  nmax <- lapply(anas_inters, sum)
  
  # Number of subjects at each interim analysis for each arm (including the control arm)
  sujet_ana <- lapply(anas_inters, function(x) {
    data.frame(nb_ana = as.integer(seq_along(x)),
               nb_patients = x)
  })
  if (methode != 4L) {
    set.seed(seed)
    on.exit(set.seed(NULL), add = TRUE) # Reset the seed to not impair the environment with the use of the function
  }
  
  if (methode == 3L) {
    
    liste_patients <- list()
    for (i in noms_ttt) {
      resp_list <- lapply(
        seq_along(anas_inters[[i]]),
        FUN = function(n) {
          n_ana <- anas_inters[[i]][n]
          mat <- rmultinom(n_sim, size = n_ana, prob = proba[[i]])
          mat <- apply(
            X = mat,
            MARGIN = 2,
            FUN = function(p) {
              efftox <- p[1]
              effnotox <- p[2]
              noefftox <- p[3]
              noeffnotox <- p[4]
              eff <- sum(mat_beta_xi[1,] * p)
              notox <- sum(mat_beta_xi[2,] * p)
              return(c(efftox, effnotox, noefftox, noeffnotox, eff, notox))
            }
          ) %>%
            t() %>%
            as.data.frame()
          colnames(mat) <- paste0(c("efftox", "effnotox", "noefftox", "noeffnotox", "eff", "notox"), "_ana", n)
          mat$n_simu <- seq_len(n_sim)
          return(mat)
        }
      )
      resp_list <- purrr::reduce(.x = resp_list,
                                 .f = dplyr::left_join,
                                 by = "n_simu")
      
      # Adding total counts at each interim analysis
      resp_list <- resp_list %>%
        tidyr::pivot_longer(
          cols = -n_simu,
          names_pattern = "^(efftox|effnotox|noefftox|noeffnotox|eff|notox)_ana(\\d+)$",
          names_to = c("critere", "nb_ana")
        ) %>%
        tidyr::pivot_wider(names_from = critere, values_from = value) %>%
        dplyr::mutate(nb_ana = as.integer(nb_ana)) %>%
        dplyr::left_join(sujet_ana[[i]], by = "nb_ana") %>%
        dplyr::group_by(n_simu) %>%
        dplyr::mutate(
          efftox = cumsum(efftox),
          effnotox = cumsum(effnotox),
          noefftox = cumsum(noefftox),
          noeffnotox = cumsum(noeffnotox),
          tot_eff = cumsum(eff),
          tot_notox = cumsum(notox),
          tot_pat = cumsum(nb_patients)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-eff,-notox,-nb_patients)
      liste_patients[[i]] <- resp_list %>%
        dplyr::mutate(ttt = i, .before = 3) %>%
        dplyr::relocate(tot_pat, .before = 8)
    }
    
    liste_patients <- do.call("rbind", liste_patients) %>%
      dplyr::arrange(n_simu, ttt, nb_ana)
    
  } else if (methode == 1) {
    
    tab_patients <- lapply(names(sujet_ana), function(x) cbind(sujet_ana[[x]], ttt = x)) %>%
      do.call(what = "rbind", arg = .)
    
    # Sample size for each treatment at each interim analysis
    liste_patients <- expand.grid(
      n_simu = seq_len(n_sim),
      nb_ana = seq_len(length(ana_inter)),
      ttt = noms_ttt
    ) %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::left_join(tab_patients, by = c("nb_ana", "ttt")) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(tot_pat = cumsum(nb_patients)) %>%
      dplyr::ungroup()
    
    # Efficacy and toxicity probabilities
    liste_patients <- liste_patients %>%
      dplyr::mutate(
        p_mult = purrr::map(ttt, ~ proba[[match(.x, names(proba))]]),
        n_mult = purrr::map2(p_mult, nb_patients, function(x, y) {
          mat <- as.double(rmultinom(n = 1, size = y, prob = x))
          names(mat) <- c("efftox", "effnotox", "noefftox", "noeffnotox")
          return(mat)
        })
      ) %>%
      tidyr::unnest_wider(n_mult) %>%
      dplyr::mutate(tot_eff = efftox + effnotox,
                    tot_notox = effnotox + noeffnotox) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(dplyr::across(efftox:tot_notox, ~ cumsum(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-nb_patients, -p_mult) %>%
      dplyr::relocate(tot_pat, .before = 9)
    
  } else if (methode == 2) {
    
    anas_inters_cum <- lapply(anas_inters, cumsum)
    liste_patients <- list()
    
    for (p in noms_ttt) {
      
      liste_patients[[p]] <- do.call("rbind",
                                     lapply(
                                       seq_len(n_sim),
                                       function(x) {
                                         tableau <- rmultinom(nmax[[p]], size = 1, prob = proba[[p]]) %>%
                                           t() %>%
                                           as.data.frame()
                                         tableau <- cbind(n_simu = x, ttt = p, tableau, tot_pat = seq_len(nmax[[p]]))
                                         tableau$nb_ana <- NA
                                         for (i in rev(anas_inters_cum[[p]])) {tableau$nb_ana <- ifelse(tableau$tot_pat <= i, match(i, anas_inters_cum[[p]]), tableau$nb_ana)}
                                         names(tableau)[3:6] <- c("efftox", "effnotox", "noefftox", "noeffnotox")
                                         return(tableau)
                                       }
                                     ))
    }
    
    liste_patients <- do.call("rbind", liste_patients) %>%
      dplyr::group_by(n_simu, ttt, nb_ana) %>%
      dplyr::summarise(dplyr::across(c(efftox:noeffnotox), ~ sum(.x)),
                       tot_pat = dplyr::last(tot_pat),
                       .groups = "drop") %>%
      dplyr::arrange(n_simu, ttt, nb_ana) %>%
      dplyr::group_by(n_simu, ttt) %>%
      dplyr::mutate(dplyr::across(c(efftox:noeffnotox), ~ cumsum(.x))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(tot_eff = efftox + effnotox,
                    tot_notox = effnotox + noeffnotox)
    
  } else if (methode == 4) {
    
    liste_patients <- GenPts(n_sim, anas_inters[[1]], matrix(unlist(proba), ncol = 4, byrow = TRUE), seed)
    liste_patients <- as.data.frame(liste_patients)
    colnames(liste_patients) <- c("ttt", "n_sim", "nb_ana", "tot_pat", "efftox", "effnotox", "noefftox", "noeffnotox", "tot_eff", "tot_notox")
    
  }
  
  return(liste_patients)
  
}



real_essai_bop <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, 
                           phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    n_eff <- data$tot_eff[data$nb_ana == i]
    n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox, 1 - prior_tox + Nb_pts - n_tox)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (n_tox + prior_tox) / (Nb_pts + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox, 1 - prior_tox + Nb_pts - n_tox) 
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox, 1 - prior_tox + Nb_pts - n_tox)
  }
  return(data)
}

real_essai_bop_borrow <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, 
                                  phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  BrasArretes <- rep(FALSE, length(unique(data$ttt)))
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    # Décomptes d'efficacités et de toxicité
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum((n_tox * BrasArretes)[-x]), numeric(1))
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_pts_autres <- vapply(seq_along(n_tox), \(x) sum((n_pts_bras * BrasArretes)[-x]), numeric(1))
    } else { # Si on n'est pas à la première analyse, il ne faut actualiser n_tox et n_pts que pour les cas où on ne s'est pas arrêté
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum((n_tox * BrasArretes)[-x]), numeric(1))
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_pts_autres <- vapply(seq_along(n_pts_bras), \(x) sum((n_pts_bras * BrasArretes)[-x]), numeric(1))
    }
    # Règles d'arrêt sur la proba a posteriori
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_pts_autres - n_tox_autres)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    if (i == 1) {
      BrasArretes <- data$arret_tox[data$nb_ana == i] == 1
    } else {
      BrasArretes[data$arret_eff[data$nb_ana == (i - 1)] == 0] <- data$arret_tox[data$nb_ana == i][data$arret_eff[data$nb_ana == (i - 1)] == 0] == 1
    }
    
  }
  return(data)
}

real_essai_bop_seq <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, 
                                  phi_eff, phi_tox, prior_eff, prior_tox) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      # n_eff_autres <- rep(0, length(n_eff))
      n_tox_autres <- rep(0, length(n_tox))
      # n_ptseff_autres <- rep(0, length(n_pts_bras))
      n_ptstox_autres <- rep(0, length(n_pts_bras))
    } else { # Si on n'est pas à la première analyse, il ne faut actualiser n_tox et n_pts que pour les cas où on ne s'est pas arrêté
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # Partage d'infos des bras arrêtés (à dose inférieure pour tox et supérieure pour eff)
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
      n_ptstox_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
      # n_eff_autres <- vapply(seq_along(n_eff), \(x) sum(n_eff * data$arret_eff[data$nb_ana == (i - 1)] * (x < seq_along(n_eff))), numeric(1))
      # n_ptseff_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras * data$arret_eff[data$nb_ana == (i - 1)] * (x < seq_along(n_eff))), numeric(1))
    }
    # Règles d'arrêt
    # PPEff <- pbeta(phi_eff, prior_eff + n_eff + n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + n_ptseff_autres - n_eff_autres)
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_ptstox_autres - n_tox_autres)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (n_tox + prior_tox + n_tox_autres) / (Nb_pts + n_ptstox_autres + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_ptstox_autres - n_tox_autres)
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_ptstox_autres - n_tox_autres)
  }
  return(data)
}

real_essai_bop_power <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, A0, 
                                 phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_pts_autres <- vapply(seq_along(n_tox), \(x) sum(n_pts_bras[-x]), numeric(1))
    } else { # Si on n'est pas à la première analyse, il ne faut actualiser n_tox et n_pts que pour les cas où on ne s'est pas arrêté
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_pts_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras[-x]), numeric(1))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    # Power prior pour partager l'information sur la toxicité
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0 * (n_pts_autres - n_tox_autres))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  return(data)
}

real_essai_bop_power_test <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, A0, Tox0, SeuilP, 
                                      phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        # Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) fisher.test(matrix(c(n_tox[x], n_pts_bras[x] - n_tox[x], n_tox[t], n_pts_bras[t] - n_tox[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_tox[-x][Pvals > SeuilP]), sum((Pvals > SeuilP) * n_pts_bras[-x])))
      }, numeric(2))
    } else { # Si on n'est pas à la première analyse, il ne faut actualiser n_tox et n_pts que pour les cas où on ne s'est pas arrêté
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # On teste si la toxicité est significativement différente de Tox0 pour partager l'information si ce n'est pas significativement différent
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        # Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) fisher.test(matrix(c(n_tox[x], n_pts_bras[x] - n_tox[x], n_tox[t], n_pts_bras[t] - n_tox[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_tox[-x][Pvals > SeuilP]), sum((Pvals > SeuilP) * n_pts_bras[-x])))
      }, numeric(2))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    # Power prior, mais sur les bras non significativement différents
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (prior_tox + n_tox + A0 * n_toxpts_autres[1, ]) / (Nb_pts + A0 * n_toxpts_autres[2, ] + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
  }
  return(data)
}

real_essai_bop_borrow_test <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, A0, Tox0, 
                                       phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        return(c(sum(n_tox[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
    } else { # Si on n'est pas à la première analyse, il ne faut actualiser n_tox et n_pts que pour les cas où on ne s'est pas arrêté
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        return(c(sum(n_tox[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    # Idem que précédemment, mais la PValue est un multiplicateur, pas juste un seuil pour le partage
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  return(data)
}

ModeleHier <- "
data {
  int<lower = 1> Nb; // Nombre de bras
  int<lower = 1> n[Nb]; // Nombre de patients par bras
  int<lower = 0> y[Nb]; // Nombre de toxicités par bras
}

parameters{
  real p_raw[Nb];
  real mu;
  real<lower = 0> sigma;
}

transformed parameters{
  real p[Nb]; // Probabilité pour chaque bras
  real logit_p[Nb]; // Prédicteur linéaire pour chaque bras
  for(j in 1:Nb){
    logit_p[j] = mu + sigma * p_raw[j]; // Non central parametrization
  }
  p = inv_logit(logit_p);
}

model{

  // prior distributions
  sigma ~ normal(0, 5);
  mu  ~ normal(0, 5);

  for (i in 1:Nb) {
    p_raw[i] ~ normal(0, 1);
    // binomial likelihood
    y[i] ~ binomial(n[i], p[i]);
  }
 
}
"
CompilHier <- stan_model(model_code = ModeleHier)

real_essai_modhier <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, 
                               phi_eff, phi_tox, prior_eff) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
    } else {
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    # Probas a posteriori calculées par MCMC sur le modèle hiérarchique
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      DonneesTox <- list(Nb = length(n_tox), n = n_pts_bras, y = n_tox)
      SampledTox <- sampling(CompilHier, 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 5000,          
                             iter = 30000,
                             thin = 5,
                             cores = 3,
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      PPred <- extract(SampledTox, pars = "p")$p
      PPTox <- colMeans(PPred > phi_tox) # On va chercher les distributions a posteriori de chaque bras pour la proportion de toxicité
      data$est_tox[data$nb_ana == i] <- mean(PPred)
      data$icinf_tox[data$nb_ana == i] <- quantile(PPred, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- quantile(PPred, probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}

ModeleExNex <- "
data {
  int<lower = 1> Nb; // Nombre de bras
  int<lower = 1> n[Nb]; // Nombre de patients par bras
  int<lower = 0> y[Nb]; // Nombre de toxicités par bras
  real mu_nex; // Moyenne du log-odds pour la partie non échangeable
  real<lower = 0> etype_nex; // Sd du log-odds pour la partie non échangeable
  real<lower = 0, upper = 1> prob_ex[Nb];
}

parameters{
  real theta[Nb]; // Log(odds) dans chaque bras
  real mu_ex; // Moyenne du log(odds) pour la distribution EX
  real<lower = 0> etype_ex; // Sd du log(odds) pour la distribution EX
}

transformed parameters{
  real p[Nb]; // Probabilité pour chaque bras
  p = inv_logit(theta); // Theta = prédicteur linéaire pour chaque bras
}

model{

  // prior distributions
  etype_ex ~ normal(0, 1);
  mu_ex  ~ normal(-1.734601, 2.616);

  for (i in 1:Nb) {
    target += log_mix(prob_ex[i],
                      normal_lpdf(theta[i] | mu_ex, etype_ex),
                      normal_lpdf(theta[i] | mu_nex, etype_nex));
    // binomial likelihood
    y[i] ~ binomial(n[i], p[i]);
  }
 
}

generated quantities{
  // Probabilité d'être EX (https://mc-stan.org/docs/stan-users-guide/finite-mixtures.html)
  real postprob_ex[Nb];
  for (i in 1:Nb) {
    postprob_ex[i] = normal_lpdf(theta[i] | mu_ex, etype_ex) + bernoulli_lpmf(0 | prob_ex[i]) - 
      log_sum_exp(normal_lpdf(theta[i] | mu_ex, etype_ex) + bernoulli_lpmf(0 | prob_ex[i]),
                  normal_lpdf(theta[i] | mu_nex, etype_nex) + bernoulli_lpmf(1 | prob_ex[i]));
  }
  postprob_ex = exp(postprob_ex);
}
"
CompilExNex <- stan_model(model_code = ModeleExNex)

ModeleLog1 <- "
data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real beta;
}
model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5); // Prior peu informatif centré sur 0
  // Likelihood
  y ~ bernoulli_logit(alpha + beta * x);
}
"
CompilLog1 <- stan_model(model_code = ModeleLog1)

ModeleLog2 <- "
data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real<lower=0> beta; // Forçage d'une pente positive (hypothèse de monotonicité)
}
model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  // Likelihood
  y ~ bernoulli_logit(alpha + beta * x);
}
"
CompilLog2 <- stan_model(model_code = ModeleLog2)

ModeleLog3 <- "
data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real beta;
}
model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(.5, 5); // Prior peu informatif centré sur plus que 0
  // Likelihood
  y ~ bernoulli_logit(alpha + beta * x);
}
"
CompilLog3 <- stan_model(model_code = ModeleLog3)

ModeleLog4 <- "
data {
  int<lower=0> N;
  vector[N] x1;
  vector[N] x2;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real beta;
  real gamma;
}
model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  gamma ~ normal(0, 5);
  // Likelihood
  y ~ bernoulli_logit(alpha + beta * x1 + gamma * x2);
}
"
CompilLog4 <- stan_model(model_code = ModeleLog4)

ModeleLog5 <- "
data {
  int<lower=0> N;
  vector[N] x1;
  vector[N] x2;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha;
  real<lower = 0> beta;
  real<lower = 0> gamma;
}
model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  gamma ~ normal(0, 5);
  // Likelihood
  y ~ bernoulli_logit(alpha + beta * x1 + gamma * x2);
}
"
CompilLog5 <- stan_model(model_code = ModeleLog5)

real_essai_bayeslog <- function(data, analyses, CPar, PPar, AnaEff, AnaTox, 
                               phi_eff, phi_tox, prior_eff, modele_log) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$dose <- as.numeric(gsub("^ttt(\\d+)$", "\\1", data$ttt)) # Doses 1/2/3/... prises dans ces simulations
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      doses_ttt <- data$dose[data$nb_ana == i]
    } else {
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      doses_ttt[VecNonArrets] <- data$dose[data$nb_ana == i][VecNonArrets]
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaEff) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    # MCMC pour la proba a posteriori seulement si nécessaire (gain de temps)
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      if (modele_log %in% c(1:3)) {
        DonneesTox <- list(N = sum(n_pts_bras),
                           y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                           x = rep(doses_ttt, n_pts_bras))
      } else {
        DonneesTox <- list(N = sum(n_pts_bras),
                           y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                           x1 = rep(c(0, rep(1, length(doses_ttt) - 1)), n_pts_bras),
                           x2 = rep(c(0, 0, rep(1, length(doses_ttt) - 2)), n_pts_bras))
      }
      SampledTox <- sampling(get(paste0("CompilLog", modele_log)), 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 5000,          
                             iter = 55000,
                             thin = 10,
                             cores = 3,
                             # J'ai dû ajouter cela car cela ne convergeait pas bien sinon pour le modèle qui restreint beta en positif
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      DistPost <- extract(SampledTox)
      if (modele_log %in% c(1:3)) {
        PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + x * DistPost$beta))))
        PPTox <- vapply(PPred, \(distrib) {mean(distrib > phi_tox)}, double(1))
      } else {
        PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + (x != 1) * DistPost$beta + (x == 3) * DistPost$gamma))))
        PPTox <- vapply(PPred, \(distrib) {mean(distrib > phi_tox)}, double(1))
      }
      data$est_tox[data$nb_ana == i] <- vapply(PPred, mean, double(1))
      data$icinf_tox[data$nb_ana == i] <- vapply(PPred, quantile, double(1), probs = .025)
      data$icsup_tox[data$nb_ana == i] <- vapply(PPred, quantile, double(1), probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% AnaTox) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


opcharac <- function(ana_inter,
                     ana_inter_tox = NULL,
                     p_n,
                     prior = NULL,
                     p_a,
                     phi = NULL,
                     mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                     CPar = .6,
                     PPar = .8,
                     methode = c("bop", "bop_borrow", "bop_power", "bop_power_test", "bop_borrow_test", "hier_tox",
                                 "bop_log1", "bop_log2", "bop_log3", "bop_log4", "bop_log5", "bop_seq"),
                     A0 = 0, SeuilP = .2, 
                     tableau_essais = NULL) {
  
  
  # Check arguments ----
  methode <- match.arg(methode)
  if (length(CPar) != 1 | length(PPar) != 1) 
    stop("CPar and PPar should be hyperparameters of the threshold determined by function deter_cutoff.", call. = FALSE)
  if (!is.matrix(mat_beta_xi) || any(dim(mat_beta_xi) != c(2, 4)))
    stop("\"mat_beta_xi\" should be a 2 by 4 matrix.", call. = FALSE)
  if (!dplyr::near(sum(p_n), 1))
    stop("Vector \"p_n\" should sum to 1.", call. = FALSE)
  if (length(p_n) != 4)
    stop("\"p_n\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (is.null(prior)) {
    prior <- p_n
  } else {
    if (length(prior) != 4)
      stop("Vector \"prior\" should be of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  }
  if (!dplyr::near(sum(p_a), 1))
    stop("Vector \"p_a\" should sum to 1.", call. = FALSE)
  if (length(p_a) != 4)
    stop("\"p_a\" should be a vector of length 4:\n1 = Pr(Eff & Tox), 2 = Pr(Eff & no Tox), 3 = Pr(no Eff & Tox) and 4 = Pr(no Eff & no Tox).", call. = FALSE)
  if (!is.numeric(ana_inter) | any(ana_inter < 0) | any(ana_inter %% 1 != 0))
    stop("\"ana_inter\" should be a vector of positive integers.", call. = FALSE)
  if (!is.null(ana_inter_tox)) {
    if (any(ana_inter_tox %% 1 != 0))
      stop("\"ana_inter_tox\" represents the number of supplementary patients at each interim analysis for toxicity and should thus be composed of integers.", call. = FALSE)
    if (sum(ana_inter) != sum(ana_inter_tox))
      stop("\"ana_inter\" and \"ana_inter_tox\" should sum to the same amount of patients.", call. = FALSE)
  }
  
  # Precompute variables ----
  
  # Reconstitution of phi. If no specified value, the values under H0 are taken.
  if (is.null(phi)) {
    phitox <- c(sum(p_n * mat_beta_xi[1, ]), 1 - sum(p_n * mat_beta_xi[2, ]))
    phi <- phitox
    phinotox <- c(sum(p_n * mat_beta_xi[1, ]), sum(p_n * mat_beta_xi[2, ]))
  } else {
    phitox <- phi
    phinotox <- c(phi[1], 1 - phi[2])
  }
  
  # Number of arms
  n_bras <- length(unique(tableau_essais$ttt))
  
  # Interim analyses
  ana_eff <- cumsum(ana_inter)
  ana_tox <- if (is.null(ana_inter_tox)) cumsum(ana_inter) else cumsum(ana_inter_tox)
  analyses <- sort(union(ana_eff, ana_tox))
  
  phi_eff <- phitox[1] # Phi for efficacy
  phi_tox <- phitox[2] # Phi for toxicity
  prior_eff <- sum(prior * mat_beta_xi[1,])
  prior_tox <- sum(prior * (1 - mat_beta_xi[2,]))
  
  # Compute the trials ----
  
  # Splitting the data to perform separately each trial
  tableau_essais <- split(tableau_essais, tableau_essais$n_simu)
  tableau_essais <- map_dfr(
    .x = tableau_essais,
    .f = function(data) {
      if (methode == "bop") {
        real_essai_bop(data, analyses, CPar, PPar, ana_eff, ana_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_borrow") {
        real_essai_bop_borrow(data, analyses, CPar, PPar, ana_eff, ana_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power") {
        real_essai_bop_power(data, analyses, CPar, PPar, ana_eff, ana_tox, A0, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power_test") {
        real_essai_bop_power_test(data, analyses, CPar, PPar, ana_eff, ana_tox, A0, phi_tox, SeuilP, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_borrow_test") {
        real_essai_bop_borrow_test(data, analyses, CPar, PPar, ana_eff, ana_tox, A0, phi_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "hier_tox") {
        real_essai_modhier(data, analyses, CPar, PPar, ana_eff, ana_tox, phi_eff, phi_tox, prior_eff)
      } else if (methode %in% c("bop_log1", "bop_log2", "bop_log3", "bop_log4", "bop_log5")) {
        ModeleAPrendre <- gsub("^bop_log(\\d)$", "\\1", methode)
        real_essai_bayeslog(data, analyses, CPar, PPar, ana_eff, ana_tox, phi_eff, phi_tox, prior_eff, ModeleAPrendre)
      } else if (methode == "bop_seq") {
        real_essai_bop_seq(data, analyses, CPar, PPar, ana_eff, ana_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      }
      
    }
  )
  
  # Decision rules
  tableau_essais <- tableau_essais %>%
    mutate(
      prior_eff = prior_eff,
      prior_tox = prior_tox,
      decision = dplyr::case_when((nb_ana != max(nb_ana)) & (arret_eff == 0) & (arret_tox == 0) ~ "Continue",
                                  (nb_ana == max(nb_ana)) & (arret_eff == 0) & (arret_tox == 0) ~ "Accept the treatment",
                                  (nb_ana != max(nb_ana)) & ((arret_eff != 0) | (arret_tox != 0)) ~ "Early stopping",
                                  (nb_ana == max(nb_ana)) & ((arret_eff != 0) | (arret_tox != 0)) ~ "Stopping"),
      decision_eff = dplyr::case_when((nb_ana != max(nb_ana)) & (arret_eff == 0) ~ "Continue",
                                      (nb_ana == max(nb_ana)) & (arret_eff == 0) ~ "Accept the treatment",
                                      (nb_ana != max(nb_ana)) & (arret_eff != 0) ~ "Early stopping (futility)",
                                      (nb_ana == max(nb_ana)) & (arret_eff != 0) ~ "Stopping (futility)"),
      decision_tox = dplyr::case_when((nb_ana != max(nb_ana)) & (arret_tox == 0) ~ "Continue",
                                      (nb_ana == max(nb_ana)) & (arret_tox == 0) ~ "Accept the treatment",
                                      (nb_ana != max(nb_ana)) & (arret_tox != 0) ~ "Early stopping (toxicity)",
                                      (nb_ana == max(nb_ana)) & (arret_tox != 0) ~ "Stopping (toxicity)")) %>%
    filter(decision != "Continue") %>%
    group_by(n_simu, ttt) %>%
    arrange(n_simu, ttt, nb_ana) %>%
    slice(1) %>%
    ungroup()
  
  
  # Output of the function ----
  return(list(
    carac_globales = data.frame(
      p_n = paste0("(", paste(p_n, collapse = "/"), ")"),
      p_a = paste0("(", paste(p_a, collapse = "/"), ")"),
      rejet_glob = tableau_essais %>%
        group_by(n_simu) %>%
        summarise(
          rejet_glob = sum(decision == "Accept the treatment") > 0,
          .groups = "drop"
        ) %>%
        summarise(rejet_glob = mean(rejet_glob)) %>%
        pull(rejet_glob)
    ),
    carac_bras = data.frame(
      summarise_decision(tableau_essais, ttt, decision, "Accept the treatment", rejet_h0) %>%
        left_join(summarise_decision(tableau_essais, ttt, decision, "Early stopping", arret_precoce), by = "ttt") %>%
        left_join(summarise_detect(tableau_essais, ttt, decision, "Stopping|stopping", arret), by = "ttt") %>%
        left_join(summarise_decision(tableau_essais, ttt, decision_eff, "Early stopping (futility)", arret_precoce_fut), by = "ttt") %>%
        left_join(summarise_decision(tableau_essais, ttt, decision_tox, "Early stopping (toxicity)", arret_precoce_tox), by = "ttt") %>%
        left_join(summarise_detect(tableau_essais, ttt, decision_eff, "futility", arret_fut), by = "ttt") %>%
        left_join(summarise_detect(tableau_essais, ttt, decision_tox, "toxicity", arret_tox), by = "ttt") %>%
        left_join(summarise_ttt(tableau_essais, ttt, tot_pat), by = "ttt") %>%
        left_join(summarise_ttt(tableau_essais, ttt, tot_eff), by = "ttt") %>%
        left_join(summarise_ttt(tableau_essais, ttt, tot_pat - tot_notox, tot_tox), by = "ttt")
    ),
    essais = tableau_essais[, c("n_simu", "ttt", "nb_ana", "est_eff", "icinf_eff", "icsup_eff", "est_tox", "icinf_tox", "icsup_tox")]
  ))
  
}
