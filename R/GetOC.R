# -------------------------------------------------------------------------------------------------- #
# Functions used to generate data, analyzed simulated trials and summarize operating characteristics #
# Autor : G. Mulier                                                                                  #
# Creation the 20/05/2024, modified the 05/11/2025                                                   #
# -------------------------------------------------------------------------------------------------- #


# Helper functions ----

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


# Generate trials ----

# Function used to generate trials with counts of efficacy and toxicity at each interim analysis supplied
# Forked from package https://github.com/GuillaumeMulier/multibrasBOP2
# Arguments :
## n_sim = number of simulated trials
## ana_inter = supplementary number of patients at each analysis for efficacy relative to precedent analysis for efficacy
## ana_inter_tox = same for toxicity (if NULL, assumed at the same time as efficacy)
## rand_ratio = if control group, randomisation ratio (n_grp / n_cont)
## multinom_cont = vector of length 4 (EffTox, EffNoTox, NoEffTox, NoEffNoTox) for multinomial distribution in control group
## multinom_ttt = list of vectors as described in multinom_cont (1 for each treatment arm)
## mat_beta_xi = matrix of 2 rows and 4 columns for contrasts for BOP2 design with first row for efficacy and second for non toxicity
## seed = seed to replicate results
## methode = chosen method for simulation of results (affects a little computation time)
# Output : data.frame of simulated trials
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
    
  } 
  
  return(liste_patients)
  
}


# Analyse trials data ----

## Simon 2 stage design + Ivanova monitoring of toxicity ----

# Compute the operating characteristics of a Simon design with fixed r1, n1, r and n
# Loop over possible values to determine optimal and minimax design
# Arguments :
## r1, r = number of responses at interim and final analysis
## n1, n = number of patients at interim and final analysis
## pu = probability of efficacy that is unacceptable
## p& = probability of efficacu that is promising
# Output : Vector with operating characteristics of such Simon's design
alpha_puiss_simon <- function(r1, n1, r, n, pu, pa) {
  
  alpha1 <- pbinom(r1, n1, pu)
  alpha2 <- sum(map_dbl(seq(r1 + 1, min(n1, r)), ~ dbinom(.x, n1, pu) * pbinom(r - .x, n - n1, pu)))
  alpha <- 1 - (alpha1 + alpha2)
  
  puissance1 <- pbinom(r1, n1, pa)
  puissance2 <- sum(map_dbl(seq(r1 + 1, min(n1, r)), ~ dbinom(.x, n1, pa) * pbinom(r - .x, n - n1, pa)))
  puissance <- 1 - (puissance1 + puissance2)
  
  EN_p0 <- n1 + (1 - alpha1) * (n - n1)
  
  return(c("r1" = r1, "n1" = n1, "r" = r, "n" = n, "pu" = pu, "pa" = pa, "alpha" = alpha, "puissance" = puissance, "EN_p0" = EN_p0, "PET_po" = alpha1))
  
}

# Stopping rules for toxicity with Ivanova's design for monitoring toxicity
# Arguments :
## prior = vector of length 2 giving the 2 parameters of prior beta distribution of toxicity
## ana_inter = number of patients at each analysis
## seuil, critere = values for the stopping rule : stop for toxicity if Pr(p_tox > seuil | Dn) > critere
# Output : data.frame with maximum number of toxicities at each analysis to continue the trial
regle_arret <- function(prior, ana_inter, seuil, critere) {
  
  map_dfr(
    .x = ana_inter,
    .f = function(x) {
      vec_x <- seq_len(x + 1) - 1
      proba <- map_dbl(vec_x, ~ 1 - pbeta(seuil, prior[1] + .x, prior[2] + x - .x))
      proba <- proba <= critere
      return(c(n = x, tox_max = sum(proba) - 1))
    }
  )
  
}

# Results of all simulated trials stored in liste_essai (generated with package multibrasBOP2 or with the forked function above)
# Arguments :
## p_n, p_a = vector of length 4 for unpromising and promising treatment (EffTox, EffNoTox, NoEffTox, NoEffNoTox)
## CaracSeuilSimon = vector of decision rules for Simon's design (r1, r, n1, n)
## CaracSeuilIva = vector of decision rules for Ivanova's design with 1 interim analysis (t1, t, n1, n)
# Output : list of
## global operating characteristics for the design
## arm-wise operating characteristics for the design
## each trials and its analysis simulated
simu_simon <- function(p_n, p_a,
                       tableau_essais, 
                       CaracSeuilSimon,
                       CaracSeuilIva) {
  
  # Stopping rules for Simon's design
  tab_eff <- data.frame(
    nb_ana = 1:2,
    tot_pat = c(CaracSeuilSimon[["n1"]], CaracSeuilSimon[["n"]]),
    seuil_eff = c(CaracSeuilSimon[["r1"]], CaracSeuilSimon[["r"]])
  )
  
  # Stopping rules for Ivanova's design
  tab_tox <- data.frame(
    nb_ana = 1:2,
    tot_pat = c(CaracSeuilIva[["n1"]], CaracSeuilIva[["n"]]),
    seuil_tox = c(CaracSeuilIva[["t1"]], CaracSeuilIva[["t"]])
  )
  
  # Apply decision rules to the simulated trials
  ## 2 analyses because of Simon's design
  tableau_essais <- left_join(tableau_essais, tab_eff, by = c("nb_ana", "tot_pat")) %>% 
    left_join(tab_tox, by = c("nb_ana", "tot_pat")) %>% 
    mutate(tot_tox = tot_pat - tot_notox,
           decision = case_when((nb_ana != 2) & tot_eff > seuil_eff & tot_tox <= seuil_tox ~ "Continue",
                                       (nb_ana == 2) & tot_eff > seuil_eff & tot_tox <= seuil_tox ~ "Accept the treatment",
                                       (nb_ana != 2) & (tot_eff <= seuil_eff | tot_tox > seuil_tox) ~ "Early stopping",
                                       (nb_ana == 2) & (tot_eff <= seuil_eff | tot_tox > seuil_tox) ~ "Stopping"),
           decision_eff = case_when((nb_ana != 2) & tot_eff > seuil_eff ~ "Continue",
                                           (nb_ana == 2) & tot_eff > seuil_eff ~ "Accept the treatment",
                                           (nb_ana != 2) & tot_eff <= seuil_eff ~ "Early stopping (futility)",
                                           (nb_ana == 2) & tot_eff <= seuil_eff ~ "Stopping (futility)"),
           decision_tox = case_when((nb_ana != 2) & tot_tox <= seuil_tox ~ "Continue",
                                           (nb_ana == 2) & tot_tox <= seuil_tox ~ "Accept the treatment",
                                           (nb_ana != 2) & tot_tox > seuil_tox ~ "Early stopping (toxicity)",
                                           (nb_ana == 2) & tot_tox > seuil_tox ~ "Stopping (toxicity)")) %>%
    filter(decision != "Continue") %>%
    arrange(n_simu, ttt, nb_ana) %>%
    slice(1, .by = c(n_simu, ttt)) %>%
    ungroup()
  tableau_essais <- bind_cols(tableau_essais, map2_dfr(tableau_essais$tot_eff, tableau_essais$tot_pat, \(x, n) {
    Test <- binom.test(x, n) # Clopper-pearson confidence interval for efficacy
    return(c("est_eff" = as.numeric(Test$estimate), "icinf_eff" = Test$conf.int[1], "icsup_eff" = Test$conf.int[2]))
  }))
  tableau_essais <- bind_cols(tableau_essais, map2_dfr(tableau_essais$tot_tox, tableau_essais$tot_pat, \(x, n) {
    ParAlpha <- CaracSeuilIva[["a"]] + x
    ParBeta <- CaracSeuilIva[["b"]] + n - x
    Estimation <- ParAlpha / (ParAlpha + ParBeta)
    IcInf <- qbeta(.025, ParAlpha, ParBeta)
    IcSup <- qbeta(.975, ParAlpha, ParBeta)
    return(c("est_tox" = Estimation, "icinf_tox" = IcInf, "icsup_tox" = IcSup)) # Beta credibility interval for Ivanova's design
  }))
  
  return(list(
    carac_globales = data.frame(
      p_n = paste0("(", paste(p_n, collapse = "/"), ")"),
      p_a = paste0("(", paste(p_a, collapse = "/"), ")"),
      nb_essais = length(unique(tableau_essais$n_simu)),
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
        left_join(summarise_ttt(tableau_essais, ttt, tot_tox), by = "ttt")
    ),
    essais = tableau_essais[, c("n_simu", "ttt", "nb_ana", "est_eff", "icinf_eff", "icsup_eff", "est_tox", "icinf_tox", "icsup_tox", "decision", "decision_eff", "decision_tox")]
  ))
  
}


## Multi-arm BOP ----

# Analyse the data from 1 trial using multi-arm BOP2
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox, 1 - prior_tox + Nb_pts - n_tox)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
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


## BOP borrow ----

# Analyse the data from 1 trial using multi-arm BOP2 by borrowing 100% information from stopped arms
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_borrow <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                  phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  BrasArretes <- rep(FALSE, length(unique(data$ttt)))
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    # Counts of efficacy and toxicity
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum((n_tox * BrasArretes)[-x]), numeric(1))
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_pts_autres <- vapply(seq_along(n_tox), \(x) sum((n_pts_bras * BrasArretes)[-x]), numeric(1))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum((n_tox * BrasArretes)[-x]), numeric(1))
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_pts_autres <- vapply(seq_along(n_pts_bras), \(x) sum((n_pts_bras * BrasArretes)[-x]), numeric(1))
    }
    # Decision rules on posterio
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_pts_autres - n_tox_autres)
    data$est_tox[data$nb_ana == i] <- (n_tox + prior_tox + n_tox_autres) / (Nb_pts + n_pts_autres + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_pts_autres - n_tox_autres) 
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_pts_autres - n_tox_autres)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
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


## Sequential borrow BOP ----

### Toxicity only ----

# Analyse the data from 1 trial using multi-arm BOP2 with sequential borrowing (if an arm is stopped for toxicity, borrow toxicity information for higher doses arms)
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_seq_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
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
      n_tox_autres <- rep(0, length(n_tox))
      n_ptstox_autres <- rep(0, length(n_pts_bras))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # borrowing for stopped arm for toxicity at inferior dose
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
      n_ptstox_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
    }
    # Stopping rules via posterior
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_ptstox_autres - n_tox_autres)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
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

### Efficacy and toxicity ----

# Analyse the data from 1 trial using multi-arm BOP2 with sequential borrowing 
# (if an arm is stopped for toxicity, borrow toxicity information for higher doses arms)
# (if an arm is stopped for futility, borrow efficacy information for lower doses arms)
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_seq_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
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
      n_eff_autres <- rep(0, length(n_eff))
      n_tox_autres <- rep(0, length(n_tox))
      n_ptseff_autres <- rep(0, length(n_pts_bras))
      n_ptstox_autres <- rep(0, length(n_pts_bras))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # Information borrowing for both efficacy and toxicity
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
      n_ptstox_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras * data$arret_tox[data$nb_ana == (i - 1)] * (x > seq_along(n_tox))), numeric(1))
      n_eff_autres <- vapply(seq_along(n_eff), \(x) sum(n_eff * data$arret_eff[data$nb_ana == (i - 1)] * (x < seq_along(n_eff))), numeric(1))
      n_ptseff_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras * data$arret_eff[data$nb_ana == (i - 1)] * (x < seq_along(n_eff))), numeric(1))
    }
    # Stopping rules via posterior
    PPEff <- pbeta(phi_eff, prior_eff + n_eff + n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + n_ptseff_autres - n_eff_autres)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff + n_eff_autres) / (Nb_pts + n_ptseff_autres + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff + n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + n_ptseff_autres - n_eff_autres) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff + n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + n_ptseff_autres - n_eff_autres)
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + n_ptstox_autres - n_tox_autres)
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
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


## Power prior BOP ----

### Toxicity only ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on toxicity
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0 = exponent of the power prior to determine borrowing
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_power_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0, 
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
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_pts_autres <- vapply(seq_along(n_tox), \(x) sum(n_pts_bras[-x]), numeric(1))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_pts_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras[-x]), numeric(1))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    # Power prior on toxicity
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0 * (n_pts_autres - n_tox_autres))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (n_tox + A0 * n_tox_autres + prior_tox) / (Nb_pts + A0 * n_pts_autres + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0 * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0 * (n_pts_autres - n_tox_autres)) 
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0 * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0 * (n_pts_autres - n_tox_autres))
  }
  return(data)
}

### Efficacy and toxicity ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on efficacy and toxicity
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0_eff, A0_tox = exponent of the power prior to determine borrowing (by default, both equals)
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_power_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_eff, A0_tox = A0_eff, 
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
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_eff_autres <- vapply(seq_along(n_eff), \(x) sum(n_eff[-x]), numeric(1))
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_pts_autres <- vapply(seq_along(n_tox), \(x) sum(n_pts_bras[-x]), numeric(1))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_tox_autres <- vapply(seq_along(n_tox), \(x) sum(n_tox[-x]), numeric(1))
      n_eff_autres <- vapply(seq_along(n_eff), \(x) sum(n_eff[-x]), numeric(1))
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_pts_autres <- vapply(seq_along(n_pts_bras), \(x) sum(n_pts_bras[-x]), numeric(1))
    }
    # Power prior
    PPEff <- pbeta(phi_eff, prior_eff + n_eff + A0_eff * n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_pts_autres - n_eff_autres))
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + A0_eff * n_eff_autres + prior_eff) / (Nb_pts + A0_eff * n_pts_autres + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff + A0_eff * n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_pts_autres - n_eff_autres)) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff + A0_eff * n_eff_autres, 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_pts_autres - n_eff_autres))
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0_tox * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_pts_autres - n_tox_autres))
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (n_tox + A0_tox * n_tox_autres + prior_tox) / (Nb_pts + A0_tox * n_pts_autres + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0_tox * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_pts_autres - n_tox_autres)) 
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0_tox * n_tox_autres, 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_pts_autres - n_tox_autres))
  }
  return(data)
}


## Power prior BOP with binomial/fisher test ----

### Toxicity only ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on toxicity (if no evidence of difference by binomial test)
# Adapted from test-then-pool strategy from Viele et al. Use of historical control data for assessing treatment effects in clinical trials. Pharm Stat. 2014 Jan-Feb;13(1):41-54. 
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0 = exponent of the power prior to determine borrowing (by default, both equals)
## SeuilP = threshold on Pvalue to assess significant difference
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_power_test_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0, SeuilP, 
                                          # Tox0,
                                          phi_eff, phi_tox, prior_eff, prior_tox) {
  # Tox0 is commented because we want to compare arm between them, not compare arm to a reference value as in historical borrowing
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
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # Test if toxicity is different between arms and borrow if no difference found
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    # Power prior, but when no significant difference at threshold SeuilP
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
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

### Efficacy and toxicity ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on toxicity and efficacy (if no evidence of difference by binomial test)
# Adapted from test-then-pool strategy from Viele et al. Use of historical control data for assessing treatment effects in clinical trials. Pharm Stat. 2014 Jan-Feb;13(1):41-54. 
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0_eff, A0_tox = exponent of the power prior to determine borrowing (by default, both equals)
## SeuilP_eff, SeuilP_tox = threshold on Pvalue to assess significant difference (by default, both equal)
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_power_test_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                             A0_eff, SeuilP_eff, A0_tox = A0_eff, SeuilP_tox = SeuilP_eff, 
                                             # Tox0, Eff0,
                                             phi_eff, phi_tox, prior_eff, prior_tox) {
  # Tox0/Eff0 is commented because we want to compare arm between them, not compare arm to a reference value as in historical borrowing
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
      n_effpts_autres <- vapply(seq_along(n_eff), \(x) {
        # Pvals <- vapply(seq_along(n_eff)[-x], \(t) binom.test(n_eff[t], n_pts_bras[t], Eff0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_eff)[-x], \(t) fisher.test(matrix(c(n_eff[x], n_pts_bras[x] - n_eff[x], n_eff[t], n_pts_bras[t] - n_eff[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_eff[-x][Pvals > SeuilP_eff]), sum((Pvals > SeuilP_eff) * n_pts_bras[-x])))
      }, numeric(2))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        # Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) fisher.test(matrix(c(n_tox[x], n_pts_bras[x] - n_tox[x], n_tox[t], n_pts_bras[t] - n_tox[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_tox[-x][Pvals > SeuilP_tox]), sum((Pvals > SeuilP_tox) * n_pts_bras[-x])))
      }, numeric(2))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      # Test if efficacy and toxicity are different between arms and borrow when no difference found
      n_effpts_autres <- vapply(seq_along(n_eff), \(x) {
        # Pvals <- vapply(seq_along(n_eff)[-x], \(t) binom.test(n_eff[t], n_pts_bras[t], Eff0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_eff)[-x], \(t) fisher.test(matrix(c(n_eff[x], n_pts_bras[x] - n_eff[x], n_eff[t], n_pts_bras[t] - n_eff[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_eff[-x][Pvals > SeuilP_eff]), sum((Pvals > SeuilP_eff) * n_pts_bras[-x])))
      }, numeric(2))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        # Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) fisher.test(matrix(c(n_tox[x], n_pts_bras[x] - n_tox[x], n_tox[t], n_pts_bras[t] - n_tox[t]), 2, 2))$p.value, numeric(1))
        return(c(sum(n_tox[-x][Pvals > SeuilP_tox]), sum((Pvals > SeuilP_tox) * n_pts_bras[-x])))
      }, numeric(2))
    }
    # Power prior for efficacy and toxicity, on arms not significantly different
    PPEff <- pbeta(phi_eff, prior_eff + n_eff + A0_eff * n_effpts_autres[1, ], 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_effpts_autres[2, ] - n_effpts_autres[1, ]))
    data$est_eff[data$nb_ana == i] <- (prior_eff + n_eff + A0_eff * n_effpts_autres[1, ]) / (Nb_pts + A0_eff * n_effpts_autres[2, ] + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff + A0_eff * n_effpts_autres[1, ], 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_effpts_autres[2, ] - n_effpts_autres[1, ]))
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff + A0_eff * n_effpts_autres[1, ], 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_effpts_autres[2, ] - n_effpts_autres[1, ]))
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    data$est_tox[data$nb_ana == i] <- (prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ]) / (Nb_pts + A0_tox * n_toxpts_autres[2, ] + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
  }
  return(data)
}


## Power prior BOP with binomial/fisher test with borrowing relative to pvalue ----

### Toxicity only ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on toxicity (dynamic borrowing with pvalue of binomial test)
# Adapted from test-then-pool strategy from Viele et al. Use of historical control data for assessing treatment effects in clinical trials. Pharm Stat. 2014 Jan-Feb;13(1):41-54. 
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0 = exponent of the power prior to determine borrowing (downscale after borrowing)
## Tox0 = reference value that we test agains in binomial test
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_borrow_test_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0, Tox0, 
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
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        return(c(sum(n_tox[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
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
    data$est_eff[data$nb_ana == i] <- (prior_eff + n_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    # Borrowing information
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$est_tox[data$nb_ana == i] <- (prior_tox + n_tox + A0 * n_toxpts_autres[1, ]) / (Nb_pts + A0 * n_toxpts_autres[2, ] + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0 * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  return(data)
}

### Efficacy and toxicity ----

# Analyse the data from 1 trial using multi-arm BOP2 with power prior on toxicity and efficacy (dynamic borrowing with pvalue of binomial test)
# Adapted from test-then-pool strategy from Viele et al. Use of historical control data for assessing treatment effects in clinical trials. Pharm Stat. 2014 Jan-Feb;13(1):41-54. 
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## A0_eff, A0_tox = exponent of the power prior to determine borrowing (downscale after borrowing)
## Eff0, Tox0 = reference value that we test agains in binomial test
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff, prior_tox = vectors of length 2 giving the 2 parameters for the beta priors for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_bop_borrow_test_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                              A0_eff, Eff0, A0_tox, Tox0, 
                                              phi_eff, phi_tox, prior_eff, prior_tox) {
  data$arret_eff <- data$arret_tox <- NA_integer_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      n_effpts_autres <- vapply(seq_along(n_eff), \(x) {
        Pvals <- vapply(seq_along(n_eff)[-x], \(t) binom.test(n_eff[t], n_pts_bras[t], Eff0)$p.value, numeric(1))
        return(c(sum(n_eff[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        return(c(sum(n_tox[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
    } else { # If not at first analysis, only update numbers when the arm is not stopped
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      n_effpts_autres <- vapply(seq_along(n_eff), \(x) {
        Pvals <- vapply(seq_along(n_eff)[-x], \(t) binom.test(n_eff[t], n_pts_bras[t], Eff0)$p.value, numeric(1))
        return(c(sum(n_eff[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
      n_toxpts_autres <- vapply(seq_along(n_tox), \(x) {
        Pvals <- vapply(seq_along(n_tox)[-x], \(t) binom.test(n_tox[t], n_pts_bras[t], Tox0)$p.value, numeric(1))
        return(c(sum(n_tox[-x] * Pvals), sum(Pvals * n_pts_bras[-x])))
      }, numeric(2))
    }
    PPEff <- pbeta(phi_eff, prior_eff + n_eff + A0_eff * n_effpts_autres[1, ], 1 - prior_eff + Nb_pts - n_eff + A0_eff * (n_effpts_autres[2, ] - n_effpts_autres[1, ]))
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    PPTox <- 1 - pbeta(phi_tox, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$est_tox[data$nb_ana == i] <- (prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ]) / (Nb_pts + A0 * n_toxpts_autres[2, ] + 1) 
    data$icinf_tox[data$nb_ana == i] <- qbeta(.025, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0 * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    data$icsup_tox[data$nb_ana == i] <- qbeta(.975, prior_tox + n_tox + A0_tox * n_toxpts_autres[1, ], 1 - prior_tox + Nb_pts - n_tox + A0_tox * (n_toxpts_autres[2, ] - n_toxpts_autres[1, ]))
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  return(data)
}


## BHM ----

# In the following, for the functions to work, it is assumed that all STAN models are stored in a list
# nammed CompiledModelsTox for toxicity and CompiledModelsEffTox for both endpoints
# See script TestPriorBop_vXX.R for usage

### Toxicity only ----

# Compile the model BHM for toxicity
# Use of non central parametrization to reduce divergences in STAN (https://ballardtj.github.io/blog/non-centered-params/)
# Arguments :
## moy_mu, sigma_mu = mean and standard deviation of normal prior for mean proportion on logit scale
## moy_sig, sigma_sig = mean and standard deviation of normal truncated prior for standard deviation of proportion on logit scale
## moy_p, sigma_p = parameters for non central parametrization, not touching in general
# Output : Compiled BHM in STAN
CompilBHM_tox <- function(moy_mu = 0, moy_sig = 0, moy_p = 0,
                          sigma_mu = 5, sigma_sig = 5, sigma_p = 1) {
  ModeleHier_t <- "
  data {
    int<lower = 1> Nb; // Number of arms
    int<lower = 1> n[Nb]; // Number of patients per arm
    int<lower = 0> y[Nb]; // Number of toxicities per arm
  }
  parameters{
    real p_raw[Nb];
    real mu;
    real<lower = 0> sigma;
  }
  transformed parameters{
    real p[Nb]; // Probability for each arm
    real logit_p[Nb]; // Linear predictor in each arm
    for(j in 1:Nb){
      logit_p[j] = mu + sigma * p_raw[j]; // Non central parametrization
    }
    p = inv_logit(logit_p);
  }
  model{
  
    // prior distributions
    sigma ~ normal(%moy_sig%, %sigma_sig%);
    mu  ~ normal(%moy_mu%, %sigma_mu%);
  
    for (i in 1:Nb) {
      p_raw[i] ~ normal(%moy_p%, %sigma_p%);
      // binomial likelihood
      y[i] ~ binomial(n[i], p[i]);
    }
  }
"
  ModeleHier_t <- gsub("%moy_mu%", moy_mu, ModeleHier_t)
  ModeleHier_t <- gsub("%moy_sig%", moy_sig, ModeleHier_t)
  ModeleHier_t <- gsub("%moy_p%", moy_p, ModeleHier_t)
  ModeleHier_t <- gsub("%sigma_mu%", sigma_mu, ModeleHier_t)
  ModeleHier_t <- gsub("%sigma_sig%", sigma_sig, ModeleHier_t)
  ModeleHier_t <- gsub("%sigma_p%", sigma_p, ModeleHier_t)
  return(stan_model(model_code = ModeleHier_t))
}

# Analyse the data from 1 trial using BHM only for toxicity and BOP for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff = vectors of length 2 giving the 2 parameters for the beta priors for efficacy
# Output : data.frame of the trial with decisions from design
real_essai_modhier_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    # Get posterior distributions via MCMC
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      DonneesTox <- list(Nb = length(n_tox), n = n_pts_bras, y = n_tox)
      SampledTox <- sampling(CompiledModelsTox[["hBOP_tox"]], 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 2000,          
                             iter = 4000,
                             thin = 1,
                             cores = 3,
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      PPred <- extract(SampledTox, pars = "p")$p
      PPTox <- colMeans(PPred > phi_tox) # posterior probability
      data$est_tox[data$nb_ana == i] <- colMeans(PPred)
      data$icinf_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}

### Efficacy and toxicity ----

# Compile the model BHM for toxicity and efficacy
# Use of non central parametrization to reduce divergences in STAN (https://ballardtj.github.io/blog/non-centered-params/)
# Arguments :
## moy_mu_eff, sigma_mu_eff, moy_mu_tox, sigma_mu_tox = mean and standard deviation of normal prior for mean proportion on logit scale
## moy_sig_eff, sigma_sig_eff, moy_sig_tox, sigma_sig_tox = mean and standard deviation of normal truncated prior for standard deviation of proportion on logit scale
## moy_p_eff, sigma_p_eff, moy_p_tox, sigma_p_tox = parameters for non central parametrization, not touching in general
# Output : Compiled BHM in STAN
CompilBHM_efftox <- function(moy_mu_eff = 0, moy_sig_eff = 0, moy_p_eff = 0, moy_mu_tox = 0, moy_sig_tox = 0, moy_p_tox = 0,
                             sigma_mu_eff = 5, sigma_sig_eff = 5, sigma_p_eff = 1, sigma_mu_tox = 5, sigma_sig_tox = 5, sigma_p_tox = 1) {
  ModeleHier_et <- "
  data {
    int<lower = 1> Nb; // Nombre de bras
    int<lower = 1> n[Nb]; // Nombre de patients par bras
    int<lower = 0> y_tox[Nb]; // Nombre de toxicités par bras
    int<lower = 0> y_eff[Nb]; // Nombre de réponses par bras
  }
  parameters{
    real p_raw_tox[Nb];
    real p_raw_eff[Nb];
    real mu_tox;
    real mu_eff;
    real<lower = 0> sigma_tox;
    real<lower = 0> sigma_eff;
  }
  transformed parameters{
    real p_tox[Nb];
    real p_eff[Nb]; // Probabilité pour chaque bras d'efficacité et de toxicité
    real logit_p_tox[Nb];
    real logit_p_eff[Nb]; // Prédicteur linéaire pour chaque bras d'efficacité et de toxicité
    for(j in 1:Nb){
      logit_p_tox[j] = mu_tox + sigma_tox * p_raw_tox[j]; // Non central parametrization
      logit_p_eff[j] = mu_eff + sigma_eff * p_raw_eff[j]; 
    }
    p_tox = inv_logit(logit_p_tox);
    p_eff = inv_logit(logit_p_eff);
  }
  model{
    // prior distributions
    sigma_tox ~ normal(%moy_sig_tox%, %sigma_sig_tox%);
    sigma_eff ~ normal(%moy_sig_eff%, %sigma_sig_eff%);
    mu_tox  ~ normal(%moy_mu_tox%, %sigma_mu_tox%);
    mu_eff  ~ normal(%moy_mu_eff%, %sigma_mu_eff%);
  
    for (i in 1:Nb) {
      p_raw_tox[i] ~ normal(%moy_p_tox%, %sigma_p_tox%);
      p_raw_eff[i] ~ normal(%moy_p_eff%, %sigma_p_eff%);
      // binomial likelihood
      y_tox[i] ~ binomial(n[i], p_tox[i]);
      y_eff[i] ~ binomial(n[i], p_eff[i]);
    }
  }
"
  ModeleHier_et <- gsub("%moy_mu_tox%", moy_mu_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%moy_sig_tox%", moy_sig_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%moy_p_tox%", moy_p_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_mu_tox%", sigma_mu_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_sig_tox%", sigma_sig_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_p_tox%", sigma_p_tox, ModeleHier_et)
  ModeleHier_et <- gsub("%moy_mu_eff%", moy_mu_eff, ModeleHier_et)
  ModeleHier_et <- gsub("%moy_sig_eff%", moy_sig_eff, ModeleHier_et)
  ModeleHier_et <- gsub("%moy_p_eff%", moy_p_eff, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_mu_eff%", sigma_mu_eff, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_sig_eff%", sigma_sig_eff, ModeleHier_et)
  ModeleHier_et <- gsub("%sigma_p_eff%", sigma_p_eff, ModeleHier_et)
  return(stan_model(model_code = ModeleHier_et))
}

# Analyse the data from 1 trial using BHM only for toxicity and efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
# Output : data.frame of the trial with decisions from design
real_essai_modhier_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                      phi_eff, phi_tox) {
  
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
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      DonneesEffTox <- list(Nb = length(n_tox), n = n_pts_bras, y_tox = n_tox, y_eff = n_eff)
      SampledEffTox <- sampling(CompiledModelsEffTox[["hBOP_efftox"]], 
                                data = DonneesEffTox,         
                                chains = 3,             
                                warmup = 2000,          
                                iter = 10000,
                                thin = 4,
                                cores = 3,
                                control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                                seed = 121221)
      PPredTox <- extract(SampledEffTox, pars = "p_tox")$p_tox
      PPTox <- colMeans(PPredTox > phi_tox) 
      data$est_tox[data$nb_ana == i] <- colMeans(PPredTox)
      data$icinf_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .975)
      PPredEff <- extract(SampledEffTox, pars = "p_eff")$p_eff
      PPEff <- colMeans(PPredEff < phi_eff) 
      data$est_eff[data$nb_ana == i] <- colMeans(PPredEff)
      data$icinf_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .025)
      data$icsup_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      PPEff <- rep(1.5, length(n_eff))
      data$est_eff[data$nb_ana == i] <- data$est_eff[data$nb_ana == (i - 1)] 
      data$icinf_eff[data$nb_ana == i] <- data$icinf_eff[data$nb_ana == (i - 1)] 
      data$icsup_eff[data$nb_ana == i] <- data$icsup_eff[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


## CBHM ----

### Calibration of parameters ----

# Calibration of the 2 parameters for CBHM
# As described in Chu Y, Yuan Y. A Bayesian basket trial design using a calibrated Bayesian hierarchical model. Clin Trials. 2018 Apr;15(2):149-158.
# Arguments :
## NSimu = number of simulated trials
## NPts = number of patients per arm in each simulated trial
## Q0, Q1 = null and alternative hypotheses
## Graine = seed to ensure reproducibility
CalibrateCBHM <- function(NSimu = 10000, NPts, Q0, Q1, Graine = 121221) {
  
  set.seed(Graine)
  on.exit(set.seed(NULL), add = TRUE)
  NBras <- length(NPts)
  
  # Treatment effective for all subgroups
  MatHomo <- do.call("cbind", lapply(seq_len(NBras), \(x) rbinom(NSimu, NPts[x], Q1)))
  N <- sum(NPts)
  SommeCols <- NPts
  H_b <- apply(MatHomo, 1, \(r) {
    suppressWarnings(Statistic <- as.numeric(chisq.test(matrix(c(r, NPts - r), nrow = 2, byrow = TRUE))$statistic))
    if (Statistic < 1) Statistic <- 1
    return(Statistic)
  })
  H_b <- median(H_b)
  
  # Heterogeneous treatments
  H_bBarreTot <- sapply(seq(1, NBras - 1), \(j) {
    MatHetero <- do.call("cbind", lapply(seq_len(NBras), \(x) {
      if (x <= j) {
        rbinom(NSimu, NPts[x], Q1)
      } else {
        rbinom(NSimu, NPts[x], Q0)
      }
    }))
    H_bBarre <- apply(MatHetero, 1, \(r) {
      suppressWarnings(Statistic <- as.numeric(chisq.test(matrix(c(r, NPts - r), nrow = 2, byrow = TRUE))$statistic))
      if (Statistic < 1) Statistic <- 1
      return(Statistic)
    })
    return(median(H_bBarre))
  })
  H_bBarreTot <- min(H_bBarreTot)
  
  # Get the parameters a and b
  SigB <- 1
  SigBBarre <- 80
  a <- log(SigB) - (log(SigBBarre) - log(SigB)) / (log(H_bBarreTot) - log(H_b)) * log(H_b)
  b <- (log(SigBBarre) - log(SigB)) / (log(H_bBarreTot) - log(H_b))
  
  return(list("a" = a, "b" = b))
  
}
# CalibrateCBHM(NSimu = 10000, NPts = c(30, 30, 30, 30), Q0 = .2, Q1 = .35, Graine = 121221)

### Toxicity only ----

# Compile the model CBHM for toxicity
# Use of non central parametrization to reduce divergences in STAN (https://ballardtj.github.io/blog/non-centered-params/)
# Arguments :
## moy_mu, sigma_mu = mean and standard deviation of normal prior for mean proportion on logit scale
## moy_p, sigma_p = parameters for non central parametrization, not touching in general
# Output : Compiled CBHM in STAN
CompilCBHM_tox <- function(moy_mu = 0, moy_p = 0,
                           sigma_mu = 5, sigma_p = 1) {
  ModeleCBHM_t <- "
  data {
    int<lower = 1> Nb; // Number of arms
    int<lower = 1> n[Nb]; // Number of patients in each arm
    int<lower = 0> y[Nb]; // Number of toxicities in each arm
    real<lower = 0> etype; // SD calibrated
  }
  
  parameters{
    real p_raw[Nb];
    real mu;
  }
  
  transformed parameters{
    real p[Nb]; // Probability in each arm
    real logit_p[Nb]; // Linear predictor in each arm
    for(j in 1:Nb){
      logit_p[j] = mu + etype * p_raw[j]; // Non central parametrization
    }
    p = inv_logit(logit_p);
  }
  
  model{
  
    // prior distributions
    mu  ~ normal(%moy_mu%, %sigma_mu%);
  
    for (i in 1:Nb) {
      p_raw[i] ~ normal(%moy_p%, %sigma_p%);
      // binomial likelihood
      y[i] ~ binomial(n[i], p[i]);
    }
   
  }
"
  ModeleCBHM_t <- gsub("%moy_mu%", moy_mu, ModeleCBHM_t)
  ModeleCBHM_t <- gsub("%moy_p%", moy_p, ModeleCBHM_t)
  ModeleCBHM_t <- gsub("%sigma_mu%", sigma_mu, ModeleCBHM_t)
  ModeleCBHM_t <- gsub("%sigma_p%", sigma_p, ModeleCBHM_t)
  return(stan_model(model_code = ModeleCBHM_t))
}

# Analyse the data from 1 trial using CBHM only for toxicity and BOP for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff = vector of 2 parameters for the beta prior on efficacy
## a, b = 2 hyperparameters of CBHM to compute variance of the model
# Output : data.frame of the trial with decisions from design
real_essai_modcbhm_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                   phi_eff, phi_tox, prior_eff,
                                   a, b) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  data$var_tox <- NA_real_
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      suppressWarnings(Statistic <- as.numeric(chisq.test(matrix(c(n_tox, n_pts_bras - n_tox), nrow = 2, byrow = TRUE))$statistic))
      if (Statistic < 1) Statistic <- 1 # Security added in the code of the article to ensure not to have too small variance or too high
      VarianceCBHM <- exp(a + b * log(Statistic))
      if (is.na(VarianceCBHM) || !is.finite(VarianceCBHM) || VarianceCBHM > 1e+4) VarianceCBHM <- 1e+4
      DonneesTox <- list(Nb = length(n_tox), n = n_pts_bras, y = n_tox, etype = sqrt(VarianceCBHM))
      SampledTox <- sampling(CompiledModelsTox[["cbhmBOP_tox"]], 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 2000,          
                             iter = 4000,
                             thin = 1,
                             cores = 3,
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      PPred <- extract(SampledTox, pars = "p")$p
      PPTox <- colMeans(PPred > phi_tox) 
      data$est_tox[data$nb_ana == i] <- colMeans(PPred)
      data$icinf_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .975)
      data$var_tox[data$nb_ana == i] <- VarianceCBHM
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      data$var_tox[data$nb_ana == i] <- data$var_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}

### Efficacy and toxicity ----

# Compile the model CBHM for toxicity and efficacy
# Use of non central parametrization to reduce divergences in STAN (https://ballardtj.github.io/blog/non-centered-params/)
# Arguments :
## moy_mu_eff, sigma_mu_eff, moy_mu_tox, sigma_mu_tox = mean and standard deviation of normal prior for mean proportion on logit scale
## moy_p_eff, sigma_p_eff, moy_p_tox, sigma_p_tox = parameters for non central parametrization, not touching in general
# Output : Compiled CBHM in STAN
CompilCBHM_efftox <- function(moy_mu_tox = 0, moy_p_tox = 0, sigma_mu_tox = 5, sigma_p_tox = 1,
                              moy_mu_eff = 0, moy_p_eff = 0, sigma_mu_eff = 5, sigma_p_eff = 1) {
  ModeleCBHM_et <- "
  data {
    int<lower = 1> Nb; // Nombre de bras
    int<lower = 1> n[Nb]; // Nombre de patients par bras
    int<lower = 0> y_tox[Nb]; // Nombre de toxicités par bras
    real<lower = 0> etype_tox; // La variance du CBHM pour la toxicité
    int<lower = 0> y_eff[Nb]; // Nombre de réponses par bras
    real<lower = 0> etype_eff; // La variance du CBHM pour l'efficacité
  }
  
  parameters{
    real p_raw_eff[Nb];
    real mu_eff;
    real p_raw_tox[Nb];
    real mu_tox;
  }
  
  transformed parameters{
    real p_eff[Nb]; // Probabilité d'efficacité pour chaque bras
    real p_tox[Nb]; // Probabilité de toxicité pour chaque bras
    real logit_p_eff[Nb]; // Prédicteur linéaire pour l'efficacité dans chaque bras
    real logit_p_tox[Nb]; // Prédicteur linéaire pour la toxicité dans chaque bras
    for(j in 1:Nb){
      logit_p_eff[j] = mu_eff + etype_eff * p_raw_eff[j]; // Non central parametrization
      logit_p_tox[j] = mu_tox + etype_tox * p_raw_tox[j]; 
    }
    p_eff = inv_logit(logit_p_eff);
    p_tox = inv_logit(logit_p_tox);
  }
  
  model{
  
    // prior distributions
    mu_eff  ~ normal(%moy_mu_eff%, %sigma_mu_eff%);
    mu_tox  ~ normal(%moy_mu_tox%, %sigma_mu_tox%);
  
    for (i in 1:Nb) {
      p_raw_eff[i] ~ normal(%moy_p_eff%, %sigma_p_eff%);
      p_raw_tox[i] ~ normal(%moy_p_tox%, %sigma_p_tox%);
      // binomial likelihood
      y_eff[i] ~ binomial(n[i], p_eff[i]);
      y_tox[i] ~ binomial(n[i], p_tox[i]);
    }
   
  }
"
  ModeleCBHM_et <- gsub("%moy_mu_eff%", moy_mu_eff, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%moy_p_eff%", moy_p_eff, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%sigma_mu_eff%", sigma_mu_eff, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%sigma_p_eff%", sigma_p_eff, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%moy_mu_tox%", moy_mu_tox, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%moy_p_tox%", moy_p_tox, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%sigma_mu_tox%", sigma_mu_tox, ModeleCBHM_et)
  ModeleCBHM_et <- gsub("%sigma_p_tox%", sigma_p_tox, ModeleCBHM_et)
  return(stan_model(model_code = ModeleCBHM_et))
}

# Analyse the data from 1 trial using CBHM only for toxicity and efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## a_eff, b_eff, a_tox, b_tox = 2 hyperparameters of CBHM to compute variance of the model (separate optimisation for efficacy and toxicity)
# Output : data.frame of the trial with decisions from design
real_essai_modcbhm_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                      phi_eff, phi_tox,
                                      a_eff, b_eff, a_tox, b_tox) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  data$var_tox <- data$var_eff <- NA_real_
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
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      suppressWarnings(StatisticEff <- as.numeric(chisq.test(matrix(c(n_eff, n_pts_bras - n_eff), nrow = 2, byrow = TRUE))$statistic))
      suppressWarnings(StatisticTox <- as.numeric(chisq.test(matrix(c(n_tox, n_pts_bras - n_tox), nrow = 2, byrow = TRUE))$statistic))
      if (StatisticEff < 1) StatisticEff <- 1 # Security added in the code of the article
      if (StatisticTox < 1) StatisticTox <- 1
      VarianceCBHMEff <- exp(a_eff + b_eff * log(StatisticEff))
      VarianceCBHMTox <- exp(a_tox + b_tox * log(StatisticTox))
      if (is.na(VarianceCBHMEff) || !is.finite(VarianceCBHMEff) || VarianceCBHMEff > 1e+4) VarianceCBHMEff <- 1e+4 # If too large variance, problems for fitting models
      if (is.na(VarianceCBHMTox) || !is.finite(VarianceCBHMTox) || VarianceCBHMTox > 1e+4) VarianceCBHMTox <- 1e+4
      DonneesEffTox <- list(Nb = length(n_tox), n = n_pts_bras, y_tox = n_tox, etype_tox = sqrt(VarianceCBHMTox), y_eff = n_eff, etype_eff = sqrt(VarianceCBHMEff))
      SampledEffTox <- sampling(CompiledModelsEffTox[["cbhmBOP_efftox"]], 
                                data = DonneesEffTox,         
                                chains = 3,             
                                warmup = 2000,          
                                iter = 4000,
                                thin = 1,
                                cores = 3,
                                control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                                seed = 121221)
      PPredTox <- extract(SampledEffTox, pars = "p_tox")$p_tox
      PPTox <- colMeans(PPredTox > phi_tox) 
      data$est_tox[data$nb_ana == i] <- colMeans(PPredTox)
      data$icinf_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .975)
      PPredEff <- extract(SampledEffTox, pars = "p_eff")$p_eff
      PPEff <- colMeans(PPredEff < phi_eff) 
      data$est_eff[data$nb_ana == i] <- colMeans(PPredEff)
      data$icinf_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .025)
      data$icsup_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .975)
      data$var_tox[data$nb_ana == i] <- VarianceCBHMTox
      data$var_eff[data$nb_ana == i] <- VarianceCBHMEff
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      PPEff <- rep(1.5, length(n_eff))
      data$est_eff[data$nb_ana == i] <- data$est_eff[data$nb_ana == (i - 1)] 
      data$icinf_eff[data$nb_ana == i] <- data$icinf_eff[data$nb_ana == (i - 1)] 
      data$icsup_eff[data$nb_ana == i] <- data$icsup_eff[data$nb_ana == (i - 1)] 
      data$var_eff[data$nb_ana == i] <- data$var_eff[data$nb_ana == (i - 1)] 
      data$var_tox[data$nb_ana == i] <- data$var_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


## EXNEX ----

### Toxicity only ----

# Compile the model EXNEX for toxicity
# Arguments :
## moy_muex, sigma_muex = mean and standard deviation of normal prior for mean proportion on logit scale in EX part of the model
## moy_sigex, sigma_sigex = mean and standard deviation of normal prior for standard deviation of proportion on logit scale in EX part of the model
## mu_nex, sigma_nex = mean and standard deviation of prior for NEX part of model
# Output : Compiled CBHM in STAN
CompilEXNEX_tox <- function(moy_muex = -1.734601, moy_sigex = 0,
                            sigma_muex = 2.616, sigma_sigex = 1,
                            mu_nex = -1.734601, sigma_nex = 2.801) {
  ModeleEXNEX_t <- "
  data {
    int<lower = 1> Nb; // Number of arms
    int<lower = 1> n[Nb]; // Number of patients in each arm
    int<lower = 0> y[Nb]; // Number of toxicities in each arm
    real<lower = 0, upper = 1> prob_ex[Nb];
  }
  parameters{
    real theta[Nb]; // Log(odds) 
    real mu_ex; // Mean of log(odds) in EX distribution
    real<lower = 0> etype_ex; // Sd of log(odds) in EX distribution
  }
  transformed parameters{
    real p[Nb]; // Probability in each arm
    p = inv_logit(theta); // Theta = linear predictor in each arm
  }
  model{
    // prior distributions
    etype_ex ~ normal(%moy_sigex%, %sigma_sigex%);
    mu_ex  ~ normal(%moy_muex%, %sigma_muex%);
    for (i in 1:Nb) {
      target += log_mix(prob_ex[i],
                        normal_lpdf(theta[i] | mu_ex, etype_ex),
                        normal_lpdf(theta[i] | %mu_nex%, %sigma_nex%));
      // binomial likelihood
      y[i] ~ binomial(n[i], p[i]);
    }
  }
  generated quantities{
    // Probability of being EX (https://mc-stan.org/docs/stan-users-guide/finite-mixtures.html)
    real postprob_ex[Nb];
    for (i in 1:Nb) {
      postprob_ex[i] = normal_lpdf(theta[i] | mu_ex, etype_ex) + bernoulli_lpmf(0 | prob_ex[i]) - 
        log_sum_exp(normal_lpdf(theta[i] | mu_ex, etype_ex) + bernoulli_lpmf(0 | prob_ex[i]),
                    normal_lpdf(theta[i] | %mu_nex%, %sigma_nex%) + bernoulli_lpmf(1 | prob_ex[i]));
    }
    postprob_ex = exp(postprob_ex);
  }
"
  ModeleEXNEX_t <- gsub("%moy_muex%", moy_muex, ModeleEXNEX_t)
  ModeleEXNEX_t <- gsub("%moy_sigex%", moy_sigex, ModeleEXNEX_t)
  ModeleEXNEX_t <- gsub("%sigma_muex%", sigma_muex, ModeleEXNEX_t)
  ModeleEXNEX_t <- gsub("%sigma_sigex%", sigma_sigex, ModeleEXNEX_t)
  ModeleEXNEX_t <- gsub("%mu_nex%", mu_nex, ModeleEXNEX_t)
  ModeleEXNEX_t <- gsub("%sigma_nex%", sigma_nex, ModeleEXNEX_t)
  return(stan_model(model_code = ModeleEXNEX_t))
}

# Analyse the data from 1 trial using EXNEX only for toxicity and BOP for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff = vector of 2 parameters for the beta prior on efficacy
## pmix = initial guess of the probability of being EX
# Output : data.frame of the trial with decisions from design
real_essai_modexnex_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                    phi_eff, phi_tox, prior_eff, p_mix) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  data$prob_ex <- NA_real_
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      DonneesTox <- list(Nb = length(n_tox), n = n_pts_bras, y = n_tox, prob_ex = p_mix)
      SampledTox <- sampling(CompiledModelsTox[["exnexBOP_tox"]], 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 2000,          
                             iter = 18000,
                             thin = 5,
                             cores = 3,
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      PPred <- extract(SampledTox, pars = "p")$p
      PPTox <- colMeans(PPred > phi_tox) # On va chercher les distributions a posteriori de chaque bras pour la proportion de toxicité
      data$est_tox[data$nb_ana == i] <- colMeans(PPred)
      data$icinf_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPred, 2, quantile, probs = .975)
      data$prob_ex[data$nb_ana == i] <- colMeans(extract(SampledTox, pars = "postprob_ex")$postprob_ex)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      data$prob_ex[data$nb_ana == i] <- data$prob_ex[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


### Efficacy and toxicity ----

# Compile the model EXNEX for toxicity and efficacy
# Arguments :
## moy_muex_tox, sigma_muex_tox, moy_muex_eff, sigma_muex_eff = mean and standard deviation of normal prior for mean proportion on logit scale in EX part of the model
## moy_sigex_tox, sigma_sigex_tox, moy_sigex_eff, sigma_sigex_eff = mean and standard deviation of normal prior for standard deviation of proportion on logit scale in EX part of the model
## mu_nex_tox, sigma_nex_tox, mu_nex_eff, sigma_nex_eff = mean and standard deviation of prior for NEX part of model
# Output : Compiled CBHM in STAN
CompilEXNEX_efftox <- function(moy_muex_tox = -1.734601, moy_sigex_tox = 0, sigma_muex_tox = 2.616, sigma_sigex_tox = 1, mu_nex_tox = -1.734601, sigma_nex_tox = 2.801,
                               moy_muex_eff = -1.734601, moy_sigex_eff = 0, sigma_muex_eff = 2.616, sigma_sigex_eff = 1, mu_nex_eff = -1.734601, sigma_nex_eff = 2.801) {
  ModeleEXNEX_et <- "
  data {
    int<lower = 1> Nb; // Number of arms
    int<lower = 1> n[Nb]; // Number of patients in each arm
    int<lower = 0> y_eff[Nb]; // Number of responses in each arm
    int<lower = 0> y_tox[Nb]; // Number of toxicities in each arm
    real<lower = 0, upper = 1> prob_ex_eff[Nb];
    real<lower = 0, upper = 1> prob_ex_tox[Nb];
  }
  parameters{
    real theta_eff[Nb]; // Log(odds) of efficacy in each arm
    real theta_tox[Nb]; // Log(odds) of toxicity in each arm
    real mu_ex_eff; // Mean of log(odds) of efficacy in EX
    real<lower = 0> etype_ex_eff; // Sd of log(odds) of efficacy in EX
    real mu_ex_tox; // Mean of log(odds) of toxicity in EX
    real<lower = 0> etype_ex_tox; // Sd of log(odds) of toxicity in EX
  }
  transformed parameters{
    real p_eff[Nb]; // Probability of response in each arm
    real p_tox[Nb]; // Probability of toxicity in each arm
    p_eff = inv_logit(theta_eff); // Theta = linear predictor in each arm
    p_tox = inv_logit(theta_tox); 
  }
  model{
    // prior distributions
    etype_ex_eff ~ normal(%moy_sigex_eff%, %sigma_sigex_eff%);
    mu_ex_eff  ~ normal(%moy_muex_eff%, %sigma_muex_eff%);
    etype_ex_tox ~ normal(%moy_sigex_tox%, %sigma_sigex_tox%);
    mu_ex_tox  ~ normal(%moy_muex_tox%, %sigma_muex_tox%);
    for (i in 1:Nb) {
      target += log_mix(prob_ex_eff[i],
                        normal_lpdf(theta_eff[i] | mu_ex_eff, etype_ex_eff),
                        normal_lpdf(theta_eff[i] | %mu_nex_eff%, %sigma_nex_eff%));
      target += log_mix(prob_ex_tox[i],
                        normal_lpdf(theta_tox[i] | mu_ex_tox, etype_ex_tox),
                        normal_lpdf(theta_tox[i] | %mu_nex_tox%, %sigma_nex_tox%));
      // binomial likelihood
      y_eff[i] ~ binomial(n[i], p_eff[i]);
      y_tox[i] ~ binomial(n[i], p_tox[i]);
    }
  }
  generated quantities{
    // Probability of being EX (https://mc-stan.org/docs/stan-users-guide/finite-mixtures.html)
    real postprob_ex_eff[Nb];
    real postprob_ex_tox[Nb];
    for (i in 1:Nb) {
      postprob_ex_eff[i] = normal_lpdf(theta_eff[i] | mu_ex_eff, etype_ex_eff) + bernoulli_lpmf(0 | prob_ex_eff[i]) - 
        log_sum_exp(normal_lpdf(theta_eff[i] | mu_ex_eff, etype_ex_eff) + bernoulli_lpmf(0 | prob_ex_eff[i]),
                    normal_lpdf(theta_eff[i] | %mu_nex_eff%, %sigma_nex_eff%) + bernoulli_lpmf(1 | prob_ex_eff[i]));
      postprob_ex_tox[i] = normal_lpdf(theta_tox[i] | mu_ex_tox, etype_ex_tox) + bernoulli_lpmf(0 | prob_ex_tox[i]) - 
        log_sum_exp(normal_lpdf(theta_tox[i] | mu_ex_tox, etype_ex_tox) + bernoulli_lpmf(0 | prob_ex_tox[i]),
                    normal_lpdf(theta_tox[i] | %mu_nex_tox%, %sigma_nex_tox%) + bernoulli_lpmf(1 | prob_ex_tox[i]));
    }
    postprob_ex_eff = exp(postprob_ex_eff);
    postprob_ex_tox = exp(postprob_ex_tox);
  }
"
  ModeleEXNEX_et <- gsub("%moy_muex_eff%", moy_muex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%moy_sigex_eff%", moy_sigex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_muex_eff%", sigma_muex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_sigex_eff%", sigma_sigex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%mu_nex_eff%", mu_nex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_nex_eff%", sigma_nex_eff, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%moy_muex_tox%", moy_muex_tox, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%moy_sigex_tox%", moy_sigex_tox, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_muex_tox%", sigma_muex_tox, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_sigex_tox%", sigma_sigex_tox, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%mu_nex_tox%", mu_nex_tox, ModeleEXNEX_et)
  ModeleEXNEX_et <- gsub("%sigma_nex_tox%", sigma_nex_tox, ModeleEXNEX_et)
  return(stan_model(model_code = ModeleEXNEX_et))
}

# Analyse the data from 1 trial using EXNEX only for toxicity and for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## pmix_eff, pmix_tox = initial guess of the probability of being EX for efficacy and toxicity
# Output : data.frame of the trial with decisions from design
real_essai_modexnex_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                       phi_eff, phi_tox, p_mix_eff, p_mix_tox = p_mix_eff) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  data$prob_ex_tox <- data$prob_ex_eff <- NA_real_
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
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      DonneesEffTox <- list(Nb = length(n_tox), n = n_pts_bras, y_eff = n_eff, prob_ex_eff = p_mix_eff, y_tox = n_tox, prob_ex_tox = p_mix_tox)
      SampledEffTox <- sampling(CompiledModelsEffTox[["exnexBOP_efftox"]], 
                                data = DonneesEffTox,         
                                chains = 3,             
                                warmup = 2000,          
                                iter = 22000,
                                thin = 5,
                                cores = 3,
                                control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                                seed = 121221)
      PPredTox <- extract(SampledEffTox, pars = "p_tox")$p_tox
      PPTox <- colMeans(PPredTox > phi_tox) 
      data$est_tox[data$nb_ana == i] <- colMeans(PPredTox)
      data$icinf_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .025)
      data$icsup_tox[data$nb_ana == i] <- apply(PPredTox, 2, quantile, probs = .975)
      PPredEff <- extract(SampledEffTox, pars = "p_eff")$p_eff
      PPEff <- colMeans(PPredEff < phi_eff) 
      data$est_eff[data$nb_ana == i] <- colMeans(PPredEff)
      data$icinf_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .025)
      data$icsup_eff[data$nb_ana == i] <- apply(PPredEff, 2, quantile, probs = .975)
      data$prob_ex_eff[data$nb_ana == i] <- colMeans(extract(SampledEffTox, pars = "postprob_ex_eff")$postprob_ex_eff)
      data$prob_ex_tox[data$nb_ana == i] <- colMeans(extract(SampledEffTox, pars = "postprob_ex_tox")$postprob_ex_tox)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      PPEff <- rep(1.5, length(n_eff))
      data$est_eff[data$nb_ana == i] <- data$est_eff[data$nb_ana == (i - 1)] 
      data$icinf_eff[data$nb_ana == i] <- data$icinf_eff[data$nb_ana == (i - 1)] 
      data$icsup_eff[data$nb_ana == i] <- data$icsup_eff[data$nb_ana == (i - 1)] 
      data$prob_ex_eff[data$nb_ana == i] <- data$prob_ex_eff[data$nb_ana == (i - 1)] 
      data$prob_ex_tox[data$nb_ana == i] <- data$prob_ex_tox[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


## Logistic regression ----

# The different logistic regression models:
## - log1: classic logistic regression
## - log2: force beta to be positive
## - log3: positive informative prior for log(OR), but still high variance
## - log4: dosis in qualitative variable
## - log5: dosis in qualitative variable and force beta and gamma positive (so that p_dose3 > p_dose2 > p_dose1)
## - log6: dosis as continuous variable with a squared term

### Toxicity only ----

# Compile the logistic model for toxicity 
# Arguments :
## mu_inter, sigma_inter = mean and standard deviation of the normal prior for intercept
## mu_coef, sigma_coef = mean and standard deviation of the normal prior for beta coefficient
## PentePos = TRUE to ensure a positive log-linear beta coefficient
## SecondCov = TRUE to add a third coefficient in the model (with same prior as "coef")
# Output : Compiled logistic regression in STAN
CompilModLog_tox <- function(mu_inter = 0, sigma_inter = 5, mu_coef = 0, sigma_coef = 5,
                             PentePos = FALSE, SecondCov = FALSE) {
  ModeleLog <- "
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
    alpha ~ normal(%mu_inter%, %sigma_inter%);
    beta ~ normal(%mu_coef%, %sigma_coef%); 
    // Likelihood
    y ~ bernoulli_logit(alpha + beta * x);
  }
"
  if (SecondCov) {
    ModeleLog <- gsub("vector\\[N\\] x;", "vector\\[N\\] x;\n    vector\\[N\\] x_2;", ModeleLog)
    ModeleLog <- gsub("real beta;", "real beta;\n    real gamma;", ModeleLog)
    ModeleLog <- gsub("%sigma_coef%);", "%sigma_coef%);\n    gamma ~ normal(%mu_coef%, %sigma_coef%);", ModeleLog)
    ModeleLog <- gsub("beta \\* x", "beta \\* x + gamma \\* x_2", ModeleLog)
  } 
  if (PentePos) {
    ModeleLog <- gsub("real beta", "real<lower = 0> beta", ModeleLog)
    ModeleLog <- gsub("real gamma", "real<lower = 0> gamma", ModeleLog)
  } 
  ModeleLog <- gsub("%mu_inter%", mu_inter, ModeleLog)
  ModeleLog <- gsub("%sigma_inter%", sigma_inter, ModeleLog)
  ModeleLog <- gsub("%mu_coef%", mu_coef, ModeleLog)
  ModeleLog <- gsub("%sigma_coef%", sigma_coef, ModeleLog)
  return(stan_model(model_code = ModeleLog))
}

# Analyse the data from 1 trial using bayesian logistic regression for toxicity and for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff = vector of 2 parameters for the beta prior on efficacy
## modele_log = Number of model as defined above
## log_dose = TRUE to incorporate log of dose in model instead of dose
# Output : data.frame of the trial with decisions from design
real_essai_bayeslog_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                    phi_eff, phi_tox, prior_eff, modele_log, log_dose) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$dose <- as.numeric(gsub("^ttt(\\d+)$", "\\1", data$ttt)) # Doses 1/2/3/... prises dans ces simulations
  # For the simulations we made the hypothesis of 3 doses with proportionality 1/2/3 so replacing the dose by 1/2/3 isn't changing anything in our case
  # But if it is not the case, one may adapt this manner to determine doses in simulations
  if (log_dose) data$dose <- log(data$dose)
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      if (modele_log %in% c(1:3, 6)) {
        DonneesTox <- list(N = sum(n_pts_bras),
                           y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                           x = rep(doses_ttt, n_pts_bras))
        if (modele_log == 6) DonneesTox[["x_2"]] <- DonneesTox[["x"]] ** 2
      } else {
        DonneesTox <- list(N = sum(n_pts_bras),
                           y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                           x = rep(c(0, rep(1, length(doses_ttt) - 1)), n_pts_bras),
                           x_2 = rep(c(0, 0, rep(1, length(doses_ttt) - 2)), n_pts_bras))
      }
      SampledTox <- sampling(CompiledModelsTox[[paste0("log", modele_log, "_tox")]], 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 2000,          
                             iter = 15000,
                             thin = 8,
                             cores = 3,
                             # For convergence issues
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      DistPost <- extract(SampledTox)
      if (modele_log %in% c(1:3)) {
        PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + x * DistPost$beta))))
      } else if (modele_log == 6) {
        PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + x * DistPost$beta + x ** 2 * DistPost$gamma))))
      } else {
        PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + (x != 1) * DistPost$beta + (x == 3) * DistPost$gamma))))
      }
      PPTox <- vapply(PPred, \(distrib) {mean(distrib > phi_tox)}, double(1))
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
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}

### Efficacy and toxicity ----

# Compile the logistic model for toxicity and BOP2 for efficacy
# Arguments :
## mu_inter_tox, sigma_inter_tox, mu_inter_eff, sigma_inter_eff = mean and standard deviation of the normal prior for intercept for toxicity and efficacy
## mu_coef_tox, sigma_coef_tox, mu_coef_eff, sigma_coef_eff = mean and standard deviation of the normal prior for beta coefficient for toxicity and efficacy
## PentePos_tox, PentePos_eff = TRUE to ensure a positive log-linear beta coefficient for toxicity and efficacy
## SecondCov_tox, SecondCov_eff = TRUE to add a third coefficient in the model (with same prior as "coef") for toxicity and efficacy
# Output : Compiled logistic regression in STAN
CompilModLog_efftox <- function(mu_inter_tox = 0, sigma_inter_tox = 5, mu_coef_tox = 0, sigma_coef_tox = 5,
                                mu_inter_eff = 0, sigma_inter_eff = 5, mu_coef_eff = 0, sigma_coef_eff = 5,
                                PentePos_tox = FALSE, SecondCov_tox = FALSE,
                                PentePos_eff = FALSE, SecondCov_eff = FALSE) {
  ModeleLog <- "
  data {
    int<lower=0> N;
    vector[N] x;
    int<lower=0,upper=1> y_tox[N];
    int<lower=0,upper=1> y_eff[N];
  }
  parameters {
    real alpha_tox;
    real beta_tox;
    real alpha_eff;
    real beta_eff;
  }
  model {
    // Priors
    alpha_tox ~ normal(%mu_inter_tox%, %sigma_inter_tox%);
    beta_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%); 
    alpha_eff ~ normal(%mu_inter_eff%, %sigma_inter_eff%);
    beta_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%); 
    // Likelihood
    y_tox ~ bernoulli_logit(alpha_tox + beta_tox * x);
    y_eff ~ bernoulli_logit(alpha_eff + beta_eff * x);
  }
"
  if (SecondCov_eff | SecondCov_tox) {
    ModeleLog <- gsub("vector\\[N\\] x;", "vector\\[N\\] x;\n    vector\\[N\\] x_2;", ModeleLog)
    if (SecondCov_tox) {
      ModeleLog <- gsub("real beta_tox;", "real beta_tox;\n    real gamma_tox;", ModeleLog)
      ModeleLog <- gsub("%sigma_coef_tox%);", "%sigma_coef_tox%);\n    gamma_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%);", ModeleLog)
      ModeleLog <- gsub("beta_tox \\* x", "beta_tox \\* x + gamma_tox \\* x_2", ModeleLog)
    }
    if (SecondCov_eff) {
      ModeleLog <- gsub("real beta_eff;", "real beta_eff;\n    real gamma_eff;", ModeleLog)
      ModeleLog <- gsub("%sigma_coef_eff%);", "%sigma_coef_eff%);\n    gamma_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%);", ModeleLog)
      ModeleLog <- gsub("beta_eff \\* x", "beta_eff \\* x + gamma_eff \\* x_2", ModeleLog)
    }
  } 
  if (PentePos_tox) {
    ModeleLog <- gsub("real beta_tox", "real<lower = 0> beta_tox", ModeleLog)
    ModeleLog <- gsub("real gamma_tox", "real<lower = 0> gamma_tox", ModeleLog)
  } 
  if (PentePos_eff) {
    ModeleLog <- gsub("real beta_eff", "real<lower = 0> beta_eff", ModeleLog)
    ModeleLog <- gsub("real gamma_eff", "real<lower = 0> gamma_eff", ModeleLog)
  } 
  ModeleLog <- gsub("%mu_inter_tox%", mu_inter_tox, ModeleLog)
  ModeleLog <- gsub("%sigma_inter_tox%", sigma_inter_tox, ModeleLog)
  ModeleLog <- gsub("%mu_coef_tox%", mu_coef_tox, ModeleLog)
  ModeleLog <- gsub("%sigma_coef_tox%", sigma_coef_tox, ModeleLog)
  ModeleLog <- gsub("%mu_inter_eff%", mu_inter_eff, ModeleLog)
  ModeleLog <- gsub("%sigma_inter_eff%", sigma_inter_eff, ModeleLog)
  ModeleLog <- gsub("%mu_coef_eff%", mu_coef_eff, ModeleLog)
  ModeleLog <- gsub("%sigma_coef_eff%", sigma_coef_eff, ModeleLog)
  return(stan_model(model_code = ModeleLog))
}

# Analyse the data from 1 trial using bayesian logistic regression for toxicity and for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## modele_log = Number of model as defined above
## log_dose = TRUE to incorporate log of dose in model instead of dose
# Output : data.frame of the trial with decisions from design
real_essai_bayeslog_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                       phi_eff, phi_tox, modele_log, log_dose) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  data$dose <- as.numeric(gsub("^ttt(\\d+)$", "\\1", data$ttt)) # Doses 1/2/3/... prises dans ces simulations
  # For the simulations we made the hypothesis of 3 doses with proportionality 1/2/3 so replacing the dose by 1/2/3 isn't changing anything in our case
  # But if it is not the case, one may adapt this manner to determine doses in simulations
  if (log_dose) data$dose <- log(data$dose)
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
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      if (modele_log %in% c(1:3, 6)) {
        DonneesEffTox <- list(N = sum(n_pts_bras),
                              y_tox = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                              y_eff = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_eff[x], n_eff[x])))),
                              x = rep(doses_ttt, n_pts_bras))
        if (modele_log == 6) DonneesEffTox[["x_2"]] <- DonneesEffTox[["x"]] ** 2
      } else {
        DonneesEffTox <- list(N = sum(n_pts_bras),
                              y_tox = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                              y_eff = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_eff[x], n_eff[x])))),
                              x = rep(c(0, rep(1, length(doses_ttt) - 1)), n_pts_bras),
                              x_2 = rep(c(0, 0, rep(1, length(doses_ttt) - 2)), n_pts_bras))
      }
      SampledEffTox <- sampling(CompiledModelsEffTox[[paste0("log", modele_log, "_efftox")]], 
                                data = DonneesEffTox,         
                                chains = 3,             
                                warmup = 2000,          
                                iter = 10000,
                                thin = 4,
                                cores = 3,
                                # Convergence issues
                                control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                                seed = 121221)
      DistPost <- extract(SampledEffTox)
      if (modele_log %in% c(1:3)) {
        PPredTox <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_tox + x * DistPost$beta_tox))))
        PPredEff <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_eff + x * DistPost$beta_eff))))
      } else if (modele_log == 6) {
        PPredTox <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_tox + x * DistPost$beta_tox + x ** 2 * DistPost$gamma_tox))))
        PPredEff <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_eff + x * DistPost$beta_eff + x ** 2 * DistPost$gamma_eff))))
      } else {
        PPredTox <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_tox + (x != 1) * DistPost$beta_tox + (x == 3) * DistPost$gamma_tox))))
        PPredEff <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_eff + (x != 1) * DistPost$beta_eff + (x == 3) * DistPost$gamma_eff))))
      }
      PPTox <- vapply(PPredTox, \(distrib) {mean(distrib > phi_tox)}, double(1))
      data$est_tox[data$nb_ana == i] <- vapply(PPredTox, mean, double(1))
      data$icinf_tox[data$nb_ana == i] <- vapply(PPredTox, quantile, double(1), probs = .025)
      data$icsup_tox[data$nb_ana == i] <- vapply(PPredTox, quantile, double(1), probs = .975)
      PPEff <- vapply(PPredEff, \(distrib) {mean(distrib < phi_eff)}, double(1))
      data$est_eff[data$nb_ana == i] <- vapply(PPredEff, mean, double(1))
      data$icinf_eff[data$nb_ana == i] <- vapply(PPredEff, quantile, double(1), probs = .025)
      data$icsup_eff[data$nb_ana == i] <- vapply(PPredEff, quantile, double(1), probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      PPEff <- rep(1.5, length(n_eff))
      data$est_eff[data$nb_ana == i] <- data$est_eff[data$nb_ana == (i - 1)] 
      data$icinf_eff[data$nb_ana == i] <- data$icinf_eff[data$nb_ana == (i - 1)] 
      data$icsup_eff[data$nb_ana == i] <- data$icsup_eff[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


## Logistic regression with CRM skeleton ----

# Compute the skeleton of the model. Allow logistic and power skeleton
# Arguments :
## Probas = assumed probability of efficacy/toxicity assumed a priori
## A0 = assumed intercept for logistic model
## B0 = assumed slope for logistic model and coefficient for power model
## MTD = number of assumed MTD dose
## Delta = delta parameter from Cheung to incorporate variability around the doses
## Model = logistic or power to select the desired model
MakeSkeleton <- function(Probas, A0 = 3, B0 = 1, MTD, Delta = NULL, Model = "logistic") {
  Model <- match.arg(Model, c("logistic", "power"))
  if (Model == "logistic") {
    if (is.null(Delta)) {
      return(-(log((1 - Probas) / Probas) + A0) / B0)
    } else {
      DoseLogit <- numeric(length(Probas))
      DoseLogit[MTD] <- (log(Probas[MTD] / (1 - Probas[MTD])) - A0) / B0
      if (MTD > 1) {
        for(i in rev(seq_len(MTD)[-1])) {
          DoseLogit[i - 1] <- (log((Probas[MTD] - Delta) / (1 - (Probas[MTD] - Delta))) - A0) * DoseLogit[i] / (log((Probas[MTD] + Delta) / (1 - (Probas[MTD] + Delta))) - A0)
        }
      }
      if (MTD < length(Probas)) {
        for (i in seq(MTD, length(Probas) - 1)) {
          DoseLogit[i + 1] <- (log((Probas[MTD] + Delta) / (1 - (Probas[MTD] + Delta))) - A0) * DoseLogit[i] / (log((Probas[MTD] - Delta) / (1 - (Probas[MTD] - Delta))) - A0)
        }
      }
      return(DoseLogit)
    }
  } else if (Model == "power") {
    return(Probas ** (exp(-B0)))
  }
}

### Toxicity only ----

# Compile the logistic model type CRM for toxicity 
# Arguments :
## mu_inter, sigma_inter = mean and standard deviation of the normal prior for intercept
## mu_coef, sigma_coef = mean and standard deviation of the normal prior for beta coefficient
## fixed_intercept = TRUE to make the intercept a constant
## PentePos = TRUE to ensure a positive log-linear beta coefficient
## SecondCov = TRUE to add a third coefficient in the model (with same prior as "coef")
# Output : Compiled logistic regression in STAN
CompilModLogCrm_tox <- function(mu_inter = 0, sigma_inter = 5, mu_coef = 0, sigma_coef = 5,
                                fixed_intercept = FALSE, PentePos = FALSE, SecondCov = FALSE) {
  if (fixed_intercept) {
    ModeleLog <- "
  data {
    int<lower=0> N;
    vector[N] x;
    int<lower=0,upper=1> y[N];
    real alpha;
  }
  parameters {
    real beta;
  }
  model {
    // Priors
    beta ~ normal(%mu_coef%, %sigma_coef%); 
    // Likelihood
    y ~ bernoulli_logit(alpha + beta * x);
  }
"
  } else {
    ModeleLog <- "
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
    alpha ~ normal(%mu_inter%, %sigma_inter%);
    beta ~ normal(%mu_coef%, %sigma_coef%); 
    // Likelihood
    y ~ bernoulli_logit(alpha + beta * x);
  }
"
  }
  if (SecondCov) {
    ModeleLog <- gsub("vector\\[N\\] x;", "vector\\[N\\] x;\n    vector\\[N\\] x_2;", ModeleLog)
    ModeleLog <- gsub("real beta;", "real beta;\n    real gamma;", ModeleLog)
    ModeleLog <- gsub("%sigma_coef%);", "%sigma_coef%);\n    gamma ~ normal(%mu_coef%, %sigma_coef%);", ModeleLog)
    ModeleLog <- gsub("beta \\* x", "beta \\* x + gamma \\* x_2", ModeleLog)
  } 
  if (PentePos) {
    ModeleLog <- gsub("real beta", "real<lower = 0> beta", ModeleLog)
    ModeleLog <- gsub("real gamma", "real<lower = 0> gamma", ModeleLog)
  } 
  ModeleLog <- gsub("%mu_inter%", mu_inter, ModeleLog)
  ModeleLog <- gsub("%sigma_inter%", sigma_inter, ModeleLog)
  ModeleLog <- gsub("%mu_coef%", mu_coef, ModeleLog)
  ModeleLog <- gsub("%sigma_coef%", sigma_coef, ModeleLog)
  return(stan_model(model_code = ModeleLog))
}

# Compile the power model (type CRM) for toxicity 
# Arguments :
## mu_coef, sigma_coef = mean and standard deviation of the normal prior for beta coefficient
## PentePos = TRUE to ensure a positive beta coefficient
# Output : Compiled power model CRM in STAN
CompilModPowCrm_tox <- function(mu_coef = 0, sigma_coef = 5, PentePos = FALSE) {
  ModelePow <- "
  data {
    int<lower=0> N;
    vector[N] x;
    int<lower=0,upper=1> y[N];
  }
  parameters {
    real beta;
  }
  model {
    // Priors
    beta ~ normal(%mu_coef%, %sigma_coef%); 
    // Likelihood
    y ~ bernoulli(exp(log(x) * exp(beta))); // For some reason on server function pow doesn't work with vectors so little tricks with exp
  }
"
  if (PentePos) {
    ModelePow <- gsub("real beta", "real<lower = 0> beta", ModelePow)
  } 
  ModelePow <- gsub("%mu_coef%", mu_coef, ModelePow)
  ModelePow <- gsub("%sigma_coef%", sigma_coef, ModelePow)
  return(stan_model(model_code = ModelePow))
}

# Analyse the data from 1 trial using bayesian logistic regression with CRM skeleton for toxicity and BOP2 for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## prior_eff = vector of 2 parameters for the beta prior on efficacy
## skeleton_tox = skeleton for toxicity (for example using MakeSkeleton function)
## fixed_intercept = TRUE for 1 parameter logistic regression
## A0 = intercept when fixed
## power_mod = TRUE for power model
# Output : data.frame of the trial with decisions from design
real_essai_bayeslogcrm_tox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                       phi_eff, phi_tox, prior_eff, skeleton_tox, fixed_intercept, A0, 
                                       power_mod = FALSE) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  # Using the crm skeleton for doses
  data$dose <- as.numeric(gsub("^ttt(\\d+)$", "\\1", data$ttt)) 
  data$dose <- skeleton_tox[data$dose]
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
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
    data$est_eff[data$nb_ana == i] <- (n_eff + prior_eff) / (Nb_pts + 1) 
    data$icinf_eff[data$nb_ana == i] <- qbeta(.025, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff) 
    data$icsup_eff[data$nb_ana == i] <- qbeta(.975, prior_eff + n_eff, 1 - prior_eff + Nb_pts - n_eff)
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      if (!power_mod) {
        if (fixed_intercept) {
          DonneesTox <- list(N = sum(n_pts_bras),
                             y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                             x = rep(doses_ttt, n_pts_bras),
                             alpha = A0)
          Modele <- "crm_fixed_tox"
        } else {
          DonneesTox <- list(N = sum(n_pts_bras),
                             y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                             x = rep(doses_ttt, n_pts_bras))
          Modele <- "crm_unfixed_tox"
        }
      } else {
        DonneesTox <- list(N = sum(n_pts_bras),
                           y = rep(rep(0:1, length(doses_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                           x = rep(doses_ttt, n_pts_bras))
        Modele <- "crmpow_tox"
      }
      SampledTox <- sampling(CompiledModelsTox[[Modele]], 
                             data = DonneesTox,         
                             chains = 3,             
                             warmup = 2000,          
                             iter = 15000,
                             thin = 8,
                             cores = 3,
                             control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                             seed = 121221)
      DistPost <- rstan::extract(SampledTox)
      if (!power_mod) {
        if (fixed_intercept) {
          PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (A0 + x * DistPost$beta))))
        } else {
          PPred <- lapply(doses_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha + x * DistPost$beta))))
        }
      } else {
        PPred <- lapply(doses_ttt, \(x) x ** exp(DistPost$beta))
      }
      PPTox <- vapply(PPred, \(distrib) {mean(distrib > phi_tox)}, double(1))
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
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}

### Efficacy and toxicity ----

# Compile the logistic model (like CRM) for toxicity  and efficacy
# Arguments :
## mu_inter_tox, sigma_inter_tox, mu_inter_eff, sigma_inter_eff = mean and standard deviation of the normal prior for intercept for toxicity and efficacy
## mu_coef_tox, sigma_coef_tox, mu_coef_eff, sigma_coef_eff = mean and standard deviation of the normal prior for beta coefficient for toxicity and efficacy
## fixed_intercept = TRUE to make the intercept a constant
## PentePos_tox, PentePos_eff = TRUE to ensure a positive log-linear beta coefficient
## SecondCov_tox, SecondCov_eff = TRUE to add a third coefficient in the model (with same prior as "coef")
# Output : Compiled logistic regression in STAN
CompilModLogCrm_efftox <- function(mu_inter_tox = 0, sigma_inter_tox = 5, mu_coef_tox = 0, sigma_coef_tox = 5,
                                   mu_inter_eff = 0, sigma_inter_eff = 5, mu_coef_eff = 0, sigma_coef_eff = 5,
                                   fixed_intercept = FALSE, PentePos_tox = FALSE, SecondCov_tox = FALSE,
                                   PentePos_eff = FALSE, SecondCov_eff = FALSE) {
  if (fixed_intercept) {
    ModeleLog <- "
  data {
    int<lower=0> N;
    vector[N] x_eff;
    vector[N] x_tox;
    int<lower=0,upper=1> y_tox[N];
    int<lower=0,upper=1> y_eff[N];
    real alpha_eff;
    real alpha_tox;
  }
  parameters {
    real beta_tox;
    real beta_eff;
  }
  model {
    // Priors
    beta_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%); 
    beta_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%); 
    // Likelihood
    y_tox ~ bernoulli_logit(alpha_tox + beta_tox * x_tox);
    y_eff ~ bernoulli_logit(alpha_eff + beta_eff * x_eff);
  }
"
  } else {
    ModeleLog <- "
  data {
    int<lower=0> N;
    vector[N] x_eff;
    vector[N] x_tox;
    int<lower=0,upper=1> y_tox[N];
    int<lower=0,upper=1> y_eff[N];
  }
  parameters {
    real alpha_tox;
    real beta_tox;
    real alpha_eff;
    real beta_eff;
  }
  model {
    // Priors
    alpha_tox ~ normal(%mu_inter_tox%, %sigma_inter_tox%);
    beta_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%); 
    alpha_eff ~ normal(%mu_inter_eff%, %sigma_inter_eff%);
    beta_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%); 
    // Likelihood
    y_tox ~ bernoulli_logit(alpha_tox + beta_tox * x_tox);
    y_eff ~ bernoulli_logit(alpha_eff + beta_eff * x_eff);
  }
"
  }
  
  if (SecondCov_eff | SecondCov_tox) {
    ModeleLog <- gsub("vector\\[N\\] x;", "vector\\[N\\] x;\n    vector\\[N\\] x_2;", ModeleLog)
    if (SecondCov_tox) {
      ModeleLog <- gsub("real beta_tox;", "real beta_tox;\n    real gamma_tox;", ModeleLog)
      ModeleLog <- gsub("%sigma_coef_tox%);", "%sigma_coef_tox%);\n    gamma_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%);", ModeleLog)
      ModeleLog <- gsub("beta_tox \\* x", "beta_tox \\* x + gamma_tox \\* x_2", ModeleLog)
    }
    if (SecondCov_eff) {
      ModeleLog <- gsub("real beta_eff;", "real beta_eff;\n    real gamma_eff;", ModeleLog)
      ModeleLog <- gsub("%sigma_coef_eff%);", "%sigma_coef_eff%);\n    gamma_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%);", ModeleLog)
      ModeleLog <- gsub("beta_eff \\* x", "beta_eff \\* x + gamma_eff \\* x_2", ModeleLog)
    }
  } 
  if (PentePos_tox) {
    ModeleLog <- gsub("real beta_tox", "real<lower = 0> beta_tox", ModeleLog)
    ModeleLog <- gsub("real gamma_tox", "real<lower = 0> gamma_tox", ModeleLog)
  } 
  if (PentePos_eff) {
    ModeleLog <- gsub("real beta_eff", "real<lower = 0> beta_eff", ModeleLog)
    ModeleLog <- gsub("real gamma_eff", "real<lower = 0> gamma_eff", ModeleLog)
  } 
  ModeleLog <- gsub("%mu_inter_tox%", mu_inter_tox, ModeleLog)
  ModeleLog <- gsub("%sigma_inter_tox%", sigma_inter_tox, ModeleLog)
  ModeleLog <- gsub("%mu_coef_tox%", mu_coef_tox, ModeleLog)
  ModeleLog <- gsub("%sigma_coef_tox%", sigma_coef_tox, ModeleLog)
  ModeleLog <- gsub("%mu_inter_eff%", mu_inter_eff, ModeleLog)
  ModeleLog <- gsub("%sigma_inter_eff%", sigma_inter_eff, ModeleLog)
  ModeleLog <- gsub("%mu_coef_eff%", mu_coef_eff, ModeleLog)
  ModeleLog <- gsub("%sigma_coef_eff%", sigma_coef_eff, ModeleLog)
  return(stan_model(model_code = ModeleLog))
}

# Compile the power model (type CRM) for toxicity and efficacy
# Arguments :
## mu_coef_tox, sigma_coef_tox, mu_coef_eff, sigma_coef_eff = mean and standard deviation of the normal prior for beta coefficient
## PentePos_tox, PentePos_eff = TRUE to ensure a positive beta coefficient
# Output : Compiled power model CRM in STAN
CompilModPowCrm_efftox <- function(mu_coef_tox = 0, sigma_coef_tox = 5, PentePos_tox = FALSE,
                                   mu_coef_eff = 0, sigma_coef_eff = 5, PentePos_eff = FALSE) {
  ModelePow <- "
  data {
    int<lower=0> N;
    vector[N] x_eff;
    vector[N] x_tox;
    int<lower=0,upper=1> y_eff[N];
    int<lower=0,upper=1> y_tox[N];
  }
  parameters {
    real beta_eff;
    real beta_tox;
  }
  model {
    // Priors
    beta_tox ~ normal(%mu_coef_tox%, %sigma_coef_tox%); 
    beta_eff ~ normal(%mu_coef_eff%, %sigma_coef_eff%); 
    // Likelihood
    y_tox ~ bernoulli(exp(log(x_tox) * exp(beta_tox)));
    y_eff ~ bernoulli(exp(log(x_eff) * exp(beta_eff)));
  }
"
  if (PentePos_tox) {
    ModelePow <- gsub("real beta_tox", "real<lower = 0> beta_tox", ModelePow)
  } 
  if (PentePos_eff) {
    ModelePow <- gsub("real beta_eff", "real<lower = 0> beta_eff", ModelePow)
  } 
  ModelePow <- gsub("%mu_coef_tox%", mu_coef_tox, ModelePow)
  ModelePow <- gsub("%sigma_coef_tox%", sigma_coef_tox, ModelePow)
  ModelePow <- gsub("%mu_coef_eff%", mu_coef_eff, ModelePow)
  ModelePow <- gsub("%sigma_coef_eff%", sigma_coef_eff, ModelePow)
  return(stan_model(model_code = ModelePow))
}

# Analyse the data from 1 trial using bayesian logistic regression with CRM skeleton for toxicity and for efficacy
# Stopping rules are those of the multi-arm BOP
# Arguments :
## data = data.frame with data from the trial simulated via gen_patients_multinom function
## analyses = number of patients in each arm at analysis
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## ana_eff_cum, ana_tox_cum = number of patients in each arm at analyses for efficacy/toxicity
## phi_eff, phi_tox = threshold values for the rule to stop : futility and toxicity thresholds
## skeleton_tox = skeleton for toxicity (for example using MakeSkeleton function)
## fixed_intercept = TRUE for 1 parameter logistic regression
## A0 = intercept when fixed
## power_mod = TRUE for power model
# Output : data.frame of the trial with decisions from design
real_essai_bayeslogcrm_efftox <- function(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, 
                                          phi_eff, phi_tox, skeleton_tox, skeleton_eff, fixed_intercept, A0,
                                          power_mod = FALSE) {
  
  data$arret_eff <- data$arret_tox <- NA_integer_
  # Using the crm skeleton for doses
  data$dose <- as.numeric(gsub("^ttt(\\d+)$", "\\1", data$ttt)) 
  data$dose_eff <- skeleton_eff[data$dose]
  data$dose_tox <- skeleton_tox[data$dose]
  data$est_eff <- data$icinf_eff <- data$icsup_eff <- NA_real_
  data$est_tox <- data$icinf_tox <- data$icsup_tox <- NA_real_
  for (i in seq_len(max(data$nb_ana))) {
    Nb_pts <- analyses[i]
    seuil <- 1 - CPar * (Nb_pts / analyses[length(analyses)]) ** PPar
    if (i == 1) {
      n_eff <- data$tot_eff[data$nb_ana == i]
      n_tox <- Nb_pts - data$tot_notox[data$nb_ana == i]
      n_pts_bras <- rep(Nb_pts, length(n_tox))
      doseseff_ttt <- data$dose_eff[data$nb_ana == i]
      dosestox_ttt <- data$dose_tox[data$nb_ana == i]
    } else {
      VecNonArrets <- data$arret_eff[data$nb_ana == (i - 1)] == 0 & data$arret_tox[data$nb_ana == (i - 1)] == 0
      n_eff[VecNonArrets] <- data$tot_eff[data$nb_ana == i][VecNonArrets]
      n_tox[VecNonArrets] <- Nb_pts - data$tot_notox[data$nb_ana == i][VecNonArrets]
      n_pts_bras[VecNonArrets] <- rep(Nb_pts, sum(VecNonArrets))
      doseseff_ttt[VecNonArrets] <- data$dose_eff[data$nb_ana == i][VecNonArrets]
      dosestox_ttt[VecNonArrets] <- data$dose_tox[data$nb_ana == i][VecNonArrets]
    }
    if ((i != 1 && any(data$arret_tox[data$nb_ana == (i - 1)] == 0 & data$arret_eff[data$nb_ana == (i - 1)] == 0)) | (i == 1)) {
      if (!power_mod) {
        if (fixed_intercept) {
          DonneesEffTox <- list(N = sum(n_pts_bras),
                                y_tox = rep(rep(0:1, length(dosestox_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                                y_eff = rep(rep(0:1, length(doseseff_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_eff[x], n_eff[x])))),
                                x_tox = rep(dosestox_ttt, n_pts_bras),
                                x_eff = rep(doseseff_ttt, n_pts_bras),
                                alpha_eff = A0,
                                alpha_tox = A0)
          Modele <- "crm_fixed_efftox"
        } else {
          DonneesEffTox <- list(N = sum(n_pts_bras),
                                y_tox = rep(rep(0:1, length(dosestox_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                                y_eff = rep(rep(0:1, length(doseseff_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_eff[x], n_eff[x])))),
                                x_tox = rep(dosestox_ttt, n_pts_bras),
                                x_eff = rep(doseseff_ttt, n_pts_bras))
          Modele <- "crm_unfixed_efftox"
        }
      } else {
        DonneesEffTox <- list(N = sum(n_pts_bras),
                              y_tox = rep(rep(0:1, length(dosestox_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_tox[x], n_tox[x])))),
                              y_eff = rep(rep(0:1, length(doseseff_ttt)), unlist(lapply(seq_along(n_pts_bras), \(x) c(n_pts_bras[x] - n_eff[x], n_eff[x])))),
                              x_tox = rep(dosestox_ttt, n_pts_bras),
                              x_eff = rep(doseseff_ttt, n_pts_bras))
        Modele <- "crmpow_efftox"
      }
      SampledEffTox <- sampling(CompiledModelsEffTox[[Modele]], 
                                data = DonneesEffTox,         
                                chains = 3,             
                                warmup = 2000,          
                                iter = 10000,
                                thin = 4,
                                cores = 3,
                                control = list(stepsize = .3, adapt_delta = .95, max_treedepth = 15),
                                seed = 121221)
      DistPost <- rstan::extract(SampledEffTox)
      if (!power_mod) {
        if (fixed_intercept) {
          PPredTox <- lapply(dosestox_ttt, \(x) 1 / (1 + exp(-1 * (A0 + x * DistPost$beta_tox))))
          PPredEff <- lapply(doseseff_ttt, \(x) 1 / (1 + exp(-1 * (A0 + x * DistPost$beta_eff))))
        } else {
          PPredTox <- lapply(dosestox_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_tox + x * DistPost$beta_tox))))
          PPredEff <- lapply(doseseff_ttt, \(x) 1 / (1 + exp(-1 * (DistPost$alpha_eff + x * DistPost$beta_eff))))
        }
      } else {
        PPredTox <- lapply(dosestox_ttt, \(x) x ** exp(DistPost$beta_tox))
        PPredEff <- lapply(doseseff_ttt, \(x) x ** exp(DistPost$beta_eff))
      }
      PPTox <- vapply(PPredTox, \(distrib) {mean(distrib > phi_tox)}, double(1))
      data$est_tox[data$nb_ana == i] <- vapply(PPredTox, mean, double(1))
      data$icinf_tox[data$nb_ana == i] <- vapply(PPredTox, quantile, double(1), probs = .025)
      data$icsup_tox[data$nb_ana == i] <- vapply(PPredTox, quantile, double(1), probs = .975)
      PPEff <- vapply(PPredEff, \(distrib) {mean(distrib < phi_eff)}, double(1))
      data$est_eff[data$nb_ana == i] <- vapply(PPredEff, mean, double(1))
      data$icinf_eff[data$nb_ana == i] <- vapply(PPredEff, quantile, double(1), probs = .025)
      data$icsup_eff[data$nb_ana == i] <- vapply(PPredEff, quantile, double(1), probs = .975)
    } else {
      PPTox <- rep(1.5, length(n_tox))
      data$est_tox[data$nb_ana == i] <- data$est_tox[data$nb_ana == (i - 1)] 
      data$icinf_tox[data$nb_ana == i] <- data$icinf_tox[data$nb_ana == (i - 1)] 
      data$icsup_tox[data$nb_ana == i] <- data$icsup_tox[data$nb_ana == (i - 1)] 
      PPEff <- rep(1.5, length(n_eff))
      data$est_eff[data$nb_ana == i] <- data$est_eff[data$nb_ana == (i - 1)] 
      data$icinf_eff[data$nb_ana == i] <- data$icinf_eff[data$nb_ana == (i - 1)] 
      data$icsup_eff[data$nb_ana == i] <- data$icsup_eff[data$nb_ana == (i - 1)] 
    }
    if (i != 1) PPTox[data$arret_tox[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_tox_cum) {
      data$arret_tox[data$nb_ana == i] <- as.integer(PPTox > seuil)
    } else {
      if (i == 1) data$arret_tox[data$nb_ana == i] <- 0 else data$arret_tox[data$nb_ana == i] <- data$arret_tox[data$nb_ana == (i - 1)]
    }
    if (i != 1) PPEff[data$arret_eff[data$nb_ana == (i - 1)] == 1] <- 1.5
    if (Nb_pts %in% ana_eff_cum) {
      data$arret_eff[data$nb_ana == i] <- as.integer(PPEff > seuil)
    } else {
      if (i == 1) data$arret_eff[data$nb_ana == i] <- 0 else data$arret_eff[data$nb_ana == i] <- data$arret_eff[data$nb_ana == (i - 1)]
    }
  }
  
  return(data)
  
}


# Get operating characteristics ----

# Apply the chosen design to the simulated trials and compute the operating characteristics
# Arguments :
## ana_inter, ana_inter_tox = number of supplement patients at each analysis for efficacy and toxicity analyses
## p_n, p_a = null and alternative hypotheses with vector of length 4 (EffTox, EffNoTox, NoEffTox, NoEffNoTox)
## prior = prior multinomial distribution (if NULL, take the null hypothesis)
## phi = threshold values for the rule to stop : futility and toxicity thresholds (vector of length 2)
## mat_beta_xi = matrix of 2 rows and 4 columns for contrasts for BOP2 design with first row for efficacy and second for non toxicity
## CPar, PPar = parameters of threshold Cn of multi-arm BOP2 design
## methode = model you want to use to make analyses
## A0_tox, A0_eff = exponents for power prior
## SeuilP_tox, SeuilP_eff = pvalue threshold for test-then-pool
## a_tox, b_tox, a_eff, b_eff = hyperparameters for CBHM model
## p_mix_tox, p_mix_eff = initial guesses for proportion of EX in EXNEX model
## log_dose = TRUE to put log dose in logistic regression
## tableau_essais = data.frame of simulated trials
# Output : list of
## global operating characteristics for the design
## arm-wise operating characteristics for the design
## each trials and its analysis simulated
opcharac <- function(ana_inter,
                     ana_inter_tox = NULL,
                     p_n,
                     prior = NULL,
                     p_a,
                     phi = NULL,
                     mat_beta_xi = matrix(c(1, 1, 0, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE),
                     CPar = .6,
                     PPar = .8,
                     methode = c("bop", "bop_borrow", "bop_seq_tox", "bop_seq_efftox", "bop_power_tox", "bop_power_efftox",
                                 "bop_power_test_tox", "bop_power_test_efftox", "bop_borrow_test_tox", "bop_borrow_test_efftox",
                                 "hier_tox", "hier_efftox", "cbhm_tox", "cbhm_efftox", "exnex_tox", "exnex_efftox",
                                 "bop_log1_tox", "bop_log2_tox", "bop_log3_tox", "bop_log4_tox", "bop_log5_tox", "bop_log6_tox",
                                 "bop_log1_efftox", "bop_log2_efftox", "bop_log3_efftox", "bop_log4_efftox", "bop_log5_efftox", "bop_log6_efftox",
                                 "crm_tox", "crm_efftox", "crmpow_tox", "crmpow_efftox"),
                     A0_tox = 0, SeuilP_tox = .2, A0_eff = 0, SeuilP_eff = .2, 
                     a_tox, b_tox, a_eff, b_eff,
                     p_mix_tox, p_mix_eff, log_dose = FALSE,
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
  ana_eff_cum <- cumsum(ana_inter)
  ana_tox_cum <- if (is.null(ana_inter_tox)) cumsum(ana_inter) else cumsum(ana_inter_tox)
  analyses <- sort(union(ana_eff_cum, ana_tox_cum))
  
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
        real_essai_bop(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_borrow") {
        real_essai_bop_borrow(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_seq_tox") {
        real_essai_bop_seq_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_seq_efftox") {
        real_essai_bop_seq_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power_tox") {
        real_essai_bop_power_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power_efftox") {
        real_essai_bop_power_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_eff, A0_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power_test_tox") {
        real_essai_bop_power_test_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_tox, SeuilP_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_power_test_efftox") {
        real_essai_bop_power_test_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_eff, SeuilP_eff, A0_tox, SeuilP_tox, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_borrow_test_tox") {
        real_essai_bop_borrow_test_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_tox, Tox0, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "bop_borrow_test_efftox") {
        real_essai_bop_borrow_test_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, A0_eff, Eff0, A0_tox, Tox0, phi_eff, phi_tox, prior_eff, prior_tox)
      } else if (methode == "hier_tox") {
        real_essai_modhier_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff)
      } else if (methode == "hier_efftox") {
        real_essai_modhier_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff)
      } else if (methode == "cbhm_tox") {
        real_essai_modcbhm_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, a_tox, b_tox)
      } else if (methode == "cbhm_efftox") {
        real_essai_modcbhm_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, a_eff, b_eff, a_tox, b_tox)
      } else if (methode == "exnex_tox") {
        real_essai_modexnex_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, p_mix_tox)
      } else if (methode == "exnex_efftox") {
        real_essai_modexnex_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, p_mix_eff, p_mix_tox)
      } else if (methode %in% c("bop_log1_tox", "bop_log2_tox", "bop_log3_tox", "bop_log4_tox", "bop_log5_tox", "bop_log6_tox")) {
        ModeleAPrendre <- gsub("^bop_log(\\d)_tox$", "\\1", methode)
        tryCatch(real_essai_bayeslog_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, ModeleAPrendre, log_dose),
                 error = function(e) NULL)
      } else if (methode %in% c("bop_log1_efftox", "bop_log2_efftox", "bop_log3_efftox", "bop_log4_efftox", "bop_log5_efftox", "bop_log6_efftox")) {
        ModeleAPrendre <- gsub("^bop_log(\\d)_efftox$", "\\1", methode)
        tryCatch(real_essai_bayeslog_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, ModeleAPrendre, log_dose),
                 error = function(e) NULL)
      } else if (methode == "crm_tox") {
        real_essai_bayeslogcrm_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, SeuilP_tox, !is.na(A0_tox), A0_tox)
      } else if (methode == "crm_efftox") {
        real_essai_bayeslogcrm_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, SeuilP_tox[[1]], SeuilP_tox[[2]], !is.na(A0_tox), A0_tox)
      } else if (methode == "crmpow_tox") {
        real_essai_bayeslogcrm_tox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, prior_eff, SeuilP_tox, !is.na(A0_tox), A0_tox, TRUE)
      } else if (methode == "crmpow_efftox")
        real_essai_bayeslogcrm_efftox(data, analyses, CPar, PPar, ana_eff_cum, ana_tox_cum, phi_eff, phi_tox, SeuilP_tox[[1]], SeuilP_tox[[2]], !is.na(A0_tox), A0_tox, TRUE)
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
      nb_essais = length(unique(tableau_essais$n_simu)),
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
    essais = tableau_essais[, c("n_simu", "ttt", "nb_ana", "est_eff", "icinf_eff", "icsup_eff", "est_tox", "icinf_tox", "icsup_tox", "decision", "decision_eff", "decision_tox")]
  ))
  
}
