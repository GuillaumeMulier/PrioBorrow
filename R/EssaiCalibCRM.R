
library(progressr)
handlers(handler_progress(
  format = ":message ==> :current/:total (:percent) :bar",
  width = 80, complete = "=", incomplete = " "
))

# CRM functions 

BinomialLik <- function(CoefBeta, Outcome, Dose, DoseLogit) {
  Probas <- 1 / (1 + exp(-A0 - CoefBeta * DoseLogit[Dose]))
  return(prod(Probas ** Outcome * (1 - Probas) ** (1 - Outcome)))
}
VecBinomialLik <- Vectorize(BinomialLik, vectorize.args = "CoefBeta")
CRM <- function(CoefBeta, Outcome, Dose, DoseLogit){
  VecBinomialLik(CoefBeta, Outcome, Dose, DoseLogit) * dnorm(CoefBeta, mean = 0, sd = sqrt(1.34))
}
CRMBeta <- function(CoefBeta, Outcome, Dose, DoseLogit){
  CoefBeta * VecBinomialLik(CoefBeta, Outcome, Dose, DoseLogit) * dnorm(CoefBeta, mean = 0, sd = sqrt(1.34))
}
DoseSuivante <- function(y, d, DoseLogit, ProbaTox) {
  Denom <- integrate(CRM, -Inf, Inf, y, d, DoseLogit)$value
  Num <- integrate(CRMBeta, -Inf, Inf, y, d, DoseLogit)$value
  BetaUpdated <- Num / Denom
  F_hat <- 1 / (1 + exp(-A0 - BetaUpdated * DoseLogit))
  next_dose <- pmin(which.min(abs(F_hat - ProbaTox)), d[length(d)]+1)
  return(list(next_dose, F_hat, BetaUpdated))
}
RealEssai <- function(Iter, TailleCohorte, DoseIni, NbCohorte, DoseLogit, ProbaTox, ScenarTox, MTD, Mess) {
  nextdose <- DoseIni
  y <- numeric(NbCohorte * TailleCohorte)
  d <- numeric(NbCohorte * TailleCohorte)
  Betas <- numeric(NbCohorte)
  for (l in seq_len(NbCohorte)) {
    d[seq((l - 1) * TailleCohorte + 1, l * TailleCohorte)] <- nextdose
    y[seq((l - 1) * TailleCohorte + 1, l * TailleCohorte)] <- rbinom(TailleCohorte, 1, ScenarTox[nextdose])
    res <- DoseSuivante(y[seq_len(l * TailleCohorte)], d[seq_len(l * TailleCohorte)], DoseLogit, ProbaTox)
    nextdose <- res[[1]]
    Betas[l] <- res[[3]]
  }
  if (Iter %% 100 == 0) {print(Mess)}
  return(nextdose == MTD)
}
PCSDelta <- function(Delta, ProbaTox, MTD, NDoses, A0, B0, 
                     Psi, NbSimu , TailleCohorte, DoseIni, NbCohorte) {
  DoseLogit <- numeric(NDoses)
  DoseLogit[MTD] <- (log(ProbaTox / (1 - ProbaTox)) - A0) / B0
  if (MTD > 1) {
    for(i in rev(seq_len(MTD)[-1])) {
      DoseLogit[i - 1] <- (log((ProbaTox - Delta) / (1 - (ProbaTox - Delta))) - A0) * DoseLogit[i] / (log((ProbaTox + Delta) / (1 - (ProbaTox + Delta))) - A0)
    }
  }
  if (MTD < NDoses) {
    for (i in seq(MTD, NDoses - 1)) {
      DoseLogit[i + 1] <- (log((ProbaTox + Delta) / (1 - (ProbaTox + Delta))) - A0) * DoseLogit[i] / (log((ProbaTox - Delta) / (1 - (ProbaTox - Delta))) - A0)
    }
  }
  VecProba <- 1 / (1 + exp(-A0 - B0 * DoseLogit))
  
  # Calibration set
  PU <- ProbaTox * Psi / (1 - ProbaTox * (1 - Psi))
  PL <- ProbaTox / (Psi + ProbaTox * (1 - Psi))
  MatScenar <- do.call("rbind", lapply(seq_len(NDoses), \(x) rep(c(PL, ProbaTox, PU), c(sum(seq_len(NDoses) < x), sum(seq_len(NDoses) == x), sum(seq_len(NDoses) > x)))))
  
  # PCS
  Res <- numeric(nrow(MatScenar))
  for (i in seq_len(nrow(MatScenar))) {
    Res[i] <- mean(vapply(seq_len(NbSimu), \(iter) {
      RealEssai(iter, TailleCohorte = TailleCohorte, DoseIni = DoseIni, NbCohorte = NbCohorte, DoseLogit = DoseLogit, ProbaTox = ProbaTox,
                ScenarTox = MatScenar[x, ], MTD = MTD, Mess = paste0("Delta = ", Delta, ", scénario n°", i, ", itération ", iter))
    }, logical(1)))
  }
  return(mean(Res))
}


set.seed(121221)
ResDelta <- data.frame(del = seq(.01, .1, .01))
ResDelta$pcs <- do.call("c", lapply(ResDelta$del, \(Delta) {
  PCSDelta(Delta = Delta, ProbaTox = .1, MTD = 2, NDoses = 5, A0 = 3, B0 = 1, 
           Psi = 2, NbSimu = 1e+3, TailleCohorte = 3, DoseIni = 1, NbCohorte = 10)
}))

crmsens(c(.05, .2, .3, .5), .2, "logistic")

prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
target <- 0.2
level <- c(3, 4, 4, 3, 3, 4, 3, 2, 2, 2)
y <- c(0, 0, 1, 0, 0, 1, 1, 0, 0, 0)
foo <- crm(prior, target, y, level, model = "logistic")

logitprior <- (log(prior / (1 - prior)) - 3) / 1
DoseSuivante(y, level, logitprior, target)



PI <- c(0.10, 0.20, 0.40, 0.50, 0.60, 0.65)
prior <- c(0.05, 0.10, 0.20, 0.35, 0.50, 0.70)
target <- 0.2
x0 <- 1
# Generate 10 replicates of CRM trial with 24 subjects
foo10 <- crmsim(PI, prior, target, 24, 3, nsim=10, mcohort=3, model = "logistic")
foo10


DoseSuivante(y, d, DoseLogit, ProbaTox)
y <- c(0, 0, 0)
d <- c(1, 1, 1)
Denom <- integrate(CRM, -Inf, Inf, y, d, DoseLogit)$value
Num <- integrate(CRMBeta, -Inf, Inf, y, d, DoseLogit)$value
BetaUpdated <- Num / Denom
F_hat <- 1 / (1 + exp(-A0 - BetaUpdated * DoseLogit))



BinomialLik(1, c(1, 1, 0), c(1, 1, 1), VecProba)
BinomialLik(15, c(1, 1, 0), c(1, 1, 1), VecProba)
VecBinomialLik(c(1, 15), c(1, 1, 0), c(1, 1, 1), VecProba)


library(dfcrm)
getprior(.01, .2, 2, 4, "logistic", 3)

halfwidth <- .01
target <- .2
nu <- 2
nlevel <- 4
intcpt <- 3
dosescaled <- prior <- rep(NA, nlevel)
b <- rep(NA, nlevel + 1)
b[1] <- -Inf
b[(nlevel + 1)] <- Inf
dosescaled[nu] <- log(target/(1 - target)) - intcpt
for (k in nu:2) {
  b[k] <- log((log((target + halfwidth)/(1 - target - 
                                           halfwidth)) - intcpt)/dosescaled[k])
  if (nu > 1) {
    dosescaled[k - 1] <- (log((target - halfwidth)/(1 - 
                                                      target + halfwidth)) - intcpt)/exp(b[k])
  }
}
if (nu < nlevel) {
  for (k in nu:(nlevel - 1)) {
    b[k + 1] <- log((log((target - halfwidth)/(1 - 
                                                 target + halfwidth)) - intcpt)/dosescaled[k])
    dosescaled[k + 1] <- (log((target + halfwidth)/(1 - 
                                                      target - halfwidth)) - intcpt)/exp(b[k + 1])
  }
}
val <- {
  1 + exp(-intcpt - dosescaled)
}^{
  -1
}