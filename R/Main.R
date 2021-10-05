##### Main Functions

#' Create an .stan model
#'
#' This is an internal function used to create a Stan model.
#' @param L0 The number of principal components kept
#' @param P The number of covariates
#' @param Lp A vector for the number of basis function for each covariate
#' @return The .stan model
#' @examples
#' cat(createStanModel(3, 5, 1:5))

createStanModel = function(L0, P, Lp = NULL) {
  scode = "data {
  int<lower=0> n;
  int<lower=0> I;
  int<lower=0> L; \n"
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode, paste("int<lower=0> L", p, ";", sep = ""), "\n")
    }
    for (p in 1:P) {
      scode = paste(scode, paste("matrix[n,L", p, "] X", p, ";", sep = ""), "\n")
    }
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("matrix[n,I] XI", l, ";", sep = ""), "\n")
  }
  scode = paste(scode,
                "vector[n] Y; \n} \n parameters { \n")
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode,
                    paste("vector[L", p, "] beta", p, ";", sep = ""), "\n",
                    paste("real<lower=0> tau", p, ";", sep = ""), "\n")
    }
  }
  scode = paste(scode, "real<lower=0> delta; \n")
  for (l in 1:L0) {
    scode = paste(scode, paste("vector[I] beta_xi", l, ";", sep = ""), "\n")
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("real<lower=0> lambda", l, ";", sep = ""), "\n")
  }
  scode = paste(scode,
                "real<lower=0> sigma; \n} \n model { \n delta ~ gamma(2, 100); \n")
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode, paste("tau", p, " ~ gamma((L", p, "+1)/2, delta/2); \n", sep = ""))
      for (l in 1:Lp[p]) {
        scode = paste(scode, paste("beta", p, "[", l, "] ~ normal(0, sqrt(tau", p, ")*sigma); \n", sep = ""))
      }
    }
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("beta_xi", l, " ~ normal(0, sqrt(lambda", l, ")); \n", sep = ""))
    scode = paste(scode, paste("lambda", l, " ~ inv_gamma(0.1, 0.1); \n", sep = ""))
  }
  scode = paste(scode, "sigma ~ inv_gamma(0.1, 0.1); \n")
  model = "Y ~ normal("
  if (P > 0) {
    for (p in 1:P) {
      model = paste(model, "X", p, "*beta", p, " + ", sep = "")
    }
  }
  for (l in 1:L0) {
    model = paste(model, "XI", l, "*beta_xi", l, sep = "")
    if (l < L0) model = paste(model, " + ", sep = "")
  }
  model = paste(model, ", sqrt(sigma)); \n", sep = "")
  scode = paste(scode, model, "}")

  return(scode)
}


#' Main WWmodel function
#'
#' This is the main function for fitting WWmodel.
#' @details
#' This is the main function used to fit a WWmodel.
#' See ?WWforecast and ?FORECASTplot for how to make forecasts using a fitted WWmodel.
#' @param modeldata The long-format data frame/table of virus concentration
#' @param ID Names of curve IDs (used to identify a unique curve)
#' @param date Name of date column
#' @param value Name of value column
#' @param covariate Names of covariates (default is NULL)
#' @param iteration A positive integer specifying the number of iterations for each chain (including burnin).
#' @param burnin The number of burnin iterations
#' @return fit The fitted MCMC model
#' @retuen Yhat A list of Yhat matrices for each site
#' @examples
#' rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
#' modeldata = DataPrep(rawdata, "N1")
#'
#' ID = c("Location", "target", "replicate")
#' date = "date"
#' value = "log10.value.raw"
#' covariate = c("dinflvol", "temp")
#' model_res = WWmodel(modeldata, ID, date, value, covariate, 5000, 2500)
#' model_res$fit

WWmodel = function(modeldata, ID, date, value, covariate = NULL, iteration, burnin) {
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  I = dim(Ymat)[1] / 2
  T = dim(Ymat)[2]
  if (!is.null(covariate)) {
    P = length(covariate)
    X = list()
    XPHI = list()
    for (p in 1:P) {
      Xp = reshape(modeldata[, names(modeldata) %in% c(ID, date, covariate[p]), with = FALSE],
                   idvar = ID, timevar = date, direction = "wide")
      Xp = Xp[order(Location, replicate)]
      Xp = Xp[, !ID, with = FALSE]
      names(Xp) = gsub(".*[.]", "", names(Xp))
      Xp = Xp[, names(Xp) %in% names(Ymat), with = FALSE]
      Xp.fit = fpca.sc(as.matrix(Xp), pve = 0.99, var = TRUE, simul = TRUE)
      Xp.fit$Yhat[2 * (1:I), ] = Xp.fit$Yhat[2 * (1:I) - 1, ]
      Xp.fit$Yhat = Xp.fit$Yhat / max(Xp.fit$Yhat)
      XPHI[[p]] = t(Xp.fit$efunctions)
      Xp = as.data.table(Xp.fit$Yhat)
      X[[p]] = Xp
    }
  }

  Ymat_mu = sweep(as.matrix(Ymat), 2, apply(Ymat, 2, function(x) mean(x, na.rm = TRUE)))
  fpca.fit = fpca.sc(as.matrix(Ymat), pve = 0.99, var = TRUE, simul = TRUE)
  fpca.fit$scores[2 * (1:I), ] = fpca.fit$scores[2 * (1:I) - 1, ]
  PHI = t(fpca.fit$efunctions)
  L0 = dim(PHI)[1]
  Lambda = sqrt(fpca.fit$evalues)
  Ymat_random = as.data.table(matrix(rep(fpca.fit$mu, 2 * I), nrow = 2 * I, byrow = TRUE))
  names(Ymat_random) = names(Ymat)
  Ymat_fix = Ymat - Ymat_random

  Y = as.numeric(t(Ymat_fix))
  missing = which(is.na(Y))
  Y = Y[-missing]
  if (!is.null(covariate)) {
    Lp = rep(NA, P)
    Xmat = list()
    for (p in 1:P) {
      Xmat[[p]] = matrix(rep(as.numeric(t(X[[p]])), dim(XPHI[[p]])[1]), ncol = dim(XPHI[[p]])[1])
      Xmat[[p]] = t(do.call("cbind", rep(list(XPHI[[p]]), I * 2))) * Xmat[[p]]
      Xmat[[p]] = Xmat[[p]][-missing, ]
      Xmat[[p]] = Xmat[[p]] / max(Xmat[[p]])
      Lp[p] = dim(Xmat[[p]])[2]
    }
    Llist = as.list(Lp)
    names(Llist) = paste("L", 1:P, sep = "")
    names(Xmat) = paste("X", 1:P, sep = "")
  }
  XImat = list()
  for (l in 1:L0) {
    XImat[[l]] = matrix(0, nrow = 2 * I * T, ncol = I)
    for (i in 1:I) XImat[[l]][((2 * i - 2) * T + 1):((2 * i) * T), i] = rep(PHI[l, ], 2)
    XImat[[l]] = XImat[[l]][-missing, ]
  }
  names(XImat) = paste("XI", 1:L0, sep = "")

  if (is.null(covariate)) {
    data = c(list(n = length(Y), I = I, L = L0),
             XImat, list(Y = Y))
    regression_model = stan_model(model_code = createStanModel(L0, 0))
  } else {
    data = c(list(n = length(Y), I = I, L = L0),
             Llist, Xmat, XImat, list(Y = Y))
    regression_model = stan_model(model_code = createStanModel(L0, P, Lp))
  }
  fit = rstan::sampling(regression_model, data = data, chains = 2, iter = iteration, refresh = 0)
  list_of_draws = extract(fit)

  after = iteration - burnin + 1
  Yhat = list()
  for (i in 1:I) {
    Yhat[[i]] = matrix(rep(unlist(Ymat_random[2 * i, ]), after), nrow = after, byrow = TRUE)
    for (l in 1:L0) {
      Yhat[[i]] = Yhat[[i]] +
        (list_of_draws[[paste("beta_xi", l, sep = "")]])[burnin:iteration, i] %*% t(PHI[l, ])
    }
    if (!is.null(covariate)) {
      for (p in 1:P) {
        Yhat[[i]] = Yhat[[i]] +
          (list_of_draws[[paste("beta", p, sep = "")]])[burnin:iteration, ] %*% XPHI[[p]] * matrix(rep(unlist(X[[1]][2 * i, ]), after), nrow = after, byrow = TRUE)
      }
    }
  }

  return(list(fit = fit, Yhat = Yhat))
}

