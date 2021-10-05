##### Forecasting Functions

#' Main forecasting function
#'
#' This is the main forecasting function for a fitted WWmodel.
#' @details
#' This function is used to make forecasts using a fitted WWmodel.
#' See ?FORECASTplot for how to plot the forecasts.
#' @param h.ahead The number of steps ahead for which prediction is required
#' @param modeldata The long-format data frame/table of virus concentration
#' @param model_res The model result from WWmodel()
#' @param ID Names of curve IDs (used to identify a unique curve)
#' @param date Name of date column
#' @param value Name of value column
#' @param covariate Names of covariates (default is NULL)
#' @param iteration A positive integer specifying the number of iterations for each chain (including burnin).
#' @param burnin The number of burnin iterations
#' @return Ypred The forecast of observed virus concentration
#' @examples
#' rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
#' modeldata = DataPrep(rawdata, "N1")
#'
#' ID = c("Location", "target", "replicate")
#' date = "date"
#' value = "log10.value.raw"
#' covariate = c("dinflvol", "temp")
#' model_res = WWmodel(modeldata, ID, date, value, covariate, 5000, 2500)
#'
#' forecast_res = WWforecast(5, modeldata, model_res, ID, date, value, covariate, 5000, 2500)

WWforecast = function(h.ahead, modeldata, model_res, ID, date, value, covariate = NULL, iteration, burnin) {
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
      X[[p]] = Xp / max(Xp)
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

  Yhat_pred = as.numeric(forecast::forecast(auto.arima(fpca.fit$mu, d = 2), h = h.ahead)$mean)
  PHIpred = t(apply(PHI, 1, function(x) as.numeric(forecast::forecast(auto.arima(x, d = 2), h = h.ahead)$mean)))
  Xpred = list()
  XPHIpred = list()
  if (!is.null(covariate)) {
    for (p in 1:P) {
      Xpred[[p]] = t(apply(X[[p]], 1, function(x) as.numeric(rep(x[length(x)], h.ahead))))
      XPHIpred[[p]] = t(apply(XPHI[[p]], 1,
                              function(x) as.numeric(forecast::forecast(auto.arima(x, d = 2), h = h.ahead)$mean)))
    }
  }
  list_of_draws = extract(model_res$fit)
  after = iteration - burnin + 1
  Ypred = list()
  for (i in 1:I) {
    Ypred[[i]] = matrix(rep(Yhat_pred, after), nrow = after, byrow = TRUE)
    for (l in 1:L0) {
      Ypred[[i]] = Ypred[[i]] +
        (list_of_draws[[paste("beta_xi", l, sep = "")]])[burnin:iteration, i] %*% t(PHIpred[l, ])
    }
    if (!is.null(covariate)) {
      for (p in 1:P) {
        Ypred[[i]] = Ypred[[i]] +
          (list_of_draws[[paste("beta", p, sep = "")]])[burnin:iteration, ] %*% XPHIpred[[p]] * matrix(rep(unlist(Xpred[[1]][2 * i, ]), after), nrow = after, byrow = TRUE)
      }
    }
  }

  return(Ypred)
}

