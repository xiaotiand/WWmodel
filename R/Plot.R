##### Plot/diagnostic Functions


plot_fitted_trend <- function(modeldata, model_res, ID, date) {

  # --- reshape input data
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))

  dat = data.frame(t(Ymat))
  names(dat) <- IDdt$Location  # <<== NOT SURE ABOUT THIS!!!!

  thedates = ymd(names(Ymat))
  dat$date = thedates
  nloc = dim(Ymat)[1] / 2
  loc = IDdt$Location[2 * (1:nloc)]

  # --- Fitted data

  Yhat = model_res$Yhat

  # --- Merging observation anf fitted trends

  datlong = dat %>%
    pivot_longer(-date, names_to = 'Location', values_to = 'value.obs')

  extract_fit_ss <- function(i, Yhat) {
    a = data.frame(t(Yhat[[i]]))
    a$date <- thedates

    a.long = a %>%
      pivot_longer(-date)

    ci = 0.95
    a.ss = a.long %>%
      group_by(date) %>%
      summarise(m = mean(value),
                qlo = quantile(value, probs = 0.5 - ci/2),
                qhi = quantile(value, probs = 0.5 + ci/2)) %>%
      mutate(Location = loc[i])
    return(a.ss)
  }

  tmp = lapply(1:nloc, extract_fit_ss, Yhat)
  df = do.call('rbind', tmp)

  dj = left_join(df, datlong, by = c('Location', 'date'))

  col.fit = 'green4'
  g = dj %>%
    ggplot(aes(x=date)) +
    geom_ribbon(aes(ymin=qlo, ymax=qhi), alpha = 0.23,
                fill = col.fit, color = col.fit, size=0.2)+
    geom_line(aes(y=m), size = 1,  color  = col.fit) +
    geom_point(aes(y = value.obs), alpha = 0.5, shape = 16) +
    facet_wrap(~Location, scales = 'free_y') +
    theme(panel.grid.minor = element_blank())+
    labs(
      title = 'Fitted trends',
      x='', y='log concentration'
    )
  return(g)
}


#' Create a plot of the sampled curves
#'
#' This function is used to create a plot of simultaneous confidence bands.
#' @param modeldata The long-format data frame/table of virus concentration
#' @param model_res The fitted model from WWmodel()
#' @param ncol The number of columns in the plot
#' @return The plot of simultaneous confidence bands.
#' @examples
#' SAMPLEplot(modeldata, model_res, 4)

SAMPLEplot = function(modeldata, model_res, ncol) {
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  I = dim(Ymat)[1] / 2

  nrow = ceiling(I / ncol)
  par(mfrow = c(nrow, ncol))
  Yhat = model_res$Yhat
  ylim = max(unlist(lapply(Yhat, max)))
  for (i in 1:I) {
    plot(Yhat[[i]][1, ], type = "l", col = "#00000010", xlab = paste(IDdt[2 * i, 1]),
         ylab = "Y", ylim = c(0, ylim))
    for(index in 2:dim(Yhat[[i]])[1]) {
      lines(Yhat[[i]][index, ], type = "l", col = "#00000010")
    }
    points(1:length(Ymat[2 * i - 1, ]), Ymat[2 * i - 1, ], col = "red", pch = 19, cex = 0.5)
    points(1:length(Ymat[2 * i, ]), Ymat[2 * i, ], col = "red", pch = 19, cex = 0.5)
  }
}


#' Create a plot of the sampled curves and forecasts
#'
#' This function is used to create a plot of simultaneous confidence bands.
#' @param modeldata The long-format data frame/table of virus concentration
#' @param targetdata The ``future'' observations to be forecasted (the default is NULL)
#' @param model_res The fitted model from WWmodel()
#' @param ncol The number of columns in the plot
#' @param legend If a lengend should be included (the default is TRUE)
#' @return The plot of simultaneous confidence bands.
#' @examples
#' rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
#' modeldata = DataPrep(rawdata, "N1")
#' targetdata = modeldata[date < as.Date("2021-06-15") & date > as.Date("2021-05-15")]
#' modeldata = modeldata[date <= as.Date("2021-05-15")]
#'
#' ID = c("Location", "target", "replicate")
#' date = "date"
#' value = "log10.value.raw"
#' covariate = c("dinflvol", "temp")
#' model_res = WWmodel(modeldata, ID, date, value, covariate, 5000, 2500)
#'
#' forecast_res = WWforecast(20, modeldata, model_res, ID, date, value, covariate, 5000, 2500)
#' FORECASTplot(modeldata, model_res, forecast_res, 4, targetdata)

FORECASTplot = function(modeldata, model_res, forecast_res, ncol, targetdata = NULL, legend = TRUE) {
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  I = dim(Ymat)[1] / 2
  if (!is.null(targetdata)) {
    Ymat2 = reshape(targetdata[, names(targetdata) %in% c(ID, date, value), with = FALSE],
                    idvar = ID, timevar = date, direction = "wide")
    Ymat2 = Ymat2[order(Location, replicate)]
    Ymat2 = Ymat2[, !ID, with = FALSE]
    names(Ymat2) = gsub(".*[.]", "", names(Ymat2))
  }

  nrow = ceiling(I / ncol)
  par(mfrow = c(nrow, ncol))
  Yhat = model_res$Yhat
  ylim = max(unlist(lapply(Yhat, max)))
  par(mfrow = c(nrow, ncol))
  for (i in 1:I) {
    l1 = length(Ymat[2 * i, ])
    l2 = dim(forecast_res[[i]])[2]
    plot(1:l1, Yhat[[i]][1, ], type = "l", col = "#00000010", xlab = paste(IDdt[2 * i, 1]),
         ylab = "Y", xlim = c(0, l1 + l2), ylim = c(0, ylim))
    lines((l1 + 1):(l1 + l2), forecast_res[[i]][1, ], type = "l", col = "#90503F10")
    for(index in 2:dim(Yhat[[i]])[1]) {
      lines(1:l1, Yhat[[i]][index, ], type = "l", col = "#00000010")
      lines((l1 + 1):(l1 + l2), forecast_res[[i]][index, ], type = "l", col = "#90503F10")
    }
    abline(v = dim(Yhat[[i]])[2], col = "red", lty = 1.5)
    points(1:l1, Ymat[2 * i - 1, ], col = "red", pch = 19, cex = 0.5)
    points(1:l1, Ymat[2 * i, ], col = "red", pch = 19, cex = 0.5)
    if (!is.null(targetdata)) {
      l2 = dim(Ymat2)[2]
      points((l1 + 1):(l1 + l2), Ymat2[2 * i - 1, ], col = "blue", pch = 19, cex = 0.5)
      points((l1 + 1):(l1 + l2), Ymat2[2 * i, ], col = "blue", pch = 19, cex = 0.5)
    }
  }
  if (legend) {
    plot(1:5, rep(4, 5), type = "l", lwd = 2, xaxt = 'n', yaxt = 'n',
         xlab = "", ylab = "", xlim = c(0, 12), ylim = c(0, 5))
    text(8, 4, "Sampled Curves")
    lines(1:5, rep(3, 5), type = "l", lwd = 2, col = "brown")
    text(8, 3, "Forecasted Curves")
    points(3, 2, col = "red", pch = 16, cex = 2)
    text(8, 2, "Observed Value (Training)")
    points(3, 1, col = "blue", pch = 16, cex = 2)
    text(8, 1, "Observed Value (Testing)")
  }
}


#' Calculate the probability of increase in virus concentration
#'
#' This function is used to calculate the probability of increase in virus concentration. A plot is returned.
#' @param modeldata The long-format data frame/table of virus concentration
#' @param model_res The fitted model from WWmodel()
#' @param lag The time lag D in working days (See the manuscript)
#' @param ncol The number of columns in the plot
#' @return The plot of probabilities.
#' @examples
#' PROBplot(modeldata, model_res, 5, 4)

PROBplot = function(modeldata, model_res, lag, ncol) {
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  I = dim(Ymat)[1] / 2
  stamps = as.Date(names(Ymat))
  iprob = list()
  Yhat = model_res$Yhat
  prob = function(xt, xt1) {
    temp = unlist(lapply(xt, "-", xt1))
    return(length(which(temp > 0)) / length(temp))
  }
  for(i in 1:I) {
    iprob[[i]] = rep(NA, dim(Yhat[[i]])[2] - lag)
    for(t in 1:length(iprob[[i]])) {
      iprob[[i]][t] = prob(10^(Yhat[[i]][, t + lag]), 10^(Yhat[[i]][, t]))
    }
  }

  nrow = ceiling(I / ncol)
  par(mfrow = c(nrow, ncol))
  for(i in 1:I) {
    plot(stamps[-(1:lag)], iprob[[i]], type = "l", col = "black", xlab = paste(IDdt[2 * i, 1]),
         ylab = "Probability", ylim = c(0, 1))
    abline(h = 0.5, col = "red")
  }
}

