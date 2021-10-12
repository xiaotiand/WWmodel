
calcProba = function(xt, xt1, trend_direction) {
  temp = unlist(lapply(xt, "-", xt1))
  s = ifelse(trend_direction == 'up', 1, -1)
  return(length(which(s*temp > 0)) / length(temp))
}

#' @title Categorize probability levels.
#' @param p Vector of probabilities.
proba_categ <- function(p) {
  res = rep(NA, length(p))
  res[p < 0.4] <- 'very low'
  res[p >= 0.4 & p <0.6] <- 'low'
  res[p >= 0.6 & p <0.8] <- 'medium'
  res[p >= 0.8] <- 'strong'
  return(res)
}


#' Calculate the probability of a trend
#'
#' @param modeldata The long-format data frame/table of virus concentration
#' @param model_res The fitted model from WWmodel()
#' @param lag The time lag D in working days (See the manuscript)
#' @param trend_direction String. Direction of the trend to calculate the probability:
#' \code{up} or \code{down}.
#' @param dataframe.format Boolean. Is the result returnedas a dataframe?
#'
#' @return
#' @export
#'
proba_trend <- function(modeldata,
                        model_res,
                        lag,
                        trend_direction,
                        ID, date,
                        dataframe.format = TRUE) {

  y = modeldata[, names(modeldata) %in% c(ID, date, value),
                with = FALSE]
  Ymat = reshape(data = y,
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location, replicate)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))

  n = dim(Ymat)[1] / 2   # TODO: change this. Assume it has always 2 replicates...
  stamps = as.Date(names(Ymat))
  iprob = list()
  Yhat = model_res$Yhat

  # Calculate probabilities across all locations
  loc.name = character(n)
  for(i in 1:n) {
    # DEBUG
    # i=1
    print(paste('Calculating probability for location #',i))
    iprob[[i]] = rep(NA, dim(Yhat[[i]])[2] - lag)
    for(t in 1:length(iprob[[i]])) {
      iprob[[i]][t] = calcProba(10^(Yhat[[i]][, t + lag]),
                                10^(Yhat[[i]][, t]),
                                trend_direction)
      loc.name[i] = IDdt$Location[2*i]  # TODO: change this. it's assuming always 2 replicates for all locations
    }
  }

  res  = iprob

  if(dataframe.format){
    df = as.data.frame(iprob)
    names(df) = loc.name
    df$date = stamps[-(1:lag)]
    res = df %>%
      pivot_longer(-date, names_to = 'wwloc', values_to = 'proba')%>%
      mutate(ctg = proba_categ(proba))

    res$ctg <- factor(x = res$ctg,
                      levels = c('very low', 'low', 'medium', 'strong'))
  }

  return(res)
}

