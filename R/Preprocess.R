##### Data Prep Functions

#' The data preparation function
#'
#' This function is used to prepare a data frame/table for the main function.
#' @param RawData The raw data for preprocessing
#' @param model_target The target virus concentration (N1 or N2)
#' @return The data frame/table for the main function
#' @examples
#' rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
#' modeldata = DataPrep(rawdata, "N1")

DataPrep = function(RawData, model_target) {
  data = RawData
  data$date = as.Date(data$date)
  data = data[labname == "NML", ]
  data = data[replicate <= 2, ]
  normdata = data[target != "PMMV" & !uqnd %in% c("UQ", "ND"), ]
  normdata = normdata[, c("Location", "date", "method", "target", "replicate", "labname", "value.raw")]
  widedata = reshape(normdata,
                     idvar = c("Location", "method", "target", "replicate", "labname"),
                     timevar = "date", direction = "wide")
  widedata = widedata[order(Location, target, replicate, method)]
  wdata = as.data.frame(matrix(NA, nrow = dim(widedata)[1] / 2, ncol = dim(widedata)[2] - 2))
  wdata[, 1:3] = unique(widedata[, c("Location", "target", "replicate")])
  names(wdata) = c("Location", "target", "replicate", names(widedata)[6:dim(widedata)[2]])
  for (i in 1:dim(wdata)[1]) {
    wdata[i, 4:dim(wdata)[2]] =
      as.numeric(apply(widedata[(2 * i - 1):(2 * i), 6:dim(widedata)[2]],
                       2, function(x) ifelse(is.na(x[2]), x[1], x[2])))
  }
  modeldata = as.data.table(melt(wdata, id.vars = c("Location", "target", "replicate"),
                                 variable.name = "date", value.name = "value.raw", na.rm = TRUE))
  modeldata$date = as.Date(gsub("value.raw.", "", modeldata$date))
  modeldata = merge(modeldata,
                    data[method == "solids",
                         c("Location", "date", "target", "replicate",
                           "dinflvol", "infltemp", "ph", "temp", "rain_t", "rain_y")],
                    by = c("Location", "target", "replicate", "date"),
                    all =  FALSE, all.x = TRUE, all.y = FALSE)
  modeldata = unique(modeldata)
  modeldata$log10.value.raw = log10(modeldata$value.raw)

  modeldata = modeldata[order(date)]
  modeldata = modeldata[target %in% model_target, ]
  return(modeldata)
}


