##### Data Prep Functions

#' The data preparation function
#'
#' This function is used to prepare a data frame/table for the main function.
#' @param RawData The raw data for preprocessing
#' @param model_target The target virus concentration (N1 or N2)
#' @param lab A list of lab sites. The default is "all". The final results will be the intersection of labs and cities selected.
#' @param city A list of cities (Edmonton, Halifax, Montreal, Toronto, and Vancouver). The default is "all".
#' The final results will be the intersection of labs and cities selected.
#' @return The data frame/table for the main function
#' @examples
#' rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
#' modeldata = DataPrep(rawdata, "N1")
#' modeldata = DataPrep(rawdata, "N1", city = c("Toronto", "Edmonton"))

DataPrep = function(RawData, model_target, lab = "all", city = "all") {
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
  if (lab != "all") {
    modeldata = modeldata[Location %in% lab, ]
  }
  if (city != "all") {
    sites = c()
    if ("Edmonton" %in% city) sites = c(sites, "EGB")
    if ("Halifax" %in% city) sites = c(sites, "HDA", "HHA", "HMC")
    if ("Montreal" %in% city) sites = c(sites, "MMN", "MMS")
    if ("Toronto" %in% city) sites = c(sites, "TAB", "THC", "THU", "TNT")
    if ("Vancouver" %in% city) sites = c(sites, "VAI", "VII", "VLG", "VLI", "VNL")
    modeldata = modeldata[Location %in% sites, ]
  }
  return(modeldata)
}


