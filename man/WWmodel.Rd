% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{WWmodel}
\alias{WWmodel}
\title{Main WWmodel function}
\usage{
WWmodel(
  modeldata,
  ID,
  date,
  value,
  covariate = NULL,
  iteration,
  burnin,
  cores = 1
)
}
\arguments{
\item{modeldata}{The long-format data frame/table of virus concentration}

\item{ID}{Names of curve IDs (used to identify a unique curve)}

\item{date}{Name of date column}

\item{value}{Name of value column}

\item{covariate}{Names of covariates (default is NULL)}

\item{iteration}{A positive integer specifying the number of iterations for each chain (including burnin).}

\item{burnin}{The number of burnin iterations}
}
\value{
fit The fitted MCMC model
}
\description{
This is the main function for fitting WWmodel.
}
\details{
This is the main function used to fit a WWmodel.
See ?WWforecast and ?FORECASTplot for how to make forecasts using a fitted WWmodel.
}
\examples{
rawdata = as.data.table(readRDS("ww-db-2021-09-10.rds"))
modeldata = DataPrep(rawdata, "N1")

ID = c("Location", "target", "replicate")
date = "date"
value = "log10.value.raw"
covariate = c("dinflvol", "temp")
model_res = WWmodel(modeldata, ID, date, value, covariate, 5000, 2500)
model_res$fit
}
