#' @rawNamespace useDynLib(tsforeign, .registration = TRUE)
#' @keywords internal
#' @import methods
#' @import tsmethods
#' @import data.table
#' @importFrom utils head tail
#' @importFrom tsaux do.call.fast box_cox check_xreg check_newxreg mape mslre mase bias mis future_dates sampling_frequency fourier_series crps tstransform
#' @importFrom stats median na.omit fitted coef quantile residuals predict as.formula rpois runif rnorm acf optim qnorm sd simulate arima.sim diffinv ts tsp tsp<- na.contiguous var alias deviance df.residual formula hat lm model.matrix pchisq pf vcov weights logLik pnorm frequency window dnorm
#' @importFrom zoo index as.zoo zoo coredata na.locf na.fill
#' @importFrom grDevices n2mfrow
#' @importFrom dlm dlmModARMA dlmModTrig dlmModSeas dlmModPoly dlm dlmFilter dlmSmooth dlmForecast
#' @importFrom graphics layout  lines  par  plot grid  legend  mtext
#' @importFrom xts xts as.xts is.xts endpoints
#' @importFrom forecast auto.arima forecast getResponse BoxCox InvBoxCox
#' @importFrom mvtnorm rmvnorm
#' @importFrom tools toTitleCase
#' @importFrom coda as.mcmc
#' @importFrom bsts AddLocalLinearTrend AddLocalLevel AddSeasonal bsts SuggestBurn AddTrig AddSemilocalLinearTrend AddAutoAr bsts.prediction.errors BstsOptions
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom Rcpp evalCpp loadModule
"_PACKAGE"
