#' Auto ARIMA specification
#'
#' @description Sets up an automatic ARIMA specification ready for estimation.
#' @param y an xts vector.
#' @param xreg an xts matrix of external regressors.
#' @param frequency frequency of y (if using a seasonal model).
#' @param seasonal whether to include a seasonal component.
#' @param seasonal_type type of seasonality (regular or trigonometric).
#' @param seasonal_harmonics number of harmonics to include in the seasonal 
#' component when seasonal_type is trigonometric.
#' @param transformation a valid transformation for y from the \dQuote{tstransform} 
#' function in the \dQuote{tsaux} package (currently box-cox or logit are available).
#' @param lambda the Box Cox lambda. If not NULL, then either a numeric value or NA 
#' denoting automatic calculation.
#' @param lower lower bound for the transformation.
#' @param upper upper bound for the transformation.
#' @param ... additional terms passed to the \code{\link{auto.arima}} function.
#' @return An object of class \dQuote{arima.spec}.
#' @note This is a wrapper to the auto.arima function from the forecast package.
#' @aliases arima_modelspec
#' @rdname arima_modelspec
#' @export
#'
#'
#'
#'
arima_modelspec = function(y, xreg = NULL, frequency = NULL, seasonal = FALSE, seasonal_type = "regular", seasonal_harmonics = NULL, 
                           transformation = "box-cox", lambda = NULL, lower = 0, upper = 1.5, ...)
{
  # 1. Check y
  if (!is.xts(y)) {
    stop("y must be an xts object")
  }
  if (any(is.na(y))) {
    stop("\nNAs found in y...not allowed.")
  }
  # 2. Check regressors
  xreg <- check_xreg(xreg, index(y))
  call <- list(frequency = frequency, seasonal = seasonal, seasonal_type = seasonal_type, transformation = transformation, 
               lambda = lambda, lower = lower, upper = upper, seasonal_harmonics = seasonal_harmonics)
  # 3. Check transformation
  y_orig <- y
  if (!is.null(lambda)) {
    if (!is.na(lambda) & lambda == 1) lambda <- NULL
  }
  transformation <- match.arg(transformation[1], c("box-cox","logit"))
  if (transformation == "logit") {
      lambda <- 1
      transform <- tstransform(method = transformation, lower = lower, upper = upper)
      transform$lambda <- 1
      transform$name <- "logit"
      transform$include_lambda <- FALSE
      transform$lower <- lower
      transform$upper <- upper
      y <- transform$transform(y)
      y <- as.numeric(y)
  } else {
      if (is.null(lambda)) {
          transform <- NULL
      } else if (is.na(lambda)) {
          include_lambda <- TRUE
          transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
          y <- transform$transform(y = y, frequency = frequency[1])
          transform$lambda <- attr(y, "lambda")
          transform$include_lambda <- TRUE
          lambda <- transform$lambda
          transform$lower <- lower
          transform$upper <- upper
          transform$name <- "box-cox"
      } else {
          include_lambda <- FALSE
          transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
          y <- transform$transform(y = y, frequency = frequency[1])
          transform$include_lambda <- FALSE
          transform$lambda <- lambda
          transform$lower <- lower
          transform$upper <- upper
          transform$name <- "box-cox"
      }
      y <- as.numeric(y)
  }
  # 5. seasonal part
  if (seasonal) {
    if (!seasonal_type == "regular" | max(frequency) > 350 | length(frequency) > 1) {
      seasonal_type <- "trigonometric"
      if (is.null(seasonal_harmonics)) {
        seasonal_harmonics <- floor(0.25*frequency)
      } else {
        if (length(seasonal_harmonics) != length(frequency)) {
          warnings("\nlength of fourier harmonics (seasonal_harmonics) not equal to length of frequency. Defaulting to 0.25*frequency.")
          seasonal_harmonics <- floor(0.25*frequency)
        }
        # Should probably check that seasonal_harmonics < 1/2xfrequency
      }
      # Create fourier terms
      seasonal_fourier <- do.call.fast(cbind, lapply(1:length(frequency), function(i){
        fourier_series(index(y), period = frequency[i], K = seasonal_harmonics[i])
      }))
    } else{
      seasonal_fourier <- NULL
      seasonal_type <- "regular"
    }
  } else{
    seasonal_fourier <- NULL
  }
  spec <- list()
  spec$target$y <- as.numeric(y)
  spec$target$y_orig <- as.numeric(y_orig)
  spec$target$index <- index(y_orig)
  spec$target$frequency <- frequency
  spec$target$sampling <- sampling_frequency(index(y_orig))
  spec$transform <- transform
  spec$arima_args <- list(...)
  spec$seasonal <- list(seasonal = seasonal, seasonal_type = seasonal_type, seasonal_fourier = seasonal_fourier, 
                        frequency = frequency, seasonal_harmonics = seasonal_harmonics)
  if (!is.null(xreg)) {
    spec$xreg$xreg <- coredata(xreg)
  } else{
    spec$xreg <- NULL
  }
  spec$call <- call
  class(spec) <- c("arima.spec","tsmodel.spec")
  return(spec)
}

#' Model Estimation
#'
#' @description Estimates a model given a specification object using
#' maximum likelihood.
#' @param object an object of class \dQuote{arima.spec} or \dQuote{bsts.spec}.
#' @param n_iter mcmc draws for the bsts model.
#' @param timeout.seconds timeout passed to the BstsOptions function.
#' @param bma.method Bayesian Model Averaging method in the presence of 
#' regressors, passed to the BstsOptions function.
#' @param trace whether to show the MCMC iterations.
#' @param ... for the bsts model, additional arguments passed to the underlying 
#' estimation functions in the BSTS package.
#' @return An object of class \dQuote{arima.estimate} or \dQuote{bsts.estimate}.
#' @aliases estimate
#' @method estimate arima.spec
#' @rdname estimate
#' @export
#'
#'
estimate.arima.spec <- function(object, ...)
{
  if (object$seasonal$seasonal) {
    if (object$seasonal$seasonal_type == "trigonometric") {
      yuse <- object$target$y
      if (!is.null(object$xreg$xreg)) {
        xreg <- cbind(object$xreg$xreg, object$seasonal$seasonal_fourier)
      } else{
        xreg <- object$seasonal$seasonal_fourier
      }
      sflag <- FALSE
    } else{
      yuse = ts(object$target$y, frequency = object$target$frequency)
      xreg <- object$xreg$xreg
      sflag <- TRUE
    }
  } else{
    yuse = ts(object$target$y, frequency = 1)
    xreg <- object$xreg$xreg
    sflag <- FALSE
  }
  args <- object$arima_args
  args$xreg <- xreg
  args$seasonal <- sflag
  args$y <- yuse
  model <- do.call(auto.arima, args = args)
  # save the order of the regressors
  if (!is.null(xreg)) {
    object$xreg$xreg_names <- colnames(xreg)
  }
  out <- list(model = model, spec = object)
  class(out) <- c("arima.estimate","tsmodel.estimate")
  return(out)
}

#' Model Fitted Values
#'
#' @description Extract the fitted values from an estimated model.
#' @param object an object of class \dQuote{arima.estimate} or 
#' \dQuote{bsts.estimate}.
#' @param raw if a transformation was used, whether to return the transformed 
#' fitted values (raw).
#' @param distribution whether to return the full distribution for the bsts 
#' model.
#' @param invdiff if a model in differences was estimated, whether to inverse 
#' the differencing when calculating the fitted values.
#' @param type for the bsts model, the filtered are the one step ahead forecast 
#' values, whilst the smoothed represent those values using the whole information 
#' set.
#' @param ... not currently used.
#' @return For the distribution choice, an object of class 
#' \dQuote{tsmodel.distribution}, else an xts vector.
#' @aliases fitted
#' @method fitted arima.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.arima.estimate = function(object, raw = FALSE, ...)
{
  transform <- object$spec$transform
  ft <- xts(fitted(object$model), object$spec$target$index)
  if (!raw & !is.null(transform)) {
    ft <- transform$inverse(ft, transform$lambda)
  }
  return(ft)
}

#' Model Residuals
#'
#' @description Extract the residual values from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param distribution whether to return the full distribution for the bsts 
#' model.
#' @param invdiff if a model in differences was estimated, whether to inverse 
#' the differencing when calculating the residual values.
#' @param standardize whether to scale the errors by the model standard 
#' deviation.
#' @param raw raw residuals are the model based values in transformed space 
#' (when either Box Cox or Logistic have been used as transformations).
#' @param type for the bsts model, the filtered are the one step ahead forecast 
#' errors, whilst the smoothed represent those values using the whole 
#' information set.
#' @param ... not currently used.
#' @return For the distribution choice, an object of class 
#' \dQuote{tsmodel.distribution}, else an xts vector.
#' @aliases residuals
#' @method residuals arima.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.arima.estimate = function(object, raw = FALSE, ...)
{
  if (raw) {
    e <- as.numeric(object$spec$target$y) - as.numeric(fitted(object, raw = raw))
  } else {
    e <- as.numeric(object$spec$target$y_orig) - as.numeric(fitted(object$model))
  }
  e <- xts(e, object$spec$target$index)
  return(e)
}

#' Model Estimation Summary
#'
#' @description Summary method.
#' @param object an object of class \dQuote{arima.spec} or \dQuote{bsts.spec}.
#' @param quantiles vector of quantiles to evaluate for each variable.
#' @param ... not currently used.
#' @return A printout of the parameter summary, model type and some model metrics. 
#' For the bsts model, the posterior distribution of parameters is fist converted 
#' into an \link{mcmc} object (coda package).
#' @aliases summary
#' @method summary arima.estimate
#' @rdname summary
#' @export
#'
#'
summary.arima.estimate <- function(object, ...)
{
  modelx <- arima_string(object$model)
  # data.table
  cat(modelx,"\n")
  cat(paste0(rep("-",nchar(modelx) + 1),collapse = ""),"\n")
  coef <- coef(object$model)
  standard_errors <- sqrt(abs(diag(object$model$var.coef)))
  t_values <- coef/standard_errors
  p_values <- 2*(1 - pnorm(abs(t_values)))
  out <- data.frame("Estimate" = as.numeric(coef), "Std.Error" = standard_errors, "t value" = t_values, "Pr(>|t|)" = p_values, row.names = names(coef), check.names = FALSE)
  print(out, digits = 4)
  cat("\nsigma :", round(sqrt(object$model$sigma2),3))
  cat("\n\n")
  mtrcs <- tsmetrics(object)
  print(data.frame(AIC = as.numeric(sprintf(mtrcs$AIC,fmt = "%.2f")),BIC = as.numeric(sprintf(mtrcs$BIC,fmt = "%.2f")),
                   AICc = as.numeric(sprintf(mtrcs$AICc,fmt = "%.2f"))), row.names = FALSE, right = T)
  cat("\n")
  print(data.frame(MAPE = as.numeric(sprintf(mtrcs$MAPE,fmt = "%.4f")), MASE = as.numeric(sprintf(mtrcs$MASE,fmt = "%.4f")),
                   MSLRE = as.numeric(sprintf(mtrcs$MSLRE,fmt = "%.4f")), BIAS = as.numeric(sprintf(mtrcs$BIAS,fmt = "%.4f"))), 
        row.names = FALSE, right = T)
  return(invisible(out))
}

#' Performance Metrics
#'
#' @description Performance metrics from an estimated or predicted model.
#' @param object An object of class \dQuote{arima.estimate}, 
#' \dQuote{arima.predict}, \dQuote{bsts.estimate} or \dQuote{bsts.predict}.
#' @param actual the actual data matched to the dates of the forecasts.
#' @param alpha the coverage level for distributional forecast metrics.
#' @param ... Optional arguments passed to the MASE function which includes 
#' \dQuote{frequency}, else will be read from the object (act as an override for 
#' instance when using fourier seasonality in xreg for arima).
#' @aliases tsmetrics
#' @method tsmetrics arima.estimate
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.arima.estimate <- function(object, ...)
{
  fx <- fitted(object)
  ac <- xts(object$spec$target$y_orig, object$spec$target$index)
  acfx <- na.omit(cbind(ac, fx))
  actual <- as.numeric(acfx[,1])
  fted <- as.numeric(acfx[,2])
  np <- length(object$model$coef)
  # + sigma
  np <- np + 1
  ny <- attr(logLik(object$model), "nobs")
  order <- object$model$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  if (m > 1 && sum(order[4:6]) > 0) {
    frequency <- m
  } else {
    frequency <- 1
  }
  m_mape <- mape(actual, fted)
  m_mase <- mase(actual, fted, original_series = actual, frequency = frequency)
  m_mslre  <- mslre(actual, fted)
  m_bias <- bias(actual, fted)
  lik  <- object$model$loglik
  aic  <- -2 * lik + 2 * np
  bic  <- -2 * lik + log(ny) * np
  aicc <- aic + 2 * np * (np + 1) / (ny - np - 1)
  data.frame("n" = ny, "no.pars" = np - 1, "LogLik" = lik, "AIC" = aic, "BIC" = bic, "AICc" = aicc, "MAPE" = m_mape, "MASE" = m_mase, 
             "MSLRE" = m_mslre, "BIAS" = m_bias)
}

#' Object Plots
#'
#' @description Plots for objects generated from the bsts or arima functions.
#' @param x an object of class \dQuote{bsts.estimate} or \dQuote{arima.estimate}.
#' @param y not used.
#' @param ... no currently used.
#' @aliases plot
#' @method plot arima.estimate
#' @rdname plot
#' @export
#'
#'
plot.arima.estimate <- function(x, y = NULL, ...)
{
  opar <- par()
  opar$cin <- NULL
  opar$cra <- NULL
  opar$csi <- NULL
  opar$cxy <- NULL
  opar$din <- NULL
  opar$page <- NULL
  modelx <- arima_string(x$model)
  f <- fitted(x)
  actual <- xts(x$spec$target$y_orig, x$spec$target$index)
  plot(as.zoo(f), type = "l", ylab = "", xlab = "", col = "black", main = modelx, cex.main = 0.9, lwd = 1.5)
  lines(as.zoo(actual), col = "brown", lwd = 2, lty = 2)
  grid()
  mtext("Fitted", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
  legend("topleft",c("Fitted","Actual"), col = c("black","brown"), bty = "n", lty = c(1,2), lwd = c(1,0.5), cex = 0.8)
  mape_val <- tsmetrics(x)$MAPE * 100
  legend("bottomright",paste0("MAPE = ",round(mape_val,3),"%"), bty = "n", cex = 0.8, inset = c(0.02,.02))
  return(invisible(NULL))
}

#' Model Prediction
#'
#' @description Prediction function for class \dQuote{arima.estimate} or 
#' \dQuote{bsts.estimate}.
#' @param object an object of class \dQuote{arima.estimate} or \dQuote{bsts.estimate}.
#' @param h the forecast horizon.
#' @param newxreg a matrix of external regressors in the forecast horizon.
#' @param nsim The number of simulations to use for generating the simulated 
#' predictive distribution. For the bsts model, this is equal to the number of 
#' MCMC samples generated during estimation, less any burn-in draws.
#' @param forc_dates an optional vector of forecast dates equal to h. If NULL will
#' use the implied periodicity of the data to generate a regular sequence of
#' dates after the last available date in the data.
#' @param init_states an optional named vector which will re-center distribution
#' of the final state from which predictions are based off. It is best to call 
#' \link{bsts_final_state} function in order to get the full matrix and make any 
#' changes to the state_means prior to submitting to the function. It is required 
#' that the full vector is provided and minimal checks other than length are performed.
#' @param posterior_means optional vector of posterior parameter means which are 
#' then use to re-center the posterior parameters. In the case of AR and regressor 
#' coefficients, these should represent the means of the non-zero values (since 
#' both of these types of parameters are based on a spike and slab prior). The 
#' posterior distribution of the parameters can be obtain by calling 
#' \link{bsts_posterior}.
#' @param innov for the bsts and arima models this is an optional vector of uniform 
#' innovations which will be translated to regular innovations using the 
#' appropriate distribution quantile function and model standard deviation. The 
#' length of this vector should be equal to nsim x horizon for the arima model 
#' and MCMC draws (less burn) x horizon for the bsts model. Burn is equal to 
#' SuggestBurn(0.1, object$model) for the bsts model.
#' @param bootstrap for the arima model whether to bootstrap the residuals.
#' @param burn optional scalar denoting the numbers of draws to burn from the 
#' posterior prior to prediction.
#' @param ... not currently used.
#' @return An object which inherits class \dQuote{tsmodel.predict} with slots for 
#' the simulated or posterior predictive distribution, the original series (as a 
#' zoo object), the original specification object and the mean forecast. The 
#' predictive distribution is inversed difference (if differencing > 0) and back 
#' transformed if lambda was not NULL in the original specification. The innov 
#' argument for bsts is unlikely to be useful for ensembling, as has been done 
#' for other models with a single source of error as a result of multiple source 
#' of errors which disables our ability to infuse the required dependence structure 
#' via this approach. 
#' @aliases predict
#' @method predict arima.estimate
#' @rdname predict
#' @export
#'
#'
predict.arima.estimate = function(object, h = NULL, newxreg = NULL, nsim = 5000, 
                                  forc_dates = NULL, bootstrap = FALSE, 
                                  innov = NULL, ...)
{
  if (is.null(newxreg)) {
    if (is.null(h)) stop("\nhorizon (h) cannot be NULL when newxreg is NULL.")
    # need to use this form for exact matching....
    if (!is.null(object$spec$xreg[["xreg", exact = TRUE]])) stop("\nmodel has regressors (xreg) but no newxreg provided.")
  } else {
    if (nrow(newxreg) != h) stop("\nnrow newxreg not equal to horizon (h)")
  }
  if (is.null(object$spec$xreg[["xreg", exact = TRUE]])) {
    newxreg <- NULL
    if (is.null(forc_dates)) {
      forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
    }
  } else {
    if (!is.null(newxreg)) {
      forc_dates <- index(newxreg)
    } else {
      if (is.null(forc_dates)) {
        forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
      }
    }
  }
  xregf <- newxreg
  if (!is.null(innov)) {
    requiredn <- h * nsim
    if (length(innov) != requiredn) {
      stop("\nlength of innov must be nsim x h")
    }
    # check that the innovations are uniform nsim (from a copula)
    if (any(innov < 0 | innov > 1 )) {
      stop("\ninnov must be >0 and <1 (uniform nsim)")
    }
    if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
    if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
    innov <- matrix(innov, h, nsim)
  }
  if (object$spec$seasonal$seasonal) {
    if (object$spec$seasonal$seasonal_type == "trigonometric") {
      seasonal_fourier <- do.call.fast(cbind, lapply(1:length(object$spec$target$frequency), function(i) {
        fourier_series(forc_dates, period = object$spec$target$frequency[i], K = object$spec$seasonal$seasonal_harmonics[i])
      }))
      xregf <- cbind(coredata(xregf), seasonal_fourier)
    }
  }
  if (!is.null(xregf)) {
    xregf <- xregf[,object$spec$xreg$xreg_names, drop = FALSE]
  } else{
    xregf <- NULL
  }
  use_drift <- is.element("drift", names(object$model$coef))
  x <- object$x <- getResponse(object$model)
  usexreg <- (!is.null(xregf) | use_drift | is.element("xreg", names(object$model)))
  if (!is.null(xregf)) {
    origxreg <- xregf <- as.matrix(xregf)
  } else {
    origxreg <- NULL
  }
  if (use_drift) {
    n <- length(x)
    if (!is.null(xregf)) {
      xregf <- cbind((1:h) + n, xregf)
      colnames(xregf)[1] = "drift"
    } else {
      xregf <- as.matrix((1:h) + n)
      colnames(xregf)[1] = "drift"
    }
  }
  if (!is.null(object$model$constant)) {
    if (object$model$constant) {
      #pred <- list(pred = rep(x[1], h), se = rep(0, h))
    } else {
      stop("Strange value of object$constant")
    }
  } else if (usexreg) {
    if (is.null(xregf)) {
      stop("No regressors provided")
    }
    object$model$call$xreg <- getxreg(object$model)
    if (NCOL(xregf) != NCOL(object$model$call$xreg)) {
      stop("Number of regressors does not match fitted model")
    }
  }
  sim <- matrix(NA, nrow = nsim, ncol = h)
  for (i in 1:nsim) {
    sim[i, ] <- simulate2_Arima(object$model, nsim = h, bootstrap = bootstrap, xreg = xregf, lambda = NULL, innov = innov[,i])
  }
  colnames(sim) = as.character(forc_dates)
  if (!is.null(object$spec$transform)) {
    sim <- matrix(object$spec$transform$inverse(as.numeric(sim), object$spec$transform$lambda), ncol = ncol(sim), nrow = nrow(sim), byrow = FALSE)
    colnames(sim) <-  as.character(forc_dates)
  }
  class(sim) <- "tsmodel.distribution"
  mean_forecast = zoo(colMeans(sim), forc_dates)
  order <- object$model$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  if (m > 1 && sum(order[4:6]) > 0) {
    frequency <- m
  } else {
    frequency <- 1
  }
  zList = list(distribution = sim, original_series = zoo(object$spec$target$y_orig, object$spec$target$index),
               h = h, mean = mean_forecast, spec = object$spec, frequency = frequency)
  class(zList) <- c("arima.predict","tsmodel.predict")
  return(zList)
}

#' @method tsmetrics arima.predict
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.arima.predict = function(object, actual, alpha = 0.1, ...)
{
  if (is.null(list(...)$frequency)) {
    frequency <- object$frequency
  } else {
    frequency <- list(...)$frequency
  }
  n <- NCOL(object$distribution)
  if (NROW(actual) != n) stop("\nactual length not equal to forecast length")
  m_mape <- mape(actual, colMeans(object$distribution))
  m_bias <- bias(actual, colMeans(object$distribution))
  m_mslre <- mslre(actual, colMeans(object$distribution))
  m_mase <- mase(actual, colMeans(object$distribution), object$original_series, frequency = frequency)
  if (!is.null(alpha)) {
    m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
  } else {
    m_mis <- as.numeric(NA)
  }
  m_crps <- crps(actual, object$distribution)
  return(data.frame("h" = n, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias, "MIS" = m_mis,"CRPS" = m_crps))
}


################################################################
# start imports from forecast
getxreg <- function(z) 
{
  if (is.element("xreg", names(z))) {
    return(z$xreg)
  }
  else if (is.element("xreg", names(z$coef))) {
    return(eval.parent(z$coef$xreg))
  }
  else if (is.element("xreg", names(z$call))) {
    return(eval.parent(z$call$xreg))
  }
  else {
    armapar <- sum(z$arma[1:4]) + is.element("intercept", names(z$coef))
    npar <- length(z$coef)
    if (npar > armapar) {
      stop("It looks like you have an xreg component but I don't know what it is.\n  Please use Arima() or auto.arima() rather than arima().")
    }
    else {
      return(NULL)
    }
  }
}

myarima.sim <- function(model, n, x, e, ...) 
{
  start.innov <- residuals(model)
  innov <- e
  data <- x
  first.nonmiss <- which(!is.na(x))[1]
  if (first.nonmiss > 1) {
    tsp.x <- tsp(x)
    start.x <- tsp.x[1] + (first.nonmiss - 1)/tsp.x[3]
    x <- window(x, start = start.x)
    start.innov <- window(start.innov, start = start.x)
  }
  model$x <- x
  n.start <- length(x)
  x <- ts(c(start.innov, innov), start = 1 - n.start, frequency = model$seasonal.period)
  flag.noadjust <- FALSE
  if (is.null(tsp(data))) {
    data <- ts(data, frequency = 1, start = 1)
  }
  if (!is.list(model)) {
    stop("'model' must be list")
  }
  if (n <= 0L) {
    stop("'n' must be strictly positive")
  }
  p <- length(model$ar)
  q <- length(model$ma)
  d <- 0
  D <- model$seasonal.difference
  m <- model$seasonal.period
  if (!is.null(ord <- model$order)) {
    if (length(ord) != 3L) {
      stop("'model$order' must be of length 3")
    }
    if (p != ord[1L]) {
      stop("inconsistent specification of 'ar' order")
    }
    if (q != ord[3L]) {
      stop("inconsistent specification of 'ma' order")
    }
    d <- ord[2L]
    if (d != round(d) || d < 0) {
      stop("number of differences must be a positive integer")
    }
  }
  if (p) {
    minroots <- min(Mod(polyroot(c(1, -model$ar))))
    if (minroots <= 1) {
      stop("'ar' part of model is not stationary")
    }
  }
  if (length(model$ma)) {
    x <- stats::filter(x, c(1, model$ma), method = "convolution", sides = 1L)
    x[seq_along(model$ma)] <- 0
  }
  len.ar <- length(model$ar)
  if (length(model$ar) && (len.ar <= length(data))) {
    if ((D != 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
      diff.data <- diff(diff.data, lag = m, differences = D)
    } else if ((D != 0) && (d == 0)) {
      diff.data <- diff(data, lag = model$seasonal.period, differences = D)
    } else if ((D == 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
    } else {
      diff.data <- data
    }
    x.new.innovations <- x[(length(start.innov) + 1):length(x)]
    x.with.data <- c(diff.data, x.new.innovations)
    for (i in (length(diff.data) + 1):length(x.with.data)) {
      lagged.x.values <- x.with.data[(i - len.ar):(i - 1)]
      ar.coefficients <- model$ar[length(model$ar):1]
      sum.multiplied.x <- sum((lagged.x.values * ar.coefficients)[abs(ar.coefficients) > .Machine$double.eps])
      x.with.data[i] <- x.with.data[i] + sum.multiplied.x
    }
    x.end <- x.with.data[(length(diff.data) + 1):length(x.with.data)]
    x <- ts(x.end, start = 1, frequency = model$seasonal.period)
    flag.noadjust <- TRUE
  }
  else if (length(model$ar)) {
    x <- stats::filter(x, model$ar, method = "recursive")
  }
  if ((d == 0) && (D == 0) && (flag.noadjust == FALSE)) {
    if (n.start >= 20) {
      xdiff <- (model$x - x[1:n.start])[n.start - (19:0)]
    }
    else {
      xdiff <- model$x - x[1:n.start]
    }
    if (all(sign(xdiff) == 1) || all(sign(xdiff) == -1)) {
      xdiff <- xdiff[length(xdiff)]
    }
    else {
      xdiff <- mean(xdiff)
    }
    x <- x + xdiff
  }
  if ((n.start > 0) && (flag.noadjust == FALSE)) {
    x <- x[-(1:n.start)]
  }
  if ((D > 0) && (d == 0)) {
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)[-(1:length.s.xi)]
  } else if ((d > 0) && (D == 0)) {
    x <- diffinv(x, differences = d, xi = data[length(data) - (d:1) + 1])[-(1:d)]
  } else if ((d > 0) && (D > 0)) {
    delta.four <- diff(data, lag = m, differences = D)
    regular.xi <- delta.four[(length(delta.four) - D):length(delta.four)]
    x <- diffinv(x, differences = d, xi = regular.xi[length(regular.xi) - (d:1) + 1])[-(1:d)]
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)
    x <- x[-(1:length.s.xi)]
  }
  x <- ts(x[1:n], frequency = frequency(data), start = tsp(data)[2] + 1/tsp(data)[3])
  return(x)
}

simulate2_Arima = function(object, nsim = length(object$x), seed = NULL, xreg = NULL, future = TRUE, bootstrap = FALSE, innov = NULL, lambda = NULL, ...){
  # Error check:
  if (object$arma[7] < 0) {
    stop("Value for seasonal difference is < 0. Must be >= 0")
  }
  else if ((sum(object$arma[c(3, 4, 7)]) > 0) && (object$arma[5] < 2)) {
    stop("Invalid value for seasonal period")
  }
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    nsim <- nrow(xreg)
  }

  ####
  # Random Seed Code
  if (is.null(innov)) {
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }
    if (is.null(seed)) {
      RNGstate <- .Random.seed
    } else {
      R.seed <- .Random.seed
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
  } else {
    bootstrap <- FALSE
    nsim <- length(innov)
  }
  ############# End Random seed code
  # Check for seasonal ARMA components and set flag accordingly. This will be used later in myarima.sim()
  flag.s.arma <- (sum(object$arma[c(3, 4)]) > 0)
  # Check for Seasonality in ARIMA model
  if (sum(object$arma[c(3, 4, 7)]) > 0) {
    # return(simulateSeasonalArima(object, nsim=nsim, seed=seed, xreg=xreg, future=future, bootstrap=bootstrap, ...))
    if (sum(object$model$phi) == 0) {
      ar <- NULL
    }
    else {
      ar <- as.double(object$model$phi)
    }
    if (sum(object$model$theta) == 0) {
      ma <- NULL
    }
    else {
      ma <- as.double(object$model$theta)
    }
    order <- c(length(ar), object$arma[6], length(ma))
    if (future) {
      model <- list(
        order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object),
        seasonal.difference = object$arma[7], seasonal.period = object$arma[5], flag.seasonal.arma = flag.s.arma,
        seasonal.order = object$arma[c(3, 7, 4)]
      )
    }
    else {
      model <- list(order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object))
    }
    flag.seasonal.diff <- (object$arma[7] > 0)
  }
  else {
    #### Non-Seasonal ARIMA specific code: Set up the model
    order <- object$arma[c(1, 6, 2)]
    if (order[1] > 0) {
      ar <- object$model$phi[1:order[1]]
    } else {
      ar <- NULL
    }
    if (order[3] > 0) {
      ma <- object$model$theta[1:order[3]]
    } else {
      ma <- NULL
    }
    if (object$arma[2] != length(ma)) {
      stop("MA length wrong")
    } else if (object$arma[1] != length(ar)) {
      stop("AR length wrong")
    }
    if (future) {
      model <- list(
        order = object$arma[c(1, 6, 2)], ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object),
        seasonal.difference = 0, flag.seasonal.arma = flag.s.arma, seasonal.order = c(0, 0, 0), seasonal.period = 1
      )
    }
    else {
      model <- list(order = object$arma[c(1, 6, 2)], ar = ar, ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object))
    }
    flag.seasonal.diff <- FALSE
    ### End non-seasonal ARIMA specific code
  }
  x <- object$x <- getResponse(object)
  if (is.null(tsp(x))) {
    x <- ts(x, frequency = 1, start = 1)
  }
  n <- length(x)
  if (bootstrap) {
    res <- na.omit(c(model$residuals) - mean(model$residuals, na.rm = TRUE))
    e <- sample(res, nsim, replace = TRUE)
  } else if (is.null(innov)) {
    e <- rnorm(nsim, 0, model$sd)
  } else if (length(innov) == nsim) {
    e <- qnorm(innov, mean = 0, sd =  model$sd)
  } else {
    stop("Length of innov must be equal to nsim")
  }
  use.drift <- is.element("drift", names(object$coef))
  usexreg <- (!is.null(xreg) | use.drift)
  xm <- oldxm <- 0
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (nrow(xreg) < nsim) {
      stop("Not enough rows in xreg")
    } else {
      xreg <- xreg[1:nsim, ,drop = FALSE]
    }
  }
  if (use.drift) {
    # Remove existing drift column
    if (NCOL(xreg) == 1) {
      xreg <- NULL
    } else {
      xreg <- xreg[, !is.element(colnames(xreg), "drift"), drop = FALSE]
    }
    # Create new drift column for historical simulation
    dft <- as.matrix(1:nsim)
    # Adapt if future simulation
    if (future) {
      dft <- dft + n
    }
    # Add to xreg
    xreg <- cbind(drift = dft, xreg)
  }
  narma <- sum(object$arma[1L:4L])
  if (length(object$coef) > narma) {
    if (names(object$coef)[narma + 1L] == "intercept") {
      xreg <- cbind(intercept = rep(1, nsim), xreg)
      object$xreg <- cbind(intercept = rep(1, n), object$xreg)
    }
    if (!is.null(xreg)) {
      xm <- if (narma == 0) {
        drop(as.matrix(xreg) %*% object$coef)
      } else {
        drop(as.matrix(xreg) %*% object$coef[-(1L:narma)])
      }
      oldxm <- if (narma == 0) {
        drop(as.matrix(object$xreg) %*% object$coef)
      } else {
        drop(as.matrix(object$xreg) %*% object$coef[-(1L:narma)])
      }
    }
  }
  if (future) {
    sim <- myarima.sim(model, nsim, x - oldxm, e = e) + xm
  } else {
    if (flag.seasonal.diff) {
      zeros <- object$arma[5] * object$arma[7]
      sim <- arima.sim(model, nsim, innov = e)
      sim <- diffinv(sim, lag = object$arma[5], differences = object$arma[7])[-(1:zeros)]
      sim <- ts(tail(sim, nsim) + xm)
    }
    else {
      if (NCOL(xm) > 1) {
        xm = rowSums(xm)
      }
      sim <- ts(tail(arima.sim(model, nsim, innov = e), nsim) + xm)
    }
    tsp(sim) <- tsp(x)
    # If model is non-stationary, then condition simulated data on first observation
    if (model$order[2] > 0 || flag.seasonal.diff) {
      sim <- sim - sim[1] + x[1]
    }
  }
  # should always be NULL, we take care of the transformation separately
  if (!is.null(lambda)) {
    warning("\narima_sim lambda is not NULL...report to developer")
    sim <- InvBoxCox(sim, lambda)
  }
  return(sim)
}
# end imports from forecast
################################################################

#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{bsts.spec} or \dQuote{arima.spec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param alpha optional numeric vector of coverage rates for which to calculate
#' the quantiles.
#' @param n_iter number of MCMC iterations.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... not currently used.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals
#' \item metrics: a summary performance table showing metrics by
#' forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @aliases tsbacktest
#' @method tsbacktest arima.spec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.arima.spec <- function(object, start = floor(NROW(object$target$y_orig)/2), 
                                  end = NROW(object$target$y_orig), h = 1, alpha = NULL, trace = FALSE, ...)
{
  data <- xts(object$target$y_orig, object$target$index)
  transform <- object$transform
  lambda <- transform$lambda
  if (!is.null(object$xreg$xreg)) {
    use_xreg <- TRUE
    xreg <- xts(object$xreg$xreg, object$target$index)
  } else {
    use_xreg <- FALSE
    xreg <- NULL
  }
  
  start_date <- index(data)[start]
  end_date <- index(data)[end - 1]
  seqdates <- index(data[paste0(start_date,"/", end_date)])
  elapsed_time <- function(idx, end_date, start_date) {
    min(h, which(end_date == idx) - which(start_date == idx))
  }
  if (!is.null(alpha)) {
    if (any(alpha <= 0)) {
      stop("\nalpha must be strictly positive")
    }
    if (any(alpha >= 1)) {
      stop("\nalpha must be less than 1")
    }
    quantiles <- as.vector(sapply(1:length(alpha), function(k) c(alpha[k]/2, 1 - alpha[k]/2)))
  }
  # setup backtest indices
  horizon <- sapply(1:length(seqdates), function(i){
    min(h, elapsed_time(index(data), index(data)[end], seqdates[i]))
  })
  if (trace) {
      prog_trace <- progressor(length(seqdates))
  }
  b <- NULL
  b %<-% future_lapply(1:length(seqdates), function(i) {
    if (trace) prog_trace()
    ytrain <- data[paste0("/", seqdates[i])]
    ix <- which(index(data) == seqdates[i])
    ytest <- data[(ix + 1):(ix + horizon[i])]
    if (use_xreg) {
      xreg_train <- xreg[index(ytrain)]
      xreg_test <- xreg[index(ytest)]
    } else {
      xreg_train <- NULL
      xreg_test <- NULL
    }
    speclist <- object$call
    speclist$y <- ytrain
    speclist$xreg <- xreg_train
    spec <- do.call(arima_modelspec, args = speclist, quote = TRUE)
    mod <- estimate(spec)
    p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(ytest))
    if (!is.null(quantiles)) {
      qp <- apply(p$distribution, 2, quantile, quantiles)
      if (length(quantiles) == 1) {
        qp <- matrix(qp, ncol = 1)
      } else{
        qp <- t(qp)
      }
      colnames(qp) <- paste0("P", round(quantiles*100,1))
    }
    out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]), 
                      "horizon" = 1:horizon[i], 
                      "size" = rep(nrow(ytrain), horizon[i]),
                      "forecast_dates" = as.character(index(ytest)), 
                      "forecast" = as.numeric(p$mean), "actual" = as.numeric(ytest))
    if (!is.null(quantiles)) out <- cbind(out, qp)
    return(out)
  }, future.packages = c("tsmethods","tsaux","xts","tsforeign","forecast"), future.seed = TRUE)
  b <- eval(b)
  b <- rbindlist(b)
  if (is.null(data_name)) data_name <- "y"
  actual <- NULL
  forecast <- NULL
  metrics <- b[,list(variable = data_name, MAPE = mape(actual, forecast), MSLRE = mslre(actual, forecast),
                     BIAS = bias(actual, forecast), n = .N), by = "horizon"]
  if (!is.null(alpha)) {
    q_names <- matrix(paste0("P", round(quantiles*100,1)), ncol = 2, byrow = TRUE)
    q <- do.call(cbind, lapply(1:length(alpha), function(i){
      b[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = "horizon"]
    }))
    q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
    colnames(q) <- paste0("MIS[",alpha,"]")
    metrics <- cbind(metrics, q)
  }
  return(list(prediction = b, metrics = metrics))
}

arima_string <- function(object)
{
  order <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 && is.element("drift", names(object$coef))) {
      result <- paste(result, "with drift        ")
    }
    else {
      result <- paste("Regression with", result, "errors")
    }
  }
  else {
    if (is.element("constant", names(object$coef)) || is.element("intercept", names(object$coef))) {
      result <- paste(result, "with non-zero mean")
    }
    else if (order[2] == 0 && order[5] == 0) {
      result <- paste(result, "with zero mean    ")
    }
    else {
      result <- paste(result, "                  ")
    }
  }
  result <- gsub("[ ]*$", "", result)
  return(result)
}
