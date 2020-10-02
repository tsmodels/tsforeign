bsts_modelspec = function(y, xreg = NULL, frequency = NULL, differences = 0, level = TRUE, slope = TRUE, damped = FALSE, seasonal = FALSE, 
                          seasonal_frequency = 4, ar = FALSE, ar_max = 1, cycle = FALSE, cycle_frequency = NULL, cycle_names = NULL, 
                          seasonal_type = "regular", lambda = NULL, lambda_lower = 0, lambda_upper = 1,
                          seasonal_harmonics = NULL, distribution = "gaussian", ...)
{
  model <- list(frequency = frequency, differences = differences, level = level, slope = slope, damped = damped, seasonal = seasonal, seasonal_frequency = seasonal_frequency, 
               ar = ar, ar_max = ar_max, cycle = cycle, cycle_frequency = cycle_frequency, cycle_names = cycle_names, seasonal_type = seasonal_type, lambda = lambda,
               seasonal_harmonics = seasonal_harmonics, distribution = distribution)
  # 1. Check y
  if (!is.xts(y)) {
    stop("y must be an xts object")
  }
  if (any(is.na(y))) {
    stop("\nNAs found in y...not allowed.")
  }
  # 2. Check regressors
  xreg <- check_xreg(xreg, index(y))
  
  valid_distributions <- c("gaussian")
  distribution <- match.arg(distribution[1], valid_distributions)

  # 3. Check transformation
  y_orig <- y
  if (!is.null(lambda)) {
    if (!is.na(lambda) & lambda == 1) lambda <- NULL
  }
  if (!is.null(lambda)) {
    transform <- box_cox(lambda = lambda, lower = lambda_lower, upper = lambda_upper)
    y <- transform$transform(y = y, frequency = frequency)
    transform$lambda <- attr(y, "lambda")
  } else{
    transform <- NULL
  }
  if (differences > 0 & distribution == "gaussian") {
    Dy <- na.omit(diff(y, differences = differences))
    if (!is.null(xreg)) {
      xreg <- xreg[index(Dy)]
      xreg <- check_xreg(xreg, index(Dy))
    }
  } else {
    differences <- 0
    Dy <- y
  }
  statespec <- c()
  b_spec <- list()
  s_spec <- NULL
  if (level & slope) {
    if (damped) {
      b_spec <- AddSemilocalLinearTrend(b_spec, as.numeric(Dy))
    } else {
      b_spec <- AddLocalLinearTrend(b_spec, as.numeric(Dy))
    }
  } else if (level & !slope) {
    b_spec <- AddLocalLevel(b_spec, as.numeric(Dy))
  }
  
  if (seasonal) {
    if (seasonal_type == "regular") {
      for (i in 1:length(seasonal_frequency)) {
        b_spec <- AddSeasonal(b_spec, as.numeric(Dy), nseasons = seasonal_frequency[i], season.duration = 1)
      }
    } else {
      if (is.null(seasonal_harmonics)) {
        seasonal_harmonics <- floor(0.25 * seasonal_frequency)
      } else {
        if (length(seasonal_harmonics) != length(seasonal_frequency)) {
          warnings("\nlength of seasonal_harmonics not equal to length of seasonal_frequency. Defaulting to 0.25*frequency.")
          seasonal_harmonics <- floor(0.25*seasonal_frequency)
        }
      }
      for (i in 1:length(seasonal_frequency)) {
        b_spec <- AddTrig(b_spec, as.numeric(Dy), period = seasonal_frequency[i], frequencies = 1:seasonal_harmonics[i])
      }
    }
  }
  
  if (cycle) {
    if (is.null(cycle_frequency)) stop("\ncycle_frequency cannot be NULL when cycle is TRUE")
    nc <- length(cycle_frequency)
    if (seasonal) {
      if (any(cycle_frequency <= min(seasonal_frequency))) {
        stop("\ncycle length is less than of equal to min(seasonal) seasonal_frequency...try again!")
      }
    }
    if (!is.null(cycle_names)) {
      if (length(cycle_names) != nc) stop("\ncycle.names length not equal to length of cycle_periods")
    } else {
      cycle_names <- paste0("Cycle_", 1:nc)
    }
    for (i in 1:nc) {
      b_spec <- AddTrig(b_spec, y = as.numeric(Dy), period = cycle_frequency[i], frequencies = 1)
     }
  }
  if (ar) {
    b_spec <- AddAutoAr(b_spec, as.numeric(Dy), lags = ar_max)
  }
  if (distribution == "gaussian") {
    s_matrices <- generate_ss_model.gaussian(level = level, slope = slope, damped = damped, seasonal = seasonal, seasonal_frequency = seasonal_frequency, 
                                             seasonal_type = seasonal_type, seasonal_harmonics = seasonal_harmonics, 
                                             cycle = cycle, cycle_frequency = cycle_frequency, ar = ar, ar_max = ar_max, xreg =  xreg)
  } else {
    # future enhancements
    stop("\nonly gaussian distribution currently implemented.")
  }
  if (!is.null(xreg)) {
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("X",1:ncol(xreg))
    }
    equation <- paste0("y~", paste(colnames(xreg), collapse = "+"), "-1")
  } else {
    equation <- paste0("y")
  }
  spec <- list()
  spec$target$Dy <- as.numeric(Dy)
  spec$target$y <- as.numeric(y)
  spec$target$y_orig <- as.numeric(y_orig)
  spec$target$index <- index(y_orig)
  spec$target$Dindex <- index(Dy)
  spec$target$frequency <- frequency
  spec$target$differences <- differences
  spec$target$sampling <- sampling_frequency(index(y_orig))
  spec$model <- model
  spec$transform <- transform
  if (!is.null(xreg)) {
    spec$xreg$xreg <- coredata(xreg)
    spec$xreg$index <- index(xreg)
    spec$xreg$xreg_names <- colnames(xreg)
  } else {
    spec$xreg <- NULL
  }
  spec$bsts_spec <- b_spec
  spec$bsts_equation <- equation
  spec$state_space <- s_matrices
  class(spec) <- c("bsts.spec","tsmodel.spec")
  return(spec)
}


estimate.bsts.spec = function(object, n_iter = 5000, timeout.seconds = Inf, bma.method = "SSVS", trace = TRUE, ...)
{
  model.opt <- BstsOptions(save.state.contributions = TRUE, save.prediction.errors = TRUE, bma.method = bma.method, 
                           oda.options = list(fallback.probability = 0.0, eigenvalue.fudge.factor = 0.01), 
                           timeout.seconds = timeout.seconds, save.full.state = TRUE)
  dat <- cbind(object$target$Dy, object$xreg$xreg)
  if (trace) {
    ping <- n_iter / 10
  } else {
    ping <- 0
  }
  colnames(dat)[1] <- "y"
  if (!is.null(object$xreg$xreg)) {
    mod <- bsts(object$bsts_equation, state.specification = object$bsts_spec, data = as.data.frame(dat), niter = n_iter, family = object$model$distribution, model.options = model.opt, ping = ping, ...)
  } else{
    mod <- bsts(dat[,"y"], state.specification = object$bsts_spec, data = as.data.frame(dat), niter = n_iter, family = object$model$distribution, model.options = model.opt, ping = ping, ...)
  }
  obj <- list(model = mod, spec = object)
  class(obj) <- c("bsts.estimate","tsmodel.estimate")
  return(obj)
}

summary.bsts.estimate <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
{
  out <- bsts_posterior(object)
  if (any(colnames(out) == "slope.oneminusar")) {
    exc <- which(colnames(out) == "slope.oneminusar")
    out <- out[ ,-exc]
  }
  coda_out <- as.mcmc(out)
  summary_out <- summary(coda_out, quantiles = quantiles)
  # check xreg and ar for inclusion probabilities
  if (object$spec$model$ar) {
    ar_out <- out[,grepl("ar[0-9]$", colnames(out)), drop = FALSE]
    ar_inclusion <- apply(ar_out, 2, function(x) sum(abs(x) > 0)/nrow(ar_out))
    summary_out$statistics <- cbind(summary_out$statistics, matrix(1, nrow = nrow(summary_out$statistics), dimnames = list(rownames(summary_out$statistics),"Prob[Include]")))
    summary_out$statistics[grepl("ar[0-9]$", colnames(out)),"Prob[Include]"] <- ar_inclusion
  }
  if (!is.null(object$spec$xreg$xreg)) {
    xreg_out <- out[,paste0("beta",1:ncol(object$spec$xreg$xreg)), drop = FALSE]
    xreg_inclusion <- apply(xreg_out, 2, function(x) sum(abs(x) > 0)/nrow(xreg_out))
    if (any(colnames(summary_out$statistics) == "Prob[Include]")) {
      summary_out$statistics[paste0("beta",1:ncol(object$spec$xreg$xreg)),"Prob[Include]"] <- xreg_inclusion
    } else {
      summary_out$statistics <- cbind(summary_out$statistics, matrix(1, nrow = nrow(summary_out$statistics), dimnames = list(rownames(summary_out$statistics),"Prob[Include]")))
      summary_out$statistics[paste0("beta",1:ncol(object$spec$xreg$xreg)),"Prob[Include]"] <- xreg_inclusion
    }
  }
  print(summary_out)
  cat("\nHarvey's Goodness of Fit Statistic:", hgof.bsts.estimate(object$model, d = 0),"\n")
  mtrcs <- tsmetrics(object)
  cat("\n")
  print(data.frame(MAPE = as.numeric(sprintf(mtrcs$MAPE,fmt = "%.4f")), MASE = as.numeric(sprintf(mtrcs$MASE,fmt = "%.4f")),
                   MSLRE = as.numeric(sprintf(mtrcs$MSLRE,fmt = "%.4f")), BIAS = as.numeric(sprintf(mtrcs$BIAS,fmt = "%.4f"))), 
        row.names = FALSE, right = T)
  return(invisible(summary_out))
}

tsmetrics.bsts.estimate <- function(object, ...)
{
  fx <- fitted(object)
  ac <- xts(object$spec$target$y_orig, object$spec$target$index)
  acfx <- na.omit(cbind(ac, fx))
  actual <- as.numeric(acfx[,1])
  fted <- as.numeric(acfx[,2])
  if (object$spec$model$seasonal) {
    frequency <- object$spec$target$frequency
  } else {
    frequency <- 1
  }
  m_mape <- mape(actual, fted)
  m_mase <- mase(actual, fted, original_series = actual, frequency = frequency)
  m_mslre  <- mslre(actual, fted)
  m_bias <- bias(actual, fted)
  ny <- NROW(object$spec$target$Dy)
  np <- unique(object$spec$state_space$coef_names)
  if (any(np == "slope.oneminusar")) np <- np[-which(np == "slope.oneminusar")]
  np <- length(np)
  # add 1 for the obs.sigma
  return(data.frame("n" = ny, "no.pars" = np + 1, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias))
}

fitted.bsts.estimate = function(object, distribution = FALSE, invdiff = TRUE, raw = FALSE, type = c("filtered","smoothed"), ...)
{
  type <- match.arg(type[1], c("filtered","smoothed"))
  burn <- SuggestBurn(0.1, object$model)
  state <- object$model$state.contributions
  sig2 <- object$model$sigma.obs^2
  if (burn > 0) {
    burn <- 1:burn
  } else {
    burn <- -(1:length(sig2))
  }
  state <- state[-burn, , , drop = FALSE]
  sig2 <- sig2[-burn]
  if (type == "filtered") {
    error <- object$model$one.step.prediction.errors[-burn, , drop = FALSE]
    out <- matrix(object$spec$target$Dy, nrow = nrow(error), ncol = ncol(error), byrow = TRUE) - error
  } else {
    out <- rowSums(aperm(state, c(1, 3, 2)), dims = 2)
  }
  if (object$spec$target$differences > 0 & invdiff) {
    if (!raw) {
      yin <- xts(object$spec$target$y_orig, object$spec$target$index)
      out <- do.call.fast(cbind, lapply(1:ncol(out), function(i) {
        coredata(d2levels_fitted_matrix(act = yin[i:(i + 1)], pred = out[,i], d = object$spec$target$differences, transform = object$spec$transform, frequency = object$spec$target$frequency))
      }))
      if (object$spec$target$differences == 1) {
        out <- cbind(matrix(as.numeric(object$spec$target$y_orig[1]), ncol = 1, nrow = nrow(out)), out)
      } else {
        out <- cbind(matrix(as.numeric(object$spec$target$y_orig[1]), ncol = 1, nrow = nrow(out)),
                     matrix(as.numeric(object$spec$target$y_orig[2]), ncol = 1, nrow = nrow(out)), out)
      }
    } else {
      yin <- xts(object$spec$target$y, object$spec$target$index)
      out <- do.call.fast(cbind, lapply(1:ncol(out), function(i) {
        coredata(d2levels_fitted_matrix(act = yin[i:(i + 1)], pred = out[,i], d = object$spec$target$differences, transform = NULL, frequency = object$spec$target$frequency))
      }))
      if (object$spec$target$differences == 1) {
        out <- cbind(matrix(as.numeric(object$spec$target$y[1]), ncol = 1, nrow = nrow(out)), out)
      } else {
        out <- cbind(matrix(as.numeric(object$spec$target$y[1]), ncol = 1, nrow = nrow(out)),
                     matrix(as.numeric(object$spec$target$y[2]), ncol = 1, nrow = nrow(out)), out)
      }
    }
  } else{
    if (!is.null(object$spec$transform) & invdiff & !raw) {
      out <- matrix(object$spec$transform$inverse(as.numeric(out), object$spec$transform$lambda), ncol = ncol(out), nrow = nrow(out), byrow = FALSE)
    }
  }
  if (object$spec$target$differences == 0 | invdiff) {
    colnames(out) <- as.character(object$spec$target$index)
    rownames(out) <- NULL
    if (!distribution) {
      out <- xts(matrix(apply(out, 2, mean), ncol = 1), object$spec$target$index)
    } else{
      class(out) <- c("bsts.distribution","tsmodel.distribution")
      attr(out,'burn') <- burn
    }
  } else{
    colnames(out) <- as.character(object$spec$target$Dindex)
    rownames(out) <- NULL
    if (!distribution) {
      out <- xts(matrix(apply(out, 2, mean), ncol = 1), object$spec$target$Dindex)
    } else{
      class(out) <- c("bsts.distribution","tsmodel.distribution")
      attr(out,'burn') <- burn
    }
  }
  #if (object$spec$distribution == "poisson") {
  #  out <- exp(out)
  #}
  return(out)
}

residuals.bsts.estimate = function(object, distribution = FALSE, invdiff = TRUE, standardize = FALSE, raw = FALSE, type = c("filtered", "smoothed"), ...)
{
  if (distribution) {
    if (standardize) {
      raw <- TRUE
      invdiff <- FALSE
      out <- fitted(object, distribution = TRUE, type = type, invdiff = invdiff, raw = raw)
      burn <- attr(out, 'burn')
      a <- matrix(as.numeric(object$spec$target$Dy), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
      out <- a - out
      sig <- object$model$sigma.obs
      if (length(burn) > 0 & burn[1] > 0) sig <- sig[-burn]
      out <- do.call.fast(cbind, lapply(1:ncol(out), function(i) out[,i]/sig))
    } else {
      if (!invdiff) {
        a <- matrix(as.numeric(object$spec$target$Dy), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
      } else {
        if (raw) {
          a <- matrix(as.numeric(object$spec$target$y), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
        } else {
          a <- matrix(as.numeric(object$spec$target$y_orig), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
        }
      }
      out <- a - out
    }
    if (!invdiff) {
      colnames(out) <- as.character(object$spec$target$Dindex)
    } else {
      colnames(out) <- as.character(object$spec$target$index)
    }
    class(out) <- c("bsts.distribution","tsmodel.distribution")
  } else{
    if (standardize) {
      burn <- SuggestBurn(0.1, object$model)
      f <- fitted(object, distribution = FALSE, invdiff = FALSE, type = type, raw = TRUE)
      a <- as.numeric(object$spec$target$Dy)
      sig <- object$model$sigma.obs
      if (burn > 0) sig <- sig[-c(1:burn)]
      sig <- mean(sig)
      out <- (a - f)/sig
    } else {
      f <- fitted(object, distribution = FALSE, invdiff = invdiff, type = type, raw = raw)
      if (!invdiff) {
        a <- as.numeric(object$spec$target$Dy)
      } else {
        if (raw) {
          a <- as.numeric(object$spec$target$y)
        } else {
          a <- as.numeric(object$spec$target$y_orig)
        }
      }
      if (!invdiff) {
        out <- xts(a, object$spec$target$Dindex) - f
      } else {
        out <- xts(a, object$spec$target$index) - f
      }
    }
    colnames(out) <- c("residuals")
  }
  return(out)
}

predict.bsts.estimate <- function(object, h = 1, newxreg  = NULL, forc_dates = NULL, last_state_means = colMeans(bsts_final_state(object)), posterior_means = NULL, innov = NULL, burn = NULL, ...)
{
  if (is.null(burn)) {
    burn <- SuggestBurn(0.1, object$model)
  }
  if (burn > 0) {
    burn <- 1:burn
    N <- length(object$model$sigma.obs[-burn])
  } else {
    N <- length(object$model$sigma.obs)
    burn <- -1 * (1:N)
  }
  requiredn <- h * N
  
  final_state <- bsts_final_state(object)
  if (!is.null(last_state_means)) {
    if (length(last_state_means) != ncol(final_state)) {
      stop("\nlast_state_means length not equal to ncol(final.state)")
    } else {
      final_state <- scale(final_state, scale = FALSE)
      final_state <- scale(final_state, center = -1 * last_state_means, scale = FALSE)
    }
  }
  final_state <- final_state[-burn,,drop = FALSE]
  if (!is.null(innov)) {
    if (length(innov) != requiredn) {
      stop(paste0("\nlength of innov must be MCMC Draws [-burn] x h :", requiredn))
    }
    if (any(innov < 0 | innov > 1 )) {
      stop("\ninnov must be >0 and <1 (uniform samples)")
    }
    if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
    if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
    innov <- t(qnorm(matrix(innov, N, h)))
  } else {
    innov <- matrix(rnorm(requiredn), h, N)
  }
  innov <- rbind(matrix(0, nrow = 1, ncol = ncol(innov)), innov)
  
  P <- bsts_posterior(object)
  P <- P[-burn, ,drop = FALSE]
  if (!is.null(posterior_means)) {
    posterior_means_names <- names(posterior_means)
    p_names <- colnames(P)
    if (is.null(posterior_means_names)) stop("\nposterior_means must be a named vector")
    # eliminate 1-slope.ar (this is a calculated value)
    if (any(posterior_means_names == "slope.oneminusar")) {
      exc <- which(posterior_means_names == "slope.oneminusar")
      posterior_means <- posterior_means[-exc]
      posterior_means_names <- names(posterior_means)
    }
    # match names and for AR and XREG retain inclusion probability by only centering non zero values
    posterior_means <- posterior_means[which(posterior_means_names %in% p_names)]
    if (is.null(posterior_means)) stop("\ncould not match any of the posterior_means names to the posterior names")
    posterior_means_names <- names(posterior_means)
    for (i in 1:length(posterior_means)) {
      use <- which(p_names == posterior_means_names[i])
      if (grepl("ar[0-9]$",posterior_means_names[i])) {
        if (any(P[,use] != 0)) {
          P[which(P[,use] != 0), use] <- (P[which(P[,use] != 0), use] - mean(P[which(P[,use] != 0), use])) + posterior_means[i]
        }
      } else if (grepl("beta[0-9]$",posterior_means_names[i])) {
        if (any(P[,use] != 0)) {
          P[which(P[,use] != 0), use] <- (P[which(P[,use] != 0), use] - mean(P[which(P[,use] != 0), use])) + posterior_means[i]
        }
      } else if (grepl("slope.ar",posterior_means_names[i])) {
        P[,use] <- (P[,use] - mean(P[,use])) + posterior_means[i]
        P[,which(p_names == "slope.oneminusar")] <- 1 - P[,use]
      } else {
        P[,use] <- (P[,use] - mean(P[,use])) + posterior_means[i]
      }
    }
  }
  if (is.null(newxreg)) {
    if (is.null(h)) stop("\nhorizon (h) cannot be NULL when newxreg is NULL.")
    if (!is.null(object$spec$xreg$xreg)) stop("\nmodel has regressors (xreg) but no newxreg provided.")
  } else {
    if (nrow(newxreg) != h) stop("\nnrow newxreg not equal to horizon (h)")
  }
  if (is.null(object$spec$xreg$xreg)) {
    newxreg <- NULL
    if (is.null(forc_dates)) {
      forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
    }
    newxreg <- matrix(0, nrow = h, ncol = 1)
    B <- matrix(0, ncol = N, nrow = 1)
    X <- NULL
  } else {
    newxreg <- check_newxreg(newxreg, xnames = object$spec$xreg$xreg_names, h = h, forc_dates = NULL)
    if (!is.null(newxreg)) {
      forc_dates <- index(newxreg)
    } else {
      if (is.null(forc_dates)) {
        forc_dates <- future_dates(tail(object$spec$target$index,1), frequency = object$spec$target$sampling, n = h)
      }
    }
    B <- t(P[,which(grepl("beta[0-9]$",colnames(P))),drop = FALSE])
    X <- t(coredata(newxreg) %*% B)
  }
  xregf <- coredata(newxreg)
  xregf <- rbind(matrix(0, nrow = 1, ncol = ncol(xregf)), xregf)
  
  Z <- object$spec$state_space$Z
  G <- object$spec$state_space$G
  Q <- object$spec$state_space$Q
  R <- object$spec$state_space$R
  init_state <- t(final_state)
  gpriors <- P[, object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "G")], drop = FALSE]
  gp_names <- colnames(gpriors)
  # for gpriors need to replace slope.mean with slope.ar and slope.ar with 1-slope.phi
  gpriors <- t(gpriors)
  qpriors <- t(P[, object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "W")], drop = FALSE])
  # square the sigmas to get variances
  qpriors <- qpriors^2
  
  vpriors <- P[,"obs.sigma"]
  
  model <- c(h, N, ncol(G))
  
  bp <- bsts_posterior_predict(model, Z, G, Q, R, init_state, xregf, B, gpriors, qpriors, vpriors, innov)
  
  if (object$spec$target$differences > 0) {
    p <- do.call(rbind, lapply(1:nrow(bp$y), function(i){
      d2levels_forecast(object$spec$target$y, bp$y[i,], d = object$spec$target$differences, transform = object$spec$transform, frequency = 1)
    }))
    meanf <- zoo(apply(p, 2, mean), forc_dates)
  } else {
    if (!is.null(object$spec$transform)) {
      p <- object$spec$transform$inverse(bp$y, object$spec$transform$lambda)
      meanf <- zoo(apply(p, 2, mean), forc_dates)
    } else{
      p <- bp$y
      meanf <- zoo(apply(p, 2, mean), forc_dates)
    }
  }
  colnames(p) <- as.character(forc_dates)
  class(p) <- "tsmodel.distribution"
  attr(p, "date_class") <- attr(object$spec$target$sampling, "date_class")
  states <- bp$states
  # names the states
  L <- list(original_series = xts(object$spec$target$y_orig, object$spec$target$index), distribution = p, mean = meanf, states = states, burn = burn, X = X, spec = object$spec)
  class(L) <- c("bsts.predict","tsmodel.predict")
  return(L)
}

tsdecompose.bsts.estimate <- function(object, ...)
{
  burn <- SuggestBurn(0.1, object$model)
  if (burn > 0) {
    burn <- 1:burn
    out <- object$model$full.state[-burn,,, drop = FALSE]
  } else {
    out <- object$model$full.state
  }
  c_names <- object$spec$state_space$component_names
  d_class <- attr(object$spec$target$sampling, "date_class")
  if (object$spec$model$level) {
    Level <- out[,which(c_names == "Level"),]
    colnames(Level) <- as.character(object$spec$target$Dindex)
    class(Level) <- "tsmodel.distribution"
    attr(Level,"date_class") <- d_class
  } else {
    Level <- NULL
  }
  if (object$spec$model$slope) {
    Slope <- out[,which(c_names == "Slope"),]
    colnames(Slope) <- as.character(object$spec$target$Dindex)
    class(Slope) <- "tsmodel.distribution"
    attr(Slope,"date_class") <- d_class
  } else {
    Slope <- NULL
  }
  if (object$spec$model$seasonal) {
    if (object$spec$model$seasonal_type == "regular") {
      Seasonal <- vector(mode = "list", length = length(object$spec$model$seasonal_frequency))
      for (i in 1:length(object$spec$model$seasonal_frequency)) {
        Seasonal[[i]] <- out[,which(c_names == paste0("Seasonal.0.",object$spec$model$seasonal_frequency[i])),]
        colnames(Seasonal[[i]]) <- as.character(object$spec$target$Dindex)
        class(Seasonal[[i]]) <- "tsmodel.distribution"
        attr(Seasonal[[i]],"date_class") <- d_class       
      }
    } else {
      Seasonal <- vector(mode = "list", length = length(object$spec$model$seasonal_frequency))
      for (i in 1:length(object$spec$model$seasonal_frequency)) {
        indx <- which(grepl(paste0("Trigonometric.",object$spec$model$seasonal_frequency[i]), c_names))
        tmp <- out[,indx,]
        indx <- seq(1, length(indx), by = 2)
        Seasonal[[i]] <- apply(tmp, 3, function(x) rowSums(x[,indx,drop = FALSE]))
        colnames(Seasonal[[i]]) <- as.character(object$spec$target$Dindex)
        class(Seasonal[[i]]) <- "tsmodel.distribution"
        attr(Seasonal[[i]],"date_class") <- d_class       
      }
    }
    names(Seasonal) <- paste0("Seasonal",object$spec$model$seasonal_frequency)
  } else {
    Seasonal <- NULL   
  }
  if (object$spec$model$cycle) {
    Cyclical <- vector(mode = "list", length = length(object$spec$model$cycle_frequency))
    for (i in 1:length(object$spec$model$cycle_frequency)) {
      indx <- which(grepl(paste0("Cyclical.",object$spec$model$cycle_frequency[i]), c_names))
      tmp <- out[,indx,]
      indx <- seq(1, length(indx), by = 2)
      Cyclical[[i]] <- apply(tmp, 3, function(x) rowSums(x[,indx,drop = FALSE]))
      colnames(Cyclical[[i]]) <- as.character(object$spec$target$Dindex)
      class(Cyclical[[i]]) <- "tsmodel.distribution"
      attr(Cyclical[[i]],"date_class") <- d_class       
    }
    names(Cyclical) <- paste0("Cyclical",object$spec$model$cycle_frequency)
  } else {
    Cyclical <- NULL
  }
  if (object$spec$model$ar) {
    indx <- which(grepl("AR.1", c_names))
    AR <- out[,indx,]
    colnames(AR) <- as.character(object$spec$target$Dindex)
    class(AR) <- "tsmodel.distribution"
    attr(AR,"date_class") <- d_class
  } else {
    AR <- NULL
  }
  if (!is.null(object$spec$xreg$xreg)) {
    B <- object$model$coefficients[-burn, , drop = FALSE]
    X <- object$spec$xreg$xreg
    X <- B %*% t(X)
    colnames(X) <- as.character(object$spec$target$Dindex)
    class(X) <- "tsmodel.distribution"
    attr(X,"date_class") <- d_class
  } else {
    X <- NULL
  }
  L <- list(Level = Level, Slope = Slope, AR = AR, X = X)
  L <- c(L, Seasonal, Cyclical)
  return(L)
}


tsdecompose.bsts.predict <- function(object, ...)
{
  c_names <- object$spec$state_space$component_names
  d_class <- attr(object$spec$target$sampling, "date_class")
  f_dates <- colnames(object$distribution)
  states <- object$states
  if (object$spec$model$level) {
    Level <- t(states[which(c_names == "Level"),,])
    colnames(Level) <- f_dates
    class(Level) <- "tsmodel.distribution"
    attr(Level,"date_class") <- d_class
  } else {
    Level <- NULL
  }
  if (object$spec$model$slope) {
    Slope <- t(states[which(c_names == "Slope"),,])
    colnames(Slope) <- f_dates
    class(Slope) <- "tsmodel.distribution"
    attr(Slope,"date_class") <- d_class
  } else {
    Slope <- NULL
  }
  if (object$spec$model$seasonal) {
    if (object$spec$model$seasonal_type == "regular") {
      Seasonal <- vector(mode = "list", length = length(object$spec$model$seasonal_frequency))
      for (i in 1:length(object$spec$model$seasonal_frequency)) {
        Seasonal[[i]] <- t(states[which(c_names == paste0("Seasonal.0.",object$spec$model$seasonal_frequency[i])),,])
        colnames(Seasonal[[i]]) <- f_dates
        class(Seasonal[[i]]) <- "tsmodel.distribution"
        attr(Seasonal[[i]],"date_class") <- d_class       
      }
    } else {
      Seasonal <- vector(mode = "list", length = length(object$spec$model$seasonal_frequency))
      for (i in 1:length(object$spec$model$seasonal_frequency)) {
        indx <- which(grepl(paste0("Trigonometric.",object$spec$model$seasonal_frequency[i]), c_names))
        tmp <- states[indx,,]
        indx <- seq(1, length(indx), by = 2)
        Seasonal[[i]] <- t(apply(tmp, 3, function(x) colSums(x[indx,,drop = FALSE])))
        colnames(Seasonal[[i]]) <- f_dates
        class(Seasonal[[i]]) <- "tsmodel.distribution"
        attr(Seasonal[[i]],"date_class") <- d_class       
      }
    }
    names(Seasonal) <- paste0("Seasonal",object$spec$model$seasonal_frequency)
  } else {
    Seasonal <- NULL   
  }
  if (object$spec$model$cycle) {
    Cyclical <- vector(mode = "list", length = length(object$spec$model$cycle_frequency))
    for (i in 1:length(object$spec$model$cycle_frequency)) {
      indx <- which(grepl(paste0("Cyclical.",object$spec$model$cycle_frequency[i]), c_names))
      tmp <- t(states[indx,,])
      indx <- seq(1, length(indx), by = 2)
      Cyclical[[i]] <- apply(tmp, 3, function(x) rowSums(x[,indx,drop = FALSE]))
      colnames(Cyclical[[i]]) <- f_dates
      class(Cyclical[[i]]) <- "tsmodel.distribution"
      attr(Cyclical[[i]],"date_class") <- d_class       
    }
    names(Cyclical) <- paste0("Cyclical",object$spec$model$cycle_frequency)
  } else {
    Cyclical <- NULL
  }
  if (object$spec$model$ar) {
    indx <- which(grepl("AR.1", c_names))
    AR <- t(states[indx,,])
    colnames(AR) <- f_dates
    class(AR) <- "tsmodel.distribution"
    attr(AR,"date_class") <- d_class
  } else {
    AR <- NULL
  }
  if (!is.null(object$spec$xreg$xreg)) {
    X <- object$X
    colnames(X) <- f_dates
    class(X) <- "tsmodel.distribution"
    attr(X,"date_class") <- d_class
  } else {
    X <- NULL
  }
  L <- list(Level = Level, Slope = Slope, AR = AR, X = X)
  L <- c(L, Seasonal, Cyclical)
  return(L)
}


plot.bsts.estimate <- function(x, y = NULL, ...)
{
  # Fitted + Residuals
  tsx <- tsdecompose(x)
  n <- sum(sapply(1:length(tsx), function(i) as.integer(!is.null(tsx[[i]]))))
  n <- n + 2
  par(mar = c(3,3,3,3))
  N <- n2mfrow(n)
  nn <- N[1]*N[2]
  layout(matrix(1:nn, n2mfrow(n)))
  plot(fitted(x, distribution = TRUE), main = "Fitted vs Actual", zero_out = FALSE)
  lines(x$spec$target$index, as.numeric(x$spec$target$y_orig), col = 2, lwd = 2)
  plot(residuals(x, distribution = TRUE, standardize = TRUE), main = "Residuals[scaled]",  zero_out = FALSE)
  if (n > 2) {
    cnames <- names(tsx)
    for (i in 1:length(tsx)) {
      if (!is.null(tsx[[i]])) plot(tsx[[i]], main = cnames[i],  zero_out = FALSE)
    }
  }
  return(invisible(tsx))
}

# Harvey's relative goodness of fit statistic
hgof.bsts.estimate <- function(object, d = 0)
{
  if (d > 0) {
    burn <- SuggestBurn(0.1, object)
    prediction_errors <- bsts.prediction.errors(object, burn = burn)$in.sample
    prediction_sse <- sum(colMeans(prediction_errors)^2)
    dy <- as.numeric(object$original.series)
    prediction_sst <- var(dy) * (length(dy) - 1)
    relative_gof = 1 - prediction_sse / prediction_sst
  } else{
    burn <- SuggestBurn(0.1, object)
    prediction_errors <- bsts.prediction.errors(object, burn = burn)$in.sample
    prediction_sse <- sum(colMeans(prediction_errors)^2)
    original_series <- as.numeric(object$original.series)
    dy <- diff(original_series)
    prediction_sst <- var(dy) * (length(dy) - 1)
    relative_gof = 1 - prediction_sse / prediction_sst
  }
  return(relative_gof)
}

tsbacktest.bsts.spec <- function(object, start = floor(NROW(object$target$y_orig)/2), end = NROW(object$target$y_orig), h = 1, alpha = NULL, 
                                 cores = 1, data_name = "y", save_output = FALSE, save_dir = "~/tmp/", trace = FALSE, n_iter = 5000, ...)
{
  if (save_output) {
    if (is.null(save_dir)) {
      stop("save.dir cannot be NULL when save.output is TRUE")
    }
    if (!dir.exists(save_dir)) {
      stop("save.dir does not exit. Create first and then resubmit")
    }
  }
  transform <- object$transform
  lambda <- transform$lambda
  frequency <- object$target$frequency
  data <- xts(object$target$y_orig, object$target$index)
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
  i <- 1
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  if (trace) {
    iterations <- length(seqdates)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
  b <- foreach(i = 1:length(seqdates), .packages = c("tsmethods","tsaux","tsforeign","xts","data.table","bsts"), .options.snow = opts, .combine = rbind) %dopar% {
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
    spec <- bsts_modelspec(ytrain, xreg = xreg_train, frequency = object$target$frequency, differences = object$target$differences, level = object$model$level, 
                           slope = object$model$slope, damped = object$model$damped, seasonal = object$model$seasonal, seasonal_frequency = object$model$seasonal_frequency, 
                           ar = object$model$ar, ar_max = object$model$ar_max, cycle = object$model$cycle, cycle_frequency = object$model$cycle_frequency, 
                           cycle_names = object$model$cycle_names, seasonal_type = object$model$seasonal_type, seasonal_harmonics = object$model$seasonal_harmonics, 
                           lambda = lambda, distribution = object$model$distribution)
    mod <- estimate(spec, n_iter = n_iter)
    p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(ytest))
    if (save_output) {
      saveRDS(mod, file = paste0(save_dir,"/model_", seqdates[i], ".rds"))
      saveRDS(p, file = paste0(save_dir,"/predict_", seqdates[i], ".rds"))
    }
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
  }
  stopCluster(cl)
  if (trace) {
    close(pb)
  }
  if (is.null(data_name)) data_name <- "y"
  actual <- NULL
  forecast <- NULL
  metrics <- b[,list(variable = data_name, MAPE = mape(actual, forecast), MSLRE = mslre(actual, forecast),
                     BIAS = bias(actual, forecast),
                     n = .N), by = "horizon"]
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
