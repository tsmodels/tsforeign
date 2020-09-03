d2levels_fitted = function(act, pred, d = 1, transform = NULL, frequency = 1)
{
  if (!is.null(transform)) {
      y_t <- transform$transform(act, transform$lambda)
      tmp <- cbind(y_t, pred)
      start <- min(which(!is.na(tmp[,2])))
      end <- nrow(tmp)
      tmp$value <- NA
      for (i in start:end) tmp[i,"value"] <- tail(diffinv(as.numeric(tmp[i,2]), differences = d, xi = as.numeric(tmp[(i - d):(i - 1), 1])), 1)
      tmp <- tmp[,"value"]
      tmp[1:d] <- y_t[1:d]
      tmp <- transform$inverse(tmp, transform$lambda)
    } else{
      tmp <- cbind(act, pred)
      start <- min(which(!is.na(tmp[,2])))
      end <- nrow(tmp)
      tmp$value <- NA
      for (i in start:end) tmp[i,"value"] <- tail(diffinv(as.numeric(tmp[i,2]), differences = d, xi = as.numeric(tmp[(i - d):(i - 1), 1])), 1)
      tmp <- tmp[,"value"]
      tmp[1:d] <- act[1:d]
  }
  return(tmp)
}


d2levels_fitted_matrix = function(act, pred, d = 1, transform = NULL, frequency = 1)
{
  if (!is.null(transform)) {
    y_t <- transform$transform(act, transform$lambda)
    if (d == 1) {
      y_t = matrix(as.numeric(y_t[1]), ncol = length(pred), nrow = 1)
    } else{
      y_t <- rbind(matrix(as.numeric(y_t[1]), ncol = length(pred), nrow = 1),
                  matrix(as.numeric(y_t[2]), ncol = length(pred), nrow = 1))
    }
    tmp <- t(diffinv(matrix(pred, nrow = 1), differences = d, xi = y_t))
    tmp <- tmp[,ncol(tmp), drop = FALSE]
    tmp <- transform$inverse(tmp, transform$lambda)
  } else {
    if (d == 1) {
      y_t <- matrix(as.numeric(act[1]), ncol = length(pred), nrow = 1)
    } else {
      y_t <- rbind(matrix(as.numeric(act[1]), ncol = length(pred), nrow = 1),
                  matrix(as.numeric(act[2]), ncol = length(pred), nrow = 1))
    }
    tmp = t(diffinv(matrix(pred, nrow = 1), differences = d, xi = y_t))
    tmp <- tmp[,ncol(tmp), drop = FALSE]
  }
  return(tmp)
}

d2levels_forecast = function(act, pred, d = 1, transform = NULL, frequency = 1)
{
  act2use <- tail(as.numeric(act), d)
  dp <- tail(diffinv(pred, differences = d, xi = act2use), length(pred))
  if (!is.null(transform)) {
    dp <- transform$inverse(dp, transform$lambda)
  }
  return(dp)
}

tsconvert.bsts.estimate <- function(object, to = "dlm", draw = "mean", burn = SuggestBurn(0.1, object$model), ...)
{
  to <- match.arg(to[1], "dlm")
  P <- bsts_posterior(object)
  if (!is.null(burn)) {
    if (burn > 0) {
      burn <- 1:burn
    } else {
      burn <- -(1:NROW(P))
    }
  }
  P <- P[-burn, , drop = FALSE]
  init_state <- bsts_full_states(object)[-burn,,1]
  Z <- object$spec$state_space$Z
  G <- object$spec$state_space$G
  Q <- object$spec$state_space$Q2
  if (draw != "mean") {
    draw <- as.integer(draw[1])
    if (!draw %in% c(1:nrow(P))) stop("\ndraw outside bounds of no.draws after burn")
    a0 <- init_state[draw,]
    P0 <- diag(1e7, ncol(init_state), ncol(init_state))
    if (any(object$spec$state_space$component_names == "Slope.Mean")) {
      ix <- which(object$spec$state_space$component_names == "Slope.Mean")
      P0[ix, ix] <- 1e-12
    }
    gpriors <- P[draw,object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "G")]]
    qpriors <- P[draw,object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "W")]]
    if (length(gpriors) > 0) {
      G[which(is.na(G))] <- gpriors
    }
    if (length(qpriors) > 0) {
      Q[which(is.na(Q))] <- qpriors
    }
    v <- object$model$sigma.obs[-burn][draw]
  } else {
    a0 <- colMeans(init_state)
    P0 <- diag(1e7, ncol(init_state), ncol(init_state))
    if (any(object$spec$state_space$component_names == "Slope.Mean")) {
      ix <- which(object$spec$state_space$component_names == "Slope.Mean")
      P0[ix, ix] <- 1e-12
    }
    gpriors <- colMeans(P[,object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "G")]])
    qpriors <- colMeans(P[,object$spec$state_space$coef_names[which(object$spec$state_space$coef_state_matrix == "W")]])
    if (length(gpriors) > 0) {
      G[which(is.na(G))] <- gpriors
    }
    if (length(qpriors) > 0) {
      Q[which(is.na(Q))] <- qpriors
    }
    v <- mean(object$model$sigma.obs[-burn])
  }
 
  if (!is.null(object$spec$xreg$xreg)) {
    X <- coredata(object$spec$xreg$xreg)
    nc <- ncol(X)
    Z <- cbind(Z, matrix(1, ncol = nc, nrow = 1))
    G <- bdiag(G, diag(1, nc, nc))
    Q <- bdiag(Q, diag(0, nc, nc))
    if (draw != "mean") {
      a0 <- c(a0, P[draw,paste0("beta",1:nc)])
    } else {
      a0 <- c(a0, colMeans(P[,paste0("beta",1:nc), drop = FALSE]))
    }
    P0 <- bdiag(P0, diag(1e-12, nc, nc))
    JFF <- cbind(matrix(0, ncol = ncol(object$spec$state_space$Z)), matrix(1:nc, nrow = 1))
    dlm_model <- dlm(FF = Z, GG = G, W = Q, V = v, m0 = a0, C0 = P0, JFF = JFF, X = X)
  } else {
    dlm_model <- dlm(FF = Z, GG = G, W = Q, V = v, m0 = a0, C0 = P0)
  }
  return(dlm_model)
}
