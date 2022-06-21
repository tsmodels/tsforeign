local_level_ss <- function(object)
{
    Z <- matrix(1, 1, 1)
    G <- matrix(1, 1, 1)
    Q <- matrix(as.numeric(NA), 1, 1)
    # Q2 is for consistency with dlm
    Q2 <- matrix(as.numeric(NA), 1, 1)
    R <- matrix(1, 1, 1)
    bsts_names <- "sigma.level"
    coef_names <- "level.sigma"
    coef_state_matrix <- "W"
    component_names <- c("Level")
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
}

local_linear_ss <- function(object)
{
    Z <- matrix(c(1,0), 1, 2)
    G <- matrix(c(1,1,0,1), 2, 2, byrow = TRUE)
    Q <- diag(as.numeric(NA), 2, 2)
    Q2 <- diag(as.numeric(NA), 2, 2)
    R <- diag(1, 2, 2)
    bsts_names <- c("sigma.trend.level","sigma.trend.slope")
    coef_names <- c("level.sigma","slope.sigma")
    coef_state_matrix <- c("W","W")
    component_names <- c("Level","Slope")
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
}

local_linear_damped_ss <- function(object)
{
    Z <- matrix(c(1,0,0), 1, 3)
    G <- matrix(c(1,1,0,0, NA, (1 - NA), 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
    Q <- diag(as.numeric(c(NA, NA)), 2, 2)
    Q2 <- diag(as.numeric(c(NA, NA, 0)), 3, 3)
    R <- rbind(diag(c(1,1), 2, 2), matrix(0, 1, 2))
    # Level, Slope, D, phi
    bsts_names <- c("trend.level.sd","trend.slope.sd","trend.slope.ar.coefficient","trend.slope.ar.coefficient")
    coef_names <- c("level.sigma","slope.sigma","slope.ar","slope.oneminusar")
    # for the second slope.ar we transform to 1-slope.ar
    coef_state_matrix <- c("W","W","G","G")
    component_names <- c("Level","Slope", "Slope.Mean")
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
}

seasonal_ss <- function(object, frequency)
{
    Z <- matrix(c(1, rep(0, frequency - 2)), nrow = 1)
    G <- matrix(0, frequency - 1, frequency - 1)
    G[row(G) == col(G) + 1] <- 1
    G[1, ] <- -1
    Q <- matrix(as.numeric(NA), 1, 1)
    Q2 <- diag(c(as.numeric(NA), rep(0, frequency - 2)), frequency - 1, frequency - 1)
    R <- matrix(c(1, rep(0, frequency - 2)), ncol = 1)
    bsts_names <- paste0("sigma.seasonal.", frequency)
    coef_names <- paste0("seasonal", frequency,".sigma")
    coef_state_matrix <- c("W")
    component_names <- c(paste0("Seasonal.0.",frequency), paste0("Seasonal.",1:(frequency - 2),".",frequency))
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
    
}

trigonometric_ss <- function(object, frequency, harmonics)
{
    if (is.integer(frequency)) {
        out <- dlmModTrig(s = frequency, q = harmonics, dW = as.numeric(NA))
    } else {
        out <- dlmModTrig(tau = frequency, q = harmonics, dW = as.numeric(NA))
    }
    Z <- out$FF
    G <- out$GG
    Q <- out$W
    Q2 <- out$W
    R <- diag(1, 2 * harmonics, 2 * harmonics)
    bsts_names <- paste0("trig.coefficient.sd.", frequency)
    coef_names <- rep(paste0("seasonal", frequency,".trig.sigma"), 2*harmonics)
    coef_state_matrix <- rep("W", 2*harmonics)
    component_names <- as.vector(sapply(1:harmonics, function(i) c(paste0("Trigonometric.", frequency,".",i,".",c(1,0)))))
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
    
}

cycle_ss <- function(object, frequency)
{
    if (is.integer(frequency)) {
        out <- dlmModTrig(s = frequency, q = 1, dW = as.numeric(NA))
    } else {
        out <- dlmModTrig(tau = frequency, q = 1, dW = as.numeric(NA))
    }
    Z <- out$FF
    G <- out$GG
    Q <- out$W
    Q2 <- out$W
    R <- diag(1, 2, 2)
    bsts_names <- paste0("trig.coefficient.sd.", frequency)
    coef_names <- rep(paste0("cycle",frequency,".sigma"),2)
    coef_state_matrix <- c("W","W")
    component_names <- paste0("Cyclical.", frequency,".","1",".",c(1,0))
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
    
}

ar_ss <- function(object, ar)
{
    out <- dlmModARMA(ar = rep(as.numeric(NA), ar))
    Z <- out$FF
    # need to transpose this to have the correct representation as bsts
    G <- t(out$GG)
    Q <- matrix(NA, 1, 1)
    Q2 <- diag(c(as.numeric(NA), rep(0, ar - 1)), ar, ar)
    if (ar > 1) {
        R <- matrix(c(1, rep(0, ar - 1)), ncol = 1)
    } else {
        R <- matrix(1, 1, 1)
    }
    bsts_names <- c(paste0("AR",ar,".coefficients"), paste0("AR",ar,".sigma"))
    coef_names <- c(paste0("ar",1:ar), paste0("ar",ar,".sigma"))
    coef_state_matrix <- c(rep("G", ar), "W")
    component_names <- paste0("AR.",1:ar)
    if (!is.null(object)) {
        Z <- cbind(object$Z, Z)
        G <- bdiag(object$G, G)
        Q <- bdiag(object$Q, Q)
        Q2 <- bdiag(object$Q2, Q2)
        R <- bdiag(object$R, R)
        bsts_names <- c(object$bsts_names, bsts_names)
        coef_names <- c(object$coef_names, coef_names)
        coef_state_matrix <- c(object$coef_state_matrix, coef_state_matrix)
        component_names <- c(object$component_names, component_names)
    }
    list(Z = Z, G = G, Q = Q, Q2 = Q2, R = R, bsts_names = bsts_names, coef_names = coef_names, coef_state_matrix = coef_state_matrix, component_names = component_names)
}

generate_ss_model.gaussian <- function(level = TRUE, slope = TRUE, damped = TRUE, seasonal = TRUE, seasonal_frequency = 4, seasonal_type = "regular", seasonal_harmonics = NULL, 
                              cycle = FALSE, cycle_frequency = 20, ar = FALSE, ar_max = 3, xreg = NULL)
{
    s_spec <- NULL
    if (level & !slope) {
        s_spec <- local_level_ss(s_spec)
    } else if (level & slope) {
        if (damped) {
            s_spec <- local_linear_damped_ss(s_spec)
        } else {
            s_spec <- local_linear_ss(s_spec)
        }
    }
    if (seasonal) {
        if (seasonal_type == "regular") {
            for (i in 1:length(seasonal_frequency)) {
                s_spec <- seasonal_ss(s_spec, seasonal_frequency[i])
            }
        } else {
            for (i in 1:length(seasonal_frequency)) {
                s_spec <- trigonometric_ss(s_spec, seasonal_frequency[i], seasonal_harmonics[i])
            }
        }
    }
    if (cycle) {
        for (i in 1:length(cycle_frequency)) {
            s_spec <- cycle_ss(s_spec, cycle_frequency[i])
        }
    }
    if (ar) {
        s_spec <- ar_ss(s_spec, ar_max)
    }
    if (!is.null(xreg)) {
        s_spec$bsts_names <- c(s_spec$bsts_names, "coefficients")
        s_spec$coef_names <- c(s_spec$coef_names, paste0("beta",1:ncol(xreg)))
    }
    return(s_spec)
}

#################################################
bdiag <- function(...)
{
    if (nargs() == 1) x <- as.list(...) else x <- list(...)
    n <- length(x)
    if (n == 0) return(NULL)
    zero_components <- sapply(x, function(y) length(y))
    if (any(zero_components == 0)) {
        exc <- which(zero_components == 0)
        x <- x[-exc]
        n <- length(x)
        if (n == 0) return(NULL)
    }
    x <- lapply(x, function(y) if (length(y)) as.matrix(y) else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2], (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}
#################################################
#' BSTS Extractor Functions
#'
#' @param object object of class \dQuote{bsts.estimate}.
#' @details The full and final states functions extract the smoothed states 
#' distribution from the BSTS estimated object, and the posterior function the
#' draws x parameters matrix.
#' @aliases bsts_full_states bsts_final_state bsts_posterior
#' @rdname bsts_extractors
#' @export
#'
bsts_full_states <- function(object)
{
    if (!inherits(object, "bsts.estimate")) stop("\nobject must inherit class bsts.estimate")
    return(object$model$full.state)
}

#' @rdname bsts_extractors
#' @export
#'
bsts_final_state <- function(object)
{
    if (!inherits(object, "bsts.estimate")) stop("\nobject must inherit class bsts.estimate")
    out <- object$model$final.state
    colnames(out) <- object$spec$state_space$component_names
    return(out)
}

#' @rdname bsts_extractors
#' @export
#'
bsts_posterior <- function(object)
{
    if (!inherits(object, "bsts.estimate")) stop("\nobject must inherit class bsts.estimate")
    out <- do.call(cbind, lapply(1:length(object$spec$state_space$bsts_names), function(i) object$model[[object$spec$state_space$bsts_names[i]]]))
    if (any(object$spec$state_space$bsts_names == "trend.slope.ar.coefficient")) {
        indx <- which(object$spec$state_space$bsts_names == "trend.slope.ar.coefficient")[1]
        out[, indx + 1] <- 1 - out[, indx + 1]
    }
    out <- cbind(object$model$sigma.obs, out)
    colnames(out) <- c("obs.sigma", unique(object$spec$state_space$coef_names))
    return(out)
}
