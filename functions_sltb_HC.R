if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

library(Matrix)
library(MASS)

# ============================================================
# Shared Laplace engine
#   - half-Cauchy prior on tau
#   - diffuse Normal(0, variance) priors on beta and gamma
# ============================================================
laplace_beta_mixed_common <- function(data,
                                      mean_formula = ~ x1 * x2 + x3 + x4,
                                      prec_formula = ~ 1,
                                      group_var = "id",
                                      y_var = "y",
                                      tau_scale = 6,
                                      b = NULL,
                                      beta_var = 1e6,
                                      gamma_var = 1e6,
                                      control_optim = list(maxit = 2000, reltol = 1e-10),
                                      tol = 1e-6,
                                      max_iter = 1000,
                                      verbose = FALSE,
                                      ridge_rel = 1e-8) {

  fail_result <- function(msg) {
    if (verbose) message("laplace_beta_mixed_common failure: ", msg)
    list(
      log_marginal = NA_real_,
      log_Cb = NA_real_,
      b = b,
      error = msg,
      lap_full = list(stage = NA_character_, err = msg),
      lap_b = list(stage = NA_character_, err = msg),
      p = NA_integer_,
      q = NA_integer_,
      G = NA_integer_
    )
  }

  if (!is.data.frame(data)) return(fail_result("data must be a data.frame"))
  if (!(group_var %in% names(data))) return(fail_result(paste0("group_var '", group_var, "' not found")))
  if (!(y_var %in% names(data))) return(fail_result(paste0("y_var '", y_var, "' not found")))

  y <- data[[y_var]]
  if (!is.numeric(y)) return(fail_result("y must be numeric"))
  if (any(!is.finite(y), na.rm = TRUE)) return(fail_result("y contains non-finite values"))

  X <- tryCatch(
    model.matrix(mean_formula, data),
    error = function(e) return(fail_result(paste0("mean_formula error: ", e$message)))
  )
  if (is.list(X) && !is.matrix(X)) return(X)

  Z <- tryCatch(
    model.matrix(prec_formula, data),
    error = function(e) return(fail_result(paste0("prec_formula error: ", e$message)))
  )
  if (is.list(Z) && !is.matrix(Z)) return(Z)

  grp <- as.factor(data[[group_var]])
  grp_idx <- as.integer(grp)
  G <- nlevels(grp)
  N <- nrow(X)
  if (length(y) != N) return(fail_result("y and design matrix rows mismatch"))

  M <- sparseMatrix(i = seq_len(N), j = grp_idx, x = 1L, dims = c(N, G))

  p <- ncol(X)
  q <- ncol(Z)

  if (is.null(b)) {
    b <- (p + 1) / N
  }
  if (!is.finite(b) || b <= 0 || b > 1) {
    return(fail_result("b must be in (0,1]"))
  }

  unpack_params <- function(par) {
    i <- 1
    beta <- par[i:(i + p - 1)]; i <- i + p
    gamma <- par[i:(i + q - 1)]; i <- i + q
    u <- par[i:(i + G - 1)]
    list(beta = beta, gamma = gamma, u = u)
  }

  obs_quantities <- function(beta, gamma, u, y_local) {
    eta <- as.numeric(X %*% beta + M %*% u)
    mu <- plogis(eta)
    dotmu <- mu * (1 - mu)
    ddotmu <- dotmu * (1 - 2 * mu)

    delt <- as.numeric(Z %*% gamma)
    phi <- exp(delt)

    y_safe <- y_local
    a <- mu * phi
    bpar <- (1 - mu) * phi

    dig_a <- digamma(a)
    dig_b <- digamma(bpar)
    dig_phi <- digamma(phi)

    tri_a <- trigamma(a)
    tri_b <- trigamma(bpar)
    tri_phi <- trigamma(phi)

    S <- -dig_a + dig_b + log(y_safe / (1 - y_safe))
    A <- phi * S
    C <- -phi^2 * (tri_a + tri_b)

    Bphi <- dig_phi - mu * dig_a - (1 - mu) * dig_b +
      mu * log(y_safe) + (1 - mu) * log1p(-y_safe)

    D_phi2 <- tri_phi - mu^2 * tri_a - (1 - mu)^2 * tri_b
    E <- S + phi * (-mu * tri_a + (1 - mu) * tri_b)
    W <- -(C * (dotmu^2) + A * ddotmu)

    list(mu = mu, phi = phi, A = A, Bphi = Bphi, C = C,
         D_phi2 = D_phi2, E = E, W = W, dotmu = dotmu)
  }

  logpost_and_grad_scaled <- function(par, tau, y_local, ll_scale = 1) {
    up <- unpack_params(par)
    beta <- up$beta
    gamma <- up$gamma
    u <- up$u

    oq <- obs_quantities(beta, gamma, u, y_local)
    mu <- oq$mu
    phi <- oq$phi
    A <- oq$A
    Bphi <- oq$Bphi
    dotmu <- oq$dotmu
    y_safe <- y_local

    a <- mu * phi
    bpar <- (1 - mu) * phi

    ll_terms <- sum(
      lgamma(phi) - lgamma(a) - lgamma(bpar) +
        (a - 1) * log(y_safe) +
        (bpar - 1) * log1p(-y_safe)
    )

    if (!is.finite(ll_terms)) {
      return(list(logpost = -Inf, grad = rep(NA_real_, p + q + G), oq = oq))
    }

    lp_u <- -0.5 * sum(u^2) / tau - (G / 2) * log(tau)
    lp_beta <- -0.5 * sum(beta^2) / beta_var - 0.5 * p * log(2 * pi * beta_var)
    lp_gamma <- -0.5 * sum(gamma^2) / gamma_var - 0.5 * q * log(2 * pi * gamma_var)

    logpost <- as.numeric(ll_scale * ll_terms + lp_u + lp_beta + lp_gamma)

    g_beta <- as.numeric(t(X) %*% (ll_scale * (A * dotmu))) - beta / beta_var
    g_gamma <- as.numeric(t(Z) %*% (ll_scale * (phi * Bphi))) - gamma / gamma_var

    obs_grad_u <- ll_scale * (A * dotmu)
    g_u <- numeric(G)
    for (gidx in seq_len(G)) {
      g_u[gidx] <- sum(obs_grad_u[grp_idx == gidx]) - u[gidx] / tau
    }

    grad <- c(g_beta, g_gamma, g_u)
    list(logpost = logpost, grad = grad, oq = oq)
  }

  logpi_tau_halfcauchy <- function(tau, s_tau) {
    if (tau <= 0) return(-Inf)
    -0.5 * log(tau) - log1p(tau / (s_tau^2))
  }

  update_tau_map_numeric <- function(u_vec, s_tau, tau_lower = 1e-12, tau_upper = 1e8) {
    f_neg <- function(tau_val) {
      if (!is.finite(tau_val) || tau_val <= 0) return(1e100)
      lp_tau <- logpi_tau_halfcauchy(tau_val, s_tau)
      lp_u <- -0.5 * sum(u_vec^2) / tau_val - (G / 2) * log(tau_val)
      -(lp_tau + lp_u)
    }

    Su2 <- 0.5 * sum(u_vec^2)
    if (is.finite(Su2) && Su2 > 0) {
      tau_upper <- max(tau_upper, Su2 * 10)
    }

    out <- tryCatch(
      optimize(f_neg, lower = tau_lower, upper = tau_upper, tol = 1e-8),
      error = function(e) {
        tryCatch(
          nlminb(start = max(1, mean(u_vec^2)), objective = f_neg),
          error = function(e2) list(par = 1)
        )
      }
    )

    tau_hat_local <- as.numeric(if (!is.null(out$minimum)) out$minimum else out$par)
    if (!is.finite(tau_hat_local) || tau_hat_local <= 0) tau_hat_local <- 1e-6
    tau_hat_local
  }

  find_map_scaled <- function(y_local, ll_scale = 1) {
    par0 <- rep(0, p + q + G)

    u0_by_grp <- tapply(y_local, grp, function(v) {
      mv <- mean(v)
      qlogis(mv)
    })
    if (length(u0_by_grp) == G && all(is.finite(u0_by_grp))) {
      par0[(p + q + 1):(p + q + G)] <- as.numeric(u0_by_grp) - mean(as.numeric(u0_by_grp))
    }

    par <- par0
    tau <- 1.0
    opt_ctrl <- modifyList(list(maxit = 2000, reltol = 1e-10), control_optim)
    iterloc <- 0
    convergedloc <- FALSE
    last_err <- NULL

    while (iterloc < max_iter) {
      iterloc <- iterloc + 1

      fn <- function(parv) {
        out <- logpost_and_grad_scaled(parv, tau, y_local, ll_scale)
        -out$logpost
      }
      gr <- function(parv) {
        out <- logpost_and_grad_scaled(parv, tau, y_local, ll_scale)
        -out$grad
      }

      opt_try <- tryCatch(
        optim(par, fn = fn, gr = gr, method = "BFGS", control = opt_ctrl),
        error = function(e) {
          last_err <<- e$message
          NULL
        }
      )

      if (is.null(opt_try) || any(!is.finite(opt_try$par))) {
        opt_try2 <- tryCatch(
          nlminb(start = par, objective = fn, gradient = gr, control = list(iter.max = 1000)),
          error = function(e) {
            last_err <<- e$message
            NULL
          }
        )
        if (is.null(opt_try2)) {
          return(list(par = par, tau = tau, iter = iterloc, converged = FALSE, err = last_err))
        } else {
          par_new <- opt_try2$par
        }
      } else {
        par_new <- opt_try$par
      }

      up <- unpack_params(par_new)
      tau_new <- update_tau_map_numeric(up$u, s_tau = tau_scale)

      dpar <- max(abs(par_new - par))
      dtau <- abs(tau_new - tau)

      par <- par_new
      tau <- tau_new

      if (dpar < tol && dtau < tol) {
        convergedloc <- TRUE
        break
      }
    }

    list(par = par, tau = tau, iter = iterloc, converged = convergedloc, err = last_err)
  }

  compute_laplace_scaled <- function(y_local, ll_scale = 1) {
    stage <- "start"

    tryCatch({
      stage <- "MAP"
      mp <- find_map_scaled(y_local, ll_scale = ll_scale)

      par_map <- mp$par
      tau_map <- mp$tau

      stage <- "unpack"
      up <- unpack_params(par_map)
      beta_hat <- up$beta
      gamma_hat <- up$gamma
      u_hat <- up$u
      tau_hat <- tau_map

      stage <- "obs_quantities"
      oq <- obs_quantities(beta_hat, gamma_hat, u_hat, y_local)

      W_vec <- oq$W * ll_scale
      dotmu <- oq$dotmu
      E_vec <- oq$E * ll_scale
      Bphi_vec <- oq$Bphi * ll_scale
      D_phi2_vec <- oq$D_phi2 * ll_scale
      phi <- oq$phi

      stage <- "hessian pieces"
      D_W <- Diagonal(x = as.numeric(W_vec))
      D_E <- Diagonal(x = as.numeric(dotmu * phi * E_vec))
      D_gamma_diag <- Diagonal(x = as.numeric(-phi^2 * D_phi2_vec - phi * Bphi_vec))

      H_bb <- as.matrix(t(X) %*% D_W %*% X) + diag(1 / beta_var, p)
      H_bg <- as.matrix(-t(X) %*% D_E %*% Z)
      H_bu <- as.matrix(t(X) %*% D_W %*% M)

      H_gb <- t(H_bg)
      H_gg <- as.matrix(t(Z) %*% D_gamma_diag %*% Z) + diag(1 / gamma_var, q)
      H_gu <- as.matrix(-t(Z) %*% D_E %*% M)

      H_ub <- t(H_bu)
      H_ug <- t(H_gu)

      sumW_by_group <- as.numeric(t(M) %*% as.numeric(W_vec))
      H_uu <- as.matrix(Diagonal(x = sumW_by_group) + Diagonal(x = rep(1 / tau_hat, G)))

      stage <- "tau hessian"
      c_tau <- tau_scale^2
      d2_logp_tau <- -sum(u_hat^2) / (tau_hat^3) +
        (G / 2) / (tau_hat^2) +
        1 / (2 * tau_hat^2) +
        1 / ((c_tau + tau_hat)^2)
      H_tautau <- -d2_logp_tau

      H_utau <- -(u_hat / (tau_hat^2))
      H_tauu <- matrix(H_utau, nrow = 1)

      stage <- "assemble H"
      top_left <- cbind(H_bb, H_bg, H_bu)
      mid_left <- cbind(H_gb, H_gg, H_gu)
      bot_left <- cbind(H_ub, H_ug, H_uu)
      H_no_tau <- rbind(top_left, mid_left, bot_left)

      ntheta <- p + q + G
      H_full <- matrix(0, nrow = ntheta + 1, ncol = ntheta + 1)
      H_full[1:ntheta, 1:ntheta] <- H_no_tau

      rows_u <- (p + q + 1):(p + q + G)
      H_full[rows_u, ntheta + 1] <- H_utau
      H_full[ntheta + 1, rows_u] <- H_tauu
      H_full[ntheta + 1, ntheta + 1] <- H_tautau

      stage <- "regularize H"
      H_sym <- (H_full + t(H_full)) / 2
      ridge <- ridge_rel * max(1, mean(abs(diag(H_sym))))
      H_sym <- H_sym + diag(ridge, nrow(H_sym))

      stage <- "invert H"
      Sigma_full <- NULL
      eig <- NULL

      chol_ok <- TRUE
      chol_try <- tryCatch({
        R <- chol(H_sym)
        Sigma_full <- chol2inv(R)
        NULL
      }, error = function(e) e$message)

      if (!is.null(chol_try)) chol_ok <- FALSE

      if (!chol_ok) {
        eig_try <- eigen(H_sym, symmetric = TRUE)
        eig_vals <- eig_try$values
        eig_vals[eig_vals < 1e-6] <- 1e-6
        eig <- list(values = eig_vals, vectors = eig_try$vectors)
        H_reg <- eig$vectors %*% diag(eig_vals) %*% t(eig$vectors)
        Sigma_full <- tryCatch(solve(H_reg), error = function(e) MASS::ginv(H_reg))
      } else {
        eig <- eigen(H_sym, symmetric = TRUE)
        eig$values[eig$values < 1e-12] <- 1e-12
        H_reg <- H_sym
      }

      stage <- "Laplace integral"
      lp_out_scaled <- logpost_and_grad_scaled(par_map, tau_map, y_local, ll_scale = ll_scale)
      logpost_scaled_no_tau <- lp_out_scaled$logpost
      logpi_tau_map <- logpi_tau_halfcauchy(tau_hat, tau_scale)
      log_numerical_at_map <- as.numeric(logpost_scaled_no_tau + logpi_tau_map)

      logdet_H <- sum(log(eig$values))
      k_full <- ntheta + 1
      log_integral <- as.numeric(log_numerical_at_map + (k_full / 2) * log(2 * pi) - 0.5 * logdet_H)

      list(
        log_integral = log_integral,
        par_map = par_map,
        tau_map = tau_map,
        H_reg = H_reg,
        Sigma_full = Sigma_full,
        eig = eig,
        stage = "ok",
        err = NULL
      )
    }, error = function(e) {
      if (verbose) message("compute_laplace_scaled internal error at stage ", stage, ": ", e$message)
      list(
        log_integral = NA_real_,
        par_map = NULL,
        tau_map = NULL,
        H_reg = NULL,
        Sigma_full = NULL,
        eig = NULL,
        stage = stage,
        err = e$message
      )
    })
  }

  lap_full <- compute_laplace_scaled(y, ll_scale = 1)
  log_marginal <- lap_full$log_integral

  lap_b <- compute_laplace_scaled(y, ll_scale = b)
  log_Cb <- lap_b$log_integral

  res <- list(
    log_marginal = log_marginal,
    log_Cb = log_Cb,
    b = b,
    lap_full = lap_full,
    lap_b = lap_b,
    p = p, q = q, G = G
  )
  class(res) <- "laplace_beta_mixed_common"
  res
}

# ============================================================
# Model formulas
# ============================================================
formula_list <- list(
  null   = as.formula("~ 1"),

  x1     = as.formula("~ x1"),
  x2     = as.formula("~ x2"),
  x3     = as.formula("~ x3"),
  x4     = as.formula("~ x4"),

  main   = as.formula("~ x1 + x2 + x3 + x4"),

  int12_3 = as.formula("~ x1 * x2 + x3"),
  int12_4 = as.formula("~ x1 * x2 + x4"),
  int13_2 = as.formula("~ x1 * x3 + x2"),
  int13_4 = as.formula("~ x1 * x3 + x4"),
  int14_2 = as.formula("~ x1 * x4 + x2"),
  int14_3 = as.formula("~ x1 * x4 + x3"),
  int23_1 = as.formula("~ x2 * x3 + x1"),
  int23_4 = as.formula("~ x2 * x3 + x4"),
  int24_1 = as.formula("~ x2 * x4 + x1"),
  int24_3 = as.formula("~ x2 * x4 + x3"),
  int34_1 = as.formula("~ x3 * x4 + x1"),
  int34_2 = as.formula("~ x3 * x4 + x2"),

  int12   = as.formula("~ x1 * x2 + x3 + x4"),
  int13   = as.formula("~ x1 * x3 + x2 + x4"),
  int14   = as.formula("~ x1 * x4 + x2 + x3"),
  int23   = as.formula("~ x2 * x3 + x1 + x4"),
  int24   = as.formula("~ x2 * x4 + x1 + x3"),
  int34   = as.formula("~ x3 * x4 + x1 + x2"),

  int123_4 = as.formula("~ x1 * x2 * x3 + x4"),
  int124_3 = as.formula("~ x1 * x2 * x4 + x3"),
  int134_2 = as.formula("~ x1 * x3 * x4 + x2"),
  int234_1 = as.formula("~ x2 * x3 * x4 + x1"),

  full   = as.formula("~ x1 * x2 * x3 * x4")
)

shared_args <- list(
  prec_formula = ~ 1,
  group_var = "id",
  y_var = "y",
  tau_scale = 6,
  beta_var = 1e6,
  gamma_var = 1e6,
  tol = 1e-8,
  max_iter = 500,
  verbose = FALSE
)
