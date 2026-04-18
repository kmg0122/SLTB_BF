if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

library(Matrix)
library(MASS)

# ============================================================
# Fixed SLT / SLTB constants
# ============================================================
s_fixed <- 1 + 10^(-8.5)
l_fixed <- 1e-9

safe_mu <- function(x) pmin(pmax(x, 1e-8), 1 - 1e-8)
safe_phi <- function(x) pmax(x, 1e-8)
safe_z <- function(z) pmin(pmax(z, 1e-12), 1 - 1e-12)

# ============================================================
# Failure helper
# ============================================================
fail_result <- function(msg, b = NA_real_, verbose = FALSE) {
  if (verbose) message("laplace_sltb_mixed_common_fast failure: ", msg)
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

# ============================================================
# SLTB helper functions
# ============================================================
sltb_logC <- function(mu, phi, s = s_fixed, l = l_fixed) {
  mu <- safe_mu(mu)
  phi <- safe_phi(phi)
  a <- mu * phi
  b <- (1 - mu) * phi

  zL <- l
  zU <- l + 1 / s

  C <- pbeta(zU, shape1 = a, shape2 = b) - pbeta(zL, shape1 = a, shape2 = b)
  log(pmax(C, .Machine$double.xmin))
}

sltb_dlogC <- function(mu, phi, s = s_fixed, l = l_fixed,
                       eps_mu = 1e-6, eps_phi = 1e-6) {
  mu <- safe_mu(mu)
  phi <- safe_phi(phi)

  mu_p <- safe_mu(mu + eps_mu)
  mu_m <- safe_mu(mu - eps_mu)
  phi_p <- safe_phi(phi + eps_phi)
  phi_m <- safe_phi(phi - eps_phi)

  dmu <- (sltb_logC(mu_p, phi, s, l) - sltb_logC(mu_m, phi, s, l)) / (mu_p - mu_m)
  dphi <- (sltb_logC(mu, phi_p, s, l) - sltb_logC(mu, phi_m, s, l)) / (phi_p - phi_m)

  c(dmu = as.numeric(dmu), dphi = as.numeric(dphi))
}

# ============================================================
# Simulation
# ============================================================
simulate_sltb_data <- function(n, G = max(5, round(n / 10)), seed = NULL,
                               beta0 = -0.5,
                               beta1 = 0.8, beta2 = -0.6, beta3 = 0.5, beta4 = 0,
                               beta12 = 0.4, beta13 = 0,
                               tau_sd = sqrt(0.20), phi = 10) {
  if (!is.null(seed)) set.seed(seed)

  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- 0.5 * x1 + rnorm(n, mean = 0, sd = 1)
  x3 <- runif(n, -1, 1)
  x4 <- rnorm(n, mean = 2, sd = 0.5)

  x1_x2 <- x1 * x2
  x1_x3 <- x1 * x3

  id <- sample(1:G, size = n, replace = TRUE)

  alpha_vec <- rnorm(G, mean = 0, sd = tau_sd)
  names(alpha_vec) <- 1:G
  u_obs <- alpha_vec[as.character(id)]

  eta <- beta0 +
    beta1 * x1 +
    beta2 * x2 +
    beta3 * x3 +
    beta4 * x4 +
    beta12 * x1_x2 +
    beta13 * x1_x3 +
    u_obs

  mu <- plogis(eta)
  a <- mu * phi
  b <- (1 - mu) * phi

  y_cont <- rbeta(n = n, shape1 = a, shape2 = b)
  y <- ifelse(y_cont < 0.01, 0,
              ifelse(y_cont > 0.99, 1, y_cont))

  df <- data.frame(
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x1_x2 = x1_x2,
    x1_x3 = x1_x3,
    id = id,
    eta = eta,
    mu = mu,
    u = u_obs
  )

  alpha_df <- data.frame(
    id = 1:G,
    alpha_true = as.numeric(alpha_vec)
  )

  list(
    data = df,
    alpha_true = alpha_df,
    alpha_vec = alpha_vec
  )
}

# ============================================================
# Shared Laplace engine
# Priors:
#   beta  : Normal(0, beta_var)
#   gamma : Normal(0, gamma_var)
#   u | tau : Normal(0, tau I)
#   tau   : half-Cauchy(scale = tau_scale)
# ============================================================
laplace_sltb_mixed_common_fast <- function(data,
                                           mean_formula = ~ x1 * x2 + x3 + x4,
                                           prec_formula = ~ 1,
                                           group_var = "id",
                                           y_var = "y",
                                           tau_scale = 6,
                                           b = 1,
                                           beta_var = 1e6,
                                           gamma_var = 1e6,
                                           control_optim = list(maxit = 2000, reltol = 1e-10),
                                           tol = 1e-6,
                                           max_iter = 1000,
                                           verbose = FALSE,
                                           ridge_rel = 1e-8) {
  if (!is.data.frame(data)) return(fail_result("data must be a data.frame", b = b, verbose = verbose))
  if (!(group_var %in% names(data))) return(fail_result(paste0("group_var '", group_var, "' not found"), b = b, verbose = verbose))
  if (!(y_var %in% names(data))) return(fail_result(paste0("y_var '", y_var, "' not found"), b = b, verbose = verbose))

  y <- data[[y_var]]
  if (!is.numeric(y)) return(fail_result("y must be numeric", b = b, verbose = verbose))
  if (any(!is.finite(y), na.rm = TRUE)) return(fail_result("y contains non-finite values", b = b, verbose = verbose))
  if (any(y < 0 | y > 1, na.rm = TRUE)) return(fail_result("y must lie in [0, 1]", b = b, verbose = verbose))

  X <- tryCatch(
    model.matrix(mean_formula, data),
    error = function(e) return(fail_result(paste0("mean_formula error: ", e$message), b = b, verbose = verbose))
  )
  if (is.list(X) && !is.matrix(X)) return(X)

  Z <- tryCatch(
    model.matrix(prec_formula, data),
    error = function(e) return(fail_result(paste0("prec_formula error: ", e$message), b = b, verbose = verbose))
  )
  if (is.list(Z) && !is.matrix(Z)) return(Z)

  grp <- as.factor(data[[group_var]])
  grp_idx <- as.integer(grp)
  G <- nlevels(grp)
  N <- nrow(X)
  if (length(y) != N) return(fail_result("y and design matrix rows mismatch", b = b, verbose = verbose))

  M <- sparseMatrix(i = seq_len(N), j = grp_idx, x = 1L, dims = c(N, G))

  p <- ncol(X)
  q <- ncol(Z)

  if (is.null(b)) b <- (p + 1) / N
  if (!is.finite(b) || b <= 0 || b > 1) {
    return(fail_result("b must be in (0,1]", b = b, verbose = verbose))
  }

  unpack_params <- function(par) {
    i <- 1
    beta <- par[i:(i + p - 1)]; i <- i + p
    gamma <- par[i:(i + q - 1)]; i <- i + q
    u <- par[i:(i + G - 1)]
    list(beta = beta, gamma = gamma, u = u)
  }

  logpi_tau_halfcauchy <- function(tau, s_tau) {
    if (tau <= 0) return(-Inf)
    -0.5 * log(tau) - log1p(tau / (s_tau^2))
  }

  obs_quantities <- function(beta, gamma, u, y_local) {
    eta <- as.numeric(X %*% beta + M %*% u)
    mu <- plogis(eta)
    dotmu <- mu * (1 - mu)
    ddotmu <- dotmu * (1 - 2 * mu)

    delt <- as.numeric(Z %*% gamma)
    phi <- exp(delt)

    z <- safe_z(y_local / s_fixed + l_fixed)
    a <- mu * phi
    bpar <- (1 - mu) * phi

    dig_a <- digamma(a)
    dig_b <- digamma(bpar)
    dig_phi <- digamma(phi)

    tri_a <- trigamma(a)
    tri_b <- trigamma(bpar)
    tri_phi <- trigamma(phi)

    S_beta <- -dig_a + dig_b + log(z) - log1p(-z)
    A_beta <- phi * S_beta
    C_beta <- -phi^2 * (tri_a + tri_b)

    Bphi_beta <- dig_phi - mu * dig_a - (1 - mu) * dig_b +
      mu * log(z) + (1 - mu) * log1p(-z)

    D_phi2_beta <- tri_phi - mu^2 * tri_a - (1 - mu)^2 * tri_b
    E_beta <- S_beta + phi * (-mu * tri_a + (1 - mu) * tri_b)
    W_beta <- -(C_beta * (dotmu^2) + A_beta * ddotmu)

    dC <- t(vapply(
      seq_along(mu),
      function(i) sltb_dlogC(mu[i], phi[i], s = s_fixed, l = l_fixed),
      numeric(2)
    ))
    colnames(dC) <- c("dmu", "dphi")

    score_mu <- A_beta - dC[, "dmu"]
    score_phi <- Bphi_beta - dC[, "dphi"]

    list(
      mu = mu,
      phi = phi,
      z = z,
      score_mu = score_mu,
      score_phi = score_phi,
      dotmu = dotmu,
      A_beta = A_beta,
      Bphi_beta = Bphi_beta,
      C_beta = C_beta,
      D_phi2_beta = D_phi2_beta,
      E_beta = E_beta,
      W_beta = W_beta
    )
  }

  logpost_and_grad_scaled <- function(par, tau, y_local, ll_scale = 1) {
    up <- unpack_params(par)
    beta <- up$beta
    gamma <- up$gamma
    u <- up$u

    oq <- obs_quantities(beta, gamma, u, y_local)
    mu <- oq$mu
    phi <- oq$phi
    z <- oq$z
    score_mu <- oq$score_mu
    score_phi <- oq$score_phi
    dotmu <- oq$dotmu

    a <- mu * phi
    bpar <- (1 - mu) * phi

    logC <- mapply(
      FUN = sltb_logC,
      mu = mu,
      phi = phi,
      MoreArgs = list(s = s_fixed, l = l_fixed)
    )

    ll_terms <- sum(
      lgamma(phi) - lgamma(a) - lgamma(bpar) +
        (a - 1) * log(z) +
        (bpar - 1) * log1p(-z) -
        logC
    ) - length(y_local) * log(s_fixed)

    if (!is.finite(ll_terms)) {
      return(list(logpost = -Inf, grad = rep(NA_real_, p + q + G), oq = oq))
    }

    lp_beta <- -0.5 * sum(beta^2) / beta_var - 0.5 * p * log(2 * pi * beta_var)
    lp_gamma <- -0.5 * sum(gamma^2) / gamma_var - 0.5 * q * log(2 * pi * gamma_var)
    lp_u <- -0.5 * sum(u^2) / tau - (G / 2) * log(tau)

    logpost <- as.numeric(ll_scale * ll_terms + lp_beta + lp_gamma + lp_u)

    g_beta <- as.numeric(t(X) %*% (ll_scale * (score_mu * dotmu))) - beta / beta_var
    g_gamma <- as.numeric(t(Z) %*% (ll_scale * (score_phi * phi))) - gamma / gamma_var

    obs_grad_u <- ll_scale * (score_mu * dotmu)
    g_u <- numeric(G)
    for (gidx in seq_len(G)) {
      g_u[gidx] <- sum(obs_grad_u[grp_idx == gidx]) - u[gidx] / tau
    }

    grad <- c(g_beta, g_gamma, g_u)
    list(logpost = logpost, grad = grad, oq = oq)
  }

  update_tau_map_numeric <- function(u_vec, s_tau, tau_lower = 1e-12, tau_upper = 1e8) {
    f_neg <- function(tau_val) {
      if (!is.finite(tau_val) || tau_val <= 0) return(1e100)
      lp_tau <- logpi_tau_halfcauchy(tau_val, s_tau)
      lp_u <- -0.5 * sum(u_vec^2) / tau_val - (length(u_vec) / 2) * log(tau_val)
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
      qlogis(safe_mu(mv))
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

    if (!convergedloc && verbose) {
      warning("Alternation (ll_scale = ", ll_scale, ") did not fully converge within max_iter. last_err = ",
              ifelse(is.null(last_err), "none", last_err))
    }

    list(par = par, tau = tau, iter = iterloc, converged = convergedloc, err = last_err)
  }

  compute_laplace_scaled <- function(y_local, ll_scale = 1) {
    stage <- "start"

    tryCatch({
      stage <- "MAP"
      mp <- find_map_scaled(y_local, ll_scale = ll_scale)
      if (!mp$converged && verbose) {
        warning("MAP did not converge for ll_scale = ", ll_scale, "; proceeding anyway.")
      }

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

      W_vec <- oq$W_beta * ll_scale
      dotmu <- oq$dotmu
      E_vec <- oq$E_beta * ll_scale
      Bphi_vec <- oq$Bphi_beta * ll_scale
      D_phi2_vec <- oq$D_phi2_beta * ll_scale
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
      H_reg <- NULL
      eig <- NULL

      chol_ok <- TRUE
      chol_try <- tryCatch({
        R <- chol(H_sym)
        Sigma_full <- chol2inv(R)
        H_reg <- H_sym
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
      }

      if (is.null(H_reg)) stop("H_reg is NULL")
      if (is.null(Sigma_full)) stop("Sigma_full is NULL")

      stage <- "Laplace integral"
      lp_out_scaled <- logpost_and_grad_scaled(par_map, tau_map, y_local, ll_scale = ll_scale)
      logpost_scaled_no_tau <- lp_out_scaled$logpost
      logpi_tau_map <- logpi_tau_halfcauchy(tau_map, tau_scale)
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
  class(res) <- "laplace_sltb_mixed_common_fast"
  res
}

# ============================================================
# Model formulas
# ============================================================
formula_list <- list(
  null       = as.formula("~ 1"),
  x1         = as.formula("~ x1"),
  x2         = as.formula("~ x2"),
  x3         = as.formula("~ x3"),
  x4         = as.formula("~ x4"),
  main       = as.formula("~ x1 + x2 + x3 + x4"),
  int12_3    = as.formula("~ x1 * x2 + x3"),
  int12_4    = as.formula("~ x1 * x2 + x4"),
  int13_2    = as.formula("~ x1 * x3 + x2"),
  int13_4    = as.formula("~ x1 * x3 + x4"),
  int14_2    = as.formula("~ x1 * x4 + x2"),
  int14_3    = as.formula("~ x1 * x4 + x3"),
  int23_1    = as.formula("~ x2 * x3 + x1"),
  int23_4    = as.formula("~ x2 * x3 + x4"),
  int24_1    = as.formula("~ x2 * x4 + x1"),
  int24_3    = as.formula("~ x2 * x4 + x3"),
  int34_1    = as.formula("~ x3 * x4 + x1"),
  int34_2    = as.formula("~ x3 * x4 + x2"),
  int12      = as.formula("~ x1 * x2 + x3 + x4"),
  int13      = as.formula("~ x1 * x3 + x2 + x4"),
  int14      = as.formula("~ x1 * x4 + x2 + x3"),
  int23      = as.formula("~ x2 * x3 + x1 + x4"),
  int24      = as.formula("~ x2 * x4 + x1 + x3"),
  int34      = as.formula("~ x3 * x4 + x1 + x2"),
  int123_4   = as.formula("~ x1 * x2 * x3 + x4"),
  int124_3   = as.formula("~ x1 * x2 * x4 + x3"),
  int134_2   = as.formula("~ x1 * x3 * x4 + x2"),
  int234_1   = as.formula("~ x2 * x3 * x4 + x1"),
  full       = as.formula("~ x1 * x2 * x3 * x4")
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
