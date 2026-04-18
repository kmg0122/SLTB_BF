rm(list = ls())

source("functions_sltb_HC.R")

# ------------------------------------------------------------
# Read arguments
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_one_condition.R n tau_name tau_sd")
}

n <- as.integer(args[1])
tau_name <- args[2]
tau_sd <- as.numeric(args[3])

R <- 100
true_model_name <- "int12_3"

base_dir <- file.path(getwd(), "mc_100rep_SLTB_halfCauchy")
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

cond_dir <- file.path(base_dir, paste0("n", n), tau_name)
dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(122)

mc_summary <- list()
selection_list <- list()
posterior_long_list <- list()
mc_best_tables <- list()
mc_theta_tables <- list()

for (r in 1:R) {
  sim_seed <- 100000 + 1000 * n + 100 * match(tau_name, c("small", "medium", "big")) + r
  
  sim <- simulate_sltb_data(
    n = n,
    seed = sim_seed,
    tau_sd = tau_sd
  )
  
  df <- sim$data
  df$id <- as.factor(df$id)
  df$x1_x2 <- df$x1 * df$x2
  df$x1_x3 <- df$x1 * df$x3
  
  n_y0 <- sum(df$y == 0, na.rm = TRUE)
  n_y1 <- sum(df$y == 1, na.rm = TRUE)
  n_y01 <- n_y0 + n_y1
  
  write.csv(df,
            file.path(cond_dir, sprintf("sim_n%d_%s_rep%03d.csv", n, tau_name, r)),
            row.names = FALSE)
  write.csv(sim$alpha_true,
            file.path(cond_dir, sprintf("alpha_n%d_%s_rep%03d.csv", n, tau_name, r)),
            row.names = FALSE)
  
  res <- tryCatch(
    {
      bf_logm <- setNames(rep(NA_real_, length(formula_list)), names(formula_list))
      
      for (nm in names(formula_list)) {
        fit <- do.call(
          laplace_sltb_mixed_common_fast,
          c(list(data = df, mean_formula = formula_list[[nm]], b = 1), shared_args)
        )
        bf_logm[nm] <- fit$log_marginal
      }
      
      bf_valid <- is.finite(bf_logm)
      if (!any(bf_valid)) stop("No valid BF fits produced.")
      
      bf_best <- max(bf_logm[bf_valid])
      bf_post <- rep(NA_real_, length(bf_logm))
      bf_post[bf_valid] <- exp(bf_logm[bf_valid] - bf_best)
      bf_post[bf_valid] <- bf_post[bf_valid] / sum(bf_post[bf_valid])
      
      bf_table <- data.frame(
        model = names(bf_logm),
        log_marginal = bf_logm,
        post_prob = bf_post,
        stringsAsFactors = FALSE
      )
      bf_table <- bf_table[order(-bf_table$post_prob), ]
      
      best_name <- bf_table$model[1]
      best_formula <- formula_list[[best_name]]
      
      best_model <- laplace_sltb_mixed_common_fast(
        data = df,
        mean_formula = best_formula,
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
      
      X_best <- model.matrix(best_formula, data = df)
      Z_best <- model.matrix(~ 1, data = df)
      
      par_map <- best_model$lap_full$par_map
      tau_map <- best_model$lap_full$tau_map
      Sigma_full <- best_model$lap_full$Sigma_full
      
      names_par <- c(
        paste0("beta_", colnames(X_best)),
        paste0("gamma_", colnames(Z_best)),
        paste0("u_", seq_len(nlevels(df$id)))
      )
      
      sd_vec <- if (!is.null(Sigma_full)) {
        sqrt(pmax(0, diag(Sigma_full)))
      } else {
        rep(NA_real_, length(par_map) + 1)
      }
      
      theta_hat <- data.frame(
        name  = c(names_par, "tau"),
        mean  = as.numeric(c(par_map, tau_map)),
        sd    = as.numeric(sd_vec),
        stringsAsFactors = FALSE
      )
      theta_hat$lower <- theta_hat$mean - 1.96 * theta_hat$sd
      theta_hat$upper <- theta_hat$mean + 1.96 * theta_hat$sd
      
      list(
        bf_table = bf_table,
        theta_hat = theta_hat,
        best_name = best_name,
        best_formula = deparse(best_formula),
        log_marginal = best_model$log_marginal,
        tau_map = tau_map,
        error = NA_character_
      )
    },
    error = function(e) {
      list(
        bf_table = NULL,
        theta_hat = NULL,
        best_name = NA_character_,
        best_formula = NA_character_,
        log_marginal = NA_real_,
        tau_map = NA_real_,
        error = e$message
      )
    }
  )
  
  if (!is.null(res$bf_table)) {
    write.csv(res$bf_table,
              file.path(cond_dir, sprintf("bf_table_n%d_%s_rep%03d.csv", n, tau_name, r)),
              row.names = FALSE)
  }
  
  if (!is.null(res$theta_hat)) {
    write.csv(res$theta_hat,
              file.path(cond_dir, sprintf("theta_hat_n%d_%s_rep%03d.csv", n, tau_name, r)),
              row.names = FALSE)
  }
  
  chosen_model <- ifelse(is.null(res$best_name), NA_character_, res$best_name)
  is_correct <- as.integer(!is.na(chosen_model) && chosen_model == true_model_name)
  
  chosen_post_prob <- NA_real_
  if (!is.null(res$bf_table) && !is.na(chosen_model)) {
    chosen_post_prob <- res$bf_table$post_prob[match(chosen_model, res$bf_table$model)]
  }
  
  mc_summary[[r]] <- data.frame(
    rep = r,
    n = n,
    tau_label = tau_name,
    tau_sd = tau_sd,
    tau_var = tau_sd^2,
    seed = sim_seed,
    n_y0 = n_y0,
    n_y1 = n_y1,
    n_y01 = n_y01,
    best_model = chosen_model,
    best_formula = ifelse(is.null(res$best_formula), NA_character_, res$best_formula),
    log_marginal = ifelse(is.null(res$log_marginal), NA_real_, res$log_marginal),
    tau_map = ifelse(is.null(res$tau_map), NA_real_, res$tau_map),
    correct = is_correct,
    chosen_post_prob = chosen_post_prob,
    error = ifelse(is.null(res$error), NA_character_, res$error),
    stringsAsFactors = FALSE
  )
  
  selection_list[[r]] <- data.frame(
    rep = r,
    n = n,
    tau_label = tau_name,
    tau_sd = tau_sd,
    seed = sim_seed,
    n_y0 = n_y0,
    n_y1 = n_y1,
    n_y01 = n_y01,
    true_model = true_model_name,
    chosen_model = chosen_model,
    correct = is_correct,
    chosen_post_prob = chosen_post_prob,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(res$bf_table)) {
    tmp_post <- res$bf_table
    tmp_post$rep <- r
    tmp_post$n <- n
    tmp_post$tau_label <- tau_name
    tmp_post$tau_sd <- tau_sd
    tmp_post$seed <- sim_seed
    tmp_post$true_model <- true_model_name
    tmp_post$chosen_model <- chosen_model
    tmp_post$correct <- is_correct
    posterior_long_list[[r]] <- tmp_post
  }
  
  mc_best_tables[[as.character(r)]] <- res$bf_table
  mc_theta_tables[[as.character(r)]] <- res$theta_hat
}

mc_summary_df <- do.call(rbind, mc_summary)
selection_df <- do.call(rbind, selection_list)
posterior_long_df <- do.call(rbind, posterior_long_list)

write.csv(mc_summary_df, file.path(cond_dir, "mc_summary_condition.csv"), row.names = FALSE)
write.csv(selection_df, file.path(cond_dir, "selection_condition.csv"), row.names = FALSE)
write.csv(posterior_long_df, file.path(cond_dir, "posterior_long_condition.csv"), row.names = FALSE)

saveRDS(mc_best_tables, file.path(cond_dir, "mc_best_tables_condition.rds"))
saveRDS(mc_theta_tables, file.path(cond_dir, "mc_theta_tables_condition.rds"))

cat("Done for n =", n, "tau =", tau_name, "\n")
