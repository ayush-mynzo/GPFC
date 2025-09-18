#####################
## Helper Function ##
#####################

##-----------------------------------------------------
#######################################################
## Build wrapper function for model parameterization ##
#######################################################

mdlParameterization_wrapper <- function(mdl_dat, m, grps = NULL) {

  # Ensure required columns exist
  stopifnot(all(c("site_id", "age", "agc", "binomial") %in% names(mdl_dat)))

  # Default groups: all binomials present
  if (is.null(grps)) {
    grps <- sort(unique(mdl_dat$binomial))
  }

  # Chapman–Richards curve for saemix
  cr_curve <- function(psi, id, x){
    t <- x[, 1]
    A <- psi[id, 1]
    k <- psi[id, 2]
    fpred <- A * (1 - exp(-k * t))^m
    return(fpred)
  }

  # Helper: fit CR for a given data.frame
  fit_cr <- function(df) {
    # Starting values (Fekedulegn et al., 1999 inspired)
    a_init <- max(df$agc, na.rm = TRUE)
    k_init <- ((max(df$agc, na.rm = TRUE) - min(df$agc, na.rm = TRUE)) /
               (max(df$age, na.rm = TRUE) - min(df$age, na.rm = TRUE))) /
               max(df$age, na.rm = TRUE)

    smx_df <<- df

    saemix_dat <- saemixData(
      name.data       = smx_df,
      name.group      = "site_id",
      name.predictors = "age",
      name.response   = "agc"
    )

    saemix_model <- saemixModel(
      model = cr_curve,
      psi0  = c(A = a_init, k = k_init)
    )

    saemix_options <- list(map = TRUE, fim = TRUE, ll.is = FALSE,
                           displayProgress = FALSE, print = FALSE)

    cr_mdl <- saemix(saemix_model, saemix_dat, saemix_options)

    data.frame(
      a   = cr_mdl@results@fixed.effects[1],
      k   = cr_mdl@results@fixed.effects[2],
      a_se = cr_mdl@results@se.fixed[1],
      k_se = cr_mdl@results@se.fixed[2]
    )
  }

  results_df <- data.frame()

  set.seed(19890402)
  seeds <- round(runif(25, 1, 1000))

  for (j in 1:25) {
    seed <- seeds[j]

    for (i in grps) {
      set.seed(seed)

      # --------- SUBSET BY BINOMIAL (requested change) ----------
      sub_dat <- dplyr::filter(mdl_dat, binomial == i)

      # If not enough data, skip safely
      if (nrow(sub_dat) < 3 || length(unique(sub_dat$site_id)) < 2) {
        df_line <- data.frame(group = i, seed = seed,
                              a = NA, a_se = NA, k = NA, k_se = NA,
                              rmse = NA, nrmse_avg = NA, nrmse_sd = NA)
        results_df <- dplyr::bind_rows(results_df, df_line)
        next
      }

      sites   <- unique(sub_dat$site_id)
      sites_n <- max(1, round(length(sites) * 0.15))

      val_sites <- sample(sites, sites_n)
      val_df <- dplyr::filter(sub_dat, site_id %in% val_sites)
      trn_df <- dplyr::filter(sub_dat, !(site_id %in% val_sites))

      agc_avg <- mean(sub_dat$agc, na.rm = TRUE)
      agc_sd  <- stats::sd(sub_dat$agc, na.rm = TRUE)

      res <- tryCatch(fit_cr(trn_df),
                      error = function(e) data.frame(a=NA,k=NA,a_se=NA,k_se=NA))

      if (!is.na(res$a)) {
        # FIX: correct parentheses in CR prediction
        h <- dplyr::select(val_df, age, agc_obs = agc) %>%
             dplyr::mutate(agc_prd = res$a * (1 - exp(-res$k * age))^m)

        rmse_val <- ModelMetrics::rmse(h$agc_obs, h$agc_prd)

        df_line <- data.frame(
          group = i, seed = seed,
          a = res$a, a_se = res$a_se, k = res$k, k_se = res$k_se,
          rmse = rmse_val,
          nrmse_avg = rmse_val / agc_avg,
          nrmse_sd  = rmse_val / agc_sd
        )
      } else {
        df_line <- data.frame(group = i, seed = seed,
                              a = NA, a_se = NA, k = NA, k_se = NA,
                              rmse = NA, nrmse_avg = NA, nrmse_sd = NA)
      }

      results_df <- dplyr::bind_rows(results_df, df_line)
    }
  }

  # Summarize RMSE values
  rmses <- results_df %>%
    dplyr::filter(!is.na(a)) %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      n = dplyr::n(),
      rmse = mean(rmse, na.rm = TRUE),
      nrmse_avg = mean(nrmse_avg, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(nrmse_avg)

  #------------------------------------------------
  # Estimate parameters using the full dataset (by binomial)
  full_params <- data.frame()

  for (i in grps) {
    message(i)

    sub_dat <- dplyr::filter(mdl_dat, binomial == i)

    res_full <- tryCatch(fit_cr(sub_dat),
                         error = function(e) data.frame(a=NA,k=NA,a_se=NA,k_se=NA))

    res2 <- data.frame(
      grp  = i,
      a    = res_full$a,
      a_se = res_full$a_se,
      k    = res_full$k,
      k_se = res_full$k_se
    )

    full_params <- dplyr::bind_rows(full_params, res2)
  }

  df_full <- full_params %>%
    dplyr::mutate(m = m) %>%
    dplyr::select(grp, a, a_se, k, k_se, m) %>%
    dplyr::left_join(rmses, by = c("grp" = "group")) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4)))

  return(df_full)
}
