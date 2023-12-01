strat_folds <- function (y, nfolds = 5) {
  foldid <- integer(length(y))
  foldids <- seq_len(nfolds)
  for (ylev in unique(y)) {
    w <- which(y == ylev)
    foldid[w] <- sample(rep_len(foldids, length(w)))
    foldids <- rev(foldids)
  }
  foldid
}

sens_spec <- function (x, threshold = NA) {
  my_param_grid <- x$cv_param_grid
  colnames(my_param_grid) <- paste('cvparam', seq_len(ncol(my_param_grid)),
                                   sep = '_')

  sensspec <- x$cv_stat %>%
    map_dfr(function (pred_probs) {
      # Compute sensitivity & specificity
      class <- if (isTRUE(is.na(threshold))) {
        t(t(pred_probs) > x$adaptive_threshold)
      } else {
        pred_probs > threshold
      }

      cbind(my_param_grid,
            tp = colSums(x$y_01 * class, na.rm = TRUE),
            tn = colSums((1 - x$y_01) * (1 - class), na.rm = TRUE)) %>%
        mutate(sens = tp / sum(x$y_01),
               spec = tn / sum(1 - x$y_01))
    }) %>%
    group_by(across(starts_with('cvparam'))) %>%
    summarize(across(everything(), list(avg = mean, sd = sd)),
              .groups = 'drop')

  colnames(sensspec)[1:ncol(my_param_grid)] <- colnames(x$cv_param_grid)

  sensspec
}


#' @param fit A glmnet fit with slot `sens_spec`, or the `sens_spec`
#'   data frame.
#' @param spec The specificity threshold.
select_best <- function (x, spec = 1, se = FALSE) {
  ss <- if (is.data.frame(x)) {
    x
  } else {
    x$sens_spec
  }

  max_sens <- if (isTRUE(se)) {
    filter(ss, spec_avg - spec_sd >= spec - .Machine$double.eps)
  } else {
    filter(ss, spec_avg >= spec - .Machine$double.eps)
  }

  max_sens %<>%
    filter(sens_avg >= max(sens_avg) - .Machine$double.eps) %>%
    slice_head(n = 1)

  if ('glmnet' %in% class(x)) {
    x$lambda_spec <- max_sens$lambda
  } else if ('gbm' %in% class(x)) {
    x$ntree_spec <- max_sens$n_trees
  }

  if (is.data.frame(x)) {
    max_sens
  } else {
    x$spec_cv_stats <- max_sens
    x
  }
}

sens_to_nr <- function (x, sd = FALSE) {
  # nr <- x$spec_cv_stats$sens_avg * sum(vital_status == 'poor')
  sens_range <- if (isTRUE(sd)) {
    with(x$spec_cv_stats,
         pmax(0, pmin(1, sens_avg + c(-1, 1) * sens_sd)))
  } else {
    rep.int(x$spec_cv_stats$sens_avg, 2)
  }

  nr <- sens_range * sum(x$y_01)

  nr <- c(floor(nr[[1]]), ceiling(nr[[2]]))
  if (abs(nr[[1]] - nr[[2]]) < .Machine$double.eps) {
    sprintf('~%d', nr[[1]])
  } else {
    paste(nr, collapse = '--')
  }
}

sens_to_pct <- function (x) {
  scales::percent(x$spec_cv_stats$sens_avg)
}


cv_roc <- function (cv_fit) {
  lambda_ind <- which.min(abs(cv_fit$lambda - cv_fit$lambda_spec))

  rocs <- cv_fit$cv_stat %>%
    map(function(cvr) {
      r <- roc(response = cv_fit$y_01, predictor = cvr[ , lambda_ind],
               algorithm = 3, direction = '<',
               quiet = TRUE)

      ss <- tibble(spec = r$specificities, sens = r$sensitivities) %>%
        group_by(spec) %>%
        summarize(sens = max(sens))

      list(auc = r$auc[[1]],
           ss = ss)
    })

  auc <- map_dbl(rocs, `[[`, 'auc')
  ss <- map_dfr(rocs, `[[`, 'ss', .id = 'repl')

  auc_label <- sprintf('AUC=%.02f [%.02f, %.02f]',
                       median(auc),
                       quantile(auc, 0.25),
                       quantile(auc, 0.75))

  ss %>%
    group_by(spec) %>%
    summarize(sens_avg = median(sens),
              sens_lw = quantile(sens, 0.25),
              sens_up = quantile(sens, 0.75),
              .groups = 'drop') %>%
    bind_rows(tibble(spec = 1.00001, sens_avg = 0, sens_lw = 0, sens_up = 0)) %>%
    bind_rows(tibble(spec = -0.00001, sens_avg = 1, sens_lw = 1, sens_up = 1)) %>%
    arrange(spec, -sens_avg) %>%
    ggplot(aes(x = -spec,
               y = sens_avg,
               ymin = sens_lw,
               ymax = sens_up)) +
    geom_step(linewidth = 1.5) +
    pammtools::geom_stepribbon(alpha = 0.3, linetype = '11') +
    geom_abline(intercept = 1, slope = 1, linetype = '22') +
    geom_label(label = auc_label,
               x = -0.35,
               y = 0.1,
               size = 3) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                       labels = function (x) sprintf("%.02f", abs(x))) +
    scale_x_continuous(breaks = seq(-1, 0, by = 0.25),
                       labels = function (x) sprintf("%.02f", abs(x))) +
    labs(x = "Specificity", y = "Sensitivity") +
    plottheme() +
    coord_fixed(xlim = c(-1, 0), ylim = c(0, 1))
}

elnet_cv <- function (x, y, alpha, nfolds, cv_repl,
                      ncores = 1) {
  library(parallel)
  requireNamespace('glmnet')
  y_01 <- as.integer(y) - 1L
  keep_obs <- which(rowSums(x^2) > .Machine$double.eps)
  y_01 <- y_01[keep_obs]
  x <- x[keep_obs, , drop = FALSE]

  x <- scale(x)

  fit_global <- glmnet::glmnet(x, y_01, family = 'binomial', alpha = alpha,
                               standardize = FALSE)

  fit_global$keep_obs <- keep_obs

  cl <- makeForkCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)

  fit_global$cv_stat <- clusterApply(
    cl = cl,
    x = seq_len(cv_repl),
    xmat = x,
    y_01 = y_01,
    nfolds = nfolds,
    fit_global = fit_global,
    fun = function (seed, y_01, xmat, nfolds, fit_global) {
      requireNamespace('glmnet')
      set.seed(seed)
      foldids <- strat_folds(y_01, nfolds = nfolds)
      folds <- split(seq_along(foldids), foldids)

      lapply(folds, function (ind) {
        train_y <- y_01[-ind]
        train_x <- xmat[-ind, , drop = FALSE]

        fit_cv <- glmnet::glmnet(train_x, train_y, family = 'binomial',
                                 standardize = FALSE,
                                 lambda = (length(train_y) / length(y_01)) *
                                   fit_global$lambda)

        predict(fit_cv, newx = xmat[ind, , drop = FALSE], type = 'response')
      }) %>%
        do.call(what = rbind) %>%
        `[`(order(unlist(folds, use.names = FALSE)), )
    })

  fit_global$adaptive_threshold <- predict(fit_global, newx = x, type = 'response') %>%
    apply(2, function (x) {
      max(x[y_01 < 1])
    })

  fit_global$cv_param_grid <- tibble(lambda = fit_global$lambda)
  fit_global$y_01 <- y_01

  fit_global
}

thresh_vs_spec <- function (cv_fit, specificities,
                            thresholds = seq(0.5, 0.95, by = 0.05)) {

  map_dfr(thresholds, function (thresh) {
    x <- sens_spec(cv_fit, thresh)

    map_dfr(specificities, function (spec) {
      spec_cv_stats <- select_best(x, spec = spec)

      if (nrow(spec_cv_stats) < 1L) {
        tibble(threshold = thresh, spec = spec)
      } else {
        spec_cv_stats %>%
          mutate(threshold = thresh, spec = spec)
      }
    })
  })
}

thresh_vs_spec_all <- function (fits, ncores = 1, ...) {
  requireNamespace('parallel')
  cl <- parallel::makeForkCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  x <- clusterApply(cl = cl, fits, fun = function (x, ...) {
    thresh_vs_spec(x, ...) %>%
      mutate(day = x$for_day)
  }, ... = ...) %>%
    bind_rows()
}
