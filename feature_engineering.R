running_mean <- function (x) {
  cumsum(replace_na(x, 0)) / cumsum(!is.na(x))
}

#' Create features from cubic regression splines fitted
#' to the individual timelines and the average timeline.
#' @param x PRx data. Needs columns `minute`, `prx_trans`, `vital_status`,
#'   and `pat_id`.
#' @param xmeanspl PRx data to estimate the mean PRx curves.
#' @param resolution At what frequency to place knots for the splines (in
#'   minutes).
#' @param spline_grid grid to evaluate splines on.
spline_features <- function (x, xmeanspl = x, resolution = 360, spline_grid) {
  requireNamespace('mgcv', quietly = TRUE)
  # We'll evaluate the splines every
  if (missing(spline_grid)) {
    spline_grid <- seq(min(x$minute), max(x$minute), by = 5)
  }
  knots <- seq(min(spline_grid), max(spline_grid), by = resolution)

  fit_spline <- function (x, grouping) {
    knots_keep <- which(knots >= min(x$minute) & knots <= max(x$minute))
    knots_pat <- knots[knots_keep]
    coefs <- numeric(length(knots))

    spline_fit <- tryCatch({
      mgcv::gam(prx_trans ~ s(minute, bs = 'cs', k = length(knots_pat)),
                data = x,
                knots = list(knots_pat))
    }, error = function (e) {
      warning(paste("Cannot fit spline to patient trajectory:",
                    as.character(e)))
      return(NULL)
    })

    if (is.null(spline_fit)) {
      return(list(smoothed = tibble(vital_status = grouping$vital_status[[1]],
                                    minute = numeric(0L),
                                    prx_smooth = numeric(0L)),
                  coefs = tibble(minute = c(0, knots[-1]), coef = coefs)))
    }

    df <- tibble(vital_status = grouping$vital_status[[1]],
                 minute = spline_grid,
                 prx_smooth = NA_real_) %>%
      filter(minute >= min(x$minute),
             minute <= max(x$minute))
    df$prx_smooth <- as.vector(mgcv::predict.gam(
      spline_fit, newdata = df,
      type = 'response',
      newdata.guaranteed = TRUE))

    coefs[c(1L, knots_keep[-1L])] <- coef(spline_fit)
    coef_tbl <- tibble(minute = c(0, knots[-1]), coef = coefs)

    list(smoothed = df, coefs = coef_tbl)
  }

  area_bw_curves <- function (x, y1, y2) {
    ydiff <- y1 - y2
    good <- which(!is.na(ydiff))
    x <- x[good]
    ydiff <- abs(ydiff[good])

    n <- length(x)
    h <- diff(x)
    area <- 0
    for (i in seq(1L, n - 2L, 2L)) {
      hph <- h[[i]] + h[[i + 1L]]
      hdh <- h[[i + 1L]] / h[[i]]
      area <- area + hph / 6 *
        ((2 - hdh) * ydiff[[i]] +
           hph ^ 2 / (h[[i]] * h[[i + 1L]]) * ydiff[[i + 1L]] +
           (2 - 1 / hdh) * ydiff[[i + 2L]])
    }

    if (isTRUE(n %% 2 == 0)) {
      hph <- h[[n - 1L]] + h[[n - 2L]]
      threehth <- 3 * h[[n - 1L]] * h[[n - 2L]]
      sixh2 <- 6 * h[[n - 2L]]
      h1sq <- h[[n - 1L]] ^ 2
      area <- area +
        (2 * h1sq + threehth) / (6 * hph) * ydiff[[n]] +
        (h1sq + threehth) / sixh2 * ydiff[[n - 1L]] -
        (h1sq * h[[n - 1L]]) / (sixh2 * hph) * ydiff[[n - 2L]]
    }
    area
  }

  mean_splines <- xmeanspl %>%
    group_by(vital_status) %>%
    group_map(fit_spline)

  pat_stats <- x %>%
    filter(!is.na(prx_trans)) %>%
    group_by(vital_status, pat_id) %>%
    group_modify(function (px, grouping) {
      pat_smooth <- px %>%
        filter(!is.na(prx_trans)) %>%
        fit_spline(grouping)

      pat_smooth$smoothed %<>%
        inner_join(px, by = c('minute')) %>%
        mutate(prx_smooth = if_else(is.na(prx_trans), NA_real_,
                                    prx_smooth))

      area <- mean_splines %>%
        map_dfr(function (msp) {
          pat_smooth$smoothed %>%
            inner_join(msp$smoothed, by = 'minute',
                       suffix = c('_pat', '_avg')) %>%
            summarize(area_abs = area_bw_curves(minute, prx_smooth_avg, prx_smooth_pat),
                      area_sq = area_abs^2) %>%
            mutate(vital_status = msp$smoothed$vital_status[[1]])
        }) %>%
        pivot_wider(names_from = 'vital_status',
                    values_from = starts_with('area'))

      pat_smooth$coefs %>%
        mutate(minute = sprintf('spc_%06d', minute)) %>%
        pivot_wider(names_from = 'minute', values_from = coef) %>%
        bind_cols(area)
    }) %>%
    ungroup()
}

prep_traj_history <- function (traj) {
  traj %>%
    group_by(pat_id) %>%
    mutate(prx_pos = pmax(0, prx_trans),
           prx_trans_d = prx_trans - lag(prx_trans),
           prx_trans_d2 = prx_trans_d - lag(prx_trans_d)) %>%
    ungroup() %>%
    mutate(day = floor(minute / (24 * 60)) + 1L) %>%
    group_by(vital_status, pat_id, day) %>%
    summarize(
      prx_avg = mean(prx_trans, na.rm = TRUE),
      prx_sd = sd(prx_trans, na.rm = TRUE),
      prx_pos_avg = mean(prx_trans[prx_trans >= 0], na.rm = TRUE),
      prx_pos_sd = sd(prx_trans[prx_trans >= 0], na.rm = TRUE),
      prxd_avg = mean(prx_trans_d, na.rm = TRUE),
      prxd_sd = sd(prx_trans_d, na.rm = TRUE),
      prxd_rms = sqrt(mean(prx_trans_d^2, na.rm = TRUE)),
      prx2dsq_avg = mean(prx_trans_d2^2, na.rm = TRUE),
      prx2dsq_sd = sd(prx_trans_d^2, na.rm = TRUE),
      .groups = 'drop_last') %>%
    arrange(pat_id, day) %>%
    mutate(across(starts_with('prx'), coalesce, 0)) %>%
    group_by(vital_status, pat_id, .add = FALSE) %>%
    mutate(across(starts_with('prx'), .fns = list(cum = running_mean))) %>%
    ungroup() %>%
    mutate(across(starts_with('prx'), ~ if_else(abs(.x) < .Machine$double.eps,
                                                NA_real_, .x)))
}

#' Prepare a summary of the history of PRx trajectories.
rearrange_history <- function (x, for_day) {
  template <- expand_grid(day = seq_len(for_day),
                          pat_id = unique(x$pat_id)) %>%
    left_join(vital_status, by = 'pat_id')

  x %<>%
    select(-vital_status) %>%
    full_join(template, by = c('day', 'pat_id')) %>%
    arrange(pat_id, day)

  current <- x %>%
    filter(day == for_day) %>%
    select(-ends_with('cum'), -day) %>%
    mutate(across(starts_with('prx'), replace_na, 0))

  past <- x %>%
    filter(day < for_day) %>%
    arrange(pat_id, day) %>%
    group_by(pat_id) %>%
    fill(ends_with('cum')) %>%
    ungroup() %>%
    select(pat_id, day, ends_with('cum')) %>%
    filter(day == for_day - 1L)

  left_join(current, past, by = 'pat_id') %>%
    mutate(across(starts_with('prx'), replace_na, 0))
}

#' Prepare a summary of the history of PRx trajectories.
rearrange_history_v2 <- function (x, for_day, cumulative = NA) {
  if (!isTRUE(is.na(cumulative))) {
    x <- if (!isTRUE(cumulative)) {
      x %>%
        select(-ends_with('cum'))
    } else {
      x %>%
        select(vital_status, pat_id, day, ends_with('cum'))
    }
  }

  x %>%
    filter(day <= for_day) %>%
    pivot_wider(id_cols = c(vital_status, pat_id),
                names_from = day,
                values_from = starts_with('prx'))
}

#' Prepare statistics of the history of PRx trajectories.
#'
#' Features created:
#'   * Running averages of of the average, SD, and RMS of the
#'     0th, 1st and 2nd derivative of PRx
#'   * The slope of these running averages from the 1st and 2nd observed days.
#'   * The slope of these running averages from the 2 days prior to `for_day`
history_stats <- function (traj, for_day, vital_status) {
  change_from <- function (x, dist) {
    (x - lag(x, dist)) / dist
  }

  slope <- function (x, day, base = 1L) {
    non_na_ind <- which(!is.na(x))
    base <- min(base, length(non_na_ind))
    base_ind <- non_na_ind[[base]]
    first_val <- x[[base_ind]]
    (x - first_val) / (day - day[[base_ind]])
  }

  template <- expand_grid(day = seq_len(max(traj$day)),
                          pat_id = unique(traj$pat_id)) %>%
    left_join(vital_status, by = 'pat_id')

  x <- traj %>%
    select(pat_id, day, ends_with('cum')) %>%
    group_by(pat_id) %>%
    mutate(across(ends_with('cum'),
                  list(chg1 = ~ slope(.x, day, 1L),
                       chg2 = ~ slope(.x, day, 2L),
                       "chg-1" = ~ change_from(.x, 1),
                       "chg-2" = ~ change_from(.x, 2)))) %>%
    ungroup() %>%
    full_join(template, by = c("pat_id", "day")) %>%
    group_by(pat_id) %>%
    fill(starts_with('prx')) %>%
    ungroup()

  if (!missing(for_day)) {
    filter(x, day == for_day | (for_day > max(day) & day == max(day))) %>%
      select(-day) %>%
      mutate(across(starts_with('prx'), replace_na, 0))
  } else {
    x
  }
}

#' Prepare statistics of the history of PRx trajectories based
#' on regression spline approximations to the PRx trajectories.
#'
#' Features created:
#'   * Area between the average trajectories and the patients' trajectories
#'   * Spline coefficients
#' @param ... additional arguments for `spline_features()`
history_stats_splines <- function (x, xmean = x, for_day, ...) {
  xmean %<>%
    filter(minute <= for_day * 60 * 24)

  x %>%
    filter(minute <= for_day * 60 * 24) %>%
    spline_features(xmean = xmean, ...) %>%
    mutate(across(starts_with('area'), replace_na, 0))
}

history_model_matrix <- function (x, interactions = 1, force_keep) {
  formula <- if (isTRUE(interactions > 1)) {
    as.formula(sprintf('~ 0 + .^%d', interactions))
  } else {
    ~ 0 + .
  }

  x %<>%
    mutate(across(c(-vital_status, -pat_id, -day), replace_na, 0))

  mm <- model.matrix(formula, select(x, -vital_status, -pat_id, -day))
  keep_cols <- (colSums(mm^2) > .Machine$double.eps)

  if (!missing(force_keep)) {
    keep_cols <- keep_cols | (colnames(mm) %in% force_keep)
  }

  list(y = x$vital_status,
       pat_id = x$pat_id,
       x = mm[, keep_cols, drop = FALSE])
}
