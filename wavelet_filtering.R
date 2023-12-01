wavelet_movingavg <- function (x, target_res) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package `zoo` is not available.")
  }
  if (!requireNamespace("wavethresh", quietly = TRUE)) {
    stop("Package `wavethresh` is not available.")
  }

  all_minutes <- tibble(minute = seq(min(x$minute),
                                     2^ceiling(log2(diff(range(x$minute)))) +
                                       min(x$minute) - 1,
                                     by = 1))

  global_wv_smooth <- x %>%
    group_by(pat_id) %>%
    group_modify(function (x, ...) {
      set.seed(123)
      x %<>% full_join(all_minutes, by = 'minute') %>%
        arrange(minute) %>%
        mutate(prx_trans_na = coalesce(prx_trans,
                                       rnorm(n(), sd = 0.01)))

      # Low-pass filter
      wlfit <- wavethresh::wd(x$prx_trans_na) %>% {
        wavethresh::threshold(., by.level = TRUE,
                              levels = 10:(wavethresh::nlevelsWT(.) - 1))
      }

      nlev <- wavethresh::nlevelsWT(wlfit)
      rough_grid_space <- 2^(nlev - target_res)

      non_na_minutes <- x %>%
        mutate(minute = minute - min(minute)) %>%
        filter(!is.na(prx_trans)) %>%
        `[[`('minute')

      wv <- tibble(minute = seq_len(2^nlev) - 1L,
                   wavelet = wavethresh::wr(wlfit),
                   wavelet_ma = zoo::rollmean(
                     wavelet,
                     k = rough_grid_space - 1L,
                     fill = NA_real_)) %>%
        pivot_longer(-minute, values_to = 'prx_trans',
                     names_to = 'smoothing') %>%
        filter(smoothing == 'wavelet' |
                 minute %% rough_grid_space == 0L &
                 !is.na(prx_trans)) %>%
        mutate(prx_trans = if_else(minute %in% non_na_minutes,
                                   prx_trans,
                                   NA_real_),
               t0 = minute,
               minute = minute + min(x$minute)) %>%
        bind_rows(transmute(x,
                            minute = minute,
                            prx_trans,
                            smoothing = 'none'))
    }) %>%
    ungroup()
}

# global_wv_smooth <- wavelet_movingavg(parsed_data, target_res = 11L)
# global_wv_smooth_10 <- wavelet_movingavg(parsed_data, target_res = 10L)
#
# global_wv_smooth %>%
#   filter(smoothing == 'wavelet_ma') %>%
#   mutate(t = t0 / min(t0)) %>%
#   filter(t <= max(t[!is.na(prx_trans)])) %>%
#   saveRDS("cache/smoothed-data-wavelet_ma.rds")
#
# global_wv_smooth_10 %>%
#   filter(smoothing == 'wavelet_ma') %>%
#   mutate(t = t0 / min(t0)) %>%
#   filter(t <= max(t[!is.na(prx_trans)])) %>%
#   saveRDS("cache/smoothed-data-wavelet_ma-10.rds")
