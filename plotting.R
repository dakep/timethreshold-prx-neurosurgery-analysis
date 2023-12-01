library(ggplot2)

COLORPAL <- c(orange = '#D55E00',
              blue = '#0072B2',
              lightblue = '#56B4E9',
              green = '#009E73',
              cyan = '#CC79A7')

plottheme <- function(base_size = 12, base_family = '') {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(plot.title = element_text(size = rel(1)),
          plot.margin = margin(0.2, 0.4, 0.5, 0.2, 'lines'),
          panel.background = element_rect(fill = 'transparent', color = NA),
          plot.background = element_rect(fill = 'transparent', color = NA),
          legend.title = element_text(size = rel(0.7)),
          legend.text = element_text(size = rel(0.7)),
          axis.title = element_text(size = rel(0.75)),
          axis.text = element_text(size = rel(0.7)),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = 0),
          panel.grid.major = element_line(
            color = 'gray30', linewidth = rel(0.5), linetype='dotted'),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(
            fill = '#ffffff', color = 'gray50', linewidth = 0.3),
          strip.text = element_text(size = rel(0.75)),
          panel.border = element_rect(color = 'gray50', linewidth = 0.3),
          legend.background = element_blank())
}

vital_scale_guide <- guide_legend(title = 'Outcome',
                                  keywidth = grid::unit(2, 'lines'))
vital_scale_col <- scale_color_manual(
  values = c(poor = COLORPAL[['orange']], good = COLORPAL[['blue']],
             combined = COLORPAL[['green']]),
  guide = vital_scale_guide)
vital_scale_fill <- scale_fill_manual(
  values = c(poor = COLORPAL[['orange']], good = COLORPAL[['blue']],
             combined = COLORPAL[['green']]),
  guide = vital_scale_guide)
vital_scale_lty <- scale_linetype_manual(
  values = c(poor = 'solid', good = '52',
             combined = '11'),
  guide = vital_scale_guide)

plot_sensspec <- function (cv_fit) {
  cv_fit$sens_spec %>%
    select(lambda, starts_with('sens'), starts_with('spec')) %>%
    pivot_longer(-lambda, names_to = c('measure', '.value'), names_sep = '_') %>%
    mutate(measure = recode(measure, sens = "Sensitivity", spec = "Specificity")) %>%
    ggplot(aes(x = lambda, y = avg,
               ymin = pmax(0, avg - sd),
               ymax = pmin(1, avg + sd),
               color = measure,
               linetype = measure,
               shape = measure)) +
    geom_ribbon(color = NA, fill = 'black', alpha = 0.3, show.legend = FALSE) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = cv_fit$lambda_spec, linetype = '22') +
    scale_color_manual(values = c(Sensitivity = '#D55E00', Specificity = '#0072B2'),
                       guide = guide_legend(title = NULL)) +
    scale_linetype_manual(values = c(Sensitivity = 'solid', Specificity = '22'),
                          guide = guide_legend(title = NULL)) +
    scale_shape_manual(values = c(Sensitivity = 16, Specificity = 15),
                       guide = guide_legend(title = NULL)) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::label_percent(),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(x = "Model complexity (lower is more complex)", y = "Prediction accuracy") +
    plottheme()
}

plot_multicut <- function (x) {
  pct_fmt <- scales::label_percent()

  x_best <- x %>%
    group_by(day, spec) %>%
    filter(sens_avg >= max(sens_avg, na.rm = TRUE)) %>%
    filter(threshold <= min(threshold)) %>%
    ungroup()

  gd <- guide_legend(title = "Specificity",
                     keywidth = unit(2.5, 'line'))

  pl <- x_best %>%
    mutate(thresh_str = sprintf("p ≥ %.02f", threshold),
           spec_str = sprintf('≥ %s', pct_fmt(spec))) %>%
    arrange(spec) %>%
    mutate(spec_str = fct_inorder(spec_str)) %>%
    ggplot(aes(x = day, y = sens_avg,
               ymin = pmax(0, sens_avg - sens_sd),
               ymax = pmin(1, sens_avg + sens_sd),
               linetype = spec_str,
               shape = spec_str,
               color = spec_str)) +
    geom_linerange(alpha = 1, linetype = 'solid') +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = unname(COLORPAL), guide = gd) +
    scale_linetype_manual(values = c('22', '42', 'solid', '12', '11'), guide = gd) +
    scale_shape_discrete(guide = gd) +
    scale_y_continuous(labels = scales::label_percent(),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(min(x_best$day), max(x_best$day), by = 2)) +
    labs(x = "Time (days)", y = "Sensitivity") +
    plottheme()

  tbl <- x_best %>%
    mutate(sens_str = sprintf("%s ±%s", pct_fmt(sens_avg),
                              pct_fmt(sens_sd))) %>%
    select(day, spec, threshold, sens_str)

  return(list(plot = pl, table = tbl))
}
