library(readxl)
library(lubridate)
library(magrittr)
library(dplyr)
library(purrr)
library(stringr)

vital_status <- read_xlsx('data/vital_status.xlsx',
                          skip = 1,
                          col_names = c('pat_nr', 'pat_id', 'vital_status', 'avg_prx',
                                        'del',
                                        'palliative_withdrawl',
                                        'non_tbi_death',
                                        'notes')) %>%
  mutate(vital_status = factor(vital_status, levels = c('Good', 'Poor'),
                               labels = c('good', 'poor')),
         palliative_withdrawl = !is.na(palliative_withdrawl),
         non_tbi_death = !is.na(non_tbi_death)) %>%
  select(-del, -avg_prx, -notes) %>%
  filter(!is.na(pat_nr))

raw_data <- c(
  list.files(file.path('data', 'Bad Outcomes'), full.names = TRUE),
  list.files(file.path('data', 'Good Outcomes'), full.names = TRUE)) %>%
  map_dfr(function (fname) {
    suppressMessages({
      read_xlsx(fname, sheet = 1, col_names = FALSE,
                skip = 1,
                col_types = c('numeric', 'date',
                              'numeric', 'skip',
                              'numeric', 'skip',
                              'numeric', 'skip')) %>%
        mutate(pat_id = paste(str_match(fname, r'(([0-9]{2})-?([A-Z]{2})\.xlsx$)')[1, 2:3],
                              collapse = '-'))
    })
  }) %>%
  left_join(select(vital_status, pat_id, pat_nr), by = 'pat_id')

first_non_na <- function (x) {
  if (length(x) == 1L) {
    return(x)
  }

  non_nas <- which(!is.na(x))
  if (length(non_nas) > 0L) {
    return(x[[non_nas[[1]]]])
  }

  x[[1L]]
}

parsed_data <- raw_data %>%
  select(pat_id,
         pat_nr,
         day = ...1,
         time = ...2,
         prx = ...3,
         cpp = ...4,
         map = ...5) %>%
  mutate(minute = 24 * 60 * (day - 1) +
           60 * hour(time) + minute(time),
         prx_trans = atanh(prx)) %>%
  group_by(pat_id, pat_nr, minute) %>%
  summarize(across(everything(), first_non_na),
            .groups = 'drop')
