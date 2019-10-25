' get_signature_set.R

Usage: get_signature_set.R -s SIGNATURES ( -n NUMBER | -m METRICS ) -o OUTPUT

Options:
    -s --signatures SIGNATURES      Path to cohort signatures
    -n --number NUMBER              Manually chosen number of mutation signatures to use as the model
    -m --metrics METRICS            Path to cohort metrics for automatic model selection
    -o --output OUTPUT              Path to output reference mutation set for cohort
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

signatures <- read_tsv(args[['signatures']])
if (! is.null(args[['metrics']])) {
    metrics <- read_tsv(args[['metrics']])

    chosen_model = metrics %>%
      spread(metric, value) %>%
      mutate(
        stability_proportion = (stability - min(stability)) / (max(stability) - min(stability)),
        reconstruction_proportion = (reconstructionError - min(reconstructionError)) / (max(reconstructionError) - min(reconstructionError)),
        combined = stability_proportion - reconstruction_proportion
      ) %>%
      filter(combined == max(combined)) %>%
      filter(row_number() == n()) %>%
      .$n_signatures
} else {
    chosen_model = as.numeric(args[['number']])
}

signatures %>%
  filter(n_signatures == chosen_model) %>%
  mutate(
    mutation_type = gsub('(.>.) (.)\\.(.)', '\\2[\\1]\\3', class),
    signature = factor(signature, levels = paste0('V', 1:length(unique(signature))))
  ) %>%
  select(
    mutation_type, signature, proportion
  ) %>%
  spread(signature, proportion) %>%
  write_tsv(args[['output']])
