## 2024-10-16
library(vegan)
library(tidyverse)

## ----------------------------------------------------------------- parameters

abundances <- c(1000, 100, 10, 1)
a_sample <- list(reads = abundances)
n_samples <- 79
n_reads <- 200
n_runs <- 10


## ------------------------------------------------------------------ functions

duplicate_sample <- function(abundance_values) {
  abundance_values %>%
    rep(x = ., n_samples)
}

amplify_sample <- function(a_sample, n_cycles) {
  if (n_cycles < 0) { stop("'n_cycles' must be >= 0") }
  if (n_cycles > 50) { stop("'n_cycles' is too large") }
  tibble(reads = a_sample) %>%
    as_tibble() %>%
    mutate(scalling_factor = runif(n = length(abundances), min = 1, max = 2),
           reads = round(reads * (scalling_factor ** n_cycles))) %>%
    select(reads) %>%
    as.data.frame()
}

sequence_sample <- function(a_sample, n_reads) {
  a_sample %>%
    as.data.frame() %>%
    vegan::rrarefy(sample = n_reads) %>%
    t() %>%
    as_tibble() %>%
    rownames_to_column(var = "species")
}

compute_stats <- function(many_samples, n_cycles) {
  many_samples %>%
    mutate(n_cycles = n_cycles) %>%
    group_by(n_cycles, species) %>%
    summarise(mean = mean(reads),
              variance = var(reads),
              .groups = "drop")
}

a_run <- function(n_cycles) {
  a_sample %>%
    duplicate_sample %>%
    purrr::modify(amplify_sample, n_cycles = n_cycles) %>%
    purrr::map_dfr(sequence_sample, n_reads = n_reads) %>%
    compute_stats(n_cycles = n_cycles)
}

many_runs <- function(n_cycles) {
  n_runs %>%
    seq() %>%
    purrr::imap_dfr(~ a_run(n_cycles = n_cycles),
                    .id = "run")
}


## ----------------------------------------------------------------------- main

many_runs(n_cycles = 0) -> metagenomics
many_runs(n_cycles = 10) -> metabarcoding

bind_rows(metagenomics, metabarcoding) %>%
  ggplot(aes(x = log(mean), y = log(variance))) +
  geom_abline(aes(intercept = 0, slope = 1))+
  geom_abline(aes(intercept = 0, slope = 2))+
  geom_point() +
  geom_point(size = 2) +
  theme_bw(base_size = 16) +
  expand_limits(x = 0, y = 0)

quit(save = "no")
