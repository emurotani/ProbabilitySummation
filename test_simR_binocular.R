library(tidyverse)
source("opi.r")
source("simR_binocular.r")
source("simDisplay.r")
source("dbTocd.r")
source("full_threshold_2.r")

set.seed(1234)

chooseOpi("SimRotterdamBinocular")
opiInitialise()

makeStim <- function(db, n) { 
    s <- list(level = dbTocd(db))
    class(s) <- "opiStaticStimulus"
    return(s)
}

FT2(est = 25, tt = c(10, 10), makeStim = makeStim, 
    fpr = 0.03, fnr = 0.00, verbose = TRUE)

crossing(tt1 = 0:3, tt2 = 0:3, n_times = 1:3) %>%
  mutate(tt = map2(tt1, tt2, ~c(.x, .y))) %>%
  mutate(FT_res = map(tt, ~ FT2(tt = .x, makeStim = makeStim))) %>%
  mutate(final = map_dbl(FT_res, ~ pluck(.x, "final")))

