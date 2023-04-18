library(tidyverse)
source("opi.r")
source("simR.r")
source("simDisplay.r")
source("dbTocd.r")
source("full_threshold_2.r")

set.seed(1234)

chooseOpi("SimRotterdam")
opiInitialise()

makeStim <- function(db, n) { 
    s <- list(level = dbTocd(db))
    class(s) <- "opiStaticStimulus"
    return(s)
}

FT2(est = 25, tt = 33, makeStim = makeStim, 
    fpr = 0.03, fnr = 0.00, verbose = TRUE)

df_retest <- crossing(tt = 0:35, n_times = 1:100) %>% 
  mutate(FT2_res = map(tt, ~ FT2(tt = .x, makeStim = makeStim))) %>% 
  mutate(final = map_dbl(FT2_res, ~ pluck(.x, "final")))

df_retest

df_retest %>% 
  ggplot(aes(x = final)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ tt) +
  labs(x = "Final estimates (dB)")

# start with random value
df_retest2 <- crossing(tt = 0:35, n_times = 1:100) %>% 
  group_by(row_number()) %>% 
  mutate(est = round(runif(1) * 10 + 20)) %>% 
  ungroup() %>% 
  select(-`row_number()`) %>% 
  mutate(FT2_res = map2(tt, est, ~ FT2(tt = .x, 
                                       makeStim = makeStim, 
                                       est = .y))) %>% 
  mutate(final = map_dbl(FT2_res, ~ pluck(.x, "final")))

df_retest2 %>% 
  ggplot(aes(x = final)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ tt) +
  labs(x = "Final estimates (dB)")


