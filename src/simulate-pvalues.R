library(limma)
library(SRP)
library(qvalue)
library(tidyverse)

#' Simulation function
srp_simu <- function(effect_size = 2, 
                     n_tests = 1000, 
                     n_effects = 200,
                     a = 1,
                     b = 1){
  
  if(n_effects > n_tests){
    stop("Number of tests less than number of effects!")
  }

  # Random normals
  z <- rnorm(n_tests)
  
  # First 20% are effects, es -- effect size
  z[1:n_effects] <- z[1:n_effects] + effect_size
  
  # Two-sided p-values
  p <- 2 * pnorm(-abs(z))
  
  # srp
  qobject <- qvalue(p)
  pi0 <- qobject$pi0
  d <- sum(qobject$qvalues < 0.05)
  td <- (1 - 0.05) * d
  
  # adjust pi0 by using beta prior around detected effects 
  pi0shrink <- (pi0 + a) / (1 + a + b)
  
  # true nonnulls using raw pi0
  tnn <- (1 - pi0) * n_tests
  
  # true nonnulls using adjusted pi0
  tnnshrink <- (1 - pi0shrink) * n_tests
  
  srp <- td / tnn
  # srp using adjusted truenonnulls
  srpshrink <- td / tnnshrink
  
  # True power
  truepower <- sum(qobject$qvalues[1:n_effects] < 0.05) / n_effects
  
  data.frame(truepower, srp, srpshrink, pi0, pi0shrink)
}

#' Run simulations at specified effect size
wrap_es <- function(x, a, b, neff) {
  replicate(100, 
            srp_simu(effect_size = x, a = a, b = b, n_effects = neff), 
            simplify = FALSE) %>% 
    bind_rows()
}

#' Effect sizes
simu <- data_frame(es = seq(1, 4, 0.2))

neff <- 200
#' pi0 beta prior parameters
a <- 2
b <- 1

#' Run simulations
simu <- simu %>% 
  mutate(srp = map(es, wrap_es, a = a, b = b, neff = neff))
simu_unnested <- simu %>% 
  unnest()

#' Plot simulation results.
#' Play with ylim: plot displays `SRP::srp(pvalues)` results when ylim cutoff is 1 and 
#' uncensored SRP formula results when cutoff is removed.
simu_unnested %>% 
  ggplot(aes(truepower, srp, color = es)) + 
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  geom_abline(slope = 1, intercept = 0.1, linetype = 2) +
  geom_abline(slope = 1, intercept = -0.1, linetype = 2) +
  ylim(0, 1) +
  xlim(0, 1) +
  labs(title = sprintf("Number of true effects %s", neff),
       x = "True power",
       y = "Multiple testing power (SRP)",
       caption = "Blue line, loess smooth. Dashed line, theoretical fit.") +
  scale_color_continuous(name = "Effect\nsize")

#' Correlation between SRP and true power
simu_unnested %>% 
  select(truepower, srp) %>% 
  filter(!is.na(srp), is.finite(srp)) %>%
  (function(x) cor(x$truepower, x$srp))

#' Classify estimates
simu_unnested %>% 
  select(truepower, srp) %>% 
  filter(!is.na(srp), is.finite(srp)) %>% 
  mutate(estimate = case_when(
    srp < truepower + 0.1 & srp > truepower - 0.1 ~ "good",
    srp > truepower + 0.1 ~ "over",
    srp < truepower - 0.1 ~ "under"
  )) %>% 
  group_by(estimate) %>% 
  summarise(n = n() / nrow(.))

# Plot of shrunken SRP
simu_unnested %>%
  select(es, truepower, srpshrink) %>% 
  ggplot(aes(truepower, srpshrink, color = es)) + 
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  geom_abline(slope = 1, intercept = 0.1, linetype = 2) +
  geom_abline(slope = 1, intercept = -0.1, linetype = 2) +
  ylim(0, 1) +
  xlim(0, 1) +
  labs(title = sprintf("SRP calculated using beta prior adjusted pi0\na = %s, b = %s. Number of effects %s", a, b, neff),
       x = "True power",
       y = "Multiple testing power (SRP)",
       caption = "Blue line, loess smooth. Dashed line, theoretical fit.") +
  scale_color_continuous(name = "Effect\nsize")

#' Classify shrunken estimates
simu_unnested %>% 
  select(truepower, srpshrink) %>% 
  mutate(estimate = case_when(
    srpshrink < truepower + 0.1 & srpshrink > truepower - 0.1 ~ "good",
    srpshrink > truepower + 0.1 ~ "over",
    srpshrink < truepower - 0.1 ~ "under"
  )) %>% 
  group_by(estimate) %>% 
  summarise(n = n() / nrow(.))

#' Correlation between shrunken SRP and true power
simu_unnested %>% 
  select(truepower, srpshrink) %>% 
  (function(x) cor(x$truepower, x$srpshrink))
