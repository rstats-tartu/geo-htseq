library(qvalue)
library(tidyverse)
#' Simulation function
srp_simu <- function(effect_size = 2,
                     a = 1,
                     b = 1,
                     n_tests = 1000, 
                     n_effects = 200){
  
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
  qobject <- try(qvalue(p))
  
  if(inherits(qobject, "try-error")) {
    return(NA)
  }
  
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
wrap_es <- function(es, a, b, neff) {
  message(sprintf("Simulating with es %s, a %s, b %s, n effects %s", es, a, b, neff))
  replicate(10, 
            srp_simu(effect_size = es, a = a, b = b, n_effects = neff), 
            simplify = FALSE) %>% 
    bind_rows()
}

#' Effect sizes
simu <- data_frame(es = seq(1, 4, 0.2))

neff <- 200
n_tests <- 1000
#' pi0 beta prior parameters
pi0mu <- (n_tests - neff) / n_tests
get_ab <- function(mu, sigma){
  x <- ((mu * (1 - mu)) / sigma ^ 2) - 1
  a <- mu * x
  b <- (1 - mu) * x
  c(a, b)
}

a <- get_ab(pi0mu, 0.1)[1]
b <- get_ab(pi0mu, 0.1)[2]

#' Run simulations
simu <- simu %>% 
  mutate(srp = map(es, wrap_es, a = a, b = b, neff = neff))
simu_unnested <- simu %>% 
  unnest()

#' ## Naive simulation results
#' 
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

#' ## Use prior on pi0
#' This is what happens when we set prior around expected proportion of true nulls (pi0).
#' Plot of shrunken SRP
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
  na.omit() %>% 
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
  na.omit %>% 
  (function(x) cor(x$truepower, x$srpshrink))

#' ## Test different prior pi0
#'
pi0mu <- seq(0.1, 0.9, by = 0.1)
shape_par <- lapply(pi0mu, get_ab, sigma = 0.1) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
colnames(shape_par) <- c("a", "b")
shape_par <- cbind(shape_par, pi0mu)

#' ### Run simulations
#' Vector of effect sizes.
es <- seq(1, 4, 0.2)
#' Merge es with pi0 shape parameters
simu <- tibble(es, shape = rep(list(shape_par), length(es))) %>% 
  unnest()
neff <- seq(100, 900, by = 100)
simu <- tibble(simu = rep(list(simu), length(neff)), neff) %>% 
  unnest()

#' Run simulations.
#+ message=FALSE
simu <- simu %>% 
  mutate(srp = pmap(list(es, a, b, neff), wrap_es))
simu_unnested <- simu %>% 
  unnest()

#' Plot simulation grid
#+ fig.height=7, fig.width=7
simu_unnested %>%
  mutate_at(vars(a, b), round, digits = 2) %>% 
  ggplot(aes(truepower, srpshrink, color = es)) + 
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  geom_abline(slope = 1, intercept = 0.1, linetype = 2) +
  geom_abline(slope = 1, intercept = -0.1, linetype = 2) +
  facet_grid(neff ~ pi0mu) +
  ylim(0, 1) +
  xlim(0, 1) +
  labs(title = "SRP calculated using beta prior adjusted pi0.",
       subtitle = "Columns, proportion of true nulls (pi0).\nRows, number of true effects.",
       x = "True power",
       y = "Multiple testing power (SRP)",
       caption = "Blue line, loess smooth. Dashed line, theoretical fit.") +
  scale_color_continuous(name = "Effect\nsize")

#' ## More simulations
neff <- 200
set.seed(515172)
#' Random normals
z <- rnorm(n_tests)
#' First 20% are effects, es -- effect size
z[1:neff] <- z[1:neff] + 3

#' Two-sided p-values
p <- 2 * pnorm(-abs(z))
hist(p)
qobject <- qvalue(p)
pi0 <- qobject$pi0
d <- sum(qobject$qvalues < 0.05)
td <- (1 - 0.05) * d

#' adjust pi0 by using beta prior around detected effects 
shape_par <- lapply(seq(0, 1, length.out = 50), get_ab, sigma = 0.1) %>% 
  do.call(rbind, .)
pi0shrink <- (pi0 + shape_par[,1]) / (1 + shape_par[,1] + shape_par[,2])

#' true nonnulls using raw pi0
tnn <- (1 - pi0) * 1000

#' true nonnulls using adjusted pi0
tnnshrink <- (1 - pi0shrink) * 1000

#' Calculate srp using raw and adjusted truenonnulls
srp <- td / tnn
srpshrink <- td / tnnshrink

# True power
truepower <- sum(qobject$qvalues[1:neff] < 0.05) / neff

data_frame(pi0shrink, srpshrink) %>% 
  filter(srpshrink != 0) %>% 
  ggplot(aes(pi0shrink, srpshrink)) +
  geom_point(color = "blue") +
  geom_line(color = "blue", size = 1) +
  geom_vline(xintercept = pi0, linetype = 2) +
  geom_hline(yintercept = srp, linetype = 2) +
  geom_hline(yintercept = truepower, color = "red") +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  labs(title = "SRP calculated over grid of beta prior adjusted pi0.",
       x = "pi0",
       y = "Multiple testing power (SRP)",
       caption = "Blue, SRP calculated over shrunken pi0.\nRed line, true power.\nDashed vertical line, unadjusted pi0.\nDashed horizontal line, raw SRP calculated using unadjusted pi0.")

