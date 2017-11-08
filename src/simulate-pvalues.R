library(limma)
library(SRP)
library(qvalue)
library(tidyverse)

#' Simulation function
srp_simu <- function(es = 2){
  
  ## Random normals
  z <- rnorm(1000)
  ## First 20% are effects, es -- effect size
  z[1:200] <- z[1:200] + es
  
  # Two-sided p-values
  p <- 2 * pnorm(-abs(z))
  
  
  
  ## True power
  truepower <- sum(qvalue(p)$qvalues[1:200] < 0.05) / 200
  ## SRP
  srp <- try(srp(p), silent = TRUE)
  
  if(inherits(srp, "try-error")){
    srp <- NA
  } else {
    srp <- srp$SRP
  }
  
  data.frame(truepower, srp)
}

#' Run simulations at specified effect size
wrap_es <- function(x) {
  replicate(100, srp_simu(x), simplify = FALSE) %>% 
    bind_rows()
}

simu <- data_frame(es = seq(1, 4, 0.2))
simu <- simu %>% 
  mutate(srp = map(es, wrap_es))
simu_unnested <- simu %>% 
  unnest() %>%
  na.omit() 

#' Plot simulation results
simu_unnested %>% 
  ggplot(aes(truepower, srp, color = es)) + 
  geom_point() +
  geom_smooth() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ylim(0, 1) + 
  xlim(0, 1) +
  labs(x = "True power",
       y = "Multiple testing power (SRP)",
       caption = "Blue line, loess smooth. Dashed line, theoretical fit.") +
  scale_color_continuous(name = "Effect\nsize")

#' Correlation between SRP and true power
cor(simu_unnested$truepower, simu_unnested$srp)

#' Classify estimates
simu_unnested %>% 
  mutate(estimate = case_when(
    srp < truepower + 0.1 & srp > truepower - 0.1 ~ "good",
    srp > truepower + 0.1 ~ "over",
    srp < truepower - 0.1 ~ "under"
  )) %>% 
  group_by(estimate) %>% 
  summarise(n = n() / nrow(.))
