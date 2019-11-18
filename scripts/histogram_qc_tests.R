library(testthat)

breaks <- 20
N <- 10000
pi0 <- 0.95

z <- runif(pi0 * N)

# uniform distribution
y <- runif((1 - pi0) * N)
x <- c(z, y)
uniform <- replicate(1000, histogram_qc(runif(N), breaks = breaks)[[2]])
test_that("is uniform", prop.test(sum(uniform == "uniform"), 1000, p = 0.95)$p.value > 0.05)

# anti-conservative, fdr adjusted, distribution
y <- replicate((1 - pi0) * N, t.test(rnorm(3, 1))$p.value)
x <- c(z, y)
anti_con <- replicate(1000, histogram_qc(c(runif(pi0 * N), replicate((1 - pi0) * N, t.test(rnorm(6, 1))$p.value)), breaks = breaks)[[2]])
test_that("is anti-conservative", prop.test(sum(anti_con == "anti-conservative"), 1000, p = 0.95)$p.value > 0.05)

# conservative, fdr adjusted, distribution
x <- p.adjust(x, method = "fdr")
con <- replicate(1000, histogram_qc(c(runif(pi0 * N), p.adjust(replicate((1 - pi0) * N, t.test(rnorm(6, 1))$p.value), method = "fdr")), breaks = breaks)[[2]])
test_that("is conservative", prop.test(sum(con == "conservative"), 1000, p = 0.95)$p.value > 0.05)
