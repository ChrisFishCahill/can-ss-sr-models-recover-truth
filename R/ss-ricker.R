# state space model self test for carl and josh, no priors
library(rstan)
library(tidybayes)
library(tidyverse)
library(cowplot)

# state-space recruitment function used in Su and Peterson
sr_model <- function(Ut = NA) {
  wt <- rnorm(n_year - k, 0, sdp) # process noise
  vt <- rnorm(n_year, 0, sdo) # observation noise
  
  # Initialize S, R, C
  S[1:k] <- R[1:k] <- ar / b # S' = R' = ln(a)/b = equilibrium S and R
  C[1:k] <- Ut[1:k] * R[1:k]
  
  # sequentially generate new recruits, catch, and spawners
  for (t in 1:(n_year - k)) {
    R[t + k] <- a * S[t] * exp(-b * S[t] + wt[t]) # truth + process noise
    C[t + k] <- Ut[t + k] * R[t + k]
    S[t + k] <- R[t + k] - C[t + k]
  }
  
  E <- S * exp(vt) # add obs. noise to S to get E
  
  out <- tibble(
    "E" = E,
    "R" = R,
    "S" = S,
    "ln_RS" = c(
      log(R[(k + 1):n_year] / S[1:(n_year - k)]),
      rep(NA, k)
    ),
    "C" = C,
    "Ut" = Ut,
    "wt" = c(rep(NA, k), wt),
    "vt" = vt,
    "year" = 1:n_year
  )
  out
}

#-------------------------------------------------------------------------------
# Some leading parameters
#-------------------------------------------------------------------------------

k <- 2 # age at maturity
n_year <- 100
ar <- b <- 1 # ln( ricker alpha), ricker b
hmsy <- 0.5 * ar - 0.07 * (ar^2) # su and peterson relationships
smsy <- hmsy / b # su and peterson relationships
a <- exp(ar)
sdp <- 0.1 # process error sd
sdo <- 0.25 # observation error sd-- note this affects catch, not calculation of R and S
E <- S <- rep(NA, n_year) # Escapement, Stock
C <- R <- rep(NA, n_year) # Catch, Recruits

# fish to low state determined by Umax
Umax <- 0.4
Ut <- rep(NA, n_year)
relU <- seq(from = 0, to = 1, by = 0.05)
Ut[1:length(relU)] <- relU
Ut[which(is.na(Ut))] <- 1
Ut <- Ut * Umax

# visualize for jim:
set.seed(1)
sim_dat <- sr_model(Ut = Ut) 

# or not using ggplot:
# plot(sim_dat$R ~ sim_dat$S, xlab = "stock", ylab = "recruits")


#-------------------------------------------------------------------------------
# Run stan, no priors:
inits <- function() {
  list(
    "ar" = ar,
    "ln_So" = rep(log(sim_dat$S[1], k)),
    "sdo" = sdo,
    "sdp" = sdp,
    "R" = sim_dat$R[(k + 1):n_year]
  )
}

stan_data <-
  list(
    "k" = k,
    "n_year" = length(sim_dat$E),
    "E" = sim_dat$E,
    "C" = sim_dat$C[1:(length(sim_dat$C)-k)]
  )
# compile the stan model
path <- "src/ss-ricker-no-priors.stan"
m <- rstan::stan_model(path, verbose = T)

fit <-
  rstan::sampling(
    m,
    data = stan_data,
    init = inits,
    pars = c("ar", "b", "sdo", "sdp", "smsy", "hmsy", "R", "S"),
    iter = 2000, warmup = 1000, chains = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15),
    verbose = TRUE
  )

#
p1 <- fit %>%
  spread_draws(R[year]) %>%
  median_qi() %>%
  ggplot(aes(x = year, y = R, ymin = .lower, ymax = .upper)) +
  geom_pointinterval()
p1
dat <- data.frame(year = 1:n_year, Rtrue = sim_dat$R)
dat <- dat[1:(nrow(dat)-k),]
p1 <- p1 + geom_point(
  data = dat, aes(
    x = year, y = Rtrue, ymin = Rtrue,
    ymax = Rtrue
  ), shape = 16, color = "red",
  size = 2
) +
  ggqfc::theme_qfc()

p1

p2 <- fit %>%
  spread_draws(S[year]) %>%
  median_qi() %>%
  ggplot(aes(x = year, y = S, ymin = .lower, ymax = .upper)) +
  geom_pointinterval()
p2
dat <- data.frame(year = 1:n_year, Strue = sim_dat$S)
p2 <- p2 + geom_point(
  data = dat, aes(
    x = year, y = Strue, ymin = Strue,
    ymax = Strue
  ), shape = 16, color = "red",
  size = 2
) +
  ggqfc::theme_qfc()
p2

p3 <- fit %>%
  gather_draws(ar, b, sdp, sdo) %>%
  median_qi() %>%
  ggplot(aes(x = .variable, ymin = .lower, ymax = .upper, y = .value)) +
  geom_pointinterval() +
  xlab("Parameter") +
  ylab("Value") +
  ggqfc::theme_qfc()
dat <- data.frame(
  .variable = c("ar", "b", "sdo", "sdp"),
  .value = c(ar, b, sdo, sdp)
)
p3 <- p3 + geom_point(
  data = dat, aes(
    x = .variable, y = .value, ymin = .value,
    ymax = .value
  ),
  shape = 16, color = "red",
  size = 2
)
p3

p <- plot_grid(p3, p1, p2, ncol = 1)
p
ggsave("plots/self-test-su-peterson-no-priors.pdf", width = 8, height = 11)
