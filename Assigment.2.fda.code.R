rm(list = ls())
library(yuima)

## 0. Data import and basic setup 

data2 <- read.csv("/Users/nikoslamprou/Desktop/MScBusinessEcon/FDA/assignments/Data2.csv",
                  header = TRUE)
BIR1    <- data2$BIR1    / 100
BIR2    <- data2$BIR2    / 100
DIR     <- data2$DIR     / 100
DIR_alt <- data2$DIR_alt / 100

# Structured note parameters
T_years <- 8          # baseline maturity in years
FV      <- 100        # face value

# Coupon rule parameters (CMS steepener style)
alpha   <- 4          # leverage on spread (BIR1 - BIR2)
beta    <- 0.5        # sensitivity to level of BIR1
F_floor <- 0.005      # floor: 0.5% minimum coupon
C_cap   <- 0.08       # cap: 8% maximum coupon

# Monte Carlo parameters
dt      <- 1/12                  # monthly time step (in years)
N_steps <- T_years * 12          # total steps for 8Y horizon
M_paths <- 5000                  # number of simulation paths

## 1. Historical descriptives & plot 

summary_hist <- data.frame(
  series = c("BIR1", "BIR2", "DIR", "DIR_alt"),
  mean   = c(mean(BIR1),    mean(BIR2),    mean(DIR),    mean(DIR_alt)),
  sd     = c(sd(BIR1),      sd(BIR2),      sd(DIR),      sd(DIR_alt)),
  min    = c(min(BIR1),     min(BIR2),     min(DIR),     min(DIR_alt)),
  max    = c(max(BIR1),     max(BIR2),     max(DIR),     max(DIR_alt))
)
summary_hist

# Figure 1: Historical interest rates
plot(BIR1, type = "l", col = "blue", ylim = c(0.01, 0.045),
     main = "Historical Interest Rates",
     ylab = "Rate (decimal)", xlab = "Time (months)")
lines(BIR2,    col = "red")
lines(DIR,     col = "darkgreen")
lines(DIR_alt, col = "orange")
legend("topright",
       legend = c("BIR1", "BIR2", "DIR", "DIR_alt"),
       col    = c("blue", "red", "darkgreen", "orange"),
       lty = 1, cex = 0.8)

## 2. Short-rate models (Vasicek & CIR) 

delta_hist <- 1/12

# Vasicek model
vas_model <- setModel(
  drift          = "kappa*(theta - r)",
  diffusion      = "sigma",
  state.variable = "r",
  solve.variable = "r"
)

# CIR model
cir_model <- setModel(
  drift          = "kappa*(theta - r)",
  diffusion      = "sigma*sqrt(r)",
  state.variable = "r",
  solve.variable = "r"
)

# Vasicek for BIR1
data_BIR1  <- setData(BIR1, delta = delta_hist)
yuima_BIR1 <- setYuima(data = data_BIR1, model = vas_model)
start_BIR1 <- list(kappa = 1, theta = mean(BIR1), sigma = sd(BIR1))
fit_BIR1   <- qmle(yuima_BIR1, start = start_BIR1)
coef_BIR1  <- coef(fit_BIR1)

# Vasicek for BIR2
data_BIR2  <- setData(BIR2, delta = delta_hist)
yuima_BIR2 <- setYuima(data = data_BIR2, model = vas_model)
start_BIR2 <- list(kappa = 1, theta = mean(BIR2), sigma = sd(BIR2))
fit_BIR2   <- qmle(yuima_BIR2, start = start_BIR2)
coef_BIR2  <- coef(fit_BIR2)

# CIR for DIR
data_DIR  <- setData(DIR, delta = delta_hist)
yuima_DIR <- setYuima(data = data_DIR, model = cir_model)
start_DIR <- list(kappa = 1, theta = mean(DIR), sigma = sd(DIR))
fit_DIR   <- qmle(yuima_DIR, start = start_DIR)
coef_DIR  <- coef(fit_DIR)

# CIR for DIR_alt
data_DIR_alt  <- setData(DIR_alt, delta = delta_hist)
yuima_DIR_alt <- setYuima(data = data_DIR_alt, model = cir_model)
start_DIR_alt <- list(kappa = 1, theta = mean(DIR_alt), sigma = sd(DIR_alt))
fit_DIR_alt   <- qmle(yuima_DIR_alt, start = start_DIR_alt)
coef_DIR_alt  <- coef(fit_DIR_alt)

# Table 2: parameters
params_table <- rbind(
  BIR1    = coef_BIR1,
  BIR2    = coef_BIR2,
  DIR     = coef_DIR,
  DIR_alt = coef_DIR_alt
)
params_table





## 3. Simulation functions 


# Vasicek: Euler-Maruyama
vasicek.sim <- function(r0, kappa, theta, sigma, dt, N, M) {
  r <- matrix(NA, nrow = N + 1, ncol = M)
  r[1, ] <- r0
  for (t in 1:N) {
    z <- rnorm(M)
    r[t + 1, ] <- r[t, ] + kappa * (theta - r[t, ]) * dt +
      sigma * sqrt(dt) * z
  }
  r
}

# CIR: Euler-Maruyama with non-negativity floor
cir.sim <- function(r0, kappa, theta, sigma, dt, N, M) {
  r <- matrix(NA, nrow = N + 1, ncol = M)
  r[1, ] <- r0
  for (t in 1:N) {
    z      <- rnorm(M)
    r_prev <- pmax(r[t, ], 0)
    r[t + 1, ] <- r_prev +
      kappa * (theta - r_prev) * dt +
      sigma * sqrt(pmax(r_prev, 0)) * sqrt(dt) * z
    r[t + 1, ] <- pmax(r[t + 1, ], 0)
  }
  r
}

# Starting values (last historical observation)
r0_BIR1    <- tail(BIR1, 1)
r0_BIR2    <- tail(BIR2, 1)
r0_DIR     <- tail(DIR, 1)
r0_DIR_alt <- tail(DIR_alt, 1)

# Unpack parameters
k1  <- coef_BIR1["kappa"];  th1  <- coef_BIR1["theta"];  s1  <- coef_BIR1["sigma"]
k2  <- coef_BIR2["kappa"];  th2  <- coef_BIR2["theta"];  s2  <- coef_BIR2["sigma"]
kD  <- coef_DIR["kappa"];    thD  <- coef_DIR["theta"];   sD  <- coef_DIR["sigma"]
kDa <- coef_DIR_alt["kappa"]; thDa <- coef_DIR_alt["theta"]; sDa <- coef_DIR_alt["sigma"]





## 4. Simulate 8Y paths (used in Task 1 & 2) 

paths_BIR1    <- vasicek.sim(r0_BIR1,    k1,  th1,  s1,  dt, N_steps, M_paths)
paths_BIR2    <- vasicek.sim(r0_BIR2,    k2,  th2,  s2,  dt, N_steps, M_paths)
paths_DIR     <- cir.sim    (r0_DIR,     kD,  thD,  sD,  dt, N_steps, M_paths)
paths_DIR_alt <- cir.sim    (r0_DIR_alt, kDa, thDa, sDa, dt, N_steps, M_paths)

# Figure 2: Fan chart for BIR1
t_grid <- (0:N_steps) * dt
plot(t_grid, paths_BIR1[, 1], type = "l",
     col = rgb(0.7, 0.7, 1, 0.3),
     ylim = c(0.01, 0.07),
     main = "Simulated BIR1 Paths (Fan Chart – first 100 paths)",
     ylab = "BIR1 (decimal)", xlab = "Years")
for (m in 2:100) {
  lines(t_grid, paths_BIR1[, m], col = rgb(0.7, 0.7, 1, 0.3))
}
lines(t_grid, rowMeans(paths_BIR1), col = "blue", lwd = 2)
abline(h = th1, col = "blue", lty = 2)
legend("topright",
       legend = c("Individual paths", "Mean path", "Long-run mean theta"),
       col = c(rgb(0.7, 0.7, 1), "blue", "blue"),
       lty = c(1, 1, 2), lwd = c(1, 2, 1), cex = 0.8)

# [optional] summary at T = 8Y
summary_sim_end <- data.frame(
  series = c("BIR1_T", "BIR2_T", "DIR_T"),
  mean   = c(mean(paths_BIR1[N_steps + 1, ]),
             mean(paths_BIR2[N_steps + 1, ]),
             mean(paths_DIR [N_steps + 1, ])),
  sd     = c(sd(paths_BIR1[N_steps + 1, ]),
             sd(paths_BIR2[N_steps + 1, ]),
             sd(paths_DIR [N_steps + 1, ]))
)
summary_sim_end




## 5. Task 1 – 8Y baseline valuation (discounting with DIR) 

# Annual coupon dates
coupon_steps <- seq(12, N_steps, by = 12)
n_coupons    <- length(coupon_steps)

# Cash-flow matrix: n_coupons x M_paths
CF <- matrix(0, nrow = n_coupons, ncol = M_paths)
for (j in 1:n_coupons) {
  step <- coupon_steps[j] + 1    # +1 because row 1 = t=0
  r1   <- paths_BIR1[step, ]
  r2   <- paths_BIR2[step, ]
  c_raw  <- alpha * (r1 - r2) + beta * r1
  c_rate <- pmin(C_cap, pmax(F_floor, c_raw))
  CF[j, ] <- c_rate * FV
}
# Add face value at final coupon
CF[n_coupons, ] <- CF[n_coupons, ] + FV

# Coupon-rate summary (excluding principal)
coupon_rates <- CF / FV
coupon_rates[n_coupons, ] <- coupon_rates[n_coupons, ] - 1

summary_coupon <- c(
  mean = mean(coupon_rates),
  sd   = sd(coupon_rates),
  min  = min(coupon_rates),
  max  = max(coupon_rates)
)
summary_coupon

# Figure 5: mean coupon per year
mean_coupon_per_year <- rowMeans(coupon_rates)
plot(1:n_coupons, mean_coupon_per_year * 100,
     type = "b", pch = 19,
     main = "Mean Simulated Coupon Rate per Year (Task 1)",
     xlab = "Year", ylab = "Mean Coupon Rate (%)")
abline(h = F_floor * 100, col = "red",  lty = 2)
abline(h = C_cap   * 100, col = "blue", lty = 2)
legend("topright",
       legend = c("Mean coupon", "Floor 0.5%", "Cap 8%"),
       col = c("black", "red", "blue"),
       lty = c(1, 2, 2), pch = c(19, NA, NA), cex = 0.8)

# Discount factors from DIR
DIR_no0 <- paths_DIR[-1, ]
cum_INT <- apply(DIR_no0 * dt, 2, cumsum)

DF <- matrix(NA, nrow = n_coupons, ncol = M_paths)
for (j in 1:n_coupons) {
  DF[j, ] <- exp(-cum_INT[coupon_steps[j], ])
}

# Path-wise PVs
PV_paths <- colSums(CF * DF)

PV_mean <- mean(PV_paths)
PV_sd   <- sd(PV_paths)
PV_se   <- PV_sd / sqrt(M_paths)
PV_CI   <- c(PV_mean - 1.96 * PV_se, PV_mean + 1.96 * PV_se)
PV_mean; PV_sd; PV_CI

# Figure 4: price distribution Task 1
hist(PV_paths, breaks = 50,
     main = "Task 1 – Simulated Prices (Baseline DIR)",
     xlab = "Present value")

# Figure 3: MC convergence
M_seq   <- seq(500, M_paths, by = 500)
PV_conv <- sapply(M_seq, function(m) mean(PV_paths[1:m]))
plot(M_seq, PV_conv, type = "l",
     main = "Monte Carlo Convergence – Task 1",
     xlab = "Number of paths (M)",
     ylab = "Running mean of PV")
abline(h = PV_mean, col = "red", lty = 2)







## 6. Task 2 – 8Y valuation with alternative discounting (DIR_alt) 

DIR_alt_no0 <- paths_DIR_alt[-1, ]
cum_INT_alt <- apply(DIR_alt_no0 * dt, 2, cumsum)

DF_alt <- matrix(NA, nrow = n_coupons, ncol = M_paths)
for (j in 1:n_coupons) {
  DF_alt[j, ] <- exp(-cum_INT_alt[coupon_steps[j], ])
}

PV_paths_alt <- colSums(CF * DF_alt)

PV_mean_alt <- mean(PV_paths_alt)
PV_sd_alt   <- sd(PV_paths_alt)
PV_se_alt   <- PV_sd_alt / sqrt(M_paths)
PV_CI_alt   <- c(PV_mean_alt - 1.96 * PV_se_alt,
                 PV_mean_alt + 1.96 * PV_se_alt)
PV_mean_alt; PV_sd_alt; PV_CI_alt

# Figure 6: price distribution Task 2
hist(PV_paths_alt, breaks = 50,
     main = "Task 2 – Simulated Prices (Alternative DIR_alt)",
     xlab = "Present value (discounted with DIR_alt)")





## 7. Task 3 – Shorter maturity 4Y with DIR 

T_years_short  <- T_years / 2           # 4 years
N_steps_short  <- T_years_short * 12    # 48 monthly steps
n_coupons_s    <- T_years_short         # 4 annual coupons
coupon_steps_s <- seq(12, N_steps_short, by = 12)

# New paths for 4Y horizon
paths_BIR1_s <- vasicek.sim(r0_BIR1, k1, th1, s1, dt, N_steps_short, M_paths)
paths_BIR2_s <- vasicek.sim(r0_BIR2, k2, th2, s2, dt, N_steps_short, M_paths)
paths_DIR_s  <- cir.sim    (r0_DIR,  kD, thD, sD, dt, N_steps_short, M_paths)

# Cash flows 4Y
CF_s <- matrix(0, nrow = n_coupons_s, ncol = M_paths)
for (j in 1:n_coupons_s) {
  step <- coupon_steps_s[j] + 1
  r1   <- paths_BIR1_s[step, ]
  r2   <- paths_BIR2_s[step, ]
  c_raw  <- alpha * (r1 - r2) + beta * r1
  c_rate <- pmin(C_cap, pmax(F_floor, c_raw))
  CF_s[j, ] <- c_rate * FV
}
CF_s[n_coupons_s, ] <- CF_s[n_coupons_s, ] + FV

# Discount factors 4Y
DIR_s_no0 <- paths_DIR_s[-1, ]
cum_INT_s <- apply(DIR_s_no0 * dt, 2, cumsum)

DF_s <- matrix(NA, nrow = n_coupons_s, ncol = M_paths)
for (j in 1:n_coupons_s) {
  DF_s[j, ] <- exp(-cum_INT_s[coupon_steps_s[j], ])
}

PV_paths_s <- colSums(CF_s * DF_s)

PV_mean_s <- mean(PV_paths_s)
PV_sd_s   <- sd(PV_paths_s)
PV_se_s   <- PV_sd_s / sqrt(M_paths)
PV_CI_s   <- c(PV_mean_s - 1.96 * PV_se_s,
               PV_mean_s + 1.96 * PV_se_s)
PV_mean_s; PV_sd_s; PV_CI_s

# Figure 7: price distribution Task 3
hist(PV_paths_s, breaks = 50,
     main = "Task 3 – Simulated Prices (Shorter 4Y Maturity)",
     xlab = "Present value (4-year note)")






## 8. Summary table across scenarios

summary_prices <- data.frame(
  scenario        = c("Task1_8Y_DIR", "Task2_8Y_DIR_alt", "Task3_4Y_DIR"),
  maturity_years  = c(T_years, T_years, T_years_short),
  discount_rate   = c("DIR", "DIR_alt", "DIR"),
  mean_PV         = c(PV_mean,     PV_mean_alt,   PV_mean_s),
  sd_PV           = c(PV_sd,       PV_sd_alt,     PV_sd_s),
  CI_low_95       = c(PV_CI[1],    PV_CI_alt[1],  PV_CI_s[1]),
  CI_high_95      = c(PV_CI[2],    PV_CI_alt[2],  PV_CI_s[2]),
  excess_over_par = c(PV_mean - FV, PV_mean_alt - FV, PV_mean_s - FV)
)
summary_prices
