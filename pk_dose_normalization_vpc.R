# ============================================================
# PK Portfolio Example: Dose Normalization + VPC Plot
# Author: Amaka 
# Description:
# This script simulates PK concentration-time data and generates:
# 1) Dose-normalized concentration plot
# 2) Visual predictive check (VPC)-style plot
# No real patient or company data are used.
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(123)

# -----------------------------
# 1. Simulate PK data
# -----------------------------

n_subjects <- 80
time_points <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
dose_levels <- c(30, 100, 200)

pk_data <- expand.grid(
  ID = 1:n_subjects,
  TIME = time_points
) %>%
  mutate(
    DOSE = sample(dose_levels, n_subjects, replace = TRUE)[ID],
    WT = rnorm(n_subjects, mean = 70, sd = 12)[ID],
    
    # Typical PK parameters
    CL = 8 * (WT / 70)^0.75 * exp(rnorm(n(), 0, 0.25)),
    V  = 60 * (WT / 70) * exp(rnorm(n(), 0, 0.20)),
    KA = 1.1,
    
    # One-compartment oral absorption model approximation
    CONC = (DOSE * KA / V) / (KA - CL / V) *
      (exp(-(CL / V) * TIME) - exp(-KA * TIME)),
    
    # Add residual variability
    CONC = CONC * exp(rnorm(n(), 0, 0.25)),
    
    # Dose-normalized concentration
    CONC_DN = CONC / DOSE
  )

# -----------------------------
# 2. Dose-normalized plot
# -----------------------------

p_dose_norm <- ggplot(
  pk_data,
  aes(
    x = TIME,
    y = CONC_DN,
    color = factor(DOSE),
    group = ID
  )
) +
  geom_line(alpha = 0.35) +
  geom_point(alpha = 0.65, size = 1.7) +
  scale_x_continuous(breaks = time_points) +
  labs(
    title = "Dose-Normalized Concentration-Time Profiles",
    x = "Time (hours)",
    y = "Dose-Normalized Concentration",
    color = "Dose (mg)"
  ) +
  theme_bw(base_size = 13)

print(p_dose_norm)

# -----------------------------
# 3. Create VPC-style summary
# -----------------------------

vpc_summary <- pk_data %>%
  group_by(TIME, DOSE) %>%
  summarise(
    P05 = quantile(CONC, 0.05, na.rm = TRUE),
    P50 = quantile(CONC, 0.50, na.rm = TRUE),
    P95 = quantile(CONC, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 4. VPC-style plot
# -----------------------------

p_vpc <- ggplot() +
  geom_ribbon(
    data = vpc_summary,
    aes(
      x = TIME,
      ymin = P05,
      ymax = P95,
      fill = factor(DOSE)
    ),
    alpha = 0.20
  ) +
  geom_line(
    data = vpc_summary,
    aes(
      x = TIME,
      y = P50,
      color = factor(DOSE)
    ),
    linewidth = 1.1
  ) +
  geom_point(
    data = pk_data,
    aes(
      x = TIME,
      y = CONC,
      color = factor(DOSE)
    ),
    alpha = 0.35,
    size = 1.4
  ) +
  labs(
    title = "Visual Predictive Check-Style Plot Using Simulated PK Data",
    x = "Time (hours)",
    y = "Concentration",
    color = "Dose (mg)",
    fill = "Dose (mg)"
  ) +
  theme_bw(base_size = 13)

print(p_vpc)

# -----------------------------
# 5. Save figures
# -----------------------------

dir.create("figures", showWarnings = FALSE)

ggsave(
  filename = "figures/dose_normalized_plot.png",
  plot = p_dose_norm,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "figures/vpc_style_plot.png",
  plot = p_vpc,
  width = 8,
  height = 5,
  dpi = 300
)
