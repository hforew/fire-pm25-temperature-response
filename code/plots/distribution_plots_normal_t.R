# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PLOTS OF COUNTRY ALPHA COEFFICIENT DISTRIBUTIONS- ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Packages ##########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(tidyverse)
library(ggplot2)
library(maps)
library(countrycode)
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Import ############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# alpha_country_meta_pierce.csv: one row per country with OLS point estimate
# and SE for alpha_c2 (slope of per-capita fire PM2.5 on GMT, µg/m³/yr per
# °C). Columns: ISO3, estimate_alpha_c2, std.error_alpha_c2

estimates_se <- read_csv(here("output", "alpha_country_meta_pierce.csv"))
head(estimates_se)

# "global" is not in alpha_country_meta_pierce.csv; omit it here
countries_to_plot <- c(
  "USA", "DEU", "RUS", "AUS", "CHN",
  "IND", "ARG", "BRA", "ETH", "NGA"
)

# df = 5 Pierce time points − 2 OLS parameters (alpha_c1, alpha_c2)
df_reg <- 3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Build density data ################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Each country's regression produces a point estimate (alpha_c2) and a
# standard error (SE). Together these define a sampling distribution — the
# range of alpha_c2 values consistent with the data given sampling
# uncertainty. We compare two distributional assumptions:
#
#   Normal:  alpha_c2 follows N(mean = estimate, variance = se^2).
#            Standard large-sample approximation. Assumes the OLS estimator
#            is approximately normally distributed, which holds when n is
#            large. With only 5 observations per country (df = 3), this
#            approximation underestimates the probability
#            of extreme values.
#
#   t (df=3): (alpha_c2_hat − alpha_c2_true) / se ~ t(df=3)
#            Exact finite-sample distribution for OLS under normally
#            distributed errors. With df = 3, the t-distribution has heavier
#            tails than the Normal, correctly reflecting the extra uncertainty
#            from estimating sigma with only 3 residual degrees of freedom.
#            As df goes to infinity, t goes to Normal.
#
# Both distributions are centered on the same point estimate
# (estimate_alpha_c2) and share the same scale (std.error_alpha_c2).
# The difference is in tail behavior.

plot_data <- estimates_se %>%
  filter(ISO3 %in% countries_to_plot) %>%
  rename(country = ISO3, estimate = estimate_alpha_c2, se = std.error_alpha_c2)

# For each country, evaluate both densities on a shared x-grid spanning
# plus minus 5 SE. This captures >99.9% of Normal mass and the bulk of
# t(df=3) mass, making the heavier tails visible without excessive whitespace.
density_df <- plot_data %>%
  pmap_dfr(function(country, estimate, se) {
    x_seq <- seq(estimate - 5 * se,  # lower: 5 SEs below point estimate
                 estimate + 5 * se,  # upper: 5 SEs above point estimate
                 length.out = 500)   # 500 points gives a smooth curve
    tibble(
      country  = country,
      estimate = estimate,  # OLS point estimate — center of both distributions
      se       = se,        # OLS standard error — spread of both distributions
      x        = x_seq,
      # Normal density: mean = estimate (center), sd = se (spread).
      # dnorm evaluates f(x) = (1/sqrt(2π*se^2))*exp(−(x−estimate)^2/(2*se^2))
      normal   = dnorm(x_seq, mean = estimate, sd = se),
      # t density: shift and scale the standard t so it is centered on
      # estimate with spread se. The standard t(df) has mean 0 and scale 1:
      #   (1) (x_seq - estimate)/se: standardizes each x_seq value (the
      #       alpha_c2 x-axis grid) onto the standard t scale, the same
      #       way a z-score centers and scales to the standard Normal
      #   (2) dt(..., df): evaluates the standard t density at each
      #       standardized value — still on the "number of SEs" scale
      #   (3) / se: Jacobian correction — stretching the x-axis by se
      #       compresses the density by 1/se so the curve integrates to 1
      # Result: t with location = estimate, scale = se, df = df_reg.
      t_dist   = dt((x_seq - estimate) / se, df = df_reg) / se
    )
  }) %>%
  # Reshape to long format so both distributions share a single "density"
  # column, enabling color-mapped overlays without manual layer duplication.
  pivot_longer(cols     = c(normal, t_dist),
               names_to = "distribution",
               values_to = "density") %>%
  mutate(distribution = factor(
    distribution,
    levels = c("normal", "t_dist"),
    labels = c("Normal", paste0("t  (df = ", df_reg, ")"))
  ))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ Distribution Plots - per country ##################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# For each country, save one combined plot overlaying both distributions:
#   dist_<iso>_pierce.png  — Normal (steelblue) and t (tomato) on one axes
#
# Both curves share the same center (estimate) and scale (se), so they
# overlay directly. The dashed vertical line marks the OLS point estimate.
# The t curve will have a lower peak and heavier tails than the Normal —
# the visual gap between them shows how much small-sample uncertainty
# (df = 3) widens the credible range beyond the Normal approximation.

dist_colors <- c("Normal"              = "steelblue",
                 "t  (df = 3)" = "tomato")

for (iso in countries_to_plot) {

  df_country <- density_df %>% filter(country == iso)
  est <- unique(df_country$estimate)  # center of both distributions
  se  <- unique(df_country$se)        # spread of both distributions
  iso_lower <- tolower(iso)

  # Both distributions plotted with color mapped to distribution type.
  # Normal: symmetric, falls off rapidly beyond +/-2-3 SE; peak at
  #   1/(se*sqrt(2*pi)).
  # t(df=3): same center and scale as Normal but heavier tails — more
  #   mass beyond +/-2 SE; lower peak because total mass (= 1) is spread
  #   over a wider range.
  p <- ggplot(df_country, aes(x = x, y = density, color = distribution)) +
    # each geom_line traces one pdf across the shared x-grid
    geom_line(linewidth = 0.8) +
    # vertical dashed line at point estimate: mode and mean of both dists.
    geom_vline(xintercept = est, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = dist_colors, name = "Distribution") +
    labs(
      title = bquote(
        .(iso) ~ ": Sampling Distribution of" ~
          alpha[c2] ~ "(Pierce)"
      ),
      # +/-1 SE contains ~68% of Normal mass; t(df=3) requires wider
      # critical values for the same coverage due to heavier tails
      subtitle = bquote(
        hat(alpha)[c2] == .(round(est, 3)) %+-%
          .(round(se, 3)) ~ "µg/m³/yr per °C"
      ),
      x = expression(alpha[c2] ~ "(µg/m³/yr per °C GMT)"),
      y = "Density"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom")

  print(p)
  ggsave(here("images/regression_alpha/alpha_pierce",
              paste0("dist_", iso_lower, "_pierce.png")),
         p, width = 7, height = 5, dpi = 300)
}

print("Distribution plots saved to images/regression_alpha/alpha_pierce/")

### THE END
