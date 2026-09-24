
library(cowplot)
library(dplyr)
library(ggplot2)
library(magick)
library(tidyr)


# Load models  --------------------------------------------------------------------

# Import the saved MCMCglmm model used for thermal-performance curves.
tpc_model_workspace <- new.env()
load(
  paste0(
    "code/MCMCglmm/TPC/output_5to8/",
    "2026-09-18quad_age_CON5to820000_100_220000.RData"
  ),
  envir = tpc_model_workspace
)
tpc_age_model <- tpc_model_workspace$optimum_fits
tpc_data <- tpc_model_workspace$dat
tpc_iterations <- tpc_model_workspace$Mynitt
tpc_burnin <- tpc_model_workspace$Myburn
tpc_thinning <- tpc_model_workspace$Mythin
rm(tpc_model_workspace)

# Import the saved heat/cold MCMCglmm model and its data.
heat_cold_data <- new.env()
load(
  "code/MCMCglmm/heat_vs_cold/output_5to8/2026-09-19heat_vs_cold_topt_5to820000_50_120000.RData",
  envir = heat_cold_data
)
# The topt run saves its analysis data under the name `repro`.
heat_cold_model <- heat_cold_data$tolerance_fits
tolerance_iterations <- heat_cold_data$Mynitt
tolerance_burnin <- heat_cold_data$Myburn
tolerance_thinning <- heat_cold_data$Mythin
rm(heat_cold_data)



# Plotting settings --------------------------------------------------------------------


female_seasons <- distinct(tpc_data, damid, year, age)
female_ages <- female_seasons |>
  group_by(damid) |>
  summarise(mean_age = mean(age), .groups = "drop")

bands <- c("≤5", "6–8", ">8")
colours <- c("≤5" = "#440154", "6–8" = "#31688E", ">8" = "#FDE725")
theme_set(theme_classic(base_size = 12))

#Use the same posterior summary for all three questions.
summarise_draws <- function(draws) {
  data.frame(
    quantity = colnames(draws),
    median = apply(draws, 2, median, na.rm = TRUE),
    lower = apply(draws, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(draws, 2, quantile, probs = 0.975, na.rm = TRUE),
    row.names = NULL
  )
}



# Prepping posterios for plotting  --------------------------------------------------------------------



# Calculate the same posterior quantities used in Figures 2–4 before the
# Results text, so inline values and plotted estimates use the same draws.
b <- do.call(rbind, lapply(tpc_age_model, `[[`, "Sol"))
ages <- sort(unique(tpc_data$age))
age_w <- (ages - 7) / 5

linear <- outer(b[, "temp"], rep(1, length(ages))) +
  outer(b[, "age_w:temp"], age_w)
quadratic <- outer(b[, "temp2"], rep(1, length(ages))) +
  outer(b[, "age_w:temp2"], age_w)
topt_age <- 20 - 5 * linear / (2 * quadratic)

temperature_range <- tpc_data |>
  group_by(age) |>
  summarise(lower = min(maxtemp_13), upper = max(maxtemp_13),
            .groups = "drop")

supported <- quadratic < 0 & is.finite(topt_age) &
  sweep(topt_age, 2, temperature_range$lower, ">") &
  sweep(topt_age, 2, temperature_range$upper, "<")

topt_age[!supported] <- NA
colnames(topt_age) <- ages

age_weights <- female_seasons |>
  mutate(band = factor(case_when(
    age <= 5 ~ "≤5", age <= 8 ~ "6–8", TRUE ~ ">8"
  ), levels = bands)) |>
  count(band, age) |>
  group_by(band) |>
  mutate(weight = n / sum(n)) |>
  ungroup()

topt_band <- matrix(NA, nrow(b), 3, dimnames = list(NULL, bands))
band_age <- setNames(numeric(3), bands)
for (band in bands) {
  weights <- filter(age_weights, .data$band == .env$band)
  topt_band[, band] <- topt_age[, as.character(weights$age), drop = FALSE] %*%
    weights$weight
  band_age[band] <- sum(weights$age * weights$weight)
}

# Compare peak output on the log-odds scale, at each draw's own optimum.
# Population is held constant and female random effects are set to zero;
# the population intercept cancels from this within-draw age contrast.
peak_log_odds_age <- outer(b[, "(Intercept)"], rep(1, length(ages))) +
  outer(b[, "age_w"], age_w) - linear^2 / (4 * quadratic)
peak_log_odds_age[!supported] <- NA
colnames(peak_log_odds_age) <- ages
peak_log_odds_band <- matrix(
  NA, nrow(b), length(bands), dimnames = list(NULL, bands)
)
for (band in bands) {
  weights <- filter(age_weights, .data$band == .env$band)
  peak_log_odds_band[, band] <-
    peak_log_odds_age[, as.character(weights$age), drop = FALSE] %*%
    weights$weight
}
peak_output_change <- summarise_draws(cbind(
  `Older minus younger` = peak_log_odds_band[, ">8"] -
    peak_log_odds_band[, "≤5"]
))

optimum_changes <- cbind(
  `Older minus younger` = topt_band[, 3] - topt_band[, 1],
  `Middle minus younger` = topt_band[, 2] - topt_band[, 1],
  `Older minus middle` = topt_band[, 3] - topt_band[, 2]
)
#summarise_draws(optimum_changes)

beta <- do.call(rbind, lapply(heat_cold_model, `[[`, "Sol"))
a <- (band_age - 7) / 5
cold <- outer(beta[, "cold"], rep(1, 3)) +
  outer(beta[, "age_w:cold"], a)
heat <- outer(beta[, "heat"], rep(1, 3)) +
  outer(beta[, "age_w:heat"], a)
colnames(cold) <- colnames(heat) <- bands

cold_change <- cold[, 3] - cold[, 1]
heat_change <- heat[, 3] - heat[, 1]
tolerance_changes <- cbind(
  Cold = cold_change, Heat = heat_change,
  `Cold minus heat` = cold_change - heat_change,
  `Cold plus heat` = cold_change + heat_change
)
#knitr::kable(summarise_draws(tolerance_changes), digits = 2)

v <- do.call(rbind, lapply(heat_cold_model, `[[`, "VCV"))
plot_ages <- sort(unique(c(seq(min(ages), max(ages), by = 0.1), band_age)))
correlations <- matrix(NA, nrow(v), length(plot_ages),
                       dimnames = list(NULL, plot_ages))
for (j in seq_along(plot_ages)) {
  x <- (plot_ages[j] - 7) / 5
  cold_variance <- v[, "cold:cold.damid"] +
    2 * x * v[, "cold:age_cold.damid"] +
    x^2 * v[, "age_cold:age_cold.damid"]
  heat_variance <- v[, "heat:heat.damid"] +
    2 * x * v[, "heat:age_heat.damid"] +
    x^2 * v[, "age_heat:age_heat.damid"]
  covariance <- v[, "cold:heat.damid"] +
    x * (v[, "cold:age_heat.damid"] + v[, "age_cold:heat.damid"]) +
    x^2 * v[, "age_cold:age_heat.damid"]
  correlations[, j] <- covariance / sqrt(cold_variance * heat_variance)
}
rho_band <- correlations[, match(band_age, plot_ages)]
correlation_change <- cbind(
  `Older minus younger` = rho_band[, 3] - rho_band[, 1]
)

#knitr::kable(summarise_draws(correlation_change), digits = 2)

pairs <- list(`≤5 to >8` = c(1, 3), `≤5 to 6–8` = c(1, 2),
              `6–8 to >8` = c(2, 3))
consistency <- list()
for (trait in c("cold", "heat")) {
  variance <- v[, paste0(trait, ":", trait, ".damid")]
  covariance <- v[, paste0(trait, ":age_", trait, ".damid")]
  age_variance <- v[, paste0("age_", trait, ":age_", trait, ".damid")]
  for (comparison in names(pairs)) {
    x <- a[pairs[[comparison]][1]]
    z <- a[pairs[[comparison]][2]]
    across_ages <- variance + (x + z) * covariance + x * z * age_variance
    first_variance <- variance + 2 * x * covariance + x^2 * age_variance
    last_variance <- variance + 2 * z * covariance + z^2 * age_variance
    draws <- cbind(across_ages / sqrt(first_variance * last_variance))
    colnames(draws) <- comparison
    consistency[[paste(trait, comparison)]] <-
      mutate(summarise_draws(draws), side = trait)
  }
}

consistency <- bind_rows(consistency) |>
  mutate(comparison = factor(quantity, levels = names(pairs)))



# Figure 2  --------------------------------------------------------------------


tpc_data <- tpc_data |>
  mutate(band = factor(case_when(
    age <= 5 ~ "≤5", age <= 8 ~ "6–8", TRUE ~ ">8"
  ), levels = bands), proportion = Neggs / (Neggs + NoEggsHalf))

bin_key <- tpc_data |>
  distinct(band, maxtemp_13) |>
  arrange(band, maxtemp_13) |>
  mutate(bin_id = row_number())

tpc_data <- left_join(tpc_data, bin_key, by = c("band", "maxtemp_13"))
raw_bins <- tpc_data |>
  group_by(bin_id, band, maxtemp_13) |>
  summarise(proportion = mean(proportion), female_seasons = n(),
            .groups = "drop")
age_counts <- tpc_data |>
  distinct(damid, year, age, band) |>
  count(age, band, name = "female_seasons")

age_histogram_panel <- ggplot(age_counts, aes(age, female_seasons, fill = band)) +
  geom_col() +
  geom_text(aes(label = female_seasons), vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_x_continuous(breaks = age_counts$age) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Female age (years)", y = "Female-seasons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))


#----------- TPC -----------#


# This calculation is repeated for every fitted chain and optimum draw.
predict_bins <- function(model, data, fixed, random) {
  keep <- round(seq(1, nrow(model$Sol), length.out = 200))
  coefficients <- model$Sol[keep, ]
  covariance <- model$VCV[keep, ]
  fixed_design <- model.matrix(fixed, data)
  random_design <- model.matrix(random, data)
  eta <- coefficients %*% t(fixed_design)
  variance <- matrix(0, nrow(eta), ncol(eta))
  for (i in seq_len(ncol(random_design))) {
    for (j in seq_len(ncol(random_design))) {
      name <- paste0(colnames(random_design)[i], ":",
                     colnames(random_design)[j], ".damid")
      variance <- variance + outer(covariance[, name],
                                   random_design[, i] * random_design[, j])
    }
  }
  variance <- variance + rowSums(covariance[
    , c("female_year", "year", "enclosure", "units")])
  normal <- statmod::gauss.quad.prob(80, dist = "normal")
  probability <- matrix(0, nrow(eta), ncol(eta))
  for (k in seq_along(normal$nodes)) {
    probability <- probability + normal$weights[k] *
      plogis(eta + normal$nodes[k] * sqrt(variance))
  }
  t(rowsum(t(probability), data$bin_id) / tabulate(data$bin_id))
}
optimum_predictions <- lapply(tpc_age_model, function(model) {
  predict_bins(model, tpc_data,
               fixed = ~ subspec + (age_w + age_b) * (temp + temp2),
               random = ~ 1 + age_c + temp + temp2)
})

# tolerance_predictions <- lapply(seq_len(n_optimum_draws), function(i) {
#   topt <- topt_profiles[selected[i], tpc_data$profile]
#   d <- tpc_data |>
#     mutate(cold = pmax(topt - maxtemp_13, 0) / 5,
#            heat = pmax(maxtemp_13 - topt, 0) / 5,
#            age_cold = age_c * cold, age_heat = age_c * heat)
#   predict_bins(tolerance_fits[[i]], d,
#     fixed = ~ subspec + (age_w + age_b) * (cold + heat),
#     random = ~ 1 + age_c + cold + heat + age_cold + age_heat)
# })

optimum_curve <- do.call(rbind, optimum_predictions)
#tolerance_curve <- do.call(rbind, tolerance_predictions)
colnames(optimum_curve) <-  bin_key$bin_id
predictions <- bind_rows(
  mutate(summarise_draws(optimum_curve), model = "Quadratic optimum model")
) |>
  mutate(bin_id = as.integer(quantity),
         model = factor(model, c("Quadratic optimum model"))) |>
  left_join(bin_key, by = "bin_id")

tpc_curve_panel <- ggplot(predictions, aes(maxtemp_13, median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = band), alpha = 0.25) +
  geom_line(linewidth = 1.2) +
  geom_line(aes(colour = band), linewidth = 0.8) +
  geom_point(data = raw_bins,
             aes(y = proportion, size = female_seasons, fill = band), shape = 21) +
  facet_grid(. ~ band) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_colour_manual(values = colours, guide = "none") +
  scale_size_area(max_size = 5, name = "Female-seasons") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Mean maximum temperature (°C)", y = "Egg-laying proportion") +
  theme(legend.position = "top") +
  guides(size = guide_legend(nrow = 1))



##------------- Topt


# Get Topt Draws used for analyses
b <- do.call(rbind, lapply(tpc_age_model, `[[`, "Sol"))
ages <- sort(unique(tpc_data$age))
age_w <- (ages - 7) / 5

linear <- outer(b[, "temp"], rep(1, length(ages))) +
  outer(b[, "age_w:temp"], age_w)
quadratic <- outer(b[, "temp2"], rep(1, length(ages))) +
  outer(b[, "age_w:temp2"], age_w)
topt_age <- 20 - 5 * linear / (2 * quadratic)

temperature_range <- tpc_data |>
  group_by(age) |>
  summarise(lower = min(maxtemp_13), upper = max(maxtemp_13),
            .groups = "drop")

supported <- quadratic < 0 & is.finite(topt_age) &
  sweep(topt_age, 2, temperature_range$lower, ">") &
  sweep(topt_age, 2, temperature_range$upper, "<")

topt_age[!supported] <- NA
colnames(topt_age) <- ages

#Summarise Topt within the three age groupings used for plotting (bands)
age_weights <- female_seasons |>
  mutate(band = factor(case_when(
    age <= 5 ~ "≤5", age <= 8 ~ "6–8", TRUE ~ ">8"
  ), levels = bands)) |>
  count(band, age) |>
  group_by(band) |>
  mutate(weight = n / sum(n)) |>
  ungroup()

topt_band <- matrix(NA, nrow(b), 3, dimnames = list(NULL, bands))
band_age <- setNames(numeric(3), bands)
for (band in bands) {
  weights <- filter(age_weights, .data$band == .env$band)
  topt_band[, band] <- topt_age[, as.character(weights$age), drop = FALSE] %*%
    weights$weight
  band_age[band] <- sum(weights$age * weights$weight)
}


optimum_age <- summarise_draws(topt_age) |>
  mutate(age = as.numeric(quantity)) |>
  left_join(age_counts, by = "age")
optimum_band <- summarise_draws(topt_band) |>
  mutate(band = factor(quantity, levels = bands))

topt_age_plot <- ggplot(optimum_age, aes(age, median)) +
  geom_line() + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(aes(fill = band, size = female_seasons), shape = 21) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_size_area(max_size = 6, name = "Female-seasons") +
  scale_x_continuous(breaks = ages) +
  labs(x = "Female age (years)", y = "Thermal optimum (°C)") +
  theme(legend.position = "top") +
  guides(size = guide_legend(nrow = 1))

topt_band_plot <- ggplot(optimum_band, aes(band, median, fill = band)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = colours, guide = "none") +
  labs(x = "Age band (years)", y = "Thermal optimum (°C)")


#----------- Grid -----------#

plot_grid(
  age_histogram_panel,
  tpc_curve_panel,
  topt_band_plot,
  topt_age_plot,
  labels = c("A", "B", "C", "D"),
  label_fontface = "bold",
  label_size = 14,
  nrow = 2,
  ncol = 2,
  rel_widths = c(0.85, 1.15)
)

# Figure 3  --------------------------------------------------------------------


slopes <- bind_rows(
  mutate(summarise_draws(cold), side = "Cold"),
  mutate(summarise_draws(heat), side = "Heat")
) |>
  mutate(band = factor(quantity, levels = bands))
changes <- summarise_draws(tolerance_changes) |>
  mutate(quantity = factor(quantity, levels = colnames(tolerance_changes)))

slope_plot <- ggplot(slopes, aes(band, median, fill = band)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(shape = 21, size = 4) + facet_wrap(~ side) +
  scale_fill_manual(values = colours, guide = "none") +
  labs(x = "Age band (years)", y = "Response (log odds per 5°C)")
change_plot <- ggplot(changes, aes(quantity, median)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_x_discrete(labels = c(
    "Cold" = "Cold", "Heat" = "Heat",
    "Cold minus heat" = "Cold minus\nheat",
    "Cold plus heat" = "Cold plus\nheat"
  )) +
  labs(x = NULL, y = "Older minus younger\n(log odds per 5°C)") +
  theme(axis.text.x = element_text(size = 9))

consistency_plot <- ggplot(consistency, aes(comparison, median)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ side, labeller = as_labeller(c(cold = "Cold", heat = "Heat"))) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Age-band comparison (years)", y = "Consistency (correlation)")

plot_grid(
  plot_grid(
    slope_plot, change_plot,
    labels = c("A", "B"), label_fontface = "bold", label_size = 14,
    nrow = 1, rel_widths = c(1, 1)
  ),
  consistency_plot,
  labels = c("", "C"),
  label_fontface = "bold", label_size = 14,
  ncol = 1, rel_heights = c(1, 1)
)

# Figure 4  --------------------------------------------------------------------



correlation_curve <- summarise_draws(correlations) |>
  mutate(age = as.numeric(quantity))

band_points <- correlation_curve[match(band_age, plot_ages), ] |>
  mutate(band = factor(bands, levels = bands))

correlation_plot <- ggplot(correlation_curve, aes(age, median)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
  geom_line() +
  geom_errorbar(data = band_points, aes(ymin = lower, ymax = upper),
                width = 0.15) +
  geom_point(data = band_points, aes(fill = band), shape = 21, size = 4) +
  scale_fill_manual(values = colours, name = "Age band") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Female age (years)", y = "Cold–heat correlation")

print(correlation_plot)
