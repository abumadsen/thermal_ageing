pacman::p_load(tidyverse, MCMCglmm, parallel, doBy)

# Load Data file from TPC curve with data and posteriors of the TPC (optimum_fits)
load("code/MCMCglmm/TPC/output_5to8/2026-09-18quad_age_CON5to820000_100_220000.RData")
outpath <- "code/MCMCglmm/heat_vs_cold/output_5to8"

Myburn = 20000
Mythin = 50
Mynitt = 100000 + Myburn
Nsamples = (Mynitt-Myburn)/Mythin

n_optimum_draws <- 100

####################################################################
#---------------------------- LOAD
####################################################################

#Extract posterios from TPC anlayses
b <- do.call(rbind, lapply(optimum_fits, `[[`, "Sol"))

#set range
temperature_range <- dat |>
  group_by(age) |>
  summarise(lower = min(maxtemp_13), upper = max(maxtemp_13),
            .groups = "drop")

#get unique prediction profiles
profiles <- dat |>
  distinct(age, mean_age, age_w, age_b) |>
  arrange(age, mean_age) |>
  left_join(temperature_range, by = "age") |>
  mutate(profile = row_number())
dat <- left_join(dat, select(profiles, age, mean_age, profile),
                 by = c("age", "mean_age"))

#estimate linear terms at the profiles in posteriors of TPC analyses
linear <- outer(b[, "temp"], rep(1, nrow(profiles))) +
  outer(b[, "age_w:temp"], profiles$age_w) +
  outer(b[, "age_b:temp"], profiles$age_b)

#estimate quadratic terms at the profiles in posteriors of TPC analyses
quadratic <- outer(b[, "temp2"], rep(1, nrow(profiles))) +
  outer(b[, "age_w:temp2"], profiles$age_w) +
  outer(b[, "age_b:temp2"], profiles$age_b)

#estimate topt in posteriors of TPC analyses
topt_profiles <- 20 - 5 * linear / (2 * quadratic)

supported <- quadratic < 0 & is.finite(topt_profiles) &
  sweep(topt_profiles, 2, profiles$lower, ">") &
  sweep(topt_profiles, 2, profiles$upper, "<")

supported_draws <- which(rowSums(!supported) == 0)

set.seed(2024)

selected <- supported_draws[sample.int(length(supported_draws),
                                       n_optimum_draws)]

####################################################################
#---------------------------- RUNNING MCMCglmm
####################################################################

tolerance_prior <- optimum_prior
tolerance_prior$G$G1 <- list(V = diag(6) / 7, nu = 7)

tolerance_fits <- lapply(seq_len(n_optimum_draws), function(i) {
  topt <- topt_profiles[selected[i], dat$profile]
  d <- dat |>
    mutate(cold = pmax(topt - maxtemp_13, 0) / 5,
           heat = pmax(maxtemp_13 - topt, 0) / 5,
           age_cold = age_c * cold, age_heat = age_c * heat)
  set.seed(6024 + i)
  MCMCglmm(
    cbind(Neggs, NoEggsHalf) ~
      subspec + (age_w + age_b) * (cold + heat),
    random = ~ us(1 + age_c + cold + heat + age_cold + age_heat):damid +
      female_year + year + enclosure,
    family = "multinomial2", data = d, prior = tolerance_prior,
    nitt = Mynitt, burnin = Myburn, thin = Mythin, verbose = FALSE
  )
})


#Save image
save.image(paste(outpath,"/", Sys.Date(),"heat_vs_cold_topt_5to8",paste(Myburn,Mythin,Mynitt, sep = "_"),".RData", sep = ""))
