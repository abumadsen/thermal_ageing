outpath <- "code/MCMCglmm/TPC/output_5to8"
datapath <- "data/intermediate"

pacman::p_load(tidyverse, MCMCglmm, parallel,dplyr,patchwork,ggplot2)

Myburn = 20000
Mythin = 100
Mynitt = 200000 + Myburn
Nsamples = (Mynitt-Myburn)/Mythin

n_chains <- 3

####################################################################
#---------------------------- LOAD & PREP
####################################################################

load(paste(datapath, "/","2026-09-04_breed1998_2024_clutch_zero_blackblue_thermal_quadanalysis_CON.RData", sep = ""))

# Input dataframe (tempCON) has the following columns

# maxtemp13 = mean temperature of the temperature bins (n = 13) used for binomial grouping of data. Temperature is here the average daily temperature in the 2-4 days before egg-laying.
# age = age of female (continous)
# damid = id of female (factorial)
# year = breeding year where observations were made (factorial) 
# enclosure = enclosure id of the enclosure where female and male resided (factorial) 
# subspec = Population of female (factorial: SAB, ZB, or Cross)
# Neggs = number of two-day intervals with an egg
# NoEggsHalf = number of two-day intervals without an egg

# If running restricted data analyses it is filtered to only contain females present as young and old:
#youngdams5 <- unique(tempCON$damid[tempCON$age <= 5])
#olddams8 <- unique(tempCON$damid[tempCON$age >= 8])
#tempCON_age5_8 <- tempCON[tempCON$damid %in% olddams8 & tempCON$damid %in% youngdams5,]

dat <- tempCON_age5_8 |>
  select(Neggs, NoEggsHalf, age, maxtemp_13, damid, year, enclosure, subspec)

#Counting female seasons per age
female_seasons <- distinct(dat, damid, year, age)

#Estimate mean age
female_ages <- female_seasons |>
  group_by(damid) |>
  summarise(mean_age = mean(age), .groups = "drop")

#To distinguish within-female ageing from differences among females observed at different ages, we decomposed age into within- and between-female components [@vandepol2006]
dat <- dat |>
  left_join(female_ages, by = "damid") |>
  mutate(
    across(c(damid, year, enclosure, subspec), factor),
    age_w = (age - mean_age) / 5,
    age_b = (mean_age - 7) / 5,
    age_c = (age - 7) / 5,
    temp = (maxtemp_13 - 20) / 5,
    temp2 = temp^2,
    female_year = interaction(damid, year, drop = TRUE)
  ) |>
  as.data.frame()


####################################################################
#---------------------------- RUNNING MCMCglmm
####################################################################

optimum_prior <- list(
  R = list(V = 1, nu = 1),
  G = list(
    G1 = list(V = diag(4) / 5, nu = 5),
    G2 = list(V = 1, nu = 1),
    G3 = list(V = 1, nu = 1),
    G4 = list(V = 1, nu = 1)
  )
)

optimum_fits <- lapply(seq_len(n_chains), function(i) {
  set.seed(2024 + i)
  MCMCglmm(
    cbind(Neggs, NoEggsHalf) ~
      subspec + (age_w + age_b) * (temp + temp2),
    random = ~ us(1 + age_c + temp + temp2):damid +
      female_year + year + enclosure,
    family = "multinomial2", data = dat, prior = optimum_prior,
    nitt = Mynitt, burnin = Myburn, thin = Mythin, verbose = FALSE
  )
})


#Save image
save.image(paste(outpath,"/", Sys.Date(),"quad_age_CON5to8",paste(Myburn,Mythin,Mynitt, sep = "_"),".RData", sep = ""))
