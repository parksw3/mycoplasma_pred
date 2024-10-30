library(dplyr)
library(tidyr)
library(lubridate)
library(rstan)
source("../script/script_incidence.R")

model <- stan_model("../stanmodel/sirs_npi_incidence.stan")

N <- nrow(data_incidence)
Npred <- 52*10
Nnonpi <- nrow(filter(data_incidence, year + week/52  < 2020 + 12/52))
npiwhich <- ceiling(((Nnonpi+1):N-Nnonpi)/4)
Nnpi <- max(npiwhich)

standata <- list(
  N=N,
  Npred=Npred,
  Nnonpi=Nnonpi,
  Nnpi=Nnpi,
  npiwhich=npiwhich,
  week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
  incidence=data_incidence$proxy,
  mu=1/80/52,
  pop=1e6,
  gamma=1/3
)

stanfit_sirs_us_incidence_npi <- sampling(model,
															 data=standata,
															 seed=104,
															 chain=4,
															 cores=4,
															 iter=4000,
															 control=list(
															 	adapt_delta=0.9
															 ))

check_hmc_diagnostics(stanfit_sirs_us_incidence_npi)
get_num_divergent(stanfit_sirs_us_incidence_npi)
get_num_max_treedepth(stanfit_sirs_us_incidence_npi)
get_low_bfmi_chains(stanfit_sirs_us_incidence_npi)
get_bfmi(stanfit_sirs_us_incidence_npi)

save("stanfit_sirs_us_incidence_npi", file="stanfit_sirs_us_incidence_npi.rda")

ss <- summary(stanfit_sirs_us_incidence_npi)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.007
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 529.6

plot(1:(N+Npred)/52, ss$summary[grepl("C\\[", rownames(ss$summary)),6], col=2)

plot(c(rep(NA, 12), zoo::rollmean(data_incidence$proxy, 12)), type="l", xlim=c(0, 1000), ylim=c(0, 50))
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6], col=2)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4], col=2)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8], col=2)

plot(data_incidence$proxy)
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("C\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("beta\\[", rownames(ss$summary)),6], ylim=c(0, 2))
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("npieff\\[", rownames(ss$summary)),6], type="l", ylim=c(0, 2))
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("S\\[", rownames(ss$summary)),6]/1e6, ylim=c(0, 1))
lines(ss$summary[grepl("S\\[", rownames(ss$summary)),4]/1e6)
lines(ss$summary[grepl("S\\[", rownames(ss$summary)),8]/1e6)
