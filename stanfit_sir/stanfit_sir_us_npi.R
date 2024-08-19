library(dplyr)
library(tidyr)
library(lubridate)
library(rstan); 

data <- read.csv("../data/Trend Mycoplasma pneumoniae Detection Rates 2024-06-30.csv")

ii <- "US"

data_gather <- data %>%
	mutate(
		date=as.Date(Week)
	) %>%
	select(-Week) %>%
	gather(key, value, -date) %>%
	mutate(
		key=factor(key, levels=c("US", "Northeast", "Midwest", "West", "South")),
		year=epiyear(date),
		week=epiweek(date),
		week=ifelse(week==53, 52, week)
	) %>%
	group_by(key, year, week) %>%
	summarize(
		value=mean(value)
	) %>%
	filter(key==ii) %>%
	head(-1)

model <- stan_model("../stanmodel/sir_npi.stan")

N <- nrow(data_gather)
Npred <- 52*10
Nnonpi <- nrow(filter(data_gather, year + week/52  < 2020 + 12/52))
npiwhich <- ceiling(((Nnonpi+1):N-Nnonpi)/4)
Nnpi <- max(npiwhich)

standata <- list(
  N=N,
  Npred=Npred,
  Nnonpi=Nnonpi,
  Nnpi=Nnpi,
  npiwhich=npiwhich,
  week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
  positivity=data_gather$value,
  mu=1/80/52,
  pop=1,
  gamma=1/3
)

stanfit_sir_us_npi <- sampling(model,
															 data=standata,
															 seed=102,
															 chain=4,
															 cores=4,
															 iter=8000,
															 control=list(adapt_delta=0.95,
															 						 max_treedepth=12))

check_hmc_diagnostics(stanfit_sir_us_npi)
get_num_divergent(stanfit_sir_us_npi)
get_num_max_treedepth(stanfit_sir_us_npi)
get_low_bfmi_chains(stanfit_sir_us_npi)
get_bfmi(stanfit_sir_us_npi)

save("stanfit_sir_us_npi", file="stanfit_sir_us_npi.rda")

ss <- summary(stanfit_sir_us_npi)

max(ss$summary[which(!is.na(ss$summary[,10])),10]) ## 1.011
min(ss$summary[which(!is.na(ss$summary[,10])),9]) ## 690

sort(ss$summary[,9])

ss$summary[grepl("sigma", rownames(ss$summary)),]

plot(1:(N+Npred)/52, ss$summary[grepl("Ifit\\[", rownames(ss$summary)),6], col=2)

plot(c(rep(NA, 12), zoo::rollmean(data_gather$value, 12)), type="l", xlim=c(0, 1000))
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),6], col=2)
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),4], col=2)
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),8], col=2)

plot(data_gather$value)
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),6])
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("Ifit\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("beta\\[", rownames(ss$summary)),6], ylim=c(0, 4))
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("beta\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("npieff\\[", rownames(ss$summary)),6], type="l", ylim=c(0, 2))
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("npieff\\[", rownames(ss$summary)),8])

plot(ss$summary[grepl("S\\[", rownames(ss$summary)),6], ylim=c(0.1, 0.3))
lines(ss$summary[grepl("S\\[", rownames(ss$summary)),4])
lines(ss$summary[grepl("S\\[", rownames(ss$summary)),8])
