library(tidyr)
library(lubridate)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(rstan)
library(egg)
load("../stanfit_sirs/stanfit_sirs_us_npi.rda")
load("../stanfit_sir/stanfit_sir_us_npi.rda")
load("../data/data_processed_mobility_us.rda")

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

ss_sirs <- summary(stanfit_sirs_us_npi)
ss_sir <- summary(stanfit_sir_us_npi)

ss_sirs$summary[rownames(ss_sirs$summary)=="tau",]
# 500.604628 = 9.62
# 227.327234--851.156062 (4.371678--16.36839)

N <- nrow(data_gather)
Npred <- 52*10

fitdata_sirs <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss_sirs$summary[grepl("Ifit\\[", rownames(ss_sirs$summary)),6],
	lwr=ss_sirs$summary[grepl("Ifit\\[", rownames(ss_sirs$summary)),4],
	upr=ss_sirs$summary[grepl("Ifit\\[", rownames(ss_sirs$summary)),8],
	lwr2=ss_sirs$summary[grepl("Ifit\\[", rownames(ss_sirs$summary)),5],
	upr2=ss_sirs$summary[grepl("Ifit\\[", rownames(ss_sirs$summary)),7]
)

npidata_sirs <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),6],
	lwr=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),4],
	upr=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),8],
	lwr2=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),5],
	upr2=ss_sirs$summary[grepl("npieff\\[", rownames(ss_sirs$summary)),7]
)

betadata_sirs <- data.frame(
	week=1:52,
	est=ss_sirs$summary[grepl("beta\\[", rownames(ss_sirs$summary)),6],
	lwr=ss_sirs$summary[grepl("beta\\[", rownames(ss_sirs$summary)),4],
	upr=ss_sirs$summary[grepl("beta\\[", rownames(ss_sirs$summary)),8],
	lwr2=ss_sirs$summary[grepl("beta\\[", rownames(ss_sirs$summary)),5],
	upr2=ss_sirs$summary[grepl("beta\\[", rownames(ss_sirs$summary)),7]
)

fitdata_sir <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss_sir$summary[grepl("Ifit\\[", rownames(ss_sir$summary)),6],
	lwr=ss_sir$summary[grepl("Ifit\\[", rownames(ss_sir$summary)),4],
	upr=ss_sir$summary[grepl("Ifit\\[", rownames(ss_sir$summary)),8],
	lwr2=ss_sir$summary[grepl("Ifit\\[", rownames(ss_sir$summary)),5],
	upr2=ss_sir$summary[grepl("Ifit\\[", rownames(ss_sir$summary)),7]
)

npidata_sir <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss_sir$summary[grepl("npieff\\[", rownames(ss_sir$summary)),6],
	lwr=ss_sir$summary[grepl("npieff\\[", rownames(ss_sir$summary)),4],
	upr=ss_sir$summary[grepl("npieff\\[", rownames(ss_sir$summary)),8],
	lwr2=ss_sir$summary[grepl("npieff\\[", rownames(ss_sir$summary)),5],
	upr2=ss_sir$summary[grepl("npieff\\[", rownames(ss_sir$summary)),7]
)

betadata_sir <- data.frame(
	week=1:52,
	est=ss_sir$summary[grepl("beta\\[", rownames(ss_sir$summary)),6],
	lwr=ss_sir$summary[grepl("beta\\[", rownames(ss_sir$summary)),4],
	upr=ss_sir$summary[grepl("beta\\[", rownames(ss_sir$summary)),8],
	lwr2=ss_sir$summary[grepl("beta\\[", rownames(ss_sir$summary)),5],
	upr2=ss_sir$summary[grepl("beta\\[", rownames(ss_sir$summary)),7]
)

g1 <- ggplot(data_gather) +
	geom_vline(xintercept=2015:2034, lty=3) +
	geom_point(aes(year+week/52, value*100), shape=1) +
	geom_line(data=fitdata_sir, aes(year+week/52, est*100), col="darkblue") +
	geom_ribbon(data=fitdata_sir, aes(year+week/52, ymin=lwr*100, ymax=upr*100), fill="darkblue", alpha=0.3) +
	geom_ribbon(data=fitdata_sir, aes(year+week/52, ymin=lwr2*100, ymax=upr2*100), fill="darkblue", alpha=0.5) +
	geom_line(data=fitdata_sirs, aes(year+week/52, est*100), col="#E02938") +
	geom_ribbon(data=fitdata_sirs, aes(year+week/52, ymin=lwr*100, ymax=upr*100), fill="#E02938", alpha=0.3) +
	geom_ribbon(data=fitdata_sirs, aes(year+week/52, ymin=lwr2*100, ymax=upr2*100), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0), limits=c(2015, 2024.5)) +
	scale_y_continuous("Test positivity (%)", limits=c(0, 3.3), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

g2 <- ggplot(NULL) +
	geom_vline(xintercept=2015:2034, lty=3) +
	geom_line(data=npidata_sir, aes(year+week/52, est), col="darkblue") +
	geom_ribbon(data=npidata_sir, aes(year+week/52, ymin=lwr, ymax=upr), fill="darkblue", alpha=0.3) +
	geom_ribbon(data=npidata_sir, aes(year+week/52, ymin=lwr2, ymax=upr2), fill="darkblue", alpha=0.5) +
	geom_line(data=npidata_sirs, aes(year+week/52, est), col="#E02938") +
	geom_ribbon(data=npidata_sirs, aes(year+week/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(data=npidata_sirs, aes(year+week/52, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0), limits=c(2015, 2024.5)) +
	scale_y_continuous("Relative transmission", expand=c(0, 0), limits=c(0, 1.8)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

g3 <- ggplot(betadata_sir) +
	geom_line(aes(week, est), col="darkblue") +
	geom_ribbon(aes(week, ymin=lwr, ymax=upr), fill="darkblue", alpha=0.3) +
	geom_ribbon(aes(week, ymin=lwr2, ymax=upr2), fill="darkblue", alpha=0.5) +
	geom_line(data=betadata_sirs, aes(week, est), col="#E02938") +
	geom_ribbon(data=betadata_sirs, aes(week, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(data=betadata_sirs, aes(week, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Week", expand=c(0, 0)) +
	scale_y_continuous("Transmission rate\n(1/week)", expand=c(0, 0), limits=c(0, 3.2)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

gcomb <- ggarrange(g1, g2, g3, nrow=3,
					labels=c("A", "B", "C"))

ggsave("figure_sirs_fit.pdf", gcomb, width=8, height=6)
