library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times", base_size = 13))
library(rstan)
library(egg)
source("../script/script_incidence.R")

load("../stanfit_sirs/stanfit_sirs_us_incidence_npi.rda")

ss <- summary(stanfit_sirs_us_incidence_npi)

N <- nrow(data_incidence)

fitdata <- data.frame(
	year=data_incidence$year,
	week=data_incidence$week,
	est=ss$summary[grepl("C\\[", rownames(ss$summary)),6][1:N],
	lwr=ss$summary[grepl("C\\[", rownames(ss$summary)),4][1:N],
	upr=ss$summary[grepl("C\\[", rownames(ss$summary)),8][1:N],
	lwr2=ss$summary[grepl("C\\[", rownames(ss$summary)),5][1:N],
	upr2=ss$summary[grepl("C\\[", rownames(ss$summary)),7][1:N]
)

npidata <- data.frame(
	year=data_incidence$year,
	week=data_incidence$week,
	est=ss$summary[grepl("npieff\\[", rownames(ss$summary)),6][1:N],
	lwr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),4][1:N],
	upr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),8][1:N],
	lwr2=ss$summary[grepl("npieff\\[", rownames(ss$summary)),5][1:N],
	upr2=ss$summary[grepl("npieff\\[", rownames(ss$summary)),7][1:N]
)

betadata <- data.frame(
	week=1:52,
	est=ss$summary[grepl("beta\\[", rownames(ss$summary)),6],
	lwr=ss$summary[grepl("beta\\[", rownames(ss$summary)),4],
	upr=ss$summary[grepl("beta\\[", rownames(ss$summary)),8],
	lwr2=ss$summary[grepl("beta\\[", rownames(ss$summary)),5],
	upr2=ss$summary[grepl("beta\\[", rownames(ss$summary)),7]
)

g1 <- ggplot(data_incidence) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_point(aes(year+week/52, proxy), shape=1) +
	geom_line(data=fitdata, aes(year+week/52, est), col="#E02938") +
	geom_ribbon(data=fitdata, aes(year+week/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(data=fitdata, aes(year+week/52, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024) +
	scale_y_continuous("Incidence proxy", limits=c(0, 8), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

g2 <- ggplot(npidata) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_hline(yintercept=1, lty=2) +
	geom_line(aes(year+week/52, est), col="#E02938") +
	geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(aes(year+week/52, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024) +
	scale_y_continuous("Relative transmission", expand=c(0, 0), limits=c(0, 1.8)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

g3 <- ggplot(betadata) +
	geom_line(aes(week, est), col="#E02938") +
	geom_ribbon(aes(week, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(aes(week, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Week", expand=c(0, 0)) +
	scale_y_continuous("Transmission rate\n(1/week)", expand=c(0, 0), limits=c(0.4, 1.5)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

gcomb <- ggarrange(g1, g2, g3, nrow=3,
									 labels=c("A", "B", "C"))

ggsave("figure1_new.pdf", gcomb, width=8, height=6)
