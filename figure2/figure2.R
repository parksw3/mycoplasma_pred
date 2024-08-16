library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times", base_size = 13))
library(rstan)
library(egg)
load("../stanfit_sirs/stanfit_sirs_us_npi.rda")
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

ss <- summary(stanfit_sirs_us_npi)

N <- nrow(data_gather)
Npred <- 52*10

fitdata <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss$summary[grepl("Ifit\\[", rownames(ss$summary)),6],
	lwr=ss$summary[grepl("Ifit\\[", rownames(ss$summary)),4],
	upr=ss$summary[grepl("Ifit\\[", rownames(ss$summary)),8],
	lwr2=ss$summary[grepl("Ifit\\[", rownames(ss$summary)),5],
	upr2=ss$summary[grepl("Ifit\\[", rownames(ss$summary)),7]
)

Sdata <- data.frame(
	year=c(2014, rep(2015:2040, each=52)[1:(N+Npred-1)]),
	week=c(52, rep(1:52, 30)[1:(N+Npred-1)]),
	est=ss$summary[grepl("S\\[", rownames(ss$summary)),6],
	lwr=ss$summary[grepl("S\\[", rownames(ss$summary)),4],
	upr=ss$summary[grepl("S\\[", rownames(ss$summary)),8],
	lwr2=ss$summary[grepl("S\\[", rownames(ss$summary)),5],
	upr2=ss$summary[grepl("S\\[", rownames(ss$summary)),7]
)

fitdata %>%
	filter(est==max(est))

g1 <- ggplot(data_gather) +
	geom_vline(xintercept=2015:2034, lty=3) +
	geom_point(aes(year+week/52, value*100), shape=1) +
	geom_line(data=fitdata, aes(year+week/52, est*100), col="#E02938") +
	geom_ribbon(data=fitdata, aes(year+week/52, ymin=lwr*100, ymax=upr*100), fill="#E02938", alpha=0.3) +
	geom_ribbon(data=fitdata, aes(year+week/52, ymin=lwr2*100, ymax=upr2*100), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0)) +
	scale_y_continuous("Test positivity (%)", limits=c(0, 12), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

fitdata$est[which.max(fitdata$est)]
fitdata$lwr[which.max(fitdata$est)]
fitdata$upr[which.max(fitdata$est)]

max(data_gather$value)

g2 <- ggplot(Sdata) +
	geom_vline(xintercept=2015:2034, lty=3) +
	geom_line(aes(year+week/52, est), col="#E02938") +
	geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(aes(year+week/52, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	scale_x_continuous("Year", expand=c(0, 0)) +
	scale_y_continuous("Proportion\nsusceptible") +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

gcomb <- ggarrange(g2, g1, nrow=2,
									 labels=c("A", "B"))

ggsave("figure2.pdf", gcomb, width=8, height=4)
