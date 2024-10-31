library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(rstan)
library(egg)
source("../script/script_incidence.R")

load("../stanfit_sirs/stanfit_sirs_us_incidence_npi.rda")
load("../data/data_processed_mobility_us.rda")

ss <- summary(stanfit_sirs_us_incidence_npi)

N <- nrow(data_incidence)

npidata <- data.frame(
	year=data_incidence$year,
	week=data_incidence$week,
	est=ss$summary[grepl("npieff\\[", rownames(ss$summary)),6][1:N],
	lwr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),4][1:N],
	upr=ss$summary[grepl("npieff\\[", rownames(ss$summary)),8][1:N],
	lwr2=ss$summary[grepl("npieff\\[", rownames(ss$summary)),5][1:N],
	upr2=ss$summary[grepl("npieff\\[", rownames(ss$summary)),7][1:N]
)

data_mob <- data_processed_mobility_us %>%
	mutate(
		week=ifelse(week==53, 52, week)
	) %>%
	group_by(year, week) %>%
	summarize(
		mean=mean(mean)
	)
	
g1 <- ggplot(npidata) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_hline(yintercept=1, lty=2) +
	geom_line(aes(year+week/52, est), col="#E02938") +
	geom_ribbon(aes(year+week/52, ymin=lwr, ymax=upr), fill="#E02938", alpha=0.3) +
	geom_ribbon(aes(year+week/52, ymin=lwr2, ymax=upr2), fill="#E02938", alpha=0.5) +
	geom_line(data=data_mob, aes(year+week/52, 1+mean/100), col="#0582CA", lwd=1.4) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024,
										 limits=c(2020, 2023)) +
	scale_y_continuous("Relative changes\nin transmission", expand=c(0, 0), limits=c(0.1, 1.8)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

cor.test(filter(npidata, year+week/52 >= 2020+8/52, year+week/52 <= 2022+41/52)$est, data_mob$mean)

npi_ccf <- ccf(filter(npidata, year+week/52 >= 2020+8/52, year+week/52 <= 2022+41/52)$est, data_mob$mean)

npi_ccf$acf[npi_ccf$lag==3]<npi_ccf$acf[npi_ccf$lag==4]
npi_ccf$acf[npi_ccf$lag==4]<npi_ccf$acf[npi_ccf$lag==5]

npidata_filter <- npidata %>%
	filter(year+week/52 >= 2020+12/52, year+week/52 <= 2022+45/52)

compdata <- data.frame(
	est=npidata_filter$est,
	mean=data_mob$mean
)

cor.test(1+compdata$mean/100, compdata$est)

ggsave("figure_mobility_new.pdf", g1, width=8, height=6)
