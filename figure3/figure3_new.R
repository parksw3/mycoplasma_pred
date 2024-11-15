library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times", base_size=16))
library(rstan)
library(egg)
source("../R/simulate_sirs.R")
source("../script/script_incidence.R")
load("../stanfit_sirs/stanfit_sirs_us_incidence_npi.rda")

ss <- summary(stanfit_sirs_us_incidence_npi)

npieff <- ss$summary[grepl("npieff\\[", rownames(ss$summary)),6]

N <- nrow(data_incidence)
Npred <- length(npieff)-N
week <- c(52, rep(1:52, 300)[1:(N+Npred-1)]) 

scalevec <- seq(0, 1, length.out=5)

scalelist <- vector('list', length(scalevec))

for (i in 1:length(scalevec)) {
	scale <- scalevec[i]
	
	npieff_scale <- 1 + (npieff-1) * scale
	
	npi1 <- data.frame(
		year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
		week=week,
		npieff=npieff_scale
	)
	
	out <- simulate_sirs(beta=ss$summary[grepl("beta\\[", rownames(ss$summary)),6],
											week=week,
											mu=1/52/80,
											delta=ss$summary[grepl("delta", rownames(ss$summary)),6],
											gamma=1/3,
											rho=ss$summary[grepl("rho", rownames(ss$summary)),6],
											npieff=npi1$npieff,
											pop=1e6,
											S0=ss$summary[grepl("S0", rownames(ss$summary)),6],
											I0=ss$summary[grepl("I0", rownames(ss$summary)),6],
											tmax=nrow(npi1)) %>%
		mutate(
			year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
			week=week,
			scale=scale,
			npieff=npieff_scale
		)
	
	scalelist[[i]] <- out
}

scaledata <- scalelist %>%
	bind_rows

durvec <- c(2020, 2021, 2022, 2023)

durlist <- vector('list', length(durvec))

for (i in 1:length(durvec)) {
	dur <- durvec[i]
	
	npi1 <- data.frame(
		year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
		week=week,
		npieff=npieff
	)
	
	npi1$npieff[npi1$year > dur] <- 1
	
	out <- simulate_sirs(beta=ss$summary[grepl("beta\\[", rownames(ss$summary)),6],
											week=week,
											mu=1/52/80,
											gamma=1/3,
											delta=ss$summary[grepl("delta", rownames(ss$summary)),6],
											rho=ss$summary[grepl("rho", rownames(ss$summary)),6],
											npieff=npi1$npieff,
											pop=1e6,
											S0=ss$summary[grepl("S0", rownames(ss$summary)),6],
											I0=ss$summary[grepl("I0", rownames(ss$summary)),6],
											tmax=nrow(npi1)) %>%
		mutate(
			year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
			week=week,
			dur=dur,
			npieff=npi1$npieff
		)
	
	durlist[[i]] <- out
}

durdata <- durlist %>%
	bind_rows

immunevec <- c(5, 10, 15, 20, 25)

immunelist <- vector('list', length(immunevec))

for (i in 1:length(immunevec)) {
	immune <- immunevec[i]
	
	npi1 <- data.frame(
		year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
		week=week,
		npieff=npieff
	)
	
	out <- simulate_sirs(beta=ss$summary[grepl("beta\\[", rownames(ss$summary)),6],
											 week=week,
											 mu=1/52/80,
											 gamma=1/3,
											 delta=1/immune/52,
											 rho=ss$summary[grepl("rho", rownames(ss$summary)),6],
											 npieff=npi1$npieff,
											 pop=1e6,
											 S0=ss$summary[grepl("S0", rownames(ss$summary)),6],
											 I0=ss$summary[grepl("I0", rownames(ss$summary)),6],
											 tmax=nrow(npi1)) %>%
		mutate(
			year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
			week=week,
			immune=immune,
			npieff=npi1$npieff
		)
	
	immunelist[[i]] <- out
}

immunedata <- immunelist %>%
	bind_rows

g1 <- ggplot(scaledata) +
	geom_vline(xintercept=2015:2025, lty=3) +
	geom_line(aes(year+week/52, npieff, col=as.factor(scale)), lwd=0.8) +
	scale_x_continuous("Year", expand=c(0, 0),
										 limits=c(2019.5, 2025.5),
										 breaks=2015:2024) +
	scale_y_continuous("Relative transmission", expand=c(0, 0), limits=c(0.3, 1.4)) +
	scale_color_manual("Relative strength of NPIs", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	guides(color=guide_legend(nrow=2, title.position="top")) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "none"
	)

g2 <- ggplot(scaledata) +
	geom_vline(xintercept=2015:2025, lty=3) +
	geom_line(aes(year+week/52, C, col=as.factor(scale)), lwd=0.8) +
	annotate("segment", x=2024, xend=2024, y=20, yend=18,
					 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
	annotate("text", x=2024, y=24, label="Beginning of\nthe observed\noutbreak",
					 family="Times") +
	scale_x_continuous("Year", expand=c(0, 0),
										 limits=c(2019.5, 2025.5),
										 breaks=2015:2024) +
	scale_y_continuous("Incidence proxy", expand=c(0, 0), limits=c(0, 30)) +
	scale_color_manual("Relative strength of NPIs", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	guides(color=guide_legend(nrow=2, title.position="top")) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "bottom"
	)

g3 <- ggplot(durdata) +
	geom_vline(xintercept=2015:2025, lty=3) +
	geom_line(aes(year+week/52, npieff, col=as.factor(dur)), lwd=0.8) +
	scale_x_continuous("Year", expand=c(0, 0),
										 limits=c(2019.5, 2025.5),
										 breaks=2015:2024) +
	scale_y_continuous("Relative transmission", expand=c(0, 0), limits=c(0.3, 1.4)) +
	scale_color_manual("Duration of NPIs", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	guides(color=guide_legend(nrow=2, title.position="top")) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "none"
	)

g4 <- ggplot(durdata) +
	geom_vline(xintercept=2015:2025, lty=3) +
	geom_line(aes(year+week/52, C, col=as.factor(dur)), lwd=0.8) +
	annotate("segment", x=2024, xend=2024, y=20, yend=18,
					 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
	annotate("text", x=2024, y=24, label="Beginning of\nthe observed\noutbreak",
					 family="Times") +
	scale_x_continuous("Year", expand=c(0, 0),
										 limits=c(2019.5, 2025.5),
										 breaks=2015:2024) +
	scale_y_continuous("Incidence proxy", expand=c(0, 0), limits=c(0, 30)) +
	scale_color_manual("Duration of NPIs", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	guides(color=guide_legend(nrow=2, title.position="top")) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "bottom"
	)

g5 <- ggplot() +
	theme(
		panel.border = element_blank()
	)
 
g6 <- ggplot(immunedata) +
	geom_vline(xintercept=2015:2025, lty=3) +
	geom_line(aes(year+week/52, C, col=as.factor(immune)), lwd=0.8) +
	annotate("segment", x=2024, xend=2024, y=20, yend=18,
					 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
	annotate("text", x=2024, y=24, label="Beginning of\nthe observed\noutbreak",
					 family="Times") +
	scale_x_continuous("Year", expand=c(0, 0),
										 limits=c(2019.5, 2025.5),
										 breaks=2015:2024) +
	scale_y_continuous("Incidence proxy", expand=c(0, 0), limits=c(0, 30)) +
	scale_color_manual("Duration of immunity\n(years)", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
	guides(color=guide_legend(nrow=2, title.position="top")) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "bottom",
	)

gcomb <- ggarrange(g1, g3, g5, g2, g4, g6, heights = c(1, 1.5),
									 labels=c("A", "C", "", "B", "D", "E"))

ggsave("figure3_new.pdf", gcomb, width=12, height=8)
