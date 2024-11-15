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

shiftvec <- 0:3

shiftlist <- vector('list', length(shiftvec))

for (i in 1:length(shiftvec)) {
	shift <- shiftvec[i]
	
	if (shift==0) {
		npieff_shift <- npieff
	} else {
		npieff_shift <- c(tail(npieff, -52*shift), rep(1, 52*shift))
	}
	
	npi1 <- data.frame(
		year=c(2014, rep(2015:2300, each=52)[1:(N+Npred-1)]),
		week=week,
		npieff=npieff_shift
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
			shift=shift,
			npieff=npieff_shift
		)
	
	shiftlist[[i]] <- out
}

shiftdata <- shiftlist %>%
	bind_rows %>%
	mutate(
		group=shift,
		group=factor(group,
								 levels=0:3,
								 labels=c("Inferred NPI effects",
								 				 "1 year earlier",
								 				 "2 years earlier",
								 				 "3 years earlier"))
	)

g1 <- ggplot(shiftdata %>% filter(C != 0)) +
	geom_vline(xintercept=(-5):6, lty=3) +
	geom_line(aes(year+week/52+shift-2020, S/20000), alpha=0.5, col="#E02938") +
	geom_line(aes(year+week/52+shift-2020, C), lwd=0.8) +
	annotate("segment", x=2024-2020, xend=2024-2020, y=20, yend=18,
					 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
	annotate("text", x=2024-2020, y=24, label="Beginning of\nthe observed\noutbreak",
					 family="Times") +
	scale_x_continuous("Time since NPI introduction (years)", expand=c(0, 0),
										 limits=c(-5, 6),
										 breaks=c(-4, -2, 0, 2, 4, 6)) +
	scale_y_continuous("Incidence proxy", expand=c(0, 0), limits=c(0, 30),
										 sec.axis = sec_axis(~.*20000/1e6, name="Susceptible proportion")) +
	facet_wrap(~group, ncol=2) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		legend.position = "bottom",
		axis.line.y.right = element_line(color="#E02938"),
		axis.ticks.y.right = element_line(color="#E02938"),
		axis.text.y.right = element_text(color="#E02938"),
		axis.title.y.right = element_text(color="#E02938")
	)

ggsave("figure3_cyclic.pdf", g1, width=12, height=8)
