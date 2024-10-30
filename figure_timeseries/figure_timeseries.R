library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(rstan)
library(egg)
source("../script/script_incidence.R")

g1 <- ggplot(data_myco_gather) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_point(aes(year+week/52, value)) +
	geom_line(aes(year+week/52, value)) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024) +
	scale_y_continuous("Test positivity (%)", limits=c(0, 3.2), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		axis.title.x = element_blank()
	)

g2 <- ggplot(data_ili_new) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_point(aes(year+week/52, ili)) +
	geom_line(aes(year+week/52, ili)) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024) +
	scale_y_continuous("Weighted ILI (%)", limits=c(0, 8), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1),
		axis.title.x = element_blank()
	)

g3 <- ggplot(data_incidence) +
	geom_vline(xintercept=2015:2024, lty=3) +
	geom_point(aes(year+week/52, proxy)) +
	geom_line(aes(year+week/52, proxy)) +
	scale_x_continuous("Year", expand=c(0, 0),
										 breaks=2015:2024) +
	scale_y_continuous("Incidence proxy", limits=c(0, 8), expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank(),
		panel.border = element_rect(linewidth=1)
	)

gcomb <- ggarrange(g1, g2, g3, nrow=3,
									 labels=c("A", "B", "C"))

ggsave("figure_timeseries.pdf", gcomb, width=8, height=6)
