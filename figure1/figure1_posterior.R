library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2); theme_set(theme_bw(base_family="Times", base_size = 13))
library(rstan)
library(egg)
load("../stanfit_sirs/stanfit_sirs_us_incidence_npi.rda")

ee <- rstan::extract(stanfit_sirs_us_incidence_npi, permuted=FALSE)

S0 <- data.frame(
	value=c(ee[,,"S0"])
)

I0 <- data.frame(
	value=c(ee[,,"I0"])
)

rho <- data.frame(
	value=c(ee[,,"rho"])
)

tau <- data.frame(
	value=c(ee[,,"tau"])
)

sigma <- data.frame(
	value=c(ee[,,"sigma"])
)

sigma_obs <- data.frame(
	value=c(ee[,,"sigma_obs"])
)

tt <- theme(
	panel.grid = element_blank(),
	panel.border = element_blank(),
	axis.line = element_line(),
	legend.position = "none"
)

g1 <- ggplot(S0) +
	stat_function(fun=function(x) dunif(x, 0, 1), col="black", lwd=0.7) +
	geom_density(aes(value), lwd=0.7, col="#E02938") +
	scale_x_continuous("Initial proportion susceptible", limits=c(0, 1)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

g2 <- ggplot(I0) +
	stat_function(fun=function(x) dnorm(x, 0, 0.01)*2, col="black", lwd=0.7) +
	geom_density(aes(value), lwd=0.7, col="#E02938") +
	scale_x_continuous("Initial proportion infected", limits=c(0, 0.02)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

g3 <- ggplot(rho) +
	stat_function(fun=function(x) dnorm(x, 0, 2)*2, col="black", lwd=0.7) +
	geom_density(aes(value), lwd=0.7, col="#E02938") +
	scale_x_continuous("Scaling parameter", limits=c(0, 0.01)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

g4 <- ggplot(tau) +
	stat_function(fun=function(x) dnorm(x, 400/52, 200/52), col="black", lwd=0.7) +
	geom_density(aes(value/52), lwd=0.7, col="#E02938") +
	scale_x_continuous("Mean duration of immunity (years)", limits=c(0, 1500/52)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

g5 <- ggplot(sigma) +
	stat_function(fun=function(x) dnorm(x, 0, 0.2)*2, col="black", lwd=0.7) +
	geom_density(aes(value), lwd=0.7, col="#E02938") +
	scale_x_continuous("Standard deviation in transmission rate", limits=c(0, 0.4)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

g6 <- ggplot(sigma_obs) +
	stat_function(fun=function(x) dnorm(x, 0, 0.5)*2, col="black", lwd=0.7) +
	geom_density(aes(value), lwd=0.7, col="#E02938") +
	scale_x_continuous("Standard deviation in residuals", limits=c(0, 2)) +
	scale_y_continuous("Probability density", expand=c(0, 0)) +
	tt

gcomb <- ggarrange(g1, g2, g3, g4, g5, g6,
					nrow=2,
					labels=c("A", "B", "C", "D", "E", "F"))

ggsave("figure1_posterior.pdf", gcomb, width=12, height=6)
