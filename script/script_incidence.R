library(vroom)
library(tidyr)
library(dplyr)
library(lubridate)

data_ili_raw <- vroom("../data/ILINet.csv", delim =",", skip=1)

data_ili_new <- data_ili_raw %>%
	filter(`REGION TYPE`=="National", YEAR+WEEK/52 >= 2015) %>%
	dplyr::select(YEAR, WEEK, `% WEIGHTED ILI`) %>%
	rename(year=YEAR, week=WEEK) %>%
	mutate(
		week=ifelse(week==53, 52, week)
	) %>%
	group_by(year, week) %>%
	summarize(
		ili=mean(as.numeric(`% WEIGHTED ILI`))
	)

data_myco_raw <- read.csv("../data/Trend Mycoplasma pneumoniae Detection Rates 2024-06-30.csv")

data_myco_gather <- data_myco_raw %>%
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
		value=mean(value)*100
	) %>%
	filter(key=="US") %>%
	head(-1)

data_incidence <- merge(data_ili_new, data_myco_gather) %>%
	mutate(
		proxy=value*ili
	)
