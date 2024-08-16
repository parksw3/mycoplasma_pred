simulate_sirs <- function(beta, 
                          week,
                          mu,
													delta,
                          gamma,
                          rho,
                          npieff,
                          pop,
                          S0,
                          I0,
                          tmax=572) {
  Svec <- Ivec <- Rvec <- rep(0, tmax)
  
  Svec[1] <- S0 * pop
  Ivec[1] <- I0 * pop
  Rvec[1] <- (1 - S0 - I0) * pop
  
  for (i in 2:tmax) {
    foi <- beta[[week[i]]] * (Ivec[i-1])/pop * npieff[i]
    
    Sout <- (1-exp(-(foi+mu))) * Svec[i-1]
    StoI <- foi/(foi+mu) * Sout
    Iout <- (1-exp(-(gamma+mu))) * Ivec[i-1]
    ItoR <- gamma/(gamma+mu) * Iout
    Rout <- (1-exp(-(delta+mu))) * Rvec[i-1]
    RtoS <- delta/(delta+mu) * Rout
    
    Svec[i] <- Svec[i-1] - Sout + mu * pop + RtoS
    Ivec[i] <- Ivec[i-1] + StoI - Iout
    Rvec[i] <- Rvec[i-1] + ItoR - Rout
  }
  
  data.frame(
    S=Svec,
    I=Ivec,
    R=Rvec,
    Ifit=Ivec*rho
  )
}

simulate_sir_stoch <- function(beta, 
												 week,
												 mu,
												 gamma,
												 rho,
												 npieff,
												 pop,
												 S0,
												 I0,
												 tmax=572) {
	Svec <- Ivec <- Rvec <- rep(0, tmax)
	
	Svec[1] <- round(S0 * pop)
	Ivec[1] <- round(I0 * pop)
	Rvec[1] <- round((1 - S0 - I0) * pop)
	
	for (i in 2:tmax) {
		foi <- beta[[week[i]]] * (Ivec[i-1])/pop * npieff[i]
		
		Sout <- rbinom(1, Svec[i-1], (1-exp(-(foi+mu))))
		StoI <- rbinom(1, Sout, foi/(foi+mu))
		Iout <- rbinom(1, Ivec[i-1], (1-exp(-(gamma+mu))))
		ItoR <- rbinom(1, Iout, gamma/(gamma+mu))
		Rout <- rbinom(1, Rvec[i-1], (1-exp(-(mu))))
		birth <- rpois(1, mu * pop)
		
		Svec[i] <- Svec[i-1] - Sout + birth
		Ivec[i] <- Ivec[i-1] + StoI - Iout
		Rvec[i] <- Rvec[i-1] + ItoR - Rout
	}
	
	data.frame(
		S=Svec,
		I=Ivec,
		R=Rvec,
		Ifit=Ivec*rho
	)
}
