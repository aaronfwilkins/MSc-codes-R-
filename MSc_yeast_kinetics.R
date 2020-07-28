## MCMC Statistical Parameter Estimation for MSc Yeast Copper Isotope Fractionation

library(deSolve)

parameters <- c(r=0.2,K=200,a1=0.2125, a2=0.2245, a3=0.2649, b1=0.2340, b2=0.2790, b3=0.2117)
state <- c(Y=1.0, Ng=42.0, Nm=1.0, Nn=1.0, Nc=1.0)
times <- seq(0, 50, by = 0.1)

kinetics <- function(times,state,parameters){
  with(as.list(c(state,parameters)),{
    # Kinetic system of equations
    dY <- r*Y *(1-(Y/K))
    dNg <- (b1^2)*Nc - (a1^2)*Ng
    dNm <- (a2^2)*Nc - (b2^2)*Nm + (Nm/(Nm+Nn+Nc))*(r*Y *(1-(Y/K)))
    dNn <- (a3^2)*Nc - (b3^2)*Nn + (Nn/(Nm+Nn+Nc))*(r*Y *(1-(Y/K)))
    dNc <- (a1^2)*Ng + (b2^2)*Nm + (b3^2)*Nn - (((a3^2)+(a2^2)+(b1^2))*Nc) + (Nc/(Nm+Nn+Nc))*(r*Y *(1-(Y/K))) 
  # Return the rate of change
    return(list(c(dY,dNg,dNm,dNn,dNc)))
    })
}
# ode integration through Fortran method 'lsoda'
out <- lsoda(y=state,times=times,func=kinetics,parms=parameters)
head(out, n=5)

reshapedout <- array(c(out[,'Y'],out[,'Ng'],out[,'Nm'],out[,'Nn'],out[,'Nc']),dim = c((50/0.1)+1,5))

matplot(out[,'time'],reshapedout, type = "l", lty=1:1,lwd=2, xlab="Time (hours)",ylab = "Concentration (ng/ml)")


