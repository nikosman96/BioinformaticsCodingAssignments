### POPM*6950 Disease Modelling 
### Deterministic Model of Ebola virus

### Loading Required Packages ---- 
#install.packages("") 
library(deSolve)
library(ggplot2)
library(gridExtra)

### ORIGINAL MODEL FROM LITERATURE ----
seihr.model <- function (t, x, params) { 
  S <- x[1] # Susceptibles as the first element of x
  E <- x[2] # Exposed as the second element of x
  I <- x[3] # Infected
  H <- x[4] # Hospitalized 
  RI <- x[5] # Removed and infectious
  RB <- x[6] # Removed and buried
  RR <- x[7] # Removed and recovered 
  
  with(as.list(params), 
       { 
         dS <- -beta1*S*I - beta2*S*RI - beta3*S*H
         dE <- beta1*S*I + beta2*S*RI + beta3*S*H - delta*E
         dI <- delta*E - gamma1*I - psi*I
         dH <- psi*I - gamma2*H
         dRI <- rho1*gamma1*I - omega*RI 
         dRB <- omega*RI + rho2*gamma2*H 
         dRR <- (1-rho1)*gamma1*I + (1-rho2)*gamma2*H
         dx <- c(dS,dE,dI,dH,dRI,dRB,dRR) # Combine results into a single vector dx
         list(dx) # Return result as a list
       }
  )
}

times <- seq(0,365,by=5) #function seq returns a sequence
xstart <-c(S=4289922,E=0,I=49/4290000,RΙ=29/4290000,RB=0,RR=0,H=0,V=0) 
params.init <- c(N= 4290000,gamma1=1/0.0542,gamma2=1/0.174
                 ,beta1=1/0.376,beta2=1/0.135,beta3=1/0.163,delta=9,omega=1/0.325,psi=1/0.5,rho1=0.98,rho2=0.88) #initial parameters


############### Simulate a model trajectory and run period pre-vax  ###############

#saveoutput for the first time points
#start at initial conditions

seihr.output = as.data.frame(lsoda(xstart, times, seihr.model, params.init))  



############### Plot  ###############

seihr.model.plot <- ggplot() + 
  geom_line(size = 1.2, data=seir.output, aes(x=time, y=S, colour = "S")) +
  geom_line(size = 1.2, data=seir.output, aes(x=time, y=E, colour = "E")) +
  geom_line(size = 1.2, data=seir.output, aes(x=time, y=I, colour = "I")) +
  geom_line(size = 1.2, data=seir.output, aes(x=time, y=R, colour = "R")) +
  labs(y="Population") + ggtitle("SEIR compartments over latent period of 0.0001 days ")