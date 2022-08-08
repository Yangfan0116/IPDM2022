###########################################################
####          Halifax Rabies Detection and Vacc       #####
###########################################################

# What's new?
# Detection code
# Vaccinate code
# Detection is based on a number of dead dogs within a week (can be adjusted)
# Vaccination is delayed to start after detection (delay can be adjusted)

library(RGeode)

#set.seed(250)

# Parameters -------------------------------------

# Probability of infection:
# What is beta made of?
# Probability of effective contact
# Probability of bite
# Probability of infection

DensityCoeff = 0.2  # If contact was based on density, how many more contacts would the dog have for every extra dog introduced to the area?
FrequencyContacts = 2  # If contact was based on frequency, what is the average number of dog contacts per dog, regardless of density?
Pbite = 0.3 # could include this is 'effective contact'
PInf = 0.5


# Duration of latent period
LPeriod_meanlog = 3.0
LPeriod_sdlog = 0.6
hist(rlnorm (1000, LPeriod_meanlog, LPeriod_sdlog), breaks = 30)


# Duration of infectious period
IPeriod_shape = 2.49
IPeriod_rate = 0.81
hist(rgamma (1000, IPeriod_shape, IPeriod_rate), breaks = 30) 

# ProbVacc
ProbVacc = 0.2


# We want to simulate 5 years:
end.time <- 1.5 * 365




# Community ---------------------------------------
n.dogs <- 220
# Create the community:
community <- data.frame(id =1:n.dogs,
                   age=round(runif(n.dogs,90,3650)),
                   state=0, 
                   latent.days = round(rlnorm (n.dogs, LPeriod_meanlog, LPeriod_sdlog)),
                   infected.days=round(rgammatr (n.dogs, IPeriod_shape, IPeriod_rate, range = c(2, 14))),
                   latent.count = 0,
                   infection.count = 0)

# We now have four compartments in the "infected" column:
# 0 = Susceptible
# 1 = latent
# 2 = Infectious
# 3 = Recovered

# Initial states and triggers-------------------------------------
# We start with one latently infected dog
community$state[1] <- 2

Vaccinate = 0 # Vaccination switch

VaccTrigger = 4 # Number of dead dogs in a week that somebody would be likely to notice
VaccTrigDay = 0 #Leave this as zero
VaccDelay = 14 # Alter to test impact of vaccine delays

# Collections -----------------------------------------

# Collect the :
age_collect <- numeric(end.time)
# Collect the number of susceptible and infected animals at each time point:
inf_collect <- data.frame(susceptible=rep(0,end.time), latent=rep(0, end.time), infected=rep(0, end.time), recovered=rep(0, end.time), vaccinated=rep(0, end.time))

DailyDead_collect = numeric(end.time)
VacTrig_collect = numeric(end.time) #rolling previous weekly sum for each day


# Model -----------------------------------------------

for (k in 1:end.time)
{
  # We update the ProbInf within the model:
  
  # Density transmission
  # risk of infection = lambda
  # lambda = beta * I
  # lambda = c * v * I/N   c = prob contact (k* N/A), v = prob infection, I/N = probability that a given contact is with an infected individual 
  # lambda = k * v * I    K is the coefficient of the slope
  
  ProbInfectionBeat_FD <- DensityCoeff * Pbite * PInf   
  
  # ProbInfection <- 1-exp(-(ProbInfectionBeat_FD * (length(community$state[community$state==2]))))
  # ProbInfection <- 1- (1- ProbInfectionBeat_FD)^ (length(community$state[community$state==2]))  
  
  
  # Frequency transmission
  # risk of infection = lambda
  # lambda = beta' * I
  # lambda = c' * v * I/N   c = prob contact (n + 0*N/A), v = prob infection, I/N = probability that a given contact is with an infected individual 
  # lambda = n * v * I/N    n is the number of contacts/time step
  
  #ProbInfection <- FrequencyContacts * Pbite * PInf * (length(community$state[community$state==2])/n.dogs)    
  ProbInfection <- 1- (1- (FrequencyContacts * Pbite * PInf))^ (length(community$state[community$state==2])/n.dogs)  
  
  #ProbInfection <- 1-exp(-((ProbInfectionBeta_FD *(length(community$state[community$state==2])) / n.dogs)))
  ProbInfection <- 1- (1- ProbInfectionBeta_FD)^ (length(community$state[community$state==2])/n.dogs)
  
  
  ### Recovery ###
  
  # Identify the dogs that will 'recover' on this day:
  PotRec       <- which(community$state==2)
  # Draw randomly which of the dogs should recover:
  NewRec       <- PotRec[which(community$infected.days[PotRec]<=community$infection.count[PotRec])]
  # If there are any recoveries, they recover:
  if(length(NewRec)>0) community$state[NewRec] <- 3
  # And update the infected days counter (at each time step):
  community$infection.count[community$state==2] <- community$infection.count[community$state==2] + 1
  
  
  ### Infectious ###
  
  # Identify the dogs that will become infectious on this day:
  PotInf       <- which(community$state==1)
  # Draw randomly which of the dogs should recover:
  NewInf       <- PotInf[which(community$latent.days[PotInf]<=community$latent.count[PotInf])]
  # If there are any recoveries, they recover:
  if(length(NewInf)>0) community$state[NewInf] <- 2
  # And update the infected days counter (at each time step):
  community$latent.count[community$state==1] <- community$latent.count[community$state==1] + 1

  
  ### Latent Infection ###
  
  # Identify the dogs that will become latently infected on this day:
  PotLat       <- which(community$state==0)
  # Draw randomly which of the cows that should be (newly) infected:
  NewLat       <- PotLat[rbinom(length(PotLat),1,prob=ProbInfection)==1]
  # If there are new infections, then we infect them:
  if(length(NewLat)>0) community$state[NewLat] <- 1
 
  
  # Add one day to the age of all the animals, for each simulated:
  community$age <- community$age + 1
  
  # Save the daily mean age of all dogs:
  age_collect[k] <- mean(community$age)

  DailyDead_collect[k] <- length(NewRec)

  if((k)>7 & Vaccinate==0) {
    VacTrig_collect[k] = sum(DailyDead_collect[(k-7):k])
  
      if(VacTrig_collect[k]>VaccTrigger) {
        Vaccinate = 1
        VaccTrigDay = k
        }
    }
  
  if ((Vaccinate) == 1 & (k>=VaccTrigDay+VaccDelay)) {
    PotVacc       <- which(community$state<=1)
    # Draw randomly which of the dogs should be vaccinated:
    NewVacc       <- PotVacc[rbinom(length(PotVacc),1,prob=ProbVacc)==1]
    # If there are any to be vaccinated they change state to 5:
    if(length(NewVacc)>0) community$state[NewVacc] <- 5
  }
      
  # Save the number of susceptible, infected and recovered:
  sus <- length(community$state[community$state==0])
  lat <- length(community$state[community$state==1])
  inf <- length(community$state[community$state==2])
  rec <- length(community$state[community$state==3])
  vac <- length(community$state[community$state==5])
  inf_collect[k,] <- c(sus, lat, inf, rec, vac)
  

}


plot(inf_collect$susceptible, type="l", lwd=2, col=1, xlab="Time", ylab="No. individuals", ylim = c(0, n.dogs))
lines(inf_collect$latent, type="l", lwd=2, col=2)
lines(inf_collect$infected, type="l", lwd=2, col=3)
lines(inf_collect$recovered, type="l", lwd=2, col=4)
lines(inf_collect$vaccinated, type="l", lwd=2, col=5)
legend("right", text.col=c(1, 2, 3, 4, 5), legend=c("Susceptible", "Latent", "Infected", "Recovered", "Vaccinated"), 
       bty="n")
abline(v = VaccTrigDay)


