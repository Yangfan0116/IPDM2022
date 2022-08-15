#set.seed(100)
n_turkeys = 1000

farm<-data.frame(id =1:n_turkeys,
                 infected=0,
                 incubationdaycount=0,
                 infecteddaycount=0,
                 totalincubationday=0,
                 totalinfectedday=0)

#initial 2800 infected
farm$infected[1] <- 2

#beta
transmissionrate <- rgamma(1,1.5,1.5) #* nrow(farm[farm$infected == 2])/n_turkeys
  
  
#INITAL SAMPLED FOR TOTAL DAYS FOR EACH ANIMAL
farm$totalincubationday<-round(rgamma(n_turkeys,3,1))
farm$totalinfectedday<- round(rgamma(n_turkeys,3,1))

end_time <- 25


inf_collect <- data.frame(susceptible=rep(0,end_time), latent=rep(0,end_time), infected=rep(0, end_time), removed=rep(0, end_time))

for (t in seq_len(end_time))
{
  
  
  # Save the number of susceptible and infected:
  sus <- length(farm$infected[farm$infected==0])
  lat <- length(farm$infected[farm$infected==1])
  inf <- length(farm$infected[farm$infected==2])
  rem <- length(farm$infected[farm$infected==3])
  inf_collect[t,1] <- sus
  inf_collect[t,2] <- lat
  inf_collect[t,3] <- inf
  inf_collect[t,4] <- rem
  
  #beta
  
  # Identify the turkey that are not infected (susceptible):
  PotInf<- which(farm$infected==0)
  # Draw randomly which of the turkey that should be (newly) latent:
  NewInf<- PotInf[rbinom(length(PotInf),1,prob=transmissionrate)==1]
  # If there are new infections, then we infect them:
  if(length(NewInf)>0) farm$infected[NewInf] <- 1
  
  
  # And update the infected days counter and incubation day (at each time step):
  farm$incubationdaycount[farm$infected==1] <- farm$incubationdaycount[farm$infected==1] + 1
  farm$infecteddaycount[farm$infected==2] <- farm$infecteddaycount[farm$infected==2] + 1
  
  # if the counts exceed the total counts then move them to the next state
  farm$infected[farm$incubationdaycount>=farm$totalincubationday] <- 2
  farm$infected[farm$infecteddaycount>=farm$totalinfectedday] <- 3

}


plot(inf_collect$susceptible, type="l", lwd=2, col=3, xlab="Time", ylab="No. individuals", ylim=c(0,n_turkeys))
lines(inf_collect$latent, type="l", lwd=2, col=2)
lines(inf_collect$infected, type="l", lwd=2, col=5)
lines(inf_collect$removed, type="l", lwd=2, col=4)
legend("right", text.col=c(3,2,5,4), legend=c("Susceptible", "Latent", "Infected", "Removed"),
       bty="n")

