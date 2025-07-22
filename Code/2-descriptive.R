### Packages ###################################################################

  library(dtms)
  library(tidyverse)
  library(writexl)  


### Load data ##################################################################

  load("Data/hrs_edited.Rda")


### Setup dtms #################################################################

  # Vector of transient states
  transient_states <- c("retired/healthy","retired/unhealthy",
                        "not working",
                        "working/healthy","working/unhealthy")

  # Note: variable step length
  hrsdtms <- dtms(transient=transient_states,
                  absorbing="dead",
                  timescale=seq(50,98,1),
                  timestep=1:3)


### Reshape and edit data some more ############################################    

  # Reshape 
  estdata <- hrs |> select(id,gender,race,education,wave,age,stateboth,weight) |> 
                    dtms_format(data=_,
                    dtms=hrsdtms,
                    idvar="id",
                    timevar="age",
                    statevar="stateboth",
                    steplength=TRUE)

  # Clean
  estdata <- dtms_clean(data=estdata,dtms=hrsdtms)
  
  # Drop "Other"
  estdata <- estdata |> filter(race!="Other")
  
  # Edit/add variables
  estdata$time2 <- estdata$time^2
  estdata$education <- as.factor(estdata$education)
  estdata$race <- as.factor(estdata$race)
  

### Subsets by gender ##########################################################

  # Main analysis
  men <- estdata |> filter(gender==1)
  women <- estdata |> filter(gender==2)
  
  # Without 2020 (note: 'wave!=14' and NOT 'wave!=15' because 'wave' here is
  # the starting wave and wave 15 is implied)
  men2020 <- men |> filter(wave!=14)
  women2020 <- women |> filter(wave!=14)


### Sample size ################################################################

  # No of transitions
  dim(men)[1]+dim(women)[1]
  
  # No of individuals
  nmen <- men |> pull(id) |> unique() |> length()
  nwomen <- women |> pull(id) |> unique() |> length()
  nmen+nwomen
  

### Descriptives ###############################################################
  
  # Sample composition: Individuals
  individuals <- estdata |> group_by(id)|> filter(time == min(time))
  
  # Average age in first wave
  mean(individuals$time)
  
  # Table with overall descriptives for individuals
  table_individuals <- data.frame(Variable=c(rep("Gender",2),
                                  rep("Race/ethnicity",3),
                                  rep("Education",3)))
  
  table_individuals$Value <- c("Men","Women",
                               "White","Black","Hispanic",
                               "Less than HS","HS","College or university")
  
  table1 <- table(individuals$gender)
  table2 <- table(factor(individuals$race,levels=c("White","Black","Hispan")))
  table3 <- table(individuals$education)
  
  table_individuals$Observations <- c(table1,table2,table3)
  table_individuals$Relative <- c(table1/sum(table1),table2/sum(table2),table3/sum(table3))
  
  write_xlsx(table_individuals, "Results/table_individuals.xlsx")
  
  # Table combined characteristics observations and by time
  estdata <- estdata |> mutate(Period=NA,
                               Period=ifelse(wave%in%2:8,1,Period),
                               Period=ifelse(wave%in%9:15,2,Period))
  
  table_observations <- estdata |> group_by(gender,Period,race,education) |> summarize(Count=n())

  table_observations$Period <- factor(table_observations$Period,levels=c(1,2),labels=c("1994-2006","2008-2020"))
  table_observations$gender <- factor(table_observations$gender,levels=c(1,2),labels=c("Men","Women"))
  table_observations$race <- factor(table_observations$race,levels=c("White","Black","Hispan"),labels=c("White","Black","Hispanic"))
  table_observations$education <- factor(table_observations$education,levels=c(0,1,2),labels=c("Low","Medium","High"))

  names(table_observations) <- c("Gender","Period","Race","Education","N")
  
  table_observations <-  arrange(table_observations,Gender,Period,Race,Education)
  totals <- table_observations |> group_by(Period,Gender) |> summarize(P=sum(N))
  
  table_observations <- table_observations |> left_join(totals)
  tabke_observations <- table_observations |> mutate(P=N/P)
  
  write_xlsx(tabke_observations, "Results/table_observations.xlsx")

  # Transitions: absolute, relative, raw transition probability
  table_transitions <- summary(estdata)
  write_xlsx(table_transitions, "Results/table_transitions.xlsx")


