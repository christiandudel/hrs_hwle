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
  estdata <- dtms_format(data=hrs,
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

  men <- estdata |> filter(gender==1)
  women <- estdata |> filter(gender==2)


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


### General settings ###########################################################

  # DTMS for prediction
  hrspredict <- dtms(transient=transient_states,
                     absorbing="dead",
                     timescale=seq(50,98,2))


### Function for analysis ######################################################

  bootfun <- function(data,dtms) {
    
    # Controls
    controls <- c("time","time2","education","race","steplength")  
    
    # Covariate values for prediction
    reused <- list(time=seq(50,98,2),
                   time2=seq(50,98,2)^2,
                   steplength=2)
    
    # Race/ethnicity and education
    low_w <- c(reused,
               list(education=factor("0",levels=c("0","1","2")),
                    race=factor("White",levels=levels(data$race))))
    
    med_w <- c(reused,
               list(education=factor("1",levels=c("0","1","2")),
                    race=factor("White",levels=levels(data$race))))
    
    hig_w <- c(reused,
               list(education=factor("2",levels=c("0","1","2")),
                    race=factor("White",levels=levels(data$race))))
    
    low_b <- c(reused,
               list(education=factor("0",levels=c("0","1","2")),
                    race=factor("Black",levels=levels(data$race))))
    
    med_b <- c(reused,
               list(education=factor("1",levels=c("0","1","2")),
                    race=factor("Black",levels=levels(data$race))))
    
    hig_b <- c(reused,
               list(education=factor("2",levels=c("0","1","2")),
                    race=factor("Black",levels=levels(data$race))))
    
    low_h <- c(reused,
               list(education=factor("0",levels=c("0","1","2")),
                    race=factor("Hispan",levels=levels(data$race))))
    
    med_h <- c(reused,
               list(education=factor("1",levels=c("0","1","2")),
                    race=factor("Hispan",levels=levels(data$race))))
    
    hig_h <- c(reused,
               list(education=factor("2",levels=c("0","1","2")),
                    race=factor("Hispan",levels=levels(data$race))))
    
    # Split data
    first <- data |> filter(wave%in%2:8)
    secon <- data |> filter(wave%in%9:15)
    
    # Model
    fit1 <- dtms_fit(data=first,controls=controls,package="mclogit")
    fit2 <- dtms_fit(data=secon,controls=controls,package="mclogit")
  
    # Predict probabilities
    probs1_low_w <- dtms_transitions(dtms=dtms,model=fit1,controls=low_w,se=F)
    probs1_med_w <- dtms_transitions(dtms=dtms,model=fit1,controls=med_w,se=F)
    probs1_hig_w <- dtms_transitions(dtms=dtms,model=fit1,controls=hig_w,se=F)
    
    probs1_low_b <- dtms_transitions(dtms=dtms,model=fit1,controls=low_b,se=F)
    probs1_med_b <- dtms_transitions(dtms=dtms,model=fit1,controls=med_b,se=F)
    probs1_hig_b <- dtms_transitions(dtms=dtms,model=fit1,controls=hig_b,se=F)
    
    probs1_low_h <- dtms_transitions(dtms=dtms,model=fit1,controls=low_h,se=F)
    probs1_med_h <- dtms_transitions(dtms=dtms,model=fit1,controls=med_h,se=F)
    probs1_hig_h <- dtms_transitions(dtms=dtms,model=fit1,controls=hig_h,se=F)
  
    probs2_low_w <- dtms_transitions(dtms=dtms,model=fit2,controls=low_w,se=F)
    probs2_med_w <- dtms_transitions(dtms=dtms,model=fit2,controls=med_w,se=F)
    probs2_hig_w <- dtms_transitions(dtms=dtms,model=fit2,controls=hig_w,se=F)
    
    probs2_low_b <- dtms_transitions(dtms=dtms,model=fit2,controls=low_b,se=F)
    probs2_med_b <- dtms_transitions(dtms=dtms,model=fit2,controls=med_b,se=F)
    probs2_hig_b <- dtms_transitions(dtms=dtms,model=fit2,controls=hig_b,se=F)
    
    probs2_low_h <- dtms_transitions(dtms=dtms,model=fit2,controls=low_h,se=F)
    probs2_med_h <- dtms_transitions(dtms=dtms,model=fit2,controls=med_h,se=F)
    probs2_hig_h <- dtms_transitions(dtms=dtms,model=fit2,controls=hig_h,se=F)
    
    # Transition matrices
    Tm1_low_w <- dtms_matrix(dtms=dtms,probs=probs1_low_w)
    Tm1_med_w <- dtms_matrix(dtms=dtms,probs=probs1_med_w)
    Tm1_hig_w <- dtms_matrix(dtms=dtms,probs=probs1_hig_w)
    
    Tm1_low_b <- dtms_matrix(dtms=dtms,probs=probs1_low_b)
    Tm1_med_b <- dtms_matrix(dtms=dtms,probs=probs1_med_b)
    Tm1_hig_b <- dtms_matrix(dtms=dtms,probs=probs1_hig_b)
    
    Tm1_low_h <- dtms_matrix(dtms=dtms,probs=probs1_low_h)
    Tm1_med_h <- dtms_matrix(dtms=dtms,probs=probs1_med_h)
    Tm1_hig_h <- dtms_matrix(dtms=dtms,probs=probs1_hig_h)
    
    Tm2_low_w <- dtms_matrix(dtms=dtms,probs=probs2_low_w)
    Tm2_med_w <- dtms_matrix(dtms=dtms,probs=probs2_med_w)
    Tm2_hig_w <- dtms_matrix(dtms=dtms,probs=probs2_hig_w)
    
    Tm2_low_b <- dtms_matrix(dtms=dtms,probs=probs2_low_b)
    Tm2_med_b <- dtms_matrix(dtms=dtms,probs=probs2_med_b)
    Tm2_hig_b <- dtms_matrix(dtms=dtms,probs=probs2_hig_b)
    
    Tm2_low_h <- dtms_matrix(dtms=dtms,probs=probs2_low_h)
    Tm2_med_h <- dtms_matrix(dtms=dtms,probs=probs2_med_h)
    Tm2_hig_h <- dtms_matrix(dtms=dtms,probs=probs2_hig_h)
    
    # Starting distribution (correct, but small sample sizes) 
    S1_low_w <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    S1_med_w <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    S1_hig_w <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    
    S1_low_b <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    S1_med_b <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    S1_hig_b <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    
    S1_low_h <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
    S1_med_h <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
    S1_hig_h <- dtms_start(dtms=dtms,data=first,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
  
    S2_low_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    S2_med_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    S2_hig_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("White",levels=levels(data$race))))
    
    S2_low_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    S2_med_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    S2_hig_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Black",levels=levels(data$race))))
    
    S2_low_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
    S2_med_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
    S2_hig_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Hispan",levels=levels(data$race))))
    
    # Codes for below:
    # period = 1 : waves 2 to 8 
    # period = 2 : waves 9 to 15
    # education = 0 : low education
    # education = 1 : medium education
    # education = 2 : high education
    # race = 0 : white
    # race = 1 : black
    # race = 2 : hispanic
    
    # 1: Expectancies
    resexp <- rbind(# Race, 0=white, 1=black, 2=hispanic
      data.frame(period=1,education=0,race=0,dtms_expectancy(dtms=dtms,matrix=Tm1_low_w,start_distr=S1_low_w))["AVERAGE",],
      data.frame(period=1,education=1,race=0,dtms_expectancy(dtms=dtms,matrix=Tm1_med_w,start_distr=S1_med_w))["AVERAGE",],
      data.frame(period=1,education=2,race=0,dtms_expectancy(dtms=dtms,matrix=Tm1_hig_w,start_distr=S1_hig_w))["AVERAGE",],
      data.frame(period=1,education=0,race=1,dtms_expectancy(dtms=dtms,matrix=Tm1_low_b,start_distr=S1_low_b))["AVERAGE",],
      data.frame(period=1,education=1,race=1,dtms_expectancy(dtms=dtms,matrix=Tm1_med_b,start_distr=S1_med_b))["AVERAGE",],
      data.frame(period=1,education=2,race=1,dtms_expectancy(dtms=dtms,matrix=Tm1_hig_b,start_distr=S1_hig_b))["AVERAGE",],
      data.frame(period=1,education=0,race=2,dtms_expectancy(dtms=dtms,matrix=Tm1_low_h,start_distr=S1_low_h))["AVERAGE",],
      data.frame(period=1,education=1,race=2,dtms_expectancy(dtms=dtms,matrix=Tm1_med_h,start_distr=S1_med_h))["AVERAGE",],
      data.frame(period=1,education=2,race=2,dtms_expectancy(dtms=dtms,matrix=Tm1_hig_h,start_distr=S1_hig_h))["AVERAGE",],
      data.frame(period=2,education=0,race=0,dtms_expectancy(dtms=dtms,matrix=Tm2_low_w,start_distr=S2_low_w))["AVERAGE",],
      data.frame(period=2,education=1,race=0,dtms_expectancy(dtms=dtms,matrix=Tm2_med_w,start_distr=S2_med_w))["AVERAGE",],
      data.frame(period=2,education=2,race=0,dtms_expectancy(dtms=dtms,matrix=Tm2_hig_w,start_distr=S2_hig_w))["AVERAGE",],
      data.frame(period=2,education=0,race=1,dtms_expectancy(dtms=dtms,matrix=Tm2_low_b,start_distr=S2_low_b))["AVERAGE",],
      data.frame(period=2,education=1,race=1,dtms_expectancy(dtms=dtms,matrix=Tm2_med_b,start_distr=S2_med_b))["AVERAGE",],
      data.frame(period=2,education=2,race=1,dtms_expectancy(dtms=dtms,matrix=Tm2_hig_b,start_distr=S2_hig_b))["AVERAGE",],
      data.frame(period=2,education=0,race=2,dtms_expectancy(dtms=dtms,matrix=Tm2_low_h,start_distr=S2_low_h))["AVERAGE",],
      data.frame(period=2,education=1,race=2,dtms_expectancy(dtms=dtms,matrix=Tm2_med_h,start_distr=S2_med_h))["AVERAGE",],
      data.frame(period=2,education=2,race=2,dtms_expectancy(dtms=dtms,matrix=Tm2_hig_h,start_distr=S2_hig_h))["AVERAGE",]
      )
  
    # Rownames for now
    rownames(resexp) <- 1:dim(resexp)[1]
    
    # Type
    resexp <- as.matrix(resexp)
  
    # Return
    return(resexp)
    
  }


### Main results ###############################################################  

  men_res <- bootfun(data=men,dtms=hrspredict)
  women_res <- bootfun(data=women,dtms=hrspredict)


### Bootstrap ##################################################################  
  
  men_boot <- dtms_boot(data=men,
                        fun=bootfun,
                        dtms=hrspredict,
                        rep=500,
                        method="block",
                        parallel=TRUE,
                        cores=10,
                        .packages=c("mclogit","VGAM","nnet","dtms","dplyr"))
  
  women_boot <- dtms_boot(data=women,
                          fun=bootfun,
                          dtms=hrspredict,
                          rep=10,
                          method="block",
                          parallel=TRUE,
                          cores=10,
                          .packages=c("mclogit","VGAM","nnet","dtms","dplyr"))
  
  # If no parallel computing:
  # men_boot <- dtms_boot(data=men,
  #                       fun=bootfun,
  #                       dtms=hrspredict,
  #                       rep=500,
  #                       method="block")
  # 
  # women_boot <- dtms_boot(data=women,
  #                         fun=bootfun,
  #                         dtms=hrspredict,
  #                         rep=500,
  #                         method="block")
  
  men_boot <- summary(men_boot)
  women_boot <- summary(women_boot)
  
  
### Save results ###############################################################

  save(list=c("men_res","women_res","men_boot","women_boot"),
       file="Results/results.Rda")
  
  # For saving as Excel
  menlist <- c(list(men_res),men_boot)
  womenlist <- c(list(women_res),women_boot)
  
  menlist <- lapply(menlist,as.data.frame)
  womenlist <- lapply(womenlist,as.data.frame)
  
  menlist[[3]]$race <- menlist[[2]]$race <- menlist[[1]]$race
  womenlist[[3]]$race <- womenlist[[2]]$race <- womenlist[[1]]$race

  write_xlsx(menlist, "Results/results_men.xlsx")
  write_xlsx(womenlist, "Results/results_women.xlsx")
  
  # If running on workstation
  if(Sys.info()["nodename"]%in%c("HYDRA01","HYDRA02","HYDRA11")) {rm(list=ls());gc()}
