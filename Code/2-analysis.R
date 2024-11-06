### Packages ###################################################################

  library(dtms)
  library(tidyverse)


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
                         statevar="state_adl",
                         steplength=TRUE)

  # Clean
  estdata <- dtms_clean(data=estdata,dtms=hrsdtms)

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


### General settings ###########################################################

  # Controls
  controls <- c("time","time2","education","race","steplength")  

  # DTMS for prediction
  hrspredict <- dtms(transient=transient_states,
                     absorbing="dead",
                     timescale=seq(50,98,2))

  # Covariate values for prediction
  reused <- list(time=seq(50,98,2),
                 time2=seq(50,98,2)^2,
                 steplength=2)

  # Race/ethnicity and education
  low_w <- c(reused,
           list(education=factor("0",levels=c("0","1","2")),
                race=factor("White",levels=levels(estdata$race))))
  
  med_w <- c(reused,
           list(education=factor("1",levels=c("0","1","2")),
                race=factor("White",levels=levels(estdata$race))))
  
  hig_w <- c(reused,
           list(education=factor("2",levels=c("0","1","2")),
                race=factor("White",levels=levels(estdata$race))))
  
  low_b <- c(reused,
             list(education=factor("0",levels=c("0","1","2")),
                  race=factor("Black",levels=levels(estdata$race))))
  
  med_b <- c(reused,
             list(education=factor("1",levels=c("0","1","2")),
                  race=factor("Black",levels=levels(estdata$race))))
  
  hig_b <- c(reused,
             list(education=factor("2",levels=c("0","1","2")),
                  race=factor("Black",levels=levels(estdata$race))))
  
  low_h <- c(reused,
             list(education=factor("0",levels=c("0","1","2")),
                  race=factor("Hispan",levels=levels(estdata$race))))
  
  med_h <- c(reused,
             list(education=factor("1",levels=c("0","1","2")),
                  race=factor("Hispan",levels=levels(estdata$race))))
  
  hig_h <- c(reused,
             list(education=factor("2",levels=c("0","1","2")),
                  race=factor("Hispan",levels=levels(estdata$race))))
  

### Function for analysis ######################################################

bootfun <- function(data,dtms) {
  
  cat(".")
  
  # Split data
  #first <- data |> filter(wave%in%1:8)
  secon <- data |> filter(wave%in%9:15)
  
  # Model
  #fit1 <- dtms_fit(data=first,controls=controls)
  fit2 <- dtms_fit(data=secon,controls=controls,package="mclogit")

  # Predict probabilities
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
  S2_low_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("White",levels=levels(estdata$race))))
  S2_med_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("White",levels=levels(estdata$race))))
  S2_hig_w <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("White",levels=levels(estdata$race))))
  
  S2_low_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Black",levels=levels(estdata$race))))
  S2_med_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Black",levels=levels(estdata$race))))
  S2_hig_b <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Black",levels=levels(estdata$race))))
  
  S2_low_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("0",levels=c("0","1","2")),race=factor("Hispan",levels=levels(estdata$race))))
  S2_med_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("1",levels=c("0","1","2")),race=factor("Hispan",levels=levels(estdata$race))))
  S2_hig_h <- dtms_start(dtms=dtms,data=secon,start_time=c(50:56),variables=list(education=factor("2",levels=c("0","1","2")),race=factor("Hispan",levels=levels(estdata$race))))
  
  # 1: Expectancies
  resexp <- rbind(
    data.frame(education=0,race="White",dtms_expectancy(dtms=dtms,matrix=Tm2_low_w,start_distr=S2_low_w))["AVERAGE",],
    data.frame(education=1,race="White",dtms_expectancy(dtms=dtms,matrix=Tm2_med_w,start_distr=S2_med_w))["AVERAGE",],
    data.frame(education=2,race="White",dtms_expectancy(dtms=dtms,matrix=Tm2_hig_w,start_distr=S2_hig_w))["AVERAGE",],
    data.frame(education=0,race="Black",dtms_expectancy(dtms=dtms,matrix=Tm2_low_b,start_distr=S2_low_b))["AVERAGE",],
    data.frame(education=1,race="Black",dtms_expectancy(dtms=dtms,matrix=Tm2_med_b,start_distr=S2_med_b))["AVERAGE",],
    data.frame(education=2,race="Black",dtms_expectancy(dtms=dtms,matrix=Tm2_hig_b,start_distr=S2_hig_b))["AVERAGE",],
    data.frame(education=0,race="Hispanic",dtms_expectancy(dtms=dtms,matrix=Tm2_low_h,start_distr=S2_low_h))["AVERAGE",],
    data.frame(education=1,race="Hispanic",dtms_expectancy(dtms=dtms,matrix=Tm2_med_h,start_distr=S2_med_h))["AVERAGE",],
    data.frame(education=2,race="Hispanic",dtms_expectancy(dtms=dtms,matrix=Tm2_hig_h,start_distr=S2_hig_h))["AVERAGE",]
    )

  # Rownames for now
  rownames(resexp) <- 1:dim(resexp)[1]

  # Return
  return(resexp)
  
}


### Main results ###############################################################  

  men_res <- bootfun(data=men,dtms=hrspredict)
  women_res <- bootfun(data=women,dtms=hrspredict)


### Save results ###############################################################

  save(list=c("men_res","women_res"),
       file="Results/first_results.Rda")
  
  library(writexl)
  
  write_xlsx(men_res, "Results/first_results_men.xlsx")
  write_xlsx(women_res, "Results/first_results_women.xlsx")
