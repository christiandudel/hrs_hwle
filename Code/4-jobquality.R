### Packages ###################################################################

  library(dtms)
  library(tidyverse)
  library(writexl)  


### Load data ##################################################################

  load("Data/hrs_edited.Rda")


### Analysis of phyiscal, stress, poverty ######################################

  # Descriptive 
  quality_wu <- hrs |> filter(stateboth=="working/unhealthy" & wave%in%9:15) |> 
    group_by(gender,race,education) |> 
    summarise(mean(physical),mean(stress),mean(poverty),mean(anybad))
  
  quality_wh <- hrs |> filter(stateboth=="working/healthy" & wave%in%9:15) |> 
    group_by(gender,race,education) |> 
    summarise(mean(physical),mean(stress),mean(poverty),mean(anybad))
  
  # Get data right for regression, only recent waves
  regdat <- hrs |> filter(wave%in%9:15 & 
                            race!="Other" &
                            stateboth %in% c("working/healthy","working/unhealthy"))
  
  # Edit variables a bit
  regdat$age2 <- regdat$age^2
  regdat$race <- as.factor(regdat$race)
  regdat$education <- as.factor(regdat$education)
  regdat$gender <- as.factor(regdat$gender)
  regdat$stateboth <- as.factor(regdat$stateboth)
  
  # Regression models
  physical_logit <- glm(physical ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                        family=binomial,
                        data=regdat)
  
  stress_logit <- glm(stress ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                      family=binomial,
                      data=regdat)
  
  poverty_logit <- glm(poverty ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                       family=binomial,
                       data=regdat)
  
  anybad_logit <- glm(anybad ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                       family=binomial,
                       data=regdat)
  
  # Random effects example, does unfortunately not converge
  # library(lme4)
  # physical_logit <- glmer(physical ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth + (1|id),
  #                         family=binomial(),
  #                         data=regdat)

  # Data frame for prediction
  regpre <- expand.grid(gender=c(1,2),
                        race=c("Black","Hispan","White"),
                        education=0:2,
                        stateboth=levels(regdat$stateboth),
                        age=63,
                        age2=63^2)
  
  regpre$race <- as.factor(regpre$race)
  regpre$education <- as.factor(regpre$education)
  regpre$gender <- as.factor(regpre$gender)
  regpre$stateboth <- as.factor(regpre$stateboth)
  
  # Predict outcomes
  regpre$physical <- predict(physical_logit,regpre,type="response")
  regpre$stress <- predict(stress_logit,regpre,type="response")
  regpre$poverty <- predict(poverty_logit,regpre,type="response")
  regpre$anybad <- predict(anybad_logit,regpre,type="response")
  
  # Format results
  tmp1 <- regpre |> filter(stateboth=="working/unhealthy")
  tmp2 <- regpre |> filter(stateboth=="working/healthy")
  riskratio <- tmp1[,7:10]/tmp2[,7:10]
  names(riskratio) <- paste0("rr_",names(riskratio))
  quality <- tmp1 |> select(gender,race,education,physical,stress,poverty,anybad)
  quality <- cbind(quality,riskratio)
  
  
### Bootstrap ##################################################################
  
  # Number of replications
  replications <- 1000
  
  # For results
  bootresults <- list()
  
  # List of individuals and rows
  ids <- unique(regdat$id)
  nids <- length(ids)
  regdat$nr <- 1:dim(regdat)[1]
  rows <- by(regdat$nr,regdat$id,function(x) return(x))
  
  # Loop for resampling
  for(rep in 1:replications) {
    
    # Output
    cat(".")
    if(rep%%50==0) cat("\n")
    
    # Get sample
    whichids <- sample(1:nids,size=nids,replace=T)
    whichrows <- unlist(as.list(rows)[whichids])
    
    # Get data
    tmpdat <- regdat[whichrows,]
    
    # Regression models
    physical_logit <- glm(physical ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                          family=binomial,
                          data=tmpdat)
    
    stress_logit <- glm(stress ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                        family=binomial,
                        data=tmpdat)
    
    poverty_logit <- glm(poverty ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                         family=binomial,
                         data=tmpdat)
    
    anybad_logit <- glm(anybad ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                        family=binomial,
                        data=tmpdat)
    
    # Predict outcomes
    regpre$physical <- predict(physical_logit,regpre,type="response")
    regpre$stress <- predict(stress_logit,regpre,type="response")
    regpre$poverty <- predict(poverty_logit,regpre,type="response")
    regpre$anybad <- predict(anybad_logit,regpre,type="response")
    
    # Format results
    tmp1 <- regpre |> filter(stateboth=="working/unhealthy")
    tmp2 <- regpre |> filter(stateboth=="working/healthy")
    tmpratio <- tmp1[,7:10]/tmp2[,7:10]
    
    # Result
    bootresults[[rep]] <- tmpratio
    
  }
  
  # Results
  rows <- nrow(riskratio)
  cols <- ncol(riskratio)
  lowratio <- matrix(data=NA,nrow=rows,ncol=cols)
  higratio <- matrix(data=NA,nrow=rows,ncol=cols)
  sigratio <- matrix(data=NA,nrow=rows,ncol=cols)
  
  for(i in 1:rows) {
    for(j in 1:cols) {
      tmp <- unlist(lapply(bootresults,function(x) x[i,j]))
      qtmp <- quantile(tmp,probs=c(0.025,0.975))
      lowratio[i,j] <- qtmp[1]
      higratio[i,j] <- qtmp[2]
      sigratio[i,j] <- as.numeric(!any(qtmp<0))
    }
  }
  
  
### Save results ###############################################################
  
  # Combine with bootstrap
  colnames(lowratio) <- c("rr_physical_95low","rr_stress_95low","rr_poverty_95low","rr_anybad_95low")
  colnames(higratio) <- c("rr_physical_95up","rr_stress_95up","rr_poverty_95up","rr_anybad_95up")
  quality <- cbind(quality,lowratio,higratio)

  # Save  
  write_xlsx(quality,"Results/quality.xlsx")
  
  
### Preliminary ################################################################
  
  # Are poor workers at higher risk of stress and strain?
  physical_pov <- glm(physical ~ poverty + poverty*gender + poverty*race + poverty*education + gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                        family=binomial,
                        data=regdat)
  
  stress_pov <- glm(stress ~ poverty + poverty*gender + poverty*race + poverty*education + gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                      family=binomial,
                    data=regdat)
  
  # Predict
  regpre_poor <- regpre_non <- regpre
  regpre_poor$poverty <- 1
  regpre_non$poverty <- 0
  
  regpre_poor$physical <- predict(physical_pov,regpre_poor,type="response")
  regpre_poor$stress <- predict(stress_pov,regpre_poor,type="response")
  
  regpre_non$physical <- predict(physical_pov,regpre_non,type="response")
  regpre_non$stress <- predict(stress_pov,regpre_non,type="response")
  
  # Quick look
  regpre_poor$physical-regpre_non$physical # Always higher for poor people
  regpre_poor$stress-regpre_non$stress # Always lower
  