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
  
  # Needs to adjust for health
  physical_logit <- glm(physical ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                        family=binomial,
                        data=regdat)
  
  stress_logit <- glm(stress ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                      family=binomial,
                      data=regdat)
  
  poverty_logit <- glm(poverty ~ gender*race + gender*education + race*education + gender*age + gender*age2 + gender*stateboth,
                       family=binomial,
                       data=regdat)
  
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
  
  # Format results
  tmp1 <- regpre |> filter(stateboth=="working/unhealthy")
  tmp2 <- regpre |> filter(stateboth=="working/healthy")
  riskratio <- tmp1[,7:9]/tmp2[,7:9]
  names(riskratio) <- paste0("rr_",names(riskratio))
  quality <- tmp1 |> select(gender,race,education,physical,stress,poverty)
  quality <- cbind(quality,riskratio)
  
  
### Save results ###############################################################
  
  write_xlsx(quality,"Results/quality.xlsx")
