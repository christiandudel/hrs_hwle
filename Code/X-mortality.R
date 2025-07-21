### Packages ###################################################################

  library(dtms)
  library(tidyverse)


### Load Packages ##############################################################

  library(tidyverse)
  library(readstata13)


### Load data ##################################################################

  # rda file, makes reloading a lot faster
  rdafile <- "Data/hrs.Rda"
  
  if(!file.exists(rdafile)) { 
    
    # Data; can be obtained from https://hrs.isr.umich.edu
    dtafile <- "Data/randhrs1992_2020v2.dta"
    
    # Load
    hrs <- read.dta13(file=dtafile,
                      convert.factors=FALSE) 
    
    # Save
    save(hrs,file=rdafile)
    
  } else load(rdafile)


### Select variables ###########################################################

  # ID, gender, death/birth year
  hrs <- hrs |> select(hhidpn,ragender,radyear,rabyear,
                       # Wave status: Response indicator (1= in wave)
                       starts_with("inw"), 
                       # Interview status (5 & 6 = dead)
                       starts_with("r")&ends_with("iwstat"),
                       # Age in years at interview month
                       starts_with("r")&ends_with("agey_e")&!contains("respagey"),
                       # Weights
                       starts_with("r")&ends_with("wtresp")
  )


### Rename vars for easier reshaping below #####################################

  # Wave status
  hrs <- hrs |> rename_with(~paste0("r",1:15,"inw"),starts_with("inw"))
  
  # Age
  hrs <- hrs |> rename_with(~paste0("r",1:15,"age"),ends_with("agey_e"))
  
  # Change format of time varying variables (not a great solution, but works)
  hrsnames <- str_split_fixed(names(hrs),"r[[:digit:]]{1,2}",2)
  hrsnames <- apply(hrsnames,1,function(x) {paste0(x,collapse="")})
  hrsnumbers <- parse_number(names(hrs))
  hrswhich <- !is.na(hrsnumbers)
  hrsnames[hrswhich] <- paste(hrsnames[hrswhich],hrsnumbers[hrswhich],sep="_")
  names(hrs) <- hrsnames


### Reshape ####################################################################

  # Get names of longitudinal vars and their ordering right 
  repvars <- grepl("_",names(hrs))   
  repvars <- names(hrs)[repvars]
  repvars <- unique(unlist(lapply(strsplit(repvars,split="_"),function(x)x[1])))
  repvars <- paste(rep(repvars, each = length(1:15)), 1:15, sep = "_")
  
  # Reshape (pivot_longer is just not intuitive to me, sorry)
  hrs <- reshape(data=as.data.frame(hrs),
                 direction="long",
                 varying=repvars,
                 sep="_",
                 idvar="hhidpn",
                 #times=1:15,
                 timevar="wave")
  
  # Sort 
  hrs <- hrs |> arrange(hhidpn,wave)

  # Drop people after death, and when not (yet) in wave
  hrs <- hrs |> filter(iwstat%in%c(1,5))


### Age ########################################################################

  # Age is missing in the year of death, add
  hrs <- hrs |> mutate(age=ifelse(iwstat==5,radyear-rabyear,age))
  
  # Age is still missing for a few people with unknown birth year and/or unknown 
  # year of death; for the latter, we impute year of death as mid-interval,
  # and generate age based on that
  hrs <- hrs |> mutate(toedit=ifelse(is.na(radyear) & !is.na(rabyear) & iwstat==5 & is.na(age),1,0),
                       radyear=case_when(
                         toedit==1 & wave==2~1993,
                         toedit==1 & wave==3~1995,
                         toedit==1 & wave==4~1997,
                         toedit==1 & wave==5~1999,
                         toedit==1 & wave==6~2001,
                         toedit==1 & wave==7~2003,
                         toedit==1 & wave==8~2005,
                         toedit==1 & wave==9~2007,
                         toedit==1 & wave==10~2009,
                         toedit==1 & wave==11~2011,
                         toedit==1 & wave==12~2013,
                         toedit==1 & wave==13~2015,
                         toedit==1 & wave==14~2017,
                         toedit==1 & wave==15~2019,
                         .default=radyear
                       ),
                       age=ifelse(iwstat==5&is.na(age),radyear-rabyear,age))
  
  # Drop if age is missing (very few, not relevant)
  hrs <- hrs |> filter(!is.na(age))


### Alive or dead ##############################################################

  # Employment (slightly more detailed/simplified)
  hrs <- hrs |> mutate(state=NA,
                       state=ifelse(iwstat==1,"alive",state), 
                       state=ifelse(iwstat==5,"dead",state))


### Limit data #################################################################

  # Limit variables
  hrs <- hrs |> select(hhidpn,ragender,wave,age,state,wtresp)
  
  # Rename
  hrs <- hrs |> rename('gender'='ragender',
                       'id'='hhidpn',
                       'weight'='wtresp')


### Setup dtms #################################################################

  # Vector of transient states
  transient_states <- c("alive")
  
  # Note: variable step length
  hrsdtms <- dtms(transient=transient_states,
                  absorbing="dead",
                  timescale=seq(50,98,1),
                  timestep=1:3)


### Reshape and edit data some more ############################################    

  # Reshape 
  estdata <- hrs |> dtms_format(data=_,
                                dtms=hrsdtms,
                                idvar="id",
                                timevar="age",
                                statevar="state",
                                steplength=TRUE)
  
  # Clean
  estdata <- dtms_clean(data=estdata,dtms=hrsdtms)
  
  # Edit/add variables
  estdata$time2 <- estdata$time^2


### Subsets by gender ##########################################################

  men2010 <- estdata |> filter(gender==1 & wave==9)
  women2010 <- estdata |> filter(gender==2 & wave==9)
  
  men2012 <- estdata |> filter(gender==1 & wave==10)
  women2012 <- estdata |> filter(gender==2 & wave==10)
  
  men2014 <- estdata |> filter(gender==1 & wave==11)
  women2014 <- estdata |> filter(gender==2 & wave==11)
  
  men2016 <- estdata |> filter(gender==1 & wave==12)
  women2016 <- estdata |> filter(gender==2 & wave==12)
  
  men2018 <- estdata |> filter(gender==1 & wave==13)
  women2018 <- estdata |> filter(gender==2 & wave==13)
  
  men2020 <- estdata |> filter(gender ==1 & wave==14)
  women2020 <- estdata |> filter(gender==2 & wave==14)

    
### General settings ###########################################################
  
  # DTMS for prediction
  hrspredict <- dtms(transient=transient_states,
                     absorbing="dead",
                     timescale=seq(50,98,2))
  
  
### Function for analysis ######################################################
  
  le <- function(data) {
    
    # Outcome
    data[,"outcome"] <- as.numeric(data[,"to"]=="dead")
    
    # Covariate values for prediction
    datapred <- data.frame(time=seq(50,98,2),
                           time2=seq(50,98,2)^2,
                           steplength=2)
    
    # Model
    fit <- glm(outcome~time+time2+steplength,data=data,weights=weight)
    
    # Predict probabilities
    probs <- predict(fit,datapred)
    
    # Area under the survival curve
    res <- sum(cumprod(c(1,1-probs)))*2

    # Return
    return(res)
    
  }
  
  
################################################################################

  men2010_res <- le(data=men2010)
  women2010_res <- le(data=women2010)
  
  men2012_res <- le(data=men2012)
  women2012_res <- le(data=women2012)
  
  men2014_res <- le(data=men2014)
  women2014_res <- le(data=women2014)
  
  men2016_res <- le(data=men2016)
  women2016_res <- le(data=women2016)
  
  men2018_res <- le(data=men2018)
  women2018_res <- le(data=women2018)
  
  men2020_res <- le(data=men2020)
  women2020_res <- le(data=women2020)

  le_men <- c(men2010_res,
              men2012_res,
              men2014_res,
              men2016_res,
              men2018_res,
              men2020_res)
  
  le_women <- c(women2010_res,
                women2012_res,
                women2014_res,
                women2016_res,
                women2018_res,
                women2020_res)

  years <- seq(2010,2020,2)  

  plot(years,le_men,type="l",xlab="Year",ylab="LE at age 50",ylim=c(24,38))
  lines(years,le_women,lty=2)
  