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

  # ID, gender, death/birth year, education
  hrs <- hrs |> select(hhidpn,ragender,radyear,rabyear,raeduc,rahispan,raracem,
                       # Wave status: Response indicator (1= in wave)
                       starts_with("inw"), 
                       # Interview status (5 & 6 = dead)
                       starts_with("r")&ends_with("iwstat"),
                       # Age in years at interview month
                       starts_with("r")&ends_with("agey_e")&!contains("respagey"),
                       # Sum of (i)ADL responses
                       starts_with("r")&ends_with("adl5a"),
                       # Labor force status
                       starts_with("r")&ends_with("lbrf")&!contains("inlbrf")
                       )


### Education/race #############################################################

  # Education, 3 levels, 0=low, 1=medium, 2=high
  hrs <- hrs |> mutate(education=case_match(raeduc,
                                            c(1,2)~0,
                                            c(3,4)~1,
                                            5~2))
  
  # Race recode
  hrs <- hrs |> mutate(race=NA) |> 
    mutate(race=ifelse(raracem%in%1 & rahispan%in%0,"White",race),
           race=ifelse(raracem%in%2 & rahispan%in%0,"Black",race),
           race=ifelse(rahispan%in%1,"Hispan",race),
           race=ifelse(raracem%in%3 & rahispan%in%0,"Other",race),
           race=ifelse(raracem%in%3 & is.na(rahispan),"Other",race))
  
  
  # Drop if education is missing (22 individuals, negligible)
  hrs <- hrs |> filter(!is.na(education) & !is.na(race))


### Rename vars for easier reshaping below #####################################

  # Wave status
  hrs <- hrs |> rename_with(~paste0("r",1:15,"inw"),starts_with("inw"))
  
  # Age
  hrs <- hrs |> rename_with(~paste0("r",1:15,"age"),ends_with("agey_e"))

  # Empty vars for reshaping later (required by reshape function)
  hrs$r1adl5a <- NA
  hrs$r1iadl5a <- NA

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
  
  # Drop if age is missing
  hrs <- hrs |> filter(!is.na(age))


### Work, disability, health (separate) ########################################
  
  # Pension age
  hrs <- hrs |> mutate(state_pension=ifelse(rabyear<=1942,65,
                                           ifelse(rabyear>1942 & rabyear<=1959,66,
                                                   ifelse(rabyear>=1960,67,NA))))
  
  
  # Employment (slightly more detailed/simplified)
  hrs <- hrs |> mutate(workstatus=case_match(
    lbrf,
    1:2~"working",
    3  ~"unemployed",
    4:5~"retired",
    6:7~"inactive"),
    workstatus=ifelse(lbrf%in%6:7&age>=state_pension,"retired",workstatus),
    worksimple=case_match(
      workstatus,
      c("unemployed","inactive")~"not working",
      .default=workstatus)) 
  
  # ADL, iADL, self-rated health (1=unhealthy/disabled)
  hrs <- hrs |> mutate(adl=case_match(adl5a,
                                      0~0,
                                      1:5~1),
                       iadl=case_match(iadl5a,
                                       0~0,
                                       1:5~1))
  
  
### Work & disability/health (combined) ########################################
  
  # Work and ADL
  hrs <- hrs |> mutate(workadl=NA,
                       workadl=ifelse(worksimple=="working" & adl==0,"working/healthy",workadl),
                       workadl=ifelse(worksimple=="working" & adl==1,"working/unhealthy",workadl),
                       workadl=ifelse(worksimple=="retired" & adl==0,"retired/healthy",workadl),
                       workadl=ifelse(worksimple=="retired" & adl==1,"retired/unhealthy",workadl),
                       workadl=ifelse(worksimple=="not working" ,"not working",workadl))
  
  # Work and iADL
  hrs <- hrs |> mutate(workiadl=NA,
                       workiadl=ifelse(worksimple=="working" & iadl==0,"working/healthy",workiadl),
                       workiadl=ifelse(worksimple=="working" & iadl==1,"working/unhealthy",workiadl),
                       workiadl=ifelse(worksimple=="retired" & iadl==0,"retired/healthy",workiadl),
                       workiadl=ifelse(worksimple=="retired" & iadl==1,"retired/unhealthy",workiadl),
                       workiadl=ifelse(worksimple=="not working" ,"not working",workiadl))

  
### State variables (including death) ##########################################
  
  # State using ADL
  hrs <- hrs |> mutate(state_adl=NA,
                       state_adl=ifelse(iwstat==1,workadl,state_adl),
                       state_adl=ifelse(iwstat==5,"dead",state_adl))
  
  # State using iADL
  hrs <- hrs |> mutate(state_iadl=NA,
                       state_iadl=ifelse(iwstat==1,workiadl,state_iadl),
                       state_iadl=ifelse(iwstat==5,"dead",state_iadl))

  
### Limit data #################################################################

  # Limit variables
  hrs <- hrs |> select(hhidpn,ragender,race,education,wave,age,state_adl,state_iadl)

  # Rename
  hrs <- hrs |> rename('gender'='ragender',
                       'id'='hhidpn')

  # Drop obs
  hrs <- na.omit(hrs)


### Saving #####################################################################

  save(hrs,file="Data/hrs_edited.Rda")

