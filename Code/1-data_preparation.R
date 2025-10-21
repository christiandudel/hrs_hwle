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
                       # Sum of mobility difficulties
                       starts_with("r")&ends_with("mobila"),
                       # Sum of large muscle difficulties
                       starts_with("r")&ends_with("lgmusa"),
                       # Labor force status
                       starts_with("r")&ends_with("lbrf")&!contains("inlbrf"),
                       # Current job requires physical effort
                       starts_with("r")&ends_with("jphys"),
                       # Current job involves lots of stress
                       starts_with("r")&ends_with("jstres"),
                       # Poverty
                       starts_with("h")&ends_with("inpov"),
                       # Weights
                       starts_with("r")&ends_with("wtresp")
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
  dim(hrs)
  hrs |> filter(is.na(education) | is.na(race)) |> count()
  hrs <- hrs |> filter(!is.na(education) & !is.na(race))


### Rename vars for easier reshaping below #####################################

  # Wave status
  hrs <- hrs |> rename_with(~paste0("r",1:15,"inw"),starts_with("inw"))
  
  # Age
  hrs <- hrs |> rename_with(~paste0("r",1:15,"age"),ends_with("agey_e"))
  
  # Poverty
  hrs <- hrs |> rename_with(~paste0("r",1:15,"inpov"),ends_with("inpov"))

  # Empty vars for reshaping later (required by reshape function)
  hrs$r1mobila <- NA
  hrs$r1lgmusa <- NA

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
  
  # Mobility, large muscle, both combined
  hrs <- hrs |> mutate(mobility=case_match(mobila,
                                      0~0,
                                      1:5~1),
                       muscle=case_match(lgmusa,
                                       0~0,
                                       1:4~1),
                       both=ifelse(mobility%in%1 | muscle %in%1, 1, 0))
  
  # Get missings right
  hrs <- hrs |> mutate(both=ifelse(is.na(mobility)|is.na(muscle), NA, both))
  
  
### Work & health (combined) ###################################################
  
  # Work and limitations
  hrs <- hrs |> mutate(workboth=NA,
                       workboth=ifelse(worksimple%in%"working" & both%in%0,"working/healthy",workboth),
                       workboth=ifelse(worksimple%in%"working" & both%in%1,"working/unhealthy",workboth),
                       workboth=ifelse(worksimple%in%"retired" & both%in%0,"retired/healthy",workboth),
                       workboth=ifelse(worksimple%in%"retired" & both%in%1,"retired/unhealthy",workboth),
                       workboth=ifelse(worksimple%in%"not working" & !is.na(both) ,"not working",workboth))
  
  # Missing person-years
  table(hrs$stateboth,useNA="always") |> prop.table()
  

### State variables (including death) ##########################################
  
  # State using limitation
  hrs <- hrs |> mutate(stateboth=NA,
                       stateboth=ifelse(iwstat==1,workboth,stateboth),
                       stateboth=ifelse(iwstat==5,"dead",stateboth))
  
  

### Recode physical, stress, poverty ###########################################
  
  # Physical (1=yes,0 =no)
  hrs <- hrs |> mutate(physical=NA,
                       physical=ifelse(jphys%in%1:2,1,physical),
                       physical=ifelse(jphys%in%3:4,0,physical))
  
  # Stress (1=yes, 0=no)
  hrs <- hrs |> mutate(stress=NA,
                       stress=ifelse(jstres%in%1:2,1,stress),
                       stress=ifelse(jstres%in%3:4,0,stress))
  
  # Poverty (1=yes, 0=no)
  hrs <- hrs |> mutate(poverty=NA,
                       poverty=ifelse(inpov%in%1,1,poverty),
                       poverty=ifelse(inpov%in%0,0,poverty))
  
  # Any (1=yes, 0=no)
  hrs <- hrs |> mutate(anybad=NA,
                       anybad=ifelse(physical%in%1 | stress%in%1 | poverty%in%1,1,anybad),
                       anybad=ifelse(physical%in%0 & stress%in%0 & poverty%in%0,0,anybad))

  
### Limit data #################################################################

  # Limit variables
  hrs <- hrs |> select(hhidpn,ragender,race,education,wave,age,stateboth,
                       stress,physical,poverty,anybad,wtresp)

  # Rename
  hrs <- hrs |> rename('gender'='ragender',
                       'id'='hhidpn',
                       'weight'='wtresp')


### Saving #####################################################################

  save(hrs,file="Data/hrs_edited.Rda")

