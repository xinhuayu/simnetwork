####################################################################
# network simulation code to explore susceptibility (Masking effect); 
# survivor bias, incomplete control and reopening
# population heterogeneity, and long tail and second wave
# Author: Xinhua Yu
# Revised date: October 25, 2020
###################################################################

rm(list=ls())
gc()
options(scipen = 999)
options(digits = 3)

library(mgcv)
library(data.table)
library(ggplot2)

rootpath ="C:/Users/xinhuayu/Desktop/COVID19/survivorbias/"

################# DATA ANALYSIS SIMULATION EPIDEMIC #################
# we can vary control measure by changing contact patterns (randomly delete links)
# contact patterns can be changed by aging population only
# we can estimate masking effect by reducing infectivity and susceptibility
# we can also only protect aging people with mask and reducing contacts
###################################################################

# use largest population
# START HERE: always read original data for each simulation,
pop = readRDS(paste0(rootpath,"sim_network_pop50_all.RDS"))
edgenet0 = readRDS(paste0(rootpath,"sim_network_net50_all.RDS"))
popsize = nrow(pop)

# fix the network connection data: broader age group
edgenet0[,c("age1","age2"):=list(
                        fifelse(age1 %in% c(1,2),1,
                          fifelse(age1 %in% c(3,4),2,
                            fifelse(age1 %in% c(5,6),3,4))),
                        fifelse(age2 %in% c(1,2),1,
                          fifelse(age2 %in% c(3,4),2,
                            fifelse(age2 %in% c(5,6),3,4))))]

# the average degree of network connections may be too high
# the original simulation try to make connections as many as possible
# randomly break 10% connections to get a bit more sparse
edgenet0 = edgenet0[,pid1:=ifelse(rbinom(nrow(edgenet0),1,0.1)==1,0,pid1)][pid1!=0,]  
# ensure uniqueness of network connections
edgenet0 = unique(edgenet0,by=c("pid1","pid2"))

# randomly assign contact weight, this assume repeat contacts every day
# most of them are 1's and a few are two or more
edgenet0[,cweight:=rpois(nrow(edgenet0),1)][,cweight:=ifelse(cweight==0,1,cweight)]

# calculating degree distributions
neta = copy(edgenet0)
setkey(neta,"pid1")
dd = neta[,degree:=.N,by="pid1"][,.SD[.N],by="pid1"]

# histogram with labels
agelab = c("0-19","20-44","45-64","65+")
names(agelab)=c("1","2","3","4")
ggplot(dd,aes(x=degree))+geom_histogram(aes(y = ..density..))+
  facet_grid(~age1,labeller=labeller(age1=agelab),scales = "free")+theme_bw() 

summary(dd$degree)
quantile(dd$degree,c(0,0.01,0.25,0.50,0.75,0.90,0.95,0.99,1))


#############################################
# outside contact matrix, eight small every ten year age groups(0-85+), 
# four large age groups: 0-19, 20-44, 45-64, 65+

agegrp1 = rep(1:8,each=8)  # as c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
agegrp2 = rep(1:8,times=8) # as c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
# Note symmetry, average contacts over large groups

# for age 0-19, 15 contacts within, 7 with age 20-39, 8 with age 40-64, and 3 with age 65+
# for age 20-44, 15 contacts within, 7 with age 0-19, 7 with age 40-64, and 5 with age 65+
# for age 45-64, 10 contacts within, 8 with age 0-19, 7 with age 20-44, and 7 with age 65+
# for age 65+, 8 contacts within, 3 with age 0-19, 5 with age 20-39, and 7 with age 40-64 
nconts = c(15,15, 7, 7, 8, 8, 3, 3,
           15,15, 7, 7, 8, 8, 3, 3,
           7, 7,15,15, 7, 7, 5, 5,
           7, 7,15,15, 7, 7, 5, 5,
           8, 8, 7, 7,10,10, 7, 7,
           8, 8, 7, 7,10,10, 7, 7,
           3, 3, 5, 5, 7, 7, 8, 8,
           3, 3, 5, 5, 7, 7, 8, 8)

contmat = data.table(cbind(agegrp1,agegrp2,nconts))

# reweight, smaller connections within each small age groups, round up to get more connections
contmat[,nconts:=round(nconts/2)]

setkey(contmat,agegrp1,agegrp2)

# tuning parameters;

# age specific infectivity and susceptibility, truncated exponential distribution
# assume independent infectivity and susceptibility, but they are likely correlated
# assume median 15% attack rate among young, and 20% attack rate among old people
# over the 7 day infectious period in household
nage4 = nrow(pop[agegrp==4,])
pop[,`:=`(infect = 1,
          suscept = ifelse(agegrp==4,rexp(popsize,1/0.3),rexp(popsize,1/0.2)))
    ][,`:=`(infect =ifelse(infect>1,0.99,infect),
        suscept =ifelse(suscept>1,0.99,suscept))]

options(scipen=999)
# check distributions, Q25%, median and Q75%
summary(pop[agegrp==2,]$suscept)
quantile(pop[agegrp==2,]$suscept,c(0,0.01,0.1,0.2,0.25,0.50,0.75,0.80,0.90,0.95,0.99,1))
summary(pop[agegrp==4,]$suscept)
quantile(pop[agegrp==4,]$suscept,c(0,0.01,0.1,0.2,0.25,0.50,0.75,0.80,0.90,0.95,0.99,1))

# rescale the susceptibility to per day for seven days
pop[,suscept:=suscept/7]
summary(pop[agegrp==2,]$suscept)
quantile(pop[agegrp==2,]$suscept,c(0,0.01,0.1,0.2,0.25,0.50,0.75,0.80,0.90,0.95,0.99,1))
summary(pop[agegrp==4,]$suscept)
quantile(pop[agegrp==4,]$suscept,c(0,0.01,0.1,0.2,0.25,0.50,0.75,0.80,0.90,0.95,0.99,1))


##################################################
######## define simulation function ##############
##################################################

# simple control measures
simepinet = function(
    lockdown = 0,  # default: 0, the strictest measure during the first main period, equivalent to double the effort
    netdense = 0,  # dense: median degree 40, moderate: median 20, sparse: median 10;
    control = 1,   # default: 1, reducing general contact control (e.g., voluntarily stay at home), time-varying probability of deleting networks
    distancing = 1, # default: 1, narrow sense of social distancing (but still have contacts with others) reduce infection risk by 20% 
    maskratio = 0.3, # default: 0.3, percent people wearing mask, also mask reducing probability of transmission by half or more
    masklong = 0,  # long term intervention after first wave, mostly mask
    agelong = 0, # extended control for aged people
    youngonly = 0, # only young people in the follow up 
    ageonly = 0,   # default: 0, limited prevention only among aged people, 0 (everybody),1 (age only within 150 days) 
    asymratio = 0.4, # default: 0.4, percent of asymptomatic cases (20-60%), assuming half of infectivity
    qtratio = 0.04,  # default: 0.04, quarantine ratio: 4% daily, 40% total over 7 days of infecting period,lockdown, double the ratio
    testing = 1,  # default: 1, testing frequency, affect the rate of case isolation and quarantine
    contweight = 0 # default: 0, using contact weight or not: poisson level of weight, that is, some contacts counts two or more times
  ) {
  
  gc()  #forced, clean up system

  ###########################################
  # run through network
  ##########################################
  
  edgenet = copy(edgenet0) # network edges
  
  # network density: 0: dense; 0.5: moderate, 0.75: sparse
  if (netdense !=0) {
    edgenet = edgenet[,pid1:=ifelse(rbinom(nrow(edgenet),1,netdense)==1,0,pid1)][pid1!=0,]  
    # ensure uniqueness of network connections
    edgenet = unique(edgenet,by=c("pid1","pid2"))
  }
  
  #rename edgenet with the primary pid
  setnames(edgenet,"pid1","pid")
  setkey(edgenet,pid)

  # create the compartments
  
  # working bucket for infectious, fluid, in and out 
  infecting = data.table(matrix(0,popsize+100,6))
  setnames(infecting,c("V1","V2","V3","V4","V5","V6"),c("pid","agegrp","timein","timeout","infectedby","asym"))
  
  # status bucket
  suscept = copy(pop) #work on susceptible, not to damage pop file
  setkey(suscept,pid)
  # add status and death indicators
  suscept[,c("status","death"):=list(0,0)]
  
  quarted = data.table(matrix(0,popsize+100,4))  # quarantined
  setnames(quarted,c("V1","V2","V3","V4"),c("pid","agegrp","timein","timeout"))
  exposed = data.table(matrix(0,popsize+100,6)) # exposed
  setnames(exposed,c("V1","V2","V3","V4","V5","V6"),c("pid","agegrp","timein","timeout","infectedby","asym"))
  infected = data.table(matrix(0,popsize+100,4)) #infected 
  setnames(infected,c("V1","V2","V3","V4"),c("pid","time","agegrp","asym"))
  removed = data.table(matrix(0,popsize+100,2))  #removed/recovered
  setnames(removed,c("V1","V2"),c("pid","time"))
  hospitalized = data.table(matrix(0,popsize+100,5)) # hospitalized
  setnames(hospitalized,c("V1","V2","V3","V4","V5"),c("pid","agegrp","time","timeout","death"))
  death = data.table(matrix(0,popsize+100,3))  # death from infection
  setnames(death,c("V1","V2","V3"),c("pid","agegrp","time"))
  
  # randomly initialize a few patients
  ninit = 20
  # start with young people
  newcase = sample(suscept[agegrp!=4,pid],ninit)
  newcase = suscept[pid %in% newcase,c("pid","agegrp")]
  
  infected[1:ninit,c("pid","agegrp","time","asym"):=list(newcase$pid,newcase$agegrp,rep(0,ninit),rep(0,ninit))]
  infecting[1:ninit,c("pid","agegrp","timein","timeout"):=list(newcase$pid,newcase$agegrp,rep(0,ninit),0 +
                                                                 round(rgamma(ninit,shape=2.8,scale=1.8)))]
  # prune those secondary new case connections (in-connect)
  edgenet[,pid2:=fifelse(pid2 %in% newcase$pid,0,pid2)]
  # flag susceptibles to be infected
  suscept = suscept[, status:=fifelse(pid %in% newcase$pid,1,0)]
  # general mortality rate
  suscept[,dprob2:=ifelse(agegrp==4,0.00014,0.00005)]
  
  # lapse = 1

  # run the lapse till fixed rounds
  for (lapse in 1:300){
    
    # susceptibles may die of other diseases, at a very low rate, depending on age
    nns = nrow(suscept)
    suscept[,death:=fifelse(death==0 & status==0 & pid!=0,rbinom(nns,1,dprob2),death)]
    dsuscept = suscept[death==1,pid]
    # remove those died
    suscept[,pid:=fifelse(death==1,0,pid)]
    # clean up all related connections if died
    edgenet[,c("pid","pid2"):=list(fifelse(pid %in% dsuscept,0,pid),fifelse(pid2 %in% dsuscept,0,pid2))]
    
    #current infecting case, 99999999 are time out cases
    infcase = infecting[pid!=0 & pid !=999999999,pid] 
    # exclude quarantined cases, they are not infecting others
    if (qtratio >0) infcase = infcase[!(infcase %in% quarted$pid)]
    
    # get the connections, but cannot connecting within cases
    worknet = unique(edgenet[pid %in% infcase,],by="pid2")
    worknet = worknet[pid2!=0,]
    
    # exclude those quarantined non-diseased people
    if (qtratio > 0) worknet = worknet[!(pid2 %in% quarted$pid),]
    
    # exclude those already infected cases
    worknet = worknet[!(pid2 %in% infected$pid),]

    # breaking connections randomly depending on timing, assessing controlling effects
    pnetbrk = 0  # probability of breaking network connections
    # time dependent, max 80% reduction during 60-90 days
    # lockdown always short, within 120 days
    if (control==1 | lockdown ==1){
      if (lapse < 15){
        pnetbrk = 0  # allow 15 days of spreading freely 
      } else if (lapse < 30){
        pnetbrk = ifelse(lockdown == 1, 0.4, 0.1)
      } else if (lapse < 60) {
        pnetbrk = ifelse(lockdown == 1, 0.8, 0.5)
      } else if (lapse < 90) {
        pnetbrk = ifelse(lockdown == 1, 0.4, 0.3)
      } else if (lapse < 120) {
        pnetbrk = ifelse(lockdown == 1, 0.2, 0.1)
      } else if (lapse < 150){
        pnetbrk = ifelse(lockdown == 1, 0.1, 0.1)
      } else pnetbrk = 0  # after 150 days, return to normal
    }
    
    nrwk = nrow(worknet)
    if (pnetbrk > 0) {
      if (ageonly == 0){
        worknet = worknet[,pid2:=ifelse(rbinom(nrwk,1,pnetbrk)==1,0,pid2)][pid2!=0,]  #everybody
      } else {
        # only among aged people
        worknet = worknet[,pid2:=ifelse(rbinom(nrwk,1,pnetbrk)==1 & age2 == 4,0,pid2)][pid2!=0,]
      } 
    }
    
    nrwk = nrow(worknet)
    # extended protection among elderly
    if (agelong == 1 & lapse >=120) {
      pnetbrk = 0.3  # old people still partially observe social distancing after reopen the society
      worknet = worknet[,pid2:=ifelse(rbinom(nrwk,1,pnetbrk)==1 & age2 == 4,0,pid2)][pid2!=0,]
    }
    
    setkey(worknet,pid)
    # recount length
    nrwk = nrow(worknet)
    
    # merge worknet with susceptible pop to get infection parameters
    tnewcase = worknet[,pid2]
    susprob = suscept[pid %in% tnewcase,c("pid","suscept")]
    setnames(susprob,"pid","pid2")
    
    # existing case parameters for infectivity
    infprob = suscept[pid %in% worknet$pid,c("pid","infect")]
    worknet = merge(worknet,infprob,on="pid",all.x=TRUE)
    
    # if infected by asymptomatic, infectivity reduce to half
    worknet = merge(worknet,infecting[,c("pid","asym")],on="pid",all.x=TRUE)
    worknet[,infect :=ifelse(asym==1,infect*0.5,infect)]
    
    # combined old case and new case, calculate probability of infection
    worknet = susprob[worknet,on="pid2"][,prb:=infect*suscept]
    
    # social distancing, lockdown
    # no social distancing measure in the first 15 days
    if ((lapse > 15) & (distancing == 1 | lockdown == 1)) {
      # later stage with lower social distancing
      if (ageonly ==0) {
        worknet[,prb:=ifelse(lapse<120, prb*0.8, prb*0.9)]
      } else {
      worknet[,prb:=ifelse(age2==4 & lapse<120,prb*0.8, prb*0.9)]
      }
    }
    
    if (agelong==1 & lapse >=120){
      # age only but with longer intervention, maintain social distancing 
      worknet[,prb:=ifelse(age2==4, prb*0.8, prb*0.9)]
    }
    
    # handling masks
    worknet[,mask:=0]
    # recount length
    nnwork = nrow(worknet)

    # for lockdown, always need some ratio of masks
    # nothing for the first 15 days, and smaller percent of masking during early period
    if (maskratio !=0){
      maskratio2 = maskratio
      if (lapse < 15) maskratio2 = 0 
      # everybody, before 60 days low mask ratio, and after 120 days, half of mask ratio if not masklong 
      else if (lapse<60) maskratio2=0.2
      else if (lapse<120) maskratio2=maskratio/2
      else if (lapse>=120){
        if (masklong == 0) maskratio2 = maskratio/2
        if (agelong == 1 | masklong == 1) maskratio2 = maskratio
      }
    
      # only among aged
      if (ageonly==1) {
        worknet[,mask := ifelse(age2==4,rbinom(nnwork,1,maskratio2),0)]
      } else {
        worknet[,mask := rbinom(nnwork,1,maskratio2)]
      }
    }  

    # adjust transmission probability for mask
    worknet[,prb:= ifelse(mask==1, prb*0.7, prb)]  # with mask, half the transmission risk
    
    # contact weight, repeated contacts 
    ssweight = function (cweight,prb){
        ss = sum(rbinom(cweight,1,prb))
        return(ifelse(ss>0,1,0))
      }
    
    if (contweight ==1){
      worknet = worknet[,prb:=ifelse(prb>0,prb,0)][,seqnum:=1:.N][,status:=ssweight(cweight,prb),by=seqnum][status==1,]
    } else {
      worknet = worknet[,prb:=ifelse(prb>0,prb,0)][,status:=rbinom(nnwork,1,prb)][status==1,]
    }
  
    # collecting new cases and old cases
    # recount length
    nnwork = nrow(worknet)
    newcase = worknet[,pid2]
    oldcase = worknet[,pid]
    nncase = length(newcase)
    newage = worknet[,age2]
    asymstatus = rbinom(nncase,1,asymratio)  #destined to be asymptomatic
    
    ########## handling exposed bucket 
    if (nncase > 0){
      # latent time in and time out of exposed, gamma distribution, 2-4 days
      timeout2 = worknet[,timeout := ifelse(age2==4,lapse + round(rgamma(nnwork,shape=2.3,scale=1.3)),
                                     lapse + round(rgamma(nnwork,shape=2.5,scale=1.6)))][,timeout]
      
      # add to exposed bucket
      nexp2 = nrow(exposed[pid!=0,])
      exposed[(nexp2+1):(nexp2+nncase),
                c("pid","agegrp","timein","timeout","infectedby","asym"):=
                  list(newcase,newage,rep(lapse,nncase),timeout2,oldcase,asymstatus)]
      
      # pruning the network, exclude those exposed cases as secondary connections, inconnection;
      edgenet[,pid2 := fifelse(pid2 %in% newcase,0,pid2)]
      
      # pruning susceptibles, exclude those exposed
      suscept[,status := fifelse(pid %in% newcase,0,1)]
    }
    
    ########Infecting
    #moving matured exposed into infecting and infected
    mtexp = exposed[timeout<=lapse & timeout != 0 & pid != 999999999,c("pid","agegrp","infectedby","asym")]
    nmcase = nrow(mtexp)
    
    # getting exposed into infecting and infected buckets
    if (nmcase > 0){
      # if new cases are old people, then longer timeout, gamma distribution of infectious period: 5 days
      timeout2 = mtexp[,timeout := ifelse(agegrp==4,lapse + round(rgamma(nmcase,shape=3.9,scale=1.8)),
                                       lapse + round(rgamma(nmcase,shape=2.8,scale=1.8)))][,timeout]
  
      # add to infected status bucket
      ninf = nrow(infected[pid!=0,])
      infected[(ninf+1):(ninf+nmcase),c("pid","agegrp","time","asym"):=list(mtexp$pid,mtexp$agegrp,rep(lapse,nmcase),mtexp$asym)]
  
      # add to infecting bucket
      ninf2 = nrow(infecting[pid!=0,])
      infecting[(ninf2+1):(ninf2+nmcase),
                c("pid","agegrp","timein","timeout","infectedby","asym"):=
                  list(mtexp$pid,mtexp$agegrp,rep(lapse,nmcase),timeout2,mtexp$infectedby,mtexp$asym)]
      
      # excluding those time out of exposed
      exposed[,pid :=ifelse(timeout>lapse | timeout == 0,pid,999999999)]
    }
    
    ######clear out infecting people:
    
    # part of infecting will be quarantined, thus not infecting in the future
    # default assume daily rate of 4% quarantined, thus overall 30% quarantined during infectious period
    # if lockdown and frequent testing, quarantine rate is doubled
#    if (testing==1) {
#      qtratio = qtratio*2
#    }
    # overall qt ratio cannot be more than 1
#    if (qtratio>0.1){
#      qtratio = 0.1
#    }
      
    if (qtratio >0){
      qtratio2 = qtratio
      nqlen = nrow(quarted[pid!=0,])
      ninf = nrow(infecting)
      if (lapse < 15) qtratio2 = 0.01
      else if (lapse<120) qtratio2 = qtratio*0.5
      if (lapse > 120 & masklong == 0) qtratio2 = qtratio*0.8 #after 120 days, 80% of quarantine speed if not longterm
      
      qt = infecting[,qtran :=rnbinom(ninf,1,qtratio2)][qtran==1,c("pid","timeout")]
      
      if (nrow(qt)>0) {
        # also 30% of their connected people
        qt2 = worknet[pid %in% qt,][,qtran:=rnbinom(ninf,1,0.3)][qtran==1,pid2]
        # for diseased people, they are not going to get out until recovered
        quarted=quarted[(nqlen+1):(nqlen+nrow(qt)),c("pid","timein","timeout"):=list(qt$pid,lapse,qt$timeout)]
        if (length(qt2)>0){
          # for nondiseased, 14 quarantine days
          quarted=quarted[(nqlen+1+nrow(qt)):(nqlen+nrow(qt)+length(qt2)),c("pid","timein","timeout"):=list(qt2,lapse,lapse+14)]
        }
      }
      # clean up quarantined
      quarted = quarted[,pid:=ifelse(timeout>lapse,pid,0)]
    }
    
    #people have a fixed probability per lapse of being hospitalized, higher among elderly
    thosp = copy(infecting[timeout>lapse & pid!=999999999 & asym==0,]) 
    if (nrow(thosp)> 0){
      thosp = thosp[,hprob:=ifelse(agegrp==4,0.06,ifelse(agegrp==3,0.04,0.025))][,hosp:=rbinom(nrow(thosp),1,hprob)][hosp==1,][,c("pid","agegrp")]
      nhosp = nrow(hospitalized[pid!=0,])
      if (length(thosp$pid)>0){
        nthosp = nrow(thosp)
        timeout2 = thosp[,timeout:=ifelse(agegrp==4,lapse + round(rgamma(nthosp,shape=5,scale=2.6)),
                                            lapse + round(rgamma(nthosp,shape=3,scale=2.3)))][,timeout]
        hospitalized[(nhosp+1):(nhosp+nthosp),c("pid","agegrp","time","timeout") :=list(thosp$pid,thosp$agegrp,rep(lapse,nthosp),timeout2)]
        # delete those hospitalized from infecting, assuming hospitalized are quarantined and not infectious any more
        infecting[,pid :=ifelse(pid %in% thosp,999999999,pid)]
      }
    }
    
    #some fixed rates of deaths from still hospitalized, higher among elderly
    tdeath = copy(hospitalized[pid!=999999999 & lapse<=timeout & death!=1,])
    if (nrow(tdeath)> 0){
      tdeath = tdeath[,dprob:=ifelse(agegrp==4,0.012,ifelse(agegrp==3,0.007,0.0035))][,death:=rbinom(nrow(tdeath),1,dprob)][death==1,][,c("pid","agegrp")]
      ndeath = nrow(death[pid!=0,])
      if (length(tdeath$pid)>0){
        ntdeath = nrow(tdeath)
        death[(ndeath+1):(ndeath+ntdeath),c("pid","agegrp","time") :=list(tdeath$pid, tdeath$agegrp,rep(lapse,ntdeath))]
        # mark those death from hospitalized as timed out 
        hospitalized[,c("death","timeout") :=list(ifelse(pid %in% tdeath,1,0),ifelse(pid %in% tdeath,lapse,timeout))]
      }
    }
    
    #for those passed timeout, then send them removed or recovery
    #we don't tally removed, as we miss those coming from hospitalized
    nrec = nrow(removed[pid!=0,])
    rec1 = copy(infecting[timeout<=lapse & timeout!=0 & pid!=999999999,][,pid])
    if (length(rec1)>0){
      removed[(nrec+1):(nrec+length(rec1)),c("pid","time") :=list(rec1,rep(lapse,length(rec1)))]
      # delete those recovered from infecting
      infecting[,pid :=ifelse(timeout>lapse | timeout == 0,pid,999999999)]
    }
  }
  
  return(list(infected=infected, hospitalized=hospitalized,death=death, suscept=suscept))
} # end of simulation

# save output
saveout = function (scene){
  
  infected = data.table(outnet$infected)
  infected = infected[pid!=0,]
  dim(infected)
  table(infected$time)
  
  saveRDS(infected,paste0(rootpath,"sim_infected_",scene,"50.RDS"))
  
  hospitalized = data.table(outnet$hospitalized)
  hospitalized = hospitalized[pid!=0,]
  dim(hospitalized)
  table(hospitalized$time)
  
  saveRDS(hospitalized,paste0(rootpath,"sim_hosp_",scene,"50.RDS"))
  
  death = data.table(outnet$death)
  death = death[pid!=0,]
  dim(death)
  table(death$time)
  saveRDS(death,paste0(rootpath,"sim_death_",scene,"50.RDS"))
}

###############FINISH RUN #####################


##########################################################
###########SIMULATING FINAL RESULTS#######################
##########################################################

# run the function under various scenarios;
# nothing;
outnet = simepinet (
  control = 0,
  netdense = 0,
  distancing = 0,   
  maskratio = 0,  
  masklong = 0,   
  ageonly = 0 ,     
  lockdown = 0,
  testing =0,
  asymratio = 0.4,  
  qtratio = 0,   
  contweight = 0) 

saveout("nothing")

# moderate intervention

outnet = simepinet (
  control = 1,
  netdense = 0,
  distancing = 1,   
  maskratio = 0.5,  
  masklong = 1,   
  ageonly = 0,     
  lockdown = 0,
  testing =0,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("moderate")


#network density: medium;
outnet = simepinet (
  control = 1,
  netdense = 0.5,
  distancing = 1,   
  maskratio = 0.5,  
  masklong = 1,   
  ageonly = 0,     
  lockdown = 0,
  testing =0,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("mod_medium_dense")

# density low;
outnet = simepinet (
  control = 1,
  netdense = 0.75,
  distancing = 1,   
  maskratio = 0.5,  
  masklong = 1,   
  ageonly = 0,     
  lockdown = 0,
  testing =0,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("mod_low_dense")

# lockdown with minimal follow up
outnet = simepinet (
  control = 1,   
  netdense = 0,
  distancing = 1,   
  maskratio = 0.5,  
  masklong = 0,   
  ageonly = 0,
  agelong = 0,     
  lockdown = 1, 
  testing = 0,
  asymratio = 0.4,  
  qtratio = 0.04,   
  contweight = 0) 
  
saveout("lockdown")

# open with high masks and distancing
outnet = simepinet (
  control = 1,   
  netdense = 0,
  distancing = 1,   
  maskratio = 0.8,  
  masklong = 1,   
  ageonly = 0 ,     
  agelong = 0 ,     
  lockdown = 1,
  testing = 1,
  asymratio = 0.4,  
  qtratio = 0.04,   
  contweight = 0) 

saveout("mask_long")


# open with distancing and high quarantine
outnet = simepinet (
  control = 1,   
  netdense = 0,
  distancing = 1,   
  maskratio = 0.5,  
  masklong = 1,   
  ageonly = 0,
  agelong = 0,
  lockdown = 1,
  testing = 1,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("mask_long_testing")

# open with more masks and distancing and high quarantine
outnet = simepinet (
  control = 1,   
  netdense = 0,
  distancing = 1,   
  maskratio = 0.8,  
  masklong = 1,   
  ageonly = 0,
  agelong = 0,
  lockdown = 1,
  testing = 1,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("mask2_long_testing")

# open with more masks and distancing and extended age protection
outnet = simepinet (
  control = 1,   
  netdense = 0,
  distancing = 1,   
  maskratio = 0.8,  
  masklong = 1,   
  ageonly = 0,
  agelong = 1,
  lockdown = 1,
  testing = 1,
  asymratio = 0.4,  
  qtratio = 0.08,   
  contweight = 0) 

saveout("mask_long_age")
 



##################################################################
##############DATA ANALYSIS: PARSING THE RESULTS #################
##################################################################

##### function for gam model #######
# good for gender, age, and gender-age specific models
####################################
myfitgam <- function(dailycase,time){
  Mixmodel <- gam(dailycase~s(time,bs="bs",k=12,m=2),family=nb,method="REML")
  
#  print(summary(Mixmodel))
  
  newd <-data.frame(time=time)
  
  # prediction + SE from NB model 
  fitTN <- predict( Mixmodel ,newd, type="response",se = TRUE )
  fitTN2 <- data.table(cbind(fitTN,newd))
  return(fitTN2)
}

# functions for checking results
checkresult = function(scene) {
  
  rm(list=c("infected","infected2","infected4","infected5"))
  
  # infected cases
  infected = data.table(readRDS(paste0(rootpath,"sim_infected_",scene,"50.RDS")))
  # drop age group 
  infected[,agegrp:=NULL]
  # merge personal information back to infected
  infected2 = pop[infected,on="pid"]
  setkey(infected2,"time")
  
  # average infectivity and susceptibility by time
  infected4 = infected2[,c("aveinfect","avesuscept"):=list(mean(infect),mean(suscept)),by=time][,.SD[.N],by=time][time!=0,]
  
  # average parameters
  print(ggplot(data=infected4[avesuscept<0.95],aes(x=time,y=avesuscept))+geom_point()+geom_smooth())

  # count by age groups
  infected5 = dcast(infected2[,c("time","agegrp")],time~agegrp)
  setnames(infected5,c("1","2","3","4"),c("age10","age30","age50","age70"))
  infected5[,c("agetotal","age4pct"):=list(age10+age30+age50+age70,age70/(age10+age30+age50+age70)),by=time]
  
  print(head(infected5))

  print(ggplot(data=infected5,aes(x=time,y=age4pct))+geom_point()+geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")))
  
  fittotal=myfitgam(infected5$agetotal,infected5$time)
  fitage10=myfitgam(infected5$age10,infected5$time)
  fitage30=myfitgam(infected5$age30,infected5$time)
  fitage50=myfitgam(infected5$age50,infected5$time)
  fitage70=myfitgam(infected5$age70,infected5$time)
  
  gfitdata <- data.frame(newdate=infected5$time,
                         fitTNA  = fittotal$fit,
                         fitTN10 = fitage10$fit,
                         fitTN30 = fitage30$fit,
                         fitTN50 = fitage50$fit,
                         fitTN70 = fitage70$fit,
                         Dailycase10=infected5$age10,
                         Dailycase30=infected5$age30,
                         Dailycase50=infected5$age50,
                         Dailycase70=infected5$age70,
                         Dailycase=infected5$agetotal)
  
  gp1=ggplot(gfitdata) + 
    geom_line(aes(newdate, fitTN10,color="0-19"),size = 1,linetype="longdash") +
    geom_line(aes(newdate, fitTN30,color="20-44"),size = 1,linetype="dotdash") +
    geom_line(aes(newdate, fitTN50,color="45-64"),size = 1,linetype="dashed") +
    geom_line(aes(newdate, fitTN70,color="65+"),size = 1,linetype="solid") +
    scale_color_manual("Age",values=c("0-19"="blue","20-44"="purple","45-64"="orange","65+"="red")) + theme_bw() +
    labs(x = "Time",y="Daily new cases")+ scale_x_continuous(breaks=seq(0,300,30)) + ggtitle(paste(scene,"Epidemic Curves")) + 
    theme(legend.position=c(.8,.8))

  print(gp1)
  
  # hospitalization
  hospitalized = data.table(readRDS(paste0(rootpath,"sim_hosp_",scene,"50.RDS")))
  hospitalized[,agegrp:=NULL]
  
  hosped2 = pop[hospitalized,on="pid"]
  setkey(hosped2,"time")
  
  # count by age agroups
  hosped5 = dcast(hosped2[,c("time","agegrp")],time~agegrp)
  setnames(hosped5,c("1","2","3","4"),c("age10","age30","age50","age70"))
  hosped5[,c("agetotal","age4pct"):=list(age10+age30+age50+age70,age70/(age10+age30+age50+age70)),by=time]
  
  print(ggplot(data=hosped5,aes(x=time,y=age4pct))+geom_point()+geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")))
  
  hfittotal=myfitgam(hosped5$agetotal,hosped5$time)
  hfitage10=myfitgam(hosped5$age10,hosped5$time)
  hfitage30=myfitgam(hosped5$age30,hosped5$time)
  hfitage50=myfitgam(hosped5$age50,hosped5$time)
  hfitage70=myfitgam(hosped5$age70,hosped5$time)
  
  hfitdata <- data.frame(newdate=hosped5$time,
                         fitTNA  = hfittotal$fit,
                         fitTN10 = hfitage10$fit,
                         fitTN30 = hfitage30$fit,
                         fitTN50 = hfitage50$fit,
                         fitTN70 = hfitage70$fit,
                         Dailycase10=hosped5$age10,
                         Dailycase30=hosped5$age30,
                         Dailycase50=hosped5$age50,
                         Dailycase70=hosped5$age70,
                         Dailycase=hosped5$agetotal)
  
  gp2=ggplot(hfitdata) + 
    geom_line(aes(newdate, fitTN10,color="0-19"),size = 1,linetype="longdash") +
    geom_line(aes(newdate, fitTN30,color="20-44"),size = 1,linetype="dotdash") +
    geom_line(aes(newdate, fitTN50,color="45-64"),size = 1,linetype="dashed") +
    geom_line(aes(newdate, fitTN70,color="65+"),size = 1,linetype="solid") +
    scale_color_manual("Age",values=c("0-19"="blue","20-44"="purple","45-64"="orange","65+"="red")) + theme_bw() +
    labs(x = "Time",y="Daily new hospitalizations")+ scale_x_continuous(breaks=seq(0,300,30)) + ggtitle(paste(scene,"Hospitalizations")) + 
    theme(legend.position=c(.8,.8))
  print(gp2)
}

checkresult("nothing")
checkresult("moderate")
checkresult("mod_medium_dense")
checkresult("mod_low_dense")
checkresult("lockdown")
checkresult("mask_long")
checkresult("mask_long_testing")
checkresult("mask2_long_testing")
checkresult("mask_long_age")


###########PLOT VARIOUS FACTORS TOGETHER#########################

# reshape data
reshapedata = function(tdt){
  # drop age group 
  tdt[,agegrp:=NULL]
  
  # merge personal information back to infected
  tdt2 = pop[tdt,on="pid"]
  setkey(tdt2,"time")

  # count by age groups
  tdt3 = dcast(tdt2[,c("time","agegrp")],time~agegrp)
  setnames(tdt3,c("1","2","3","4"),c("age10","age30","age50","age70"))
  tdt3[,c("agetotal","age4pct"):=list(age10+age30+age50+age70,age70/(age10+age30+age50+age70)),by=time]
  return(tdt3)
} 


# combine data;
combinedata = function(dep){
  
  dt1 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_nothing50.RDS")))
  dt2 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_moderate50.RDS")))
  dt3 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mod_medium_dense50.RDS")))
  dt4 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mod_low_dense50.RDS")))

  dt1a = reshapedata(dt1)
  dt2a = reshapedata(dt2)
  dt3a = reshapedata(dt3)
  dt4a = reshapedata(dt4)

  # not all groups have the same length of time
  tdt1a = dt1a[,c("time","agetotal")]
  setnames(tdt1a,"agetotal","nothing")
  tdt2a = dt2a[,c("time","agetotal")]
  setnames(tdt2a,"agetotal","moderate")
  tdt3a = dt3a[,c("time","agetotal")]
  setnames(tdt3a,"agetotal","medium_dense")
  tdt4a = dt4a[,c("time","agetotal")]
  setnames(tdt4a,"agetotal","low_dense")

  setkey(tdt1a,"time")
  setkey(tdt2a,"time")
  setkey(tdt3a,"time")
  setkey(tdt4a,"time")

  # multiple merge/join by time
  totaldt = tdt1a[tdt2a,][tdt3a,][tdt4a,]
  
  # age 70
  # not all groups have the same length of time
  tdt1b = dt1a[,c("time","age70")]
  setnames(tdt1b,"age70","nothing")
  tdt2b = dt2a[,c("time","age70")]
  setnames(tdt2b,"age70","moderate")
  tdt3b = dt3a[,c("time","age70")]
  setnames(tdt3b,"age70","medium_dense")
  tdt4b = dt4a[,c("time","age70")]
  setnames(tdt4b,"age70","low_dense")

  setkey(tdt1b,"time")
  setkey(tdt2b,"time")
  setkey(tdt3b,"time")
  setkey(tdt4b,"time")
  
  # multiple merge/join by time
  age70dt = tdt1b[tdt2b,][tdt3b,][tdt4b,]

  return(list(totaldt=totaldt, age70dt=age70dt))
}

#combining results and plot

# infected cases
infected = combinedata("infected")

totalinfect = data.table(infected$totaldt)
age70infect = data.table(infected$age70dt)

print(head(totalinfect))
print(head(age70infect))

# obtain fitted values
fit1=myfitgam(totalinfect$nothing,totalinfect$time)
fit2=myfitgam(totalinfect$moderate,totalinfect$time)
fit3=myfitgam(totalinfect$medium_dense,totalinfect$time)
fit4=myfitgam(totalinfect$low_dense,totalinfect$time)

gfitdata <- data.frame(newdate=totalinfect$time,
                       fitTN1 = fit1$fit,
                       fitTN2 = fit2$fit,
                       fitTN3 = fit3$fit,
                       fitTN4 = fit4$fit,
                       Dailycase1=totalinfect$nothing,
                       Dailycase2=totalinfect$moderate,
                       Dailycase3=totalinfect$medium_dense,
                       Dailycase4=totalinfect$low_dense)

gp1=ggplot(gfitdata) + 
  geom_line(aes(newdate, fitTN1,color="No intervention"),size = 1,linetype="longdash") +
  geom_line(aes(newdate, fitTN2,color="Moderate intervention"),size = 1,linetype="dotdash") +
  geom_line(aes(newdate, fitTN3,color="Medium density"),size = 1,linetype="dashed") +
  geom_line(aes(newdate, fitTN4,color="Low density"),size = 1,linetype="solid") +
  scale_color_manual("Network",values=c("No intervention"="red","Moderate intervention"="blue","Medium density"="orange","Low density"="green")) + 
  theme_bw() + labs(x = "Time",y="Predicted daily new cases")+ scale_x_continuous(breaks=seq(0,300,15)) + ggtitle(paste("Network density","all age groups")) + 
  theme(legend.position=c(.8,.8))

print(gp1)
  
# obtain fitted values among aged people
fit1=myfitgam(age70infect$nothing,age70infect$time)
fit2=myfitgam(age70infect$moderate,age70infect$time)
fit3=myfitgam(age70infect$med_dense,age70infect$time)
fit4=myfitgam(age70infect$low_dense,age70infect$time)

agfitdata <- data.frame(newdate=age70infect$time,
                       fitTN1 = fit1$fit,
                       fitTN2 = fit2$fit,
                       fitTN3 = fit3$fit,
                       fitTN4 = fit4$fit,
                       Dailycase1=age70infect$nothing,
                       Dailycase2=age70infect$moderate,
                       Dailycase3=age70infect$medium_dense,
                       Dailycase4=age70infect$low_dense)

gp2=ggplot(agfitdata) + 
  geom_line(aes(newdate, fitTN1,color="No intervention"),size = 1,linetype="longdash") +
  geom_line(aes(newdate, fitTN2,color="Moderate intervention"),size = 1,linetype="dashed") +
  geom_line(aes(newdate, fitTN3,color="Medium density"),size = 1,linetype="dotted") +
  geom_line(aes(newdate, fitTN4,color="Low density"),size = 1,linetype="solid") +
  scale_color_manual("Network",values=c("No intervention"="red","Moderate intervention"="red","medium density"="orange","Low density"="green")) + 
  theme_bw() + labs(x = "Time",y="Predicted daily New Cases")+ scale_x_continuous(breaks=seq(0,300,15)) + ggtitle(paste("Infected","Age 60 or older")) + 
  theme(legend.position=c(.8,.8))

print(gp2)


# other interventions/factors;
# combine data;
combinedata = function(dep){
  
  dt1 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mask_long_age50.RDS")))
  dt2 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mask_long50.RDS")))
  dt3 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mask_long_testing50.RDS")))
  dt4 = data.table(readRDS(paste0(rootpath,"sim_",dep,"_mask2_long_testing50.RDS")))
  
  dt1a = reshapedata(dt1)
  dt2a = reshapedata(dt2)
  dt3a = reshapedata(dt3)
  dt4a = reshapedata(dt4)
  
  # not all groups have the same length of time
  tdt1a = dt1a[,c("time","agetotal")]
  setnames(tdt1a,"agetotal","mask_long_age")
  tdt2a = dt2a[,c("time","agetotal")]
  setnames(tdt2a,"agetotal","mask_long")
  tdt3a = dt3a[,c("time","agetotal")]
  setnames(tdt3a,"agetotal","mask_long_testing")
  tdt4a = dt4a[,c("time","agetotal")]
  setnames(tdt4a,"agetotal","mask2_long_testing")
  
  setkey(tdt1a,"time")
  setkey(tdt2a,"time")
  setkey(tdt3a,"time")
  setkey(tdt4a,"time")
  
  # multiple merge/join by time
  totaldt = tdt1a[tdt2a,][tdt3a,][tdt4a,]
  
  # age 70
  # not all groups have the same length of time
  tdt1b = dt1a[,c("time","age70")]
  setnames(tdt1b,"age70","mask_long_age")
  tdt2b = dt2a[,c("time","age70")]
  setnames(tdt2b,"age70","mask_long")
  tdt3b = dt3a[,c("time","age70")]
  setnames(tdt3b,"age70","mask_long_testing")
  tdt4b = dt4a[,c("time","age70")]
  setnames(tdt4b,"age70","mask2_long_testing")
  
  setkey(tdt1b,"time")
  setkey(tdt2b,"time")
  setkey(tdt3b,"time")
  setkey(tdt4b,"time")
  
  # multiple merge/join by time
  age70dt = tdt1b[tdt2b,][tdt3b,][tdt4b,]
  
  return(list(totaldt=totaldt, age70dt=age70dt))
}

#combining results and plot

# infected cases
infected = combinedata("infected")

totalinfect = data.table(infected$totaldt)
age70infect = data.table(infected$age70dt)

print(head(totalinfect))
print(head(age70infect))

# obtain fitted values
fit1=myfitgam(totalinfect$mask_long_age,totalinfect$time)
fit2=myfitgam(totalinfect$mask_long,totalinfect$time)
fit3=myfitgam(totalinfect$mask_long_testing,totalinfect$time)
fit4=myfitgam(totalinfect$mask2_long_testing,totalinfect$time)

gfitdata <- data.frame(newdate=totalinfect$time,
                       fitTN1 = fit1$fit,
                       fitTN2 = fit2$fit,
                       fitTN3 = fit3$fit,
                       fitTN4 = fit4$fit,
                       Dailycase1=totalinfect$mask_long_age,
                       Dailycase2=totalinfect$mask_long,
                       Dailycase3=totalinfect$mask_long_testing,
                       Dailycase4=totalinfect$mask2_long_testing)

gp1=ggplot(gfitdata) + 
  geom_line(aes(newdate, fitTN1,color="Targeted intervention after lockdown"),size = 1,linetype="longdash") +
  geom_line(aes(newdate, fitTN2,color="More mask after lockdown"),size = 1,linetype="dotdash") +
  geom_line(aes(newdate, fitTN3,color="More testing after lockdown"),size = 1,linetype="dashed") +
  geom_line(aes(newdate, fitTN4,color="More mask and testing after lockdown"),size = 1,linetype="solid") +
  scale_color_manual("Intervention",values=c("Targeted intervention after lockdown"="red","More mask after lockdown"="blue","More testing after lockdown"="orange","More mask and testing after lockdown"="green")) + 
  theme_bw() + labs(x = "Time",y="Predicted daily new cases")+ scale_x_continuous(breaks=seq(0,300,15)) + ggtitle(paste("Various interventions,","all age groups")) + 
  theme(legend.position=c(.8,.8))

print(gp1)

# obtain fitted values among aged people
fit1=myfitgam(age70infect$mask_long_age,age70infect$time)
fit2=myfitgam(age70infect$mask_long,age70infect$time)
fit3=myfitgam(age70infect$mask_long_testing,age70infect$time)
fit4=myfitgam(age70infect$mask2_long_testing,age70infect$time)

agfitdata <- data.frame(newdate=age70infect$time,
                        fitTN1 = fit1$fit,
                        fitTN2 = fit2$fit,
                        fitTN3 = fit3$fit,
                        fitTN4 = fit4$fit,
                        Dailycase1=age70infect$mask_long_age,
                        Dailycase2=age70infect$mask_long,
                        Dailycase3=age70infect$mask_long_testing,
                        Dailycase4=age70infect$mask2_long_testing)

gp2=ggplot(agfitdata) + 
  geom_line(aes(newdate, fitTN1,color="Targeted intervention after lockdown"),size = 1,linetype="longdash") +
  geom_line(aes(newdate, fitTN2,color="More mask after lockdown"),size = 1,linetype="dotdash") +
  geom_line(aes(newdate, fitTN3,color="More testing after lockdown"),size = 1,linetype="dashed") +
  geom_line(aes(newdate, fitTN4,color="More mask and testing after lockdown"),size = 1,linetype="solid") +
  scale_color_manual("Intervention",values=c("Targeted intervention after lockdown"="red","More mask after lockdown"="blue","More testing after lockdown"="orange","More mask and testing after lockdown"="green")) + 
  theme_bw() + labs(x = "Time",y="Predicted daily new cases")+ scale_x_continuous(breaks=seq(0,300,15)) + ggtitle(paste("Various intervention","elderly")) + 
  theme(legend.position=c(.8,.8))

print(gp2)





#############################
# time varying R (Cori A, AJE paper)
library(incidence)
library(EpiEstim)

# Cori A. instanteneous R, AJE 2013, Epidemic 2019
# simple parametric for serial interval distribution is enough;
# other methods tried, but also need assumptions ;
# possible: before and after stay home rule;

#genearating time with mean 7.5, sd 3.4, NEJM Li et al.
# mean 4.7, sd: 2.9, Nishiura et al 2020
# much higher R, maybe unrealiastic in some region;
#model assumes gamma distribution, always positive

# by different age group
epic<-as.vector(infected5$agetotal)

metro_r<- estimate_R(epic, 
                     method="parametric_si",
                     config = make_config(list(mean_si = 4.7, std_si = 2.9)))
summary(metro_r)
metro_r$dates<-infected5$time
plot(metro_r, what = c("incid"),options_I=list(col=c("blue")))
plot(metro_r, what = c("R"))



epic<-as.vector(infected5$age70)

metro_r<- estimate_R(epic, 
                     method="parametric_si",
                     config = make_config(list(mean_si = 4.7, std_si = 2.9)))
summary(metro_r)
metro_r$dates<-infected5$time
plot(metro_r, what = c("incid"),options_I=list(col=c("blue")))
plot(metro_r, what = c("R"))

# metro_r$R

