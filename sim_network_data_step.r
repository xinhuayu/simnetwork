####################################################################
# network simulation code to explore susceptibility (Masking effect); 
# survivor bias, incomplete control and reopening
# population heterogeneity, and long tail and second wave
# Author: Xinhua Yu
# Revised date: July 6, 2020
###################################################################

rm(list=ls())
gc()
options(scipen = 999)
options(digits = 3)

library(data.table)

############DATA STEP: BUILD NETWORKS, RUN ONCE ONLY ##############
# build a small world network with random within
# bottom up approach
###################################################################

# also try parallel computing
library(iterators)
library(foreach)
library(doParallel)

# initialize
cl<-makeCluster(16)  # check core numbers depending on your computer
registerDoParallel(cl)


#############################################
# outside contact matrix, eighteen small every five year age groups(0-85+), 
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

# skeleton function for loop convenience to create small data 
# flab for saved file label to avoid overwriting old files
floop = function(pop,contmat,ii,flab) {
  
  ############################################
  # build age stratified network connections, within age and cross age groups
  ############################################
  gc()
  
  setkey(pop,pid)
  
  ################### STEP 1: within age group #######################
  # function for building random network within each age group 
  withingrpnet = function(dage,contmat,inc){
    # add one (empty) index to avoid accidental overwriting;
    inc = inc + 1
    popage = pop[magegrp==dage,]
    ncontage = contmat[agegrp1==dage & agegrp2 == dage,nconts]
    npopage = nrow(popage)
    # calculate probability of connections between any people, based on contact matrix
    probcontact = (npopage*ncontage)/(choose(npopage,2))
    
    # collect connections;
    for(i in 1:npopage){
      for(j in (i+1):npopage){
        if (rbinom(1,1,probcontact) == 1) {
          # repeat the connections for both parties: ids and alters
          # so we only need to search the primary id (pid1) later
          tedgenet[inc,c("pid1","pid2","type","age1","age2"):=list(popage[i,pid],popage[j,pid],1,dage,dage)]
          tedgenet[inc+1,c("pid1","pid2","type","age1","age2"):=list(popage[j,pid],popage[i,pid],1,dage,dage)]
          inc = inc+2
        }
      }
    }
  }
  
  # loops to create networks within the age group
  # %dopar% with parallel, %do% only with serial
  # need .multicombine and .packages to send them to each R worker 
  edgenet1 = foreach( i=1:8,.combine=function(...) rbindlist(list(...)),.multicombine=TRUE,.packages="data.table") %dopar% {
    # initialize a huge network connection data
    tedgenet = data.table(matrix(0,20000000,3))
    setnames(tedgenet,c("V1","V2","V3"),c("pid1","pid2","type"))
    edglen =  1
    withingrpnet(i,contmat,edglen)
    tedgenet[pid1!=0,]
  }
  
  ###############STEP 2: between age group####################
  # function for building network across age groups
  btwgrpnet = function(cage,dage,contmat,inc){
    inc = inc + 1;
    popage1 = pop[magegrp==cage,]
    popage2 = pop[magegrp==dage,]
    npage1 =nrow(popage1)
    npage2 =nrow(popage2)
    ncontage = contmat[agegrp1==cage & agegrp2==dage,nconts]
    
    # calculate probability of connections, ensure sufficient connections for the smaller group
    probcontact = (min(npage1,npage2)*ncontage)/(npage1*npage2)
    
    for(i in 1:npage1){
      for(j in 1:npage2){
        if (rbinom(1,1,probcontact) == 1) {
          # two copies with flipped connections for ids and alters
          tedgenet[inc,c("pid1","pid2","type","age1","age2"):=list(popage1[i,pid],popage2[j,pid],2,cage,dage)]
          tedgenet[inc+1,c("pid1","pid2","type","age1","age2"):=list(popage2[j,pid],popage1[i,pid],2,dage,cage)]
          inc = inc+2
        }
      }
    }
  }
  
  # generate between group connections;
  # %dopar% with parallel, %do% only with serial
  edgenet2 = foreach( i=1:7,.combine=function(...) rbindlist(list(...)),.multicombine=TRUE,.packages="data.table") %dopar% {
    # initialize a huge network connection data
    tedgenet = data.table(matrix(0,50000000,3))
    setnames(tedgenet,c("V1","V2","V3"),c("pid1","pid2","type"))
    for (j in (i+1):8){
      edglen = nrow(tedgenet[pid1!=0,])+1
      btwgrpnet(i,j,contmat,edglen)
    }
    tedgenet[pid1!=0,]
  }
  
  ##### merge within and between age group networks, also delete zero connections
  edgenet = data.table::rbindlist(list(edgenet1,edgenet2))
  edgenet[pid1!=0,]
  
  nrow(edgenet)
  table(edgenet$type)
  
  ##### make sure change file name tag to avoid accidental overwriting
  saveRDS(edgenet,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,flab,".RDS"))
  saveRDS(pop,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_pop",ii,flab,".RDS"))
  
  # check contact distributions by age group
  table(edgenet$age1,edgenet$age2)
  table(pop$magegrp)  # need to divide by two in the within group

} # end of skeleton loop function


##############STEP 0: create small communities#######################
# first create small batch-run for small communities (every 5000);
for (ii in 1:100) {
  popsize = 5000
  pop = data.table(matrix(0,popsize,7))
  setnames(pop,c("V1","V2","V3","V4","V5","V6","V7"),c("pid","rannum","gender","ncontacts","agegrp","infect","suscept"))
  
  # randomly generate age and gender groups, initialize data
  # https://www.census.gov/data/tables/2019/demo/age-and-sex/2019-age-sex-composition.html
  # eighteen age groups, but will group them into four at the end;
  pop[,`:=`(
    pid = seq_along(1:nrow(pop)),
    rannum = runif(nrow(pop)),
    gender = rbinom(nrow(pop),1,0.5)+1,
    ncontacts = 0)][,magegrp:=fifelse(rannum < 0.123,1,
                              fifelse(rannum < 0.251,2,
                              fifelse(rannum < 0.389,3,
                              fifelse(rannum < 0.583,4,
                              fifelse(rannum < 0.709,5,
                              fifelse(rannum < 0.837,6,
                              fifelse(rannum < 0.935,7,8)))))))][,rannum:=NULL]
  
  # combine into larger groups used in later simulations
  pop[,agegrp:=fifelse(magegrp %in% c(1,2),1,
               fifelse(magegrp %in% c(3,4),2,
               fifelse(magegrp %in% c(5,6),3,4)))]
  
  prop.table(table(pop$magegrp))
  prop.table(table(pop$agegrp))
  prop.table(table(pop$gender))
  
  popsize = nrow(pop)
  # batch run
  floop(pop,contmat,ii,"a")
}

# fix pid codes in the saved data;
for (ii in 1:100) {
  edgeneta = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,"a.RDS")))
  popa = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_pop",ii,"a.RDS")))
  popa[,pid := ii*1000000+pid][,community:=ii]
  edgeneta[,c("pid1","pid2") := list(ii*1000000+pid1,ii*1000000+pid2)]
  saveRDS(edgeneta,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_networka",ii,"a.RDS"))
  saveRDS(popa,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",ii,"a.RDS"))
}  

# check degreee distribution for each community
neta = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_networka1a.RDS")))
setkey(neta,"pid1")
dd = neta[,degree:=.N,by="pid1"][,.SD[.N],by="pid1"]
# degree distributions
hist(dd$degree,breaks=50)
summary(dd$degree)
quantile(dd$degree,c(0,0.01,0.25,0.50,0.75,0.90,0.95,0.99,1))

# create the final large population;
popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",1,"a.RDS")))
pop = copy(popb1)
for (ii in 2:100) {
  popb2 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",ii,"a.RDS")))
  pop = rbindlist(list(pop,popb2))
}

saveRDS(pop,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_pop50_all.RDS"))

###################################

##############STEP 3: STEPWEDGE-CASE: from small to medium and large community #####################
# a random sample to have long distance (cross-community) links
# like small world network, may create some people with a lot of links
#
# first link within every 25,000 big community
# these population will be those having larger number of connections (in and outside communities)
# random poisson distribution of degrees;

pop=data.table()

for (ii in 0:19){
  for (jj in 1:5){
    popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",jj+ii*5,"a.RDS")))
    if (jj==1) pop = popb1[runif(nrow(popb1),0,1)<0.2,]
    else  pop = data.table(rbindlist(list(pop,popb1[runif(nrow(popb1),0,1)<0.2,])))
  }
  popsize=nrow(pop)
  floop(pop,contmat,ii,"bb")
}  

# link within every 50,000 big community
pop=data.table()

for (ii in 0:9){
  for (jj in 1:10){
    popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",jj+ii*10,"a.RDS")))
    if (jj==1) pop = popb1[runif(nrow(popb1),0,1)<0.15,]
    else  pop = data.table(rbindlist(list(pop,popb1[runif(nrow(popb1),0,1)<0.15,])))
  }
  popsize=nrow(pop)
  
  floop(pop,contmat,ii,"cc")
}  

# link within every 100,000 big community
pop=data.table()

for (ii in 0:4){
  for (jj in 1:20){
    popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",jj+ii*20,"a.RDS")))
    if (jj==1) pop = popb1[runif(nrow(popb1),0,1)<0.1,]
    else  pop = data.table(rbindlist(list(pop,popb1[runif(nrow(popb1),0,1)<0.1,])))
  }
  popsize=nrow(pop)
  
  floop(pop,contmat,ii,"dd")
}  

# link within every 250,000 big community
pop=data.table()

for (ii in 0:1){
  for (jj in 1:40){
    popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",jj+ii*50,"a.RDS")))
    if (jj==1) pop = popb1[runif(nrow(popb1),0,1)<0.1,]
    else  pop = data.table(rbindlist(list(pop,popb1[runif(nrow(popb1),0,1)<0.1,])))
  }
  popsize=nrow(pop)
  
  floop(pop,contmat,ii,"ee")
}  

# one last random link within the big region, 500,000;
# give everybody some random chances of connecting with somebody unknown;
popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",1,"a.RDS")))
pop = popb1[runif(nrow(popb1),0,1)<0.05,]
for (ii in 2:100) {
  popb2 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",ii,"a.RDS")))
  pop = rbindlist(list(pop,popb2[runif(nrow(popb2),0,1)<0.05,]))
}

popsize=nrow(pop)

floop(pop,contmat,555,"ff")

# do twice for the large network
popb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",1,"a.RDS")))
pop = popb1[runif(nrow(popb1),0,1)<0.05,]
for (ii in 2:100) {
  popb2 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_popa",ii,"a.RDS")))
  pop = rbindlist(list(pop,popb2[runif(nrow(popb2),0,1)<0.05,]))
}

popsize=nrow(pop)

floop(pop,contmat,555,"gg")

###### pool all connections together;
# individual 5000 pop
netb1 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_networka",1,"a.RDS")))
neta = copy(netb1)
for (ii in 2:100){
  netb2 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_networka",ii,"a.RDS")))
  neta = rbindlist(list(neta,netb2))
}
#25,000 pop
for (ii in 0:19){
  netb3 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,"bb.RDS")))
  netb3[,type:=3]
  neta = rbindlist(list(neta,netb3))
}
# 50,000 pop
for (ii in 0:9){
  netb4 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,"cc.RDS")))
  netb4[,type:=4]
  neta = rbindlist(list(neta,netb4))
}
# 100,000 pop
for (ii in 0:4){
  netb5 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,"dd.RDS")))
  netb5[,type:=5]
  neta = rbindlist(list(neta,netb5))
}
# 250,000 pop
for (ii in 0:1){
  netb6 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",ii,"ee.RDS")))
  netb6[,type:=6]
  neta = rbindlist(list(neta,netb6))
}

# last two big connections
netb7 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",555,"ff.RDS")))
netb7[,type:=7]
neta = rbindlist(list(neta,netb7))

netb8 = data.table(readRDS(paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network",555,"gg.RDS")))
netb8[,type:=7]
neta = rbindlist(list(neta,netb8))

neta= unique(neta,by=c("pid1","pid2"))

saveRDS(neta,paste0("C:/Users/xinhuayu/Desktop/COVID19/sim_network_net50_all.RDS"))

# degree distribution
neta = data.table(readRDS("C:/Users/xinhuayu/Desktop/COVID19/sim_network_net50_all.RDS"))
setkey(neta,"pid1")
dd = neta[,degree:=.N,by="pid1"][,.SD[.N],by="pid1"]

summary(dd$degree)
quantile(dd$degree,c(0,0.01,0.25,0.50,0.75,0.90,0.95,0.99,1))


#################################################################################
stopCluster(cl)

