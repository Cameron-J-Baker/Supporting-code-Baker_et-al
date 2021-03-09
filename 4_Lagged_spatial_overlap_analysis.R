# Author: Cameron J baker
# 03-12-2020
# Email: cameron.baker@uqconnect.edu.au
# Code for the analyses contained within the manuscript "Membership of crocodile social environments is dictated by male 
# philopatry" that is submitted for review to Behavioral Ecology
# Aim : To examine the long term stability of HR overlaps between conspecifics acoss years using the Lagged spatial 
# overlap rate

# Load the required R packages
library(asnipe)
library(data.table)
library(tidyverse)
library(lubridate)
library(plyr)
library(foreach)
library(doParallel)
library(rlist)
library(tibble)
library(ggplot2)
library(ggpubr)

##########################################################################################################################
# Import the required datasets ###########################################################################################

# The movement tactics of each crocodile
croc_categories <- fread("Data/Crocodile_movement_tactics.csv")

#====================================================================================================================#
# Impoty the monthly overlap matrices
#====================================================================================================================#
files <- list.files(path ="Output/Overlap_Matrices/", full.names = T) # create a vector list with the names of all files

files <- files[13:113] # exclude all data before the 2010-08 due to the ammount of missing months present for this period

monthly <- lapply(files, fread) # Import all files using lapply

# Convert each of the files imported into matrices and move the TRANSMITTERIDs into the rownames
for(i in 1:length(monthly)){
  monthly[[i]] <- as.data.frame(monthly[[i]]) # convert to data.frame
  colnames(monthly[[i]]) = monthly[[i]][1, ] # convert the first row back to column names
  monthly[[i]] = monthly[[i]][-1, ] # remove the first row with the names  
  rownames(monthly[[i]]) <- c() # remove the current row names
  monthly[[i]] <- tibble::column_to_rownames(monthly[[i]], var = "NA") # convert iD back into row.names
  # Remove the individuals that are not included within the study
  monthly[[i]] <- monthly[[i]][rownames(monthly[[i]]) %in% ID, ]
  monthly[[i]] <- monthly[[i]][ , colnames(monthly[[i]]) %in% ID]
  
  monthly[[i]] <- as.matrix(monthly[[i]])}

# To ensure that things are kept in the order that they were imported,create the month list from the names of each file
timepoints <- substr(files, 25, 31)

########################################################################################################################
# Import the required functions ########################################################################################

source("R_functions/LARprep.R")

########################################################################################################################
# Run the LSOR analysis on the observed datasets

# Resident-resident LSOR
prepdataR <- LARprep(monthly, years = timepoints, class = c("Resident")) 
groupsdataR <- prepdataR[[1]] # extract the groups
time_pointsR <- prepdataR[[2]] # extract the time series of the data
lagged_ratesRprep <- LAR(groupsdataR,time_pointsR,31, identities = categories$TRANSMITTERID)
lagged_ratesRprep <- as.data.frame(lagged_ratesRprep)
lagged_ratesRprep$Class <- "Resident-Resident"

# Nomadic-Nomadic LSOR
prepdataN <- LARprep(monthly, years = timepoints, class = c("Nomadic")) 
groupsdataN <- prepdataN[[1]] # extract the groups
time_pointsN <- prepdataN[[2]] # extract the time series of the data
lagged_ratesNprep <- LAR(groupsdataN,time_pointsN,31, identities = categories$TRANSMITTERID)
lagged_ratesNprep <- as.data.frame(lagged_ratesNprep)
lagged_ratesNprep$Class <- "Nomadic-Nomadic"

# Resident-Nomadic LSOR
prepdataRN <- LARprep(monthly, years = timepoints, class = c("Resident","Nomadic")) 
groupsdataRN <- prepdataRN[[1]] # extract the groups
time_pointsRN <- prepdataRN[[2]] # extract the time series of the data
lagged_ratesRNprep <- LAR(groupsdataRN,time_pointsRN,31, identities = categories$TRANSMITTERID)
lagged_ratesRNprep <- as.data.frame(lagged_ratesRNprep)
lagged_ratesRNprep$Class <- "Resident-Nomadic"

# Female-Female LSOR
prepdataF <- LARprep(monthly, years = timepoints, class = c("F")) 
groupsdataF <- prepdataF[[1]] # extract the groups
time_pointsF <- prepdataF[[2]] # extract the time series of the data
lagged_ratesFprep <- LAR(groupsdataF,time_pointsF,31, identities = categories$TRANSMITTERID)
lagged_ratesFprep <- as.data.frame(lagged_ratesFprep)
lagged_ratesFprep$Class <- "Female-Female"

# Nomadic-Nomadic LSOR
prepdataRF <- LARprep(monthly, years = timepoints, class = c("Resident", "F")) 
groupsdataRF <- prepdataRF[[1]] # extract the groups
time_pointsRF <- prepdataRF[[2]] # extract the time series of the data
lagged_ratesRFprep <- LAR(groupsdataRF,time_pointsRF,31, identities = categories$TRANSMITTERID)
lagged_ratesRFprep <- as.data.frame(lagged_ratesRFprep)
lagged_ratesRFprep$Class <- "Resident-Female"

# Nomadic-Female LSOR
prepdataNF <- LARprep(monthly, years = timepoints, class = c("Nomadic", "F")) 
groupsdataNF <- prepdataNF[[1]] # extract the groups
time_pointsNF <- prepdataNF[[2]] # extract the time series of the data
lagged_ratesNFprep <- LAR(groupsdataNF,time_pointsNF,31, identities = categories$TRANSMITTERID)
lagged_ratesNFprep <- as.data.frame(lagged_ratesNFprep)
lagged_ratesNFprep$Class <- "Nomadic-Female"

lag_ratesprep <- rbind(lagged_ratesRprep,lagged_ratesNprep,lagged_ratesRNprep,lagged_ratesFprep,lagged_ratesRFprep,lagged_ratesNFprep)

# Visulise the LSOR through time
ggplot()+
  geom_line(data = lag_ratesprep, aes(x = exp(lag_ratesprep[,1]), y = lag_ratesprep[,2], colour = Class))+
  ylim(0, 1) + ylab("Lagged spatial overlap rate") + xlab("Lag (1 month)")+
  labs(colour = "Movement strategy")+
  scale_x_continuous(breaks = seq(0,121, by = 12))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))

##############################################################################################################################
# Determine mean and SD of the expected random overlaps
##############################################################################################################################
# Note : This section requires random overlap matrices. For this, please consult the R scripts 1_Determine_monthly_HR_overlap.R
# and 3_Spatially_Explict_random_walk_model for the code required to generate these matrices. 
# The code provided below details how the null LSOR was calculated.
# Before use of the code below, minor changes need to be made to LINE 139 to reflect your own data formating

# First step, import the simulated overlaps between all of the individuals
# First import each simulations monthly overlaps into R in separate lists
simulated_overlap <- list() # create an empty list in which to store each of the simulated overlaps within

for(t in 1:100){
  files <- list.files(path = paste("Data/Overlap_Matrices/Simulations/Monthly/Overlap_Matrices/Simulation", t, sep = "-"),
                      full.names = T) # create a vector list with the names of all files
  #-----------------------------#
  # Remove all of the monthly overlaps before 2010-08
  if(t <= 9){
    files2 <- substr(files, 73,79 )
  } else if (t <= 99){
    files2 <- substr(files, 74,80 )
  } else {
    files2 <- substr(files, 75,81 )
  }
  files2 <- files2[files2 >= "2010-08"] # extract all of the months after this threshold date
  
  files <- files[as.numeric(length(files)-length(files2)+1):length(files)] # this line subsets files to only the months equal to  or after 2010-08
  
  simulate <- lapply(files, fread) # Import all files using lapply
  
  #-------------------------------------------------------------------------------------------------------#
  # Convert each of the imported data.frames into matrices
  for(i in 1:length(simulate)){
    simulate[[i]] <- as.data.frame(simulate[[i]]) # convert to data.frame
    colnames(simulate[[i]]) = simulate[[i]][1, ] # convert the first row back to column names
    simulate[[i]] = simulate[[i]][-1, ] # remove the first row with the names  
    rownames(simulate[[i]]) <- c() # remove the current row names
    simulate[[i]] <- tibble::column_to_rownames(simulate[[i]], var = "NA") # convert iD back into row.names
    simulate[[i]] <- as.matrix(simulate[[i]]) # convert to matrix
  }
  #assign(paste("Simulation",t, sep = "_"), simulate) # reassign the name of the simulate list to the dynamic name I would like
  # Place each of these imported simulations into another list
  # Put the simulated overlaps into another list along with the time series they came in with
  simulated <- list(simulate, files2)
  
  # place all of the simulations and their time lists into their own list
  simulated_overlap[[t]] <- simulated
}


# Extraction test
sim1 <- simulated_overlap[[100]]
sim1time <- sim1[[2]]
sim1overlaps <- sim1[[1]]

#============================================================================================================#
# Create a function to create the time series and run LAR prep of each simulation - this will allow paraellel processing to occur
#============================================================================================================#

simLAR <- function(i , population = T, Class = c("Resident")){
  # First; Extract the overlap matrices list and underlying movement data from their respective lists
  simulation <- simulated_overlap[[i]]
  Simoverlap <- simulation[[1]]
  
  timepoints <- simulation[[2]] # extract the time points of the simulation of interest
  timepoints <- sort(timepoints) # make sure they are in the correct temporal order before conducting LARprep
  
  
  # Run the LAR prep function to get the data ready
  prepdata <- LARprep(Simoverlap, years = timepoints, class = Class) 
  
  if(population == T){ 
    # Calculate the population level LSOR
    lagged_rates <- as.data.frame(LAR(prepdata[[1]],prepdata[[2]],30, identities = categories$TRANSMITTERID))
    names(lagged_rates) <- c("X", "Y") # Rename the columns 
    lagged_rates$Simulation <- i
    lagged_rates$Step <- seq(1,nrow(lagged_rates))
    #lagged_rates$Class <- Class
  } else{
    # Calculate the dyad specific LSOR
    lagged_rates <- LRA(prepdata[[1]],prepdata[[2]],30, identities = prepdata[[3]], output_style = 2)
    lagged_rates$Simulation <- i
  }
  
  return(lagged_rates)
}

####=============================================================================================================#####
# The below code calculates the null LSOR for each movement tactic and creates the final plots that were used in 
# Figure 4 of the manuscript

# Do low movement males first
registerDoParallel(10) # set the number of cores to use

lmigR <- foreach(i = 1:(length(simulated_overlap)), 
                 .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                               "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = "Resident")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORR <- data.table(Reduce(rbind, lmigR))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsR <- simLSORR[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]
meanOverlapsR$sdmin <- ifelse(meanOverlapsR$Mean - meanOverlapsR$Sd < 0, 0, meanOverlapsR$Mean - meanOverlapsR$Sd)
meanOverlapsR$sdmax <- meanOverlapsR$Mean + meanOverlapsR$Sd

p1 <- ggplot()+
  #geom_point(data = lagged_ratesL, aes(x = exp(lagged_ratesL[,1]), y = lagged_ratesL[,2]), col = "blue")+
  geom_line(data = lagged_ratesRprep, aes(x = (exp(lagged_ratesRprep[,1]))/12, y = lagged_ratesRprep[,2]))+
  geom_line(data = meanOverlapsR, aes(x = (exp(X))/12, y = Mean), col = "red",  linetype = "dashed")+
  geom_ribbon(data = meanOverlapsR, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = sdmax), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("") + xlab("")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  #geom_hline(yintercept = 0.047, linetype = "dotted", col = "red")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p1

# Do low movement males first
registerDoParallel(10) # set the number of cores to use

lmigN <- foreach(i = 1:(length(simulated_overlap)), 
                 .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                               "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = "Nomadic")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORN <- data.table(Reduce(rbind, lmigN))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsN <- simLSORN[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]
meanOverlapsN$sdmin <- ifelse(meanOverlapsN$Mean - meanOverlapsN$Sd < 0, 0, meanOverlapsN$Mean - meanOverlapsN$Sd)
meanOverlapsN$sdmax <- meanOverlapsN$Mean + meanOverlapsN$Sd

# remove the annoying gaps
meanOverlapsN[is.nan(meanOverlapsN)] <- 0
meanOverlapsN$Mean[is.nan(meanOverlapsN$Mean)]<-0
lagged_ratesNprep[,2][is.nan(lagged_ratesNprep[,2])]<-0

p2 <- ggplot()+
  geom_line(data = lagged_ratesNprep, aes(x = (exp(lagged_ratesNprep[,1]))/12, y = lagged_ratesHprep[,2]))+
  geom_line(data = meanOverlapsN, aes(x = (exp(X))/12, y = Mean), col = "red",  linetype = "dashed")+
  geom_ribbon(data = meanOverlapsN, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = sdmax), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("Lagged spatial overlap rate") + xlab("")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p2

# Do resident to nomadic individuals 
registerDoParallel(10) # set the number of cores to use

lmigRN <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                                "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = c("Resident", "Nomadic"))

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORRN <- data.table(Reduce(rbind, lmigRN))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsRN <- simLSORRN[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]
meanOverlapsRN$sdmin <- ifelse(meanOverlapsRN$Mean - meanOverlapsRN$Sd < 0, 0, meanOverlapsRN$Mean - meanOverlapsRN$Sd)
meanOverlapsRN$sdmax <- meanOverlapsRN$Mean + meanOverlapsRN$Sd

p3 <- ggplot()+
  geom_line(data = lagged_ratesRNprep, aes(x = (exp(lagged_ratesRNprep[,1]))/12, y = lagged_ratesRNprep[,2]))+
  geom_line(data = meanOverlapsRN, aes(x = (exp(X))/12, y = Mean), col = "red", linetype = "dashed")+
  geom_ribbon(data = meanOverlapsRN, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = sdmax), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p3

# Do females 
registerDoParallel(10) # set the number of cores to use

lmigF <- foreach(i = 1:(length(simulated_overlap)), 
                 .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                               "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = "F")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORF <- data.table(Reduce(rbind, lmigF))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsF <- simLSORF[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]

meanOverlapsF$sdmin <- ifelse(meanOverlapsF$Mean - meanOverlapsF$Sd < 0, 0, meanOverlapsF$Mean - meanOverlapsF$Sd)
meanOverlapsF$sdmax <- meanOverlapsF$Mean + meanOverlapsF$Sd

p4 <- ggplot()+
  geom_line(data = lagged_ratesFprep, aes(x = (exp(lagged_ratesFprep[,1]))/12, y = lagged_ratesFprep[,2]))+
  geom_line(data = meanOverlapsF, aes(x = (exp(X))/12, y = Mean), col = "red", linetype = "dashed")+
  geom_ribbon(data = meanOverlapsF, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = Mean + Sd), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("") + xlab("")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p4

# Do resident to females 
registerDoParallel(10) # set the number of cores to use

lmigRF <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                                "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = c("Resident","F"))

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORRF <- data.table(Reduce(rbind, lmigRF))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsRF <- simLSORRF[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]
meanOverlapsRF$sdmin <- ifelse(meanOverlapsRF$Mean - meanOverlapsRF$Sd < 0, 0, meanOverlapsRF$Mean - meanOverlapsRF$Sd)
meanOverlapsRF$sdmax <- meanOverlapsRF$Mean + meanOverlapsRF$Sd

p5 <- ggplot()+
  geom_line(data = lagged_ratesRFprep, aes(x = (exp(lagged_ratesRFprep[,1]))/12, y = lagged_ratesRFprep[,2]))+
  geom_line(data = meanOverlapsRF, aes(x = (exp(X))/12, y = Mean), col = "red", linetype = "dashed")+
  geom_ribbon(data = meanOverlapsRF, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = sdmax), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("") + xlab("")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p5

# Do nomadic to females 
registerDoParallel(10) # set the number of cores to use

lmigNF <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                                "asnipe","rlist","tibble","tidyverse")) %dopar%
  simLAR(i, population = T, Class = c("Nomadic", "F"))

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# reduce the list into a single data table
simLSORNF <- data.table(Reduce(rbind, lmigNF))

# Calculate the mean LSOR, sd and CI for each time step present
# Dyad specific
#meanOverlaps <- simLSOR[,.(Mean = mean(Y), Sd = sd(Y), CI = qnorm(0.975)*sd(RATE)/sqrt(10)), by = Simulation]
# Population
meanOverlapsNF <- simLSORNF[,.(Mean = mean(Y, na.rm = T), Sd = sd(Y, na.rm = T), CI = qnorm(0.975)*sd(Y)/sqrt(10)), by = X]
meanOverlapsNF$sdmin <- ifelse(meanOverlapsNF$Mean - meanOverlapsNF$Sd < 0, 0, meanOverlapsNF$Mean - meanOverlapsNF$Sd)
meanOverlapsNF$sdmax <- meanOverlapsNF$Mean + meanOverlapsNF$Sd

p6 <- ggplot()+
  geom_line(data = lagged_ratesNFprep, aes(x = (exp(lagged_ratesNFprep[,1]))/12, y = lagged_ratesNFprep[,2]))+
  geom_line(data = meanOverlapsNF, aes(x = (exp(X))/12, y = Mean), col = "red", linetype = "dashed")+
  geom_ribbon(data = meanOverlapsNF, aes(x = (exp(X))/12 ,ymin = sdmin, ymax = sdmax), fill = "grey50", alpha = 0.4)+
  ylim(0, 1) + ylab("") + xlab("Lag (1 year)")+
  scale_x_continuous(breaks = seq(0,5, by = 1), limits = c(0,5))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))
p6
#-----------------------------------------------------------------------------------------------------------#
# Compile the final figure together
figure <- ggarrange(p1,p4,p2,p5,p3,p6,
                    labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
                    ncol = 2, nrow = 3)
figure

# Final modifications to the figure were performed in Adobe Illustrator
