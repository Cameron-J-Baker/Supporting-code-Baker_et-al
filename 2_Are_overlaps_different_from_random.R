# Author: Cameron J baker
# 03-12-2020
# Email: cameron.baker@uqconnect.edu.au
# Code for the analyses contained within the manuscript "Membership of crocodile social environments is dictated by male 
# philopatry" that is submitted for review to Behavioral Ecology
# Aim : To examine the whether the social environment of estuarine crocodiles is different from random chance
# This script examines the monthly stability of HR overlaps and the repeatability of these overlaps

# Load the required R packages
library(data.table)
library(lubridate)
library(foreach)
library(rptR)
library(tibble)
library(plyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(brms)

#### 1 : Examine the stability of the HR overlaps between months ######################################################

##########################################################################################################################
# Import the required datasets ###########################################################################################

# The movement tactics of each crocodile
croc_categories <- fread("Data/Crocodile_movement_tactics.csv")
### Repeat of the code in 1_Determine-monthly_HR_overlap to determine a list of all the different combinations of 
# movement tactic
ID <- unique(croc_categories$TRANSMITTERID) # create a vector of all the crocodile IDs

# To create a list of all the possible combinations, we first created a matrix with each individual and their movement
# tactic
moveCat<- function(i, Crocodile_movement_tactics){ 
  tt <- list()
  
  indiv <- Crocodile_movement_tactics[TRANSMITTERID == ID[i]]
  tt <- paste(indiv$Class,Crocodile_movement_tactics$Class, sep = "")
}

movMatrix <- foreach(i = 1:nrow(Crocodile_movement_tactics), .combine = "cbind") %do% moveCat(i, categories)
rownames(movMatrix) <- as.character(ID) # rename the columns to their respective CrocIDs
colnames(movMatrix) <- as.character(ID) # rename the rows to their respective CrocIDs
# Reduce to only one version of each column
movMatrix[movMatrix == "FemaleResident"] <- "ResidentFemale"
movMatrix[movMatrix == "FemaleNomadic"] <- "NomadicFemale"
movMatrix[movMatrix == "NomadicResident"] <- "ResidentNomadic"

# Convert the size matrix to a list
moveClasses <- na.omit(as.data.table(as.table(movMatrix, na="", row.names=T, col.names=T))) # converts from matrix to table
names(moveClasses) <- c("ID1", "ID2","Class") # Rename the columns 

#===============================================================================================================================#
# Import the monthly overlap matrices #################################################################################
files <- list.files(path ="Output/Overlap_Matrices/", full.names = T) # create a vector list with the names of all files

monthly <- lapply(files, fread) # Import all files using lapply
# To ensure that things are kept in the order that they were imported,create the month list from the names of each file
month <- substr(files,25,31)

# Convert each of the files imported into matrices and then to a data.frame of the individuals present
o <- 1
for(i in 1:length(monthly)){
  monthly[[i]] <- as.data.frame(monthly[[i]]) # convert to data.frame
  colnames(monthly[[i]]) = monthly[[i]][1, ] # convert the first row back to column names
  monthly[[i]] = monthly[[i]][-1, ] # remove the first row with the names  
  rownames(monthly[[i]]) <- c() # remove the current row names
  monthly[[i]] <- tibble::column_to_rownames(monthly[[i]], var = "NA") # convert iD back into row.names
  # Remove the individuals that are not included within the study
  monthly[[i]] <- monthly[[i]][rownames(monthly[[i]]) %in% ID, ]
  monthly[[i]] <- monthly[[i]][ , colnames(monthly[[i]]) %in% ID]
  
  monthly[[i]] <- as.matrix(monthly[[i]])
  
  # set the upper half and diagonal of the matrix to NA
  monthly[[i]][upper.tri(monthly[[i]], diag = T)] <- NA 
  # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
  monthly[[i]] <- na.omit(as.data.table(as.table(monthly[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
  names(monthly[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns 
  monthly[[i]]$Month <- month[o] # add a column for the month this was observed occurring in
  o <- o + 1
}

overlaps <- rbindlist(monthly) # Merge each of the matrices in the list into a single dataframe

head(overlaps)

# Next step: remove any dyads that have an overlap less than 0.01 - including them would skew the social environments of each
# individual
overlaps <- overlaps[HRoverlap >= 0.01]
# Join the two dataframes together
overlaps <- plyr::join(overlaps, moveClasses, by = c("ID1", "ID2"))
head(overlaps)

##### Import the function to determine the monthly stability of HR overlaps ###########################################
source("R_functions/socialStability.R")

### For each combination of movement tactics, determine the mean and SD of the number of conspecifics that individuals
# maintain overlaps with between months

# First do it for the overall combination of all tactics
Socialenvironment <- socialStability(overlaps, categories = croc_categories, Class = "Overall")
mean(Socialenvironment$Prop_last_month, na.rm = T)
sd(Socialenvironment$Prop_last_month, na.rm = T)

# Resident-Resident
ResRes <- socialStability(overlaps, categories = croc_categories,Class = "ResRes")
mean(ResRes$Prop_last_month, na.rm = T)
sd(ResRes$Prop_last_month, na.rm = T)

# Nomadic-Nomadic
NomNom <- socialStability(overlaps, categories = croc_categories,Class = "NomNom")
mean(NomNom$Prop_last_month, na.rm = T)
sd(NomNom$Prop_last_month, na.rm = T)

# Resident - Nomadic
ResNom <- socialStability(overlaps, categories = croc_categories,Class = "ResNom")
mean(ResNom$Prop_last_month, na.rm = T)
sd(ResNom$Prop_last_month, na.rm = T)

# Female - Female
FemFem <- socialStability(overlaps, categories = croc_categories,Class = "FemFem")
mean(FemFem$Prop_last_month, na.rm = T)
sd(FemFem$Prop_last_month, na.rm = T)

# Resident - Female
ResFem <- socialStability(overlaps, categories = croc_categories,Class = "ResFem")
mean(ResFem$Prop_last_month, na.rm = T)
sd(ResFem$Prop_last_month, na.rm = T)

# Nomadic - Female
NomFem <- socialStability(overlaps, categories = croc_categories,Class = "NomFem")
mean(NomFem$Prop_last_month, na.rm = T)
sd(NomFem$Prop_last_month, na.rm = T)
# Test to see if the function works

#######################################################################################################################
# 2: Examine the repeatability of the proportion of HR overlap between months ############################################
#######################################################################################################################

#====================================================================================================================#
# Impoty the monthly overlap matrices
#====================================================================================================================#
files <- list.files(path ="Data/Overlap_Matrices/lcUD_Monthly", full.names = T) # create a vector list with the names of all files
#files <- files[-119] # remove the last item as it is another folder within the directory
files <- files[13:113] # only import the files from 2010-08 onwards

monthly <- lapply(files, fread) # Import all files using lapply

# Convert each of the files imported into matrices and move the TRANSMITTERIDs into the rownames
for(i in 1:length(monthly)){
  monthly[[i]] <- as.data.frame(monthly[[i]]) # convert to data.frame
  colnames(monthly[[i]]) = monthly[[i]][1, ] # convert the first row back to column names
  monthly[[i]] = monthly[[i]][-1, ] # remove the first row with the names  
  rownames(monthly[[i]]) <- c() # remove the current row names
  monthly[[i]] <- tibble::column_to_rownames(monthly[[i]], var = "NA") # convert iD back into row.names
  
  monthly[[i]] <- as.matrix(monthly[[i]])}

# retrieve the month and sequence that the monthly file is in by using substring to retrieve the names from
# the files list
timepoints <- substr(files,25,31)

#===============================================================================================================#
# Import the movement, sex and size metrics
#===============================================================================================================#
croc_categories <- fread("Data/Crocodile_movement_tactics.csv")

# Import the required function to prepare the matrix data to examine the repeatability of HR overlaps
source("R_functions/rptPrep.R")

###############################################################################################################
# Prepare the raw overlap matrices into the format required to complete the repeatability analyses
# First step, create two lists with all the low movement and femaels
Res <- unique((subset(categories, Class == "Resident"))$TRANSMITTERID) # create a list of the low/female individuals
Nom <- unique((subset(categories, Class == "Nomadic"))$TRANSMITTERID) # create a list of the Nom/female individuals
Fem <- unique((subset(categories, Class == "Female"))$TRANSMITTERID) # create a list of the female individuals

# Second, create a copy of the monthly for both low and Nom movement individuals
RROverlap <- monthly
NNOverlap <- monthly
FFOverlap <- monthly
RFOverlap <- monthly
NFOverlap <- monthly
RNOverlap <- monthly

# Using the lists created above subset out the low and Nom movement individuals from their respective categories
for(i in 1:length(RROverlap)){
  RROverlap[[i]] <- RROverlap[[i]][rownames(RROverlap[[i]]) %in% Res, ]
  RROverlap[[i]] <- RROverlap[[i]][ , colnames(RROverlap[[i]]) %in% Res]
}
for(i in 1:length(NNOverlap)){
  NNOverlap[[i]] <- NNOverlap[[i]][rownames(NNOverlap[[i]]) %in% Nom, ]
  NNOverlap[[i]] <- NNOverlap[[i]][ , colnames(NNOverlap[[i]]) %in% Nom]
}
for(i in 1:length(FFOverlap)){
  FFOverlap[[i]] <- FFOverlap[[i]][rownames(FFOverlap[[i]]) %in% Fem, ]
  FFOverlap[[i]] <- FFOverlap[[i]][ , colnames(FFOverlap[[i]]) %in% Fem]
}
for(i in 1:length(RFOverlap)){
  RFOverlap[[i]] <- RFOverlap[[i]][rownames(RFOverlap[[i]]) %in% Res, ]
  RFOverlap[[i]] <- RFOverlap[[i]][ , colnames(RFOverlap[[i]]) %in% Fem]
}
for(i in 1:length(NFOverlap)){
  NFOverlap[[i]] <- NFOverlap[[i]][rownames(NFOverlap[[i]]) %in% Nom, ]
  NFOverlap[[i]] <- NFOverlap[[i]][ , colnames(NFOverlap[[i]]) %in% Fem]
}
for(i in 1:length(RNOverlap)){
  RNOverlap[[i]] <- RNOverlap[[i]][rownames(RNOverlap[[i]]) %in% Res, ]
  RNOverlap[[i]] <- RNOverlap[[i]][ , colnames(RNOverlap[[i]]) %in% Nom]
}

# Then run the list of matrices through the rptPrep function to create an overlap dataframe
RROverlap <- rptPrep(RROverlap, years = timepoints, dyad = T)
NNOverlap <- rptPrep(NNOverlap, years = timepoints, dyad = T)
FFOverlap <- rptPrep(FFOverlap, years = timepoints, dyad = T)
RFOverlap <- rptPrep(RFOverlap, years = timepoints, dyad = T)
NFOverlap <- rptPrep(NFOverlap, years = timepoints, dyad = T)
RNOverlap <- rptPrep(RNOverlap, years = timepoints, dyad = T)

#### Determine the observed level of repeatability ####################################################################
# Resident-resident males
RROverlap <- plyr::join(RROverlap, size, by = "ID") # Join the overlap and size datasets
RROverlap$Month <- substr(RROverlap$Month, 6,7)

RR.brm <- brm(HRoverlap ~ 1+ (1 | ID), data = RROverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
RR.brm <- add_criterion(RR.brm, "waic")
summary(RR.brm)
var.animal_id <- posterior_samples(RR.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(RR.brm)$"sigma"^2
rep1 <- mean(var.animal_id/(var.animal_id + var.res))

# Nomadic-nomadic males
NNOverlap <- plyr::join(NNOverlap, size, by = "ID") # Join the overlap and size datasets
NN.brm <- brm(HRoverlap ~ (1 | ID), data = NNOverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
NN.brm <- add_criterion(NN.brm, "waic")
summary(NN.brm)
var.animal_id <- posterior_samples(NN.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(NN.brm)$"sigma"^2
rep2 <- mean(var.animal_id/(var.animal_id + var.res))

# Resident-nomadic males
RNOverlap <- plyr::join(RNOverlap, size, by = "ID") # Join the overlap and size datasets
RN.brm <- brm(HRoverlap ~ (1 | ID), data = RNOverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
RN.brm <- add_criterion(RN.brm, "waic")
summary(RN.brm)
var.animal_id <- posterior_samples(RN.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(RN.brm)$"sigma"^2
rep3 <- mean(var.animal_id/(var.animal_id + var.res))

# Female-fenmale 
FFOverlap <- plyr::join(FFOverlap, size, by = "ID") # Join the overlap and size datasets
FF.brm <- brm(HRoverlap ~ (1 | ID), data = FFOverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
FF.brm <- add_criterion(FF.brm, "waic")
summary(FF.brm)
var.animal_id <- posterior_samples(FF.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(FF.brm)$"sigma"^2
rep4 <- mean(var.animal_id/(var.animal_id + var.res))

# Resident males -females
RFOverlap <- plyr::join(RFOverlap, size, by = "ID") # Join the overlap and size datasets
RF.brm <- brm(HRoverlap ~ (1 | ID), data = RFOverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
RF.brm <- add_criterion(RF.brm, "waic")
summary(RF.brm)
var.animal_id <- posterior_samples(RF.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(RF.brm)$"sigma"^2
rep5 <- mean(var.animal_id/(var.animal_id + var.res))

# Nomadic males-females
NFOverlap <- plyr::join(NFOverlap, size, by = "ID") # Join the overlap and size datasets
NF.brm <- brm(HRoverlap ~ (1 | ID), data = NFOverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
              inits = "random", cores = my.cores, seed = 12345)
NF.brm <- add_criterion(NF.brm, "waic")
summary(NF.brm)
var.animal_id <- posterior_samples(NF.brm)$"sd_ID__Intercept"^2
var.res <- posterior_samples(NF.brm)$"sigma"^2
rep6 <- mean(var.animal_id/(var.animal_id + var.res)) 

#==============================================================================================================#
# Create a dataframe with the observed results
Observed <- data.table(Simulation = "Observed",
                       Movement_strats = c("LL", "HH", "LH", "FF", "LF", "HF"),
                       ID_rpt = c(rep1$R[2,1],rep2$R[2,1],rep3$R[2,1],rep4$R[2,1],rep5$R[2,1],rep6$R[2,1]),
                       Size_rpt = c(rep1$R[2,2],rep2$R[2,2],rep3$R[2,2],rep4$R[2,2],rep5$R[2,2],rep6$R[2,2]))

Observed <- data.table(Simulation = "Observed",
                       Movement_strats = c("LL", "HH", "LH", "FF", "LF", "HF"),
                       ID_rpt = c(rep1,rep2,rep3,rep4,rep5,rep6))

write.csv(Observed,"Data/Observed_repeatability.csv")

###############################################################################################################
# Determine the ICC random distribution and then determine the significance of the above repeatabilities
###############################################################################################################
# Note : This section requires random overlap matrices. For this, please consult the R scripts 1_Determine_monthly_HR_overlap.R
# and 3_Spatially_Explict_random_walk_model for the code required to generate these matrices. 
# The code provided below details how the random ICCs were calculated and how the P-values for were generated. Before
# use of the code below, minor changes need to be made to LINE 310 to reflect your own data formating

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

###############################################################################################################
# The below function calculates the random ICC distribution from the 100 simulations for each specific category

simICC <- function(i, Class = c("RR","NN","RN","FF","RF","NF")){
  
  # First the function susbets out the specific simulation of interest
  simulation <- simulated_overlap[[i]]
  Simoverlap <- simulation[[1]] # create a list with just the simulated overlaps
  timepoints <- simulation[[2]] # extract the time points of the simulation of interest
  timepoints <- sort(timepoints) # make sure they are in the correct temporal order before conducting LARprep
  
  #----------------------------------------------------------------------------------------------------------#
  # Subset the simulated dataset down to the movement category combination of choice
  if(Class == "RR"){
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Res, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Res]} 
  } else if(Class == "NN"){
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Nom, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Nom]} 
  } else if(Class == "RN"){
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Res, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Nom]} 
  } else if(Class == "FF"){
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Fem, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Fem]} 
  } else if(Class == "RF"){
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Res, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Fem]} 
  } else {
    for(t in 1:length(Simoverlap)){
      Simoverlap[[t]] <- Simoverlap[[t]][rownames(Simoverlap[[t]]) %in% Nom, ]
      Simoverlap[[t]] <- Simoverlap[[t]][ , colnames(Simoverlap[[t]]) %in% Fem]} 
  }
  Simoverlap <- rptPrep(Simoverlap, years = timepoints, dyad = T) # run the rptPrep function to prepare the data
  
  # Run the rpt function to calculate the repeatability
  Simoverlap <- plyr::join(Simoverlap, size, by = "ID") # Join the overlap and size datasets
  m1.brm <- brm(HRoverlap ~ (1 | ID), data = Simoverlap, warmup = 500, iter = 3000, thin = 2, chains = 2,
                inits = "random", cores = 1, seed = 12345)
  m1.brm <- add_criterion(m1.brm, "waic")
  var.animal_id <- posterior_samples(m1.brm)$"sd_ID__Intercept"^2
  var.res <- posterior_samples(m1.brm)$"sigma"^2
  reps <- mean(var.animal_id/(var.animal_id + var.res))
  
  
  # Extract the information that is required
  output <- data.table(Simulation = i,
                       Movement_strats = Class,
                       ID_rpt = reps)
  
}

testing <- simICC(1, Class = "NF")

#==============================================================================================================#
# Now that the function works, calculate the expected repeatabilities for each movement strategy combinations

#-------------------------------------------------------------------------------------------------------------#
# Resident-resident
registerDoParallel(10) # set the number of cores to use

lmigRR <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "RR")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
simRptRR <- data.table(Reduce(rbind, lmigRR))

#-------------------------------------------------------------------------------------------------------------#
# nomadic-nomadic
registerDoParallel(10) # set the number of cores to use

lmigNN <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "NN")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
simRptNN <- data.table(Reduce(rbind, lmigNN))
View(simRptNN)
#-------------------------------------------------------------------------------------------------------------#
# Resident-nomadic
registerDoParallel(10) # set the number of cores to use

lmigRN <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "RN")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
simRptRN <- data.table(Reduce(rbind, lmigRN))

#-------------------------------------------------------------------------------------------------------------#
# female-female
registerDoParallel(10) # set the number of cores to use

lmigFF <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "FF")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
simRptFF <- data.table(Reduce(rbind, lmigFF))

#-------------------------------------------------------------------------------------------------------------#
# Resident-female
registerDoParallel(10) # set the number of cores to use

lmigRF <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "RF")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
simRptRF <- data.table(Reduce(rbind, lmigRF))

#-------------------------------------------------------------------------------------------------------------#
# nomadic-female
registerDoParallel(10) # set the number of cores to use

lmigNF <- foreach(i = 1:(length(simulated_overlap)), 
                  .packages = c("data.table","rptR", "brms")) %dopar%
  simICC(i, Class = "NF")

stopImplicitCluster() # stop the cluster to return memory and resources to other systems
# reduce the list into a single data table
lmigNF <- list()

for(i in 1:100){
  lmigNF[[i]] <- simICC(i, Class = "NF")
}

simRptNF <- data.table(Reduce(rbind, lmigNF))

#-------------------------------------------------------------------------------------------------------------#
# Combine all of the repeatability estimates into a single data frame

allRptdist <- rbind(simRptRR,simRptNN,simRptRN,simRptFF,simRptRF,simRptNF)

# write a copy of this data frame to save the data and not have to run these functions multiple times
write.csv(allRptdist, "Data/Simulated/Repeatabilities_bays_prop.csv")

###############################################################################################################
# Determine the significance of each of the repeatability values

# P - value calculated as the proportion of times that the expected repeatabilities are above the observed
Observed$Pvalue <- NA # create an empty column to save the P-value results to
Observed$ExpectedMin <- NA
Observed$ExpectedMax <- NA
for(i in 1:6){
  OBS <- Observed[i,] # subset the observed data to only a single row
  expected <- allRptdist[Movement_strats == OBS$Movement_strats]
  expected$Prop <- ifelse(expected$ID_rpt > OBS$ID_rpt, 1, 0)
  Observed$Pvalue[i] <- sum(expected$Prop)/100
  Observed$ExpectedMin[i] <- min(expected$ID_rpt)
  Observed$ExpectedMax[i] <- max(expected$ID_rpt)
}
