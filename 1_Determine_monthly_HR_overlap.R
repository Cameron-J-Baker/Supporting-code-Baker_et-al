# Author: Cameron J baker
# 03-12-2020
# Email: cameron.baker@uqconnect.edu.au
# Code for the analyses contained within the manuscript "Membership of crocodile social environments is dictated by male 
# philopatry" that is submitted for review to the Proceedings of the Royal Society B
# Aim : To determine the monthly proportion of spatial overlap between conspecifics and then how this is influenced by
# proximity to the mating season and movement tactic. 

# Load the required R packages
library(VTrack)
library(tidyverse)
library(sf)
library(raster)
library(mapview)
library(data.table)
library(rgdal)
library(rgeos)
library(foreach)
library(doParallel)
library(wicket)
library(ggmap)
library(spatstat)
library(FSA)


##########################################################################################################################
# Import the required datasets ###########################################################################################

# The raw acoustic detections
crocs <- fread("Data/Wenlock_array_crocodile_2019.csv")

# The movement tactics of each crocodile
croc_categories <- fread("Data/Crocodile_movement_tactics.csv")

# The location of each of the acoustic hydrophones
hydro_loc <- fread("Data/Acoustic_hydrophone_locations.csv")

# A value to convert WGS data to the UTM projection for the Wenlock and Ducie river systems
UTM <- CRS("+init=epsg:32754") 

# The line file of the Wenlock and Ducie river systems
Wen_River <- readOGR("Data/Line_file","Wenlock_Ducie_linefile") # Load the Wenlock spatiallines  data
Wen_River <- spTransform(Wen_River, UTM) # project to UTM
plot(Wen_River) # Visualise the line file

# Import the cost.raster for calculating the least cost HR estimates
cost.raster <- readRDS("Data/cost_raster.rds") # Import the cost raster 
plot(cost.raster)

########################################################################################################################
# Source the required functions ########################################################################################
source("R_functions/COAV2.R")
source("R_functions/SnapCOAtoRiver.R")
source("R_functions/lcDistance.R")
source("R_functions/lcUDoverlap.R")
source("R_functions/monthlyoverlap.R")

########################################################################################################################
# Determine the monthly HR overlap between conspecifics ################################################################

# First step is to convert the raw acoustic detection data into centers of activity
coa_crocs <- COAV2(tagdata = crocs, hydro_loc, Vtrack = TRUE, id = "TRANSMITTERID", timestep = 360)
coa_crocs <- data.frame(Reduce(rbind, coa_crocs))
coa_crocs <- SnapCOAtoRiver(coa_crocs, Wen_River, max_dist = 5000) # snap the COAs to the river

# Write a copy of the output to avoid having to replicate this in later scripts
write.csv("Outputs/Crocodile_COA.csv", row.names = F)

#======================================================================================================================#
# Prepare the data for use in the monthlyoverlap function to determine the monthly HR overlap of conspecifics ##########

# set the DATETIME column to be a posixct object
coa_crocs$DATETIME <- as.POSIXct(coa_crocs$DATETIME, format="%Y-%m-%d %H:%M:%S", tz = "Australia/Brisbane")

# Exclude all detections from before 2010-08 to prevent potential biases in VI estimation due to low sample size at the 
# start of the study duration
coa_crocs <- coa_crocs[DATETIME >= "2010-08-01"]

# create a list of the years individuals have been tracked for
coa_crocs$Year <- paste(year(coa_crocs$DATETIME), substr(coa_crocs$DATETIME, 6,7), sep = "-")

years <- sort(unique(coa_crocs$Year)) # extract a list of each of the unique months that crocodiles were detected for and sort it from smallest to largest

#======================================================================================================================#
# Determine the monthly HR overlap between conspecifics using the monthlyoverlap function ##############################

# TO decrease the run time, paraellel processing using the 'foreach' package has been implemented. Alternatively, the
# function could also be run using parLapply, lapply or a for loop
registerDoParallel(10) # set the number of cores to use

lmig <- foreach(i = 1:(length(years)), 
                .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
                              "adehabitatLT", "adehabitatHR", "rgeos", "rgdal","VTrack",
                              "geosphere","raster","gdistance","tidyverse")) %dopar%
  monthlyoverlap(i, coa_crocs, loc, cost.raster = cost.raster,  h = 300)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

# The output of this will be a list containing the HR overlap matrices of each of the months run. The matrices will be
# in order of the month that they were within the years vector, i.e., years[1] = lmig[1]

#---------------------------------------------------------------------------------------------------------------------#
# For each individual determine the mean HR overlap and number of conspecifics per month. This has been calculated for
# each combination of movement tactic to examine how the movement tactic of conspecifics influences the social 
# environment of an individual

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

#-----------------------------------------------------------------------#
# Create a lists of the Resident and high movement males IDs
ID <- unique(croc_categories$TRANSMITTERID) # All individuals
Resident <- unique((subset(croc_categories, Class == "Resident"))$TRANSMITTERID) # create a list of the Resident individuals
High <- unique((subset(croc_categories, Class == "High"))$TRANSMITTERID) # create a list of the Nomadic individuals
Females <- unique((subset(croc_categories, Class == "Female"))$TRANSMITTERID)

meanON <- list()
timepoints <- years
for(i in 1:length(lmig)) { # Do the following for each month
  
  monthDat <- lmig[[i]] # Extract the month of interest
  # Remove the individuals not found in croc_categories as they are not of interest for this study
  monthDat <- monthDat[rownames(monthDat) %in% ID, ]
  monthDat <- monthDat[ , colnames(monthDat) %in% ID]
  # Convert the matrix to a list
  # set the upper half and diagonal of the matrix to NA
  monthDat[upper.tri(monthDat, diag = T)] <- NA 
  # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
  monthDat <- na.omit(as.data.table(as.table(monthDat, na="", row.names=T, col.names=T))) # converts from matrix to table
  names(monthDat) <- c("ID1", "ID2","HRoverlap") # Rename the columns
  monthDat <- as.data.table( plyr::join(monthDat, moveClasses, by = c("ID1", "ID2")))
  monthDat <- monthDat[HRoverlap >= 0.01] # remove any dyds that are less than 1% to avoid including dyads that don't overlap
  # First go through and create subsets of the overlap data to only include the groups of interest
  Resident_overlap <- monthDat[Class == "ResidentResident"] # extract just the Resident movement males to all others
  Nomadic_overlap <- monthDat[Class == "NomadicNomadic"] # extract just the Nomadic movement males to all others
  ResidentFemList <- monthDat[Class == "ResidentF"]
  NomadicFemList <- monthDat[Class == "NomadicF"]
  femFemlist <- monthDat[Class == "FF"]
  ResidentNomadic_overlap <- monthDat[Class == "ResidentNomadic"]
  
  #-------------------------------------------------------------------------------------------------------------#
  # Next for each of the male individuals present, determine their mean overlap and number of overlapping individuals
  # with other males/females within each group 
  overlap_Resident <- list()
  for(u in 1:length(Resident)){
    indiv <- Resident_overlap[ID1 == Resident[u]| ID2 == Resident[u]]
    indivCat <- categories[TRANSMITTERID == Resident[u]]
    
    overlap_Resident[[u]] <- data.table(ID = Resident[u],
                                   Size = indivCat$Size,
                                   HRoverlap = mean(indiv$HRoverlap),
                                   No.Indivs = nrow(indiv),
                                   Class = "Resident",
                                   Sex = "MM",
                                   Month = timepoints[i])
  }
  ResidentMGroups <- data.table(Reduce(rbind, overlap_Resident))
  
  overlap_Nomadic <- list()
  for(t in 1:length(Nomadic)){
    indiv <- Nomadic_overlap[ID1 == Nomadic[t]| ID2 == Nomadic[t]]
    indivCat <- categories[TRANSMITTERID == Nomadic[t]]
    
    overlap_Nomadic[[t]] <- data.table(ID = Nomadic[t],
                                    Size = indivCat$Size,
                                    HRoverlap = mean(indiv$HRoverlap),
                                    No.Indivs = nrow(indiv),
                                    Class = "Nomadic",
                                    Sex = "MM",
                                    Month = timepoints[i])
  }
  NomadicMGroups <- data.table(Reduce(rbind, overlap_Nomadic))
  
  overlap_ResidentNomadic <- list()
  for(t in 1:length(Nomadic)){
    indiv <- ResidentNomadic_overlap[ID1 == Nomadic[t]| ID2 == Nomadic[t]]
    indivCat <- categories[TRANSMITTERID == Nomadic[t]]
    
    overlap_ResidentNomadic[[t]] <- data.table(ID = Nomadic[t],
                                       Size = indivCat$Size,
                                       HRoverlap = mean(indiv$HRoverlap),
                                       No.Indivs = nrow(indiv),
                                       Class = "ResidentNomadic",
                                       Sex = "MM",
                                       Month = timepoints[i])
  }
  ResidentNomadicMGroups <- data.table(Reduce(rbind, overlap_ResidentNomadic))
  
  overlap_ResidentF <- list()
  for(e in 1:length(Resident)){
    indiv <- ResidentFemList[ID1 == Resident[e]| ID2 == Resident[e]]
    indivCat <- categories[TRANSMITTERID == Resident[e]]
    
    overlap_ResidentF[[e]] <- data.table(ID = Resident[e],
                                    Size = indivCat$Size,
                                    HRoverlap = mean(indiv$HRoverlap),
                                    No.Indivs = nrow(indiv),
                                    Class = "ResidentF",
                                    Sex = "MF",
                                    Month = timepoints[i])
  }
  ResidentFGroups <- data.table(Reduce(rbind, overlap_ResidentF))
  
  overlap_NomadicF <- list()
  for(q in 1:length(Nomadic)){
    indiv <- NomadicFemList[ID1 == Nomadic[q]| ID2 == Nomadic[q]]
    indivCat <- categories[TRANSMITTERID == Nomadic[q]]
    
    overlap_NomadicF[[q]] <- data.table(ID = Nomadic[q],
                                     Size = indivCat$Size,
                                     HRoverlap = mean(indiv$HRoverlap),
                                     No.Indivs = nrow(indiv),
                                     Class = "NomadicF",
                                     Sex = "MF",
                                     Month = timepoints[i])
  }
  NomadicFGroups <- data.table(Reduce(rbind, overlap_NomadicF))
  
  overlap_F <- list()
  for(o in 1:length(Females)){
    indiv <- femFemlist[ID1 == Females[o]| ID2 == Females[o]]
    indivCat <- categories[TRANSMITTERID == Females[o]]
    
    overlap_F[[o]] <- data.table(ID = Females[o],
                                 Size = indivCat$Size,
                                 HRoverlap = mean(indiv$HRoverlap),
                                 No.Indivs = nrow(indiv),
                                 Class = "F",
                                 Sex = "FF",
                                 Month = timepoints[i])
  }
  FGroups <- data.table(Reduce(rbind, overlap_F))
  
  # Combine the above data.frames into a single one
  meanON[[i]] <- rbind(ResidentMGroups, NomadicMGroups,ResidentNomadicMGroups, ResidentFGroups, NomadicFGroups,FGroups)
}

crocOverlap_neighbours <- data.table(Reduce(rbind, meanON))

crocOverlap_neighbours <- na.omit(crocOverlap_neighbours)

# Write a copy to avoid having to rerun all of the above code again
write.csv(crocOverlap_neighbours, "Data/Output/monthly_overlap.csv")

#######################################################################################################################
# What is the influence of movement tactic and proximity to the mating season on HR overlap ##########################
#######################################################################################################################

# Import the monthly HR overlap for each individual
crocOverlap_neighbours <- fread("Data/Output/monthly_overlap.csv")

# Add a season column. August - November defined as the breeding season
crocOverlap_neighbours$Season <- ifelse(as.numeric(substr(crocOverlap_neighbours$Month, 6,7)) >= 8 & as.numeric(substr(crocOverlap_neighbours$Month, 6,7)) < 12, "Breeding", "Non-breeding")

# To account for how the different sex combinations may influence the amount of overlap between conspecifics, subset
# out all of the different possibilities
MMoverlap_neighbours <- crocOverlap_neighbours[Sex == "MM"]
MFoverlap_neighbours <- crocOverlap_neighbours[Sex == "MF"]
FFoverlap_neighbours <- crocOverlap_neighbours[Sex == "FF"]

###=================================================================================================================###
# For each of the different sex combinations, determine the mean proportion of HR overlap that an individual has with
# conspecifics for each of the two seasons

# Do the male - male combinations first

#-----------------------------------------------------------------------------#
# Subset out the Resident Nomadic combinations to run them independently to avoid it getting lumped in with the low individuals values
MMResNomB <- MMoverlap_neighboursB[Class == "ResidentNomadic"]
MMResNomNB <- MMoverlap_neighboursNB[Class == "ResidentNomadic"]
MMoverlap_neighboursB <- MMoverlap_neighboursB[Class != "ResidentNomadic"]
MMoverlap_neighboursNB <- MMoverlap_neighboursNB[Class != "ResidentNomadic"]
# Determine the mean overlap with males per month
MMoverlapB <- MMoverlap_neighboursB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                       Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                       Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
MMoverlapNB <- MMoverlap_neighboursNB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                         Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                         Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
MMoverlapRNB <- MMResNomB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                              Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                              Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
MMoverlapRNNB <- MMResNomNB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]

MMoverlap <- rbind(MMoverlapB,MMoverlapNB,MMoverlapRNB,MMoverlapRNNB) # combine the two

###-----------------------------------------------------------------------------------------------------------------###
# Male - female combinations

MFoverlap_neighboursB <- MFoverlap_neighbours[Season == "Breeding"]
MFoverlap_neighboursNB <- MFoverlap_neighbours[Season == "Non-breeding"]

# Determine the mean overlap with males per month
MFoverlapB <- MFoverlap_neighboursB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                       Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                       Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
MFoverlapNB <- MFoverlap_neighboursNB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                         Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                         Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
MFoverlap <- rbind(MFoverlapB,MFoverlapNB) # combine the two

###-----------------------------------------------------------------------------------------------------------------###
# Finally do the female - female combinations
FFoverlap_neighboursB <- FFoverlap_neighbours[Season == "Breeding"]
FFoverlap_neighboursNB <- FFoverlap_neighbours[Season == "Non-breeding"]

# Determine the mean overlap with males per month
FFoverlapB <- FFoverlap_neighboursB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                       Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                       Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
FFoverlapNB <- FFoverlap_neighboursNB[,.(Mean.overlap = mean(HRoverlap, na.rm = T), SE.overlap = plotrix::std.error(HRoverlap, na.rm = T), 
                                         Mean.indivs= mean(No.Indivs, na.rm = T), SE.indivs = plotrix::std.error(No.Indivs, na.rm = T),
                                         Class = Class[1], Season = Season[1], Length = mean(Length)), by = TRANSMITTERID]
FFoverlap <- rbind(FFoverlapB,FFoverlapNB) # combine the two

###----------------------------------------------------------------------------------------------------------------###
# Combine each of these datasets together to perform the analysis on
Overlaps <- rbind(MMoverlap, MFoverlap, FFoverlap)

# Set the order of class for plotting the data
Overlaps$Class <- factor(Overlaps$Class,levels = 
                           c("FemaleFemale","ResidentFemale","NomadicFemale","Resident","Nomadic","ResidentNomadic"))
# rename the classes for creating the plot of the manuscript
Overlaps$Class <- ifelse(Overlaps$Class == "FemaleFemale", "Female-Female",
                         ifelse(Overlaps$Class == "ResidentFemale", "Resident-Female",
                                ifelse(Overlaps$Class == "NomadicFemale", "Nomadic-Female",
                                       ifelse(Overlaps$Class == "Resident", "Resident-Resident",
                                              ifelse(Overlaps$Class == "Nomadic", "Nomadic-Nomadic", "Resident-Nomadic")))))
# Finally set the order of the different movement tactics to create the plot
Overlaps$Class <- factor(Overlaps$Class,levels = 
                           c("Nomadic-Female","Resident-Nomadic", "Female-Female","Resident-Female","Nomadic-Nomadic","Resident-Resident"))
# create the plot for Figure 3 of the manuscript. Final edits to the Figure where then made in Adobe Illustrator 
p1 <- ggplot()+
  geom_boxplot(data = Overlaps, aes(x = Class, y = Mean.overlap, fill = Season))+
  ylab("Proportion of home range overlap (%)")+
  xlab("Crocodile movement strategy")+
  #ylim(0,16)+
  scale_fill_manual(values = c("#4b7a2b","#c9dbb4"))+
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

# Check the homogentiy of the data
car::leveneTest(Mean.overlap~Class*Season, Overlaps) # Non-homogenity in the data, need to switch to non-parametric test

#------------------------------------------------------W
# To determine if there is an impact of class within season, sepearate Kruskal-Wallis tests were performed for each season
# breeding
breeding <- Overlaps[Season == "Breeding"]
kruskal.test(Mean.overlap~Class, data = breeding)
dunnTest(Mean.overlap~Class,
         data=breeding,
         method="bh")

#non-breeding
non_breeding <- Overlaps[Season == "Non-breeding"]
kruskal.test(Mean.overlap~Class, data = non_breeding)
dunnTest(Mean.overlap~Class,
         data=non_breeding,
         method="bh")

###===================================================####
# Now within classes, is there an impact of season?
#resid-resid
RR <- Overlaps[Class == "Resident-Resident"]
kruskal.test(Mean.overlap~Season, data = RR)

#nomad-nomad
NN <- Overlaps[Class == "Nomadic-Nomadic"]
kruskal.test(Mean.overlap~Season, data = NN)

#resid-Nomad
RN <- Overlaps[Class == "Resident-Nomadic"]
kruskal.test(Mean.overlap~Season, data = RN)

#female-female
FF <- Overlaps[Class == "Female-Female"]
kruskal.test(Mean.overlap~Season, data = FF)

#resid-resid
RF <- Overlaps[Class == "Resident-Female"]
kruskal.test(Mean.overlap~Season, data = RF)

#nomadic-female
NF <- Overlaps[Class == "Nomadic-Female"]
kruskal.test(Mean.overlap~Season, data = NF)

#===========================================#
# Does male movement strat influence overlap with females?
kruskal.test(Mean.overlap~Class, data = MFoverlap)


