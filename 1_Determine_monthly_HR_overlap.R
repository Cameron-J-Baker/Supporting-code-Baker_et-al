# Author: Cameron J baker
# 03-12-2020
# Email: cameron.baker@uqconnect.edu.au
# Code for the analyses contained within the manuscript "Membership of crocodile social environments is dictated by male 
# philopatry" that is submitted for review to Behavioral Ecology
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
library(effects)
library(lme4)
library(emmeans)
library(performance)


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
# The output from the monthlyoverlap function are also saved to file to allow the below analyses to be completed without
# re-running the code above

#---------------------------------------------------------------------------------------------------------------------#
# Import all of the months into a list and then make each data.frame into a matrix
files <- list.files(path ="Output/lcUD-files", full.names = T) # create a vector list with the names of all files
files <- files[13:113] # only import the files from 2010-08 onwards
lmig <- lapply(files, fread) # Import all files using lapply

# Convert each of the files imported into matrices and move the TRANSMITTERIDs into the rownames and then create a dataframe with the IDs and thier overlaps
for(i in 1:length(lmig)){
  lmig[[i]] <- as.data.frame(lmig[[i]]) # convert to data.frame
  colnames(lmig[[i]]) = lmig[[i]][1, ] # convert the first row back to column names
  lmig[[i]] = lmig[[i]][-1, ] # remove the first row with the names  
  rownames(lmig[[i]]) <- c() # remove the current row names
  lmig[[i]] <- tibble::column_to_rownames(lmig[[i]], var = "NA") # convert iD back into row.names
  
  lmig[[i]] <- as.matrix(lmig[[i]])
  # set the upper half and diagonal of the matrix to NA
  lmig[[i]][upper.tri(lmig[[i]], diag = T)] <- NA 
  # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
  lmig[[i]] <- na.omit(as.data.table(as.table(lmig[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
  names(lmig[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns
  lmig[[i]]$Month <- substr(files[i], 36, 42)
  lmig[[i]] <- lmig[[i]][HRoverlap >= 0.01] # remove any dyds that are less than 1% to avoid including dyads that don't overlap
  #lmig[[i]] <- as.data.table( plyr::join(lmig[[i]], moveClasses, by = c("ID1", "ID2")))
  }

lmig2 <- plyr::compact(lmig)

#View(lmig2[[1]])

# Now convert the HR overlap matrices into a dataframe with the IDs of each individual
monthly_overlaps <- rbindlist(lmig2)

# Next, import the movement tactic of each animal and then attach it to hte HR overlaps for each one of the dyads
categories <- fread("Data/Metrics/Croc_categories_updated_V2.csv") # import the category of each individual

# Next create two dataframes, one for ID1 and then another for ID2 so that the movement tactics of each indivdual can be assigned
ID1_cat <- data.table(ID1 = categories$TRANSMITTERID,
                      ID1_tactic = categories$Class)
ID2_cat <- data.table(ID2 = categories$TRANSMITTERID,
                      ID2_tactic = categories$Class)

# join the two datasets to the overlaps dataframe
monthly_overlaps <- monthly_overlaps |>
  plyr::join(ID1_cat, by = "ID1") |>
  plyr::join(ID2_cat, by = "ID2")

# As there are some indivduals included who do not have a movement tactic assigned, use NA omit to remove them from the data set
monthly_overlaps <- na.omit(monthly_overlaps)

# Next, create the dyad and dyad tactic column
monthly_overlaps$Dyads <- paste0(monthly_overlaps$ID1, "-", monthly_overlaps$ID2)
monthly_overlaps$Tactic <- paste0(monthly_overlaps$ID1_tactic, "-", monthly_overlaps$ID2_tactic)

# Check how many different types of dyad tactics are present
unique(monthly_overlaps$Tactic)

# Recode the values so that only Low-Low, Low-High, Low-F, High-High, High-F and F-F are present
monthly_overlaps$Tactic <- ifelse(monthly_overlaps$Tactic == "Nomadic-Resident", "Resident-Nomadic", 
                                  ifelse(monthly_overlaps$Tactic == "Female-Resident", "Resident-Female", 
                                         ifelse(monthly_overlaps$Tactic == "Female-Nomadic", "Nomadic-Female", monthly_overlaps$Tactic)))
# Then remove the unncessary columns
monthly_overlaps$ID1_tactic <- monthly_overlaps$ID2_tactic <- NULL

#------------------------------------------------------------------#
# Assign each overlap to the month and then season (mating/non-mating) season they occurred
monthly_overlaps$Year <- substr(monthly_overlaps$Month, 1,4) # create a column for year so it can potentially be used as a random effect
# then redo the month column to only have the month of interest present
monthly_overlaps$Month <- ifelse(substr(monthly_overlaps$Month, 6,6) == 1,substr(monthly_overlaps$Month, 6,7), substr(monthly_overlaps$Month, 7,7) )
# Now create the season column. August - Novemeber is the Mating season
monthly_overlaps$Season <- ifelse(monthly_overlaps$Month < 8 & monthly_overlaps$Month > 11, "Non-mating", "Mating")

######################################################################################################################################################################
# The data is now ready to model

# First make sure that each variable is in the correct format
str(monthly_overlaps)

# Make each of the variables into factors
monthly_overlaps$Month <- as.factor(monthly_overlaps$Month)
monthly_overlaps$Dyads <- as.factor(monthly_overlaps$Dyads)
monthly_overlaps$Tactic <- as.factor(monthly_overlaps$Tactic)
monthly_overlaps$Season <- as.factor(monthly_overlaps$Season)

# Visualise the trends within the raw data
ggplot()+
  geom_boxplot(data = monthly_overlaps, aes(x = Tactic, y = HRoverlap, colour = Season))
ggplot()+
  geom_boxplot(data = monthly_overlaps, aes(x = Season, y = HRoverlap))

#----------------------------------# 
# create the model
model1 <- lmer(sqrt(HRoverlap) ~ Tactic * Season + (1|Dyads), data = monthly_overlaps)

check_model(model1) # check that the model meets all of the assumptions of linear regression

car::Anova(model1, test.statistic = "F") # Use the Anova function to determine what the significant values are

summary(model2)


emmeans(model1, list(pairwise ~ Tactic:Season), adjust = "tukey") # run a Tukey post hoc test using the estimated marginal means to see where the significant differences are

lsmeans::lsmeans(model1, pairwise~Tactic:Season, adjust="tukey") # just to validate things, also run a Tukey post hoc using the least square means as well

preditions <- as.data.table(predictorEffect("Tactic", model1)) # get the predictions and 95% confidence intervals from the model
preditions$Tactic <- factor(preditions$Tactic,levels = 
                       c("Nomadic-Female","Resident-Nomadic", "Female-Female","Resident-Female","Nomadic-Nomadic","Resident-Resident"))
monthly_overlaps$Tactic <- factor(monthly_overlaps$Tactic,levels = 
                                    c("Nomadic-Female","Resident-Nomadic", "Female-Female","Resident-Female","Nomadic-Nomadic","Resident-Resident"))

# Create Figure 3 of the manuscript
ggplot()+
  geom_point(data = monthly_overlaps, aes(x = Tactic, y = HRoverlap, shape = Season), 
             position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3), alpha = 0.7, colour = "grey")+
  geom_point(data = preditions, aes(x = Tactic, y = (fit)^2, shape = Season), size = 3, position = position_dodge(width = 0.6))+
  geom_errorbar(data = preditions, aes(x= Tactic, ymin= (lower)^2 , ymax= (upper)^2, shape = Season),# add in the standard deviation as error bars
                width=0.1, position = position_dodge(width = 0.6)) +
  ylab("Proportion of home range overlap (%)") + xlab("Crocodile movement tactic")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        legend.position = c(0.1,0.95))
