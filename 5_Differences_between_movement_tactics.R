# Author: Cameron J baker
# 03-12-2020
# Email: cameron.baker@uqconnect.edu.au
# Code for the analyses contained within the manuscript "Membership of crocodile social environments is dictated by male 
# philopatry" that is submitted for review to the Proceedings of the Royal Society B
# Aim : To compare the mean monthly distance travelled, home range size and proportion of home range overlap with the
# previous month between each of the three movement tactics

# Load the required R packages
library(data.table)
library(ggplot2)
library(ggpubr)
library(foreach)
library(doParallel)
library(sp)
library(rgeos)
library(rgdal)
library(raster)
library(lubridate)
library(gdistance)
library(maptools)
library(lme4)
library(performance)
library(effects)

#######################################################################################################################
#### 1 : Examine how the mean monthly distance moved varies between each movement tactic
#######################################################################################################################

# First import the COA of the population, as this is what will be simulated by the model
coa_crocs <- fread("Outputs/Crocodile_COA.csv")
# set the DATETIME column to be a posixct object
coa_crocs$DATETIME <- as.POSIXct(coa_crocs$DATETIME, format="%Y-%m-%d %H:%M:%S", tz = "Australia/Brisbane")

# The movement tactics of each crocodile
croc_categories <- fread("Data/Crocodile_movement_tactics.csv")

# Subset out the individuals from the raw data who do not have the required data to calculate home ranges
coa_crocs <- subset(coa_crocs, (coa_crocs$TRANSMITTERID %in% croc_categories$TRANSMITTERID))#remove once the lcUDs have been completed

# create a vector with the IDs of each individual
myids <- unique(croc_categories$TRANSMITTERID)

# Define UTM to project data into for calculation of home ranges
UTM <- CRS("+init=epsg:32754") #a value to convert data to the UTM projection


# Import the raster of the Wenlock river to calculate distances along
Wen_Raster <- raster("Data/WenlockRaster.asc")
plot(Wen_Raster)
# Convert the raster to a transition object - required for the river distance algorythm
tr <- transition(Wen_Raster, transitionFunction=mean, directions=8) # Create a Transition object from the raster

#####===========================================================================================================######
# Import the required functions

source("R_functions/consecutiveDistances.R")

#-------------------------------------------------------------------------------------------------------------#
# Run the above function through for each individual
registerDoParallel(10) # set the number of cores to use

lmig <- foreach(i = 1:(length(myids)), 
                .packages = c("data.table","gdistance","plyr", "sp")) %dopar%
  consecDistTravelled(i, coa_crocs, tr)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

CrocDistances<- data.table(Reduce(rbind, lmig))

write.csv(CrocDistances, "Outputs/distances_moved.csv", row.names = F) # write a copy to avoid having to run that again!
#=============================================================================================================#
# Combine the distances dataset to the categories one to assign each of the crocodiles to their categories

CrocDistances <- plyr::join(CrocDistances, croc_categories, by = "TRANSMITTERID")

# Create a new column to assign the month and year of the distances moved
CrocDistances$Month <- paste(year(CrocDistances$DATETIME), substr(CrocDistances$DATETIME, 6,7), sep = "-")

# Calculate the total distance travelled per month for each individual
monthlyDists <- list()

for(i in 1:length(myids)){
  indiv <- CrocDistances[TRANSMITTERID == myids[i]]
  
  monthlyDists[[i]] <-  indiv[,.(Total.dist = sum(rDISTANCE, na.rm = T), Class = Class[1], 
                                 TRANSMITTERID = myids[i]), by = Month]
}

monthlyDists <- data.table(Reduce(rbind, monthlyDists))

# conduct a linear model on the influence of movement tactic on the distance moved per month
# first visualise the data
ggplot()+
  geom_boxplot(data = monthlyDists, aes(x = Class, y = Total.dist))

# check the normality of the distances traveled
hist(sqrt(monthlyDists$Total.dist + 0.01))

model1 <- lmer(sqrt(Total.dist+0.01) ~ Class + (1|TRANSMITTERID), data = monthlyDists)

check_model(model1)

car::Anova(model1, test.statistic = "F")

lsmeans::lsmeans(model1, pairwise~Class, adjust="tukey") 

# Get the model predictions and then plot them to create panel one of supplementary figure 2
Dist_travelled <- as.data.table(predictorEffect("Class", model1))

#=============================================================================================================#
# Create the plot to visulise what is occurring

p1 <- ggplot()+
  geom_point(data = monthlyDists, aes(x = Class, y = Total.dist), position = "jitter", alpha = 0.3, colour = "grey50")+
  geom_point(data = Dist_travelled, aes(x = Class, y = (fit)^2))+
  geom_errorbar(data = Dist_travelled, aes(x= Class, ymin= (lower)^2, ymax= (upper)^2),
                width=.2,                    # Width of the error bars
                position=position_dodge(.0))+
  xlab("")+
  ylab("Monthly distance travelled (km)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.84,0.9))
p1


#######################################################################################################################
# 2 : examine the mean monthly home range size of each movement tactic ################################################
#====================================================================================================================================#
# First extract the monthly HR size of each individual from the UDs gerenated in 1_Determine_HR_overlap

# Import the UDs for each month
files <- list.files(path ="Outputs/lcUD-files/", full.names = T) # create a vector list with the names of all files
lmig <- lapply(files, readRDS) # Import all files using lapply

# To ensure that things are kept in the order that they were imported,create the month list from the names of each file
months <- substr(files, 19,25)

monthlyHR <- list()

for(i in 1:length(months)){
  # Extract the month of interest
  lmigM <- lmig[[i]]
  # Extract the IDs of the individuals present
  IDs <- names(lmigM)
  lcUD95 <- NULL
  lcUD50 <- NULL
  for(p in 1:length(IDs)){
    indiv <- lmigM[[p]]
    lcUD95[p] <- indiv$kernel.areas[,2]
    lcUD50[p] <- indiv$kernel.areas[,1]
  }
  monthlyHR[[i]] <- data.table(TRANSMITTERID = IDs,
                               lcUD95 = lcUD95/1000000 , # convert to km2
                               lcUD50 = lcUD50/1000000 , # convert to km2
                               Month = months[i])
}

monthlyHR <- data.table(Reduce(rbind, monthlyHR))

write.csv(monthlyHR, "Data/lc_HR/monthlyHRsize_20191204.csv", row.names = F)

# Join the home range size of each individual with their movement tactic
monthlyHR <- plyr::join(monthlyHR, croc_categories, by = "TRANSMITTERID")

# Remove all of the individual who are not included in the analysis that could not be classified into a movement class
monthlyHR <- na.omit(monthlyHR)

monthlyHR <- monthlyHR[Month >= "2010-08"]

### Run the analysis but using a linear mixed effect model instead of a Kuralis wallis test below
# first visualise the data
ggplot()+
  geom_boxplot(data = monthlyHR, aes(x = Class, y = lcUD95))

# check the normality of the distances traveled
hist(log(monthlyHR$lcUD95))

model2 <- lmer(log(lcUD95) ~ Class + (1|TRANSMITTERID), data = monthlyHR)

check_model(model2)

car::Anova(model2, test.statistic = "F")
anova(model2)

lsmeans::lsmeans(model2, pairwise~Class, adjust="tukey") 

# Get the model predictions and then plot them to create panel two of supplementary figure 2
HR_size <- as.data.table(predictorEffect("Class", model2))

#=============================================================================================================#
# Create the plot to visulise what is occurring

p2 <- ggplot()+
  geom_point(data = monthlyHR, aes(x = Class, y = lcUD95), position = "jitter", alpha = 0.3, colour = "grey50")+
  geom_point(data = HR_size, aes(x = Class, y = exp(fit)))+
  geom_errorbar(data = HR_size, aes(x= Class, ymin= exp(lower), ymax= exp(upper)),
                width=.2,                    # Width of the error bars
                position=position_dodge(.0))+
  xlab("Movement tactic")+
  ylab("Monthly home range size " ~(km^2))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = c(0.84,0.9))
p2

###---------------------------------------------------------------------------------------------------------------####
# Create Figure 2 of the manuscript by combining the last two ggplots created

figure <- ggarrange(p1,p2,
                    labels = c("a)", "b)"),
                    ncol = 1, nrow = 2)
#figure <- annotate_figure(figure, bottom = text_grob("Month"))
figure

# Final adjustments to the figures were completed in Adobe Illustrator

#######################################################################################################################
##### 3 : Does movement tactic influence an individuals proportion of HR overlap with their previous HR ###############

# This section used the list of lcUD estimates and the Crocodile movement tactic dataframe that were imported previously

#=====================================================================================================================#
# First, create a dataframe that lists all of the months that each individual have a HR generated for
monthIDs <- list()
for(i in 1:length(lmig)){
  
  monthIDs[[i]] <- data.table(ID = names(lmig[[i]]),
                              Month = substr(files[i],40,46),
                              lmigPos = i)
  
}

crocUDmonths <- data.table(Reduce(rbind, monthIDs))

myIDs <- unique(crocUDmonths$ID)

#####-------------------------------------------------------------------------------------------------------------####
# Import the functions required
source("R_functions/lcUDoverlap.R")

#=============================================================================================================#
# Create a function to subset out the UDs of individuals and then determines their mean and SD of their overlap with
# themselves between months

indivOverlap <- function(t){
  #crocID <- myIDs[t]
  # First step, subset crocUDmonths down to the individual that is being focused on
  indiv <- crocUDmonths[ID == myIDs[t]]
  
  # Next create a new list and then use a for loop to place all of the monthly KUDs of each crocodile into it
  indiv_lmig <- list()
  
  for(i in 1: nrow(indiv)){
    
    month <- lmig[[indiv$lmigPos[i]]] # subset lmig to the month of interest
    
    IDs <- names(month) # work out the names of the individuals present
    for(p in 1:length(IDs)){ # go through the list of IDs present and only pull out the data of the animal of interest
      if(IDs[p] == myIDs[t]){
        indiv_lmig[[i]] <- month[[p]] 
      }
    }
  }
  names(indiv_lmig) <- indiv$Month # Name each of the points in the list to their corresponding month 
  
  #----------------------------------------------------------------------------------------------------------#
  # Now that a list file has been created containing only an individuals UDs, run it through the overlap function
  # to calculate the overlap matrix
  
  indiv_overlap <- lcUDoverlap(indiv_lmig, percent = 95, method = "VI")
  
  delta <- row(indiv_overlap) - col(indiv_overlap) # create another matrix the same size as the indiv_overlap one, but with
  # it counting out from the diagonal
  indiv_overlap[delta < 1 | delta > 1] <- NA # set all values except those directly next to the diagonal to NA
  
  indiv_overlap <- na.omit(as.data.table(as.table(indiv_overlap, na="", row.names=T, col.names=T))) # converts from matrix to table
  names(indiv_overlap) <- c("Month1", "Month2","HRoverlap") # Rename the columns
  
  # calculate the mean and sd of HR with themselves and then return the values
  output <- data.table(TRANSMITTERID = myIDs[t],
                       meanHRoverlap = mean(indiv_overlap$HRoverlap),
                       sdHRoverlap = sd(indiv_overlap$HRoverlap))
  return(output)
}

# test the above function
#indivOverlap(17)

######################################################################################################################
# Run the above function in paraellel to calculate the mean HR values for each individual with themselves
#registerDoParallel(5) # set the number of cores to use

#HRoverlaps <- foreach(i = 1:(length(myIDs)), 
#               .packages = c("data.table","lubridate","plyr", "sp", "maptools", 
#                              "adehabitatLT", "adehabitatHR", "rgeos", "rgdal","VTrack",
#                             "geosphere","raster","gdistance","tidyverse")) %dopar%
# indivOverlap(i)

#stopImplicitCluster() 

#CrocMetrics <- data.frame(Reduce(rbind, HRoverlaps))

HRoverlaps <- list()

for(e in 1:length(myIDs)){
  HRoverlaps[[e]] <- indivOverlap(e)
}
CrocMetrics <- data.frame(Reduce(rbind, HRoverlaps))

#====================================================================================================================#
# Merge crocMetrics with categories then plot what differences are present or not. 

CrocMetrics <- plyr::join(categories, CrocMetrics, by = "TRANSMITTERID")

p3 <- ggplot()+
  geom_boxplot(data = CrocMetrics, aes(x = Class, y = meanHRoverlap))+
  ylim(0, 1)+
  xlab("Movement strategy")+
  ylab("Proportion of home range overlap (%)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        #axis.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position='none')

p3

hist(CrocMetrics$HRoverlap)

model3 <- lmer(HRoverlap ~ Class + (1|TRANSMITTERID), data = CrocMetrics)

check_model(model3)

car::Anova(model3, test.statistic = "F")
