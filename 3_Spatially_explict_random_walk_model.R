# Author: Cameron Baker
# 03/02/2020
# cameron.baker@uqconnect.edu.au
# PhD Capter 1: Spatial structure of the estuarine crocodile
# Aim: Below is the code for the spatially explict random model that was used to simulate the movements of individuals 
# to determine if the observed HR overlaps between conspecifics were different from chance. 

# Load the required packages
library(data.table)
library(sp)
library(rgeos)
library(rgdal)
library(gdistance)
library(foreach)
library(doParallel)

#########################################################################################################################
# Import the required data

# First import the COA of the population, as this is what will be simulated by the model
coa_crocs <- fread("Outputs/Crocodile_COA.csv")

# set the DATETIME column to be a posixct object
coa_crocs$DATETIME <- as.POSIXct(coa_crocs$DATETIME, format="%Y-%m-%d %H:%M:%S", tz = "Australia/Brisbane")

# create a list of the IDs for each individual
myid <- unique(coa_crocs$TRANSMITTERID)

# Define UTM to project data into for calculation of home ranges
UTM <- CRS("+init=epsg:32754") #a value to convert data to the UTM projection

# The location of each of the acoustic hydrophones
hydro_loc <- fread("Data/Acoustic_hydrophone_locations.csv")
#--------------------------------------------------------------------------------------------------------------------#
# Prepare the river lines data for the COA snap to lines below
#--------------------------------------------------------------------------------------------------------------------#
# firstly clip the Cape york lines file to the Wenlock river
Wen_Riverclip <- readOGR("Data/Line_file","Wenlock_Ducie_linefile") # Load the Wenlock spatiallines  data
Wen_Riverclip <- spTransform(Wen_Riverclip,UTM) # project to UTM
#plot(Wen_Riverclip)

# Import the raster of the Wenlock river to calculate distances along
Wen_Raster <- raster("Data/WenlockRaster.asc")
#plot(Wen_Raster)
# Convert the raster to a transition object - required for the river distance algorythm
tr <- transition(Wen_Raster, transitionFunction=mean, directions=8) # Create a Transition object from the raster

#######################################################################################################################
# Import the required functions #######################################################################################
source("R_functions/centroidDistances.R")
source("R_functions/consecutiveDistances.R")

#######################################################################################################################
# Next, for each inidivdual calculate their mean and SD distance travelled between successive detections, along with
# the probability that they either go upstream, downstream or remain in the same position

moveSummary <- function(i, COAs, tr){
  
  #---------------------------------------------------#
  # i : which position along the myid vector to choose
  # COAs : The output of the COA and SnapCOAtoRiver functions
  # tr : The transition object of the study area
  
  # first subset down to a specific individual
  indiv <- COAs[TRANSMITTERID == myid[i]]
  indiv$DATETIME <- as.POSIXct(indiv$DATETIME)
  indiv <- indiv[Longitude.coa.snap >= 141.8688] 
  indiv <- indiv[Latitude.coa.snap <= -11.92442]
  # Calculate the consecutive distance between each COA
  cosecdist <- consecDistTravelled(indiv, tr)
  
  #Create a dataframe with the longitude and latitude data for a polygon around most downstream vr2W (PM Cullen AR)
  AMTD0pos <- data.frame(LONGITUDE = c(141.91287), 
                         LATITUDE = c(-11.95887))
  coordinates(AMTD0pos) <- ~ LONGITUDE + LATITUDE #convert it to a spatial points dataframe
  proj4string(AMTD0pos) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #give it a proj4string to WGS84
  AMTD0posUTM <- spTransform(AMTD0pos,UTM) #project to UTM
  
  distmouth <- centroidDist(indiv, cent = AMTD0posUTM)
  
  # generate a new column in distmouth which determines if an individual has gone either up or down stream. Up will be indicated
  # by a positive value, down by a negative value. 
  t <- 1
  for(p in 1:nrow(distmouth)){
    distmouth$Direction[t] <- ifelse(t == 1, 0, distmouth[t,10] - distmouth[t-1,10]) 
    t <- t + 1
  }
  
  up <- subset(distmouth, Direction > 0)
  down <- subset(distmouth, Direction < 0)
  stay <- subset(distmouth, Direction == 0)
  
  output <- data.frame(TRANSMITTERID = myid[i],
                       meanStep = mean(cosecdist$rDISTANCE),
                       sdStep = sd(cosecdist$rDISTANCE),
                       propUp = nrow(up)/nrow(indiv),
                       propDown = nrow(down)/nrow(indiv),
                       propStay = nrow(stay)/nrow(indiv),
                       StartDate =format(min(indiv$DATETIME), "%Y-%m-%d %H:%M:%S"),
                       EndDate = format(max(indiv$DATETIME), "%Y-%m-%d %H:%M:%S"),
                       sLongitude = indiv$Longitude.coa.snap[1], # retrieve the first long and lat observed
                       sLatitude = indiv$Latitude.coa.snap[1])
  
  return(output)
}

#-----------------------------------------------------------------------------------------------------------#
# Run the function in paraellel for all individuals and save the output

registerDoParallel(10) # set the number of cores to use

lmig <- foreach(i = 1:(length(myid)), 
                .packages = c("data.table","gdistance","plyr", "sp")) %dopar%
  moveSummary(i, coa_crocs, tr)

stopImplicitCluster() # stop the cluster to return memory and resources to other systems

CrocMetrics <- rbindlist(lmig)

write.csv(CrocMetrics,"Data/Rand_walk_parameters_lcUD.csv", row.names = F) 

#######################################################################################################################
# Create the two functions that form the spatially explict random walk model ##########################################

# This function generates the random movements points as part of the randomwalk model
nextPoint <- function(X = X, direction, steps, indiv_stretch, AMTD0posUTM = AMTD0posUTM){
  for(t in 1:length(direction)){
    
    #print(paste("Simulation", p))
    if(direction[t] == 2){ 
      X[t + 1,] <- X[t,]
      #p <- p + 1
      #OBS[t,] <- 0 # save a vector with the step lengths selected to allow debugging as required
      
    } else {
      
      # Extract the previous point
      previous <- X[t,]
      coordinates(previous) <- c("x", "y") # convert the data to a spatial points format
      proj4string(previous) <-"+init=epsg:32754" # give it a proj4sting to UTM
      # Take a random step length from the generated normal distibution of step lengths
      # Repeat this step until the buffer surrounding a point intersects the river at least once. This is to prevent
      # the model getting stuck within small subsets of the river that are present when the area used is subset out
      q <- 1
      repeat{
        # This part takes a random step length from the distribution
        # If the step length choosen though is not enough to overlap with a part of the river, it will take the sample
        # again from a restricted subset of step lengths greater than the previous one. This is to stop individuals 
        # becomming stuck on small subsets of creeks that overlap their HR area
        if(q == 1){ # if this is the first repeat run as normal randomly selecting a point
          n.step <- n.step2 <- sample(steps, 1) # take a random step length sample while saving a copy encase it's needed
        } else { # if an individuals last step attempt did not overlap anywhere, do it again this time with a larger step to prevent the individual being stuck
          lessstep <- subset(steps, steps > n.step2 ) # subset steps to only ones bigger than the previous
          
          if(length(lessstep) != 0){
            n.step <- n.step2 <- sample(lessstep, 1)
          } else {
            # If an individuals last step attempt was it's maximum, and it has overshot it's HR area, resample 
            # from steps and start again. Will only occur for weird low moving animals like the Trans individuals
            n.step <- n.step2 <- sample(steps, 1)
          }
          
        }
        
        #OBS[t,] <- n.step # save a vector with the step lengths selected to allow debugging as required
        # using this length generate a buffer surrounding the previously generated point
        buffered.point <- gBuffer(previous, width = n.step) # create a buffer surrounding the previous point
        proj4string(buffered.point) <- "+init=epsg:32754" #give it a proj4string to WGS84
        buffered.point <- spTransform(buffered.point,UTM) # project to UTM
        buffered.point <- as(buffered.point, "SpatialLines") # convert to a spatial lines file to extract the overlapping points
        #plot(buffered.point, add = T)
        
        # using the intersect function, extract the points that the buffer intersects the river.
        inter <- gIntersection(indiv_stretch, buffered.point, byid = T) # Obtain all of the points that the buffered point intersects the line
        q <- q + 1
        if(length(inter) != 0){
          break
        }
      }
      
      inter <- as.data.frame(coordinates(inter)) # create a data.frame with the intersect points present
      
      inter <- rbind(inter, as.data.frame(previous)) # At the bottom add in the previous point
      ####------------------------------------------------------------------------------------------------------#####
      # Calculate the distance that each point is from the mouth of the Wenlock river
      ####------------------------------------------------------------------------------------------------------#####
      # Create a duplicate copy and transform it into a spatial points object to calculate the distance from river mouth
      inter2 <- inter
      coordinates(inter2) <- c("x", "y") # convert the data to a spatial points format
      proj4string(inter2) <-"+init=epsg:32754" # give it a proj4sting to UTM
      SP <- SpatialPoints(inter2,proj4string=UTM) # retrive only the spatial points
      
      # create empty vector to store the results into
      rDISTANCE <- rep(0,nrow(inter)-1)
      mouthDISTANCE <- rep(0,nrow(inter)-1)
      crowDirection <- rep(0,nrow(inter)-1)
      iCount <- 1
      # For each line in track file
      while (iCount <= nrow(inter))
        #while (iCount <= 20)
        
      {
        SP1 <- AMTD0posUTM
        SP2 <- inter2[iCount,]
        SP3 <- inter2[length(inter2),]     
        #plot(SP2,add=T,col=4)
        # Calculate the distance from the mouth of the river to determine the direction of travel (up or downstream)
        sPath2 <- shortestPath(tr, SP1, SP2, output="SpatialLines") 
        mouthdist <- SpatialLinesLengths(sPath2)
        mouthDISTANCE[iCount] <- mouthdist 
        # Also calculate the straightline distance from the mouth encase an individual ends up a side creek that runs
        # longer than the river, with the next point being selected onto the river
        crowDirection[iCount] <- spDists(SP1,SP2)
        
        # Calculate the distance from the last detection to choose the correct position
        sPath3 <- shortestPath(tr, SP3, SP2, output="SpatialLines") 
        riverdist <- SpatialLinesLengths(sPath3)
        
        rDISTANCE[iCount] <- riverdist
        
        
        iCount <- iCount + 1
      }
      
      inter$Direction <- mouthDISTANCE - mouthDISTANCE[length(mouthDISTANCE)] # Assigns the direction (either up or downstream) of each point from the previous
      inter$crowDirection <- crowDirection - crowDirection[length(crowDirection)] # assigns the direction travelled based on the straight line distance from the mouth
      # this is a backup to identify when individuals are going up or down stream encase an individual only goes a single direction from the river distances (i.e. up a side creek that is longer than the corresponding points on a river it may end up at)
      inter$Distance <- rDISTANCE # the distance each point is from the first
      
      ####------------------------------------------------------------------------------------------------------#####
      # Assign next point of the sequence
      ####------------------------------------------------------------------------------------------------------#####
      if(direction[t] == 0){ # if direction at time t is equal to 0 than an individual is going downstream
        
        inter_down <- subset(inter, Direction < 0) #subset out all points less than 0 indicating downstream
        
        if(nrow(inter_down) == 0) { # if there are no points downstream
          # If there are no points downstream using the river distance, examine again using crow distance
          inter_down <- subset(inter, crowDirection < 0) #subset out all points less than 0 indicating downstream using crowdistance
        }
        #---------------------------------------------#
        if(nrow(inter_down) == 1) { # If there is only one intersecting point assign that one to the next position
          X[t + 1,] <- inter_down[1,1:2] # assign the point
          # p <- p + 1
          #coordinates(inter_down) <- ~ x + y #convert it to a spatial points dataframe
          #plot(inter_down, col = "red", add = T) # plot the next point
          #---------------------------------------------#
        } else if(nrow(inter_down ) > 1){ #if there are two or more though, determine which one is closest to the step length choosen and assign it
          # Calculate the difference of the remaining points from the selected step length
          inter_down$select <- inter_down$Distance - n.step
          # Subset out the point which has the minimum difference between the distance travelled and step size selected
          inter_down <- data.table::setDT(inter_down)[, .SD[which.min(select)]] 
          
          X[t + 1,] <- inter_down[1,1:2] # assign the point
          #  p <- p + 1
          #coordinates(inter_down) <- ~ x + y #convert it to a spatial points dataframe
          #plot(inter_down, col = "red", add = T) # plot the next point
          #---------------------------------------------#
        } else {
          # if it is still 0 then the individual has reached the end of its HR and is to go in the opposite direction
          inter_up <- subset(inter, Direction > 0) # subset out all points upstream of the last
          
          if (nrow(inter_up) == 0){ # if this is less than 0 some how resubset it using the crowDistacne
            inter_up <- subset(inter, crowDirection > 0) #subset out all points greater than 0 indicating downstream using crowdistance
          }
          
          if(nrow(inter_up) == 1) { # If there is only one intersecting point assign that one to the next position
            X[t + 1,] <- inter_up[1,1:2] # assign the point
            #   p <- p + 1
            #coordinates(inter_up) <- ~ x + y #convert it to a spatial points dataframe
            #plot(inter_up, col = "orange", add = T) # plot the next point
            
          } else { #if there are two or more though, determine which one is closest to the step length choosen and assign it
            # Calculate the difference of the remaining points from the selected step length
            inter_up$select <- inter_up$Distance - n.step
            # Subset out the point which has the minimum difference between the distance travelled and step size selected
            inter_up <- data.table::setDT(inter_up)[, .SD[which.min(select)]] 
            
            X[t + 1,] <- inter_up[1,1:2] # assign the point
            #   p <- p + 1
            #coordinates(inter_up) <- ~ x + y #convert it to a spatial points dataframe
            #plot(inter_up, col = "orange", add = T) # plot the next point
          }
          
        }
        
      } else { # the individual moves upstream
        inter_up <- subset(inter, Direction > 0) #subset out all points greater than 0 indicating upstream
        
        if(nrow(inter_up) == 0) { # if there are no points downstream
          # If there are no points downstream using the river distance, examine again using crow distance
          inter_up <- subset(inter, crowDirection > 0) #subset out all points greater than 0 indicating upstream using crowdistance
        }
        #---------------------------------------------#
        if(nrow(inter_up) == 1) { # If there is only one intersecting point assign that one to the next position
          X[t + 1,] <- inter_up[1,1:2] # assign the point
          # p <- p + 1
          #coordinates(inter_up) <- ~ x + y #convert it to a spatial points dataframe
          #plot(inter_up, col = "red", add = T) # plot the next point
          #---------------------------------------------#
        } else if (nrow(inter_up) > 1){ #if there are two or more though, determine which one is closest to the step length choosen and assign it
          # Calculate the difference of the remaining points from the selected step length
          inter_up$select <- inter_up$Distance - n.step
          # Subset out the point which has the minimum difference between the distance travelled and step size selected
          inter_up <- data.table::setDT(inter_up)[, .SD[which.min(select)]] 
          
          X[t + 1,] <- inter_up[1,1:2] # assign the point
          # p <- p + 1
          #coordinates(inter_up) <- ~ x + y #convert it to a spatial points dataframe
          #plot(inter_up, col = "red", add = T) # plot the next point
          #---------------------------------------------#
        } else { # if there are no points upstream of the current one, then the indiv has reached that end of their HR and need to go in the opposite direction
          
          inter_down <- subset(inter, Direction < 0) #subset out all points less than 0 indicating downstream
          
          if (nrow(inter_down) == 0){ # if this is less than 0 some how resubset it using the crowDistacne
            inter_down <- subset(inter, crowDirection < 0) #subset out all points greater than 0 indicating downstream using crowdistance
          }
          
          if(nrow(inter_down) == 1) { # If there is only one intersecting point assign that one to the next position
            X[t + 1,] <- inter_down[1,1:2] # assign the point
            # p <- p + 1
            #coordinates(inter_down) <- ~ x + y #convert it to a spatial points dataframe
            #plot(inter_down, col = "orange", add = T) # plot the next point
            
          } else { #if there are two or more though, determine which one is closest to the step length choosen and assign it
            # Calculate the difference of the remaining points from the selected step length
            inter_down$select <- inter_down$Distance - n.step
            # Subset out the point which has the minimum difference between the distance travelled and step size selected
            inter_down <- data.table::setDT(inter_down)[, .SD[which.min(select)]] 
            
            X[t + 1,] <- inter_down[1,1:2] # assign the point
            # p <- p + 1
            #coordinates(inter_down) <- ~ x + y #convert it to a spatial points dataframe
            #plot(inter_down, col = "orange", add = T) # plot the next point
          } 
        }
      }
    }
  }
  return(X)
}

# The spatially explict random walk function.
randWalk_COA <- function(i, observed, moveMet, linefile,   
                         simulations = 1){
  
  ###-----------------------------------------------------------------------------##
  # i : which position along the myid vector to choose which individual to simulate
  # observed : The observed COA data
  # moveMet : The output of the moveSummary function containing the movement parameters of each individual
  # linefile : A line file of the river system
  # simulations : The number of simulations that you would like to run
  
  #####=========================================================================================================######
  # Initial data preparation and model parameterization
  #####=========================================================================================================######
  # First subset the observed data to the specific individual of interest
  observed <- as.data.table(observed)
  HRparameters <- as.data.table(HRparameters)
  moveMet <- as.data.table(moveMet)
  observed <- observed[TRANSMITTERID == myid[i]]
  observed$TRANSMITTERID <- factor(observed$TRANSMITTERID)
  observed2 <- observed # create a copy of the data to project in order to generate an individuals HR
  HRpara <- HRparameters[TRANSMITTERID == myid[i]] # Subset out the individuals HR parameters
  moveMet <- moveMet[TRANSMITTERID == myid[i]]
  
  # Create a UTM parameter to convert the observations to to calculate distances
  UTM <- CRS("+init=epsg:32754")
  
  # if the simulation is to be unconstrained, the individual is free to move within the full extent of the study area
  indiv_stretch <- linefile
  
  # If a full simulation is not required, only simulate the number of observed detections an individual was tracked
  n <- nrow(observed) # create an n vector to get the length of a time simulations must be completed for
  #---------------------------------------------------------------------------------------------------------#
  
  # Next create a dataframe with the longitude and latitude data for a point at the most downstream vr2W (PM Cullen AR)
  AMTD0pos <- data.frame(LONGITUDE = c(141.91287), 
                         LATITUDE = c(-11.95887))
  coordinates(AMTD0pos) <- ~ LONGITUDE + LATITUDE #convert it to a spatial points dataframe
  proj4string(AMTD0pos) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #give it a proj4string to WGS84
  AMTD0posUTM <- spTransform(AMTD0pos,UTM) #project to UTM
  #-----------------------------------------------------------------------------------------------------------#
  # Paramiterise the step length and direction that the individual will move, along with the vectors in which to save the movements to
  A <- 1
  # Run the random walk function for specified number of simulations required
  for(q in 1:simulations){
    print(paste(myid[i],"Simulation", A))
    # Generate the step length distribution - multiple bu 1000 to convert from km to m
    steps <- rnorm(1000, moveMet$meanStep, moveMet$sdStep) * 1000 # creates a random normal distribution based on the mean and SD of step lengths
    steps <- subset(steps, steps >=0 ) # as an individual can not have a negative step length, remove any values less than 0
    
    # Generate the movement sequence up and down the river based on an individuals probability of moving each direction
    # 0 = downstream, 1 = upstream, 2 = remains at same position
    direction <- sample(c("0","1","2"), size = n-1, replace = T, # n - 1 as the first position will be randomly selected
                        prob = c(moveMet$propUp, moveMet$propDown, moveMet$propStay))
    
    # Create the vectors in which to store the movements and directions of each step
    X <- data.frame(matrix(NA, nrow = n, ncol = 2)) # dataframe to store the simulated lat longs
    names(X) <- c("x", "y")
    OBS <- matrix(NA,n,2) # matrix for the steps lengths and directions
    
    #####=========================================================================================================######
    #Generate the starting point of the sequence
    #####=========================================================================================================######
    start <- data.frame(x = moveMet$sLongitude, y = moveMet$sLatitude)
    coordinates(start) <- ~ x + y #convert it to a spatial points dataframe
    proj4string(start) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #give it a proj4string to WGS84
    start <- spTransform(start, UTM) # project to UTM
    #plot(start, col = "blue", add = T) # add the starting point to the plot as a red dot
    
    X[1,] <- coordinates(start) # Add the coordinates of the starting point to the simulated lat long data.frame
    
    #####=========================================================================================================######
    # Generate the rest of the sequence
    #####=========================================================================================================######
    
    X <- nextPoint(X = X, direction, steps, indiv_stretch, AMTD0posUTM)
    
    # If only an observed simulation has been conducted paste the simulated points to the observed dates and export
    sim.seq <- data.frame(DATETIME = observed$DATETIME,
                          x = X[,1],
                          y = X[,2])
    # Export the simulation
    write.csv(sim.seq, paste("Outputs/Simulated/", observed$TRANSMITTERID[1], "-simulation-", A, ".csv", sep = ""), row.names = F)
    A <- A + 1 
  } 
  return(sim.seq)
}

# Example use
simulated <- list() # create an empty list
for(i in 1:length(myid)){
  simulated[[i]] <- randWalk_COA(i, coa_crocs, moveMet = CrocMetrics, linefile =  Wen_Riverclip, 
                                 simulations = 1)
}

simulated <- rbindlist(simulated)

# The output from the randwalk_COA function can then be utilised instead of the observed COAs in the 
# 1_Determine_monthly_HR_overlap R script to gernerate the random HR overlaps. 

# For the manuscript this process was performed using the high performance computing facilaties available at the 
# University of Queensland. Specifically, with the assisatance of Dr. David Green, the director of the Research Computing
# Center (https://rcc.uq.edu.au/), the above code was run on the Tinaroo HPC cluster using an embedded Nimrod project,
# which is an embedded RCC specific workflow for parametric compution. 