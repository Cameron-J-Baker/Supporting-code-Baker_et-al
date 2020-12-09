# This function prepares the raw overlap matrices into the format required when using the LAR function from the asnipe
#  R package
LARprep <- function(lmig, years, class = c("Total","Resident","Nomadic","Female")){ #years is the list of unique monthly IDs from the crocs data. Will be used to add column to each matrix to assign when they occurred
  
  ####========================================================================================================####
  # lmig : The list of HR overlap matrices
  # years : A vector with all of the time points in the same order as the matrices within the list
  # class : Which movement tactics would you like to examine. Can do any combination of the three movement tactics
  
  
  # First for each of the matrices present, go through and remove all excess individuals and set overlaps to either 0 or 1
  # then convert each matrix to a data.frame for setting the time point for each individual.
  o <- 1
  for(i in 1 : length(lmig)){
    
    lmig[[i]][is.na(lmig[[i]])] <- 0 # replace all NA values with 0
    
    # Next, set the indviduals who are overlapping with each other as 1 to demonstrate they are apart of a group. If a
    # pair of individuals are 
    lmig[[i]][lmig[[i]] < 0.01] <- 0 # any value less than 0.05 will be change to 0
    lmig[[i]][lmig[[i]] >= 0.01] <- 1 # any value greater than or equal to 0.05 will be change to 0
    
    # Convert the matrix to a data.frame and then make row.names into a column named TRANSMITTERID
    lmig[[i]] <- as.data.frame(lmig[[i]]) # convert from a matrix to a data.frame
    lmig[[i]] <- lmig[[i]][rowSums(lmig[[i]][, -1])>0, ] # remove any rows that only have 0 present i.e., individuals not found that year
    
    if(nrow(lmig[[i]]) == 0){
      #lmig <- list.remove(lmig, i)
      #lmig[[i]] <- NULL # if after doing above that a month is found to have 0 overlaps, remove from the list
    } else { # else assign a column to include an individuals ID and one for what month this is
      lmig[[i]] <- rownames_to_column(lmig[[i]], var="TRANSMITTERID") # sets the rownames to a column
      lmig[[i]]$Month <- years[o]
      o <- o +1
    }
  }
  
  overlaps <- as.data.frame(Reduce(rbind, lmig)) # Merge each of the matrices in the list into a single dataframe
  overlaps$Month <- as.Date(paste(overlaps$Month, "01", sep = "-"))
  
  #----------------------------------------------------------------------------------------------------------------#
  # Create the time series file so that each individual as their own unique 0 point and time series for when they're present
  myids <- unique(overlaps$TRANSMITTERID)
  
  IDorder <- list() # create an empty list to place each individual into, essentially reordering the detections by individuals
  time_list <- list() # create an empty list to place each individuals time points into
  
  for(t in 1: length(myids)){
    
    if(length(class) == 1){
      if(class == "Total"){
        # IF "Total" is choosen, simply run the function as originally intended to retrieve the data for all
        # potential dyads
        IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
        # Create the time series data for the analysis
        time_between <- 0
        r <- 1
        while(r <= nrow(IDorder[[t]])){
          # create the time series for the individual from 0
          time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
          r <- r + 1
        }
        time_list[[t]] <- time_between
        IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
        IDorder[[t]]$Month <- NULL # remvoe the month column
        
      } else {
        # If we are not after the total, susbet out the group of interest
        IDs <- unique((subset(categories, Class == class))$TRANSMITTERID)
        if(is.element(myids[t], IDs) == T){
          #----------------------------------------------------------------------#
          # subset out the focal individuals data
          IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
          # Create the time series data for the analysis
          time_between <- 0
          r <- 1
          while(r <= nrow(IDorder[[t]])){
            # create the time series for the individual from 0
            time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
            r <- r + 1
          }
          time_list[[t]] <- time_between
          IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
          IDorder[[t]]$Month <- NULL # remvoe the month column
          
          for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
            if(is.element(names(IDorder[[t]][q]), IDs) == F){ # test to see if the individual for that column is in the list
              IDorder[[t]][,q] <- 0 # if not set their values to 0
            }
          }
        }
      }
    }  else {
      # Subset out the two lists of interest
      IDs1 <- unique((subset(categories, Class == class[1]))$TRANSMITTERID)
      IDs2 <- unique((subset(categories, Class == class[2]))$TRANSMITTERID)
      #----------------------------------------------------------------------#
      if(is.element(myids[t], IDs1) == T | is.element(myids[t], IDs2) == T){
        # subset out the focal individuals data
        IDorder[[t]] <- subset(overlaps, TRANSMITTERID == myids[t]) # subset out a specific individual and place the individual into the list in position t
        # Create the time series data for the analysis
        time_between <- 0
        r <- 1
        while(r <= nrow(IDorder[[t]])){
          # create the time series for the individual from 0
          time_between[[r]] <- as.numeric(difftime(IDorder[[t]]$Month[r],paste0(years[1],"-01")))
          r <- r + 1
        }
        time_list[[t]] <- time_between
        IDorder[[t]]$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
        IDorder[[t]]$Month <- NULL # remvoe the month column
        #----------------------------------------------------------------#
        # Determine which list the individual of interest belongs to
        # This is important as we want this individuals matrix to reflect how it spatially overlaps with only
        # individuals in the other strategy to it
        if(is.element(myids[t], IDs1) == F){
          # If the focal individual is not in the first list, add it's ID to the first list to create the list
          # of IDs for which to keep the overlaps (1) within the matrix
          IDs1[length(IDs1)+1] <- myids[t] # add the ID of interest to the previous column
          for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
            if(is.element(names(IDorder[[t]][q]), IDs1) == F){ # test to see if the individual for that column is in the ID list to keep
              IDorder[[t]][,q] <- 0 # if not set their values to 0
            }
          }
        } else {
          # If that individual is in that list, add it's ID to the other list and then filter the data
          IDs2[length(IDs2)+1] <- myids[t] # add the ID of interest to the previous column
          for(q in 1:length(names(IDorder[[t]]))){ # for all of the columns present
            if(is.element(names(IDorder[[t]][q]), IDs2) == F){ # test to see if the individual for that column is in the ID list to keep
              IDorder[[t]][,q] <- 0 # if not set their values to 0
            }
          }
        }
      }
    }
    
    
    #time_list[[t]] <- time_between
    
    #time_list[[t]] <- seq(0, by = 365, length.out = nrow(IDorder[[t]])) # create the time sequence that they are present
  }
  
  
  IDorder <- data.frame(Reduce(rbind, IDorder)) # convert the list to a data.table
  IDlist <- unique(IDorder$TRANSMITTERID) # create a list of the ID's in the order they're in the data frame
  #IDorder$TRANSMITTERID <- NULL # remove the TRANSMITTERID column
  #IDorder$Month <- NULL
  
  IDorder <- as.matrix(IDorder) # then convert ID order into a matrix to be used in the LAR function
  
  # Combine the list of vectors for the time data
  time_list <- Reduce(c,time_list)
  
  # combine the times vector and the overlap matrix into a list to be output from the function
  output <- list(IDorder, time_list, IDlist)
  
  return(output) 
  # Add a new string of vectors to the times list to indicate this time point
  #timepoint <- nrow(lmig[[i]])
  #time_list[[i]] <- rep(time_vec[i], timepoint)
}
