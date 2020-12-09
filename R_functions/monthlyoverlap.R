#---------------------------------------------------------------------------------------------------------------------#
# This function calculates the lcUD for each individual present within a specific month and then calculates the
# amount of overlap between each individual. While doing this, it saves a copy of the lcUD and the resulting
# HR overlap matrix to drive to save having to run this function multiple times.
# PRIOR TO USE EITHER THE DESTINATION THAT THE FUNCTION IS SAVING THE OUTPUTS TO SHOULD BE CHANGED TO SUIT YOUR OWN 
# FILE DESITINATIONS OR DISABLED BY EITHER DELETING OR # OUT

monthlyoverlap <- function(i, coa_crocs, loc, cost.raster = cost.raster,  h = 200){
  
  ##########################################################################################
  # i : The point along the years vector to determine which month to perform the function on
  # coa_crocs : COA object produced from running the COA or COAV2 functions followed by the SnapCOAtoRiver function
  # loc : A dataframe containing the name and location (lat/long) of acoustic hydrophone stations
  # cost.raster : The cost raster of the study system
  # h : The smoothing parameter for the least cost utilisation distribution
  ##########################################################################################
  
  # First subset down to the specific year that will be examined
  crocyear <- coa_crocs[Year == years[i]]
  
  # create a list of the unique crocodile IDs present
  myid <- unique(crocyear$TRANSMITTERID)
  
  # Remove any individuals that have less than 10 centers of activity for a specific month or if they have only
  # been detected at 
  
  for(j in 1:length(myid)){ 
    
    indiv <- subset(crocyear, TRANSMITTERID == myid[j])
    # Remove any detections that are outside of the river system
    indiv <- indiv[Longitude.coa.snap >= 141.8688] 
    indiv <- indiv[Latitude.coa.snap <= -11.92442]
    if (nrow(indiv) < 10){
      crocyear <- subset(crocyear, TRANSMITTERID != myid[j])
    } else if (length(unique(indiv$Latitude.coa.snap)) == 1 & length(unique(indiv$Longitude.coa.snap)) == 1){
      crocyear <- subset(crocyear, TRANSMITTERID != myid[j])
    } else if (length(unique(indiv$Latitude.coa.snap)) <= 3 | length(unique(round(indiv$Longitude.coa.snap,4))) <= 3) { # if an indiv has less than 3 unique positions remove it
      crocyear <- subset(crocyear, TRANSMITTERID != myid[j])
    }
  }
  
  if(nrow(crocyear) > 0 & length(unique(crocyear$TRANSMITTERID)) > 1){ # if a month has no data skip it
    
    if(length(unique(crocyear$Latitude.coa)) != 1 & length(unique(crocyear$Longitude.coa)) != 1){
      # Convert COA_crocs$TRASMITTERID to a numeric then back to a factor to reset the number of individuals present
      crocyear$TRANSMITTERID <- as.integer(as.character(crocyear$TRANSMITTERID))
      
      myid2 <- unique(crocyear$TRANSMITTERID)
      
      # Create a list to save the HR estimates into
      lcUDs <- list()
      
      # For each individual left, run the lcdistance function to determine their lcUDs and save the results
      # into the empty list
      for(e in 1:length(myid2)){
        # Subset down to a specific individual to examine
        indiv <- crocyear[TRANSMITTERID == myid2[e]]
        
        lcUDs[[e]] <- lcDistance(tag_data = indiv,  ## station information, tagdata and taginfo data all in one ATTdata object
                                 receiver_loc = loc,
                                 cost = cost.raster, ## cost raster for your study site. if NULL it finds coastline data from OSM server
                                 cost2 = cost.ras,
                                 tr = trCost,
                                 ll_epsg = 4326,     ## EPSG code for the raw data (in lat/long)
                                 utm_epsg = 32754,    ## EPSG code for the Projected CRS for your study site (in meters)
                                 timestep = 360,      ## timestep in minutes for COA estimation (see COA() function for details of timestep)
                                 h = h,            ## smoothing parameter for UD estimate
                                 UDgrid = 50,
                                 cost.res = 50,
                                 directions = 16) 
        
      }
      # Add the names of each individual to their respective point on the list
      names(lcUDs) <- myid2
      
      # Write a copy of the lcUDs list to file to allow the home ranges of all individuals to be called in easily
      # to speed up future analyses
      write_rds(lcUDs, paste0("Output/lcUD-files/",years[i],".rds"))
      #-------------------------------------------------------------------#
      # Run the overlap analysis
      HRoverlap <- lcUDoverlap(lcUDs, percent = 95, method = "VI", UTM)
      
      #names(dimnames(HRoverlap)) <- c("ID", "")
      
      HRoverlap2 <- as.data.frame(HRoverlap) # convert from a matrix to a data.frame
      HRoverlap2 <- rownames_to_column(HRoverlap2, var="ID") # sets the rownames to a column
      
      # Combine both of the data.frames above while summing the Value columns together
      All_crocs <- data.frame( ID = unique(coa_crocs$TRANSMITTERID),
                               Value = 1)
      
      Year_crocs <- data.frame( ID = myid2,
                                value = 1)
      
      croc_list <- plyr::join(All_crocs, Year_crocs, by = "ID") # join the two lists together
      
      croc_list <- subset(croc_list, is.na(croc_list$value)) # subset the data so only the rows with NAs in the second list remain
      
      # For all of the crocodiles not present/HRs could not be determined for this year add a row and column of NA values
      for(h in 1 : nrow(croc_list)){
        
        # Add a new row of NA values for this individual
        HRoverlap2[nrow(HRoverlap2)+1, 1] <- croc_list[h, 1] # can use this to create a new column with the ID of each croc not included and then everything else NA
        
        HRoverlap2$new <- NA # create a new column
        
        colnames(HRoverlap2)[colnames(HRoverlap2)=="new"] <- croc_list[h, 1] # rename the new column to the individuals name
      }
      
      # Reorder the HR overlap data.frame so that the rows are in ascending order
      
      HRoverlap3 <- HRoverlap2 %>% arrange(ID) # rearrange the rows so that they are in ascending order
      
      HRoverlap3 <- HRoverlap3[ , order(names(HRoverlap3))] # place the columns into ascending order
      
      HRoverlap3 <- column_to_rownames(HRoverlap3, var = "ID") # convert ID back into row.names
      
      HRoverlap3 <- data.matrix(HRoverlap3, rownames.force = "ID") # convert back into a matrix
      
      write.csv(HRoverlap3, file=paste0("Output/Overlap_Matrices/",years[i], '.csv')) # wrtie a copy to be used later for other anylses
      
      return(HRoverlap3)
    }
  }
}

# The output of this function is a symmetrical matrix of the HR overlap between conspecifics for a specific month
