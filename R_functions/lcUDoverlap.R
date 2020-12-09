#The below function determines the proportion of overlap present between least coast utilisation distributions. The 
# function will eventually be incorporated in the VTrack R package

lcUDoverlap <- function(lcUDs, percent = 95, method = c("VI", "UDOI", "AO"), UTM =CRS("+init=epsg:32754")){
  ###############################################
  # lcUDs : A list object containing the outputs for multiple individuals from the lcdistance package
  # percent : The level of HR overlap that you would like to examine. Default is 95%
  # method : Define the method to be used to determine the proportion of HR overlap present. Currently only the Volumne
  #          of intersection (VI) and utilisation distribution overlap index (UDOI). The absolute overlap (AO) is not functional
  # UTM : UTM epsg code of the study area. Deafult is currently for the Wenlock and Ducie river systems
  ##########################################################
  # Extract the UD.rasters from each of the outputs and standardise the observed values to between 0 and 1
  UDlist <- list()
  
  for(i in 1:length(lcUDs)){
    # Extract the sum of all values present in the raster
    UDras <- lcUDs[[i]]$UD.raster
    UDras[values(UDras) > percent] <- NA # select the 95% UD
    values(UDras) <- 100 - values(UDras) # invert the UD values so that areas used more frequently have higher values than at the edge of the HR
    UDras[is.na(values(UDras))] <- 0 # set all NA's to 0
    values(UDras) <- values(UDras)/sum(values(UDras), na.rm = T) # standardise the values to sum to 1
    
    UDlist[[i]] <- UDras # add each raster progressively into a list
  }
  #names(UDlist) <- names(lcUDs) # reassign the name of each raster
  
  res <- matrix(0, ncol = length(lcUDs), nrow = length(lcUDs)) # create the matrix to save the results into
  
  for(i in 1:length(lcUDs)){
    for(j in 1:i){
      
      if(method == "VI"){
        vi <- values(UDlist[[i]]) # retrieve the values of the first individual
        vj <- values(UDlist[[j]]) # retrieve the values of the second individual
        res[i, j] <- res[j, i] <- sum(pmin(vi,vj)) 
      }
      
      if(method == "UDOI"){
        vi <- values(UDlist[[i]]) # retrieve the values of the first individual
        vj <- values(UDlist[[j]]) # retrieve the values of the second individual
        # determine the number of cells that the individuals overlap together
        ai <- vi # create a copy of vj
        aj <- vj # create a copy of vj
        ai[ai != 0] <- 1 #set all values that are not 0 to 1
        aj[aj != 0] <- 1
        ak <- sum(ai*aj) # retrieves the number of cells that individuals are overlapping
        # assign the values
        res[i, j] <- res[j, i] <- ak * sum(vi * vj) 
      }
      
      if(method == "AO"){
        vi <- rasterToPolygons(UDlist[[i]], dissolve = T)
        vi <- spTransform(vi,UTM) #convert it to UTM
        vj <- rasterToPolygons(UDlist[[j]], dissolve = T)
        vj <- spTransform(vj,UTM) #convert it to UTM
        
        # Compute the area of overlap
        res[i, j] <- res[j, i] <- gArea(gIntersection(vi,vj))/gArea(gUnion(vi,vj))
      }
      
    }
  }
  rownames(res) <- names(lcUDs) # add the ID labels to the rows
  colnames(res) <- names(lcUDs) # add the ID labels to the columns
  return(res)
}
