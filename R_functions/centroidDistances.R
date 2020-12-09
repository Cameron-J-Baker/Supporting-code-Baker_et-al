# Function to calculate the distance an individual is from it's home range centroid
# currently set to run with already UTM transformed data
centroidDist <- function(sdata, cent) # set the values it is searching for
{
  # sdata : A data frame outputted from the COAV2 and SnapCOAtoRiver functions
  # cent : The lat long points that you want to measure from
  
  xy <- sdata
  coordinates(xy) <- c("Longitude.coa.snap", "Latitude.coa.snap") # convert the data to a spatial points format
  proj4string(xy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #give it a proj4string to WGS84
  xy <- spTransform(xy,UTM) # project to UTM
  
  SP <- SpatialPoints(xy,proj4string=UTM) # retrive only the spatial points
  
  # Initialise variables
  iCount <- 1
  rDISTANCE <- rep(0,length(xy)-1)
  sDISTANCE <- rep(0,length(xy)-1)
  
  # For each line in track file
  while (iCount <= nrow(xy))
    #while (iCount <= 20)
    
  {
    SP1 <- cent
    SP2 <- SP[iCount,]
    #plot(SP1,add=T,col=3)      
    #plot(SP2,add=T,col=4)
    
    crowdist <- spDists(SP1,SP2)
    sPath2 <- shortestPath(tr, SP1, SP2, output="SpatialLines")
    riverdist <- SpatialLinesLengths(sPath2)
    
    sDISTANCE[iCount] <- crowdist      
    rDISTANCE[iCount] <- riverdist 
    
    iCount <- iCount + 1
  }
  sdata$rDISTANCE <- rDISTANCE/1000 #divide by 1000 to convert from m to km
  sdata$sDISTANCE <- sDISTANCE/1000 #divide by 1000 to convert from m to km
  return(sdata)
}
