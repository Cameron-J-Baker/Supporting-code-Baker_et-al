# This function prepares the data into a format suitable for performing the repeatability analysis
rptPrep <- function(lmig, years, dyad = T){
  
  #---------------------------------------------------------------------------#
  # Lmig : A list of overlap matrcies
  # years : a vector with the month and year of each of the corresponding enteries in the list
  # dyad : if T, determines dyad specific values
  
  o <- 1 # create a value to add the correct years to a matrix
  for(i in 1:length(lmig)){
    # set the upper half and diagonal of the matrix to NA
    lmig[[i]][upper.tri(lmig[[i]], diag = T)] <- NA 
    # This code converts the matrix to a data.frame with only the lower triangle of the matrix excluding the diagonal
    lmig[[i]] <- na.omit(as.data.table(as.table(lmig[[i]], na="", row.names=T, col.names=T))) # converts from matrix to table
    names(lmig[[i]]) <- c("ID1", "ID2","HRoverlap") # Rename the columns 
    lmig[[i]]$Month <- years[o] # add a column for the month this was observed occurring in
    lmig[[i]]$HRpresence = ifelse(lmig[[i]]$HRoverlap > 0.001, 1, 0)
    o <- o +1
  }
  
  overlaps <- as.data.table(Reduce(rbind, lmig)) # Merge each of the matrices in the list into a single dataframe
  
  overlaps$ID <- paste(overlaps$ID1, overlaps$ID2, sep = "-")
  
  if(dyad == F){
    indivMonth <- list()
    for(q in 1:length(years)){
      month <- overlaps[Month == years[q]] # subset down to a specific month of interst
      IDv <- unique(c(month$ID1, month$ID2))
      
      indivO <- list()
      for(t in 1:length(ID)){
        indiv <- month[ID1 == IDv[t]| ID2 == IDv[t]]
        
        indivO[[t]] <- data.table(ID = IDv[t],
                                  HRoverlap = mean(indiv$HRoverlap),
                                  No.Indivs = nrow(indiv),
                                  Month = years[q]) 
        
      }
      indivMonth[[q]] <- as.data.table(Reduce(rbind, indivO))
    }
    
    overlaps <- as.data.table(Reduce(rbind, indivMonth))
  }
  
  # set any overlaps less than 0.01 to 0
  overlaps$HRoverlap <- ifelse(overlaps$HRoverlap < 0.01, 0, overlaps$HRoverlap)
  return(overlaps)
}