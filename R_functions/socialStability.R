# The function below determines whether individuals are maintaining overlaps with the same conspecifics for each month


socialStability <- function(overlaps, categories = croc_categories,
                            Class = c("ResRes", "NomNom", "ResFem", "NomFem", "FemFem", "ResNom", "Overall")){
  ############################################################################################
  # The inputs for this function are:
  # overlaps : A dataframe with the monthly HR proportion of each dyad
  # categories : A dataframe with the ID and movement tactic of each individual
  # Class : The movememnt tactic combination that you would like to examine
  
  #============================================================================#
  # First step; subset the data down to the group of interest and create a vector with the IDs to examine
  if(Class == "ResRes"){
    overlaps1 <- overlaps[Class == "ResidentResident"] # extract just the low movement males to all others
    IDs <- unique((subset(categories, Class == "Resident"))$TRANSMITTERID) # create a list of the low individuals
    
  } else if(Class == "NomNom"){
    overlaps1 <- overlaps[Class == "NomadicNomadic"] # extract just the high movement males to all others
    IDs <- unique((subset(categories, Class == "Nomadic"))$TRANSMITTERID) # create a list of the high individuals
    
  } else if(Class == "ResFem"){
    overlaps1 <- overlaps[Class == "ResidentFemale"]
    IDs <- unique((subset(categories, Class == "Low"))$TRANSMITTERID) # create a list of the low individuals
    
  } else if(Class == "NomFem"){
    overlaps1 <- overlaps[Class == "NomadicFemale"]
    IDs <- unique((subset(categories, Class == "High"))$TRANSMITTERID) # create a list of the high individuals
    
  } else if(Class == "ResNom"){
    overlaps1 <- overlaps[Class == "ResidentNomadic"]
    IDs <- unique((subset(categories, Class == "High"))$TRANSMITTERID)
    
  } else if(Class == "FemFem") {
    overlaps1 <- overlaps[Class == "FemaleFemale"]
    IDs <- unique((subset(categories, Class == "Female"))$TRANSMITTERID)
    
  } else {
    overlaps1 <- overlaps
    IDs <- unique(categories$TRANSMITTERID)
  }
  
  #==============================================================================#
  # Begin the main part of the function by subsetting down to specific individuals of interest
  social <- list() # create a list to store the outputs from below into
  # t <- 1
  for(i in 1:length(IDs)){
    indiv <- overlaps1[ID1 == IDs[i]| ID2 == IDs[i]]
    
    if(nrow(indiv) >= 1){
      # Determine the number of months they have overlaps for
      months <- unique(indiv$Month)
      
      indiv_list <- list() # create the list to save the results found here into
      for(u in 1:length(months)){
        if(u == 1){ # for the first run assign the month ID and prop present to NA
          indiv_list[[u]] <- data.table(TRANSMITTERID = IDs[i],
                                        Month = months[u],
                                        Class = Class,
                                        Prop_last_month = NA)
          
          # Then create the list of IDs that they are overlapping with this month
          indiv_month <- indiv[Month == months[u]] # subset down to the month of interest
          pastID <- unique(c(indiv_month$ID1, indiv_month$ID2)) # make a vector with the ID list
          pastID <- pastID[pastID != IDs[i]] # remove the individual of interests ID to stop any potential biases for the next step
          
        } else {
          indiv_month <- indiv[Month == months[u]] # subset down to the month of interest
          indivspast <- indiv_month[indiv_month$ID1 %in% pastID | indiv_month$ID2 %in% pastID] # subset again this time looking to see if individuals are present from last month
          
          indiv_list[[u]] <- data.table(TRANSMITTERID = IDs[i],
                                        Month = months[u],
                                        Class = Class,
                                        Prop_last_month = (nrow(indivspast)/length(pastID))) # proportion of individuals that are present this month from last month
          
          # Redo pastID so that it can then be examined for next month
          pastID <- unique(c(indiv_month$ID1, indiv_month$ID2)) # make a vector with the ID list
          pastID <- pastID[pastID != IDs[i]] # remove the individual of interests ID to stop any potential biases for the next step
        }
      }
      social[[i]] <- data.table(Reduce(rbind, indiv_list)) 
    }
    #t <- t +1 
  }
  
  social <- data.table(Reduce(rbind, social))
  print("everything is finished")
  return(social)
}
