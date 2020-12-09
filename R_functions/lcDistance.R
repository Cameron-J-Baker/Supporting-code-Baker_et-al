# The below function was orignially created by Vinay Udawyer and calculates the least cost utilisation distribution.
# Modifcations were made by Cameron Baker to allow the function to utilize VTrack format data along with a data frame
# detailing the location of acoustic hydrophone stations. This version currently only runs on R version 3.6.1 and will
# be incoroporated into the VTrack R package when it has been modified to work with later versions of R
# NOTE: THE CURRENT IMPLEMENTATION OF THIS FUNCTION IS NOT YET EFFIENCENT AND CAN HAVE A CONSIDERABLE RUN TIME DEPENDING
# ON THE NUMBER OF DETECTIONS THAT AN INDIVIDUAL HAS

lcDistance <- function(tag_data, receiver_loc, cost = NULL,  utm_epsg, ll_epsg = 4326, timestep = 60, h = 100, UDextent = 4, UDgrid = 100, cost.res = 100, directions = 16, ID = "test", ...){
  ####################################################################################
  ## tag_data       COA object produced from running the COA or COAV2 functions followed by the SnapCOAtoRiver function
  ## receiver_loc   Dataframe containing the location (lat/long) of the aoucsitc hydrophone stations
  ## cost           cost raster, with land pixel values of 1000 and sea as 1, if none provided
  ##                this is extracted and calculated using polygon downloaded from 'osmdata' package
  ## utm_epsg       EPSG code for CRS object with projected coordinate reference system (units m)
  ## ll_epsg        EPSG code for CRS object with unprojected coordinate reference system (units deg) [default 4326]
  ## timestep       timestep for Center of activity positions (see COA() function in VTrack)
  ## h              smoothing parameter for kernel density estimator in meters (default 100m)
  ## UDextent       a value controlling the extent of the grid used for the estimation of kernel density (default = 4) 
  ## UDgrid         a number giving the size of the grid on which the kernel density should be estimated (in meters; default 100)
  ## cost.res       resolution of cost layer (in degrees) calculated if not provided [default 100]
  ## directions     number of directional axes to consider when calculating 
  ##                least cost path (options: 4, 8, 16 [default])
  ####################################################################################
  
  ## load required libraries and set up CRSs
  sapply(c("lubridate","sf","dplyr","raster","gdistance","spatstat","maptools","rgeos","VTrack"), require, character.only=TRUE, warn.conflicts = F, quietly = T)
  ll<-CRS("+init=epsg:4326")
  utm<-CRS(paste0("+init=epsg:", 32754))
  
  combdata <-
    indiv  %>% 
    rename(Longitude = Longitude.coa.snap,
           Latitude = Latitude.coa.snap)
  
  statinfo <-
    receiver_loc %>% 
    st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326)
  
  ## Setup input tagdata; convert to spatial object and transform to projection
  spdata_ll <- 
    combdata %>%
    mutate(step = 0:(nrow(.)-1)) %>% 
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
  
  spdata_utm <-
    spdata_ll %>% 
    st_transform(crs = 32754)
  
  ## Extract landmass to calculate cost raster if not provided
  if(is.null(cost)){
    message("- No cost layer provided. Downloading coastline data from Open Street Map server")
    require("osmdata")
    tryCatch({
      polydat <-
        opq(bbox = (extent(statinfo) + 0.05)[c(1,3,2,4)]) %>% 
        add_osm_feature(key = 'admin_level', value = '6') %>% 
        osmdata_sf
      
      poly <- 
        polydat$osm_multipolygons %>% 
        st_transform(utm_epsg) %>%
        st_crop(extent(statinfo %>% st_transform(utm_epsg)) + 5000)
      
      cost.ras<-rasterize(poly, raster(extent(statinfo %>% st_transform(utm_epsg)) + 5000, res = cost.res), 1000)
      cost.ras[is.na(values(cost.ras))] <- 1
      projection(cost.ras) <- utm
      
    }, error=function(e){message("Error: no land in sight!\nConsider adding your own cost layer")})
  }else{
    cost.in_utm <- raster::projectRaster(cost, crs=utm, method = "ngb")
    cost.ras <- resample(cost.in_utm, raster(extent(cost.in_utm), res = cost.res), method = "ngb")
    projection(cost.ras) <- utm
  }
  
  ## Produce transition matrices, and correct for distortion
  tryCatch({
    message("- Calculating transition matrices for least cost path estimation")
    trCost <- transition(1/cost.ras, mean, directions = directions)
    trCost <- geoCorrection(trCost, type = "c")
  }, error=function(e){message("Error in calculating Transition layer")})
  
  ## Construct shortest path trajectories between sequence of detection steps
  outdat <- 
    spdata_ll %>% 
    as_Spatial %>% 
    as_tibble %>% 
    rename(Longitude = coords.x1, Latitude = coords.x2) %>% 
    mutate(Distance_m = NA) %>% 
    dplyr::select(-step)
  
  tryCatch({
    message("- Constructing least cost path trajectories between consecutive detections")
    traj<-list()
    for(i in 1:max(spdata_utm$step)){
      if(i %in% 1){pb <- txtProgressBar(min=1, max=max(spdata_utm$step), style=3)}
      origin <- spdata_utm %>% filter(step %in% (i-1)) %>% as_Spatial
      goal <- spdata_utm %>% filter(step %in% i) %>% as_Spatial
      traj[[i]]<-shortestPath(trCost, origin, goal, output = "SpatialLines")
      outdat$Distance_m[i+1]<-as.numeric(costDistance(trCost, origin, goal))
      setTxtProgressBar(pb, i)
    }
    trajectory_utm <- do.call(rbind, traj)
    trajectory_ll <- st_as_sf(trajectory_utm) %>% st_transform(ll_epsg)
  }, error=function(e){message("\nError in calculating least cost path trajectories")})
  
  # Write a copy of the least cost trajectory for each individual to speed up analysis for calculating H
  #write_rds(trajectory_utm, paste("Data/lc_HR/lc_trajectories/",ID,".rds", sep = ""))
  ## Calculate Kernel Density from trajectory
  UDwin <- owin(xrange = (extent(cost.in_utm) + diff(extent(cost.in_utm)[1:2])/UDextent)[1:2],
                yrange = (extent(cost.in_utm) + diff(extent(cost.in_utm)[3:4])/UDextent)[3:4])
  dimyx <- c(diff(UDwin$yrange)/UDgrid,
             diff(UDwin$xrange)/UDgrid)
  
  # point UD
  tryCatch({ 
    message("\n- Estimating Kernel Density from COA positions")
    coa.data <- 
      combdata %>% 
      st_as_sf(coords=c("Longitude","Latitude"), crs=4326) %>% 
      st_transform(32754) %>% as_Spatial() %>% as_tibble() %>% 
      rename(Longitude = coords.x1, Latitude = coords.x2)
    suppressWarnings(pt.ppp <- ppp(x = coa.data$Longitude, y = coa.data$Latitude, window = UDwin))
    suppressWarnings(UDras.pt <- raster(density(pt.ppp, sigma = h*2, method = "C", dimyx = dimyx)))
    values(UDras.pt) <- abs(values(UDras.pt)/max(values(UDras.pt)) - 1) * 100
    projection(UDras.pt) <- utm
  }, error=function(e){message("Error in estimating UD from COA positions")})
  
  # trajectory UD
  tryCatch({
    message("- Estimating Kernel Density from least cost path trajectory")
    traj.psp <- as.psp(trajectory_utm, window = UDwin)
    UDras.traj <- raster(density(traj.psp, sigma = h, method = "C", dimyx = dimyx))
    values(UDras.traj) <- abs(values(UDras.traj)/max(values(UDras.traj)) - 1) * 100
    projection(UDras.traj) <- utm
  }, error=function(e){message("Error in estimating UD for trajectories")})
  
  # Set the NA values in each raster to 0
  # This way if the individual is only found in a single cell and subsequently has no trajectory UD only the point one is used below
  values(UDras.pt)[is.na(values(UDras.pt))] <- 0
  values(UDras.traj)[is.na(values(UDras.traj))] <- 0
  
  # Merge trajectory UD with point UD and remove land areas
  #tryCatch({ 
  if (sum(values(UDras.traj)) == 0){ # If the individual is only active within a single grid cell, create the UD off just that one point
    UDras <- mean(UDras.pt)
    mask.ras <- resample(cost.ras, UDras, method = "ngb")
    values(mask.ras)[values(mask.ras) > 1] <- NA
    UDras.m_utm <- raster::mask(UDras, mask = mask.ras)
    UDras.m_ll <- projectRaster(UDras.m_utm, crs = ll)
  } else {
    UDras <- mean(UDras.traj, UDras.pt)
    mask.ras <- resample(cost.ras, UDras, method = "ngb")
    values(mask.ras)[values(mask.ras) > 1] <- NA
    UDras.m_utm <- raster::mask(UDras, mask = mask.ras)
    UDras.m_ll <- projectRaster(UDras.m_utm, crs = ll)
  }
  #}, error=function(e){message("Error in estimating UD")})
  
  # calculate UD areas
  #tryCatch({
  message("- Estimating UD area for 50% and 95% contours")
  UD50 <- UDras.m_utm < 50; UD50[values(UD50) %in% 0] <- NA
  UD50.area <- gArea(rasterToPolygons(UD50, dissolve=T))
  
  UD95 <- UDras.m_utm < 95; UD95[values(UD95) %in% 0] <- NA
  UD95.area <- gArea(rasterToPolygons(UD95, dissolve=T))
  
  #}, error=function(e){message("Error in estimating UD areas")})
  
  ## return list output with step distances and spatial trajectory file
  out<-list(tagdata = outdat, 
            coa.data = coa.data,
            kernel.areas = data.frame(UD50_m2 = UD50.area, UD95_m2 = UD95.area),
            lc.traj = trajectory_ll, 
            UD.raster = UDras.m_ll, 
            cost.raster = projectRaster(cost.ras, crs = ll, method = "ngb"))
  
  return(out)
}
