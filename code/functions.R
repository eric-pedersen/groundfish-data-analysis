
#function to calculate geometric mean abundance
CalcZeroInfGeomDens = function(x){
  prob_obs= sum(x>0)/length(x)
  geom_mean_dens = ifelse(prob_obs>0,exp(mean(log(x[x>0]))),0)
  return(geom_mean_dens*prob_obs)
}


community.sync.window <- function (d, twin=5, nrands=999) {
  years = as.numeric(rownames(d))
  starts=seq(from=1, to=length(years)-twin+1, by=1)
  ends=starts+twin-1
  mid=starts+(twin-1)/2
  
  results=matrix(nrow=length(starts), ncol=6, NA)
  for (i in 1:length(starts)) {
    locs=which(rownames(d) %in% years[starts[i]:ends[i]])
    c = community.sync(d[locs,], nrands=nrands)
    results[i,] = c(years[starts[i]], years[ends[i]],years[mid[i]], c$meancorr, c$obs, c$pval)
  }
  colnames(results)=c("start", "end","mid", "meancorr", "sync", "pval")
  results=as.data.frame(results)
  return(results)
}


CalcGeoDistMat = function(lat, long){
  # This function calculates the geodesic distance matrix between all points in a set of 
  # lat/long coordinates.
  stopifnot(is.numeric(lat)&is.numeric(long))
  n_data= length(lat)
  lat_data = expand.grid(lat,lat)
  long_data = expand.grid(long,long)
  dist_calc = geoDist(lat1=lat_data[,1],lat2=lat_data[,2],
                      lon1=long_data[,1],lon2=long_data[,2],NAOK=T)
  dist_mat = matrix(dist_calc,nrow=n_data)
  return(dist_mat)
}

CalcSpaceDepthDists = function(x,use_geodesic){
  #Calculates in a single function both depth and distance dissimilarity between polygons
  if(use_geodesic){
    spatial_dist = CalcGeoDistMat(x$LAT_DEC,x$LONG_DEC)
    spatial_dist = spatial_dist[lower.tri(spatial_dist,diag=F)]
  }else{
    spatial_dist = as.vector(dist(cbind(x$LAT_DEC,x$LONG_DEC)))
  }
  depth_dist = as.vector(dist(x$Depth))
  return(data.frame(space_dist=spatial_dist, depth_dist=depth_dist))
}

EstimateDistR2 = function(x){
  # Estimates the partial R^2 explained of species dissimilarity
  # by differences in distance and depth (after transformations)
  distance_scl = x$distance_scl
  boot_data = x
  boot_data = subset(boot_data, !is.na(distance_scl))
  distance_scl = distance_scl[!is.na(distance_scl)]
  var_decomp = varpart(Y=distance_scl, ~space_dist_scl,~depth_dist_scl,
                       data=boot_data)
  return(var_decomp$part$indfract$Adj.R.square[c(1,3)])
}

CalculateRSquareValues = function(x){
  #Wrapper function for EstimateDistR2 that also calculates
  #jackknife se values for the estimates
  var_frac = EstimateDistR2(x)
  if(boot_errors){
    unique_sites = unique(c(x$site_1,x$site_2))
    n_sites=length(unique_sites)
    jack_R_sqr = matrix(rep(0,times=n_sites*2),ncol=2)
    jack_mean = rep(0,times=n_sites)
    jack_var =rep(0,times=n_sites)
    for(i in 1:n_sites){
      indexes = x$site_1!=i&x$site_2!=i
      current_data = x[indexes,]
      jack_mean[i] = mean(current_data$distance_scl,na.rm=T)
      jack_var[i] =var(current_data$distance_scl,na.rm=T)
      jack_R_sqr[i,]  =  EstimateDistR2(current_data)
    }
    jack_R_sqr_est = colMeans(jack_R_sqr)
    jack_mean_est = mean(jack_mean)
    jack_var_est = mean(jack_var)
    jack_R_sqr_se = sqrt((n_sites-1)*(colMeans(jack_R_sqr^2)-jack_R_sqr_est^2))
    jack_mean_se = sqrt((n_sites-1)*var(jack_mean))
    jack_var_se = sqrt((n_sites-1)*var(jack_var))
  }
  model_slopes = coef(lm(distance_scl~space_dist_scl+depth_dist_scl,
                         data=x))
  output_data = data.frame(coef = c("Mean community distance","variance",
                                    "Spatial Distance","Depth"),
                           R2_value = c(mean(x$distance_scl,na.rm=T),
                                        var(x$distance_scl,na.rm=T), 
                                        var_frac),
                           slope_value = c(mean(x$distance_scl,na.rm=T),
                                           var(x$distance_scl,na.rm=T),
                                           model_slopes[2:3]))
  output_data$coef= factor(output_data$coef, 
                           levels= c("Mean community distance",
                                     "variance","Spatial Distance",
                                     "Depth"))
  if(boot_errors){
    
    R2_value = c(jack_mean_est,jack_var_est, jack_R_sqr_est)
    output_data$se = c(jack_mean_se,jack_mean_se,jack_R_sqr_se)
  }
  return(output_data)
}

R2NameLabeller =function(variable, value){
  #Function to aid in plotting. Parses various expressions into caption labels.
  dist_labels= list("Mean community distance"=expression("Mean scaled community distance"), 
                    "variance" = expression("Variance of community distance"),
                    "Spatial Distance"= expression(paste('Partial ',R^2, ' explained by distance')),
                    "Depth"= expression(paste('Partial ',R^2, ' explained by depth differential'))
  )
  R2_labels = list(
    "R2_value" = "Total community",
    "R2_top4"="four most abundant",
    "R2_nontop"="Remaining species"
  )
  if(variable=="coef"){
    return(dist_labels[value])
  }else{
    return(R2_labels[value])
  }
}



space_trans = function(x) {
  #Returns centered, scaled, log-10 transformation of spatial distances
  return(scale(log10(x)))
}

depth_trans = function(x) {
  #Returns centered, scaled depth dissimilarities
  return(scale(x,center=T,scale=T))
}

comm_trans = function(x, logt=T,scale_to_max=T,com_dist="euclidean", pad_val =0.01,
                      center_val=F,scale_val = F) {
  # a modified logistic transformation for species dissimilarity. Transforms a value
  # ranging from 0 to some maximum value into a range from -inf to +inf. 
  if(com_dist=="euclidean"){
    max_val = sqrt(2)
  }else if(com_dist=="bray"){
    max_val = 1 
  }else stop("community distance ",com_dist, " is not a recognized distance currently") 
  new_x = x
  if(scale_to_max) new_x = (new_x+pad_val)/(max_val+pad_val-new_x)
  if(logt) new_x =log(new_x)
  new_x= scale(new_x,center=center_val,scale=scale_val)
  return(new_x)
}



PlotMultipleGgplotObjs <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

