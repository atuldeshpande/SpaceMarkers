#' @import spatstat
#' @import moranfast
#' @import pracma
#import description end
0



## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_kernel_outliers_for_sensitivity <- function(pattern, locs, pattern_threshold = 0.15, sigma = 10, kernelthreshold = 2, method = "Pattern_Threshold", outlier = "positive")
{
  numReps <- 100
  allwin <- spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)), yrange = c(min(locs$y),max(locs$y)))
  X <-spatstat.geom::ppp(x = locs$x, y = locs$y, window = allwin, marks = pattern)
  
  Kact = spatstat.explore::Smooth(X, at = "points", sigma = sigma, leaveoneout = T)
  Kvar <- Kmean <- Kact*0
  for (i in seq(1,numReps)){
      Xr = X
      spatstat.geom::marks(Xr) = spatstat.geom::marks(X)[sample(length(spatstat.geom::marks(X)))]
      K_i = spatstat.explore::Smooth(Xr, at = "points", sigma = sigma, leaveoneout = T)
      Kmean <- Kmean + K_i
      Kvar <- Kvar + K_i^2
  }
  mKvec <- sum(Kmean)/numReps/length(Kmean)
  sKvec <- sum(Kvar)/(numReps*length(Kmean) - 1)
  
  #Karr = sapply(seq(1,100), function(i) {Xr = X; spatstat.geom::marks(Xr) = spatstat.geom::marks(X)[pracma::randperm(1:length(spatstat.geom::marks(X)))]; temp = spatstat.explore::Smooth(Xr, at = "points", sigma = sigma, leaveoneout = T); return(temp)})
      
  #Kvec = unlist(Karr)
  #mKvec = mean(Kvec)
  #sKvec = sd(Kvec)
  
  Kact<-(Kact-mKvec)/sKvec
  return(Kact)
}

getOptimalSigmaThresh <- function(pattern, locs, sigmaVec, threshVec){
  allwin <- spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)), yrange = c(min(locs$y),max(locs$y)))
  X <-spatstat.geom::ppp(x = locs$x, y = locs$y, window = allwin, marks = pattern)
  mor_min <- Inf
  sigma_min <- 0
    for (i in 1:length(sigmaVec)) {
      K_i <- spatstat.explore::Smooth(X, at = "points", sigma = sigmaVec[i], leaveoneout = T)
      #mor_i <- ape::Moran.I(spatstat.geom::marks(X)-K_i, visium.dist.inv))
      mor_i <- moranfast::moranfast(spatstat.geom::marks(X)-K_i, locs$x, locs$y)
      if (abs(mor_i[[1]]-mor_i[[2]]) < mor_min) {
          sigma_min <- i
          mor_min <- abs(mor_i[[1]]-mor_i[[2]])
      }
    } 
  sigmaOpt1_ind <- sigma_min
  if (sigmaOpt1_ind>1&&sigmaOpt1_ind<length(sigmaVec)){
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind-1],sigmaVec[sigmaOpt1_ind+1],(sigmaVec[sigmaOpt1_ind+1] - sigmaVec[sigmaOpt1_ind-1])/10)
  }else if (sigmaOpt1_ind==1){
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind],sigmaVec[sigmaOpt1_ind+1],(sigmaVec[sigmaOpt1_ind+1] - sigmaVec[sigmaOpt1_ind])/10)
  }else{
    smallSigmaVec <- seq(sigmaVec[sigmaOpt1_ind-1],sigmaVec[sigmaOpt1_ind],(sigmaVec[sigmaOpt1_ind] - sigmaVec[sigmaOpt1_ind-1])/10)
  }
  
  mor_min <- Inf
  sigma_min <- 0
  for (i in 1:length(smallSigmaVec)) {
      K_i <- spatstat.explore::Smooth(X, at = "points", sigma = smallSigmaVec[i], leaveoneout = T)
      #mor_i <- ape::Moran.I(spatstat.geom::marks(X)-K_i, visium.dist.inv))
      mor_i <- moranfast::moranfast(spatstat.geom::marks(X)-K_i, locs$x, locs$y)
      if (abs(mor_i[[1]]-mor_i[[2]]) < mor_min) {
          sigma_min <- i
          mor_min <- abs(mor_i[[1]]-mor_i[[2]])
      }
  } 
  sigmaOpt1_ind <- sigma_min
  Kact2 = find_kernel_outliers_for_sensitivity(pattern = pattern, locs = locs, sigma = smallSigmaVec[sigmaOpt1_ind], method = "Kernel2", kernelthreshold = 0, outlier = "positive")
  K_i <- spatstat.explore::Smooth(X, at = "points", sigma = smallSigmaVec[sigmaOpt1_ind], leaveoneout = T)
  inds2<- (kronecker(matrix(1,1,length(threshVec)),Kact2)>threshVec)*1
  mor_min <- Inf
  mor_2_ind <- 0
  for (i in 1:length(threshVec)){
      inds2 <- which(Kact2>threshVec[i]) 
      mor_2 <- moranfast::moranfast(spatstat.geom::marks(X)[inds2]-K_i[inds2],locs$x[inds2],locs$y[inds2])
      if ((mor_2[[1]]-mor_2[[2]]) < mor_min){
          mor_min <- (mor_2[[1]]-mor_2[[2]])
          mor_2_ind <- i
      }
  }
  threshOpt1_ind <- mor_2_ind
  return(data.frame(sigmaOpt = smallSigmaVec[sigmaOpt1_ind], threshOpt = threshVec[threshOpt1_ind]))
}
#===================
#' getSpatialParameters
#' Calculate ...
#'
#' This function calculates ...
#'
#' @export
#'
#' @param spatialPatterns 	...
#'
#'

getSpatialParameters <- function(spatialPatterns){
  good_gene_threshold <- 3;  
  sigmaRes <- max(floor(min(diff(range(spatialPatterns$x)),diff(range(spatialPatterns$y)))/250),1)
  sigmaVec <- seq(2,40*sigmaRes,sigmaRes)
  threshVec <- seq(1,3,0.1)
  
  patternList <- colnames(spatialPatterns)[startsWith(colnames(spatialPatterns),"Pattern_")]
  optParams<-sapply(patternList, function(i) unlist(getOptimalSigmaThresh(pattern = spatialPatterns[,i], locs = data.frame(x = spatialPatterns$x, y = spatialPatterns$y), sigmaVec, threshVec)))
  return(optParams)
}
