#' @import spatstat
#' @importFrom ape where
#' @import pracma
#import description end
0

## author: Atul Deshpande
## email: adeshpande@jhu.edu

find_kernel_outliers_for_sensitivity <- function(pattern,locs,
                                                 pattern_threshold=0.15,
                                                 sigma = 10,kernelthreshold = 2,
                                                 method = "Pattern_Threshold",
                                                 outlier = "positive",...)
{
    allwin<-spatstat.geom::owin(xrange = c(min(locs$x),max(locs$x)),
                                yrange=c(min(locs$y),max(locs$y)))
    X<-spatstat.geom::ppp(x=locs$x,y=locs$y,window=allwin,marks=pattern)
    Kact<-spatstat.explore::Smooth(X,at ="points",sigma=sigma,leaveoneout=TRUE)
    Karr<-vapply(seq(1,100), function(i){Xr<-X;
    spatstat.geom::marks(Xr)<-spatstat.geom::marks(X)[pracma::randperm(
        seq_len(length(spatstat.geom::marks(X))))];
    temp<-spatstat.explore::Smooth(Xr,at="points",sigma=sigma,leaveoneout=TRUE);
    return(temp)}, numeric(length(Kact)))
    Kvec <- unlist(Karr)
    mKvec <- mean(Kvec)
    sKvec <- sd(Kvec)
    Kact<-(Kact-mKvec)/sKvec
    return(Kact)
}
### Add sigVec and threshVec
getOptimalSigmaThresh <- function(pattern, locs,sigVec,threshVec,...){
    visium.dist <- as.matrix(dist(locs))
    visium.dist.inv <-1/visium.dist
    diag(visium.dist.inv) <- 0
    allwin<-spatstat.geom::owin(xrange=c(min(locs$x),max(locs$x)),
                                yrange=c(min(locs$y),max(locs$y)))
    X<-spatstat.geom::ppp(x=locs$x,y=locs$y,window=allwin,marks=pattern)
    Ks<-vapply(sigVec,function(i) spatstat.explore::Smooth(X,at="points",
                                                           sigma=i,
                                                           leaveoneout=TRUE),
               numeric(X$n))
    mor_1<-vapply(seq_len(length(sigVec)),function(i) unlist(ape::Moran.I(
        spatstat.geom::marks(X)-Ks[,i], visium.dist.inv)),numeric(4))
    sigOpt1_ind <- which.min(abs(unlist(mor_1[1,])-unlist(mor_1[2,])))
    if (sigOpt1_ind>1&&sigOpt1_ind<length(sigVec)){
        smallsigVec<-seq(sigVec[sigOpt1_ind-1],sigVec[sigOpt1_ind+1],
                         (sigVec[sigOpt1_ind+1]-sigVec[sigOpt1_ind-1])/10)
    }else if (sigOpt1_ind==1){
        smallsigVec <- seq(sigVec[sigOpt1_ind],sigVec[sigOpt1_ind+1],
                           (sigVec[sigOpt1_ind+1] - sigVec[sigOpt1_ind])/10)
    }else{
        smallsigVec <- seq(sigVec[sigOpt1_ind-1],sigVec[sigOpt1_ind],
                           (sigVec[sigOpt1_ind] - sigVec[sigOpt1_ind-1])/10)
    }
    smallKs<-vapply(smallsigVec,function(i) spatstat.explore::Smooth(
        X,at="points", sigma = i, leaveoneout = TRUE), numeric(X$n))
    smallmor_2<-vapply(seq_len(length(smallsigVec)), function(i) unlist(
        ape::Moran.I(spatstat.geom::marks(X)-smallKs[,i],visium.dist.inv)),
        numeric(4))
    sigOpt1_ind <- which.min(abs(unlist(smallmor_2[1,])-unlist(smallmor_2[2,])))
    Kact2<-find_kernel_outliers_for_sensitivity(pattern=pattern,locs=locs,
                                                sigma=smallsigVec[sigOpt1_ind],
                                                method="Kernel2",
                                                kernelthreshold=0,
                                                outlier="positive")
    inds2<-(kronecker(matrix(1,1,length(threshVec)),Kact2)>threshVec)*1
    mor_2_ind<-vapply(seq_len(length(threshVec)),function(j){
        visium.dist.inv_1<- visium.dist.inv;visium.dist.inv_1[inds2[,j]==0]<-0;
        unlist(ape::Moran.I(spatstat.geom::marks(X)-smallKs[,sigOpt1_ind],
                            visium.dist.inv_1))}, numeric(4))
    threshOpt1_ind<-which.min(abs(unlist(mor_2_ind[1,])-unlist(mor_2_ind[2,])))
    return(data.frame(
        sigmaOpt=smallsigVec[sigOpt1_ind],threshOpt=threshVec[threshOpt1_ind]))
}

find_pattern_hotspots <- function(
        spPatterns, params = NULL,
        patternName = "Pattern_1",
        outlier = "positive",...){
    
    if (is.null(params)){
        if (!exists('sigVec')){
            sigmaRes <- max(floor(min(diff(range(spPatterns$x)),
                                      diff(range(spPatterns$y)))/250),1)
            sigVec <- seq(2,40*sigmaRes,sigmaRes)
        }
        if (!exists('threshVec')){
            threshVec <- seq(1,3,0.1)
        }
        ### Added sigVec and threshVec
        params <- getOptimalSigmaThresh(spPatterns[,patternName], locs = spPatterns[,c("x","y")],sigVec,threshVec)
    }   
    sigmaPair <- params[[1]]
    kernelthreshold <- params[[2]]
    
    allwin <- spatstat.geom::owin(
        xrange = c(min(spPatterns$x),max(spPatterns$x)),yrange =c(
            min(spPatterns$y),max(spPatterns$y)))
    patternVector <- as.matrix(spPatterns[,patternName])
    X <-spatstat.geom::ppp(
        x=spPatterns$x,y = spPatterns$y,window = allwin,marks = patternVector)
    Kact1 <- spatstat.explore::Smooth(
        X, at = "points", sigma = sigmaPair[1],leaveoneout = TRUE)
    Karr1 <- vapply(seq(1,100),function(i){
        Xr<-X;
        spatstat.geom::marks(Xr) <- spatstat.geom::marks(X)[pracma::randperm(
            seq_len(length(spatstat.geom::marks(
                X))))];temp <- spatstat.explore::Smooth(
                    Xr, at="points", sigma = sigmaPair[1],leaveoneout=TRUE); 
                return(temp)}, numeric(length(Kact1)))
    Karr1 <- unlist(Karr1)
    mKvec <- mean(Karr1)
    sKvec <- sd(Karr1)
    upthresh <- mKvec+kernelthreshold*sKvec
    lothresh <- mKvec-kernelthreshold*sKvec
    if (outlier == "positive"){
        ind1 <- which(Kact1 > upthresh[1])
    }
    else if (outlier == "two.sided")
    {
        ind1 <- which((Kact1 > upthresh)|(Kact1 < lothresh))
    }
    region <- array(NA, length(Kact1))
    region[ind1] <- patternName
    return(region)
}
