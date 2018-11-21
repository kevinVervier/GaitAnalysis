#-------------------------#
# Kinect data analysis V1 #
#-------------------------#

# This package aims at:
# * processing raw Kinect data (start/end trimming, glitch detection, rotation to the sagittal plan)
# * Extract Fourier components from the recorded signal
# see Demo at the very end of the script

if(!require(randomForest)){
  install.packages('randomForest')
  require(randomForest)
}

###############################################################################################################################
# Main function from raw data to Fourier matrix
# Inputs:
#         - file: .csv file name that need to be processed
#         - (optional) nfreq: number of frequencies in Fourier (default: 10)
#         - (optional) bodysite: condisered body sites (default: 1 to 25 --> all)
#         - (optional) rmSteady: boolean for performing steady start/end removal (default = TRUE)
#         - (optional) rmOdd: boolean for performing odd start/end removal (default = TRUE)
#         - (optional) rmGlitch: boolean for performing glitch removal (default = TRUE)
#         - (optional) rotated: boolean for performing rotation (default = TRUE)
#         - (optional) output.file: .csv file where scores will be output (if not, display in standard output)

mainProcess <- function(file.std = NULL, file.ht = NULL, nfreq = 10, bodysite = 1:25, rmSteady = TRUE, rmOdd = TRUE, rmGlitch = TRUE, rotated = TRUE, verbose = FALSE, output.file = NULL){

  ################################
  # read standard walking data
  if(is.null(file.std)){
    tmp = strsplit(file.ht,'HT')[[1]][1]
    file.std = paste(tmp,'W.csv',sep='')
  }
  if(file.exists(file.std)){
    x = loadData(file.std)
    cat('Standard Walk: Data loaded...\n')
    if(verbose) dataViz(x)
    # remove steady start/end
    if(rmSteady){ 
      x = removeSteady(x)
      cat('Standard Walk: Steady parts removed...\n')
      if(verbose) dataViz(x)
    }
    # remove odd start/end
    if(rmOdd){ 
      x = removeOdd(x)
      cat('Standard Walk: Odd start/end removed...\n')
      if(verbose) dataViz(x)
    }
    # remove glitches
    if(rmGlitch){
      x = removeGlitch(x)
      cat('Standard Walk: Glitches removed...\n')
      if(verbose) dataViz(x)
    }
    # rotate data
    if(rotated){
      x = rotateData(x)
      cat('Standard Walk: Data rotated...\n')
      if(verbose) dataViz(x)
    }
    # apply Fourier
    fourier.mat.std = applyFourier(x,nfreq,bodysite)
    cat('Standard Walk: Fourier analysis done...\n')
    # derive prediction for standard walking
    load('../output/RF-std-model_and_scales.Rdata')
    s2 = matrix(fourier.mat.std[feat_imp],nrow=1)
    s2=scale(s2,attr(s,"scaled:center"),attr(s,"scaled:scale"))
    s2[which(is.na(s2))] = 0
    s2[which(s2 == Inf)] = 0
    x=predict(m,s2,type='prob')
    preds.std = x[,2]
  }else{
    preds.std = 0.5
  }
  
  ################################
  # read corresponding heel-toe data
  if(is.null(file.ht)){
    tmp = strsplit(file.std,'W')[[1]][1]
    file.ht = paste(tmp,'HT.csv',sep='')
  }
  if(file.exists(file.ht)){
    x = loadData(file.ht)
    cat('HT walk: Data loaded...\n')
    if(verbose) dataViz(x)
    # remove steady start/end
    if(rmSteady){ 
      x = removeSteady(x)
      cat('HT walk: Steady parts removed...\n')
      if(verbose) dataViz(x)
    }
    # remove odd start/end
    if(rmOdd){ 
      x = removeOdd(x)
      cat('HT walk: Odd start/end removed...\n')
      if(verbose) dataViz(x)
    }
    # remove glitches
    if(rmGlitch){
      x = removeGlitch(x)
      cat('HT walk: Glitches removed...\n')
      if(verbose) dataViz(x)
    }
    # rotate data
    if(rotated){
      x = rotateData(x)
      cat('HT walk: Data rotated...\n')
      if(verbose) dataViz(x)
    }
    # apply Fourier
    fourier.mat.ht = applyFourier(x,nfreq,bodysite)
    cat('HT walk: Fourier analysis done...\n')
    #output
    graphics.off()
    # derive prediction for heel-toe
    load('../output/RF-ht-model_and_scales.Rdata')
    s2 = matrix(fourier.mat.ht[feat_imp],nrow=1)
    s2=scale(s2,attr(s,"scaled:center"),attr(s,"scaled:scale"))
    s2[which(is.na(s2))] = 0
    s2[which(s2 == Inf)] = 0
    x=predict(m,s2,type='prob')
    preds.ht = x[,2]
  }else{
    preds.ht = 0.5
  }
  
  # ensemble model predicted score
  load('../output/hybrid_model.Rdata')
  fitted.results <- predict(m,newdata=data.frame('w.score'=preds.std,'ht.score' = preds.ht),type='response')
  
  if(is.null(output.file)){
    cat('Standard walk score:',preds.std,'\n')
    cat('Heel-toe walk score:',preds.ht,'\n')
    cat('Hybrid score:',fitted.results,'\n')
  }else{
    #store output
    OUTPUT = c(paste(strsplit(file.std,'_')[[1]][-length(strsplit(file.std,'_')[[1]])],collapse='_'),preds.std,preds.ht,fitted.results)
    names(OUTPUT) = c('individual','standard_walk_score','heel_toe_score','hybrid_score')
    write.table(OUTPUT,file=output.file,quote = FALSE,col.names = FALSE,sep=',')
  }
  

}


##################################
# Function for pre-processing data
# Warning: only consider .csv files as inputs

loadData <- function(file){
  
  #set default header with body points and coordinates
  bodyPoints = c('SpineBase','SpineMid','Neck','Head','ShoulderLeft','ElbowLeft','WristLeft','HandLeft','ShoulderRight','ElbowRight','WristRight','HandRight',
                 'HipLeft','KneeLeft','AnkleLeft','FootLeft','HipRight','KneeRight','AnkleRight','FootRight','SpineShoulder','HandTipLeft','ThumbLeft','HandTipRight','ThumbRight')
  coordinates = c('X','Y','Z')  
  # combine body points and coordinates for the final header
  header = as.vector(t(outer(bodyPoints,coordinates,paste, sep=".")))
  
  #init output
  x = NULL
  #read file
  tmp = read.csv2(file = file,header=F,stringsAsFactors = F)
  #check for two-lines header
  if(tmp[1,1] == 'SpineBase'){ 
    x = tmp[-(1:2),-76] #remove it if yes and also remove the 76th column which is NA...
  }else{
    x = tmp[,-76] #remove the 76th column = NA
  }
  
  #need to convert the matrix content into numeric values
  x=sapply(x,as.numeric)
  return(x)
}

##################################################################
# Function allwing to visualize data for all bodyPoints and 3 axis

dataViz <- function(x){
  #################
  # define colors #
  
  #plot X-axis values for 25 body points
  COL = rep(0,25)
  LTY = rep(1,25)
  #spinal+head
  COL[c(1:4,21)] = colorRampPalette(c("darkgreen", "green"))(5)
  #arms
  COL[c(5:8)] = colorRampPalette(c("orchid4", "orchid1"))(4)
  LTY[c(5:8)] = 2
  COL[c(9:12)] = colorRampPalette(c("orchid4", "orchid1"))(4)
  LTY[c(9:12)] = 3
  #legs
  COL[c(13:16)] = colorRampPalette(c("dodgerblue4", "dodgerblue1"))(4)
  LTY[c(13:16)] = 2
  COL[c(17:20)] = colorRampPalette(c("dodgerblue4", "dodgerblue1"))(4)
  LTY[c(17:20)] = 3
  #hands
  COL[c(22:23)] = c("firebrick","firebrick3")
  LTY[c(22:23)] = 2
  COL[c(24:25)] = c("firebrick","firebrick3")
  LTY[c(24:25)] = 3
  
  par(mfrow = c(2,2))
  #x-axis
  plot(x[,1],type='l',ylim=range(x[,3*(1:25)-2]),col=COL[1],ylab="x-axis value",xlab="Time")
  for(j in 2:25){
    lines(x[,3*j-2],col=COL[j],lty=LTY[j])
  }
  #y-axis
  plot(x[,2],type='l',ylim=range(x[,3*(1:25)-1]),col=COL[1],ylab="y-axis value",xlab="Time")
  for(j in 2:25){
    lines(x[,3*j-1],col=COL[j],lty=LTY[j])
  }
  #z-axis
  plot(x[,3],type='l',ylim=range(x[,3*(1:25)]),col=COL[1],ylab="z-axis value",xlab="Time")
  for(j in 2:25){
    lines(x[,3*j],col=COL[j],lty=LTY[j])
  }
  #legend
  plot(1:10,type = "n", axes=F, xlab="", ylab="")
  legend(1, 5,ncol=2,col = c("green","orchid1","dodgerblue","firebrick","white","black","black","white"), lty = c(0,0,0,0,0,2,3,0),legend = c("Spinal+Head","Arms","Legs","Hands","","Left","Right",""),pch=c(16,16,16,16,NA,NA,NA,NA),lwd=c(0,0,0,0,0,2,2,0),box.lwd = 0,box.col = "white",bg = "white")
  
  scan()
}

######################################################
# Function computing the first derivative along Y-axis

computeDerivative<-function(x){
  derivative = rep(0,length(x)-1)
  derivative = 0.5*(x[2:length(x)]-x[1:(length(x)-1)])
  return(derivative)
}

########################################################
# Function returning the index of the first moving point
removeStart<-function(y,thresh=10^-3){
  i = 1
  tmp=y[i]
  while(abs(tmp)<thresh & i < length(y)){
    i = i+1
    tmp=y[i]
  }
  return(i)
}

########################################################
# Function returning the index of the last moving point
removeEnd<-function(y,thresh=10^-3){
  i = 1
  tmp=rev(y)[i]
  while(abs(tmp)<thresh & i < length(y)){
    i = i+1
    tmp=rev(y)[i]
  }
  return(length(y)-i)
}


##############################################
# Function trimming the steady start and end
removeSteady <- function(x){
  
  #compute the first derivative along Y-axis
  y = apply(x,2,computeDerivative)
  # trim steady start and end for each record
  start=removeStart(y[,2])   #process by removing steady start
  end=removeEnd(y[,2])   #process by removing steady end
  if(start < end){ 
    x.proc= x[c(start:end),] #keep only moving time
  }else{
    x.proc = NULL
  }
  return(x.proc)
}

##################################################
# Function detecting odd starts (too much movement)
# based on median value +/- 2*sd
pruneOddMvtStart <- function(x,thresh=2){
  med = apply(x,2,mean)
  inf = med -thresh*apply(x,2,sd)
  max = med + thresh*apply(x,2,sd)
  
  i = 1
  tmp=x[i,]
  while( any(tmp < inf) | any(tmp > max)){
    i = i+1
    tmp=x[i,]
  }
  return(i+1)
}

################################################
# Function detecting odd ends (too much movement)
# based on median value +/- 2*sd
pruneOddMvtEnd <- function(x,thresh = 2){
  med = apply(x,2,mean)
  inf = med -thresh*apply(x,2,sd)
  max = med + thresh*apply(x,2,sd)
  
  i = nrow(x)
  tmp=x[i,]
  while( any(tmp < inf) | any(tmp > max)){
    i = i-1
    tmp=x[i,]
  }
  return(i-1)
}

#####################################################
# Function detecting and removing odd starts and ends
removeOdd <- function(x){
  if(is.null(dim(x))) return(NULL)
  # trim start and end
  start=pruneOddMvtStart(x,thresh=2) #detect odd starts
  end=pruneOddMvtEnd(x,thresh=2)     #detect odd ends
  x.proc = x[c(start:end),] #remove
  
  return(x.proc)
}

##########################################
# Function detecting and removing glitches 

removeGlitch <- function(x){
  if(is.null(dim(x))) return(NULL)
  #compute the first derivative along Y-axis
  y = apply(x,2,computeDerivative)
  #look for the point where, on average, each bodypoint has the highest variation
  tmp = max(apply(as.matrix(y,ncol=ncol(y)),1,function(i) mean(abs(i))))
  #detect glitch using a percentile/outlier
  idx = (tmp > 0.15) # empirical value found during February 2016 study
  # the output will be almost identical to the input
  x.proc = x
  #correct detected glitches
  if(idx){
    tmp = y
    j = which.max(apply(tmp,1,function(i) mean(abs(i)))) + 1
    if(j > 0.5*nrow(x.proc)){                     # remove the smaller part of the signal
      x.proc = x.proc[-c(j:nrow(x.proc)),]
    }else{
      x.proc = x.proc[-c(1:j),]
    } 
  }
  return(x.proc)
}

#########################################################################################################
# get center of mass direction (SpineBase point, instead of HipCenter that does not exist in our data)
getDirectionZ <-function(x){
  if(nrow(x) > 2){
    start = apply(x[1:(round(0.1*nrow(x))),1:3],2,mean) 
  }else{
    start = apply(x[1:2,1:3],2,mean) 
  }
  return(-start/sqrt(sum(start^2)))
}

################################################
#get y-axis given the direction SpineBase-Neck
getDirectionY <-function(x){
  spineBase = apply(x[,1:3],2,mean) 
  neck = apply(x[,7:9],2,mean)
  return((neck-spineBase)/sqrt(sum((neck-spineBase)^2)))
}

############################################################
# cross product function for 3D vectors
# Compute the vector cross product between x and y, and return the components indexed by i.
CrossProduct3D <- function(x, y, i=1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)
  
  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1
  
  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

#######################################################
#get x-axis based on y and z axes crossproduct
getDirectionX <- function(y,z){
  return(CrossProduct3D(y,z,i=1:3))
}

########################################################
# Function rotating data in the sagittal plan
rotateData <- function(x){
  if(is.null(dim(x))) return(NULL)
  #compute direction along z-axis
  z_directions = getDirectionZ(x)
  #compute direction along y-axis
  y_directions = getDirectionY(x)
  #compute direction along x-axis
  x_directions = getDirectionX(y_directions,z_directions)
  
  #apply rotation to the data
  x.proc = list()
  #bind x,y,z directions in a rotation matrix
  rotate = cbind(x_directions, y_directions, z_directions)
  tmp = NULL
  for(j in 1:(ncol(x)/3)){
    # rotate data
    tmp = cbind(tmp, x[,3*j-c(2:0)]%*%rotate)
  }
  x.proc = tmp
  return(x.proc)
}

###################################################################################3
# Fourier analysis of records signal
# Function computing the frequency domain for nfreq and bodysite(s)
computeFreq <- function(x,nfreq = 10,bodysite = c(1:25)){
  #remove trend (constant component) ~ bias in signal
  tmp = apply(x[,sort(3*bodysite-rep(c(0,1,2),length(bodysite)))],2,function(y)lm(y ~ c(1:length(y))))
  detrended.trajectory <- lapply(tmp,function(y)y$residuals)
  #apply Fourier analysis
  X.k <- lapply(detrended.trajectory,fft)
  #get the power of each component (Module)
  mod = lapply(X.k,function(y)Mod(y[2:(nfreq+1)])^2)
  #output matrix: row = bodysite, col = 3*nfreq
  X = t(matrix(unlist(mod),ncol=length(bodysite),nrow=nfreq*3))
  return(X)
}

########################################################################
# Function converting temporal data into frequency data using Fourier 
applyFourier <- function(x,nfreq = 10,bodysite = c(1:25)){
  if(is.null(dim(x))) return(rep(0,nfreq * length(bodysite) * 3))
  #set default body points
  bodyPoints = c('SpineBase','SpineMid','Neck','Head','ShoulderLeft','ElbowLeft','WristLeft','HandLeft','ShoulderRight','ElbowRight','WristRight','HandRight',
                 'HipLeft','KneeLeft','AnkleLeft','FootLeft','HipRight','KneeRight','AnkleRight','FootRight','SpineShoulder','HandTipLeft','ThumbLeft','HandTipRight','ThumbRight')[bodysite]
  
  # compute frequency component for each record
  data = computeFreq(x,nfreq=nfreq,bodysite=bodysite)
  data = as.vector(t(data))
  #compute the new feature names
  tmp2 = outer(outer(c('X','Y','Z'),paste("Freq",1:nfreq),paste,sep='.'),bodyPoints[bodyPoints != '' & !is.na(bodyPoints)],paste,sep='.')
  tmp3 = NULL
  for(i in 1:dim(tmp2)[3]) tmp3 = c(tmp3, as.vector(t(tmp2[,,i])))
  names(data) = tmp3
  return(data)
}


###################################################################
#                DEMO (single file)                               #
#data = mainProcess(file.std = '../data/144/144_1_W.csv',verbose = FALSE) 

