################################################################################
# Updated version of the R code for the analysis in:
#
#   "Attributable risk from distributed lag models"
#   Antonio Gasparrini & Michela Leone
#   BMC Medical Research Methodology 2014
#   http://www.ag-myresearch.com/bmcmrm2014.html
#
# Update: 14 March 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

library(dlnm) ; library(splines) ; library(foreign) ; library(tsModel)

# CHECK VERSION OF THE PACKAGE
  if(packageVersion("dlnm")<"2.2.0")
    stop("update dlnm package to version >= 2.2.0")

# LOAD THE DATA
lndn <- read.csv("london.csv")

################################################################################
# DERIVE THE CROSS-BASIS

# KNOTS FOR EXPOSURE-RESPONSE FUNCTION
vk <- equalknots(lndn$tmean,fun="bs",df=4,degree=2)

# KNOTS FOR THE LAG-RESPONSE FUNCTION
maxlag <- 25
lk <- logknots(maxlag,3)

# CENTERING VALUE (AND PERCENTILE)
cen <- 20
sum(lndn$tmean<cen,na.rm=T)/sum(!is.na(lndn$tmean))

# COMPUTE THE CROSS-BASIS
cb <- crossbasis(lndn$tmean, lag=maxlag, argvar=list(fun="bs",degree=2,
  knots=vk), arglag=list(knots=lk))

# SUMMARY
summary(cb)

# COMPUTE 1ST AND 99TH PERCENTILES
(perc <- quantile(lndn$tmean,c(0.01,0.99)))

################################################################################
# RUN THE MODEL AND OBTAIN PREDICTIONS

# RUN THE MODEL
model <- glm(death~cb+ns(time,10*14)+dow,family=quasipoisson(),lndn)

# PREDICTION FOR 3D AND OVERALL CUMULATIVE
cp1 <- crosspred(cb,model,at=seq(min(lndn$tmean,na.rm=T),
  max(lndn$tmean,na.rm=T),length=25),cen=cen)
cp1cold <- crosspred(cb,model,at=seq(min(lndn$tmean,na.rm=T),
  cen,length=20),cen=cen)
cp1hot <- crosspred(cb,model,at=seq(cen,max(lndn$tmean,na.rm=T),
  length=20),cen=cen)

# PREDICTIONS AT 1ST AND 99TH PERCENTILES
cp2 <- crosspred(cb,model,at=perc,bylag=0.2,cen=cen)

# PREDICTION FOR LAG-RESPONSE CURVES
cp3 <- crosspred(cb,model,at=c(25,27,29),bylag=0.2,cen=cen)

#
