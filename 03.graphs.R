################################################################################
# Updated version of the R code for the analysis in:
#
#   "Attributable risk from distributed lag models"
#   Antonio Gasparrini & Michela Leone
#   BMC Medical Research Methodology 2014
#   http://www.ag-myresearch.com/bmcmrm2014.html
#
# Update: 17 March 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

################################################################################
# FIGURE 1

fcurve <- function(x) dchisq(x,4)*7

plot(0:7,seq(-0.5,1.5,length=8),type="n",frame.plot=F,ylim=c(-0.5,1.5),
  axes=F,xlab="Time (Lags)",ylab="Effect",xlim=c(-0.5,7.5),mgp=c(2.5,1,0))
#rect(1.5,-0.5,2.5,1.5,col=grey(0.9),border=NA)
abline(v=1,lty=2)
axis(1,at=0:8-0.5,labels=F)
axis(1,at=1:6,tick=F,labels=c(expression(bold(t)),expression(t+1),
  expression(t+2),expression(...),expression(...),expression(t+L)))
axis(2,at=-1:3*0.5,labels=c(NA,0,NA,NA,NA))
lines(0:50/10+1,fcurve(0:50/10*2.5+0.5),lwd=1.5)
points(1:6,fcurve(0:5*2.5+0.5),pch=21,bg="green3",cex=1.6)
abline(h=0)
title("Forward perspective")

plot(0:7,seq(-0.5,1.5,length=8),type="n",frame.plot=F,ylim=c(-0.5,1.5),
  axes=F,xlab="Time (Lags)",ylab="Effect",xlim=c(-0.5,7.5),mgp=c(2.5,1,0))
#rect(1.5,-0.5,2.5,1.5,col=grey(0.9),border=NA)
abline(v=6,lty=2)
axis(1,at=0:8-0.5,labels=F)
axis(1,at=6:1,tick=F,labels=c(expression(bold(t)),expression(t-1),
  expression(t-2),expression(...),expression(...),expression(t-L)))
axis(2,at=-1:3*0.5,labels=c(NA,0,NA,NA,NA))
for(i in 0:4*10) lines(i:50/10+1,fcurve(0:(50-i)/10*2.5+0.5),lwd=1.5)
points(rep(6,5),fcurve(1:5*2.5+0.5),pch=22,bg="yellow",cex=1.6)
abline(h=0)
title("Backward perspective")

################################################################################
# FIGURE 2

plot(cp1,xlab="Temperature (C)",zlim=c(0.85,1.45),xlim=c(-5,31),
  zlab="RR",ltheta=170,phi=35,lphi=30)
title("Exposure-lag-response surface")
mtext("London 1993-2006",cex=0.75)

plot(cp1cold,"overall",col=4,ylab="RR",xlab="Temperature",xlim=c(-5,35),
  ylim=c(-0.5,3.5),axes=F,lwd=1.5)
lines(cp1hot,"overall",ci="area",col=2,lwd=1.5)
axis(1,at=-1:7*5)
axis(2,at=c(1:7*0.5))
title("Overall cumulative association and temperature distribution")
mtext("London 1993-2006",cex=0.75)
par(new=T)
hist(lndn$tmean,xlim=c(-5,35),ylim=c(0,1200),axes=F,ann=F,col=grey(0.95),breaks=30)
abline(v=quantile(lndn$tmean,c(0.01,0.99)),lty=2)
abline(v=cen,lty=3)
axis(4,at=0:4*100)
mtext("Freq",4,line=2.5,at=200,cex=0.8)


################################################################################
# FIGURE 3 (NB: DIFFERENT THAN ORIGINAL AS ROME DATA ARE NOT AVAILABLE)

plot(cp3,var=29,ci="n",col="red3",ylab="RR",ylim=c(.95,1.2),xlab="Lag",lwd=2)
lines(cp3,var=27,col="red",lty=2,lwd=2)
lines(cp3,var=25,col="orange2",lty=4,lwd=2)
title("Lag-response association at different temperatures")
mtext("London 1993-2006",cex=0.7)
legend("topright",c(expression(paste(25*degree,"C",sep="")),
  expression(paste(27*degree,"C",sep="")),
  expression(paste(29*degree,"C",sep=""))),lty=c(4,2,1),lwd=2,
  col=c("orange2","red","red3"),inset=0.05)

attrback <- attrdl(lndn$tmean,cb,lndn$death,model,tot=F,type="an",cen=cen,
  range=c(cen,100))
attrforw <- attrdl(lndn$tmean,cb,lndn$death,model,tot=F,type="an",dir="forw",
  cen=cen,range=c(cen,100))

sub <- lndn$year==1995 & lndn$month%in%7:9
plot(as.Date(lndn$date[sub]),attrback[sub],ylim=c(-50,100),ylab="",xlab="Date",
  axes=F,frame.plot=F,cex=0.8,pch=22,bg="yellow")
points(as.Date(lndn$date[sub]),attrforw[sub],cex=0.8,pch=21,bg="green3")
axis(1,at=as.Date(paste(1995,c(7:9,9),c(1,1,30),sep="-")),labels=F)
axis(1,at=as.Date(paste(1995,7:9,15,sep="-")),tick=F,
  labels=c("July","Aug","Sept"))
axis(2,at=-1:20*5)
abline(h=0)
title("Daily attributable deaths and temperature")
mtext("London July-Sept 1995",cex=0.7)
mtext("Attributable deaths",2,line=2.5,at=40)
legend("topright",c("Forward","Backward"),pch=21:22,
  pt.bg=c("green3","yellow"),inset=0.05)
par(new=T)
plot(as.Date(lndn$date[sub]),lndn$tmean[sub],type="l",ann=F,axes=F,ylim=c(10,90))
abline(h=cen,lty=3)
axis(4,at=3:6*5)
mtext("Degree C",4,line=2.5,at=22.5,cex=0.8)

#
