args<-commandArgs(trailingOnly=TRUE)
z<-read.csv(args[1],check.names=FALSE);stem<-args[2];nr<-as.integer(args[3]);nc<-as.integer(args[4])
xr<-as.numeric(args[5:6]);yr<-as.numeric(args[7:8]);vmax<-as.numeric(args[9]);xl<-args[10];yl<-args[11];ttl<-args[12]
pal<-hcl.colors(128,'Cividis');xe<-seq(xr[1],xr[2],length.out=151);ye<-seq(yr[1],yr[2],length.out=151)
draw<-function(){
 mm<-matrix(seq_len(nr*nc),nr,nc,byrow=TRUE);mm<-cbind(mm,rep(nr*nc+1,nr))
 layout(mm,widths=c(rep(1,nc),.38));par(oma=c(3.4,3.7,2.8,1),family='sans',cex=1)
 for(i in seq_len(nr*nc)){
  par(mar=c(2.3,2.6,1.9,.45),mgp=c(1.3,.3,0),tcl=-.18,xaxs='i',yaxs='i')
  m<-matrix(NA_real_,150,150);zz<-z[z$panel==i,];m[cbind(zz$ix,zz$iy)]<-log10(zz$events)
  image(xe,ye,m,col=pal,zlim=c(0,log10(vmax)),xlab='',ylab='',axes=FALSE,useRaster=TRUE)
  vv<-c(-1000,0,1000,100000);at<-asinh(vv/150);labs<-c('-1000','0','1000','100000')
  if((i-1)%%nc==0){ok<-at>=yr[1]&at<=yr[2];axis(2,at=at[ok],labels=labs[ok],las=1,cex.axis=.9,lwd=.5)}
  if(i>(nr-1)*nc){ok<-at>=xr[1]&at<=xr[2];axis(1,at=at[ok],labels=labs[ok],cex.axis=.9,lwd=.5)}
  box(lwd=.4,col='#64727c');mtext(unique(zz$title),3,line=.45,adj=0,cex=.9)
 }
 par(mar=c(2.3,.05,1.9,1.9),mgp=c(1.3,.3,0),xaxs='i',yaxs='i')
 sc<-seq(0,log10(vmax),length.out=128)
 image(c(0,1),sc,matrix(rep(sc,each=2),2,128),col=pal,xlab='',ylab='',axes=FALSE,useRaster=TRUE)
 br<-seq.int(0,floor(log10(vmax)));axis(4,at=br,labels=format(10^br,scientific=FALSE,trim=TRUE),las=1,cex.axis=.9,lwd=.5)
 mtext('Events/bin',4,line=1.8,cex=.9)
 mtext(xl,1,outer=TRUE,line=1.6,cex=1);mtext(yl,2,outer=TRUE,line=2.3,cex=1)
 mtext(ttl,3,outer=TRUE,line=1,adj=0,font=2,cex=1)
}
hh<-if(nr==4)6.6 else 6.3
pdf(paste0(stem,'.pdf'),width=6.15,height=hh,pointsize=9,family='Helvetica',useDingbats=FALSE);draw();dev.off()
cat('All source events remain represented in the exported histogram bins. No gates or tests applied.\n')
