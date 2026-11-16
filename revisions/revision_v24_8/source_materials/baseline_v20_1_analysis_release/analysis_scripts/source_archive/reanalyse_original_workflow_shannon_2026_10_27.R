# New exploratory review models; original tables/labels are never modified.
suppressPackageStartupMessages({library(lmerTest);library(ggplot2)})
o <- 'review_response/03_analysis/Prakash_original_workflow_review_2026-10-27'
m <- read.csv(file.path(o,'data/source_metadata_original_codes.csv'),stringsAsFactors=FALSE)
put <- function(d,n) write.csv(d,file.path(o,'results',paste0(n,'.csv')),row.names=FALSE)
alltests <- list();diagnostics <- list();vals <- list()
for(tag in c('baseline_39','sensitivity_41')) {
 x <- as.matrix(read.csv(file.path(o,'data',paste0(tag,'_source_counts.csv')),row.names=1,check.names=FALSE))
 stopifnot(identical(colnames(x),m$SourceCode),all(colSums(x)>0),all(x>=0))
 p <- sweep(x,2,colSums(x),'/');lp <- matrix(0,nrow(p),ncol(p));lp[p>0]<-log(p[p>0])
 z <- m;z$Shannon <- -colSums(p*lp)
 z$Participant<-factor(z$Participant);z$SpecimenKey<-factor(z$SpecimenKey)
 z$Location<-factor(z$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'))
 z$Fraction<-factor(z$Fraction,levels=c('PRE','IgA','IgM','IgG'))
 z$Phase<-factor(z$Phase,levels=c('Proliferative','Secretory'))
 put(z,paste0(tag,'_Shannon_per_source_library'));vals[[tag]]<-z
 for(whichfit in c('all149_nested','all149_nested_phase_adjusted','proliferative_nested','omit99_unsuffixed_endo_nested','omit99_suffix2_endo_nested')) {
  d <- z
  if(whichfit=='proliferative_nested')d<-droplevels(z[z$Phase=='Proliferative',])
  if(whichfit=='omit99_unsuffixed_endo_nested')d<-droplevels(z[!z$SourceCode%in%c('99EA','99EG','99EM','99EP'),])
  if(whichfit=='omit99_suffix2_endo_nested')d<-droplevels(z[!z$SourceCode%in%c('99EA2','99EG2','99EM2','99EP2'),])
  form <- Shannon~Location*Fraction+(1|Participant)+(1|SpecimenKey)
  if(whichfit=='all149_nested_phase_adjusted')form<-update(form,.~.+Phase)
  fit <- lmer(form,d,REML=TRUE)
  name<-paste(tag,whichfit,sep='_');saveRDS(fit,file.path(o,'results',paste0(name,'.rds')))
  a<-anova(fit,ddf='Satterthwaite');a<-data.frame(term=rownames(a),a,row.names=NULL,check.names=FALSE)
  a$model<-name;a$libraries<-nrow(d);a$women<-nlevels(droplevels(d$Participant));alltests[[name]]<-a
  put(data.frame(term=rownames(summary(fit)$coefficients),summary(fit)$coefficients,row.names=NULL,check.names=FALSE),paste0(name,'_coefficients'))
  put(as.data.frame(VarCorr(fit)),paste0(name,'_variance_components'))
  diagnostics[[name]]<-data.frame(model=name,libraries=nrow(d),women=nlevels(droplevels(d$Participant)),specimens=nlevels(droplevels(d$SpecimenKey)),singular=isSingular(fit),convergence=paste(fit@optinfo$conv$lme4$messages,collapse=';'),residual_SD=sigma(fit),minimum_fitted=min(fitted(fit)),maximum_fitted=max(fitted(fit)))
  pdf(file.path(o,'results',paste0(name,'_diagnostics.pdf')),width=8,height=8)
  par(mfrow=c(2,2));plot(fitted(fit),residuals(fit),xlab='Fitted Shannon',ylab='Residual');abline(h=0,lty=2)
  qqnorm(residuals(fit));qqline(residuals(fit));hist(residuals(fit),main='Residual distribution',xlab='Residual')
  plot(d$Shannon,fitted(fit),xlab='Observed Shannon',ylab='Fitted Shannon');abline(0,1,lty=2);dev.off()
 }
}
put(do.call(rbind,alltests),'Exploratory_Shannon_tests_all_specifications')
put(do.call(rbind,diagnostics),'Model_diagnostics_register')
z<-vals[['baseline_39']]
g<-ggplot(z,aes(Fraction,Shannon,group=SpecimenKey))+geom_line(alpha=.3,color='grey40')+geom_point(aes(color=Phase),size=1.6)+facet_wrap(~Location,nrow=1)+theme_bw(base_size=11)+labs(x='Library fraction',y="Shannon diversity (39-taxon baseline)",caption=sprintf('%d libraries; %d women; %d workbook preparation keys. Both participant-99 endometrial preparations retained.\nExploratory report-linked revision analysis; not a replacement of the submitted figure.',nrow(z),length(unique(z$Participant)),length(unique(z$SpecimenKey))))
ggsave(file.path(o,'results/Shannon_original_source_baseline39_REVIEW.png'),g,width=11,height=4.5,dpi=300)
ggsave(file.path(o,'results/Shannon_original_source_baseline39_REVIEW.pdf'),g,width=11,height=4.5)
capture.output(sessionInfo(),file=file.path(o,'R_sessionInfo.txt'))
print(do.call(rbind,alltests)[,c('model','term','Pr(>F)')])
