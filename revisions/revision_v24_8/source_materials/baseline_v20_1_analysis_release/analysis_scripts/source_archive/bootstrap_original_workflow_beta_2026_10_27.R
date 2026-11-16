# Woman-cluster percentile intervals for each separately defined focused model.
suppressPackageStartupMessages(library(vegan))
o<-'review_response/03_analysis/Prakash_original_workflow_review_2026-10-27'
b<-read.csv(file.path(o,'results/Original_source_focused_beta_tests.csv'),check.names=FALSE)
options(contrasts=c('contr.sum','contr.poly'));B<-1000L;rows<-list();draws<-list();pointchecks<-list()
fit_r2<-function(D,z,context){
 z<-droplevels(z);for(n in c('Participant','SpecimenKey','Location','Fraction'))z[[n]]<-factor(z[[n]])
 n<-nrow(D);a<-D^2;G<--.5*(a-rowMeans(a)-matrix(colMeans(a),n,n,byrow=TRUE)+mean(a))
 if(startsWith(context,'PRE_site')){cform<-~Participant;fform<-~Participant+Location}
 else if(context=='Fraction_complete'){cform<-~SpecimenKey;fform<-~SpecimenKey+Fraction}
 else{cform<-~SpecimenKey+Fraction;fform<-~SpecimenKey+Fraction+Location:Fraction}
 xc<-model.matrix(cform,z);xf<-model.matrix(fform,z);qc<-qr(xc);qf<-qr(xf)
 proj<-function(q){Q<-qr.Q(q)[,seq_len(q$rank),drop=FALSE];sum(Q*(G%*%Q))}
 total<-sum(diag(G));cs<-proj(qc);fs<-proj(qf);ss<-fs-cs
 c(R2_total=ss/total,R2_partial=ss/(total-cs),rank=qf$rank-qc$rank)
}
for(i in seq_len(nrow(b))){
 name<-b$analysis[i];obj<-readRDS(file.path(o,'results',paste0(name,'.rds')));z<-obj$manifest
 x<-as.matrix(read.csv(file.path(o,'data',paste0(b$taxon_arm[i],'_source_counts.csv')),row.names=1,check.names=FALSE));rel<-sweep(x,2,colSums(x),'/')
 D<-as.matrix(vegdist(t(rel[,z$SourceCode,drop=FALSE]),'bray'));context<-b$context[i];pt<-fit_r2(D,z,context)
 err<-max(abs(pt[1:2]-c(b$R2_total[i],b$R2_partial[i])));stopifnot(err<1e-10,pt[3]==obj$model$CCA$rank)
 pointchecks[[name]]<-data.frame(analysis=name,max_QR_vs_vegan_R2_error=err,rank=pt[3])
 people<-unique(z$Participant);set.seed(20261030L);bt<-matrix(NA_real_,B,3,dimnames=list(NULL,names(pt)))
 for(j in seq_len(B)){
  selected<-sample(people,length(people),replace=TRUE);ii<-integer();zz<-list()
  for(k in seq_along(selected)){
   ix<-which(z$Participant==selected[k]);s<-z[ix,,drop=FALSE];s$Participant<-paste0('draw',k);s$SpecimenKey<-paste0('draw',k,'_',s$SpecimenKey);ii<-c(ii,ix);zz[[k]]<-s
  }
  zz<-do.call(rbind,zz)
  bt[j,]<-tryCatch(fit_r2(D[ii,ii,drop=FALSE],zz,context),error=function(e)rep(NA_real_,3))
 }
 ok<-is.finite(bt[,1])&is.finite(bt[,2])&bt[,3]==pt[3]
 ci<-if(sum(ok)>1)apply(bt[ok,1:2,drop=FALSE],2,quantile,c(.025,.975),names=FALSE) else matrix(NA_real_,2,2)
 rows[[name]]<-data.frame(analysis=name,R2_total=pt[1],total_CI95_low=ci[1,1],total_CI95_high=ci[2,1],R2_partial=pt[2],partial_CI95_low=ci[1,2],partial_CI95_high=ci[2,2],valid_draws=sum(ok),attempted_draws=B,rank_changing_draws=sum(is.finite(bt[,3])&bt[,3]!=pt[3]),numerical_failures=sum(!is.finite(bt[,1])|!is.finite(bt[,2])),less_than_95percent_valid=sum(ok)<.95*B,interval='Pointwise 1000 woman-cluster bootstrap percentile interval; cluster copies uniquely identified; rank-changing draws excluded; exploratory, not a comparison of different-model R2 values.')
 draws[[name]]<-bt;cat(name,':',sum(ok),'/',B,'valid\n')
}
write.csv(do.call(rbind,rows),file.path(o,'results/Original_source_beta_R2_cluster_intervals.csv'),row.names=FALSE)
write.csv(do.call(rbind,pointchecks),file.path(o,'results/Independent_QR_vs_vegan_R2_checks.csv'),row.names=FALSE)
saveRDS(draws,file.path(o,'results/Original_source_beta_R2_bootstrap_draws.rds'))
