suppressPackageStartupMessages(library(vegan))
o<-'review_response/03_analysis/Prakash_original_workflow_review_2026-10-27'
b<-read.csv(file.path(o,'results/Original_source_focused_beta_tests.csv'));out<-list()
for(i in seq_len(nrow(b))){
 name<-b$analysis[i];ob<-readRDS(file.path(o,'results',paste0(name,'.rds')));z<-ob$manifest;P<-ob$permutations
 stopifnot(nrow(P)==9999,ncol(P)==nrow(z),!anyDuplicated(z$SourceCode))
 stopifnot(all(apply(P,1,function(ii)identical(sort(ii),seq_len(nrow(z)))&&all(z$Participant[ii]==z$Participant))))
 if(startsWith(b$context[i],'PRE_site'))stopifnot(all(z$Fraction=='PRE'),!anyDuplicated(paste(z$Participant,z$Location))) else{
  stopifnot(all(table(z$SpecimenKey)==4),all(apply(P,1,function(ii)all(z$SpecimenKey[ii]==z$SpecimenKey))))
  for(person in levels(z$Participant))for(f in levels(z$Fraction)){
   ix<-which(z$Participant==person&z$Fraction==f);stopifnot(all(apply(P[,ix,drop=FALSE],1,function(ii)length(unique(z$Fraction[ii]))==1)))
  }
 }
 write.csv(z,file.path(o,'results',paste0(name,'_manifest.csv')),row.names=FALSE)
 out[[name]]<-data.frame(analysis=name,all_permutation_rows_checked=nrow(P),between_woman_exchange=FALSE,source_rows=nrow(z),design_check='PASS')
}
write.csv(do.call(rbind,out),file.path(o,'results/Independent_saved_permutation_checks.csv'),row.names=FALSE)
# Requested joint factorial effect summaries, descriptive only. A single unrestricted
# row permutation within women is not silently used to test all nesting levels.
for(tag in c('baseline_39','sensitivity_41')){
 x<-as.matrix(read.csv(file.path(o,'data',paste0(tag,'_source_counts.csv')),row.names=1,check.names=FALSE));rel<-sweep(x,2,colSums(x),'/')
 z<-read.csv(file.path(o,'data/source_metadata_original_codes.csv'));for(n in c('Participant','Location','Fraction'))z[[n]]<-factor(z[[n]])
 y<-t(rel[,z$SourceCode,drop=FALSE])
 for(order in c('site_first','fraction_first')){
  a<-if(order=='site_first')adonis2(y~Location*Fraction,data=z,method='bray',permutations=0,by='terms') else adonis2(y~Fraction*Location,data=z,method='bray',permutations=0,by='terms')
  # Remove placeholder permutation p values if emitted with permutations=0.
  d<-data.frame(term=rownames(a),a,row.names=NULL,check.names=FALSE);d<-d[,!grepl('Pr\\(',names(d)),drop=FALSE]
  d$interpretation<-'Descriptive sequential partition only; order-dependent; no permutation inference or R2-ranking test.'
  write.csv(d,file.path(o,'results',paste0(tag,'_Joint_factorial_DESCRIPTIVE_',order,'.csv')),row.names=FALSE)
 }
}
cat('PASS: all saved restricted permutation rows checked; joint tables explicitly descriptive.\n')
