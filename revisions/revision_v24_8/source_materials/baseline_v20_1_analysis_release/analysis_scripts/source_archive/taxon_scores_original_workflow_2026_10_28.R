# Matched source-preparation scores with equal woman weighting.
suppressPackageStartupMessages(library(ggplot2))
o<-'review_response/03_analysis/Taxon_scores_PRE_review_2026-10-28'
m<-read.csv(file.path(o,'data/source_metadata_with_vaginal_category.csv'),stringsAsFactors=FALSE)
tx<-read.csv(file.path(o,'data/original_taxonomy.csv'),stringsAsFactors=FALSE)
put<-function(x,n)write.csv(x,file.path(o,'results',paste0(n,'.csv')),row.names=FALSE)
B<-5000L;allout<-list();poolout<-list();womenout<-list();bootout<-list();pairout<-list();missing<-list()
pairkeys<-unique(m$SpecimenKey)
for(sp in pairkeys){
 a<-m[m$SpecimenKey==sp&m$Fraction=='IgA',];b<-m[m$SpecimenKey==sp&m$Fraction=='IgM',]
 if(nrow(a)==1&&nrow(b)==1)pairout[[sp]]<-data.frame(SpecimenKey=sp,Participant=a$Participant,Location=a$Location,Phase=a$Phase,IgA_SourceCode=a$SourceCode,IgM_SourceCode=b$SourceCode)
 else missing[[sp]]<-data.frame(SpecimenKey=sp,Participant=m$Participant[match(sp,m$SpecimenKey)],Location=m$Location[match(sp,m$SpecimenKey)],IgA_available=nrow(a)==1,IgM_available=nrow(b)==1,reason='Both IgA and IgM source libraries required; missing fractions not imputed.')
}
pairs<-do.call(rbind,pairout);put(pairs,'Matched_IgA_IgM_preparation_manifest');put(do.call(rbind,missing),'Preparations_without_matched_IgA_IgM')
for(arm in c('baseline_39','sensitivity_41')){
 x<-as.matrix(read.csv(file.path(o,'data',paste0(arm,'_source_counts.csv')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$SourceCode))
 genus<-tx$Genus[match(rownames(x),tx$Taxon)];stopifnot(!anyNA(genus),all(nchar(genus)>0))
 for(rank in c('provisional_species','genus_aggregation')){
  mat<-if(rank=='provisional_species')x else rowsum(x,genus,reorder=TRUE);rel<-sweep(mat,2,colSums(mat),'/');taxa<-rownames(rel)
  for(scope in c('All','Proliferative')){
   pp<-if(scope=='All')pairs else pairs[pairs$Phase=='Proliferative',];allpeople<-sort(unique(pp$Participant));siteM<-list()
   for(site in c('VAGINA','CERVIX','ENDOMETRIUM')){
    z<-pp[pp$Location==site,];people<-sort(unique(z$Participant));n<-length(people);stopifnot(n>1)
    A<-t(sapply(people,function(p)rowMeans(rel[,z$IgA_SourceCode[z$Participant==p],drop=FALSE])))
    M<-t(sapply(people,function(p)rowMeans(rel[,z$IgM_SourceCode[z$Participant==p],drop=FALSE])))
    rownames(A)<-rownames(M)<-people;colnames(A)<-colnames(M)<-taxa
    # A woman with two paired preparations contributes the mean relative abundance
    # per fraction, then receives the same final weight as every other woman.
    muA<-colMeans(A);muM<-colMeans(M);den<-muA+muM;sc<-muA/den;sc[den==0]<-NA_real_
    set.seed(20261028L+match(site,c('VAGINA','CERVIX','ENDOMETRIUM'))+100L*(scope=='Proliferative'))
    W<-t(rmultinom(B,size=n,prob=rep(1/n,n)));ba<-W%*%A;bm<-W%*%M;bs<-ba/(ba+bm);bs[ba+bm==0]<-NA_real_
    ci<-apply(bs,2,function(v)if(sum(is.finite(v))>1)quantile(v[is.finite(v)],c(.025,.975),names=FALSE) else c(NA_real_,NA_real_))
    detected<-colSums(A+M>0);ws<-A/(A+M);ws[A+M==0]<-NA_real_
    alt<-apply(ws,2,function(v)if(any(is.finite(v)))mean(v[is.finite(v)]) else NA_real_)
    name<-paste(arm,rank,scope,site,sep='_');siteM[[site]]<-list(A=A,M=M)
    allout[[name]]<-data.frame(taxon_arm=arm,rank=rank,scope=scope,Location=site,Taxon=taxa,n_matched_women=n,n_paired_preparations=nrow(z),women_detected_either=detected,women_detected_IgA=colSums(A>0),women_detected_IgM=colSums(M>0),mean_IgA_relative_abundance=muA,mean_IgM_relative_abundance=muM,Score_ratio_of_matched_woman_means=sc,CI95_low=ci[1,],CI95_high=ci[2,],bootstrap_valid=colSums(is.finite(bs)),bootstrap_attempted=B,bootstrap_valid_fraction=colMeans(is.finite(bs)),mean_woman_ratio_conditional_detection=alt,descriptive_flag=ifelse(detected<3,'Fewer than 3 detected women; sparse descriptive estimate','At least 3 detected women; not an experimental validation threshold'),interpretation='Relative sequencing representation in matched IgA/IgM libraries; not per-cell coating, antibody affinity or percent of bacteria coated. Endpoint 0 or 1 does not establish exclusive binding.')
    bootout[[name]]<-list(scores=bs,draw_weights=W,women=people,taxa=taxa)
    womenout[[name]]<-do.call(rbind,lapply(seq_along(people),function(i)data.frame(taxon_arm=arm,rank=rank,scope=scope,Location=site,Participant=people[i],Taxon=taxa,IgA_relative_mean=A[i,],IgM_relative_mean=M[i,],n_preparations=sum(z$Participant==people[i]))))
   }
   # Secondary equal-site score requires all three site scores. Resample women
   # jointly across sites so within-woman site dependence is preserved.
   set.seed(20261031L+100L*(scope=='Proliferative'));W<-t(rmultinom(B,size=length(allpeople),prob=rep(1/length(allpeople),length(allpeople))))
   arr<-array(NA_real_,dim=c(B,length(taxa),3));points<-matrix(NA_real_,length(taxa),3)
   for(j in seq_along(siteM)){
    s<-siteM[[j]];A<-M<-matrix(0,length(allpeople),length(taxa));ix<-match(rownames(s$A),allpeople);A[ix,]<-s$A;M[ix,]<-s$M
    points[,j]<-colSums(A)/(colSums(A)+colSums(M));sa<-W%*%A;sm<-W%*%M;arr[,,j]<-sa/(sa+sm)
   }
   bs<-apply(arr,c(1,2),function(v)if(all(is.finite(v)))mean(v) else NA_real_)
   pt<-apply(points,1,function(v)if(all(is.finite(v)))mean(v) else NA_real_)
   ci<-apply(bs,2,function(v)if(sum(is.finite(v))>1)quantile(v[is.finite(v)],c(.025,.975),names=FALSE) else c(NA_real_,NA_real_))
   name<-paste(arm,rank,scope,'Equal_three_site_secondary',sep='_')
   poolout[[name]]<-data.frame(taxon_arm=arm,rank=rank,scope=scope,Taxon=taxa,n_sites_estimable=rowSums(is.finite(points)),equal_three_site_secondary_score=pt,CI95_low=ci[1,],CI95_high=ci[2,],bootstrap_valid=colSums(is.finite(bs)),bootstrap_attempted=B,interpretation='Secondary equal-site mean; undefined unless all three site scores estimable. Joint woman bootstrap; pointwise conditional intervals, not independent validation.')
   bootout[[name]]<-list(scores=bs,draw_weights=W,women=allpeople,taxa=taxa)
  }
 }
}
put(do.call(rbind,allout),'Site_specific_matched_IgA_IgM_representation_scores');put(do.call(rbind,poolout),'Equal_three_site_secondary_scores');put(do.call(rbind,womenout),'Woman_site_IgA_IgM_relative_abundances')
saveRDS(bootout,file.path(o,'results/Score_bootstrap_draws_and_weights.rds'))
res<-do.call(rbind,allout);key<-c('Sneathia vaginalis','Gardnerella vaginalis','Lactobacillus iners','Lactobacillus crispatus','Lactobacillus jensenii','Lactobacillus gasseri','Prevotella buccalis','Fannyhessea vaginae')
z<-res[res$taxon_arm=='baseline_39'&res$rank=='provisional_species'&res$scope=='All'&res$Taxon%in%key,];z$Taxon<-factor(z$Taxon,levels=rev(key));z$Location<-factor(z$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'))
g<-ggplot(z,aes(Score_ratio_of_matched_woman_means,Taxon))+geom_vline(xintercept=.5,lty=2,color='grey60')+geom_errorbar(aes(xmin=CI95_low,xmax=CI95_high),orientation='y',width=.2)+geom_point(aes(size=women_detected_either),color='#2166AC')+facet_wrap(~Location,nrow=1)+scale_x_continuous(limits=c(0,1))+theme_bw(base_size=11)+labs(x='Matched IgA / (IgA + IgM) representation score',y=NULL,size='Detected women',caption='Provisional taxonomic assignments; 95% pointwise woman-bootstrap intervals. Sparse and undefined scores are flagged in tables.\nScores describe relative library representation, not antibody coating per bacterial cell.')
ggsave(file.path(o,'results/Matched_site_scores_key_taxa_REVIEW.png'),g,width=12,height=6,dpi=300);ggsave(file.path(o,'results/Matched_site_scores_key_taxa_REVIEW.pdf'),g,width=12,height=6)
capture.output(sessionInfo(),file=file.path(o,'R_score_sessionInfo.txt'));print(table(pairs$Location));cat('Site and secondary scores saved with detection and bootstrap denominators.\n')
