# Separate exploratory questions; these tests are not one joint PERMANOVA table.
suppressPackageStartupMessages({library(vegan);library(permute);library(lmerTest);library(ggplot2)})
o<-'review_response/03_analysis/Prakash_original_workflow_review_2026-10-27'
m<-read.csv(file.path(o,'data/source_metadata_original_codes.csv'),stringsAsFactors=FALSE)
put<-function(x,n)write.csv(x,file.path(o,'results',paste0(n,'.csv')),row.names=FALSE)
tidy<-function(x)data.frame(term=rownames(x),x,check.names=FALSE,row.names=NULL)
N<-9999L;B<-5000L;options(contrasts=c('contr.sum','contr.poly'))
beta<-list();disp<-list();checks<-list();matched_models<-list();contrasts_out<-list();ord<-list()
prep<-function(z){for(v in c('Participant','SpecimenKey','Location','Fraction'))z[[v]]<-factor(z[[v]]);z$Fraction<-factor(z$Fraction,levels=c('PRE','IgA','IgM','IgG'));droplevels(z)}
fraction_perms<-function(z,seed){
 set.seed(seed);lf<-levels(z$Fraction);P<-matrix(NA_integer_,N,nrow(z))
 for(i in seq_len(N))for(person in levels(z$Participant)){
  sh<-sample(lf)
  for(sp in unique(z$SpecimenKey[z$Participant==person])){
   ix<-which(z$SpecimenKey==sp);stopifnot(length(ix)==4)
   P[i,ix]<-ix[match(sh[match(z$Fraction[ix],lf)],z$Fraction[ix])]
  }
 }
 stopifnot(all(apply(P,1,function(i)identical(sort(i),seq_len(nrow(z)))&&all(z$SpecimenKey[i]==z$SpecimenKey))))
 for(person in levels(z$Participant))for(frac in lf){ix<-which(z$Participant==person&z$Fraction==frac);stopifnot(all(apply(P[,ix,drop=FALSE],1,function(i)length(unique(z$Fraction[i]))==1)))}
 P
}
for(tag in c('baseline_39','sensitivity_41')){
 x<-as.matrix(read.csv(file.path(o,'data',paste0(tag,'_source_counts.csv')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$SourceCode))
 rel<-sweep(x,2,colSums(x),'/');D<-as.matrix(vegdist(t(rel),'bray'));put(data.frame(SourceCode=rownames(D),D,check.names=FALSE),paste0(tag,'_BrayCurtis_all149'))
 pairs<-list();missing<-list()
 for(i in which(m$Fraction!='PRE')){
  pr<-which(m$SpecimenKey==m$SpecimenKey[i]&m$Fraction=='PRE')
  if(length(pr)!=1){missing[[length(missing)+1]]<-m[i,];next}
  pairs[[length(pairs)+1]]<-data.frame(SourceCode=m$SourceCode[i],PRE=m$SourceCode[pr],Participant=m$Participant[i],SpecimenKey=m$SpecimenKey[i],Location=m$Location[i],Fraction=m$Fraction[i],Phase=m$Phase[i],Distance=D[m$SourceCode[i],m$SourceCode[pr]])
 }
 paired<-do.call(rbind,pairs);miss<-do.call(rbind,missing);put(paired,paste0(tag,'_Matched_own_PRE'));put(miss,paste0(tag,'_Missing_own_PRE'))
 # Unconstrained coordinates are descriptive; Lingoes correction affects coordinates,
 # not the untransformed Bray-Curtis distances used in tests or paired distances.
 pc<-wcmdscale(as.dist(D),k=2,eig=TRUE,add='lingoes');pco<-cbind(m,PCoA1=pc$points[,1],PCoA2=pc$points[,2]);put(pco,paste0(tag,'_PCoA_Lingoes_coordinates'))
 put(data.frame(axis=seq_along(pc$eig),eigenvalue=pc$eig),paste0(tag,'_PCoA_Lingoes_eigenvalues'))
 ord[[tag]]<-pco
 for(scope in c('All','Proliferative')){
  z<-prep(if(scope=='All')m else m[m$Phase=='Proliferative',]);complete<-names(which(table(z$SpecimenKey)==4));mc<-prep(z[z$SpecimenKey%in%complete,]);Pfrac<-fraction_perms(mc,20261027L)
  contexts<-list(Fraction_complete=mc,Interaction_complete=mc)
  pre<-z[z$Fraction=='PRE',]
  contexts$PRE_site_omit99_unsuffixed<-prep(pre[!pre$SourceCode%in%'99EP',])
  contexts$PRE_site_omit99_suffix2<-prep(pre[!pre$SourceCode%in%'99EP2',])
  for(context in names(contexts)){
   zz<-contexts[[context]];y<-t(rel[,zz$SourceCode,drop=FALSE]);ispre<-startsWith(context,'PRE_site')
   if(ispre){
    stopifnot(!anyDuplicated(paste(zz$Participant,zz$Location)))
    set.seed(20261028L);P<-shuffleSet(nrow(zz),N,control=how(blocks=zz$Participant))
    stopifnot(all(apply(P,1,function(ii)all(zz$Participant[ii]==zz$Participant))))
    fit<-dbrda(y~Location+Condition(Participant),zz,distance='bray')
   }else{
    P<-Pfrac
    if(context=='Fraction_complete')fit<-dbrda(y~Fraction+Condition(SpecimenKey),zz,distance='bray') else{
     mm<-model.matrix(~Location*Fraction,zz);Xi<-mm[,attr(mm,'assign')==3,drop=FALSE]
     fit<-dbrda(y~Xi+Condition(SpecimenKey+Fraction),zz,distance='bray');stopifnot(fit$CCA$rank==6)
    }
   }
   name<-paste(tag,scope,context,sep='_');a<-tidy(anova(fit,permutations=P,model='reduced'))
   b<-a[1,,drop=FALSE];b$analysis<-name;b$taxon_arm<-tag;b$scope<-scope;b$context<-context;b$n_libraries<-nrow(zz);b$n_women<-nlevels(zz$Participant);b$n_specimens<-nlevels(zz$SpecimenKey)
   b$R2_total<-fit$CCA$tot.chi/fit$tot.chi;b$R2_partial<-fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi);b$residual_df<-a$Df[2];b$permutations<-nrow(P)
   b$interpretation<-'Different focused nulls/subsets; R2 values are not directly ranked or statistically compared.';beta[[name]]<-b
   saveRDS(list(model=fit,permutations=P,manifest=zz),file.path(o,'results',paste0(name,'.rds')))
   checks[[name]]<-data.frame(analysis=name,participant_block_verified=TRUE,specimen_fraction_block_verified=!ispre,synchronized_within_woman=!ispre,permutations=nrow(P),unique_permutation_rows=nrow(unique(P)),minimum_group_n=min(table(if(ispre)zz$Location else zz$Fraction)))
   if(!startsWith(context,'Interaction')){
    grp<-if(ispre)zz$Location else zz$Fraction;dd<-vegdist(y,'bray');bd<-betadisper(dd,grp,type='median',bias.adjust=TRUE);dt<-tidy(permutest(bd,permutations=P)$tab)
    ddout<-dt[1,,drop=FALSE];ddout$analysis<-name;ddout$residual_df<-dt$Df[2];ddout$context<-context;ddout$taxon_arm<-tag;ddout$scope<-scope;disp[[name]]<-ddout
    put(data.frame(SourceCode=zz$SourceCode,group=grp,distance_to_median=bd$distances),paste0(name,'_dispersion_distances'))
    saveRDS(bd,file.path(o,'results',paste0(name,'_dispersion_model.rds')))
   }
  }
  pp<-prep(if(scope=='All')paired else paired[paired$Phase=='Proliferative',])
  fit<-lmer(Distance~Location*Fraction+(1|Participant)+(1|SpecimenKey),pp,REML=TRUE)
  name<-paste(tag,scope,'Matched_PRE_diagnostic',sep='_');put(tidy(anova(fit)),paste0(name,'_ANOVA'));saveRDS(fit,file.path(o,'results',paste0(name,'.rds')))
  matched_models[[name]]<-data.frame(analysis=name,n_distances=nrow(pp),women=nlevels(pp$Participant),specimens=nlevels(pp$SpecimenKey),singular=isSingular(fit),convergence=paste(fit@optinfo$conv$lme4$messages,collapse=';'),min_fitted=min(fitted(fit)),max_fitted=max(fitted(fit)),use='Exploratory Gaussian diagnostic; bounded shared-PRE response, small n and source-key assumptions require caution.')
  pdf(file.path(o,'results',paste0(name,'_residuals.pdf')),width=8,height=4);par(mfrow=c(1,2));plot(fitted(fit),residuals(fit),xlab='Fitted distance',ylab='Residual');abline(h=0,lty=2);qqnorm(residuals(fit));qqline(residuals(fit));dev.off()
  set.seed(20261029L)
  for(site in levels(pp$Location))for(k in list(c('IgA','IgM'),c('IgA','IgG'),c('IgM','IgG'))){
   aa<-pp[pp$Location==site&pp$Fraction==k[1],c('Participant','SpecimenKey','Distance')];bb<-pp[pp$Location==site&pp$Fraction==k[2],c('Participant','SpecimenKey','Distance')]
   both<-merge(aa,bb,by=c('Participant','SpecimenKey'));if(!nrow(both))next
   both$difference<-both$Distance.x-both$Distance.y
   # Average repeated preparations within each woman before resampling women.
   bywoman<-aggregate(difference~Participant,both,mean);v<-bywoman$difference
   boot<-replicate(B,mean(sample(v,length(v),replace=TRUE)));ci<-quantile(boot,c(.025,.975),names=FALSE)
   contrasts_out[[length(contrasts_out)+1]]<-data.frame(taxon_arm=tag,scope=scope,Location=site,contrast=paste(k,collapse=' minus '),women=length(v),paired_preparations=nrow(both),mean_difference_equal_weight_women=mean(v),CI95_low=ci[1],CI95_high=ci[2],interval='Pointwise 5000 woman-bootstrap percentile interval; within-woman preparation differences averaged; no simultaneous coverage claim.')
  }
 }
}
b<-do.call(rbind,beta);b$p_BH_four_focused_tests<-ave(b$`Pr(>F)`,interaction(b$taxon_arm,b$scope),FUN=function(p)p.adjust(p,'BH'));put(b,'Original_source_focused_beta_tests')
d<-do.call(rbind,disp);d$p_BH_three_dispersion_tests<-ave(d$`Pr(>F)`,interaction(d$taxon_arm,d$scope),FUN=function(p)p.adjust(p,'BH'));put(d,'Original_source_focused_dispersion_tests')
put(do.call(rbind,checks),'Original_source_permutation_design_checks');put(do.call(rbind,matched_models),'Original_source_matched_model_checks');put(do.call(rbind,contrasts_out),'Original_source_matched_PRE_woman_weighted_contrasts')
g<-ggplot(ord[['baseline_39']],aes(PCoA1,PCoA2,color=Fraction,shape=Location))+geom_point(size=2,alpha=.8)+theme_bw(base_size=11)+labs(x='PCoA axis 1 (Lingoes-corrected coordinates)',y='PCoA axis 2',caption='Descriptive ordination of 149 source-code libraries; 15 women. Tests use original Bray-Curtis distances.\nPoint separation is not a test of immune causation or anatomical ascension.')
ggsave(file.path(o,'results/PCoA_original_source_baseline39_REVIEW.png'),g,width=9,height=6,dpi=300)
p<-read.csv(file.path(o,'results/baseline_39_Matched_own_PRE.csv'));p$Fraction<-factor(p$Fraction,levels=c('IgA','IgM','IgG'))
g<-ggplot(p,aes(Fraction,Distance,group=SpecimenKey,color=Phase))+geom_line(alpha=.35)+geom_point(size=1.5)+facet_wrap(~Location,nrow=1)+theme_bw(base_size=11)+labs(x='Sorted fraction',y='Bray-Curtis distance to own PRE',caption='Matched within recorded source preparation; unmatched libraries listed separately.\nParticipant-99 preparations remain repeated observations from one woman.')
ggsave(file.path(o,'results/Matched_PRE_original_source_baseline39_REVIEW.png'),g,width=11,height=4.5,dpi=300)
capture.output(sessionInfo(),file=file.path(o,'Beta_PRE_R_sessionInfo.txt'))
print(b[,c('analysis','n_libraries','n_women','R2_total','R2_partial','Pr(>F)','p_BH_four_focused_tests')])
