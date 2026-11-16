# Additional exploratory joint factorial dbRDA/PERMANOVA for Reviewer 2 Q15.
# The full model is common to all marginal contrasts. Residual permutations
# preserve women, whole anatomical units and synchronized fraction structure.
# Approximate inference assumes residual exchangeability under this hierarchy;
# known unequal dispersion prevents a centroid-only biological interpretation.
suppressPackageStartupMessages({library(vegan);library(permute)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))normalizePath(args[1]) else normalizePath('.')
o<-if(length(args)>1)args[2] else file.path(root,'review_response/Q15_joint_inference_2026-11-05')
N<-if(length(args)>2)as.integer(args[3]) else 9999L
B<-if(length(args)>3)as.integer(args[4]) else 1000L
dir.create(o,recursive=TRUE,showWarnings=FALSE)
for(d in c('inputs','results','diagnostics'))dir.create(file.path(o,d),showWarnings=FALSE)
i<-if(file.exists(file.path(root,'inputs/Metadata_original_codes.csv')))file.path(root,'inputs') else file.path(root,'review_response/Author_review_v19_2026-11-05/04_data_and_code/bundle/source_filter_inputs')
for(f in c('Metadata_original_codes.csv','Counts_39_original_codes.csv','Counts_41_original_codes.csv')){
 stopifnot(file.copy(file.path(i,f),file.path(o,'inputs',f),overwrite=FALSE))
}
m<-read.csv(file.path(i,'Metadata_original_codes.csv'),stringsAsFactors=FALSE,colClasses=c(Participant='character',SourceSuffix='character'))
stopifnot(nrow(m)==149L,length(unique(m$Participant))==15L,length(unique(m$SpecimenKey))==40L)
spl<-split(m,m$SpecimenKey)
complete<-names(Filter(function(z)nrow(z)==4L && setequal(z$Fraction,c('PRE','IgA','IgM','IgG')),spl))
mc<-m[m$SpecimenKey%in%complete,]
stopifnot(nrow(mc)==124L,length(unique(mc$Participant))==14L,length(unique(mc$SpecimenKey))==31L)
write.csv(mc,file.path(o,'results/Complete_preparation_124_library_manifest.csv'),row.names=FALSE)
write.csv(m[!m$SpecimenKey%in%complete,],file.path(o,'results/Excluded_incomplete_preparation_libraries.csv'),row.names=FALSE)
options(contrasts=c('contr.sum','contr.poly'))
fractions<-c('PRE','IgA','IgM','IgG')
counts<-list()
for(n in c(39,41)){
 x<-as.matrix(read.csv(file.path(i,paste0('Counts_',n,'_original_codes.csv')),row.names=1,check.names=FALSE))
 stopifnot(nrow(x)==n,identical(colnames(x),m$SourceCode),all(colSums(x)>0))
 counts[[as.character(n)]]<-sweep(x,2,colSums(x),'/')
}
profile_data<-function(n,choice){
 q<-mc
 if(choice=='keep99_unsuffixed')q<-q[!(q$Participant=='99'&q$Location=='ENDOMETRIUM'&q$SourceSuffix=='2'),]
 if(choice=='keep99_suffix2')q<-q[!(q$Participant=='99'&q$Location=='ENDOMETRIUM'&q$SourceSuffix!='2'),]
 q$AnatomicalUnit<-paste(q$Participant,q$Location,sep='_')
 key<-paste(q$AnatomicalUnit,q$Fraction,sep='_')
 by<-split(seq_len(nrow(q)),key)
 y<-t(vapply(by,function(ix)rowMeans(counts[[as.character(n)]][,q$SourceCode[ix],drop=FALSE]),numeric(n)))
 zz<-do.call(rbind,lapply(names(by),function(k){
  z<-q[by[[k]],];stopifnot(length(unique(z$Participant))==1L,length(unique(z$Location))==1L,length(unique(z$Fraction))==1L)
  data.frame(ProfileID=k,Participant=z$Participant[1],Location=z$Location[1],Fraction=z$Fraction[1],AnatomicalUnit=z$AnatomicalUnit[1],source_preparations=length(unique(z$SpecimenKey)),source_libraries=paste(z$SourceCode,collapse=';'),stringsAsFactors=FALSE)
 }))
 rownames(zz)<-zz$ProfileID
 zz$Participant<-factor(zz$Participant);zz$Location<-factor(zz$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'))
 zz$Fraction<-factor(zz$Fraction,levels=fractions);zz$AnatomicalUnit<-factor(zz$AnatomicalUnit)
 stopifnot(nrow(y)==120L,nrow(zz)==120L,nlevels(zz$Participant)==14L,nlevels(zz$AnatomicalUnit)==30L,max(abs(rowSums(y)-1))<1e-12,identical(rownames(y),zz$ProfileID),all(table(zz$AnatomicalUnit,zz$Fraction)==1L))
 list(y=y,z=zz,underlying_libraries=nrow(q),underlying_preparations=length(unique(q$SpecimenKey)))
}
base<-profile_data(39,'mean_both')
z<-base$z;n<-nrow(z)
set.seed(20261105L)
P<-matrix(NA_integer_,N,n)
for(b in seq_len(N))for(person in levels(z$Participant)){
 units<-unique(as.character(z$AnatomicalUnit[z$Participant==person]))
 shuffled_units<-sample(units,length(units));shuffled_fractions<-sample(fractions)
 for(j in seq_along(units)){
  target<-which(z$AnatomicalUnit==units[j]);source<-which(z$AnatomicalUnit==shuffled_units[j])
  P[b,target]<-source[match(shuffled_fractions[match(z$Fraction[target],fractions)],z$Fraction[source])]
 }
}
structure_checks<-list()
relation<-list(woman=outer(as.character(z$Participant),as.character(z$Participant),'=='),anatomical_unit=outer(as.character(z$AnatomicalUnit),as.character(z$AnatomicalUnit),'=='),woman_fraction=outer(paste(z$Participant,z$Fraction),paste(z$Participant,z$Fraction),'=='))
for(k in names(relation)){
 rel<-relation[[k]];pass<-vapply(seq_len(N),function(b)identical(rel[P[b,],P[b,]],rel),TRUE)
 stopifnot(all(pass));structure_checks[[k]]<-data.frame(relation=k,permutations_checked=N,violations=sum(!pass))
}
stopifnot(all(vapply(seq_len(N),function(b)identical(sort(P[b,]),seq_len(n))&&all(z$Participant[P[b,]]==z$Participant),TRUE)))
write.csv(P,file.path(o,'diagnostics/Hierarchical_permutation_indices_1based.csv'),row.names=FALSE)
saveRDS(P,file.path(o,'diagnostics/Hierarchical_permutations.rds'))
write.csv(do.call(rbind,structure_checks),file.path(o,'diagnostics/Permutation_relation_invariants.csv'),row.names=FALSE)

# Bootstrap independent women, relabelling repeated draws as separate clusters.
set.seed(20261106L);people<-levels(z$Participant)
draws<-matrix(sample(seq_along(people),B*length(people),replace=TRUE),B,length(people))
write.csv(draws,file.path(o,'diagnostics/Woman_bootstrap_indices_1based.csv'),row.names=FALSE)
write.csv(data.frame(index=seq_along(people),Participant=people),file.path(o,'diagnostics/Bootstrap_woman_key.csv'),row.names=FALSE)
H<-function(x){q<-qr(x);qq<-qr.Q(q)[,seq_len(q$rank),drop=FALSE];tcrossprod(qq)}
inertia<-function(G,H)sum(G*H)
matrix_parts<-function(zz){
 X<-model.matrix(~Location*Fraction,zz);a<-attr(X,'assign')
 Z<-model.matrix(~Participant,zz)
 full<-cbind(Z,X[,a>0,drop=FALSE]);hf<-H(full)
 list(X=X,assign=a,Z=Z,full=full,Hfull=hf,rank=qr(full)$rank)
}
table_out<-list();bootstrap_out<-list();dispersion_out<-list();all_null<-list();fit_summaries<-list();k<-0
arms<-data.frame(tag=c('primary39_mean_both','sensitivity41_mean_both','sensitivity39_keep99_unsuffixed','sensitivity39_keep99_suffix2'),taxa=c(39,41,39,39),choice=c('mean_both','mean_both','keep99_unsuffixed','keep99_suffix2'),stringsAsFactors=FALSE)
for(ii in seq_len(nrow(arms))){
 arm<-arms[ii,];a<-profile_data(arm$taxa,arm$choice);y<-a$y;zz<-a$z
 stopifnot(identical(as.character(zz$ProfileID),as.character(z$ProfileID)))
 tag<-arm$tag;message('Running ',tag,'; ',N,' residual permutations per marginal contrast')
 write.csv(zz,file.path(o,'results',paste0(tag,'_profile_manifest.csv')),row.names=FALSE)
 write.csv(data.frame(ProfileID=rownames(y),y,check.names=FALSE),file.path(o,'results',paste0(tag,'_relative_abundance_profiles.csv')),row.names=FALSE)
 D<-as.matrix(vegdist(y,'bray'));J<-diag(n)-matrix(1/n,n,n);G<- -0.5*J%*%(D^2)%*%J
 write.csv(D,file.path(o,'diagnostics',paste0(tag,'_BrayCurtis.csv')),row.names=FALSE)
 parts<-matrix_parts(zz);stopifnot(parts$rank==25L)
 write.csv(parts$full,file.path(o,'diagnostics',paste0(tag,'_full_design.csv')),row.names=FALSE)
 fit_joint<-dbrda(y~Location*Fraction+Condition(Participant),zz,distance='bray')
 stopifnot(abs(fit_joint$tot.chi-sum(diag(G)))<1e-9)
 fullSS<-inertia(G,parts$Hfull-parts$Z%*%solve(crossprod(parts$Z),t(parts$Z)))
 # The partial-model full inertia must agree with the common matrix projection.
 stopifnot(abs(fullSS-fit_joint$CCA$tot.chi)<1e-8)
 marginal<-list();fits<-list()
 for(term in seq_len(3)){
  termname<-c('Site','Fraction','Site_by_fraction')[term]
  Xt<-parts$X[,parts$assign==term,drop=FALSE]
  Zn<-cbind(parts$Z[,-1,drop=FALSE],parts$X[,parts$assign>0 & parts$assign!=term,drop=FALSE])
  fit<-dbrda(y~Xt+Condition(Zn),distance='bray')
  red<-cbind(1,Zn);hr<-H(red);ht<-parts$Hfull-hr
  ss<-inertia(G,ht);resSS<-inertia(G,diag(n)-parts$Hfull)
  df<-c(2,3,6)[term];rdf<-n-parts$rank
  stopifnot(abs(ss-fit$CCA$tot.chi)<1e-8,abs(resSS-fit$CA$tot.chi)<1e-8,qr(red)$rank==parts$rank-df,rdf==95L)
  F<-ss/df/(resSS/rdf)
  # Explicit fixed-nuisance Freedman-Lane engine. Re-addition of fitted
  # nuisance values cancels in both orthogonal projection statistics.
  # Do not substitute vegan's partial-model permutation implementation:
  # its nuisance handling need not equal this fixed-design residual scheme.
  Rred<-diag(n)-hr;Rfull<-diag(n)-parts$Hfull
  Gres<-Rred%*%G%*%Rred
  fnull<-vapply(seq_len(N),function(b){
   gp<-Gres[P[b,],P[b,]]
   inertia(gp,ht)/df/(inertia(gp,Rfull)/rdf)
  },0.)
  stopifnot(all(is.finite(fnull)))
  p<-(1+sum(fnull>=F-sqrt(.Machine$double.eps)))/(N+1)
  test<-list(F.0=F,F.perm=matrix(fnull,ncol=1),df=c(df,rdf),model='Explicit fixed-nuisance Freedman-Lane residual permutations; custom Gower projection implementation')
  write.csv(red,file.path(o,'diagnostics',paste0(tag,'_',termname,'_reduced_design.csv')),row.names=FALSE)
  all_null[[paste(tag,termname,sep='__')]]<-as.numeric(test$F.perm[,1])
  marginal[[term]]<-data.frame(analysis=tag,taxon_membership=arm$taxa,HUT099_handling=arm$choice,term=termname,df=df,residual_df=rdf,marginal_SS=ss,pseudo_F=F,term_R2_total=ss/fit$tot.chi,partial_R2=ss/(ss+resSS),p_raw=p,n_profiles=n,n_women=14L,n_anatomical_units=30L,underlying_libraries=a$underlying_libraries,underlying_preparations=a$underlying_preparations,permutations=N,permutation_seed=20261105L,stringsAsFactors=FALSE)
  fits[[termname]]<-list(model=fit,test=test)
 }
 mm<-do.call(rbind,marginal);mm$p_Holm_three_terms<-p.adjust(mm$p_raw,'holm');mm$p_BH_three_terms<-p.adjust(mm$p_raw,'BH')
 # Primary family: three term tests, Holm. BH is supplied transparently as a
 # secondary correction only, not selected after seeing which is significant.
 table_out[[tag]]<-mm
 saveRDS(list(joint=fit_joint,marginal=fits,manifest=zz),file.path(o,'diagnostics',paste0(tag,'_models_and_tests.rds')))
 eig<-eigen(G,symmetric=TRUE,only.values=TRUE)$values
 fit_summaries[[tag]]<-data.frame(analysis=tag,total_SS=fit_joint$tot.chi,participant_SS=fit_joint$pCCA$tot.chi,joint_effect_SS=fit_joint$CCA$tot.chi,joint_effect_df=fit_joint$CCA$qrank,residual_SS=fit_joint$CA$tot.chi,residual_df=n-parts$rank,joint_model_R2_total=fit_joint$CCA$tot.chi/fit_joint$tot.chi,negative_eigenvalues=sum(eig< -1e-10),negative_eigenvalue_mass=sum(abs(eig[eig<0]))/sum(abs(eig)),R2_definition='Each marginal term SS / total profile-level inertia; marginal term R2 values need not sum to model R2',stringsAsFactors=FALSE)

 bvals<-matrix(NA_real_,B,3);reason<-rep('',B)
 for(b in seq_len(B)){
  ix<-integer();newperson<-character()
  for(j in seq_along(people)){
   ids<-which(zz$Participant==people[draws[b,j]])
   ix<-c(ix,ids);newperson<-c(newperson,rep(paste0('draw',j),length(ids)))
  }
  zb<-zz[ix,];zb$Participant<-factor(newperson)
  if(nlevels(droplevels(zb$Location))!=3L){reason[b]<-'Missing site level in cluster draw';next}
  pp<-matrix_parts(zb)
  if(pp$rank!=nlevels(zb$Participant)+11L){reason[b]<-'Factorial design rank loss in cluster draw';next}
  nn<-length(ix);jj<-diag(nn)-matrix(1/nn,nn,nn)
  gg<- -0.5*jj%*%(D[ix,ix]^2)%*%jj;total<-sum(diag(gg))
  if(!is.finite(total)||total<=0){reason[b]<-'Invalid total inertia';next}
  for(t in seq_len(3)){
   rr<-cbind(pp$Z,pp$X[,pp$assign>0 & pp$assign!=t,drop=FALSE])
   bvals[b,t]<-inertia(gg,pp$Hfull-H(rr))/total
  }
 }
 colnames(bvals)<-c('Site','Fraction','Site_by_fraction')
 write.csv(data.frame(draw=seq_len(B),bvals,exclusion_reason=reason),file.path(o,'results',paste0(tag,'_R2_woman_bootstrap_draws.csv')),row.names=FALSE)
 for(t in seq_len(3)){
  good<-is.finite(bvals[,t]);stopifnot(sum(good)>B*.8)
  ci<-quantile(bvals[good,t],c(.025,.975),names=FALSE,type=7)
  mm$R2_CI95_low[t]<-ci[1];mm$R2_CI95_high[t]<-ci[2];mm$valid_bootstrap_draws[t]<-sum(good);mm$attempted_bootstrap_draws[t]<-B
 }
 table_out[[tag]]<-mm
 for(group in c('Location','Fraction')){
  bd<-betadisper(as.dist(D),zz[[group]],type='median',bias.adjust=TRUE)
  dt<-permutest(bd,permutations=P)
  dispersion_out[[paste(tag,group,sep='_')]]<-data.frame(analysis=tag,grouping=group,df=dt$tab$Df[1],residual_df=dt$tab$Df[2],F=dt$tab$F[1],p_raw=dt$tab$`Pr(>F)`[1],method='Median-based betadisper with small-group bias adjustment; same hierarchical restrictions',stringsAsFactors=FALSE)
  write.csv(data.frame(ProfileID=zz$ProfileID,group=zz[[group]],distance_to_spatial_median=bd$distances),file.path(o,'diagnostics',paste0(tag,'_',group,'_dispersion_distances.csv')),row.names=FALSE)
 }
 message('Finished ',tag)
}
tab<-do.call(rbind,table_out);rownames(tab)<-NULL
write.csv(tab,file.path(o,'results/Joint_factorial_marginal_tests_and_intervals.csv'),row.names=FALSE)
write.csv(tab[tab$analysis=='primary39_mean_both',],file.path(o,'results/Table_S2R_8_joint_factorial_primary.csv'),row.names=FALSE)
write.csv(tab[tab$analysis!='primary39_mean_both',],file.path(o,'results/Table_S2R_9_joint_factorial_sensitivities.csv'),row.names=FALSE)
disp<-do.call(rbind,dispersion_out);rownames(disp)<-NULL
disp$p_Holm_two_dispersion_tests<-ave(disp$p_raw,disp$analysis,FUN=function(x)p.adjust(x,'holm'))
write.csv(disp,file.path(o,'results/Table_S2R_10_joint_subset_dispersion.csv'),row.names=FALSE)
write.csv(do.call(rbind,fit_summaries),file.path(o,'results/Joint_model_inertia_and_metric_diagnostics.csv'),row.names=FALSE)
write.csv(data.frame(permutation=seq_len(N),all_null,check.names=FALSE),file.path(o,'diagnostics/All_marginal_null_F_statistics.csv'),row.names=FALSE)
write.csv(data.frame(arm=arms$tag,status='Completed',assumption='Hierarchically exchangeable reduced-model residuals; approximate Freedman-Lane inference, not an exact randomization of an assigned treatment',profile_rule='Normalize each complete source library then equally average repeated preparations within woman/site/fraction; 120 profiles from up to 124 libraries',interpretation='Observational community association; significant dispersion prevents a centroid-only or antibody-causal interpretation',stringsAsFactors=FALSE),file.path(o,'results/Analysis_scope_and_assumptions.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(o,'R_sessionInfo.txt'))
print(tab[,c('analysis','term','pseudo_F','term_R2_total','p_raw','p_Holm_three_terms','R2_CI95_low','R2_CI95_high')])
print(disp)
cat('COMPLETED: supplementary joint factorial analysis; see explicit assumptions and independent verification before manuscript use.\n')
