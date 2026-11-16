# Fresh R --vanilla refits from frozen source-coded counts. No source mutations.
# Saved permutation draws are replayed to compare the same Monte Carlo results.
suppressPackageStartupMessages({library(lmerTest);library(vegan)})
args<-commandArgs(TRUE);root<-if(length(args))normalizePath(args[1]) else normalizePath('.')
V<-file.path(root,'bundle')
old<-file.path(root,'bundle/frozen_model_replay')
O<-file.path(root,'reproduced/fresh_model_check');dir.create(O,recursive=TRUE,showWarnings=FALSE)
put<-function(x,n)write.csv(x,file.path(O,paste0(n,'.csv')),row.names=FALSE)
m<-read.csv(file.path(V,'source_inputs/source_metadata_original_codes.csv'))
reference<-read.csv(file.path(old,'results/Exploratory_Shannon_tests_all_specifications.csv'),check.names=FALSE)
tests<-list();diagnostics<-list();j<-0
for(arm in c('baseline_39','sensitivity_41')){
 cp<-if(arm=='baseline_39')file.path(V,'source_inputs/baseline_39_source_counts.csv') else file.path(V,'taxon_analysis/data/sensitivity_41_source_counts.csv')
 X<-as.matrix(read.csv(cp,row.names=1,check.names=FALSE));stopifnot(identical(colnames(X),m$SourceCode),all(colSums(X)>0))
 rel<-sweep(X,2,colSums(X),'/');lp<-matrix(0,nrow(rel),ncol(rel));lp[rel>0]<-log(rel[rel>0]);z<-m;z$Shannon<--colSums(rel*lp)
 # vegan additionally renormalizes a matrix of relative abundances. Its tiny
 # roundoff differences can perturb finite-difference Satterthwaite derivatives.
 # Check that implementation independently, then refit the recorded exact H values.
 stopifnot(max(abs(z$Shannon-vegan::diversity(t(rel))))<1e-12)
 z$Participant<-factor(z$Participant);z$SpecimenKey<-factor(z$SpecimenKey)
 z$Location<-factor(z$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'));z$Fraction<-factor(z$Fraction,levels=c('PRE','IgA','IgM','IgG'));z$Phase<-factor(z$Phase,levels=c('Proliferative','Secretory'))
 for(spec in c('all149_nested','all149_nested_phase_adjusted','proliferative_nested','omit99_unsuffixed_endo_nested','omit99_suffix2_endo_nested')){
  d<-z
  if(spec=='proliferative_nested')d<-droplevels(z[z$Phase=='Proliferative',])
  if(spec=='omit99_unsuffixed_endo_nested')d<-droplevels(z[!z$SourceCode%in%c('99EA','99EG','99EM','99EP'),])
  if(spec=='omit99_suffix2_endo_nested')d<-droplevels(z[!z$SourceCode%in%c('99EA2','99EG2','99EM2','99EP2'),])
  form<-Shannon~Location*Fraction+(1|Participant)+(1|SpecimenKey)
  if(spec=='all149_nested_phase_adjusted')form<-update(form,.~.+Phase)
  fit<-lmerTest::lmer(form,data=d,REML=TRUE);a<-as.data.frame(anova(fit,ddf='Satterthwaite'))
  name<-paste(arm,spec,sep='_');a$term<-rownames(a);a$model<-name
  ref<-reference[reference$model==name,];ref<-ref[match(a$term,ref$term),];stopifnot(!anyNA(ref$term))
  for(column in c('NumDF','DenDF','F value','Pr(>F)')){
   difference<-abs(a[[column]]-ref[[column]])
   tests[[length(tests)+1]]<-data.frame(model=name,term=a$term,quantity=column,refit=a[[column]],saved=ref[[column]],absolute_difference=difference)
   if(any(difference>=1e-6)){put(do.call(rbind,tests),'Shannon_partial_comparison_on_failure');print(tail(tests,1));stop('Fresh refit differs from saved output')}
  }
  diagnostics[[length(diagnostics)+1]]<-data.frame(model=name,libraries=nrow(d),women=nlevels(d$Participant),preparations=nlevels(d$SpecimenKey),singular=lme4::isSingular(fit),convergence=paste(fit@optinfo$conv$lme4$messages,collapse=';'))
  if(arm=='baseline_39'&&spec%in%c('all149_nested','all149_nested_phase_adjusted','proliferative_nested')){
   png(file.path(O,paste0(name,'_diagnostics.png')),width=1200,height=1000,res=130)
   par(mfrow=c(2,2),mar=c(4.5,4.5,2.5,1))
   plot(fitted(fit),residuals(fit),xlab='Fitted Shannon',ylab='Residual',pch=16,col='#386985');abline(h=0,lty=2)
   qqnorm(residuals(fit),main='Residual normal quantiles',pch=16,col='#386985');qqline(residuals(fit))
   hist(residuals(fit),main='Residual distribution',xlab='Residual',col='#c4d7e2',border='white')
   plot(d$Shannon,fitted(fit),xlab='Observed Shannon',ylab='Fitted Shannon',pch=16,col='#386985');abline(0,1,lty=2)
   dev.off()
  }
 }
}
put(do.call(rbind,tests),'Shannon_fresh_refit_vs_saved');put(do.call(rbind,diagnostics),'Shannon_fresh_refit_diagnostics')

X<-as.matrix(read.csv(file.path(V,'source_inputs/baseline_39_source_counts.csv'),row.names=1,check.names=FALSE));rel<-sweep(X,2,colSums(X),'/')
reference<-read.csv(file.path(V,'statistics/S2R3_focused_beta.csv'),check.names=FALSE)
dr<-read.csv(file.path(V,'statistics/S2R4_restricted_dispersion.csv'),check.names=FALSE)
out<-list();dout<-list();perm_checks<-list()
options(contrasts=c('contr.sum','contr.poly'))
for(i in seq_len(nrow(reference))){
 r<-reference[i,];obj<-readRDS(file.path(old,'results',paste0(r$analysis,'.rds')));zz<-obj$manifest;P<-obj$permutations
 stopifnot(nrow(P)==9999,all(apply(P,1,function(ii)all(zz$Participant[ii]==zz$Participant))))
 ispre<-startsWith(r$context,'PRE_site');y<-t(rel[,zz$SourceCode,drop=FALSE])
 if(ispre)fit<-vegan::dbrda(y~Location+Condition(Participant),data=zz,distance='bray') else if(r$context=='Fraction_complete')fit<-vegan::dbrda(y~Fraction+Condition(SpecimenKey),data=zz,distance='bray') else{
  mm<-model.matrix(~Location*Fraction,zz);Xi<-mm[,attr(mm,'assign')==3,drop=FALSE]
  fit<-vegan::dbrda(y~Xi+Condition(SpecimenKey+Fraction),data=zz,distance='bray')
 }
 a<-as.data.frame(anova(fit,permutations=P,model='reduced'));R2<-fit$CCA$tot.chi/fit$tot.chi
 vals<-c(F=a$F[1],p=a$`Pr(>F)`[1],R2_total=R2);refs<-c(F=r$F,p=r$`Pr(>F)`,R2_total=r$R2_total)
 stopifnot(max(abs(vals-refs))<1e-10)
 out[[i]]<-data.frame(analysis=r$analysis,quantity=names(vals),refit=vals,saved=refs,absolute_difference=abs(vals-refs))
 if(!ispre){
  stopifnot(all(apply(P,1,function(ii)all(zz$SpecimenKey[ii]==zz$SpecimenKey))))
  for(person in levels(zz$Participant))for(frac in levels(zz$Fraction)){
   ix<-which(zz$Participant==person&zz$Fraction==frac)
   stopifnot(all(apply(P[,ix,drop=FALSE],1,function(ii)length(unique(zz$Fraction[ii]))==1)))
  }
 }
 perm_checks[[i]]<-data.frame(analysis=r$analysis,libraries=nrow(zz),women=nlevels(zz$Participant),permutation_rows=nrow(P),participant_restriction_verified=TRUE,specimen_and_synchronized_fraction_restrictions_verified=if(ispre)NA else TRUE)
 if(r$context!='Interaction_complete'){
  group<-if(ispre)zz$Location else zz$Fraction
  bd<-vegan::betadisper(vegan::vegdist(y,'bray'),group,type='median',bias.adjust=TRUE)
  da<-vegan::permutest(bd,permutations=P)$tab;ref<-dr[dr$analysis==r$analysis,]
  stopifnot(nrow(ref)==1,abs(da$`Pr(>F)`[1]-ref$`Pr(>F)`)<1e-12,abs(da$F[1]-ref$F)<1e-10)
  dout[[length(dout)+1]]<-data.frame(analysis=r$analysis,F_refit=da$F[1],F_saved=ref$F,p_refit=da$`Pr(>F)`[1],p_saved=ref$`Pr(>F)`)
 }
}
put(do.call(rbind,out),'Focused_beta_fresh_refit_vs_saved');put(do.call(rbind,dout),'Dispersion_fresh_refit_vs_saved');put(do.call(rbind,perm_checks),'Permutation_restriction_checks')
capture.output(sessionInfo(),file=file.path(O,'R_sessionInfo.txt'))
writeLines(c('PASS: ten Shannon model refits; four focused beta model refits; three dispersion refits.',
 'All refits use the frozen source-coded counts and recorded model specifications in a fresh R --vanilla session.',
 'The beta and dispersion checks replay saved 9999-row permutation matrices; they do not constitute an independent randomization-design validation or new biological evidence.',
 'This does not reconstruct the raw FASTQ pipeline, classify taxa anew or validate DAPI.'),file.path(O,'REFIT_RESULT.txt'))
cat('PASS: 10 Shannon, 4 focused beta and 3 dispersion models reproduced in fresh R session.\n')
