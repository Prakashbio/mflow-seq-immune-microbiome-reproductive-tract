suppressPackageStartupMessages(library(lmerTest))
o<-'review_response/03_analysis/Prakash_original_workflow_review_2026-10-27'
z<-read.csv(file.path(o,'results/baseline_39_Shannon_per_source_library.csv'))
z<-z[z$Phase=='Proliferative',]
subsets<-list(omit99_unsuffixed_endo=z[!z$SourceCode%in%c('99EA','99EG','99EM','99EP'),],omit99_suffix2_endo=z[!z$SourceCode%in%c('99EA2','99EG2','99EM2','99EP2'),])
for(p in unique(z$Participant))subsets[[paste0('leave_participant_',p,'_out')]]<-z[z$Participant!=p,]
out<-list()
for(n in names(subsets)) {
 d<-subsets[[n]];d$Participant<-factor(d$Participant);d$SpecimenKey<-factor(d$SpecimenKey)
 d$Location<-factor(d$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'));d$Fraction<-factor(d$Fraction,levels=c('PRE','IgA','IgM','IgG'))
 f<-lmer(Shannon~Location*Fraction+(1|Participant)+(1|SpecimenKey),d,REML=TRUE)
 a<-anova(f);out[[n]]<-data.frame(sensitivity=n,term=rownames(a),p_nominal=a$`Pr(>F)`,p_BH_three_omnibus=p.adjust(a$`Pr(>F)`,'BH'),libraries=nrow(d),women=nlevels(d$Participant),singular=isSingular(f),convergence=paste(f@optinfo$conv$lme4$messages,collapse=';'))
 saveRDS(f,file.path(o,'results',paste0('proliferative39_',n,'.rds')))
}
ans<-do.call(rbind,out);write.csv(ans,file.path(o,'results/Proliferative39_influence_checks.csv'),row.names=FALSE)
print(ans[ans$term=='Location:Fraction',])
