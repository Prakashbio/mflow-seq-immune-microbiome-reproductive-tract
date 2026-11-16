# Reviewer-scoped plots from named-source audited outputs. No label editing.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(TRUE);out<-if(length(args))args[1] else 'review_response/taxon_integrated_working_v3_2026-10-28'
src<-if(length(args)>1)args[2] else 'review_response/03_analysis/Taxon_scores_PRE_review_2026-10-28'
beta<-if(length(args)>2)args[3] else 'review_response/diversity_integrated_working_v2_2026-10-28/statistics'
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,'taxon_statistics'),showWarnings=FALSE)
put<-function(z,n)write.csv(z,file.path(out,'taxon_statistics',paste0(n,'.csv')),row.names=FALSE)
s<-read.csv(file.path(src,'results/Site_specific_matched_IgA_IgM_representation_scores.csv'))
w<-read.csv(file.path(src,'results/Woman_site_IgA_IgM_relative_abundances.csv'))
bs<-readRDS(file.path(src,'results/Score_bootstrap_draws_and_weights.rds'))
sites<-c('VAGINA','CERVIX','ENDOMETRIUM');site_names<-c('Vagina','Cervix','Endometrium')
theme<-theme_bw(base_size=14)+theme(legend.position='bottom',panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey94'),plot.caption=element_text(size=11,hjust=0))
savefig<-function(p,n,width,height){ggsave(file.path(out,'figures',paste0(n,'.png')),p,width=width,height=height,dpi=300,bg='white');ggsave(file.path(out,'figures',paste0(n,'.pdf')),p,width=width,height=height,bg='white')}
format_site<-function(z){z$Location<-factor(z$Location,levels=sites,labels=site_names);z}
scoreplot<-function(z,caption){
 z<-format_site(z);z$Taxon<-factor(z$Taxon,levels=rev(sort(unique(z$Taxon))));z$sparse<-factor(ifelse(z$women_detected_either<3,'Fewer than 3 women detected','At least 3 women detected'))
 ggplot(z,aes(Score_ratio_of_matched_woman_means,Taxon))+geom_vline(xintercept=.5,lty=2,colour='grey65')+geom_errorbar(aes(xmin=CI95_low,xmax=CI95_high),orientation='y',width=.15,linewidth=.45)+geom_point(aes(size=women_detected_either,shape=sparse),colour='#235f84')+scale_shape_manual(values=c(16,2))+facet_wrap(~Location,nrow=1)+scale_x_continuous(breaks=c(0,.5,1))+coord_cartesian(xlim=c(0,1))+theme+guides(shape=guide_legend(nrow=1),size=guide_legend(nrow=1))+labs(x='IgA / (IgA + IgM) relative sequencing representation',y=NULL,size='Women detected',shape=NULL,caption=caption)
}
key<-c('Lactobacillus crispatus','Lactobacillus iners','Lactobacillus jensenii','Lactobacillus gasseri','Sneathia vaginalis','Gardnerella vaginalis','Prevotella buccalis','Fannyhessea vaginae')
z<-s[s$taxon_arm=='baseline_39'&s$rank=='provisional_species'&s$scope=='All'&s$Taxon%in%key,]
savefig(scoreplot(z,'Eight previously discussed provisional assignments; matched preparations and equal weight to women.\nPointwise 95% bootstrap intervals; zero-width endpoint intervals do not establish exclusive binding.'),'Figure5_matched_site_scores',12,6.8)
z<-s[s$taxon_arm=='baseline_39'&s$rank=='genus_aggregation'&s$scope=='All',]
savefig(scoreplot(z,'Genus aggregation of the 39-taxon baseline. Matched IgA/IgM sets: V15, C9 and E12 women.\nPointwise conditional intervals; sparse detections and undefined draws are tabulated. Aggregation does not validate classification.'),'FigureS10_genus_representation',12.5,10.3)
# Additional descriptive abundance differences use the SAME saved woman draws.
diffs<-list()
for(site in sites){
 name<-paste('baseline_39_provisional_species_All',site,sep='_');obj<-bs[[name]];A<-M<-matrix(0,length(obj$women),length(obj$taxa))
 z<-w[w$taxon_arm=='baseline_39'&w$rank=='provisional_species'&w$scope=='All'&w$Location==site,]
 for(i in seq_along(obj$women)){zz<-z[z$Participant==obj$women[i],];zz<-zz[match(obj$taxa,zz$Taxon),];A[i,]<-zz$IgA_relative_mean;M[i,]<-zz$IgM_relative_mean}
 delta<-100*(A-M);boot<-obj$draw_weights%*%delta/length(obj$women);ci<-apply(boot,2,quantile,c(.025,.975),names=FALSE)
 diffs[[site]]<-data.frame(Location=site,Taxon=obj$taxa,women=length(obj$women),women_detected=colSums(A+M>0),mean_difference_pp=colMeans(delta),CI95_low=ci[1,],CI95_high=ci[2,],bootstrap_draws=nrow(boot))
 weights<-as.data.frame(obj$draw_weights);names(weights)<-paste0('woman_',obj$women);put(weights,paste0('Absolute_difference_bootstrap_weights_',site))
}
dd<-do.call(rbind,diffs);put(dd,'Matched_mean_abundance_differences_pp');dd<-format_site(dd);dd$Taxon<-factor(dd$Taxon,levels=rev(sort(unique(dd$Taxon))))
p<-ggplot(dd,aes(mean_difference_pp,Taxon))+geom_vline(xintercept=0,lty=2,colour='grey60')+geom_errorbar(aes(xmin=CI95_low,xmax=CI95_high),orientation='y',width=.18,linewidth=.4)+geom_point(aes(size=women_detected),colour='#235f84')+facet_wrap(~Location,nrow=1)+theme+theme(axis.text.y=element_text(size=13))+labs(x='Mean IgA minus IgM relative abundance (percentage points)',y=NULL,size='Women detected',caption='All 39 provisional assignments. Preparation pairs are averaged within woman before estimating site means.\nBars are pointwise 95% intervals from 5,000 saved woman-bootstrap draws; no multiple-comparison significance claim.')
savefig(p,'Figure8_matched_abundance_differences',13,12.5)
# PRE groups: select only one source preparation per woman/site, without
# treating repeated endometrial preparations as extra women.
v<-read.csv(file.path(src,'results/PRE_taxon_community_per_woman_values.csv'))
tests<-read.csv(file.path(src,'results/PRE_taxon_community_exact_tests.csv'))
keys<-c('Lactobacillus crispatus','Lactobacillus iners','Sneathia vaginalis','Gardnerella leopoldii','Prevotella oris','Parvimonas micra')
z<-v[v$taxon_arm=='baseline_39'&v$rank=='provisional_species'&v$selection=='retain99EP2'&v$Taxon%in%keys,]
stopifnot(nrow(z)==length(keys)*36);put(z,'Figure6_plotted_PRE_values')
z<-format_site(z);z$Taxon<-factor(z$Taxon,levels=keys,labels=c('L. crispatus','L. iners','S. vaginalis','G. leopoldii','P. oris','P. micra'));z$group<-factor(ifelse(z$Category=='LACTOBACILLUS_DESCRIPTION','L-desc','D-desc'),levels=c('L-desc','D-desc'))
med<-aggregate(100*z$relative_abundance,by=list(Taxon=z$Taxon,Location=z$Location,group=z$group),FUN=median);names(med)[4]<-'median_pct';put(med,'Figure6_PRE_medians')
p<-ggplot(z,aes(group,relative_abundance*100,colour=group))+geom_point(position=position_jitter(width=.1,height=0,seed=20261028),size=1.4,alpha=.8)+geom_errorbar(data=med,aes(x=group,ymin=median_pct,ymax=median_pct),inherit.aes=FALSE,width=.5,colour='black',linewidth=.3)+facet_grid(Taxon~Location)+scale_y_sqrt(breaks=c(0,1,4,16,64,100))+scale_colour_manual(values=c('#30698a','#a44448'))+theme+theme(legend.position='none',strip.text.y=element_text(size=12))+labs(x='Original vaginal-community description',y='PRE relative abundance (%) on square-root scale',caption='Each point is one woman at one site; black bars mark medians. No provisional species survives within-site BH.\nL-desc/D-desc are original composition descriptions, not clinical diagnoses. Both participant-99 PRE choices are tabulated.')
savefig(p,'Figure6_PRE_group_descriptives',11.5,11.5)
t<-tests[tests$taxon_arm=='baseline_39'&tests$rank=='provisional_species'&tests$selection=='retain99EP2'&is.finite(tests$p_exact_label_permutation),];t<-format_site(t)
p<-ggplot(t,aes(100*median_difference_D_minus_L,-log10(p_exact_label_permutation)))+geom_point(colour='grey40',size=1.8)+facet_wrap(~Location,nrow=1)+theme+labs(x='Median PRE abundance difference, D-desc minus L-desc (percentage points)',y='Minus log10 exact nominal p',caption='All estimable provisional species in the 39-taxon baseline. None has within-site BH q below 0.05.\nNominal p is displayed on the y-axis; no point is labelled as an established enrichment finding.')
savefig(p,'FigureS5_PRE_exact_nominal_p',11.5,4.8)
# Show the individual paired relative library abundances underlying the
# Lactobacillus score discussion rather than a community-state-type inference.
z<-w[w$taxon_arm=='baseline_39'&w$rank=='provisional_species'&w$scope=='All'&w$Taxon%in%key[1:4],]
zz<-rbind(data.frame(z[,c('Location','Taxon','Participant')],Fraction='IgA',RA=100*z$IgA_relative_mean),data.frame(z[,c('Location','Taxon','Participant')],Fraction='IgM',RA=100*z$IgM_relative_mean));zz<-format_site(zz);zz$Taxon<-factor(zz$Taxon,levels=key[1:4]);zz$Fraction<-factor(zz$Fraction,levels=c('IgA','IgM'))
put(zz,'FigureS12_Lactobacillus_plotted_values')
p<-ggplot(zz,aes(Fraction,RA,group=Participant))+geom_line(colour='grey65',alpha=.5,linewidth=.3)+geom_point(aes(colour=Fraction),size=1.3)+facet_grid(Taxon~Location)+theme+scale_colour_manual(values=c('#30698a','#a44448'))+theme(legend.position='none',strip.text.y=element_text(size=12))+labs(x='Matched sorted-library fraction',y='Relative abundance (%)',caption='Four provisional Lactobacillus assignments; lines pair IgA/IgM library abundances within woman/site.\nRepeated preparations are averaged within woman. These values do not measure antibody coating per bacterial cell.')
savefig(p,'FigureS12_Lactobacillus_paired_values',11.5,9)
c<-read.csv(file.path(beta,'Original_source_matched_PRE_contrasts_with_exact_tests.csv'));c<-c[c$taxon_arm=='baseline_39'&c$scope=='All',];c<-format_site(c)
p<-ggplot(c,aes(mean_difference_equal_weight_women,contrast))+geom_vline(xintercept=0,lty=2,colour='grey60')+geom_errorbar(aes(xmin=CI95_low,xmax=CI95_high),orientation='y',width=.15)+geom_point(colour='#235f84',size=2)+facet_wrap(~Location,nrow=1)+theme+labs(x='Difference between sorted-to-own-PRE distances',y=NULL,caption='Contrasts average repeated-preparation differences within woman; pointwise 95% woman-bootstrap intervals.\nNo contrast survives BH across nine exact sign-flip tests. Interval exclusion of zero is not multiplicity-adjusted significance.')
savefig(p,'FigureS9_matched_PRE_contrast_intervals',12,3.8)
capture.output(sessionInfo(),file=file.path(out,'figures/Taxon_plot_R_sessionInfo.txt'))
cat('Seven PNG/PDF figure pairs and matched abundance-difference data saved.\n')
