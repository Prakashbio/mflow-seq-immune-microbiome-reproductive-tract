# Publication-layout candidates from immutable, previously audited results.
# No model is refitted, no gates are moved, and no microscopy image is edited.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);O<-normalizePath(if(length(args))args[1] else 'review_response/Figure_design_review_v22_2026-11-10')
D<-file.path(O,'source_data');F<-file.path(O,'figures');A<-file.path(O,'audit')
captured<-list();envs<-list()
run_recorded<-function(name,argv){
 e<-new.env(parent=globalenv());e$commandArgs<-function(...)argv
 e$ggsave<-function(filename,plot=last_plot(),...){captured[[tools::file_path_sans_ext(basename(filename))]]<<-plot;invisible(NULL)}
 e$png<-e$pdf<-e$dev.off<-e$draw<-function(...)invisible(NULL)
 e$write.csv<-e$writeLines<-e$capture.output<-function(...)invisible(NULL)
 expressions<-parse(file.path(O,'scripts',name))
 for(ex in expressions){
  if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&identical(ex[[2]],as.name('draw')))next
  eval(ex,envir=e)
 }
 envs[[name]]<<-e
 invisible(e)
}
run_recorded('verify_and_plot_PRE_v7_2026_10_29.R',file.path(D,'pre_context'))
run_recorded('plot_diversity_integration_2026_10_28.R',c(D,file.path(D,'statistics')))
run_recorded('plot_taxon_integration_v3_2026_10_28.R',c(D,file.path(D,'taxon_analysis'),file.path(D,'statistics')))
run_recorded('verify_plot_legacy_supplement_v6_2026_10_29.R',D)
run_recorded('plot_verify_source_filter_v5_2026_10_28.R',D)
run_recorded('plot_verify_storage_review_v4_2026_10_28.R',D)
run_recorded('update_storage_figure_notes_v14_2026_11_03.R',D)

# Typography is specified at a 180-mm double-column output width.
style<-theme_classic(base_size=9.5,base_family='Helvetica')+theme(
 text=element_text(colour='#17212B'),axis.text=element_text(size=8,colour='#28343E'),
 axis.title=element_text(size=9),axis.line=element_line(linewidth=.3,colour='#65717A'),
 axis.ticks=element_line(linewidth=.25,colour='#65717A'),
 panel.grid.major.y=element_line(colour='#E9EDF0',linewidth=.25),
 strip.background=element_rect(fill='#F0F3F5',colour=NA),strip.text=element_text(face='bold',size=9),
 panel.spacing=unit(3,'mm'),legend.position='bottom',legend.title=element_text(size=8.5),
 legend.text=element_text(size=8),legend.key.size=unit(3.5,'mm'),
 plot.title=element_text(face='bold',size=10,hjust=0),plot.subtitle=element_text(size=8),
 plot.caption=element_blank(),plot.tag=element_text(face='bold',size=12),
 plot.margin=margin(5,6,5,5))
polish<-function(p){
 if(inherits(p,'patchwork'))return(p&style)
 for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,'GeomPath')||inherits(p$layers[[i]]$geom,'GeomLine')){
  p$layers[[i]]$aes_params$alpha<-.65;p$layers[[i]]$aes_params$linewidth<-.45
 }
 p+style
}
checks<-list();serial<-0
signature<-function(p){
 if(inherits(p,'patchwork'))return(list(lapply(p$patches$plots,signature),signature(structure(p,class=setdiff(class(p),'patchwork')))))
 z<-ggplot_build(p)$data
 lapply(z,function(d)d[,intersect(names(d),c('x','y','xmin','xmax','ymin','ymax','xend','yend','group','PANEL','label')),drop=FALSE])
}
saveplot<-function(p,name,h,width=7.0866,baseline=NULL){
 if(!is.null(baseline)){
  same<-isTRUE(all.equal(signature(baseline),signature(p),check.attributes=TRUE,tolerance=0))
  stopifnot(same);checks[[name]]<<-data.frame(figure=name,check='All plotted coordinates, groups and labels before/after styling',pass=same)
 }
 ggsave(file.path(F,paste0(name,'.png')),p,width=width,height=h,dpi=600,bg='white',limitsize=FALSE)
 ggsave(file.path(F,paste0(name,'.pdf')),p,width=width,height=h,device=grDevices::pdf,useDingbats=FALSE,bg='white',limitsize=FALSE)
 ggsave(file.path(F,paste0(name,'.svg')),p,width=width,height=h,device=svglite::svglite,bg='white',limitsize=FALSE)
 cat('Saved',name,'\n')
}
take<-function(key,name,h){
 p<-captured[[key]];stopifnot(!is.null(p));q<-polish(p)
 if(name=='Supplementary_Figure_S4'){
  e<-envs[['plot_verify_source_filter_v5_2026_10_28.R']]
  bottom<-polish(e$pb)+facet_wrap(~stage,nrow=1,labeller=as_labeller(setNames(c('268 to 41 assignments','41 to 39 assignments'),levels(e$long$stage))))
  q<-polish(e$pa)/bottom+plot_layout(heights=c(1,1.05))
 }
 if(name=='Supplementary_Figure_S9'){
  e<-envs[['plot_verify_storage_review_v4_2026_10_28.R']]
  top<-polish(e$p7a)
  top$layers[[3]]$aes_params$angle<-90;top$layers[[3]]$aes_params$size<-2.8
  top<-top+coord_cartesian(ylim=c(0,120))
  bottom<-polish(e$p7b)+theme(axis.text.x=element_text(angle=55,hjust=1,size=8))
  fill_scale<-bottom$scales$get_scales('fill');fill_scale$breaks<-c(0,50,100)
  q<-top/bottom+plot_layout(heights=c(1,1.15))+plot_annotation(tag_levels='A')
 }
 if(name=='Supplementary_Figure_S10'){
  q$layers[[3]]$mapping$hjust<-aes(hjust=ifelse(Storage=='Fresh',0,1))$hjust
  pct1<-sub('.*\\(([^%]+)%.*','\\1',p$labels$x);pct2<-sub('.*\\(([^%]+)%.*','\\1',p$labels$y)
  q<-q+labs(x=paste0('PRE PCoA 1 (',pct1,'%)'),y=paste0('PRE PCoA 2 (',pct2,'%)'))
 }
 if(name=='Supplementary_Figure_S11')q<-q+theme(legend.position='none')
 if(name=='Supplementary_Figure_S12')q<-q+scale_x_discrete(labels=c('V','C','E'))
 if(name=='Supplementary_Figure_S13')q<-q+scale_colour_manual(values=c(Proliferative='#D55E00',Secretory='#0072B2'))
 if(name=='Supplementary_Figure_S14')q<-q+theme(panel.spacing.x=unit(6,'mm'))
 if(name=='Supplementary_Figure_S15')q<-q+facet_grid(Taxon~Location,labeller=labeller(Taxon=c('Lactobacillus crispatus'='L. crispatus','Lactobacillus iners'='L. iners','Lactobacillus jensenii'='L. jensenii','Lactobacillus gasseri'='L. gasseri')))+theme(legend.position='none')
 if(name=='Supplementary_Figure_S16')q<-q+scale_shape_manual(values=c(16,2),labels=c('3+ women detected','Fewer than 3 women detected'))+theme(legend.box='vertical')
 if(name=='Supplementary_Figure_S17')q<-q+theme(legend.position='none')
 saveplot(q,name,h,baseline=p);q
}

# Main 2: preserve the composition and Shannon displays, rebalance their areas.
e<-envs[['verify_and_plot_PRE_v7_2026_10_29.R']]
p2a<-polish(e$ga)+labs(title='A  PRE genus composition',x='Participant')+guides(fill=guide_legend(nrow=2,byrow=TRUE))+theme(legend.text=element_text(size=8),axis.text.x=element_text(size=8))
p2b<-polish(e$gb)+labs(title='B  PRE diversity within women',y='Shannon entropy')+theme(legend.position='none')
stopifnot(isTRUE(all.equal(signature(e$ga),signature(p2a),tolerance=0)),isTRUE(all.equal(signature(e$gb),signature(p2b),tolerance=0)))
saveplot(p2a/p2b+plot_layout(heights=c(2.1,1)),'Main_Figure_2',7.4)

# Main 3: A-C use the same reconstructed gate values; D uses the already saved
# paired IgM-minus-IgA contrasts. No new bootstrap or significance test is run.
r<-read.csv(file.path(D,'FACS_reconstruction/data/Rebuilt_primary_42_woman_site_records.csv'))
s<-read.csv(file.path(D,'FACS_reconstruction/data/Rebuilt_Table_2_and_pointwise_CIs.csv'))
con<-read.csv(file.path(D,'FACS_reconstruction/data/FACS_paired_IgM_minus_IgA_v12.csv'))
con<-con[con$subset=='All',]
sites<-c('VAGINA','CERVIX','ENDOMETRIUM');sn<-c('Vagina','Cervix','Endometrium')
endpoints<-c('IgA_only','IgM_only','Double_positive','IgG_positive');gatecols<-c('#D97727','#7B59A6','#697781','#16826F')
long<-do.call(rbind,lapply(endpoints,function(v)data.frame(Participant=r$Participant,Location=r$Location,endpoint=v,percent=r[[paste0(v,'_percent')]])))
long$endpoint<-factor(long$endpoint,levels=endpoints);s$endpoint<-factor(s$endpoint,levels=endpoints)
stopifnot(nrow(long)==168,nrow(s)==12,nrow(con)==3)
ps<-list()
for(i in 1:3){
 z<-long[long$Location==sites[i],];m<-s[s$Location==sites[i],]
 ps[[i]]<-ggplot(z,aes(endpoint,percent))+geom_point(aes(colour=endpoint),position=position_jitter(width=.13,height=0,seed=20261110),size=1.5,alpha=.78)+
 geom_errorbar(data=m,aes(y=mean_percent,ymin=ci95_low_percent,ymax=ci95_high_percent),width=.18,linewidth=.6)+
 geom_point(data=m,aes(y=mean_percent),shape=18,size=2.8)+scale_colour_manual(values=gatecols,guide='none')+
 scale_x_discrete(labels=c('IgA only','IgM only','Double +','IgG +'))+scale_y_continuous(limits=c(0,50),breaks=seq(0,50,10),expand=expansion(mult=c(.02,.04)))+
 style+labs(title=paste0(LETTERS[i],'  ',sn[i],' (',length(unique(z$Participant)),' women)'),x=NULL,y='Ig-gate events (% of DAPI-negative parent)')+
 theme(axis.text.x=element_text(size=8),axis.title.y=element_text(size=8.5))
}
con$Location<-factor(con$Location,levels=rev(sites),labels=rev(sn))
ps[[4]]<-ggplot(con,aes(mean_M_minus_A_pp,Location))+geom_vline(xintercept=0,linetype=2,colour='#A0A7AD',linewidth=.35)+geom_errorbar(aes(xmin=ci95_low_pp,xmax=ci95_high_pp),orientation='y',width=.12,linewidth=.7,colour='#435B70')+geom_point(size=2.5,shape=18,colour='#435B70')+style+labs(title='D  Paired IgM minus IgA',x='Mean difference (percentage points)',y=NULL)+theme(panel.grid.major.y=element_blank(),axis.title.x=element_text(size=8))
saveplot(wrap_plots(ps,ncol=2),'Main_Figure_3',6.2)
write.csv(long,file.path(A,'Main_Figure_3_points.csv'),row.names=FALSE)
write.csv(s,file.path(A,'Main_Figure_3_intervals.csv'),row.names=FALSE)
write.csv(con,file.path(A,'Main_Figure_3_paired_contrasts.csv'),row.names=FALSE)

# Main 4 retains the preparation-linked Shannon display; visual density alone
# is not a reason to add mathematically related indices to the main result.
p4<-captured[['Figure2_Shannon_preparation_pairs']]
p4new<-polish(p4)+scale_colour_manual(values=c(Proliferative='#D55E00',Secretory='#0072B2'))+theme(legend.title=element_text(size=8),legend.text=element_text(size=8))
saveplot(p4new,'Main_Figure_4',3.6,baseline=p4)

# Main 5 candidate combines descriptive ordination with the already calculated
# sorted-to-own-PRE distances expressly requested by the reviewer.
p5<-polish(captured[['Figure3_PCoA_participant_pairs']])+labs(title='A  Community ordination')+guides(colour=guide_legend(nrow=3,override.aes=list(shape=16,alpha=1)),shape=guide_legend(nrow=1))+theme(legend.text=element_text(size=8),legend.key.width=unit(2.5,'mm'),legend.box='vertical')
p5b<-polish(captured[['FigureS4_matched_own_PRE']])+scale_colour_manual(values=c(Proliferative='#D55E00',Secretory='#0072B2'))+labs(title='B  Matched distance to own PRE',y='Bray-Curtis distance to own PRE')+theme(legend.title=element_text(size=8))
saveplot(p5/p5b+plot_layout(heights=c(1.4,1)),'Main_Figure_5',7.2)

# Main 6 candidate adds the corresponding stored mean abundances below the
# interval display, distinguishing relative representation from abundance.
p6base<-captured[['Figure5_matched_site_scores']]
p6<-polish(p6base)+labs(title='A  Matched relative representation')+scale_shape_manual(values=c(16,2),labels=c('3+ women detected','Fewer than 3 women detected'))+theme(axis.text.y=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.box='vertical')
scores<-read.csv(file.path(D,'taxon_statistics/Site_specific_matched_IgA_IgM_representation_scores.csv'))
keys<-c('Fannyhessea vaginae','Gardnerella vaginalis','Lactobacillus crispatus','Lactobacillus gasseri','Lactobacillus iners','Lactobacillus jensenii','Prevotella buccalis','Sneathia vaginalis')
z<-scores[scores$taxon_arm=='baseline_39'&scores$rank=='provisional_species'&scores$scope=='All'&scores$Taxon%in%keys,]
stopifnot(nrow(z)==24)
h<-rbind(data.frame(z[,c('Taxon','Location')],Fraction='IgA',pct=100*z$mean_IgA_relative_abundance),data.frame(z[,c('Taxon','Location')],Fraction='IgM',pct=100*z$mean_IgM_relative_abundance))
h$Taxon<-factor(h$Taxon,levels=rev(keys));h$Location<-factor(h$Location,levels=sites,labels=sn);h$Fraction<-factor(h$Fraction,levels=c('IgA','IgM'))
ph<-ggplot(h,aes(Fraction,Taxon,fill=pct))+geom_tile(colour='white',linewidth=.5)+geom_text(aes(label=sprintf('%.1f',pct)),size=2.9,colour='black')+facet_wrap(~Location,nrow=1)+scale_fill_gradient(low='#F6F1DD',high='#85B6D2',limits=c(0,100),breaks=c(0,50,100),name='Mean abundance (%)')+style+theme(axis.text.y=element_text(size=8),panel.grid=element_blank())+labs(title='B  Corresponding mean relative abundance',x=NULL,y=NULL)
saveplot(p6/ph+plot_layout(heights=c(1.45,1)),'Main_Figure_6',7.8)
write.csv(h,file.path(A,'Main_Figure_6_abundance_panel.csv'),row.names=FALSE)

# Existing supplementary analyses: every plotted coordinate/group/label is
# checked before/after the typography and spacing changes.
mapping<-list(
 c('FigureS14_filtering_reconstruction','Supplementary_Figure_S4','6.1'),
 c('FigureS6_PRE_descriptive_rank_correlations','Supplementary_Figure_S5','6.5'),
 c('FigureS8_PRE_participant_accumulation','Supplementary_Figure_S6','3.6'),
 c('FigureS11_PRE_and_retained_count_descriptions','Supplementary_Figure_S7','6.8'),
 c('FigureS2_storage_alpha_descriptive','Supplementary_Figure_S8','7.4'),
 c('Figure7_descriptive_storage','Supplementary_Figure_S9','7.1'),
 c('Supplementary_Figure_S10','Supplementary_Figure_S10','3.8'),
 c('Supplementary_Figure_S11','Supplementary_Figure_S11','5.4'),
 c('FigureS3_Shannon_woman_site_pairs','Supplementary_Figure_S12','4.1'),
 c('FigureS4_matched_own_PRE','Supplementary_Figure_S13','3.7'),
 c('FigureS9_matched_PRE_contrast_intervals','Supplementary_Figure_S14','3.4'),
 c('FigureS12_Lactobacillus_paired_values','Supplementary_Figure_S15','6.4'),
 c('FigureS10_genus_representation','Supplementary_Figure_S16','8.2'),
 c('Figure6_PRE_group_descriptives','Supplementary_Figure_S17','8.4'))
for(m in mapping){
 if(m[2]=='Supplementary_Figure_S5'){
  e<-envs[['verify_plot_legacy_supplement_v6_2026_10_29.R']]
  hp<-lapply(e$plots,function(p)polish(p)+theme(panel.grid=element_blank())+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),breaks=c(-1,0,1),na.value='#BBBBBB',name='Rank correlation')+guides(fill=guide_colourbar(barwidth=unit(28,'mm'))))
  key<-e$legend+theme_void()+theme(plot.margin=margin(5,0,0,0))
  p<-(hp[[1]]+hp[[2]]+hp[[3]]+key)+plot_layout(ncol=2,guides='collect')&theme(legend.position='bottom')
  saveplot(p,m[2],as.numeric(m[3]),baseline=captured[[m[1]]])
 }else take(m[1],m[2],as.numeric(m[3]))
}

# Optional S18: a new presentation of already saved four-index measurements.
# No new inferential tests. Both preparations of woman 99 remain displayed.
metrics<-read.csv(file.path(D,'storage_statistics/Storage_library_metrics.csv'))
metrics<-metrics[metrics$taxon_arm=='baseline_39',]
phase<-read.csv(file.path(D,'statistics/baseline_39_Shannon_per_source_library.csv'))
metrics$Phase<-phase$Phase[match(metrics$SourceCode,phase$SourceCode)]
stopifnot(nrow(metrics)==149,!anyNA(metrics$Phase))
v<-c('Shannon','Simpson_1_minus_D','Richness','Pielou_evenness')
z<-do.call(rbind,lapply(v,function(k)data.frame(metrics[,c('SourceCode','Participant','SpecimenKey','Location','Fraction','Phase')],metric=k,value=metrics[[k]])))
z$metric<-factor(z$metric,levels=v,labels=c('Shannon entropy','Simpson diversity (1-D)','Observed richness','Pielou evenness'))
z$Location<-factor(z$Location,levels=sites,labels=sn);z$Fraction<-factor(z$Fraction,levels=c('PRE','IgA','IgM','IgG'))
stopifnot(nrow(z)==596,sum(is.na(z$value))==6)
p18<-ggplot(z,aes(Fraction,value,group=SpecimenKey))+geom_line(colour='#A7B0B7',alpha=.6,linewidth=.35,na.rm=TRUE)+geom_point(aes(colour=Phase),size=1.1,alpha=.85,na.rm=TRUE)+facet_grid(metric~Location,scales='free_y',switch='y')+scale_colour_manual(values=c(Proliferative='#D55E00',Secretory='#0072B2'))+style+theme(strip.placement='outside',strip.text.y.left=element_text(angle=90,size=8),axis.text.x=element_text(size=8))+labs(x='Library fraction',y=NULL,colour='Recorded cycle phase')
saveplot(p18,'OPTIONAL_S18_Alpha_metrics_by_site_and_fraction',8.3)
write.csv(z,file.path(A,'Optional_S18_plotted_values.csv'),row.names=FALSE)
# Verify the preserved components in the new multi-panel combinations.
composites<-list(Main2_A=list(envs[['verify_and_plot_PRE_v7_2026_10_29.R']]$ga,p2a),Main2_B=list(envs[['verify_and_plot_PRE_v7_2026_10_29.R']]$gb,p2b),Main5_A=list(captured[['Figure3_PCoA_participant_pairs']],p5),Main5_B=list(captured[['FigureS4_matched_own_PRE']],p5b),Main6_A=list(p6base,p6))
for(key in names(composites)){
 same<-isTRUE(all.equal(signature(composites[[key]][[1]]),signature(composites[[key]][[2]]),tolerance=0))
 stopifnot(same);checks[[key]]<-data.frame(figure=key,check='All original plotted coordinates, groups and labels in composite component',pass=same)
}
saveRDS(composites,file.path(A,'Composite_plot_components.rds'))
write.csv(do.call(rbind,checks),file.path(A,'Theme_only_coordinate_preservation.csv'),row.names=FALSE)
saveRDS(captured,file.path(A,'Recorded_plot_objects.rds'))
capture.output(sessionInfo(),file=file.path(A,'R_sessionInfo.txt'))
capture.output(warnings(),file=file.path(A,'Warnings.txt'))
cat('PASS: figure layout candidates completed. No sealed document or package modified.\n')
