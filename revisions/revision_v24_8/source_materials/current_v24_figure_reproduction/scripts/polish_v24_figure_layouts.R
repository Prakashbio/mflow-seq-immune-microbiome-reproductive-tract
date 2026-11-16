# Layout changes only. Frozen results are copied to temporary storage before use.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);src<-normalizePath(args[1]);out<-args[2];stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
for(n in c('figures','audit'))dir.create(file.path(out,n))
work<-tempfile('mflow_v24_capture_');dir.create(work)
for(n in c('source_data','scripts'))stopifnot(file.copy(file.path(src,n),work,recursive=TRUE))
for(n in c('figures','audit'))dir.create(file.path(work,n))
cap<-new.env(parent=globalenv());figs<-list();cap$commandArgs<-function(...)work
cap$ggsave<-function(filename,plot=last_plot(),...){figs[[tools::file_path_sans_ext(basename(filename))]]<<-plot;invisible(NULL)}
cap$write.csv<-cap$writeLines<-cap$saveRDS<-cap$capture.output<-function(...)invisible(NULL)
for(ex in parse(file.path(work,'scripts/redesign_v22_figures.R')))eval(ex,envir=cap)
mapping<-c('4'='4','5'='5','6'='6','7'='7','12'='8','13'='UNUSED','14'='10','15'='11','16'='12','17'='13','8'='14','9'='15','10'='16','11'='17')
old<-list()
for(i in c(3,4,6))old[[paste0('Main_Figure_',i)]]<-figs[[paste0('Main_Figure_',i)]]
for(i in names(mapping))if(mapping[[i]]!='UNUSED')old[[paste0('Supplementary_Figure_S',mapping[[i]])]]<-figs[[paste0('Supplementary_Figure_S',i)]]
old[['Supplementary_Figure_S9']]<-figs[['OPTIONAL_S18_Alpha_metrics_by_site_and_fraction']]
stopifnot(length(old)==17,all(vapply(old,function(x)!is.null(x),logical(1))))
# At the retained 6.15-inch Word width, 9.25 native points print at about 8 pt.
big<-theme(axis.text=element_text(size=9.25),axis.title=element_text(size=9.5),
 strip.text=element_text(size=9.5),legend.text=element_text(size=9.25),legend.title=element_text(size=9.25,face='bold'),
 plot.title=element_text(size=10.5,face='bold'),legend.key.size=unit(3.5,'mm'),
 legend.margin=margin(1,1,1,1),legend.box.margin=margin(1,0,1,0),legend.spacing.y=unit(1,'mm'))
pol<-function(p){if(inherits(p,'patchwork'))p&big else p+big}
new<-lapply(old,pol)
# The fourth S5 component is a text-only taxon key, so it must stay axis-free.
new[['Supplementary_Figure_S5']][[4]]<-new[['Supplementary_Figure_S5']][[4]]+theme_void()+theme(plot.margin=margin(5,0,0,0))
abbr<-function(x)sub('^Lactobacillus ','L. ',x)
compact_detection<-function(p)p+guides(size=guide_legend(order=1,nrow=1,title.position='left'),shape=guide_legend(order=2,nrow=1,title.position='left',override.aes=list(size=2.5)))+
 theme(legend.box='vertical',legend.spacing.y=unit(.5,'mm'),legend.key.width=unit(3,'mm'))

# Main 3: wrap long units, keeping every event value, summary and CI.
ps<-lapply(cap$ps,pol)
for(i in 1:3)ps[[i]]<-ps[[i]]+labs(y='Ig-gate events\n(% of DAPI-negative parent)')
new[['Main_Figure_3']]<-wrap_plots(ps,ncol=2)
new[['Main_Figure_4']]<-new[['Main_Figure_4']]+guides(colour=guide_legend(nrow=1,title.position='left'))

# Main 6: reduce repeated long genus labels and compact detection keys.
p6<-compact_detection(pol(cap$p6))+scale_y_discrete(labels=abbr)+theme(axis.text.y=element_text(size=9.25,face='italic'))
ph<-pol(cap$ph)+scale_y_discrete(labels=abbr)+theme(axis.text.y=element_text(size=9.25,face='italic'))
ph$layers[[2]]$aes_params$size<-3.4
ph<-ph+guides(fill=guide_colourbar(title.position='top',title.hjust=.5,barwidth=unit(30,'mm'),barheight=unit(3,'mm')))
new[['Main_Figure_6']]<-p6/ph+plot_layout(heights=c(1.45,1))

# Four long analytical selection rules fit in two clearly ordered rows.
new[['Supplementary_Figure_S4']][[1]]<-new[['Supplementary_Figure_S4']][[1]]+guides(fill=guide_legend(nrow=2,byrow=TRUE))
# Both accumulation keys retain all levels; separate horizontal rows avoid overflow.
new[['Supplementary_Figure_S6']]<-new[['Supplementary_Figure_S6']]+guides(colour=guide_legend(order=1,nrow=1,title.position='left'),linetype=guide_legend(order=2,nrow=1,title.position='left'))+theme(legend.box='vertical')

# S7: consistent site ordering and a single horizontal key in each relevant panel.
leg<-cap$envs[['verify_plot_legacy_supplement_v6_2026_10_29.R']]
pa<-pol(cap$polish(leg$pA))+labs(title='A  PRE trajectories\n6 women')
pb<-pol(cap$polish(leg$pB))+labs(title='B  Retained counts\n149 libraries')+guides(colour=guide_legend(nrow=1,byrow=TRUE))
pc<-pol(cap$polish(leg$pC))+labs(title='C  PRE distances\nare descriptive')+guides(fill=guide_legend(nrow=1))
pd<-pol(cap$polish(leg$pD))+labs(title='D  Count-diversity\ndisplay')+guides(colour=guide_legend(nrow=1,byrow=TRUE))
for(pname in c('pb','pd')){
 pp<-get(pname);sc<-pp$scales$get_scales('colour');sc$breaks<-c('VAGINA','CERVIX','ENDOMETRIUM');assign(pname,pp)
}
new[['Supplementary_Figure_S7']]<-(pa+pb)/(pc+pd)+plot_layout(guides='keep')

# S8: all 15 participant codes remain in one horizontal row, with the original colours.
new[['Supplementary_Figure_S8']]<-new[['Supplementary_Figure_S8']]+guides(colour=guide_legend(nrow=1,byrow=TRUE,title.position='left',override.aes=list(linetype=0,alpha=1)))+theme(legend.key.width=unit(2.1,'mm'),legend.spacing.x=unit(.6,'mm'))
new[['Supplementary_Figure_S9']]<-new[['Supplementary_Figure_S9']]+guides(colour=guide_legend(nrow=1,title.position='left'))+theme(strip.text.y.left=element_text(size=9.25))
new[['Supplementary_Figure_S12']]<-compact_detection(new[['Supplementary_Figure_S12']])+theme(axis.text.y=element_text(size=9.25))

# Storage figures: enlarge metric names and simplify the crowded heatmap axes.
new[['Supplementary_Figure_S14']]<-new[['Supplementary_Figure_S14']]+guides(colour=guide_legend(nrow=1,title.position='left'))+theme(strip.text.y=element_text(size=9.25))
storage<-cap$envs[['plot_verify_storage_review_v4_2026_10_28.R']]
st<-pol(cap$polish(storage$p7a))+coord_cartesian(ylim=c(0,120))+labs(title='A  Lactobacillus abundance')
st$layers[[3]]$aes_params$angle<-0;st$layers[[3]]$aes_params$size<-3.3
count_labels<-st$layers[[3]]$data;old_count_labels<-count_labels$label
stopifnot(all(grepl('^n=[0-9]+$',old_count_labels)))
count_labels$label<-sub('^n=','',old_count_labels);st$layers[[3]]$data<-count_labels
write.csv(data.frame(old_label=old_count_labels,new_label=count_labels$label,identical_count=as.integer(sub('^n=','',old_count_labels))==as.integer(count_labels$label)),file.path(out,'audit/S15_count_label_values.csv'),row.names=FALSE)
short_storage<-function(x)gsub('\n',' ',gsub('Frozen','Fz',gsub('Fresh','Fr',x,fixed=TRUE),fixed=TRUE),fixed=TRUE)
sb<-pol(cap$polish(storage$p7b))+labs(title='B  Mean relative abundance',x='Fraction and storage (Fr, fresh; Fz, frozen)')+
 scale_y_discrete(labels=abbr)+scale_x_discrete(labels=short_storage)+theme(axis.text.x=element_text(angle=45,hjust=1,size=9.25),axis.text.y=element_text(size=9.25,face='italic'))
storage_fill_scale<-sb$scales$get_scales('fill');storage_fill_scale$breaks<-c(0,50,100)
sb<-sb+guides(fill=guide_colourbar(title.position='top',title.hjust=.5,barwidth=unit(30,'mm'),barheight=unit(3,'mm')))
new[['Supplementary_Figure_S15']]<-st/sb+plot_layout(heights=c(1,1.15))

# Preserve the exact data-layer signatures, including all point/interval/box/
# tile coordinates, groups, label values, colours, fill, symbols and line types.
signature<-function(p,normalize_counts=FALSE){
 if(inherits(p,'patchwork'))return(list(lapply(p$patches$plots,signature,normalize_counts=normalize_counts),signature(structure(p,class=setdiff(class(p),'patchwork')),normalize_counts=normalize_counts)))
 lapply(ggplot_build(p)$data,function(z){if(normalize_counts&&'label'%in%names(z))z$label<-sub('^n=','',z$label);z[,intersect(names(z),c('x','y','xmin','xmax','ymin','ymax','xend','yend','group','PANEL','label','colour','fill','shape','linetype','linewidth','alpha')),drop=FALSE]})
}
heights<-c(Main_Figure_3=6.2,Main_Figure_4=3.6,Main_Figure_6=7.8,
 Supplementary_Figure_S4=6.1,Supplementary_Figure_S5=6.5,Supplementary_Figure_S6=4,
 Supplementary_Figure_S7=7.0,Supplementary_Figure_S8=4.1,Supplementary_Figure_S9=8.3,
 Supplementary_Figure_S10=3.4,Supplementary_Figure_S11=6.4,Supplementary_Figure_S12=8.2,
 Supplementary_Figure_S13=8.4,Supplementary_Figure_S14=7.4,Supplementary_Figure_S15=7.5,
 Supplementary_Figure_S16=3.8,Supplementary_Figure_S17=5.4)
checks<-list()
for(name in names(new)){
 a<-signature(old[[name]],normalize_counts=name=='Supplementary_Figure_S15');b<-signature(new[[name]],normalize_counts=name=='Supplementary_Figure_S15')
 same<-isTRUE(all.equal(a,b,tolerance=0,check.attributes=TRUE));if(!same){print(all.equal(a,b,tolerance=0));stop(name,' data signature changed')}
 checks[[name]]<-data.frame(figure=name,numerical_and_aesthetic_layer_signature_identical=same)
 saveRDS(b,file.path(out,'audit',paste0(name,'_verified_data_signature.rds')))
 for(ext in c('png','pdf','svg')){
  kw<-list(filename=file.path(out,'figures',paste0(name,'.',ext)),plot=new[[name]],width=7.0866,height=heights[[name]],dpi=600,bg='white',limitsize=FALSE)
  if(ext=='pdf'){kw$device<-grDevices::pdf;kw$useDingbats<-FALSE};if(ext=='svg')kw$device<-svglite::svglite
  do.call(ggsave,kw)
 }
 cat('Preserved and saved',name,'\n')
}
write.csv(do.call(rbind,checks),file.path(out,'audit/All_figure_layer_preservation.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'audit/R_sessionInfo.txt'))
writeLines(c('PASS','17 layouts checked against frozen plotted values','No points, intervals, summaries, colours, fill, shapes or line data changed','Taxon abbreviations and storage-axis abbreviations change display labels only','No original FCS, microscopy or numerical source files changed','No model refitting or new statistical tests'),file.path(out,'audit/Presentation_preservation.txt'))
unlink(work,recursive=TRUE);cat('PASS: all 17 figure layouts saved.\n')
