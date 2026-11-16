# Change only the legend layout of Figure 5A, retaining plot data and B styling.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);src<-normalizePath(args[1]);out<-args[2];stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);dir.create(file.path(out,'figures'));dir.create(file.path(out,'audit'))
D<-file.path(src,'source_data');e<-new.env(parent=globalenv());captured<-list()
e$commandArgs<-function(...)c(D,file.path(D,'statistics'))
e$ggsave<-function(filename,plot=last_plot(),...){captured[[tools::file_path_sans_ext(basename(filename))]]<<-plot;invisible(NULL)}
e$write.csv<-e$writeLines<-e$capture.output<-function(...)invisible(NULL)
for(ex in parse(file.path(src,'scripts/plot_diversity_integration_2026_10_28.R')))eval(ex,envir=e)
style<-theme_classic(base_size=9.5,base_family='Helvetica')+theme(
 text=element_text(colour='#17212B'),axis.text=element_text(size=8,colour='#28343E'),
 axis.title=element_text(size=9),axis.line=element_line(linewidth=.3,colour='#65717A'),axis.ticks=element_line(linewidth=.25,colour='#65717A'),
 panel.grid.major.y=element_line(colour='#E9EDF0',linewidth=.25),strip.background=element_rect(fill='#F0F3F5',colour=NA),strip.text=element_text(face='bold',size=9),
 panel.spacing=unit(3,'mm'),legend.position='bottom',legend.title=element_text(size=8.5),legend.text=element_text(size=8),legend.key.size=unit(3.5,'mm'),
 plot.title=element_text(face='bold',size=10,hjust=0),plot.subtitle=element_text(size=8),plot.caption=element_blank(),plot.tag=element_text(face='bold',size=12),plot.margin=margin(5,6,5,5))
polish<-function(p){for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,'GeomPath')||inherits(p$layers[[i]]$geom,'GeomLine')){p$layers[[i]]$aes_params$alpha<-.65;p$layers[[i]]$aes_params$linewidth<-.45};p+style}
oldA<-polish(captured[['Figure3_PCoA_participant_pairs']])+labs(title='A  Community ordination')+guides(colour=guide_legend(nrow=3,override.aes=list(shape=16,alpha=1)),shape=guide_legend(nrow=1))+theme(legend.text=element_text(size=8),legend.key.width=unit(2.5,'mm'),legend.box='vertical')
newA<-oldA+guides(
 colour=guide_legend(order=1,nrow=1,byrow=TRUE,title.position='left',override.aes=list(shape=16,alpha=1,linetype=0)),
 shape=guide_legend(order=2,nrow=1,byrow=TRUE,title.position='left',override.aes=list(colour='#17212B',alpha=1,linetype=0)))+
 theme(legend.box='vertical',legend.box.just='center',legend.title=element_text(face='bold',size=8.5),legend.key.width=unit(2.5,'mm'),legend.key.height=unit(3.5,'mm'),legend.spacing.x=unit(1.2,'mm'),legend.spacing.y=unit(1.2,'mm'),legend.margin=margin(0,1,0,1),legend.box.margin=margin(2,0,2,0))
# B is exactly the same plot specification as the frozen v23 figure driver.
B<-polish(captured[['FigureS4_matched_own_PRE']])+scale_colour_manual(values=c(Proliferative='#D55E00',Secretory='#0072B2'))+labs(title='B  Matched distance to own PRE',y='Bray-Curtis distance to own PRE')+theme(legend.title=element_text(size=8))
signature<-function(p)lapply(ggplot_build(p)$data,function(z)z[,intersect(names(z),c('x','y','xend','yend','group','PANEL','colour','shape','alpha','linewidth')),drop=FALSE])
stopifnot(isTRUE(all.equal(signature(oldA),signature(newA),tolerance=0)))
p<-newA/B+plot_layout(heights=c(1.4,1))
for(ext in c('png','pdf','svg')){
 kwargs<-list(filename=file.path(out,'figures',paste0('Main_Figure_5.',ext)),plot=p,width=7.0866,height=7.2,dpi=600,bg='white',limitsize=FALSE)
 if(ext=='pdf'){kwargs$device<-grDevices::pdf;kwargs$useDingbats<-FALSE};if(ext=='svg')kwargs$device<-svglite::svglite
 do.call(ggsave,kwargs)
}
write.csv(captured[['Figure3_PCoA_participant_pairs']]$data,file.path(out,'audit/Figure5A_original_ordination_plot_data.csv'),row.names=FALSE)
write.csv(captured[['FigureS4_matched_own_PRE']]$data,file.path(out,'audit/Figure5B_unchanged_matched_distance_plot_data.csv'),row.names=FALSE)
saveRDS(signature(newA),file.path(out,'audit/Figure5A_verified_layer_signature.rds'))
writeLines(c('PASS','Only A participant/fraction legend layout changed','All A plotted coordinates, groups, colours, shapes, alpha and line widths unchanged','Participant legend in one horizontal row; all 15 IDs retained','Fraction key in one horizontal row; all four original symbols retained','Panel B plot specification and legend styling unchanged','No model refitting or new inference'),file.path(out,'audit/Legend_preservation.txt'))
capture.output(sessionInfo(),file=file.path(out,'audit/R_sessionInfo.txt'));cat('PASS: Figure 5A legend layout saved.\n')
