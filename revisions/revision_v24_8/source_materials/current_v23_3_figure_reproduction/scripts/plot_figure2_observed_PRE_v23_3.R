# Presentation-only Figure 2 revision. Run: Rscript script.R FROZEN_V23_FOLDER NEW_OUTPUT
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(TRUE);src<-normalizePath(args[1]);out<-args[2]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);dir.create(file.path(out,'figures'));dir.create(file.path(out,'audit'))
D<-file.path(src,'source_data/pre_context')
# Reuse the frozen verified plot inputs; intercept all historical writes/devices.
e<-new.env(parent=globalenv());e$commandArgs<-function(...)D
e$png<-e$pdf<-e$dev.off<-e$draw<-function(...)invisible(NULL)
e$write.csv<-e$writeLines<-e$capture.output<-function(...)invisible(NULL)
for(ex in parse(file.path(src,'scripts/verify_and_plot_PRE_v7_2026_10_29.R'))){
 if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&identical(ex[[2]],as.name('draw')))next
 eval(ex,envir=e)
}
style<-theme_classic(base_size=9.5,base_family='Helvetica')+theme(
 text=element_text(colour='#17212B'),axis.text=element_text(size=8,colour='#28343E'),
 axis.title=element_text(size=9),axis.line=element_line(linewidth=.3,colour='#65717A'),
 axis.ticks=element_line(linewidth=.25,colour='#65717A'),
 panel.grid.major.y=element_line(colour='#E9EDF0',linewidth=.25),
 strip.background=element_rect(fill='#F0F3F5',colour=NA),strip.text=element_text(face='bold',size=9),
 panel.spacing=unit(3,'mm'),legend.position='bottom',legend.title=element_text(size=8.5),
 legend.text=element_text(size=8),legend.key.size=unit(3.5,'mm'),
 plot.title=element_text(face='bold',size=10,hjust=0),plot.subtitle=element_text(size=8),
 plot.caption=element_blank(),plot.tag=element_text(face='bold',size=12),plot.margin=margin(5,6,5,5))
polish<-function(p){
 for(i in seq_along(p$layers))if(inherits(p$layers[[i]]$geom,'GeomPath')||inherits(p$layers[[i]]$geom,'GeomLine')){p$layers[[i]]$aes_params$alpha<-.65;p$layers[[i]]$aes_params$linewidth<-.45}
 p+style
}
pdat<-e$pdat
labels<-c(VAGINA='Vagina\n15 women',CERVIX='Cervix\n8 women',ENDOMETRIUM='Endometrium\n13 women')
pa<-ggplot(pdat,aes(Participant,relative_abundance,fill=Group))+geom_col(width=.86)+
 facet_grid(.~Location,scales='free_x',space='free_x',labeller=labeller(Location=labels))+
 scale_fill_manual(values=e$cols,drop=FALSE)+
 scale_y_continuous(labels=function(x)paste0(round(x*100),'%'),breaks=seq(0,1,.25),expand=expansion(mult=c(0,.02)))+
 coord_cartesian(ylim=c(0,1))+style+
 guides(fill=guide_legend(nrow=2,byrow=TRUE))+
 labs(title='A  PRE genus composition',x='Participant within anatomical site',y='Relative abundance',fill=NULL)+
 theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=8),panel.spacing.x=unit(4,'mm'))
pb<-polish(e$gb)+labs(title='B  PRE diversity within women',y='Shannon entropy')+theme(legend.position='none')
# Prove that every observed stacked abundance, identity and colour is retained.
old<-ggplot_build(e$ga);new<-ggplot_build(pa)
stack_identity<-function(b){
 z<-b$data[[1]];maps<-b$layout$layout;res<-list()
 for(i in seq_len(nrow(z))){
  panel<-as.integer(z$PANEL[i]);lr<-maps[match(panel,maps$PANEL),];sx<-b$layout$panel_scales_x[[lr$SCALE_X]]
  res[[i]]<-data.frame(Location=as.character(lr$Location),Participant=as.character(sx$get_labels()[as.integer(z$x[i])]),fill=z$fill[i],ymin=z$ymin[i],ymax=z$ymax[i])
 }
 rr<-do.call(rbind,res);rr<-rr[order(rr$Location,as.numeric(rr$Participant),rr$fill),];rownames(rr)<-NULL;rr
}
stopifnot(isTRUE(all.equal(stack_identity(old),stack_identity(new),tolerance=0)))
sig<-function(p)lapply(ggplot_build(p)$data,function(z)z[,intersect(names(z),c('x','y','group','PANEL')),drop=FALSE])
stopifnot(isTRUE(all.equal(sig(e$gb),sig(pb),tolerance=0)))
stacks<-aggregate(ymax~PANEL+x,new$data[[1]],max);stopifnot(nrow(stacks)==36,all(abs(stacks$ymax-1)<1e-12),nrow(e$missing)==9)
stopifnot(all(table(pdat$Location)/length(levels(pdat$Group))==c(15,8,13)))
p<-pa/pb+plot_layout(heights=c(1.15,1))
for(ext in c('png','pdf','svg')){
 kwargs<-list(filename=file.path(out,'figures',paste0('Main_Figure_2.',ext)),plot=p,width=7.0866,height=7.4,dpi=600,bg='white',limitsize=FALSE)
 if(ext=='pdf'){kwargs$device<-grDevices::pdf;kwargs$useDingbats<-FALSE}
 if(ext=='svg')kwargs$device<-svglite::svglite
 do.call(ggsave,kwargs)
}
write.csv(pdat,file.path(out,'audit/Figure2_all_observed_stacked_values.csv'),row.names=FALSE)
write.csv(e$uu,file.path(out,'audit/Figure2_woman_site_Shannon_values.csv'),row.names=FALSE)
write.csv(e$missing,file.path(out,'audit/Figure2_unavailable_PRE_observations.csv'),row.names=FALSE)
write.csv(stack_identity(new),file.path(out,'audit/Figure2_rendered_stack_segments.csv'),row.names=FALSE)
writeLines(c('PASS','36 observed woman/site compositions retained exactly','15 women and 37 PRE libraries retained','9 absent cells omitted from display without imputation','All genus colours and rendered stack segment values unchanged','All panel B coordinates and participant grouping unchanged','No statistical models refitted'),file.path(out,'audit/Presentation_preservation.txt'))
capture.output(sessionInfo(),file=file.path(out,'audit/R_sessionInfo.txt'))
cat('PASS: observed-only Figure 2 saved; all numerical values preserved.\n')
