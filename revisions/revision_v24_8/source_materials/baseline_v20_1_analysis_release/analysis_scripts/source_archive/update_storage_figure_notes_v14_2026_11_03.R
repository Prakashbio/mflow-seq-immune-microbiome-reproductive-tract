# Update only obsolete figure-note wording after confirmed same-biopsy pairing.
# Read saved source coordinates and measured values; perform no new analysis.
suppressPackageStartupMessages(library(ggplot2))
a<-commandArgs(TRUE);stopifnot(length(a)==1)
p<-a[1];d<-file.path(p,'storage_statistics');f<-file.path(p,'figure_note_update_v14');dir.create(f,showWarnings=FALSE)
cols<-c(Fresh='#D97B31',Frozen='#3A79B8')
theme<-theme_bw(base_size=14)+theme(legend.position='bottom',panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey94'),plot.caption=element_text(size=11,hjust=0))
co<-read.csv(file.path(d,'FigureS15_PRE_PCoA_coordinates.csv'))
co$Location<-factor(co$Location,levels=c('VAGINA','CERVIX','ENDOMETRIUM'),labels=c('Vagina','Cervix','Endometrium'))
co$Fraction<-factor(co$Fraction,levels=c('PRE','IgA','IgM','IgG'))
co$Storage<-factor(co$Storage,levels=c('FRESH','FROZEN'),labels=c('Fresh','Frozen'))
pa<-read.csv(file.path(d,'FigureS15_PRE_PCoA_parameters.csv'))
ev<-read.csv(file.path(d,'FigureS15_PRE_PCoA_eigenvalues.csv'))
pc<-100*ev$eigenvalue[1:2]/pa$positive_eigenvalue_sum
pt99<-co[co$Participant==99&co$Location=='Endometrium',]
p10<-ggplot(co,aes(Axis1,Axis2,colour=Storage))+geom_point(size=2.2,alpha=.85)+geom_path(data=pt99,aes(group=Participant),linetype=2,colour='grey50',linewidth=.5)+geom_text(data=pt99,aes(label=SourceCode),size=3,nudge_y=.025,show.legend=FALSE)+facet_wrap(~Location,nrow=1)+scale_colour_manual(values=cols)+theme+labs(x=sprintf('PRE PCoA1 (%.1f%% of positive corrected eigenvalues)',pc[1]),y=sprintf('PRE PCoA2 (%.1f%%)',pc[2]),colour='Recorded storage',caption="One common Bray-Curtis PCoA of 37 PRE libraries, with Lingoes correction, displayed on common axes across sites.\nThe dashed line joins paired fresh and cryopreserved fragments from woman 99's same biopsy. No storage test or confidence ellipse.")
pair<-read.csv(file.path(d,'FigureS19_participant99_plotted_values.csv'))
pair$Fraction<-factor(pair$Fraction,levels=c('PRE','IgA','IgM','IgG'))
pair$Storage<-factor(pair$Storage,levels=c('Fresh','Frozen'))
pair$metric<-factor(pair$metric,levels=c('Shannon H','Lactobacillus abundance (%)'))
p11<-ggplot(pair,aes(Storage,value,group=1))+geom_line(colour='grey50',linetype=2,linewidth=.5)+geom_point(aes(colour=Storage),size=2.5)+facet_grid(metric~Fraction,scales='free_y')+scale_colour_manual(values=cols)+theme+theme(legend.position='none',strip.text.y=element_text(size=11))+labs(x='Recorded preparation condition',y=NULL,caption='Participant 99 endometrium only: one woman, two preparations, four available fractions per preparation.\nLines join paired fresh and cryopreserved fragments of the same biopsy. Fractions are not independent storage replicates. No equivalence test.')
for(n in c('S10','S11')){
 plot<-if(n=='S10')p10 else p11;ht<-if(n=='S10')5.2 else 6.8
 ggsave(file.path(f,paste0('Supplementary_Figure_',n,'.png')),plot,width=12,height=ht,dpi=300,bg='white')
 ggsave(file.path(f,paste0('Supplementary_Figure_',n,'.pdf')),plot,width=12,height=ht,bg='white')
}
capture.output(sessionInfo(),file=file.path(f,'R_sessionInfo.txt'))
writeLines('Only two storage figure notes are updated. All saved coordinates, measured source values and input statistics are read-only. No microscopy image is opened or edited.',file.path(f,'SCOPE.txt'))
