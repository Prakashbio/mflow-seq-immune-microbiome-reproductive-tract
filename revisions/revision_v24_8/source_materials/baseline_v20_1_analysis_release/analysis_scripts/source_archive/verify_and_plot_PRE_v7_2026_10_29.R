args<-commandArgs(trailingOnly=TRUE);o<-if(length(args))args[1] else 'review_response/03_analysis/PRE_baseline_review_2026-10-28'
suppressPackageStartupMessages({library(vegan);library(ggplot2);library(grid)})
m<-read.csv(file.path(o,'data/source_metadata_original_codes.csv'));tx<-read.csv(file.path(o,'data/original_taxonomy.csv'));z<-m[m$Fraction=='PRE',]
lib<-read.csv(file.path(o,'results/PRE_library_metrics.csv'));u<-read.csv(file.path(o,'results/PRE_woman_site_metrics.csv'));sp<-read.csv(file.path(o,'results/PRE_woman_site_taxon_relative.csv'));gen<-read.csv(file.path(o,'results/PRE_woman_site_genus_relative.csv'));sums<-read.csv(file.path(o,'results/PRE_site_descriptive_summary.csv'));sens<-read.csv(file.path(o,'results/PRE_participant99_sensitivity.csv'))
eq<-function(a,b)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-11,check.attributes=FALSE)))
for(arm in c('baseline_39','sensitivity_41')){
 x<-as.matrix(read.csv(file.path(o,'data',paste0(arm,'_source_counts.csv')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$SourceCode))
 mat<-x[,z$SourceCode];rel<-sweep(mat,2,colSums(mat),'/');H<-diversity(t(mat),index='shannon')
 ll<-lib[lib$taxon_arm==arm,];ll<-ll[match(z$SourceCode,ll$SourceCode),];eq(H,ll$Shannon);eq(colSums(mat),ll$sum_retained_taxon_counts)
 uu<-u[u$taxon_arm==arm,]
 for(i in seq_len(nrow(uu))){
  row<-uu[i,];ids<-z$SourceCode[z$Participant==row$Participant&z$Location==row$Location];avg<-rowMeans(rel[,ids,drop=FALSE])
  eq(mean(H[ids]),row$mean_preparation_Shannon);eq(diversity(avg),row$Shannon_of_mean_profile);eq(length(ids),row$n_PRE_libraries)
  ss<-sp[sp$taxon_arm==arm&sp$Participant==row$Participant&sp$Location==row$Location,];eq(avg,ss$relative_abundance[match(names(avg),ss$Taxon)])
  g<-rowsum(matrix(avg,ncol=1),tx$Genus[match(rownames(x),tx$Taxon)]);gg<-gen[gen$taxon_arm==arm&gen$Participant==row$Participant&gen$Location==row$Location,];eq(g[,1],gg$relative_abundance[match(rownames(g),gg$Genus)])
 }
 for(site in c('VAGINA','CERVIX','ENDOMETRIUM')){
  v<-uu[uu$Location==site,];ss<-sums[sums$taxon_arm==arm&sums$Location==site,]
  eq(median(v$mean_preparation_Shannon),ss$Shannon_median);eq(quantile(v$mean_preparation_Shannon,c(.25,.75)),c(ss$Shannon_Q1,ss$Shannon_Q3));eq(nrow(v),ss$n_women);eq(sum(v$n_PRE_libraries),ss$n_PRE_libraries)
 }
 for(choice in c('both_preparations_equal_woman','retain99EP2','retain99EP')){
  ids<-z$SourceCode[z$Location=='ENDOMETRIUM'];if(choice=='retain99EP2')ids<-setdiff(ids,'99EP');if(choice=='retain99EP')ids<-setdiff(ids,'99EP2')
  zz<-z[match(ids,z$SourceCode),];v<-tapply(H[ids],zz$Participant,mean);ss<-sens[sens$taxon_arm==arm&sens$choice==choice,]
  eq(c(median(v),quantile(v,c(.25,.75)),mean(v)),c(ss$Shannon_median,ss$Shannon_Q1,ss$Shannon_Q3,ss$Shannon_mean))
 }
}
writeLines('PASS: independent vegan/R verification of all74 PRE library metrics,72 woman/site rows,all species/genus profiles,6 site summaries and6 participant99 sensitivities. No hypothesis tests performed.',file.path(o,'results/Independent_R_verification.txt'))

people<-sort(unique(z$Participant));sites<-c('VAGINA','CERVIX','ENDOMETRIUM');g<-gen[gen$taxon_arm=='baseline_39',];uu<-u[u$taxon_arm=='baseline_39',]
# Rank genera by their total across equally weighted observed woman/site profiles.
top<-names(sort(tapply(g$relative_abundance,g$Genus,sum),decreasing=TRUE))[1:6]
g$Group<-ifelse(g$Genus%in%top,g$Genus,'Other retained genera');pdat<-aggregate(relative_abundance~Participant+Location+Group,g,sum)
pdat$Group<-factor(pdat$Group,levels=c(top,'Other retained genera'));pdat$Participant<-factor(pdat$Participant,levels=people);pdat$Location<-factor(pdat$Location,levels=sites)
missing<-read.csv(file.path(o,'results/PRE_missing_woman_site_cells.csv'));missing$Participant<-factor(missing$Participant,levels=people);missing$Location<-factor(missing$Location,levels=sites)
site_labels<-c('VAGINA'='Vagina','CERVIX'='Cervix','ENDOMETRIUM'='Endometrium')
cols<-setNames(c('#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB'),levels(pdat$Group))
theme<-theme_bw(base_size=11)+theme(axis.text=element_text(size=10.5),legend.text=element_text(size=10.5),legend.position='bottom',strip.text=element_text(size=10.5),panel.spacing.y=unit(.18,'in'),plot.title=element_text(size=12,face='bold'))
ga<-ggplot(pdat,aes(Participant,relative_abundance,fill=Group))+geom_col(width=.86)+geom_rect(data=missing,aes(xmin=as.numeric(Participant)-.43,xmax=as.numeric(Participant)+.43,ymin=0,ymax=1),inherit.aes=FALSE,fill='grey95',color='grey60')+geom_text(data=missing,aes(Participant,.5,label='NA'),inherit.aes=FALSE,size=3.8,angle=90,color='grey35')+facet_grid(Location~.,labeller=labeller(Location=site_labels))+scale_fill_manual(values=cols,drop=FALSE)+scale_y_continuous(labels=function(x)paste0(round(x*100),'%'),expand=expansion(mult=c(0,.02)))+coord_cartesian(ylim=c(0,1))+theme+guides(fill=guide_legend(nrow=3,byrow=TRUE))+labs(title='A  PRE composition',x='Original participant code',y='Relative abundance',fill=NULL)
# Check the rendered stack inputs as well as the source profiles. Coordinate zoom
# avoids dropping near-one cumulative sums because of floating-point rounding.
plotted<-ggplot_build(ga)$data[[1]];stopifnot(all(is.finite(plotted$ymax)))
stacks<-aggregate(ymax~PANEL+x,plotted,max);eq(stacks$ymax,rep(1,36));stopifnot(nrow(stacks)==36)
uu$Location<-factor(uu$Location,levels=sites)
gb<-ggplot(uu,aes(Location,mean_preparation_Shannon,group=Participant))+geom_line(color='grey65',linewidth=.4,alpha=.8)+geom_point(aes(color=Location),size=2.5)+scale_color_manual(values=c('#4477AA','#228833','#AA3377'),labels=c('Vagina','Cervix','Endometrium'))+scale_x_discrete(labels=c('Vagina','Cervix','Endometrium'))+theme+theme(legend.position='none')+labs(title='B  PRE diversity',x=NULL,y='Shannon H')
draw<-function(){grid.newpage();pushViewport(viewport(layout=grid.layout(2,1,heights=unit(c(.72,.28),'npc'))));print(ga,vp=viewport(layout.pos.row=1));print(gb,vp=viewport(layout.pos.row=2));popViewport()}
png(file.path(o,'results/PRE_baseline_context_REVIEW.png'),width=2840,height=2720,res=400);draw();dev.off()
pdf(file.path(o,'results/PRE_baseline_context_REVIEW.pdf'),width=7.1,height=6.8);draw();dev.off()
write.csv(data.frame(legend_order=seq_along(top),Genus=top,selection_basis='Total over36 equally weighted observed woman/site profiles in39-taxon baseline'),file.path(o,'results/PRE_plot_genus_selection.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(o,'R_sessionInfo.txt'));cat('Independent verification and descriptive PRE figure saved.\n')
