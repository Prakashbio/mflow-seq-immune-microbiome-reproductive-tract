args<-commandArgs(TRUE);O<-if(length(args))args[1] else 'review_response/supplement_audit_working_v6_2026-10-29'
I<-file.path(O,'supplement_audit_inputs');D<-file.path(O,'supplement_audit_statistics');F<-file.path(O,'figures');dir.create(F,showWarnings=FALSE,recursive=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
rd<-function(n)read.csv(file.path(D,paste0(n,'.csv')),check.names=FALSE,stringsAsFactors=FALSE)
m<-read.csv(file.path(I,'Metadata_original_codes.csv'),colClasses='character',check.names=FALSE)
sites<-c('VAGINA','CERVIX','ENDOMETRIUM');colors<-c(VAGINA='#0072B2',CERVIX='#009E73',ENDOMETRIUM='#CC79A7')
ps<-rd('PRE_equal_woman_profiles');hs<-rd('PRE_participant_Shannon');li<-rd('Library_assignment_count_metrics');co<-rd('PRE_correlation_all_pairs');ac<-rd('PRE_exact_accumulation');di<-rd('PRE_cross_site_distances');key<-rd('S6_display_taxon_key')
same<-function(x,y)stopifnot(length(x)==length(y),all(is.na(x)==is.na(y)),all(abs(x[!is.na(x)]-y[!is.na(y)])<1e-10))
checks<-c(profiles=0,library_metrics=0,correlations=0,accumulation=0,distances=0)
for(arm in c(39,41)){
 x<-as.matrix(read.csv(file.path(I,paste0('Counts_',arm,'_original_codes.csv')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$SourceCode));rel<-sweep(x,2,colSums(x),'/')
 sh<-apply(rel,2,function(v)-sum(v[v>0]*log(v[v>0])))
 for(j in seq_len(nrow(m))){z<-li[li$arm==arm&li$SourceCode==m$SourceCode[j],];same(z$retained_assignment_counts,sum(x[,j]));same(z$Shannon,sh[j]);same(z$observed_assignments,sum(x[,j]>0));checks['library_metrics']<-checks['library_metrics']+1}
 pre<-m[m$Fraction=='PRE',];profiles<-list()
 for(s in sites){
  people<-sort(unique(pre$Participant[pre$Location==s]));mats<-sapply(people,function(p){ids<-pre$SourceCode[pre$Location==s&pre$Participant==p];rowMeans(rel[,ids,drop=FALSE])});colnames(mats)<-people;profiles[[s]]<-mats
  for(p in people){z<-ps[ps$arm==arm&ps$Location==s&ps$Participant==p,];same(z$mean_relative_abundance,mats[match(z$Taxon,rownames(mats)),p]);hh<-hs[hs$arm==arm&hs$Location==s&hs$Participant==p,];ids<-pre$SourceCode[pre$Location==s&pre$Participant==p];same(hh$mean_preparation_Shannon,mean(sh[ids]));checks['profiles']<-checks['profiles']+nrow(z)}
  zz<-co[co$arm==arm&co$Location==s,];expected<-suppressWarnings(cor(t(mats),method='spearman'));constant<-apply(mats,1,function(v)length(unique(v))<2);expected[constant,]<-NA_real_;expected[,constant]<-NA_real_
  same(zz$rho,expected[cbind(match(zz$Taxon1,rownames(expected)),match(zz$Taxon2,colnames(expected)))]);checks['correlations']<-checks['correlations']+nrow(zz)
  incidence<-rowSums(mats>0);n<-ncol(mats)
  for(k in seq_len(n)){z<-ac[ac$arm==arm&ac$Location==s&ac$women_accumulated==k,];same(z$expected_observed_assignments,sum(1-choose(n-incidence,k)/choose(n,k)));checks['accumulation']<-checks['accumulation']+1}
 }
 z<-di[di$arm==arm,]
 for(j in seq_len(nrow(z))){q<-z[j,];u<-profiles[[q$Site1]][,as.character(q$Participant1)];v<-profiles[[q$Site2]][,as.character(q$Participant2)];same(q$Bray_Curtis,sum(abs(u-v))/sum(u+v));checks['distances']<-checks['distances']+1}
}
# Exhaustive subset check of the exact accumulation formula, independently of
# the incidence formula, for the small cervical PRE dataset in both arms.
for(arm in c(39,41)){
 z<-ps[ps$arm==arm&ps$Location=='CERVIX',];mat<-xtabs(mean_relative_abundance~Taxon+Participant,z)>0;n<-ncol(mat)
 for(k in seq_len(n)){
  actual<-mean(combn(seq_len(n),k,FUN=function(j)sum(rowSums(mat[,j,drop=FALSE])>0)))
  same(actual,ac$expected_observed_assignments[ac$arm==arm&ac$Location=='CERVIX'&ac$women_accumulated==k])
 }
}
ranked<-aggregate(mean_relative_abundance~Location+Taxon,ps[ps$arm==39,],mean);ranked<-aggregate(mean_relative_abundance~Taxon,ranked,mean);ranked<-ranked[order(-ranked$mean_relative_abundance,ranked$Taxon),];stopifnot(identical(key$Taxon,head(ranked$Taxon,12)))
writeLines(c('PASS: independent R reconstruction from frozen original-code counts.',paste(names(checks),checks,sep=': '),'Cervical exact accumulation independently checked by exhaustive woman subsets in both memberships.','Rank correlations have no inferential p values; constant inputs remain missing.','No FACS-event counts, raw-read QC status or no-depth-confounding claim is produced.'),file.path(D,'Independent_R_supplement_audit_verification.txt'))
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=11),legend.position='bottom'))
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,device='pdf',useDingbats=FALSE,bg='white')}
plots<-list()
for(s in sites){
 z<-co[co$arm==39&co$Location==s&co$Taxon1%in%key$Taxon&co$Taxon2%in%key$Taxon,];z$x<-match(z$Taxon2,key$Taxon);z$y<-match(z$Taxon1,key$Taxon)
 n<-unique(z$n_women);plots[[s]]<-ggplot(z,aes(x,y,fill=rho))+geom_tile()+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),na.value='#BBBBBB',name='Rank correlation')+scale_x_continuous(breaks=1:12,expand=c(0,0))+scale_y_reverse(breaks=1:12,expand=c(0,0))+coord_equal()+labs(title=paste0(tools::toTitleCase(tolower(s)),' PRE: ',n,' women'),x=NULL,y=NULL)+theme(panel.grid=element_blank(),axis.text=element_text(size=8),plot.title=element_text(size=10))
}
key$label<-paste0(key$display_ID,'  ',key$Taxon)
legend<-ggplot(key,aes(0,13-display_ID,label=label))+geom_text(hjust=0,size=3.4,fontface='italic')+annotate('text',0,13,label='Displayed assignments',hjust=0,fontface='bold',size=3.5)+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(.3,13.5))+theme_void()+theme(plot.margin=margin(5,0,0,0))
p6<-(plots[[1]]+plots[[2]]+plots[[3]]+legend)+plot_layout(ncol=2,guides='collect')&theme(legend.position='bottom')
savefig(p6,'FigureS6_PRE_descriptive_rank_correlations',7.1,6.5)
ac$Location<-factor(ac$Location,levels=sites);ac$membership<-factor(ac$arm,levels=c(39,41),labels=c('39 assignments','41 assignments'))
p8<-ggplot(ac,aes(women_accumulated,expected_observed_assignments,color=Location,linetype=membership))+geom_line(linewidth=.8)+geom_point(size=1.4)+scale_color_manual(values=colors,labels=c('Vagina','Cervix','Endometrium'))+scale_linetype_manual(values=c('solid','dashed'))+scale_x_continuous(breaks=1:15)+labs(title='Observed PRE assignments accumulate across women',x='Women included at the stated site',y='Expected observed assignments',color='Site',linetype='Membership')+theme(legend.box='vertical',legend.margin=margin(0,0,0,0))
savefig(p8,'FigureS8_PRE_participant_accumulation',7.1,4.3)
base<-hs[hs$arm==39,];all3<-Reduce(intersect,lapply(sites,function(s)base$Participant[base$Location==s]));base<-base[base$Participant%in%all3,];base$Location<-factor(base$Location,levels=sites)
pA<-ggplot(base,aes(Location,mean_preparation_Shannon,group=Participant,color=factor(Participant)))+geom_line(linewidth=.45,show.legend=FALSE)+geom_point(size=1.7,show.legend=FALSE)+scale_color_viridis_d(end=.85)+labs(title=paste0('A  PRE trajectories: ',length(all3),' women'),x=NULL,y='Mean preparation Shannon')+scale_x_discrete(labels=c('Vagina','Cervix','Endometrium'))+theme(legend.position='none',axis.text.x=element_text(size=8))
l<-li[li$arm==39,];l$Fraction<-factor(l$Fraction,levels=c('PRE','IgA','IgM','IgG'));l$Location<-factor(l$Location,levels=sites)
pB<-ggplot(l,aes(Fraction,retained_assignment_counts,color=Location))+geom_boxplot(aes(group=Fraction),color='#444444',fill=NA,outlier.shape=NA,width=.6)+geom_point(position=position_jitter(width=.13,height=0,seed=20261029),size=1.4,alpha=.65)+scale_color_manual(values=colors,labels=c('Vagina','Cervix','Endometrium'))+scale_y_log10()+labs(title='B  Retained counts: 149 libraries',x=NULL,y='39-assignment sum (log scale)',color=NULL)+theme(legend.position='none')
d<-di[di$arm==39,];d$site_pair<-paste(substr(d$Site1,1,1),substr(d$Site2,1,1),sep='-');d$site_pair<-factor(d$site_pair,levels=c('V-C','V-E','C-E'));d$pair_type<-factor(d$pair_type,levels=c('Within woman','Between women'))
pC<-ggplot(d,aes(site_pair,Bray_Curtis,fill=pair_type))+geom_boxplot(width=.65,outlier.size=.6)+scale_fill_manual(values=c('#56B4E9','#E69F00'))+labs(title='C  PRE distances are descriptive',x='Site pair (V vagina; C cervix; E endometrium)',y='Bray-Curtis on relative abundance',fill=NULL)+theme(axis.title.x=element_text(size=8),legend.text=element_text(size=8),legend.key.width=unit(.3,'cm'),legend.spacing.x=unit(.1,'cm'))
pD<-ggplot(l,aes(retained_assignment_counts,Shannon,color=Location))+geom_point(size=1.5,alpha=.7)+scale_x_log10()+scale_color_manual(values=colors,labels=c('Vagina','Cervix','Endometrium'))+labs(title='D  Count-diversity display',x='39-assignment sum (log scale)',y='Shannon diversity',color=NULL)+guides(color=guide_legend(nrow=2))+theme(legend.text=element_text(size=8),legend.key.width=unit(.3,'cm'))
p11<-(pA+pB)/(pC+pD)+plot_layout(guides='keep')&theme(plot.title=element_text(size=10),axis.text=element_text(size=10.2),axis.title=element_text(size=10.2),legend.text=element_text(size=10.2))
savefig(p11,'FigureS11_PRE_and_retained_count_descriptions',7.1,7.1)
capture.output(sessionInfo(),file=file.path(D,'R_sessionInfo.txt'))
cat('Independent reconstruction PASS; three PNG/PDF figure pairs saved.\n')
