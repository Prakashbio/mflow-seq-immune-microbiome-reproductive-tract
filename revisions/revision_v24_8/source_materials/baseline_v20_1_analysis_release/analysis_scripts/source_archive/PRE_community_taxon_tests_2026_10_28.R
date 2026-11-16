# One original PRE library per woman/site; each HUT099 choice kept separately.
suppressPackageStartupMessages(library(ggplot2))
o<-'review_response/03_analysis/Taxon_scores_PRE_review_2026-10-28'
m<-read.csv(file.path(o,'data/source_metadata_with_vaginal_category.csv'),stringsAsFactors=FALSE);tx<-read.csv(file.path(o,'data/original_taxonomy.csv'),stringsAsFactors=FALSE)
put<-function(x,n)write.csv(x,file.path(o,'results',paste0(n,'.csv')),row.names=FALSE)
B<-5000L;rows<-list();units<-list();prevalues<-list();allocations<-list()
for(arm in c('baseline_39','sensitivity_41')){
 x<-as.matrix(read.csv(file.path(o,'data',paste0(arm,'_source_counts.csv')),row.names=1,check.names=FALSE));stopifnot(identical(colnames(x),m$SourceCode))
 genus<-tx$Genus[match(rownames(x),tx$Taxon)];stopifnot(!anyNA(genus))
 for(rank in c('provisional_species','genus_aggregation')){
  mat<-if(rank=='provisional_species')x else rowsum(x,genus,reorder=TRUE);rel<-sweep(mat,2,colSums(mat),'/')
  for(choice in c('retain99EP2','retain99EP')){
   drop<-if(choice=='retain99EP2')'99EP' else '99EP2';z<-m[m$Fraction=='PRE'&m$SourceCode!=drop,]
   stopifnot(nrow(z)==36,!anyDuplicated(paste(z$Participant,z$Location)))
   units[[paste(arm,rank,choice)]]<-cbind(taxon_arm=arm,rank=rank,selection=choice,z)
   for(site in c('VAGINA','CERVIX','ENDOMETRIUM')){
    zz<-z[z$Location==site,];ids<-zz$SourceCode;A<-zz$VaginalCategory=='VAGINAL_DYSBIOSIS_DESCRIPTION';nD<-sum(A);nL<-sum(!A);n<-nrow(zz);stopifnot(nD>1,nL>1)
    C<-combn(n,nD);W<-matrix(0,ncol(C),n);W[cbind(rep(seq_len(ncol(C)),each=nD),as.vector(C))]<-1
    allocations[[paste(site,choice)]]<-list(SourceCode=ids,observed_dysbiosis=A,assignment_matrix=W)
    set.seed(20261101L+match(site,c('VAGINA','CERVIX','ENDOMETRIUM')))
    bd<-matrix(sample.int(nD,B*nD,replace=TRUE),B,nD);bl<-matrix(sample.int(nL,B*nL,replace=TRUE),B,nL)
    draws<-cbind(bd,bl);colnames(draws)<-c(paste0('D',seq_len(nD)),paste0('L',seq_len(nL)))
    write.csv(draws,file.path(o,'results',paste0('PRE_bootstrap_indices_',site,'_',choice,'.csv')),row.names=FALSE)
    for(tax in rownames(rel)){
     values<-rel[tax,ids];a<-values[A];b<-values[!A];rr<-rank(values,ties.method='average');U<-sum(rr[A])-nD*(nD+1)/2
     constant<-all(values==values[1]);permU<-as.vector(W%*%rr)-nD*(nD+1)/2
     p<-if(constant)NA_real_ else mean(abs(permU-nD*nL/2)>=abs(U-nD*nL/2)-1e-12)
     boot<-apply(matrix(a[bd],B,nD),1,median)-apply(matrix(b[bl],B,nL),1,median);ci<-quantile(boot,c(.025,.975),names=FALSE)
     approx<-if(constant)NA_real_ else suppressWarnings(wilcox.test(a,b,exact=FALSE,correct=TRUE)$p.value)
     rows[[length(rows)+1]]<-data.frame(taxon_arm=arm,rank=rank,selection=choice,Location=site,Taxon=tax,n_dysbiosis_description=nD,n_Lactobacillus_description=nL,detected_D=sum(a>0),detected_L=sum(b>0),median_D=median(a),median_L=median(b),median_difference_D_minus_L=median(a)-median(b),CI95_low=ci[1],CI95_high=ci[2],rank_biserial=if(constant)NA_real_ else 2*U/(nD*nL)-1,U_D=U,exact_assignments=nrow(W),p_exact_label_permutation=p,p_asymptotic_Wilcoxon_sensitivity=approx,status=if(constant)'Non-estimable constant observed values' else 'Exploratory exact tied-rank label-permutation test',sparse_detection=sum(a>0)<3||sum(b>0)<3,interpretation='Vaginal-community descriptions, not clinical diagnoses; source-linked observational PRE comparison; no evidence of ascension or independent validation.')
     prevalues[[length(prevalues)+1]]<-data.frame(taxon_arm=arm,rank=rank,selection=choice,Location=site,Taxon=tax,SourceCode=ids,Participant=zz$Participant,Category=zz$VaginalCategory,Phase=zz$Phase,relative_abundance=values)
    }
   }
  }
 }
}
res<-do.call(rbind,rows);res$p_BH_within_site<-ave(res$p_exact_label_permutation,interaction(res$taxon_arm,res$rank,res$selection,res$Location),FUN=function(p)p.adjust(p,'BH'))
res$p_BH_all_three_sites<-ave(res$p_exact_label_permutation,interaction(res$taxon_arm,res$rank,res$selection),FUN=function(p)p.adjust(p,'BH'))
put(res,'PRE_taxon_community_exact_tests');put(do.call(rbind,prevalues),'PRE_taxon_community_per_woman_values');put(do.call(rbind,units),'PRE_community_selected_library_manifests')
saveRDS(allocations,file.path(o,'results/PRE_exact_label_assignments.rds'))
counts<-unique(res[,c('selection','Location','n_dysbiosis_description','n_Lactobacillus_description','exact_assignments')]);put(counts,'PRE_community_group_denominators')
phase<-unique(m[,c('Participant','VaginalCategory','Phase')]);put(as.data.frame(table(phase$VaginalCategory,phase$Phase)),'Community_category_by_phase_counts')
g<-ggplot(res[res$taxon_arm=='baseline_39'&res$rank=='provisional_species'&res$selection=='retain99EP2'&is.finite(res$p_exact_label_permutation),],aes(median_difference_D_minus_L*100,-log10(p_exact_label_permutation)))+geom_point(aes(color=p_BH_within_site<.05),size=1.8)+facet_wrap(~Location,nrow=1)+scale_color_manual(values=c('FALSE'='grey45','TRUE'='#B2182B'))+theme_bw(base_size=11)+labs(x='Median PRE abundance difference (percentage points)',y='Minus log10 exact nominal p',color='BH q below 0.05',caption='Exploratory PRE-only comparisons of original vaginal-community descriptions. One observation per woman/site.\nTwo participant-99 endometrial choices are reported separately; groups are not clinical diagnoses.')
ggsave(file.path(o,'results/PRE_community_taxon_comparisons_REVIEW.png'),g,width=11,height=4.5,dpi=300)
capture.output(sessionInfo(),file=file.path(o,'R_PRE_sessionInfo.txt'));print(counts);print(aggregate(p_BH_within_site~taxon_arm+rank+selection,res,min))
