# Reconstruct saved bootstrap scores from exported woman-level values and weights.
o<-'review_response/03_analysis/Taxon_scores_PRE_review_2026-10-28'
b<-readRDS(file.path(o,'results/Score_bootstrap_draws_and_weights.rds'))
w<-read.csv(file.path(o,'results/Woman_site_IgA_IgM_relative_abundances.csv'))
s<-read.csv(file.path(o,'results/Site_specific_matched_IgA_IgM_representation_scores.csv'))
ps<-read.csv(file.path(o,'results/Equal_three_site_secondary_scores.csv'))
eq<-function(a,b)stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-10,check.attributes=FALSE)))
n<-0L
for(arm in c('baseline_39','sensitivity_41'))for(rank in c('provisional_species','genus_aggregation'))for(scope in c('All','Proliferative')){
 mats<-list()
 for(site in c('VAGINA','CERVIX','ENDOMETRIUM')){
  name<-paste(arm,rank,scope,site,sep='_');v<-b[[name]];z<-w[w$taxon_arm==arm&w$rank==rank&w$scope==scope&w$Location==site,]
  A<-M<-matrix(0,length(v$women),length(v$taxa));for(i in seq_along(v$women)){
   zz<-z[z$Participant==v$women[i],];zz<-zz[match(v$taxa,zz$Taxon),];A[i,]<-zz$IgA_relative_mean;M[i,]<-zz$IgM_relative_mean
  }
  sa<-v$draw_weights%*%A;sm<-v$draw_weights%*%M;re<-sa/(sa+sm);re[sa+sm==0]<-NA_real_
  eq(re,v$scores);stopifnot(all(rowSums(v$draw_weights)==length(v$women)))
  zz<-s[s$taxon_arm==arm&s$rank==rank&s$scope==scope&s$Location==site,];zz<-zz[match(v$taxa,zz$Taxon),]
  ci<-apply(re,2,function(x)if(sum(is.finite(x))>1)quantile(x[is.finite(x)],c(.025,.975),names=FALSE) else c(NA,NA))
  eq(ci[1,],zz$CI95_low);eq(ci[2,],zz$CI95_high);eq(colSums(is.finite(re)),zz$bootstrap_valid)
  mats[[site]]<-list(A=A,M=M,women=v$women);n<-n+length(v$taxa)
 }
 v<-b[[paste(arm,rank,scope,'Equal_three_site_secondary',sep='_')]];arr<-array(NA_real_,c(nrow(v$scores),length(v$taxa),3));pt<-matrix(NA_real_,length(v$taxa),3)
 for(j in seq_along(mats)){
  z<-mats[[j]];A<-M<-matrix(0,length(v$women),length(v$taxa));ix<-match(z$women,v$women);A[ix,]<-z$A;M[ix,]<-z$M
  sa<-v$draw_weights%*%A;sm<-v$draw_weights%*%M;arr[,,j]<-sa/(sa+sm);pt[,j]<-colSums(A)/(colSums(A)+colSums(M))
 }
 re<-apply(arr,c(1,2),function(x)if(all(is.finite(x)))mean(x) else NA_real_);eq(re,v$scores)
 zz<-ps[ps$taxon_arm==arm&ps$rank==rank&ps$scope==scope,];zz<-zz[match(v$taxa,zz$Taxon),]
 ci<-apply(re,2,function(x)if(sum(is.finite(x))>1)quantile(x[is.finite(x)],c(.025,.975),names=FALSE) else c(NA,NA))
 eq(ci[1,],zz$CI95_low);eq(ci[2,],zz$CI95_high);eq(colSums(is.finite(re)),zz$bootstrap_valid)
 eq(apply(pt,1,function(x)if(all(is.finite(x)))mean(x) else NA_real_),zz$equal_three_site_secondary_score)
}
writeLines(paste('PASS: reconstructed all saved score bootstrap draws and percentile intervals;',n,'site score rows and all secondary joint-site rows.'),file.path(o,'results/Bootstrap_verification.txt'))
cat('Bootstrap score and interval verification PASS\n')
