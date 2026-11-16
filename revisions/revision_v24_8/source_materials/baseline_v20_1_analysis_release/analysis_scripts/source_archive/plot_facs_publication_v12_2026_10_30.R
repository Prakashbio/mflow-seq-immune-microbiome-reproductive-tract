args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('review_response/facs_integrated_working_v12_2026-10-30')
library(ggplot2)
d <- read.csv(file.path(root,'FACS_reconstruction/data/Rebuilt_primary_42_woman_site_records.csv'))
s <- read.csv(file.path(root,'FACS_reconstruction/data/Rebuilt_Table_2_and_pointwise_CIs.csv'))
ep <- c('IgA_only','IgM_only','Double_positive','IgG_positive')
labs <- c('IgA','IgM','IgA/IgM\ndouble+','IgG')
long <- do.call(rbind,lapply(seq_along(ep),function(i) data.frame(Woman=d$Participant,Location=d$Location,endpoint=ep[i],percent=d[[paste0(ep[i],'_percent')]])))
sites <- c('VAGINA','CERVIX','ENDOMETRIUM')
sl <- c('Vagina\nn = 16 women','Cervix\nn = 10 women','Endometrium\nn = 16 women')
long$Location <- factor(long$Location,levels=sites,labels=sl);s$Location <- factor(s$Location,levels=sites,labels=sl)
long$endpoint <- factor(long$endpoint,levels=ep,labels=labs);s$endpoint <- factor(s$endpoint,levels=ep,labels=labs)
p <- ggplot(long,aes(endpoint,percent))+
 geom_point(aes(colour=endpoint),position=position_jitter(width=.13,height=0,seed=20261030),size=1.3,alpha=.7)+
 geom_errorbar(data=s,aes(ymin=ci95_low_percent,ymax=ci95_high_percent,y=mean_percent),width=.16,linewidth=.45,colour='black')+
 geom_point(data=s,aes(y=mean_percent),shape=95,size=5,colour='black')+
 facet_wrap(~Location,nrow=1)+scale_colour_manual(values=c('#D78037','#7658A2','#666666','#3A8A64'),guide='none')+
 scale_y_continuous(limits=c(0,50),breaks=seq(0,50,10),expand=expansion(mult=c(.015,.02)))+
 labs(x=NULL,y='Ig-gated events (% of\nDAPI-negative parent)')+
 theme_classic(base_size=9.5,base_family='sans')+
 theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=10),axis.text.x=element_text(size=8.5),axis.title.y=element_text(size=9),panel.spacing=grid::unit(.8,'lines'))
out <- file.path(root,'figure');dir.create(out,showWarnings=FALSE)
ggsave(file.path(out,'Main_Figure_3_v12.png'),p,width=6.5,height=3.1,dpi=600,bg='white')
ggsave(file.path(out,'Main_Figure_3_v12.svg'),p,width=6.5,height=3.1,device=svglite::svglite,bg='white')
capture.output(sessionInfo(),file=file.path(out,'FACS_plot_session_info_v12.txt'))
