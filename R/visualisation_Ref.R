##Potential ggplot templates
g = ggplot()+theme_classic()+theme(axis.title.x=element_text(size=10),
                                   axis.title.y=element_text(size=10),
                                   axis.text.x=element_text(size=10),
                                   axis.text.y=element_text(size=10),
                                   plot.title = element_text(face="bold",size = 12)
)
g1 = ggplot()+theme_classic()+theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none")


##Data exploration
plot.stats.metadata = function(metadata,data.clusterid="cluster.id",split.by="Tissue"){
  data.plot=metadata[,c(data.clusterid,split.by)]
  colnames(data.plot)=c("cluster.id","split.by")
  data.plot=as.data.frame(data.plot)
  plot = ggplot(data.plot,aes(x=cluster.id))+
    geom_bar(stat="count",aes(fill=cluster.id))+
    #facet_grid(cols=vars(Tissue),scales = "free")+
    facet_wrap(vars(split.by),scales="free")+
    geom_text(stat="count",aes(label=stat(count)),vjust=-1)+
    theme_bw()+labs(y="Number of cells")+
    theme(legend.position = "")+
    theme(axis.text.x = element_text(angle = 90,hjust=1))
  return(plot)
}


##Plots to visualise various steps of ScRB pipeline 
## Normalisation - sctransform 
gexp.dist.compare=function(data_list,gene,only_nonzero=T){ 
  g = ggplot2::ggplot()+ggplot2::theme_classic()+ggplot2::theme(plot.margin=unit(c(0,0,0,0),"cm"))+ 
    ggplot2::theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "top") 
  
    plot_data=lapply(names(data_list),function(x){ 
    dat=as.data.frame(t(data_list[[x]][gene,])) 
    if(only_nonzero){ 
      dat=dat[dat[,1]>0,] 
    } 
    var=rep(x,length(dat)) 
    return(cbind(counts=as.numeric(dat[,1]),var=var)) 
  }) 
  dat=do.call(rbind,plot_data) 
  dat=data.frame(counts=as.numeric(dat[,1]),var=dat[,2]) 
  
  g3=g+geom_histogram(data=dat,aes(counts,fill=var),position = "dodge",stat = "density")+ 
  labs(title=gene)+xlim(c(0,max(dat$counts))) 
  return(g3) 
  
} 


##Sampling
show.sampling<-function(gene,data_unsampled,data_sampled,metadata,metadata.cellid="cell.name"){
  clusterids=setdiff(colnames(metadata),metadata.cellid)
  rownames(metadata)=metadata[,metadata.cellid]
  data_
  data_plot_s=as.data.frame(cbind(cell=names(data_sampled[[gene]]),gexp=data_sampled[[gene]],metadata[names(data_sampled[[gene]]),]))
  data_plot_us=as.data.frame(cbind(cell=names(data_unsampled[[gene]]),gexp=data_unsampled[[gene]],metadata[names(data_unsampled[[gene]]),]))
  p1=ggplot(data_plot_us,aes(gexp))+geom_bar(stat="count",width=0.1)+xlim(c(range(data_unsampled[[gene]])))
  p2=ggplot(data_plot_s,aes(gexp))+geom_bar(stat="count",width = 0.1)+xlim(c(range(data_unsampled[[gene]])))
  grid.arrange(p1,p2)
}

##Outliers and Scaling


##threshold distributions ---- written for old sampling !!
g_ecdf = function(data,metadata,th_m,gene){
  ggplot=ggplot()+theme_classic()+theme(axis.title.x=element_text(size=10),
                                            axis.title.y=element_text(size=10),
                                            axis.text.x=element_text(size=10),
                                            axis.text.y=element_text(size=10),
                                            plot.title = element_text(face="bold",size = 12)
  )
  th_l = mean(th_m[[gene]][,1])
  th_h = mean(th_m[[gene]][,2])
  dat = as.data.frame(data[[gene]])
  #samp_set = rep(sets,length(dat[,1]))
  samp_set=metadata$Tissue[match(rownames(dat),metadata$cell.name)]
  dat=cbind(dat,samp_set)
  lim=(range(dat[,1])[2]-range(dat[,1])[1])/4
  g2 = ggplot+geom_histogram(data=dat,aes(dat[,1],col=samp_set))+
    xlab("")+
    ylab("")+
    theme(plot.margin=unit(c(0,0,0,0),"cm"))+
    theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position = "none")
  g3 = ggplot+geom_histogram(aes(th_m[[gene]][,1]),fill="blue")+
    geom_histogram(aes(th_m[[gene]][,2]),fill="red")+
    xlab("")+
    ylab("")+
    theme(plot.margin=unit(c(0,0,0,0),"cm"))+
    theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank())
  g1=ggplot+stat_ecdf(data=dat,aes(dat[,1]))+
    ylim(0,2.5)+
    annotation_custom(grob=ggplotGrob(g3)
                      ,ymin=1,ymax=1.8
                      ,xmin=-Inf,xmax=Inf)+
    annotation_custom(grob=ggplotGrob(g2)
                      ,ymin=1.8,ymax=2.5
                      ,xmin=-Inf,xmax=Inf)+
    labs(title=paste(gene,"Background"))+
    xlab("bg_norm_data")+
    ylab("")+
    geom_vline(xintercept=th_l,col="blue")+
    geom_vline(xintercept = th_h,col="red")
  return(g1)
}
