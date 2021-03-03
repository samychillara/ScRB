library(ggplot2)
library(gridExtra)
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
    theme(axis.text.x = element_text(angle = 45,hjust=1))
  return(plot)
}
##input: genes X cells
plot.gene.stats.dist=function(data_mat,range=c(0,10)){
  genes=rownames(data_mat)
  counts=rowSums(data_mat)
  data_plot=data.frame(genes=genes,counts=counts)
  g = ggplot()+theme_classic()+theme(axis.title.x=element_text(size=10),
                                     axis.title.y=element_text(size=10),
                                     axis.text.x=element_text(size=10),
                                     axis.text.y=element_text(size=10),
                                     plot.title = element_text(face="bold",size = 12)
  )
  g+geom_histogram(data=data_plot,aes(y=counts),bins = binsize,stat="count")
}


##gene expression distributions


gexp.dist.compare=function(data_list,gene,only_nonzero=T){
  
  g = ggplot2::ggplot()+ggplot2::theme_classic()+ggplot2::theme(plot.margin=unit(c(0,0,0,0),"cm"))+
    
    ggplot2::theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "top")
  
  plot_data=lapply(names(data_list),function(x){
    
    dat=as.data.frame(t(data_list[[x]][gene,]))
    
    if(only_nonzero){
      
      dat=dat[dat[,1]>0,]
      
    }
    
    var=rep(x,length(dat))
    
    return(cbind(counts=as.numeric(dat),var=var))
    
  })
  
  dat=do.call(rbind,plot_data)
  
  dat=data.frame(counts=as.numeric(dat[,1]),var=dat[,2])
  
  g3=g+geom_density(data=dat,aes(counts,fill=var,after_stat(ndensity)),position = "identity",alpha=0.5)+
    
    labs(title=gene)+xlim(c(0,max(dat$counts)))
  
  return(g3)
  
}


##After discreting viewing gene expression in clusters
dotplot.scrb=function(data_dis,metadata=NULL,genes,split.clusterid=NULL,combine.clusterid=NULL,pick.cluster=NULL,title=NULL){

  if(is.list(data_dis)){
    data_dis_list=lapply(data_dis,function(x){
      dat=x$counts[intersect(genes,rownames(x$counts)),]
      return(list(counts=dat,x$metadata))
    })
  }else{
    data_dis=data_dis[intersect(genes,rownames(data_dis)),]
    if(is.null(split.clusterid)){
      data_dis_list=data_dis
    }else{
    data_dis_list=Split.Data(data_dis,metadata,gen.all.combos = F,use.clusterid = split.clusterid,pick.clusters = pick.cluster,return.data = "Q")
    }
  }

  genes_sum=lapply(data_dis_list,function(x){
    if(is.null(dim(x$counts))){
      return(NULL)
    }else {
      gene.activity.summary(x$counts)
    }})

  genes_sum=genes_sum[names(unlist(sapply(genes_sum,function(x){dim(x)[2]})))]
  genes_sum=do.call(rbind,lapply(names(genes_sum),function(x){
    data.frame(genes_sum[[x]],cluster=x,gene=rownames(genes_sum[[x]]))
  }))

  g = ggplot2::ggplot()+ggplot2::theme_classic()+ggplot2::theme(plot.margin=unit(c(0,0,0,0),"cm"))+

    ggplot2::theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "top")

  g1=g+geom_point(data=genes_sum,aes(x=cluster,y=gene,size=pct_active,col=pct_unknown))+ggtitle(title)

  return(g1)
}

##Ratsor plots of discretised data

data_prep<-function(data_list,genes_vec){
  if(is.list(data_list)){
    data_plot=do.call(rbind,lapply(data_list,function(x){
      reshape2::melt(x[,intersect(as.vector(colnames(x)),genes_vec)])
    }))
  }else{

  }

  data_types=do.call(rbind,lapply(names(data_list),function(x){
    cbind(cell.name=rownames(data_list[[x]]),cell.type=rep(x,dim(data_list[[x]])[1]))
  }))
  data_plot=merge(data_plot,data_types,by.x="Var1",by.y="cell.name")
  return(data_plot)
}

plot_gexp_rastor<-function(data_list,genes_vec,by.celltype=F
                           ,multiTissue=F,cell.types=NULL,fixed=F,cols=c(red,blue)){
  if(multiTissue){
    if(is.null(cell.types)==FALSE){
      data_list=subset.cellType(data_list,cell.types)
    }
    data_plot=lapply(data_list,function(x){data_prep(x,genes_vec)})
    tissue=do.call(c,lapply(names(data_plot),function(x){rep(x,dim(data_plot[[x]])[1])}))
    data_plot=do.call(rbind,data_plot)
    data_plot=cbind(data_plot,tissue)

  }else{
    if(is.null(cell.types)==FALSE){
      select=do.call(c,lapply(cell.types,function(y){
        pick=grep(y,names(data_list))
      }))
      data_list=data_list[select]
    }

    data_plot=data_prep(data_list,genes_vec)
  }

  if(by.celltype){
    xVAR=data_plot$cell.type
  }else{
    xVAR=data_plot$Var1
  }
  if(multiTissue){
    plot=ggplot2::ggplot(data_plot,ggplot2::aes(xVAR,Var2,fill=tissue,alpha=value))+ggplot2::geom_raster()+
      ggplot2::theme_minimal()+
      ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme(axis.text.x = ggplot2::element_blank())+
      ggplot2::facet_grid(cols=ggplot2::vars(cell.type),scales = "free")+
      ggplot2::scale_alpha_continuous(range=c(0,1))
  }else{
    plot=ggplot2::ggplot(data_plot,ggplot2::aes(xVAR,Var2,fill=value))+ggplot2::geom_raster()+
      ggplot2::theme_minimal()+
      ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme(axis.text.x = ggplot2::element_blank())+
      ggplot2::facet_grid(cols=ggplot2::vars(cell.type),scales = "free")
  }

  return(plot)
}

##Exploring results of marker analysis
#
# for(j in tissues){
#   pdf(paste0(j,".pdf"),height=15,width=10)
#   for(i in 1:(length(compare[[j]]))){
#     data=compare[[j]][[i]]
#     cell_type=names(compare[[j]])[i]
#     ct_n1=ct_n[[j]]
#
#     genes_ct=unlist(strsplit(as.vector(data$mGenes_exp),","))
#     gene_info=do.call(rbind,lapply(unique(genes_ct),function(x){
#       mat1=data[grep(x,data$mGenes_exp),c("tool","cell.ann")]
#       rank=data[grep(x,data$mGenes_exp),"rank"]
#       mat2=apply(mat1,1,function(y){paste(y,collapse = "-")})
#       return(cbind(gene=rep(x,dim(mat1)[1]),mat1,weight=gene_weight[x],var=mat2,rank))
#     }))
#
#     g=ggplot(gene_info,aes(var,gene,alpha=as.numeric(as.vector(weight))))+
#       geom_tile(aes(fill=cell.ann))+geom_point(aes(shape=as.factor(rank),col=tool),size=4)+
#       scale_x_discrete(labels=NULL)+scale_alpha_continuous(range=c(0.3,1))+
#       scale_color_manual(values=c("blue","darkred","darkgreen","black"))+
#       guides(col=guide_legend("Tool"),
#              fill=guide_legend("Cell Type"),
#              shape=guide_legend("Rank"),
#              alpha=guide_legend("weight"))+
#       labs(title=paste(cell_type,as.vector(ct_n1[cell_type])))+
#       theme_bw()
#
#     #plots1[[i]]=ggplot(data,aes(tool,cell.ann))+geom_point(aes(col=test_score,shape=factor(rank)),size=10)+theme_classic()+
#     #theme(axis.text = element_text(size=16),title=element_text(size=20))+labs(title=paste(cell_type,as.vector(ct_n1[cell_type])))
#
#     # by_ct=split(gene_info,as.factor(gene_info$cell.ann))
#     # by_ct=do.call(rbind,(lapply(names(by_ct),function(x){
#     #   diff=setdiff(markers_all[[x]],unique(by_ct[[x]]$gene))
#     #   mat=cbind(gene=diff,tool="none",cell.ann=x,weight=gene_weight[diff])
#     #   rbind(by_ct[[x]],mat)
#     # })))
#
#     # #  plots=list()
#     #   for(k in unique(gene_info$tool)){
#     #
#     #     plots[[k]]=ggplot(subset(gene_info,tool==k),aes(cell.ann,gene))+geom_tile(aes(fill=weight))
#     #     plots2[[j]][[i]]=grid.arrange(grobs=plots,ncol=3,top=cell_type)
#     #   }
#     #   }
#
#     #   ggsave(paste0(j,".jpg"),p,height=15,width=45,unit="in")
#
#   }
#   #plots2[[j]]=grid.arrange(grobs=plots1,ncol=3)
#
#   dev.off()
# }
#
##input = output of gene.activity.summary
plot.dis.gene.stats=function(data_list){
  data_plot=lapply(names(data_list),function(x){
    var=rep(x,dim(data_list[[x]])[1])
    return(data.frame(data_list[[x]],var=var))
  })
  if(is.list(data_plot)){
    data_plot=do.call(rbind,data_plot)
  }
  g = ggplot2::ggplot()+ggplot2::theme_classic()+ggplot2::theme(plot.margin=unit(c(0,0,0,0),"cm"))+xlab("gene activity in cells(percent)")+ylab("proportion of genes")+
    ggplot2::theme(legend.position = "top")
  g3=g+geom_histogram(data=data_plot,aes(pct_active,fill=var),position = "dodge",stat = "density")+  labs(title="active")
  g4=g+geom_histogram(data=data_plot,aes(pct_inactive,fill=var),position = "dodge",stat = "density")+  labs(title="inactive")
  g5=g+geom_histogram(data=data_plot,aes(pct_unknown,fill=var),position = "dodge",stat = "density")+  labs(title="unknown")
  return(grid.arrange(grobs=list(g3,g4,g5),nrows=1))
}
