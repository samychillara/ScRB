#' Gene expression level computed from raw data
#' input- raw counts data matrix: genes x cells; threshold vector: a vector of integer values which are bins for classifying counts data; percentage of min cells required to express a given gexp value
#' output- named vector: names-gene name, value-maximum gene expression.
#' A gene is highly expressed if there exists at least min.cells percent of cell population that have above th.gexp UMI counts. Otherwise gene is lowly expressed
#' @export
Data.gexp.stats=function(data_counts,th.gexp_vec,min.cells){
 gene.exp=sapply(rownames(data_counts),function(x){
   min.cells=length(data_counts[x,])*(min.cells/100)
   vec=as.numeric(data_counts[x,])
   if(max(vec)==0){
     return(0)
   }else if(max(table(vec))<min.cells){
       return(0)
     }else {
   tab=table(cut(vec,th.gexp_vec,include.lowest = T))
   names(tab)=as.character(th.gexp_vec)[-1]
   tab=tab[is.na(tab)==F]
   classify=max(as.numeric(names(tab[tab>min.cells])))
   return(as.numeric(classify))
   }
 })
  gene.exp=reshape2::melt(gene.exp)
  return(gene.exp)
}

#' Descriptive statistics of gene exp distribution
#' @export
gexp.dstats=function(data_counts,genes){
  data=data_counts[genes,]
  mean=apply(data,1,mean)
  sd=apply(data,1,sd)
  return(list(mean=mean,sd=sd))
}
########################################################################################################
##Utility function
#' To select appropriate cutoff values for number of cells given gexp data
#' input- data: gene expression data; vector of integer values with number of cells; vector of integer values with respective cutoff values (in percentage)
#' output- interger - cutoff value (in percentage)

#' @export
cutoff.select<-function(data,ncells_vec=c(50,250,500,1000),cutoffpct_vec=c(40,30,20,10)){
  ncells=dim(data)[2]
  if(ncells>max(ncells_vec)){
    i=which(ncells_vec==max(ncells_vec))
    cutoff=cutoffpct_vec[i]
  }else if(ncells<min(ncells_vec)){
    i=which(ncells_vec==min(ncells_vec))
    cutoff=cutoffpct_vec[i]
  }else {
  i=which(ncells_vec==min(ncells_vec[ncells_vec>=ncells]))
  cutoff=cutoffpct_vec[i]
  }
  return(cutoff)
}
########################################################################################################
##For discretised data
#' Summary of discretised data: information for every gene, list of active, inactive, unknown genes
#' input- discretised data matrix: genes x cells ; cutoff (default selected by cutoff.select)
#' output- list: 1.geneinfo: matrix genes X gene info; 2.genelist: vectors with names of active, inactive, unknown genes
#' @export
gene.activity.summary=function(data_dis,cutoff=cutoff.select(data_dis),return.sum=c("all","genelist","geneinfo")){
  data_mat=data_dis
  ncells=dim(data_mat)[2]
  data_mat=apply(data_mat,1,function(x){
    ncells_0=length(x[which(x==0)])
    ncells_1=length(x[which(x==1)])
    ncells_0.5=length(x[which(x==0.5)])
    pct_ncells0=round((ncells_0/ncells)*100,2)
    pct_ncells1=round((ncells_1/ncells)*100,2)
    pct_ncells0.5=round((ncells_0.5/ncells)*100,2)
    return(cbind(ncells=ncells,active=ncells_1,inactive=ncells_0,unknown=ncells_0.5,pct_ncells1,pct_ncells0,pct_ncells0.5))
  })
  data_mat=t(data_mat)
  colnames(data_mat)=c("sample_size","active","inactive","unknown","pct_active","pct_inactive","pct_unknown")

  active.genes=data_mat[,"pct_active"][data_mat[,"pct_active"]>=cutoff]
  inactive.genes=data_mat[,"pct_inactive"][data_mat[,"pct_inactive"]>=cutoff]
  unknown.genes=data_mat[,"pct_unknown"][data_mat[,"pct_unknown"]>cutoff]

  genelist=list(active.genes=active.genes,inactive.genes=inactive.genes,unknown.genes=unknown.genes)

  if(return.sum=="all"){
    return(list(geneinfo=data_mat,genelist=genelist))
  }else if(return.sum=="geneinfo"){
    return(geneinfo)
  }else{
    return(genelist)
  }
}

###########################################################################################################
### For mulitple discretised data sets/combined discretised data set with metadata
## compare active, inactive, unknown gene lists
###########################################################################################################
#' Identifying Unique genes sets (active/inactive/unknown)
#' input - list of output of genes.activity.summary for multiple data sets
#' output - list of unique active/inactive/unknown genes for each data set
#' @export
Compare.GeneActivity<-function(data_list,all.comparisons=F){
  # gene.act.sum=lapply(data_list,function(x){
  #   if(is.null(cutoff)){
  #     cutoff=cutoff.select(x)
  #   }else {
  #   cutoff=cutoff
  #   }
  #   return(gene.activity.summary(x,cutoff=cutoff,return.sum = "all"))
  # })
  gene.act.sum=data_list
  genes.compare.active=Genelist.compare(lapply(gene.act.sum,function(x){
    x$genelist$active.genes
  }),all.pairwise = all.comparisons)
  genes.compare.inactive=Genelist.compare(lapply(gene.act.sum,function(x){
    x$genelist$inactive.genes
  }),all.pairwise = all.comparisons)
  genes.compare.unknown=Genelist.compare(lapply(gene.act.sum,function(x){
    x$genelist$unknown.genes
  }),all.pairwise = all.comparisons)

  return(list(active.genes=genes.compare.active,inactive.genes=genes.compare.inactive,unknown.genes=genes.compare.unknown))
}

###to compare given genelists to identify sets of unique and common genes

Genelist.compare<-function(genelist,all.pairwise=T){
  genes_exp_all=Reduce(intersect,lapply(genelist,names))
  genes_Uexp_byeach=lapply(c(1:length(genelist)),function(x){
    set=setdiff(c(1:length(genelist)),x)
    genes_diff=setdiff(names(genelist[[x]]),do.call(c,lapply(genelist[set],names)))
    return(genes_diff)
  })
  names(genes_Uexp_byeach)=names(genelist)
  if(all.pairwise){
    genes_diff_all=lapply(names(genelist),function(x){
      set=setdiff(names(genelist),x)
      sapply(set,function(y){
        setdiff(names(genelist[[x]]),names(genelist[[y]]))
      })
    })
    names(genes_diff_all)=names(genelist)
    return(list(GenesExpByAll=genes_exp_all,
                UniqueGenesExpB=genes_Uexp_byeach,
                PairwiseGeneExp=genes_diff_all))
  }else{
    return(list(GenesExpByAll=genes_exp_all,
                UniqueGenesExpB=genes_Uexp_byeach))
  }

}

############################################################################################################
##Identifying DE genes between data sets
##data_mat = output of gene.activity.summary$geneinfo
############################################################################################################

#' Identifying differentially activated genes (DAG)
#' input- list of gene information (output of gene.activity.summary$geneinfo) for the two data sets that are compared
#' @export
DAG<-function(data_list,cutoff=30,write.files=T,fileid="myQuery"){
  common_genes=Reduce(intersect,lapply(data_list,rownames))
      data_list=lapply(data_list,function(x){x[common_genes,]})

  diff.act_inact= sapply(data_list,function(x){
    diff.pct=x[,"pct_active"]-x[,"pct_inactive"]
  })

  diff.betweenData = apply(diff.act_inact,1,function(x){x[1]-x[2]})

  genes.de=sapply(names(diff.betweenData),function(x){
    if(abs(diff.betweenData[x])>cutoff){
      de=T
      if(diff.betweenData[x]>0){
        stat="active"
      }else{
        stat="inactive"
      }
    }else if(abs(data_list[[1]][x,"pct_unknown"]-data_list[[2]][x,"pct_unknown"])>cutoff){
      de=T
      stat="unknown"
    }else {
      de=F
      stat="not DA"
    }
    return(c(de,stat))
  })
  genes.de=t(genes.de)

  genes.de=cbind(genes.de,diff.betweenData)
  colnames(genes.de)=c("DA","up/down","diff.pct")

  return(genes.de)
}


# ## Analysing gene exp dist
# Genes_exp_cutoffs=function(data_dis,cutoffs,th_cutoffs,pick){
#   th_cutoffs=th_cutoffs[order(cutoffs,decreasing = F)]
#   cutoffs=as.character(cutoffs[order(cutoffs,decreasing = F)])
#   names(th_cutoffs)=cutoffs
#   lapply(data_dis,function(y){
#     th=0.1
#     ncells=dim(y)[1]
#     for(m in cutoffs){
#       if(ncells<=as.numeric(m)){
#         th=th_cutoffs[m]
#         break
#       }
#     }
#     genes_exp=Genes_exp(y,pick=pick,pick.all = F,threshold = th)
#     return(genes_exp)
#   })
#
# }
