#' Data Normalisation using SCTransform
#' input- raw data matrix - gene counts x cells
#' ouput- normalised (default including log transfrom) and corrected data matrix: genes x cells
#' @export
DataNorm<-function(data,log_transform=c(TRUE,FALSE)){
  data=as.matrix(data)
  dat_norm<-sctransform::vst(data,
                             return_corrected_umi = TRUE
                             ,n_genes=NULL
                             ,min_cells=5
  )
  if(log_transform){
    return(log1p(dat_norm$umi_corrected))
  }else{
    return(dat_norm$umi_corrected)
  }
}

#' data normalisation by log transformation
#' input- raw data matrix - gene counts x cells, scaling factor (default=10000)
#' output- log normalised matrix : genes x cells
#' @export
LogNorm<-function(data_mat,s.f=10000){
  norm=apply(data_mat,2,function(x){
    total_counts=sum(x)
    if(total_counts==0){
      return(log1p(x))
    }else{
      return(log1p((x/total_counts)*s.f))
    }
  })
  return(norm)
}



###internal function
sample.data<-function(cells_dat,cellid,nlevels,ncells.sample){
  set_samp=cellid
  for(i in nlevels:1){
    if(i!=nlevels){
      cell.id.split=sapply(set_samp,function(x){paste(strsplit(x,"[.]")[[1]][c(1:i)],collapse = ".")})
      set_cellid=cell.id.split
      names(cell.id.split)=names(set_samp)
      set=split(cell.id.split,as.factor(cell.id.split))
      replace.cells=F
    }else{
      cell.id.split=set_samp
      set=split(set_samp,as.factor(cell.id.split))
      replace.cells=T
    }
    set_cells=unlist(sapply(set,function(x){
      if(i==nlevels){
        pick=intersect(cells_dat,names(x))
      }else{
        pick=names(x)
      }
      if(length(pick)>=1&&length(x)>=ncells.sample){
        sub_cell= sample(pick,ncells.sample,replace=replace.cells)
      }else if(i==nlevels&&length(pick)>=1){
        sub_cell=pick
      }else{
        sub_cell=NA
      }
      return(sub_cell)
    }))
    if(is.null(length(set_cells))){
      set_cells="NA"
      break
    }else{
      set_cells=set_cells[which(is.na(set_cells)==FALSE)]
    }
    set_samp=set_samp[set_cells]
  }
  return(names(set_samp))
}

## internal function
##Filter cells with non zero gene expression
##input: gene expression matrix: genes X UMIs
##ouput: list of genes with names of cells (UMIs) that have non zero gene exp values

cells_exp<-function(data_matrix,return.values=F,min.cells.gexp){
  data_nonzero_list=apply(data_matrix,1,function(x){
    nonzero=x[x!=0]
    if(return.values){
      return(nonzero)
    }else{
      return(names(nonzero))
    }
  })
  len=sapply(data_nonzero_list,length)
  return(data_nonzero_list[len>min.cells.gexp])
}

##Background selection:
#' sampling cells for reference data construction
#'input- cell_clusterid_mat: matrix with columns --"cell.name(UMI)"--same as colnames of counts matrix,"cluster.id","Tissue"; counts_norm: normalised umi counts matrix-- genes X UMIs; ncells.sample: number of cells to sample; clusterids_vec: vector with elements corresponding column names in cell_clusterid_mat with the correct heirarchiel order of clusterids (eg. Tissue, cell type, cell subtype)
#' output- list: for every gene, sampled gene expression vector
#' @export
Background_sampling<-function(counts_norm_mat,cell_clusterid_mat,ncells.sample=50,clusterids_vec=c("Tissue","clusterid.1","clusterid.2"),min.cells.gexp=10){
  #cellidXclusterid
  nlevels=length(clusterids_vec)
  cellid<-apply(cell_clusterid_mat,1,function(x){paste(x[clusterids_vec],collapse = ".")})
  names(cellid)=cell_clusterid_mat$cell.name

  cells_sample=cells_exp(counts_norm_mat,min.cells.gexp=min.cells.gexp)

  genes=names(cells_sample)

  bar=utils::txtProgressBar(min=0,max=length(genes),style = 3)
  ids_samples_byGene=foreach(n=1:length(genes))%do%{
    utils::setTxtProgressBar(bar,n)
      set_samp=sample.data(cells_sample[[n]],cellid,nlevels,ncells.sample)

    return(set_samp)
  }
  close(bar)
  names(ids_samples_byGene)=genes
  genes_exclude=names(ids_samples_byGene)[which(is.null(sapply(ids_samples_byGene,length))==TRUE)]
  pick=setdiff(genes,genes_exclude)
  ids_samples_byGene=ids_samples_byGene[pick]
  new_samples_byGene=lapply(names(ids_samples_byGene),function(x){
    counts_norm_mat[x,ids_samples_byGene[[x]]]
  })
  names(new_samples_byGene)=names(ids_samples_byGene)
  return(new_samples_byGene)
}


## internal function
##no gene expressed
geneNotExp<-function(data){
  data=as.matrix(data)
  counts = rowSums(data)
  return(names(counts[counts==0]))
}

## internal function
##genes with low count
genes_lowCount<-function(data_v,count){
  res=length(data_v)<count
  return(res[res])
}


#' data prep: removing genes with no expression, removing zero gexp values, scaling with maximum gene expression, removing outliers
#' input- normalised gexp matrix(genes x cells)
#' output- list with precessed data, scaling factors, outliers, genes not expressed and genes with low counts
#' @export
Data_prep1<-function(data, filter.genesNotExpressed=T, remove.zeros=T, scale.WithMaxGexp=T, remove.outliers=T){
  if(filter.genesNotExpressed){
    genes_not_expressed<-no_gene_exp(data)
    data_new<-subset(data,rowSums(as.matrix(data))>0)
  }else{
    data_new=data
    genes_not_expressed=NA}
  genes=names(data_new)
  if(remove.zeros){
    data_new<-apply(data_new,1,function(x){keep<-which(x!=0)
    dat=x[keep]
    return(dat)
    })
  }
  if(remove.outliers){
    genes=names(data_new)
    data_new=sapply(data_new,function(x){y=as.numeric(x)
    names(y)=names(x)
    return(y)})
    outliers<-lapply(data_new,function(x){
      boxplot.stats(x)$out
    })
    num_outliers<-sapply(outliers,length)[sapply(outliers,length)!=0]
    num_outliers<-num_outliers[order(num_outliers,decreasing = TRUE)]
    data_new<-lapply(1:length(data_new),function(x){
      data_new[[x]][setdiff(names(data_new[[x]]),names(outliers[[x]]))]
    })
    names(data_new)=genes
  }
  if(scale.WithMaxGexp){
    data_new=sapply(data_new,as.numeric)
    sf<-sapply(data_new,max)
    data_new<-lapply(1:length(sf),function(x){
      data_new[[x]]=data_new[[x]]/sf[x]
    })
  }else{sf=cbind(names(data_new),rep(1,length(data_new)))}

  names(data_new)=genes
  return(list(data_processed=data_new,genes_not_expressed=genes_not_expressed,scaling_factors=sf,outliers=num_outliers))
}


#' data prep: scaling with maximum gene expression, removing outliers
#' input: output of background sampling - list of genes with corresponding gexp
#' output: list of scaled, filtered data, scaling factors, outliers
#' @export
Data_prep2<-function(data,scale.WithMaxGexp=T,remove.outliers=T){
  genes=names(data)
  data_new=lapply(data,function(x){y=as.numeric(x)
  names(y)=names(x)
  return(y)})
  if(remove.outliers){
    outliers<-lapply(data_new,function(x){
      boxplot.stats(x)$out
    })
    num_outliers<-sapply(outliers,length)[sapply(outliers,length)!=0]
    num_outliers<-num_outliers[order(num_outliers,decreasing = TRUE)]
    data_new<-lapply(1:length(data_new),function(x){
      data_new[[x]][setdiff(names(data_new[[x]]),names(outliers[[x]]))]
    })
    names(data_new)=genes
  }else{
    num_outliers=NULL
  }
  if(scale.WithMaxGexp){
    data_new=sapply(data_new,as.numeric)
    sf<-sapply(data_new,max)
    data_new<-lapply(1:length(sf),function(x){
      data_new[[x]]=data_new[[x]]/sf[x]
    })
  }else{sf=cbind(names(data_new),rep(1,length(data_new)))}

  names(data_new)=genes
  return(list(data_processed=data_new,scaling_factors=sf,outliers=num_outliers))
}



#' Scaling gene expression data with given scaling factors
#' input- scaling factors: matrix - rownames=gene symbols, column-value of scaling factor; data matrix: genes X cells
#' output- scaled gene expression matrix: genes X cells
#' @export
scaleMaxGexp<-function(data_mat,scalingfactors,genes){
  bar=utils::txtProgressBar(min=0,max=length(genes),style = 3)
  data_s=foreach(n=1:length(genes))%do%{
    utils::setTxtProgressBar(bar,n)
    data=data_mat[genes[n],]/scalingfactors[genes[n],]
    return(data)
  }
  close(bar)
  names(data_s)=genes
  data_s=do.call(rbind,data_s)
  return(data_s)
}

##internal function
##obtain lower and upper th values
.optLowFun <- function(x,empcdf,maxVal){
  return(abs(x* (1-empcdf(x))))
}

##internal function
.optHighFun <- function(x,empcdf,maxVal){
  return(abs(empcdf(x) * (maxVal-x)))
}

##internal function
minimizeRectangle <- function(empCDF,maxVal,minVal){
  lowThr <- optimize(f = .optLowFun, interval = c(minVal,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum

  highThr<- optimize(f = .optHighFun, interval = c(minVal,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum

  return(c(lowThr,highThr))
}

##internal function
##defining lower gexp value to be lower threshold and higher gexp value to be higher
swap_th<-function(th_dist){
  mat=apply(th_dist,1,function(x){
    if(x[1]>x[2]){
      temp=x[1]
      x[1]=x[2]
      x[2]=temp
      return(c(x[1],x[2]))
    }else{
      return(c(x[1],x[2]))
    }
  }
  )
  return(t(mat))
}


#' Generate threshold distributions per gene
#' input- vector of gene expression values
#' output- matrix: numbootstrap x (lower,upper thresholds)
#' @export
createThresholdDist <- function(gexp, numBootstrapSamples = 1000,swap=TRUE){
  gexp=as.numeric(gexp)
  thrs <- do.call("rbind",lapply(seq(1:numBootstrapSamples),function(x){
    samp <- base::sample(gexp,length(gexp),replace = TRUE)
    if(length(unique(samp))==1){
      th<-c(0,unique(samp))
    }else{
      th<-minimizeRectangle(ecdf(samp),max(samp,na.rm=T),min(samp,na.rm=T))
    }
    return(th)
  }))
  colnames(thrs) <- c("Lower","Upper")
  if(swap){
    thrs=swap_th(thrs)
  }
  return(thrs)
}


#' compute p-val and q vals
#' input- threshold list; query gene expression matrix
#' output- list: q values for gene being active, q values for gene being inactive
#' @export
compute_pval<-function(th_list,gexp_mat,genes,p_correct=T){
  print("computing p values from threshold distributions")
  bar=utils::txtProgressBar(min=0,max=length(genes),style = 3)
  p_vals=foreach(n=1:length(genes))%do%{
    utils::setTxtProgressBar(bar,n)
    x=genes[n]
    p_vals_exp=1-(ecdf(th_list[[x]][,"Upper"])(gexp_mat[x,]))
    p_vals_nexp=ecdf(th_list[[x]][,"Lower"])(gexp_mat[x,])
    return(list(p_vals_exp=p_vals_exp,p_vals_nexp=p_vals_nexp))
  }
  close(bar)
  if(p_correct){
    print("Correcting for multiple testing: computing q values")
    sig_vals=lapply(p_vals,function(x){
      q_vals_exp=stats::p.adjust(x$p_vals_exp,method = "BH")
      q_vals_nexp=stats::p.adjust(x$p_vals_nexp,method = "BH")
      return(list(sig_exp=q_vals_exp,sig_nexp=q_vals_nexp))
    })
  }else{
    sig_vals=list(sig_gexp=p_vals$p_vals_exp,sig_nexp=p_vals$p_vals_nexp)
  }
  names(sig_vals)=genes
  sig_exp=do.call(rbind,lapply(sig_vals,function(x){x$sig_exp}))
  sig_nexp=do.call(rbind,lapply(sig_vals,function(x){x$sig_nexp}))
  colnames(sig_exp)=colnames(gexp_mat)
  colnames(sig_nexp)=colnames(gexp_mat)
  return(list(sig_exp=sig_exp,sig_nexp=sig_nexp))
}

#' Classify gene activity
#' input- output of compute_pval, thresholds
#' output- list: discretised data matrix: genes x cells, significance values
#' @export
classify<-function(sig_vals,th,th_i){
  print("Discretising gene expression values")
  discrete_exp<-matrix(NaN,nrow=dim(sig_vals[[1]])[1],ncol=dim(sig_vals[[1]])[2])
  rownames(discrete_exp)=rownames(sig_vals[[1]])
  colnames(discrete_exp)=colnames(sig_vals[[1]])
  bar=utils::txtProgressBar(min=0,max=dim(sig_vals[[1]])[1],style = 3)
  for(i in 1:dim(discrete_exp)[1]){
    utils::setTxtProgressBar(bar,i)
    for(j in 1:dim(discrete_exp)[2]){
      if(is.nan(sig_vals[["sig_exp"]][i,j])){
        discrete_exp[i,j]=0
      }
      if(is.nan(sig_vals[["sig_nexp"]][i,j])){
        discrete_exp[i,j]=0
      }else if(sig_vals[["sig_exp"]][i,j]<=th){
        discrete_exp[i,j]=1
      }else if(sig_vals[["sig_nexp"]][i,j]<=th){
        discrete_exp[i,j]=0
      }else if(sig_vals[["sig_nexp"]][i,j]>=1-th_i&&sig_vals[["sig_exp"]][i,j]>=th_i){
        discrete_exp[i,j]=0.5
      }
    }
  }
  close(bar)
  return(discrete_exp)
}



