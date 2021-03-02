###Utility functions

#' Split given data by given split parameter, return query or background data as required. Default: returns data that can be queried
#' input: clusterids: matrix with first column containing cell names corresponding to gexp data matrix col names
#' @export
Split.Data<-function(data_mat,clusterids,use.clusterid="cluster.id",pick.clusters=NULL,return.onlyGxp=F,return.data="Q"){
  return.data=match.arg(return.data,c("Q","R","all"),several.ok = F)
  cluster_dat<-split(clusterids,as.factor(clusterids[,use.clusterid]))
  if(is.null(pick.clusters)==FALSE){
    cluster_dat=cluster_dat[pick.clusters]
  }
  data_list<-lapply(cluster_dat,function(x){
    cells_query=as.vector(x[,1])
    cells_bg = setdiff(as.vector(clusterids[,1]),cells_query)
    if(return.onlyGxp){
      query=data_mat[,which(colnames(data_mat)%in%cells_query)]
      bg=data_mat[,which(colnames(data_mat)%in%cells_bg)]
    }else{
    query=list(gexp=data_mat[,which(colnames(data_mat)%in%cells_query)],metadata=x)
    bg=list(gexp=data_mat[,which(colnames(data_mat)%in%cells_bg)],metadata=clusterids[which(clusterids[,1]%in%cells_bg),])
    }
    if(return.data=="all"){
      return(list(query=query,Background=bg))
    }else if(return.data=="Q"){
      return(query)
    }else if(return.data=="R"){
      return(bg)
    }
  })
  return(data_list)
}


#'Combine data lists - each element is a matrix: genesXcells -
#'Output: combines data into one matrix - genesXcells
#' @export
combine.subset<-function(data_list,names.set=names(data_list),return.list=F,intersect=F,subset.genes=NULL){
  data_set=data_list[names.set]
  genes_data_set=lapply(data_set,function(x){
    if(is.list(x)){
      mat=sapply(x,rownames)
    }else{
      mat=rownames(x)
    }
    if(is.list(mat)){
      return(Reduce(intersect,mat))
    }else{
      return(mat)
    }
  })

  genes_common_set=Reduce(intersect,genes_data_set)
  if(is.null(subset.genes)==F){
    genes_common_set=intersect(genes_common_set,subset.genes)
  }
  data_set=lapply(data_set,function(x){
    if(is.list(x)){
      list=lapply(x,function(y){
        y[,genes_common_set]
      })
      if(length(list)>1){
        return(do.call(rbind,list))
      }else{
        return(list[[1]])
      }
    }else{
      list=x[,genes_common_set]
      return(list)
    }
  })
  if(return.list){
    if(intersect){
      return(data_set)
    }else{
      data_set=do.call(c,data_list[names.set])
      return(data_set)
    }
  }else{
    return(do.call(rbind,data_set))
  }
}

#'Read a file such that each column is a list of genes including NAs
#'Output: list with each element corresponding to the column in
#' @export
read.genelist<-function(file,na.rm=T){
  data <- readr::read_csv(file,col_types = cols())
  data_list<-as.list(data)
  if(na.rm){
    data_list=lapply(data_list,function(x){x[which(is.na(x)==F)]})
  }
  return(data_list)
}


#'Read a file such that each column is a list of genes including NAs
#'Output: list with each element corresponding to the column in
#' @export
write.genelist<-function(vec_list,filename,rownames=T){
  table=fill.combine(vec_list)
  write.csv(table,filename,row.names = rownames)
}

fill.combine<-function(list,sort=T){
  if(is.data.frame(list[[1]])){
    max=max(sapply(list,function(x){dim(x)[1]}))
  }else{
    max=max(sapply(list,length))
  }
  list_new=lapply(list,function(x){
    if(is.data.frame(x)){
      diff=max-dim(x)[1]
      if(diff!=0){
        new=apply(x,2,function(y){c(y,rep(" ",diff))})
        return(new)
      }else{
        return(x)
      }
    } else{
      diff=max-length(x)
      if(diff!=0){
        new=c(x,rep(" ",diff))
        return(new)
      }else{
        return(x)
      }
    }
  })
  new_table=do.call(cbind,list_new)
  if(sort){
    new_table=new_table[,sort(colnames(new_table))]
    return(new_table)
  } else{
    return(new_table)
  }
}
