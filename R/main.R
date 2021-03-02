##Obtain threshold distributions from background data
##GiveN: ScRNAseq data: input: UMI counts - genesXsamples, cellclusterids- matrix: samples(cell.name)Xids(cluster.id)XTissue
#library(foreach)
#library(parallel)
#library(progressr)
#library(doFuture)

#' Data Normalisation by ScTransform
#' @export
scrb_norm<-function(data_counts,log_transform = T,data=c("Ref","query"),file.id="mydata",saveData=T,return=T){
  tictoc::tic("SCT normalisation")
  print("Normalising data")
  data_norm<-DataNorm(data_counts,log_transform = log_transform)
  tictoc::toc()
  print("saving data")
  if(saveData){
    if(data=="Ref"){
    dir_name=paste0("ScRbRef_",file.id)
    }else if (data=="query"){
      dir_name=paste0("ScRbQuery_",file.id)
    }
    dir.create(dir_name)
    saveRDS(data_norm,paste0(dir_name,"/",file.id,"_NormData.rds"))
  }
  if(return){
    return(as.matrix(data_norm))
  }
}


#' Sampling, scaling and removing outliers
#' @export
PrepRefData<-function(data_norm,metadata,clusterid_vec=c("Tissue","clusterid.1","clusterid.2"),file.id="mydata",saveData=T,return=T,
                      sample_data=T,sample_n=50,data_prep=2){
  print("preparing data to generate gene threshold distributions")
  tictoc::tic("prepared data to generate gene threshold distributions")
  data_bg=list()

  data_bg[["Norm_data"]]=data_norm
  data_bg[["Metadata"]]=metadata

  if(sample_data){
    print("Sampling cells for every gene to construct reference data per gene")
    tictoc::tic("Background Sampling")
    data_bg[["Sampled_data"]]=Background_sampling(data_bg[["Norm_data"]],data_bg[["Metadata"]],ncells.sample = sample_n,clusterids_vec =clusterid_vec)
    tictoc::toc()
    if(saveData){
      print("saving data")
      dir_name=paste0("ScRbRef_",file.id)
      if(dir.exists(dir_name)==FALSE){
      dir.create(dir_name)
      }
      saveRDS(data_bg[["Sampled_data"]],paste0(dir_name,"/",file.id,"_PrepData.rds"))
    }
  }
  if(is.null(data_prep)==F){
    print("Scaling sampled data")
    if(data_prep==1){
      tictoc::tic("removing zeros from the data, removing genes not expressed, removing outliers, scaling gexp")
      data_bg[["filtered"]]<-Data_prep1(data_bg$Sampled_data)
      tictoc::toc()
    }else if(data_prep==2){
      print("removing outliers, scaling gexp")
      tictoc::tic("Scaled data")
      data_bg[["filtered"]]<-Data_prep2(data_bg$Sampled_data)
      tictoc::toc()
    }
    if(saveData){
      print("saving data")
      dir_name=paste0("ScRbRef_",file.id)
      if(dir.exists(dir_name)==FALSE){
        dir.create(dir_name)
      }
      saveRDS(data_bg,paste0(dir_name,"/",file.id,"_PrepData.rds"))
    }
    print(paste("saving scaling factors to file"))
    write.table(data_bg$filtered$scaling_factors,file=paste0(dir_name,"/",file.id,"_sf.txt"),col.names = F)
  }
  tictoc::toc()
  if(return){
    return(data_bg$filtered$data_processed)
  }
}

#' Computing threshold distributions
#' @export
GetThreshDist<-function(data,parallelize=F,ncores=detectCores()-2,file.id,return=T,saveData=T){
  check=names(which(sapply(data,length)<1))
  if(length(check)>0){
    data=data[setdiff(names(data),check)]
  }
  th_list<-list()

  print(paste("computing threshold distributions for",length(data),"genes"))
  tictoc::tic("computedthresholds")
  if(parallelize){
    cl = parallel::makeCluster(ncores)
    registerDoParallel(cl)
    handlers("progress")
    registerDoFuture()
    plan(multisession)

    nmodules=c(1:ncores)*as.integer(length(data)/ncores)
    nmodules[ncores]=nmodules[ncores]+(length(data)%%ncores)
    start=c(1,nmodules[1:ncores-1]+1)
    progressr::with_progress({
      bar=progressr::progressor(along=1:length(data))
      th_list =foreach::foreach(j=1:ncores)%:%
        options(future.rng.onMisuse="ignore")
      foreach::foreach(n=start[j]:nmodules[j])%dopar%{
        len=start[j]:nmodules[j]
        bar(sprintf("n=%g",n))
        Sys.sleep(0.1)
        createThresholdDist(data[[n]],numBootstrapSamples = 1000)
      }
    })
    th_list=do.call(c,th_list)
    names(th_list)=names(data)
    tictoc::toc()
    parallel::stopCluster(cl)
  }else{
    bar=utils::txtProgressBar(min=0,max=length(data),style = 3)
    th_list = foreach::foreach(n=1:length(data))%do%{
      utils::setTxtProgressBar(bar,n)
      createThresholdDist(data[[n]],numBootstrapSamples = 1000)
    }
    close(bar)
    names(th_list)=names(data)
    tictoc::toc()
  }
  if(saveData){
    print("saving data")
    dir_name=paste0("ScRbRef_",file.id)
    if(dir.exists(dir_name)==FALSE){
      dir.create(dir_name)
    }
    saveRDS(th_list,paste0(dir_name,"/",file.id,"_th_dist.rds"))
  }
  if(return){
    return(th_list)
  }
}


####################################################################################################################333
############ QUERY ############

#' Normalising, Scaling query data
#' @export
PrepQueryData<-function(query_dat,Ref_data,file.id="query",normalise=T,scale_data=T,saveData=T,return.data=T){
  scaling.factors.file = list.files(Ref_data,pattern = "sf",full.names = T)
  dir_name=paste0("ScRbQuery_",file.id)
  if(dir.exists(dir_name)==FALSE){
    dir.create(dir_name)
  }
  log.file=paste0(dir_name,"/",file.id,"_log.txt")
  cat("Log file",file=log.file,sep="\n",append=T)
  cat(paste("Find data in",paste0(getwd(),"/",dir_name)),file=log.file,sep="\n",append = T)
  cat("Preparing data: Normalisation and Scaling",file=log.file,sep="\n",append = T)
  tictoc::tic.clearlog()
  ##normalisation
  if(normalise){
    print("Normalising data: ScTransform")
    tictoc::tic("Data Normalisation")
    query_Normdata=scrb_norm(query_dat,file.id = file.id,data="query",saveData = T,return = T)
    tictoc::toc(log=T)
  }else{
    query_Normdata=query_dat
  }
  ##scaling
  sf=read.delim(scaling.factors.file,sep=" ",row.names = 1,header = F)
  if(scale_data){
    genes_common<-intersect(rownames(query_Normdata),rownames(sf))
    print(paste("Number of genes in the query with corresponding reference distributions: ",length(genes_common)))
    print("Scaling query data with scaling factors obtained from the background")
    tictoc::tic("Scaling data")
    query_scaled_data=scaleMaxGexp(query_Normdata,sf,genes_common)
    tictoc::toc(log=T)
    if(saveData){
      print("saving scaled data")
      saveRDS(query_scaled_data,file=paste0(dir_name,"/",file.id,"_ScaledData.rds"))
    }
  }
  log=tictoc::tic.log(format = T)
  tictoc::tic.clearlog()
  cat("Time log",file=log.file,sep="\n",append = T)
  cat(unlist(log),file=log.file,sep="\n",append = T)
  cat("End Time log",file=log.file,sep="\n",append = T)
  if(return.data){
  return(query_scaled_data)
  }
}

#' computing q values and discretising query data
#' @export
Discretise_data<-function(query_ScaledData,Ref_data,th=0.5,th_i=0.5,file.id="myquery",return.data=T,saveData=T){
  dir_name=paste0("ScRbQuery_",file.id)
  if(dir.exists(dir_name)==FALSE){
    dir.create(dir_name)
  }

  log.file=paste0(dir_name,"/",file.id,"_log.txt")

  cat("Discretising data",file=log.file,sep="\n",append = T)
  tictoc::tic.clearlog()

  threshold.dist.file=list.files(Ref_data,pattern = "_dist",full.names = T)
  threshold.dist=readRDS(threshold.dist.file)
  ##computing pvals
  genes_common<-intersect(rownames(query_ScaledData),names(threshold.dist))
  tictoc::tic("computed q values")
  query_pvals=compute_pval(threshold.dist,query_ScaledData,genes_common)
  tictoc::toc(log=T)
  ##discretising gene expression
  tictoc::tic("Discretised data")
  query_discretised = classify(query_pvals,th,th_i)
  tictoc::toc(log=T)
  data=list(discretised_data=query_discretised,significance=query_pvals)
  if(saveData){
    print("saving discretised data")
    saveRDS(data,file=paste0(dir_name,"/",file.id,"_DiscretisedData.rds"))
  }
  log=tictoc::tic.log(format = T)
  tictoc::tic.clearlog()
  cat("Time log",file=log.file,sep="\n",append = T)
  cat(unlist(log),file=log.file,sep="\n",append = T)
  cat("End Time log",file=log.file,sep="\n",append = T)
  if(return.data){
  return(data)
  }
}



