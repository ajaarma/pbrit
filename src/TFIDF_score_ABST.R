######################################################################
# TFIDF calculation for all annotation sources. 
# It employs usage of "SNOW" package for parallelizing the computation
#
# Author: Ajay Anand Kumar
#         ajay.kumar@uantwerpen.be
#
#                       
#
######################################################################
#args = (commandArgs(TRUE))
library("Matrix")
#load("go_sparse_mat")
#source("process_sparse.R")
library("snow")
library(hash)
#load("col_hash_obj")

########################################
#
#  LOADING SPARSE MATRIX from TEXT FILE
########################################
process_sparse <- function(sparse_data_arg,sparse_row_arg,sparse_col_arg)
{
  sparse_matr = read.table(sparse_data_arg,sep="\t",header=FALSE,quote="")
  sparse_row = read.table(sparse_row_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")
  sparse_col = read.table(sparse_col_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")

  sz_rownames = dim(sparse_row)[1]
  sz_colnames = dim(sparse_col)[1]

  if(max(sparse_matr[,1]) > sz_rownames[1])
  {
      print("row names don't fit the matrix. Too few row names")
  }
  
  if(max(sparse_matr[,2]) > sz_colnames[1])
  {
     print("col names don't fit the matrix. Too few col names")
  }

  data <- sparseMatrix(i = sparse_matr[,1], j = sparse_matr[,2],
                       dims = c(sz_rownames[1], sz_colnames[1]))

  rownames(data) = sparse_row[,1]
  colnames(data) = sparse_col[,1]
  out <- list(sparse_data = data, all_colnames = sparse_col, all_rownames = sparse_row)
  return(out)
}

###################################################
#
# Computing TFIDF using parallel SNOW package
#
#
###################################################
myFunc <- function(query_list)
{
  library(hash)
  library("Matrix")

  counter = query_list[[5]]
  if(counter !=20)
  { 
     cat("First few: ",query_list[[1]][1:4],"\n",file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/console.txt",append=TRUE)
     row.ind = query_list[[1]]
  }
  else
  {
    row.ind = query_list[[1]]
  }
  mat_new = query_list[[2]]
  col_list = query_list[[3]]
  row_list = query_list[[4]]
  counter = query_list[[5]]
  #cat("after mat new\n",file="console.txt",append=TRUE)
  #row.name = query_list[[2]][query_list[[1]]]
  #cat(length(query_list[[1]]),"\n",file="console.txt",append=TRUE)
  #mat = sparseMatrix(row.name,query_list[[3]],dims = c(length(query_list[[1]]),length(query_list[[3]])))
  #cat("After sparse mat",file="console.txt",append=TRUE)
  #mat = query_list[[2]]
  nrow = dim(query_list[[2]])[1]
  col_list = query_list[[3]]
  for(ele in 1:length(row.ind))
  {
    cat("COUNTER:",ele,"\t",counter,"\n",file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/console.txt",append=TRUE)
    #ind = which(mat_new[ele,] !=0)
    row_col_ind = row.ind[ele]
    
    ind = row_list[[row_col_ind]]
    #termFreq = 1/length(ind)    
    for(i in ind)
    {
       
       termFreq = 1+log(mat_new[ele,i])
       a1 = mat_new[ele,i]
       #termFreq = mat_new[ele,i]/length(ind)
       #pst_str = noquote(paste("\'",i,"\'",sep=""))
       #idf = log(nrow/(as.numeric(values(col_hash,keys=i))+1))
       #cat("COL hash values: ",counter,"\t",ele,"\t",col_list[[i]],"\n",file="/home/aakumar/WORK_LOCAL/MAPPINGS/ENSG_PUBMED/console.txt",append=TRUE)
       idf = log(nrow/as.numeric(col_list[[i]]+1))
       #cat("TFIDF value:",tfidf,"\n",file="console.txt",append=TRUE)
       mat_new[ele,i]=idf*termFreq
       tfidf = idf*termFreq
       #cat("TFIDF value:",tfidf,"\t",as.numeric(col_list[[i]]+1),a1,"\n",file="/home/aakumar/WORK_LOCAL/PUBMED_TEXT/RUN_2/console1.txt",append=TRUE)
       #cat("MAT:",ele,"\t",i,"\t",counter,"\n",file="/home/aakumar/WORK_LOCAL/PUBMED_TEXT/RUN_2/console2.txt",append=TRUE)
       
    }
    normal = sqrt(sum(mat_new[ele,]^2))
    if(normal !=0)
    {
      #cat("Normal value:",normal,"\n",file="console.txt",append=TRUE)
      mat_new[ele,] = mat_new[ele,]/normal
    }
  }
  return(mat_new)
}

################################################################################
#
# parColSel : Function to do parallel computation for the the columns of sparse 
#             matrix.
#
# #####################################################

parColSel = function(query_ind_arg)
{
  library("hash")
  library("Matrix")
  seq.num = query_ind_arg[[1]]
  mat = query_ind_arg[[2]]#[,query_ind_arg[[1]]]
  ncol = dim(mat)[2]
  counter = query_ind_arg[[3]]
  col_list = list()
  for(i in 1:ncol)
  {
     cat("COL: ",counter,"\t",i,"\t",seq.num[i],"\n",file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/console.txt",append=TRUE)
     ind = as.numeric(which(mat[,i]!=0))
     col_list[[seq.num[i]]]=length(ind)    
  }  
  return(col_list)
}

parRowSel = function(query_ind_arg)
{
  library("hash")
  library("Matrix")
  seq.num = query_ind_arg[[1]]
  mat = query_ind_arg[[2]]#[,qyery_ind_arg[[1]]]
  counter = query_ind_arg[[3]]
  nrow = dim(mat)[1]
  row_list = list()
  for(i in 1:nrow)
  {
   cat("ROW: ",counter,"\t",i,"\t",seq.num[i],"\n",file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/console.txt",append=TRUE)
   ind = as.numeric(which(mat[i,]!=0))
   row_list[[seq.num[i]]] = ind   
  }
  return(row_list)
}
######################################################################
#
# Initialisin matrix for preparing the sample 
# 
#  colMat_init: Divinding the column matrix for columns of the matrix
#
#
######################################################################


colMat_init = function(mat_obj,num_node)
{
  print("inside col Mat")
  mat = mat_obj
  ncol = dim(mat)[[2]] 
  seq.num = seq(1,ncol,round(ncol/num_node))
  print(seq.num)
  query.ind = list()
  for(i in 1:length(seq.num))
  {
    print(i)
    if(i != length(seq.num))
    {
     col_ind = c(seq.num[i]:(seq.num[i+1]-1))
     query.ind[[i]] = list(col_ind,mat[,col_ind],i)
    }
    else
    {
      col_ind = c(seq.num[i]:ncol)
      query.ind[[i]] = list(col_ind,mat[,col_ind],i)
    }
  }
  save(query.ind,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/query_ind_col") 
  tryCatch({
  #print("before cl")
   cl = makeCluster(length(query.ind),type="SOCK")
   result = clusterApply(cl,query.ind,parColSel)
   },interrupt = function(ex){print(ex)},
   error = function(ex){print(ex)})
  stopCluster(cl)
  save(result,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/result_col")
  mat.final = list()
  mat.final = result[[1]]
  #save(result,file="MATRICES/col_hash_result_mat")
  for(i in 2:length(query.ind))
  {
   print(i)
   A = result[[i]]
   mat.final = c(mat.final,A[!sapply(A,is.null)])#result[[i]])
  }  
  col_list = mat.final
  save(col_list,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/col_list_result")
  return(mat.final)
}

RowMat_init = function(mat_obj,num_node)
{
  print("Inside Row mat")
  mat = mat_obj
  nrow = dim(mat)[1]
  print(nrow)
  seq.num = seq(1,nrow,round(nrow/num_node))
  print(seq.num)
  query.ind = list()
  for(i in 1:length(seq.num))
  {
    print(i)
    if(i != length(seq.num))
    {
      row_ind = c(seq.num[i]:(seq.num[i+1]-1))
      query.ind[[i]] = list(row_ind,mat[row_ind,],i)
    }
    else
    {
      row_ind = c(seq.num[i]:nrow)
      query.ind[[i]] = list(row_ind,mat[row_ind,],i)
    }
  }
  print("Begining Cluster in Row Mat")   
  tryCatch({
   cl = makeCluster(length(query.ind),type="SOCK")
   result = clusterApply(cl,query.ind,parRowSel)
   },interrupt = function(ex){print(ex)},
   error = function(ex){print("Inside Row Mat",ex)})
  stopCluster(cl)
  save(result,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/result_row")
  mat.final = list()
  mat.final = result[[1]]
  #save(result,file="MATRICES/col_hash_result_mat")
  for(i in 2:length(query.ind))
  {
   print(i)
   A = result[[i]]
   mat.final = c(mat.final,A[!sapply(A,is.null)])#result[[i]])
  }
  row_list = mat.final  
  save(row_list,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/row_list_result")
  return(mat.final)
}


#########################################################
#
# Initializing The Matrix for SNOW package run
#
#
#########################################################


matrix_init <- function(mat_obj, num_node, mode)#, rowInd=rowInd)
{
  mat = mat_obj[[1]]
  if(mode ==1)
  {
     mat <- mat_obj[[1]]
     colnames(mat) = gsub(" ","",mat_obj[[2]][,1],fixed=TRUE)
     rownames(mat) = gsub(" ","",rownames(mat),fixed=TRUE)
     #mat.old = mat
     mat = mat*1
     nrow = dim(mat)[1]
     ncol = dim(mat)[2]
  }
  else
  {
     mat = sparseMatrix(length(colnames(mat)),length(colnames(mat)))
     rownames(mat) = colnames(mat)
     colnames(mat) = colnames(mat)
     mat = mat*1
     nrow = dim(mat)[1]
     ncol = dim(mat)[2]

  }
    col_list = list()
    col_list = colMat_init(mat,num_node)
    #load("/home/aakumar/WORK_LOCAL/PUBMED_TEXT/RUN_2/MATRICES/col_list_result")
    #col_hash = mat.final
    #col_list=mat.final
    #mat.final=c()
    row_list = list()
    row_list = RowMat_init(mat,num_node)
    #col_list = mat.final
    #mat.final = c()
    #load("/home/aakumar/WORK_LOCAL/PUBMED_TEXT/RUN_2/MATRICES/row_list_result")
    #row_list = mat.final
    #mat.final = c()
    
  #for(ele in 1:ncol)
  #{
    #cat(ele,"\n",file="/home/aakumar/WORK_LOCAL/MAPPINGS/ENSG_PUBMED/console.txt",append=TRUE)
    #if (mode ==1)
    #{
      #ind = as.numeric(which(mat[,ele]!=0))
      #col_hash[ele]=length(ind)
    #}
    #else
    #{
      #ind = which(mat[,ele]!=0)
      #col_hash[ele] = as.numeric(ind)
    #}
  #}
  query.ind = list()
  #seq.num = seq(1,nrow,round(nrow/30))
  seq.num = seq(1,nrow,round(nrow/num_node))
  print(nrow)
  for(i in 1:length(seq.num))
  {
    print(i)
    if(i != length(seq.num))
    {
       row_ind = c(seq.num[i]:(seq.num[i+1]-1))
       query.ind[[i]] = list(row_ind,mat[row_ind,],col_list,row_list,i)
    }
    else
    {
      row_ind = c(seq.num[i]:nrow)
      query.ind[[i]] = list(row_ind,mat[row_ind,],col_list,row_list,i)
    }
  }
  #save(query.ind,file="/home/aakumar/WORK_LOCAL/PUBMED_TEXT/RUN_5/MATRICES/query_ind_mat")
  return(query.ind)
}


TFIDFWeight <- function(mat_obj)
{
    mat1 = mat_obj[[1]]
    mat = as.matrix(mat1)
    ND = dim(mat)[1]
    NW = dim(mat)[2]
    for(i in 1:ND){
             cat("ND\t",i,"\n",file="console.txt",append=TRUE )
             ind = which(mat[i,] !=0)
             for(e in ind){
                    mat[i,e] = 1+log(mat[i,e])
             }
    }
   
    for(i in 1:NW){
             cat("NW",i,"\n",file="console.txt",append=TRUE)
             mat[,i] = mat[,i]*log2(ND/sum(mat[,i]>0,na.rm=T)+1)
    }

    for(i in 1:ND){
             cat("Normal ND",i,"\n",file="console.txt",append=TRUE)
             mat[i,which(is.na(mat[i,])==TRUE)] = 0
             mat[i,] = mat[i,]/sqrt(sum(mat[i,]^2,na.rm=T))
    }
   
    mat[is.nan(mat)] = 0
    mat
}


#############################################
#
# Making the SNOW Cluster 
#
############################################


make_cluster <- function(query.ind,mode)
{
  cl = makeCluster(length(query.ind),type="SOCK")
  #cl = makeCluster(15,type="SOCK")
  cat("after clusater")
  tryCatch({
  #print("before cl")
  if(mode==1)
  {
   result = clusterApply(cl,query.ind,myFunc)
  }
  else
  {
   print("Inside else, make cluster")
   result = clusterApply(cl,query.ind,myInfoGain)
  }
   },interrupt = function(ex){print(ex)},
   error = function(ex){print(ex)})  
  stopCluster(cl)
  mat.final = result[[1]]
  save(result,file="/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES/CL_result_mat")
  for(i in 2:length(query.ind))
  {
   mat.final = rBind(mat.final,result[[i]])

  }
  return(mat.final)
}

hist_plot = function(obs_mat,ind)
{
   ind1 = which(obs_mat[ind,]!=0)
   
   y = as.vector(obs_mat[ind,ind1])
 
   x = c(1:length(y))
   #plot(x,y,xlab=colnames(obs_mat),ylab="TFIDF values")
   ##plot(x,y,xlab=colnames(mat.PB.final)[which(mat.PB.final[9665,]!=0)],ylab="TFIDF Values",col="red")
   plot(x,y,xlab=colnames(obs_mat)[ind1],ylab="TFIDF Values",col="red")
   text(x,y,colnames(obs_mat)[ind1],srt=90,cex=0.8,pos=1,xpd=TRUE)
   ##text(x,y,colnames(mat.PB.final)[which(mat.PB.final[9665,]!=0)],srt=90,cex=0.8,pos=1,xpd=TRUE)

}

mat_sparse = "/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/INPUT/HGNC_sparse.txt"#args[[1]]
mat_row = "/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/INPUT/HGNC_rownames.txt"#args[[2]]
mat_col = "/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/INPUT/HGNC_colnames.txt"#args[[3]]
mat_dir = "/home/shared_data_medgen_aorta/pbrit_data/V3/PUB/MATRICES"#args[[4]]
anno_type = "PUB"#args[[5]]
cat("Loading the Sparse MAtrix\n")

anno_sparse_data = c()

anno_sparse_data = process_sparse(mat_sparse, mat_row, mat_col)
save(anno_sparse_data,file=paste(mat_dir,"/",anno_type,"_sparse_data",sep=""))

cat("Initializing the Matrices and Indicies\n")
query.IND = matrix_init(anno_sparse_data,21,1)
cat("Saving the Index Matrices\n")
save(query.IND, file=paste(mat_dir,"/","query_ind_mat",sep=""))

cat("Computing the TF-IDF matrix across the cluster\n")
mat.PB.TFIDF = make_cluster(query.IND,1)
save(mat.PB.TFIDF,file=paste(mat_dir,"/","mat.PB.TFIDF",sep=""))

cat("Computing the cosine similarity matrix\n")
PB.sim.mat = mat.PB.TFIDF %*% t(mat.PB.TFIDF)
save(PB.sim.mat,file=paste(mat_dir,"/","PB.sim.mat",sep=""))


