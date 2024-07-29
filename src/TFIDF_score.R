######################################################################
# TFIDF calculation for all annotation sources. 
# It employs usage of "SNOW" package for parallelizing the computation
#
# Author: Ajay Anand Kumar
#         aakumar1707@gmail.com
#
# Copyright: Plant Stress Bioinformatics Group                       
#
######################################################################
args = (commandArgs(TRUE))

library("Matrix")
#load("go_sparse_mat")
#source("process_sparse.R")
library("snow")
library(hash)
#load("col_hash_obj")


process_sparse <- function(sparse_data_arg,sparse_row_arg,sparse_col_arg){
	sparse_matr = read.table(sparse_data_arg,sep="\t",header=FALSE,quote="")
  	sparse_row = read.table(sparse_row_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")
  	sparse_col = read.table(sparse_col_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")

  	sz_rownames = dim(sparse_row)[1]
  	sz_colnames = dim(sparse_col)[1]

  	if(max(sparse_matr[,1]) > sz_rownames[1]){
		print("row names don't fit the matrix. Too few row names")
  	}
  
  	if(max(sparse_matr[,2]) > sz_colnames[1]){
		print("col names don't fit the matrix. Too few col names")
  	}

  	data <- sparseMatrix(i = sparse_matr[,1], j = sparse_matr[,2],
                       dims = c(sz_rownames[1], sz_colnames[1]))

  	rownames(data) = sparse_row[,2]
  	colnames(data) = sparse_col[,1]
  	out <- list(sparse_data = data, all_colnames = sparse_col, all_rownames = sparse_row)
  	return(out)
	}


process_sparsePPI <- function(sparse_data_arg,sparse_row_arg,sparse_col_arg){

    sparse_matr = read.table(sparse_data_arg,sep="\t",header=FALSE,quote="")
    sparse_row = read.table(sparse_row_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")
    sparse_col = read.table(sparse_col_arg,sep="\t",header=FALSE,quote="",comment.char="",colClasses="character")

    sz_rownames = dim(sparse_row)[1]
    sz_colnames = dim(sparse_col)[1]

    if(max(sparse_matr[,1]) > sz_rownames[1])            {
        print("row names don't fit the matrix. Too few row names")
    }

    if(max(sparse_matr[,2]) > sz_colnames[1]){
        print("col names don't fit the matrix. Too few col names")
    }

    data <- sparseMatrix(i = sparse_matr[,1], j = sparse_matr[,2],x=sparse_matr[,4],
                                        dims = c(sz_rownames[1], sz_colnames[1]))

    data[is.na(data)] <- 0
    rownames(data) = sparse_row[,2]
    colnames(data) = sparse_col[,1]
    out <- list(sparse_data = data, all_colnames = sparse_col, all_rownames = sparse_row)
    return(out)
}

myFunc <- function(query_list){
	library(hash)
  	library("Matrix")
  	row.ind = query_list[[1]]
  	mat_new = query_list[[2]][query_list[[1]],]
  	#cat("after mat new\n",file="console.txt",append=TRUE)
  	#row.name = query_list[[2]][query_list[[1]]]
  	#cat(length(query_list[[1]]),"\n",file="console.txt",append=TRUE)
  	#mat = sparseMatrix(row.name,query_list[[3]],dims = c(length(query_list[[1]]),length(query_list[[3]])))
  	#cat("After sparse mat",file="console.txt",append=TRUE)
  	#mat = query_list[[2]]
  	nrow = dim(query_list[[2]])[1]
  	col_hash = query_list[[3]]
  	for(ele in 1:length(row.ind)){
		cat(ele,"\n",file="console.txt",append=TRUE)
    	ind = which(mat_new[ele,] !=0)
    		#termFreq = 1/length(ind)
    
    	for(i in ind){
       		termFreq = 1+log(mat_new[ele,i])
            #termFreq = mat_new[ele,i]/length(ind)
            idf = log(nrow/(as.numeric(values(col_hash,keys=i))+1))
            #cat("TFIDF value:",tfidf,"\n",file="console.txt",append=TRUE)
            mat_new[ele,i]=idf*termFreq
        }
        normal = sqrt(sum(mat_new[ele,]^2))
        if(normal !=0){
            #cat("Normal value:",normal,"\n",file="console.txt",append=TRUE)
            mat_new[ele,] = mat_new[ele,]/normal
        }
    }
    return(mat_new)
}

matrix_init <- function(mat_obj, num_node, mode){#, rowInd=rowInd){
    mat = mat_obj[[1]]
    if(mode ==1){

        mat <- mat_obj[[1]]
        colnames(mat) = gsub(" ","",mat_obj[[2]][,1],fixed=TRUE)
        rownames(mat) = gsub(" ","",rownames(mat),fixed=TRUE)
        #mat.old = mat
        mat = mat*1
        nrow = dim(mat)[1]
        ncol = dim(mat)[2]
    }
    else{
        mat = sparseMatrix(length(colnames(mat)),length(colnames(mat)))
        rownames(mat) = colnames(mat)
        colnames(mat) = colnames(mat)
        mat = mat*1
        nrow = dim(mat)[1]
        ncol = dim(mat)[2]
    }
    col_hash = hash()
  
    for(ele in 1:ncol){
        print(ele)
        if (mode ==1){
            ind = as.numeric(which(mat[,ele]!=0))
            col_hash[ele]=length(ind)
        }
        else{
            ind = which(mat[,ele]!=0)
            col_hash[ele] = as.numeric(ind)
        }
    }
    query.ind = list()
  #seq.num = seq(1,nrow,round(nrow/30))
    seq.num = seq(1,nrow,round(nrow/num_node))
    print(seq.num)
    #print(nrow)
    for(i in 1:length(seq.num)){
        if(i != length(seq.num)){
            query.ind[[i]] = list(c(seq.num[i]:(seq.num[i+1]-1)),mat,col_hash)
        }
        else {
            query.ind[[i]] = list(c(seq.num[i]:nrow),mat,col_hash)
        }
    }
    return(query.ind)
}

make_cluster <- function(query.ind,mode,path){

    cl = makeCluster(length(query.ind),type="SOCK")
    #cl = makeCluster(15,type="SOCK")
    cat("after clusater")
    tryCatch({
    #print("before cl")
    if(mode==1){
        result = clusterApply(cl,query.ind,myFunc)
    }
    else{
        print("Inside else, make cluster")
        result = clusterApply(cl,query.ind,myInfoGain)
    }
    },interrupt = function(ex){print(ex)},
    error = function(ex){print(ex)})  
    stopCluster(cl)
    mat.final = result[[1]]
    save(result,file=paste(path,"/CL_result_mat",sep=""))
    for(i in 2:length(query.ind)){
        mat.final = rBind(mat.final,result[[i]])

    }
    return(mat.final)
}

hist_plot = function(obs_mat,ind){
    ind1 = which(obs_mat[ind,]!=0)
    y = as.vector(obs_mat[ind,ind1])
    x = c(1:length(y))
    #plot(x,y,xlab=colnames(obs_mat),ylab="TFIDF values")
    ##plot(x,y,xlab=colnames(mat.PB.final)[which(mat.PB.final[9665,]!=0)],ylab="TFIDF Values",col="red")
    plot(x,y,xlab=colnames(obs_mat)[ind1],ylab="TFIDF Values",col="red")
    text(x,y,colnames(obs_mat)[ind1],srt=90,cex=0.8,pos=1,xpd=TRUE)
    ##text(x,y,colnames(mat.PB.final)[which(mat.PB.final[9665,]!=0)],srt=90,cex=0.8,pos=1,xpd=TRUE)

}

    
mat_sparse = args[[1]] #"/home/shared_data_medgen_aorta/pbrit_data/V3/GO/INPUT/HGNC_sparse.txt"
mat_row = args[[2]] #"/home/shared_data_medgen_aorta/pbrit_data/V3/GO/INPUT/HGNC_rownames.txt"
mat_col = args[[3]] #"/home/shared_data_medgen_aorta/pbrit_data/V3/GO/INPUT/HGNC_colnames.txt"
mat_dir = args[[4]]
anno_type = args[[5]]


cat("Loading the Sparse Matrix\n")

anno_sparse_data = c()
if (anno_type=="PPI"){
    anno_sparse_data = process_sparsePPI(mat_sparse, mat_row, mat_col)
}else{
    anno_sparse_data = process_sparse(mat_sparse, mat_row, mat_col)
}

#save(GO_sparse_data,file="/home/shared_data_medgen_aorta/pbrit_data/V3/GO/MATRICES/GO_sparse_data")
save(anno_sparse_data,file=paste(mat_dir,"/",anno_type,"_sparse_data",sep=""))
#query.IND = matrix_init(GO_sparse_data,15,1)
cat("Initializing the Matrices and Indices\n")
query.IND = matrix_init(anno_sparse_data,16,1)
cat("Saving the Index Matrices\n")
save(query.IND, file=paste(mat_dir,"/","query_ind_mat",sep=""))

mat.TFIDF = c()
sim.mat = c()

cat("Computing the TF-IDF matrix across the cluster\n")
mat.TFIDF = make_cluster(query.IND,1,mat_dir)
cat("Computing the cosine similarity matrix\n")
sim.mat = mat.TFIDF %*% t(mat.TFIDF)


if (anno_type=="GO"){
    cat("Inside GO Matrix\n")
    mat.GO.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    GO.sim.mat = sim.mat
    save(mat.GO.TFIDF,file=paste(mat_dir,"/","mat.GO.TFIDF",sep=""))
    save(GO.sim.mat,file=paste(mat_dir,"/","GO.sim.mat",sep=""))
}

if (anno_type=="HPO"){
    cat("Inside HPO Matrix\n")
    mat.HP.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    HP.sim.mat = sim.mat
    save(mat.HP.TFIDF,file=paste(mat_dir,"/","mat.HP.TFIDF",sep=""))
    save(HP.sim.mat,file=paste(mat_dir,"/","HP.sim.mat",sep=""))
}
if (anno_type=="DO"){
    cat("Inside DO Matrix\n")
    mat.DO.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    DO.sim.mat = sim.mat
    save(mat.DO.TFIDF,file=paste(mat_dir,"/","mat.DO.TFIDF",sep=""))
    save(DO.sim.mat,file=paste(mat_dir,"/","DO.sim.mat",sep=""))
}
if (anno_type=="MPO"){
    cat("Inside MPO Matrix\n")
    mat.MP.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    MP.sim.mat = sim.mat
    save(mat.MP.TFIDF,file=paste(mat_dir,"/","mat.MP.TFIDF",sep=""))
    save(MP.sim.mat,file=paste(mat_dir,"/","MP.sim.mat",sep=""))
}
if (anno_type=="GAD"){
    cat("Inside GAD Matrix\n")
    mat.GD.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    GD.sim.mat = sim.mat
    save(mat.GD.TFIDF,file=paste(mat_dir,"/","mat.GD.TFIDF",sep=""))
    save(GD.sim.mat,file=paste(mat_dir,"/","GD.sim.mat",sep=""))
}

if (anno_type=="HUGE"){
    cat("Inside HuGE Matrix\n")
    mat.HD.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    HD.sim.mat = sim.mat
    save(mat.HD.TFIDF,file=paste(mat_dir,"/","mat.HD.TFIDF",sep=""))
    save(HD.sim.mat,file=paste(mat_dir,"/","HD.sim.mat",sep=""))
}

if (anno_type=="PATH"){
    cat("Inside Pathway Matrix\n")
    mat.PY.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    PY.sim.mat = sim.mat
    save(mat.PY.TFIDF,file=paste(mat_dir,"/","mat.PY.TFIDF",sep=""))
    save(PY.sim.mat,file=paste(mat_dir,"/","PY.sim.mat",sep=""))
}

if (anno_type=="PPI"){
    cat("Inside PPI Matrix\n")
    mat.PP.TFIDF = mat.TFIDF #make_cluster(query.IND,1)
    PP.sim.mat = sim.mat
    save(mat.PP.TFIDF,file=paste(mat_dir,"/","mat.PP.TFIDF",sep=""))
    save(PP.sim.mat,file=paste(mat_dir,"/","PP.sim.mat",sep=""))
}
