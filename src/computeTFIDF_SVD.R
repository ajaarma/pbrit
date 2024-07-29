library("Matrix")
library("irlba")
require('optparse')
require('methods')

option_list = list(
                make_option(c('-d','--dataDir'),type='character',default=NULL,
                            help='pBRIT data directory',metavar='character'),
                make_option(c('-a','--annoType'),type='character',default=NULL,
                            help='annotation source',metavar='character'),
                make_option(c('-o','--outDir'),type='character',default=NULL,
                            help='output SVD directory',metavar='character')
                  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dataDir)){
    print_help(opt_parser)
    stop('At least one argument should be provided',call.= FALSE)
}


singValD <- function(tfidf_arg,anno_source,out_svd_dir){
	
    obj = load(tfidf_arg)
	X = get(obj)#tfidf_arg[ind]
    message('-- Processing: ',anno_source)
    message('-- Begin estimating SVD PCs')
    if (anno_source == 'PB'|| anno_source=='GO') {

        row_ind = seq(1,dim(X)[1],by=round(dim(X)[1]/5))
        message('-- Splitted rwo intervals: ',row_ind)
        mat_list = list()

        # Split the Pubmed-TFIDF Matrix by row in equal set of 5 for ease 
        # computing zero center mean
        for(i in 1:length(row_ind)){
            if (i != length(row_ind)) {
                message('-- Spliting by row number: ',row_ind[i],' ',row_ind[i+1])
                X.sub = as.matrix(X[row_ind[i]:(row_ind[i+1]-1),])
                X.tmp = apply(X.sub,1,function(x) x - mean(x))
                X.tmp = t(X.tmp)
                mat_list[[i]] = X.tmp
            }else {
                message('-- Spliting by row number: ',row_ind[i],' ',dim(X)[1])
                X.sub = as.matrix(X[row_ind[i]:dim(X)[1],])
                X.tmp = apply(X.sub,1,function(x) x - mean(x))
                X.tmp = t(X.tmp)
                mat_list[[i]] = X.tmp
            }
        }
        
        message(' -- Rbinding the list of matrices')
        X.zmean = do.call(rbind,mat_list)

        message('-- Using irlba package')
        start_time = Sys.time() 
        X.pc = irlba(X.zmean,nv=100,tol=1e-10,maxit=5000)
        end_time = Sys.time()
        time_diff = end_time - start_time
        message('-- Time taken for computing SVD: ',time_diff)
        message('-- Projecting data to PCs')
        X.proj = X %*% X.pc$v

    }else{
        
        X = as.matrix(X)
        X.zmean = t(apply(X,1,function(x) x-mean(x)))
        
        message('-- Using irlba package')
        X.pc = irlba(X.zmean,nv=100,tol=1e-10,maxit=5000)
        message('-- Projecting PCs SVD')
        #X.proj = X %*% X.pc$v[,1:propVarIndex]
        X.proj = X %*% X.pc$v
    }

	message('-- Computing unit norm')
    X.proj.norm = t(apply(X.proj,1,function(x) x/sqrt(sum(x^2))))
    message('-- Computing dot product')
	X.sim = X.proj.norm %*% t(X.proj.norm)
	save(X.pc, file=paste(out_svd_dir,'/',anno_source,".PC.mat",sep=""))
	save(X.sim,file=paste(out_svd_dir,'/',anno_source,'.SVD.sim.mat',sep=''))
}

#Get input directory arguments
mat_path = opt$dataDir
anno_source = opt$annoType
out_svd_dir = opt$outDir

#Get ENSEMBL mapped TFIDF matrices for individual annotation sources
singValD(mat_path,anno_source,out_svd_dir)


#a2 = sort(HP.sim.mat[which(rownames(HP.sim.mat)=="TGFBR2"),],decreasing=TRUE)
#a1 = sort(X.sim[which(rownames(X.sim)=="TGFBR2"),],decreasing=TRUE)

#plot(sv$d^2/sum(sv$d^2),type="b",pch=16,xlab="Singular Values",ylab="Variance explained")
