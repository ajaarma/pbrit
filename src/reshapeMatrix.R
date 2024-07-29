require('Matrix')
require('optparse')
require('methods')
require('R.utils')

option_list = list(
                    make_option(c('-d','--dataDir'),type='character',default=NULL,
                        help='pBRIT ensembl mapped data directory - TFIDF/SVD',
                        metavar='character'),
                    make_option(c('-m','--method'),type='character',default=NULL,
                        help='analysis method - TFIDF/SVD',
                        metavar='character'),
                    make_option(c('-o','--outDir'),type='character',default=NULL,
                        help='output composite matrices - TFIDF/SVD',
                        metavar='character')

                  )


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$dataDir)){
    print_help(opt_parser)
    stop('At least one argument should be provided',call.= FALSE)
}

reshapeCompMat <- function(compMatFile,outFile,matType) {
    
    'Function to extract similarity score from the composite matrices X and Y'

    obj = load(compMatFile)
    compMat = get(obj)

    row_names = rownames(compMat)
    col_names = colnames(compMat)
    counter_file = paste(dirname(dirname(dirname(dirname(compMatFile)))),'/QOUT/counter_',
                                                      matType,'.txt',sep='')
    #out_file = paste(outDir,'/composite_',matType,'.txt',sep='')

    if(file.exists(counter_file)){
        file.remove(counter_file)
    }
    if(file.exists(outFile)){
        file.remove(outFile)
    }

    for(i in 1:length(row_names)) {
        cat(i,'\n',file=counter_file,append=TRUE)
        for(j in i:length(col_names)){
            cat(row_names[i],'\t',col_names[j],'\t',compMat[i,j],'\n',file=outFile,
                 append=TRUE)
            }
        }
    gzip(outFile,overwrite=T,remove=T)
}

reshapeTfidfMat <- function(tfidfMat,outFile,matType) {
    'Function to extract TFIDF score from the TFIDF Matrices if the 
     Annotation sources ' 

    obj = load(tfidfMat)
    tfidfMat = get(obj)

    row_names = rownames(tfidfMat)
    col_names = colnames(tfidfMat)

    for(i in 1:length(row_names)) {
        tfidf_mat_vals = tfidfMat[i,tfidfMat[i,]!=0]

        for(j in 1:length(tfidf_mat_vals)){
            cat(row_names[i],'\t',names(tfidf_mat_vals)[j],'\t',tfidf_mat_vals[j],
                '\n',file=outFile, append=TRUE)
        }
    }
    
    gzip(outFile,overwrite=T,remove=T)

}

compMat = opt$dataDir
method  = opt$method
outDir  = opt$outDir

if (grepl('ANNO',method)) {
    reshapeTfidfMat(compMat,outDir,method)
}else {
    reshapeCompMat(compMat,outDir,method)
}
