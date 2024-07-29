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

rm.zeros <- function(mat){
    # Remove zero rows/columns in the matrix
    
    i1 = which(rowSums(mat)>0)
    i2 = which(colSums(mat)>0)
    return(mat[i1,i2])
}

loadAnnoMat <- function(annoMat,annoType){
    # Function to load individual annotation matrices
    # Turn all negative values to zero
    # Remove all zero rows in the matrix
    
    obj = load(annoMat)
    X.sim = get(obj)
    message('-- Processing: ',annoType)
    
    X.sim[X.sim<0] <- 0
    X.sim = rm.zeros(X.sim)
    return(X.sim)
}

extractSimScore <- function(compMat,outDir,matType) {
    'Function to extract similarity score from the composite matrices X and Y'
    row_names = rownames(compMat)
    col_names = colnames(compMat)
    counter_file = paste(dirname(dirname(dirname(outDir))),'/QOUT/counter_',
                         matType,'.txt',sep='')
    out_file = paste(outDir,'/composite_',matType,'.txt',sep='')
    
    if(file.exists(counter_file)){
        file.remove(counter_file)
    }
    if(file.exists(out_file)){
        file.remove(out_file)
    }

    for(i in 1:length(row_names)) {
        cat(i,'\n',file=counter_file,append=TRUE)
        for(j in i:length(col_names)){
            cat(row_names[i],'\t',col_names[j],'\t',compMat[i,j],'\n',file=out_file,
                append=TRUE)
        }
    }
    gzip(out_file,overwrite=T,remove=T)

}



getCompositeMatrixXY <- function(matDir,method,outDir){
    # Function to load annotation similarity matrix and
    # compute the composite matrices

    matDirList = list.files(matDir)

    annoX_list = c('PP','PB','GO','BL','MP','PY')
    annoY_list = c('HP','HD','GD','DO')


    for (annoType in c(annoX_list,annoY_list)) {
        
        if (method == 'TFIDF') {
            annoMat = paste(matDir,'/',method,'/',annoType,'.ENSG.sim.mat',sep='')
            
        }else{
            annoMat = paste(matDir,'/',method,'/',annoType,'.SVD.sim.mat',sep='')
        }
       
        message(' -- The annotation matrix is: ',annoMat)
        if(annoType == 'PP') {
            PP.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'PB') {
            PB.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'GO') {
            GO.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'BL') {
            annoMat = paste(dirname(dirname(matDir)),'/SEQ/MATRICES/',
                                            annoType,'.ENSG.sim.mat',sep='')
            BL.sim.mat = loadAnnoMat(annoMat,annoType)           
            message(' -- The annotation matrix is: ',annoMat)

        } 
        if(annoType == 'MP') {
            MP.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'PY') {
            PY.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'HP') {
            HP.sim.mat = loadAnnoMat(annoMat,annoType)            
        } 
        if(annoType == 'HD') {
            HD.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'GD') {
            GD.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        if(annoType == 'DO') {
            DO.sim.mat = loadAnnoMat(annoMat,annoType)            
        }
        message(' -- The annotation matrix is: ',annoMat)

    }

    ####### Creating composite Matrices X and Y #####################

    # Get combined row names of X and Y representative annotation 
    # matrices
    gene.names.X = unique(c(rownames(PP.sim.mat),rownames(GO.sim.mat),
                            rownames(BL.sim.mat),rownames(PB.sim.mat),
                            rownames(PY.sim.mat),rownames(MP.sim.mat)
                            )
                         )
    gene.names.Y = unique(c(rownames(HP.sim.mat),rownames(HD.sim.mat),
                            rownames(GD.sim.mat),rownames(DO.sim.mat)
                            )
                         )
    
    # Unique names of Ensembl genes from all the combined matrices
    master_list = unique(c(gene.names.X,gene.names.Y))
    message(' -- After creating unique master list')

    # Create empty matrix
    composite_X = matrix(0,length(master_list),length(master_list),
                       dimnames=list(master_list,master_list)
                        )

    composite_Y=matrix(0,length(master_list),length(master_list),
                       dimnames=list(master_list,master_list)
                       )

    message(" -- After initializing Composite Matrices")

    composite_X[rownames(GO.sim.mat),
                colnames(GO.sim.mat)] = composite_X[rownames(GO.sim.mat),
                                                    colnames(GO.sim.mat)]+
                                                    as.matrix(GO.sim.mat)
    GO.sim.mat = c()
    message(" -- after adding GO similarity matrix to composite matrix")

    composite_X[rownames(PP.sim.mat),
                colnames(PP.sim.mat)] = composite_X[rownames(PP.sim.mat),
                                                    colnames(PP.sim.mat)]+
                                                    as.matrix(PP.sim.mat)
    PP.sim.mat = c()
    message(" -- after adding PP similarity matrix to composite matrix") 
    
    composite_X[rownames(BL.sim.mat),
                colnames(BL.sim.mat)] = composite_X[rownames(BL.sim.mat),
                                                    colnames(BL.sim.mat)]+
                                                    as.matrix(BL.sim.mat)
    BL.sim.mat = c()
    message(" -- after adding BL similarity matrix to composite matrix")

    composite_X[rownames(PY.sim.mat), 
                colnames(PY.sim.mat)] = composite_X[rownames(PY.sim.mat),
                                                    colnames(PY.sim.mat)]+
                                                    as.matrix(PY.sim.mat)
    PY.sim.mat = c()
    message(" -- after adding PY similarity matrix to composite matrix")

    composite_X[rownames(PB.sim.mat), 
                colnames(PB.sim.mat)] = composite_X[rownames(PB.sim.mat),
                                                    colnames(PB.sim.mat)]+
                                                    as.matrix(PB.sim.mat)
    PB.sim.mat = c()
    message(" -- after adding PB similarity matrix to composite matrix")
               
    composite_X[rownames(MP.sim.mat), 
                colnames(MP.sim.mat)] = composite_X[rownames(MP.sim.mat),
                                                    colnames(MP.sim.mat)]+
                                                    as.matrix(MP.sim.mat)
    MP.sim.mat = c()
    message(" -- after adding MP similarity matrix to composite matrix")

    ####### Composite Matrix for Y ##############
    composite_Y[rownames(HP.sim.mat), 
                colnames(HP.sim.mat)] = composite_Y[rownames(HP.sim.mat),
                                                    colnames(HP.sim.mat)]+
                                                    as.matrix(HP.sim.mat)
    HP.sim.mat = c()
    message(" -- after adding HP similarity matrix to composite matrix Y")

    composite_Y[rownames(GD.sim.mat), 
                colnames(GD.sim.mat)] = composite_Y[rownames(GD.sim.mat),
                                                    colnames(GD.sim.mat)]+
                                                    as.matrix(GD.sim.mat)
    GD.sim.mat = c()
    message(" -- after adding GD similarity matrix to composite matrix Y")


    composite_Y[rownames(DO.sim.mat), 
                colnames(DO.sim.mat)] = composite_Y[rownames(DO.sim.mat),
                                                    colnames(DO.sim.mat)]+
                                                    as.matrix(DO.sim.mat)
    DO.sim.mat = c()
    message(" -- after adding DO similarity matrix to composite matrix Y")

    composite_Y[rownames(HD.sim.mat), 
                colnames(HD.sim.mat)] = composite_Y[rownames(HD.sim.mat),
                                                    colnames(HD.sim.mat)]+
                                                    as.matrix(HD.sim.mat)
    HD.sim.mat = c()
    message(" -- after adding HD similarity matrix to composite matrix Y")

    #####################################################################

    ########## Averaging the composite Matrix for X and Y ###############
    composite_X = composite_X/6
    composite_Y = composite_Y/4
    message(" -- after averaging composite Matrix X and Y")

    composite_X[is.na(composite_X)] <- 0
    composite_Y[is.na(composite_Y)] <- 0
    message(" -- after assigning all negative values=>0 to X and Y")

    diag(composite_X) <- 1
    diag(composite_Y) <- 1
    message(" -- after diagonalising composite Matrix X and Y to 1")

    ############## Saving the averaged out composite Matrices X and Y #########
    save(composite_X,file=paste(outDir,'/composite_X',sep=''))
    message(" -- after saving composite matrix X")

    #extractSimScore(composite_X,outDir,'X')
    message(' -- Extracting the Sim score of composite-X matrix in 
                                                compressed Text file')

    save(composite_Y,file=paste(outDir,'/composite_Y',sep=''))
    message(" -- after saving composite Matrix Y")
    #extractSimScore(composite_Y,outDir,'Y')
    message(' -- Extracting the Sim score of composite-Y matrix in 
                                                compressed Text file')
    
}


matDir = opt$dataDir
method = opt$method
outDir = paste(opt$outDir,'/COMP/',method,sep='')
system(paste0('mkdir -p ',outDir))

compMatList = getCompositeMatrixXY(matDir,method,outDir)

