######################################################################
# Loads all the TFIDF matrix and maps HGNC to ENSEMBL Ids 
#
# Author: Ajay Anand Kumar
#         aakumar1707@gmail.com
#
#
######################################################################
#args = (commandArgs(TRUE))

library("Matrix")
require('methods')
library('optparse')

###############################################
# 
# Map HGNC to ENSG for input TF-IDF Matrix
#
# 
#
#
###############################################
`%ni%` <- Negate(`%in%`)

checkIn <- function(x,hgnc_id){

    check_list = c()
    for (ele in x){
         vals = hgnc_id %in% gsub(' ','',strsplit(ele,',')[[1]])
         check_list = c(check_list,vals)
    }
    return(check_list)
}


getEnsgMap <- function(hgnc.mat.TFIDF,hgnc_ensg_file) {
    
    data_hgnc_ensg = read.table(hgnc_ensg_file,sep='\t',stringsAsFactor=F,na.strings='')
    data_hgnc_ensg[,1] = gsub(' ','',data_hgnc_ensg[,1])
    data_hgnc_ensg[,2] = gsub(' ','',data_hgnc_ensg[,2])

    rownames(hgnc.mat.TFIDF) = toupper(rownames(hgnc.mat.TFIDF))
    ensg_list = c()
    # Know the HGNC gene ids of TFIDF mat present in Master HGNC-ENSEMBL Map
    #hgnc_list = as.data.frame(rownames(hgnc.mat.TFIDF));
    #colnames(hgnc_list) = 'hgnc_id'
   

    hgnc_list_in = rownames(hgnc.mat.TFIDF)[rownames(hgnc.mat.TFIDF) %in% data_hgnc_ensg[,1]]
    hgnc_list_nin = rownames(hgnc.mat.TFIDF)[rownames(hgnc.mat.TFIDF) %ni% data_hgnc_ensg[,1]]
    print(length(hgnc_list_in))
    print(length(unique(hgnc_list_in)))

    # Create an empty Sparse Matrix which is subset of Orignal HGNC-TFIDF Matrix.
    ensg.mat.TFIDF = Matrix(0,nrow=length(hgnc_list_in),ncol=dim(hgnc.mat.TFIDF)[2],sparse=TRUE)
    rownames(ensg.mat.TFIDF) = hgnc_list_in
    colnames(ensg.mat.TFIDF) = colnames(hgnc.mat.TFIDF)
   
    # Copy the HGNC-TFIDF score values for each of the HGNC ids and 
    # keep a track of corresponding ensembl IDs as well
    ensg_list = c()
    for (i in 1:dim(ensg.mat.TFIDF)[1]) {

        hgnc_id = rownames(ensg.mat.TFIDF)[i]
        ind_1 = which(rownames(hgnc.mat.TFIDF)==hgnc_id)
        ind_2 = which(data_hgnc_ensg[,1]==hgnc_id)
       
        if (length(ind_1)!=0){
            vals = as.vector(hgnc.mat.TFIDF[ind_1,])
            #cat('hgnc-id: ',i,'\t',hgnc_id,'\n')
            #cat('indices: ',ind_1,'\t',ind_2,'\n')
            #cat(length(vals),'\t',dim(ensg.mat.TFIDF)[2],'\n')
            ensg.mat.TFIDF[i,] = vals
            ensg_id = as.vector(data_hgnc_ensg[ind_2,2])[1]
            ensg_list = c(ensg_list,ensg_id)
        }else{
            message('Non ID: ',hgnc_id)
        }
    }
    cat('Length hgnc: ',dim(hgnc.mat.TFIDF)[1],'\n')
    cat('Length hgnc_in: ',length(hgnc_list_in),'\n')
    cat('Length ensg_list: ',length(ensg_list),'\n')
    rownames(ensg.mat.TFIDF) = ensg_list
    return(ensg.mat.TFIDF)
}

# Create command line options
option_list = list(
            make_option(c("-i", "--inpMap"), type="character", default=NULL,
                help="Input HGNC-ENSG Map file", metavar="character"),
            make_option(c("-p", "--ppiMap"), type="character", default=NULL,
                help="Input UniprotID-ENSG Map file", metavar="character"),
            make_option(c("-o", "--outDir"), type="character", default=NULL,
                help="Output data directory", metavar="character")
            )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpMap)){
          print_help(opt_parser)
  stop("At least one argument must be input Normalized BLAST score file", call.=FALSE)
}

hgnc_ensg_file = opt$inpMap
ppi_ensg_file = opt$ppiMap
outDir = opt$outDir

anno_list = c('GO','HPO','DO','MPO','GAD','HUGE','PATH','PPI')
#anno_list = c('GAD','HUGE','PATH','PPI')
#anno_list = c('GAD')

# Processing each of the annotation sources of pbrit internal database
# except Pubmed and BLAST as the HGNC genes are already mapped to Ensembl IDs


# Check if the Mapping-Report file already exists
report_file = paste(outDir,'MappingReport.txt',sep='')
if(!file.exists(report_file)){
        system(paste('rm ',report_file,sep=''))
    }
map_mat_dir = paste(outDir,'/MAP/MATRICES/TFIDF/',sep='')
system(paste0('mkdir -p ',map_mat_dir,sep=''))

for (anno_type in anno_list) {

    mat_dir = paste(outDir,'/',anno_type,'/MATRICES/',sep='')
    #map_mat_dir = paste(outDir,'/MAP/MATRICES/TFIDF/',sep='')
    #system(paste0('mkdir -p ',map_mat_dir,sep=''))

    # Check if the Mapping-Report file already exists
    if (anno_type=="GO"){
        cat("\nInside GO Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.GO.TFIDF',sep=''))
        ensg.GO.TFIDF = getEnsgMap(mat.GO.TFIDF,hgnc_ensg_file)
        GO.sim.mat = ensg.GO.TFIDF %*% t(ensg.GO.TFIDF)
        save(ensg.GO.TFIDF,file=paste(map_mat_dir,"/","ensg.GO.TFIDF",sep=""))
        save(GO.sim.mat,file=paste(map_mat_dir,"/","GO.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.GO.TFIDF),'\n',
                                                file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.GO.TFIDF),'\n',
                                                file=report_file,append=TRUE)
    }
    if (anno_type=="HPO"){
        cat("\nInside HPO Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.HP.TFIDF',sep=''))
        ensg.HP.TFIDF = getEnsgMap(mat.HP.TFIDF,hgnc_ensg_file)
        HP.sim.mat = ensg.HP.TFIDF %*% t(ensg.HP.TFIDF)
        save(ensg.HP.TFIDF,file=paste(map_mat_dir,"/","ensg.HP.TFIDF",sep=""))
        save(HP.sim.mat,file=paste(map_mat_dir,"/","HP.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.HP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.HP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }
    if (anno_type=="DO"){
        cat("\nInside DO Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.DO.TFIDF',sep=''))
        ensg.DO.TFIDF = getEnsgMap(mat.DO.TFIDF,hgnc_ensg_file) 
        DO.sim.mat = ensg.DO.TFIDF %*% t(ensg.DO.TFIDF)
        save(ensg.DO.TFIDF,file=paste(map_mat_dir,"/","ensg.DO.TFIDF",sep=""))
        save(DO.sim.mat,file=paste(map_mat_dir,"/","DO.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.DO.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.DO.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }
    if (anno_type=="MPO"){
        cat("\nInside MPO Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.MP.TFIDF',sep=''))
        ensg.MP.TFIDF = getEnsgMap(mat.MP.TFIDF,hgnc_ensg_file) 
        MP.sim.mat = ensg.MP.TFIDF %*% t(ensg.MP.TFIDF)
        save(ensg.MP.TFIDF,file=paste(map_mat_dir,"/","ensg.MP.TFIDF",sep=""))
        save(MP.sim.mat,file=paste(map_mat_dir,"/","MP.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.MP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.MP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }

    if (anno_type=="GAD"){
        cat("\nInside GAD Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.GD.TFIDF',sep=''))
        ensg.GD.TFIDF = getEnsgMap(mat.GD.TFIDF,hgnc_ensg_file) 
        GD.sim.mat = ensg.GD.TFIDF %*% t(ensg.GD.TFIDF)
        save(ensg.GD.TFIDF,file=paste(map_mat_dir,"/","ensg.GD.TFIDF",sep=""))
        save(GD.sim.mat,file=paste(map_mat_dir,"/","GD.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.GD.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.GD.TFIDF),'\n',
                                                    file=report_file,append=TRUE)

    }
    if (anno_type=="HUGE"){
        cat("\nInside HuGE Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.HD.TFIDF',sep=''))
        ensg.HD.TFIDF = getEnsgMap(mat.HD.TFIDF,hgnc_ensg_file) 
        HD.sim.mat = ensg.HD.TFIDF %*% t(ensg.HD.TFIDF)
        save(ensg.HD.TFIDF,file=paste(map_mat_dir,"/","ensg.HD.TFIDF",sep=""))
        save(HD.sim.mat,file=paste(map_mat_dir,"/","HD.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.HD.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.HD.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }
    if (anno_type=="PATH"){
        cat("\nInside Pathway Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.PY.TFIDF',sep=''))
        ensg.PY.TFIDF = getEnsgMap(mat.PY.TFIDF,hgnc_ensg_file) 
        PY.sim.mat = ensg.PY.TFIDF %*% t(ensg.PY.TFIDF)
        save(ensg.PY.TFIDF,file=paste(map_mat_dir,"/","ensg.PY.TFIDF",sep=""))
        save(PY.sim.mat,file=paste(map_mat_dir,"/","PY.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.PY.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.PY.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }
    if (anno_type=="PPI"){
        cat("\nInside PPI Matrix\n",file=report_file,append=TRUE)
        load(paste(mat_dir,'mat.PP.TFIDF',sep=''))
        ensg.PP.TFIDF = getEnsgMap(mat.PP.TFIDF,ppi_ensg_file) 
        PP.sim.mat = ensg.PP.TFIDF %*% t(ensg.PP.TFIDF)
        save(ensg.PP.TFIDF,file=paste(map_mat_dir,"/","ensg.PP.TFIDF",sep=""))
        save(PP.sim.mat,file=paste(map_mat_dir,"/","PP.ENSG.sim.mat",sep=""))
        cat('Dimension original- HGNC: ',dim(mat.PP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
        cat('Dimension mapping - ENSMBL: ',dim(ensg.PP.TFIDF),'\n',
                                                    file=report_file,append=TRUE)
    }
}
