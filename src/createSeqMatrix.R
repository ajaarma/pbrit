#### Script to create Sequence similarity Matrix #########
#### 
#### 
###/#####################################################

require('optparse')
require('Matrix')
require('reshape2')
require('data.table')

option_list = list(
           make_option(c("-i", "--inpBit"), type="character", default=NULL,
                help="Input Normalized Bit Score", metavar="character"),
           make_option(c("-o", "--outDir"), type="character", default=NULL,
                help="Output Matrices directory", metavar="character"),
           make_option(c("-m", "--outMat"), type="character", default=NULL,
                help="Output Matrix File", metavar="character")
           )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpBit)){
      print_help(opt_parser)
  stop("At least one argument must be input Normalized BLAST score file", call.=FALSE)
}

print(opt)

mat_norm_data = fread(opt$inpBit,header=F,stringsAsFactors=F)
matFrame = data.frame(mat_norm_data)

#stop()
outDir = opt$outDir
outMat = opt$outMat

matFrame[,2] = gsub(' ','',matFrame[,2])
matFrame[,3] = gsub(' ','',matFrame[,3])

#matFrame = matFrame[,c(2:4)]

blastMat = acast(matFrame,V2~V3,value.var='V4')
blastMat[is.na(blastMat)] <- 0

blastMat = blastMat*1
message(' -- Finished creating the BLAST Matrix')

save(blastMat,file=paste(outDir,'/MATRICES/blastMatFull',sep=''))
message(' -- Finished saving the BLAST Full Matrix')

blastMatUnitNorm = t(apply(
                        blastMat,1,function(x){
                                tmp = sqrt(sum(x^2))
                                y = x/tmp
                                return(y)
                                }
                          )
                    )
message(' -- Finished creating the BLAST Unit Norm Matrix')

save(blastMatUnitNorm,file=paste(outDir,'/MATRICES/blastMatUnitNorm',sep=''))
message(' -- Finished creating the BLAST Matrix')

blastMatUnitNormSparse = Matrix(blastMatUnitNorm,sparse=TRUE)
message(' -- Finished creating unit Norm Sparse BLAST Matrix')

BL.sim.mat = blastMatUnitNormSparse %*% t(blastMatUnitNormSparse)
save(BL.sim.mat,file=outMat)
message(' -- Finished computing dot product similarity Matrix')


