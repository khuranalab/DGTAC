args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
    stop("Usage: Rscript --vanilla batchCorrect.R input output", call.=FALSE)
} else {
    print(args)
}


library(sva)
pwd <- Sys.getenv("PWD")
print(pwd)
################################################################################
## Combat to match ATAC-seq data

# batch correct among all tcga samples
data <- read.csv(args[1],stringsAsFactors = FALSE, row.names = 1,check.names=FALSE)
batchInfo <- read.csv(paste0(pwd,'/.tmp/batchInfo.csv'), row.names = 1,check.names=FALSE)


matrixData <- as.matrix(data)
rm(data)
mod = model.matrix(~as.factor(cancerType), data=batchInfo)

# reference-batch version, with covariates
combat_edata <- ComBat(dat=matrixData, batch=batchInfo$batch, mod=mod, par.prior=TRUE, ref.batch=1)
write.csv(combat_edata, args[2], quote = FALSE)
