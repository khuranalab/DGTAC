
library(sva)

pwd <- Sys.getenv("PWD")
print(pwd)
################################################################################
## Combat to match ATAC-seq data

# batch correct among all tcga samples
data <- read.csv(paste0(pwd,'/.tmp/expTpmMerged.csv'),stringsAsFactors = FALSE, row.names = 1,check.names=FALSE)
batchInfo <- read.csv(paste0(pwd,'/.tmp/batchInfo.csv'), row.names = 1,check.names=FALSE)


matrixData <- as.matrix(data)
rm(data)
mod = model.matrix(~as.factor(cancerType), data=batchInfo)

# reference-batch version, with covariates
combat_edata <- ComBat(dat=matrixData, batch=batchInfo$batch, mod=mod, par.prior=TRUE, ref.batch=1)
write.csv(combat_edata, paste0(pwd,'/.tmp/expTpmAdjusted.csv'), quote = FALSE)
