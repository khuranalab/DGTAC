args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
    stop("Usage: Rscript --vanilla batchCorrect.R input output", call.=FALSE)
} else {
    print(args)
}


library(sva)
library(preprocessCore)

pwd <- Sys.getenv("PWD")
print(pwd)
################################################################################
## Combat to match ATAC-seq data

# batch correct among all tcga samples
data <- read.csv(args[1],stringsAsFactors = FALSE, row.names = 1,check.names=FALSE)
batchInfo <- read.csv(paste0(pwd,'/.tmp/batchInfo.csv'), row.names = 1,check.names=FALSE)


sampleBin <- data[,row.names(subset(batchInfo, batch != 1))]
sampleBinNorm <- normalize.quantiles(as.matrix(sampleBin),copy=TRUE)
sampleBinNorm_df <- data.frame(sampleBinNorm)
rownames(sampleBinNorm_df) <- rownames(sampleBin)
colnames(sampleBinNorm_df) <- colnames(sampleBin)

data[,row.names(subset(batchInfo, batch != 1))] <- sampleBinNorm_df

correctPeaks <- rownames(sampleBinNorm_df)[rowSums(sampleBinNorm_df) != 0]
length(correctPeaks)

matrixData <- as.matrix(data)
rm(data)
mod = model.matrix(~as.factor(cancerType), data=batchInfo)

# reference-batch version, with covariates
combat_edata <- ComBat(dat=matrixData[correctPeaks,], batch=batchInfo$batch, mod=mod, par.prior=T, ref.batch=1)
matrixData[correctPeaks,] <- combat_edata

# combat_edata[,row.names(subset(batchInfo, batch == 2))]
write.csv(matrixData, args[2], quote = FALSE)
