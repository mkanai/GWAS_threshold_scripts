library(data.table)
library(doParallel)
# library(foreach)
# library(pforeach)
library(Rcpp)

setwd(Sys.getenv("workingdir"))
sourceCpp(paste(Sys.getenv("scriptdir"), "Rmeta_analysis.cpp", sep = "/"))

args <- commandArgs(trailingOnly = T)
if (length(args) < 3) {
    stop("Invalid Arguments.")
}

chrs <- c(1:23, 25)
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
sex <- c(".males", ".females")

maf <- args[1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])

nSNP <- as.numeric(read.table(sprintf("analysis/%s/nSNP.txt", maf), colClasses = "numeric")[1,])
nSNP[25] <- nSNP[24]

SNPs.which <- as.list(NA)
invisible(foreach(chr = chrs) %do% {
              SNPs.which[[chr]] <- read.table(sprintf("analysis/%s/chr%d.%s.ALL.bim.which", maf, chr, maf), header = T, colClasses = "numeric")
          })

# cl <- makeCluster(ifelse(detectCores() > length(pops), length(pops), detectCores()))
cl <- makeCluster(6)
registerDoParallel(cl)

input_fmt <- "analysis/%s/%s/chr%d/chr%d.%s.%s.%d.assoc.logistic.P.gz"
ret <- foreach (i = start:end, .combine = "rbind") %do% {
           pvals <- foreach (chr = chrs, .combine = "rbind", .packages = c("foreach", "data.table")) %dopar% {
               foreach (p = pops, .combine = "cbind") %do% {
                   input <- sprintf(input_fmt, maf, p, chr, chr, maf, p, i)
                   df <- as.matrix(fread(sprintf("gunzip -c %s", input), header = T, colClasses = "numeric"))
                   # df <- as.matrix(read.table(gzfile(input), header = T, colClasses = "numeric"))
                   mat <- matrix(NA, nSNP[chr], ncol(df))
                   mat[na.omit(SNPs.which[[chr]][[p]]),] <- df
                   return (mat)
               }
           }
           force(pvals)
           minP <- getMinP_meta(pvals)
           force(minP)
           return (minP)
       }

stopCluster(cl)

pops <- c(pops, "ALL")
colnames(ret) <- c(pops, paste(pops, "GC", sep = "_"), "ALL_GC2", paste(pops, "lambda", sep = "_"))

output <- sprintf("analysis/%s/ALL/genome.%s.ALL.%d-%d.assoc.logistic.P.gz", maf, maf, start, end)
write.table(ret, gzfile(output), sep = "\t", quote = F, row.names = F)

warnings()
