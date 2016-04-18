library(data.table)
library(doParallel)
library(Rcpp)

setwd(Sys.getenv("workingdir"))
sourceCpp(paste(Sys.getenv("scriptdir"), "Rmeta_analysis.cpp", sep = "/"))

args <- commandArgs(trailingOnly = T)

if (length(args) < 3) {
    stop("Invalid Arguments.")
}

pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
sex <- c("males", "females")

maf <- args[1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])

cl <- makeCluster(2)
registerDoParallel(cl)

input_fmt <- "analysis/%s/%s/chr23/chr23.%s.%s.%s.%d.assoc.logistic.P.gz"
invisible(foreach (i = start:end) %do% {
    foreach (p = pops) %do% {
        output <- sprintf("analysis/%s/%s/chr23/chr23.%s.%s.%d.assoc.logistic.P.gz", maf, p, maf, p, i)
        pvals.sex <- foreach (s = sex, .combine = "cbind", .packages = "data.table") %dopar% {
            input <- sprintf(input_fmt, maf, p, maf, p, s, i)
            return (as.matrix(fread(sprintf("gunzip -c %s", input), header = T, colClasses = "numeric")))
            # read.table(gzfile(input), header = T, colClasses = "numeric")
        }
        force(pvals.sex)
        ret <- getMetaP(pvals.sex)
        force(ret)
        colnames(ret) <- c("BETA", "SE", "P")
        write.table(ret, gzfile(output), sep = "\t", quote = F, row.names = F)
    }
})

stopCluster(cl)
