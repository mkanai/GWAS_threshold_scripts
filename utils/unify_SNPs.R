library(foreach)
library(pforeach)
library(data.table)

setwd(Sys.getenv("workingdir"))

args <- commandArgs(trailingOnly = T)
if (length(args) < 1) {
    stop("Invalid Arguments.")
}

chrs <- c(1:23, 25)
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
maf <- args[1]
input_fmt <- sprintf("1KGP/plink/%s/%%s/chr%%d.%s.%%s.bim", maf, maf)
output_fmt <- sprintf("analysis/%s/chr%%d.%s.ALL.bim.which", maf, maf)

nSNP <- foreach(chr = chrs) %do% {
            SNPs <- as.list(NA)
            common <- foreach(p = pops, .combine = "rbind") %do% {
                          input <- sprintf(input_fmt, p, chr, p)
                          dt <- fread(input)
                          SNPs[[p]] <- dt$V2
                          return(dt)
                      }
            setkey(common, V2)
            common.uniq <- unique(common)
            setkey(common.uniq, V4)

            N <- max(sapply(SNPs, function(x){length(x)}))

            SNPs.which <- pforeach(p = pops, .combine = "cbind")({
                              cm <- chmatch(SNPs[[p]], common.uniq$V2)
                              return (c(cm, rep(NA, N - length(cm))))
                          })
            colnames(SNPs.which) <- pops
            write.table(SNPs.which, sprintf(output_fmt, chr), quote = F, row.names = F, sep = "\t")
            return (nrow(common.uniq))
}

write.table(data.frame(nSNP), sprintf("analysis/%s/nSNP.txt", maf), quote = F, col.names = F, row.names = F, sep = "\t")
warnings()
