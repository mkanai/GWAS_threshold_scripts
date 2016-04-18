library(ggplot2)
library(gridExtra)
library(reshape2)
library(WRShd)

width <- 8
height <- 8
dpi <- 300

cores <- 20

args <- commandArgs(trailingOnly = T)

maf <- args[1]
start <- as.numeric(args[2])
end <- as.numeric(args[3])
s <- ifelse(length(args) > 3, args[4], "")

input <- sprintf("analysis/%s/ALL%s/genome.%s.ALL%s.%d-%d.assoc.logistic.P.gz", maf, s, maf, s, start, end)
df <- read.table(gzfile(input), header = T, colClasses = "numeric")
npop <- (ncol(df) - 1) / 3


pvals.log10 <- -log10(df[1:(npop*2 + 1)])
molten <- melt(pvals.log10[, 1:npop], variable.name = "population")
molten.GC <- melt(pvals.log10[, (npop+1):ncol(pvals.log10)], variable.name = "population")

molten.95 <- melt(data.frame(t(apply(pvals.log10[, 1:npop], 2, function(x) {hd(x, q = .95, cores = cores)}))), value.name = "q95", variable.name = "population")
molten.95.ci <- cbind(as.factor(colnames(pvals.log10))[1:npop], data.frame(t(apply(pvals.log10[, 1:npop], 2, function(x) {hdpb(x, q = .95, cores = cores)$ci}))))
colnames(molten.95.ci) <- c("population", "low", "high")
molten.GC.95 <- melt(data.frame(t(apply(pvals.log10[, (npop+1):ncol(pvals.log10)], 2, function(x) {hd(x, q = .95, cores = cores)}))), value.name = "q95", variable.name = "population")
molten.GC.95.ci <- cbind(as.factor(colnames(pvals.log10))[(npop+1):ncol(pvals.log10)], data.frame(t(apply(pvals.log10[, (npop+1):ncol(pvals.log10)], 2, function(x) {hdpb(x, q = .95, cores = cores)$ci}))))
colnames(molten.GC.95.ci) <- c("population", "low", "high")

plt <- ggplot(molten, aes(x = value, fill = population, color = population)) +
           stat_density(aes(ymax = ..density.., ymin = ..density..)) +
           xlim(0, 10) + xlab("-log10 P") + #ylim(0, 1) +
           theme(legend.position = "bottom") +
           facet_grid(. ~ population) +
           geom_vline(x = c(-log10(5e-8))) +
           geom_vline(data = molten.95, aes(xintercept = q95, color = population), linetype = "dashed") +
           geom_segment(data = molten.95.ci, aes(x = low, xend = low, y = .25, yend = .75, color = population), alpha = .8) +
           geom_segment(data = molten.95.ci, aes(x = high, xend = high, y = .25, yend = .75, color = population), alpha = .8) +
           geom_segment(data = molten.95.ci, aes(x = low, xend = high, y = .5, yend = .5, color = population), alpha = .8) +
           coord_flip() #+
#           annotation_custom(tableGrob(molten.95), xmin = 0, xmax = 2.5, ymin = 0.75, ymax = 1)

plt.GC <- ggplot(molten.GC, aes(x = value, fill = population, color = population)) +
              stat_density(aes(ymax = ..density.., ymin = ..density..)) +
              xlim(0, 10) + xlab("-log10 P") + #ylim(0, 1) +
              theme(legend.position = "bottom") +
              facet_grid(. ~ population) +
              geom_vline(x = c(-log10(5e-8))) +
              geom_vline(data = molten.GC.95, aes(xintercept = q95, color = population), linetype = "dashed") +
              geom_segment(data = molten.GC.95.ci, aes(x = low, xend = low, y = .25, yend = .75, color = population), alpha = .8) +
              geom_segment(data = molten.GC.95.ci, aes(x = high, xend = high, y = .25, yend = .75, color = population), alpha = .8) +
              geom_segment(data = molten.GC.95.ci, aes(x = low, xend = high, y = .5, yend = .5, color = population), alpha = .8) +
              coord_flip() #+
#              annotation_custom(tableGrob(molten.GC.95), xmin = 0, xmax = 3, ymin = 0.75, ymax = 1)


lambda <- df[(npop*2 + 2):ncol(df)]
molten <- melt(lambda, variable.name = "population")
levels(molten$population) <- gsub("_lambda", "", levels(molten$population))
plt.lambda <- ggplot(molten, aes(x = value, y = ..density.., fill = population, color = population)) +
                  geom_density(alpha = .2) + xlab("lambda") +
                  theme(legend.position = "bottom")

output_fmt <- sprintf("analysis/%s/ALL%s/genome.%s.ALL%s.%d-%d.VOLCANO%%s.%%s", maf, s, maf, s, start, end)

for (ext in c("pdf", "png")) {
    ggsave(sprintf(output_fmt, "", ext), plt, width = width, height = height, dpi = dpi)
    ggsave(sprintf(output_fmt, ".GC", ext), plt.GC, width = width, height = height, dpi = dpi)
    ggsave(sprintf(output_fmt, ".lambda", ext), plt.lambda, width = width, height = height, dpi = dpi)
}

write.table(rbind(cbind(molten.95, molten.95.ci[,2:3]), cbind(molten.GC.95, molten.GC.95.ci[,2:3])), sprintf(output_fmt, ".q95", "tsv"), quote = F, row.names = F, sep = "\t")

