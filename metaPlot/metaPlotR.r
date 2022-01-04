library (ggplot2)

args <- commandArgs(trailingOnly=TRUE)
DIST <- args[1]
OUT <- args[2]
m6a.dist <- read.delim (DIST, header = T)

# Determine longest length transcript for each gene
trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
colnames(temp) <- c("gene_name", "gid", "trx_len") 
temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
temp.df <- temp[!duplicated(temp$gene_name),]

# limit m6a data to one transcript per gene (longest)
m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]

utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
## Determine scale factor
utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]

# rescale 5'UTR and 3'UTR
library("scales")
utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))


# Combine and plot
## Histogram
m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
pout <- qplot(m6a.metagene.coord, geom="density") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
ggsave(pout, file=OUT,width=16,height=12)
