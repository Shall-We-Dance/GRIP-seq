library (ggplot2)

args <- commandArgs(trailingOnly=TRUE)
DIST1 <- args[1]
DIST2 <- args[2]
DIST3 <- args[3]
DIST4 <- args[4]
OUT <- args[5]
m6a.dist <- read.delim (DIST1, header = T)

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
m6a1.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)




m6a.dist <- read.delim (DIST2, header = T)

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
m6a2.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)


m6a.dist <- read.delim (DIST3, header = T)

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
m6a3.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)



m6a.dist <- read.delim (DIST4, header = T)

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
m6a4.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)



metagene.cord <- c(m6a1.metagene.coord, m6a2.metagene.coord, m6a3.metagene.coord, m6a4.metagene.coord)
mod <- c(rep("YTH-WT", length(m6a1.metagene.coord)), 
        rep("YTH-397SFY_1", length(m6a2.metagene.coord)),
        rep("YTH-397SFY_2", length(m6a3.metagene.coord))
        rep("YTH-397SFY_3", length(m6a4.metagene.coord))) 
df <- data.frame(metagene.cord, mod)

ggplot(df) + geom_density(aes(x = metagene.cord, colour = mod)) + xlim(0, 3) + 
  theme_bw() + geom_vline(xintercept = 1:2, col = "grey")

ggsave(pout, file=OUT,width=16,height=12)
