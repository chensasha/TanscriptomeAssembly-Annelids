library(tximport)
library(readr)

tx2gene <- read.csv("PATH/TO/Trinity.fasta.gene_trans_map.tsv", sep = "\t")
txi.salmon <- tximport("PATH/TO/quant.sf", type = "salmon", tx2gene = tx2gene)
data <- txi.salmon[txi.salmon$abundance >= 1.0]
write.csv(data, "file.csv")
