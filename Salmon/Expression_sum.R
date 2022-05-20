library(tximport)
library(readr)
library(ggplot2)

tx2gene <- read.csv("/Users/aleksandr/Downloads/Gene_expression/Trinity.fasta.gene_trans_map_Arenicola.tsv", sep = "\t")
txi.salmon <- tximport("/Users/aleksandr/Downloads/Gene_expression/quant_Arenicola.sf", type = "salmon", tx2gene = tx2gene)
data_A <- txi.salmon[txi.salmon$abundance >= 1.0]
write.csv(data_A, "Arenicola_genes.csv")

tx2gene <- read.csv("/Users/aleksandr/Downloads/Gene_expression/Trinity.fasta.gene_trans_map_Pygospio.tsv", sep = "\t")
txi.salmon <- tximport("/Users/aleksandr/Downloads/Gene_expression/quant_Pygospio.sf", type = "salmon", tx2gene = tx2gene)
data_P <- txi.salmon[txi.salmon$abundance >= 1.0]
write.csv(data_P, "Pygospio_genes.csv")


data <- data.frame(matrix(unlist(txi.salmon$abundance), nrow=length(txi.salmon$abundance), byrow=TRUE),stringsAsFactors=FALSE)
ggplot(data, aes(x=matrix.unlist.txi)) + geom_histogram()
