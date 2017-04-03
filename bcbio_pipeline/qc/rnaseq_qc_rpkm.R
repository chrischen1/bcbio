# Clustering on parental lines in batch2 and basal cases
rpkm_batch2 <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/20170115/processed_data/rpkm_batch2.tsv", row.names=1, stringsAsFactors=FALSE)
RNAseq.rpkm <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/RNAseq_BrCaPorfiling_201604/processed_data/rpkm_basal.tsv", stringsAsFactors=FALSE,row.names = 1)

rpkm_batch2_basal <- rpkm_batch2[,grep('parent|PR',colnames(rpkm_batch2))]
rpkm_basal <- RNAseq.rpkm[-1,]
colnames(rpkm_basal) <- gsub('_NA','',paste(colnames(RNAseq.rpkm),RNAseq.rpkm[1,],sep = '_'))

genes_i <- intersect(rownames(rpkm_basal),rownames(rpkm_batch2_basal))
rpkm_all <- cbind(rpkm_batch2_basal[genes_i,-6],rpkm_basal[genes_i,])
rpkm_all <- apply(rpkm_all,2,as.numeric)
rpkm_all <- log(rpkm_all+1,2)

cor_m <- cor(rpkm_all,method = 'spearman')
dissimilarity <- 1 - cor_m
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation of RPKM", xlab="")

