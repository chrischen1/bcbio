library(dplyr)
library(reshape2)
library(biomaRt)

combined <- read.delim("/n/scratch2/cc400/rnaseq_bcbio/rnaseq_batch2/samples/final/2017-01-21_samples/combined.sf",as.is = T)
tx2gene <- read.csv("~/Dropbox/_ChrisProject/workspace/RNAseq/processed_data/batch2/tx2gene.csv", header=FALSE, stringsAsFactors=FALSE)
ercc_mapping <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/ercc/cms_095047.txt", stringsAsFactors=FALSE)
# RNAseq_sample <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/RNAseq_sample.tsv", stringsAsFactors=FALSE)


#Counts for transcripts to RPKM for genes
tpm2rpkm <- function(combined,tx2gene,ercc_mapping = NULL){
  gene_mapping <- cbind('transcript'= c(tx2gene$V1,ercc_mapping$GenBank),'gene' = c(tx2gene$V2,ercc_mapping$ERCC_ID))
  genes <- gene_mapping[,2]
  names(genes) <- gene_mapping[,1]
  lib_size <- data.frame('numreads'=combined$numreads,'sample'=combined$sample)
  x <- lib_size %>% group_by(sample) %>% summarise_each(funs(sum))
  scale_factor <- x$numreads/1000000
  names(scale_factor) <- x$sample
  
  combined$RPM <- combined$numreads/scale_factor[combined$sample]
  combined$RPKM <- combined$RPM/(combined$effectiveLength/1000)
  combined$gene <- genes[combined$id]
  
  rpkm_combined <- data.frame('sample'=combined$sample,'gene'=combined$gene,'RPKM'=combined$RPKM)
  rpkm_combined_gene <- rpkm_combined %>% group_by(sample,gene)%>% summarise_each(funs(sum))
  
  rpkm_raw <- acast(rpkm_combined_gene,gene~sample)
}


rpkm_raw <- tpm2rpkm(combined,tx2gene,ercc_mapping)
colnames(rpkm_raw) <- gsub('.+_S(\\d+)_R.+','\\1',colnames(rpkm_raw))

# transform ensumbl id to gene symbol
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(rpkm_raw), mart = mart)
ercc <- grep('^ERCC',rownames(rpkm_raw),value = T)
ens_mapping <- results[results$hgnc_symbol != '' & !(results$hgnc_symbol %in% results$hgnc_symbol[duplicated(results$hgnc_symbol)]) ,]
ens_mapping2 <- rbind(ens_mapping,cbind('ensembl_gene_id'=ercc,'hgnc_symbol'=ercc))
# write.csv(ens_mapping2,'~/Dropbox/_ChrisProject/workspace/bcbio/annotations/ens_mapping.csv',row.names = F)

# ens_mapping2 <- read.csv('~/Dropbox/_ChrisProject/workspace/bcbio/annotations/ens_mapping.csv',as.is = T)
gene_symbols <- ens_mapping2$hgnc_symbol
names(gene_symbols) <- ens_mapping2$ensembl_gene_id

rpkm <- rpkm_raw[rownames(rpkm_raw) %in% ens_mapping2$ensembl_gene_id,]
rownames(rpkm) <- gene_symbols[rownames(rpkm)]

write.table(rpkm,'rpkm.tsv',sep = '\t',quote = F)

# add column infomation
# cells <- RNAseq_sample$cellines
# tags <- RNAseq_sample$rep_tag
# names(cells) <- names(tags) <- RNAseq_sample$sample_id
# rpkm_lincs <- rbind('Cellline'=cells[colnames(rpkm)],'tag'=tags[colnames(rpkm)],rpkm)
# 
# RNAseq_rpkm <- NULL
# data_rpkm <- rpkm_lincs[-(1:2),]
# for(i in unique(RNAseq_sample$cellines)){
#   selected_samples <- RNAseq_sample$sample_id[RNAseq_sample$cellines==i]
#   if(length(selected_samples)==1){
#     if(is.null(RNAseq_rpkm)){
#       RNAseq_rpkm <- as.numeric(data_rpkm[,selected_samples])
#     }else{
#       RNAseq_rpkm <- cbind(RNAseq_rpkm,as.numeric(data_rpkm[,selected_samples]))
#     }
#   }else{
#     avg_col <- apply(data_rpkm[,selected_samples],1,function(x)mean(as.numeric(x)))
#     if(is.null(RNAseq_rpkm)){
#       RNAseq_rpkm <- avg_col
#     }else{
#       RNAseq_rpkm <- cbind(RNAseq_rpkm,avg_col)
#     }
#   }
# }
# RNAseq_rpkm <- log(RNAseq_rpkm+1,2)
# RNAseq_rpkm2 <- rbind(c('id',unique(RNAseq_sample$cellines)),cbind(rownames(RNAseq_rpkm),RNAseq_rpkm))
# rownames(rpkm_lincs)[1]='id'
# write.table(rpkm_lincs,'~/Dropbox/_ChrisProject/workspace/bcbio/processed_data/basal/rpkm_lincs.tsv',sep = '\t',quote = F,col.names =F)
# write.table(RNAseq_rpkm2,'~/Dropbox/_ChrisProject/workspace/bcbio/processed_data/basal/RNAseq-rpkm.tsv',sep = '\t',quote = F,col.names =F,row.names = F)

#QC
# rpkm_old <- read.delim("~/rpkm.tsv", row.names=1, stringsAsFactors=FALSE)
# rpkm_new <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/processed_data/basal/RNAseq-rpkm.tsv", row.names=1, stringsAsFactors=FALSE)
# gene_i <- intersect(rownames(rpkm_old),rownames(rpkm_new))
# 
# par(mfrow=c(1,2))
# smoothScatter(log(as.numeric(unlist(rpkm_new[gene_i,]))+0.1),log(as.numeric(unlist(rpkm_old[gene_i,colnames(rpkm_new)]))+0.1),xlab = 'bcbio result',ylab='DIY result')
# 
# plot(density(log(as.numeric(unlist(rpkm_old[gene_i,]))+0.1)),main='Log RPKM')
# lines(density(log(as.numeric(unlist(rpkm_new[gene_i,colnames(rpkm_old)]))+0.1)),col='red')
# legend('topright',c('bcbio result','DIY result'),col=c(2,1),lty=1)
# 



