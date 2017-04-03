library(biomaRt)
library(RUVSeq)
cnt_ref <- read.delim("~/Downloads/140331_VSDexport.txt", row.names=1, stringsAsFactors=FALSE)
cnt_exp <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/processed_data/counts.tsv", row.names=1, stringsAsFactors=FALSE)

cell_info <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/cellines_info/celllines_info_Genentech.txt",stringsAsFactors = F)
cell_info_exp <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/annotation/RNAseq_sample.tsv",stringsAsFactors = F)

cell_exp <- cell_info_exp$cellines
colnames(cnt_ref) <- toupper(gsub('[\\.|-| ]','',colnames(cnt_ref)))
cell_intersect=intersect(cell_exp,colnames(cnt_ref))
cnt_ref_intersect <- cnt_ref[,cell_intersect]

col_name_exp <- paste(cell_info_exp$cellines,cell_info_exp$rep_tag)
names(col_name_exp) <- cell_info_exp$sample_id

colnames(cnt_exp) <- gsub('[(\\.bio\\.replic)|(\\tech\\.replic)]','',gsub('X(.+)_R','\\1',colnames(cnt_exp)))
colnames(cnt_exp) <- col_name_exp[colnames(cnt_exp)]

# processing rows:
cnt_exp2 <- cnt_exp[apply(cnt_exp,1,max)>0,]
cnt_ref2 <- cnt_ref_intersect[apply(cnt_ref_intersect,1,max)>0,]
# convert to ens:
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ref_gene <- getBM(attributes=c('ensembl_gene_id','entrezgene'),filters = 'entrezgene', values = rownames(cnt_ref2), mart = ensembl)
dup_gene <- ref_gene$ensembl_gene_id[duplicated(ref_gene$ensembl_gene_id)]
ref_gene <- ref_gene[!ref_gene$ensembl_gene_id %in% dup_gene,]
cnt_ref3 <- cnt_ref2[as.character(ref_gene$entrezgene),]
rownames(cnt_ref3) <- ref_gene$ensembl_gene_id

intersect_gene <- intersect(rownames(cnt_ref3),rownames(cnt_exp2))
cnt_exp_final <- cnt_exp2[intersect_gene,]
cnt_ref_final <- cnt_ref3[intersect_gene,]

set <- newSeqExpressionSet(as.matrix(cnt_exp2),phenoData = data.frame(colnames(cnt_exp2), row.names=colnames(cnt_exp2)))
spikes <- grep('^ERCC',rownames(cnt_exp2),value = T)
set1 <- RUVg(set, spikes, k=1)
cnt_exp_final_ruv <- normCounts(set1)[intersect_gene,]

colnames(cnt_ref_final) <- paste(colnames(cnt_ref_final),'ref',sep = '_')

cnt_final <- cbind(cnt_exp_final,cnt_ref_final)
cor_mat <- cor(cnt_final,method = 'spearman')
dissimilarity <- 1 - cor_mat
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - Correlation", xlab="")

cnt_final_ruv <- cbind(cnt_exp_final_ruv,cnt_ref_final)
cor_mat_ruv <- cor(cnt_final_ruv,method = 'spearman')
dissimilarity_ruv <- 1 - cor_mat_ruv
distance_ruv <- as.dist(dissimilarity_ruv)
plot(hclust(distance_ruv),main="Dissimilarity = 1 - Correlation, after RUV", xlab="")


