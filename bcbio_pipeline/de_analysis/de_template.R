source('~/Dropbox/_ChrisProject/workspace/bcbio/src/api_rnaseq.R')
count_file = "~/Dropbox/_ChrisProject/workspace/RNAseq/processed_data/batch2/combined.counts"
mata_file  = "~/Dropbox/_ChrisProject/workspace/RNAseq/processed_data/batch2/20161220RNAseqCDKiPR.tsv"
out_folder = '~/Dropbox/_ChrisProject/workspace/RNAseq/processed_data/batch2/'

# data cleaning
counts <- read.delim(count_file,as.is = T,row.names = 1)
colnames(counts) <- gsub('.+_S(\\d+)_R.+','\\1',colnames(counts))

# read meta data and change it to starndard format
sample_annotations <- read.delim(mata_file, stringsAsFactors=FALSE)
sample_annotations <- sample_annotations[sample_annotations$Treatment != 'N/A',]
trt <- strsplit(sample_annotations$Treatment,split = ' ')
samples <- data.frame(cbind('well'=sample_annotations$Sample,'CellLine'=sample_annotations$Cell.line,'DrugName'=unlist(lapply(trt,function(x)x[2])),
                 'Conc'=unlist(lapply(trt,function(x)x[1])),'Time'=gsub(' hours','',sample_annotations$Time.point)),stringsAsFactors = F)
samples$DrugName[is.na(samples$DrugName)] <- '-'
samples$Conc[samples$Conc=='ctrl'] <- '0.0'
samples$Conc <- as.numeric(samples$Conc)
samples$ctrl <- samples$Conc==0
samples$Time[samples$ctrl] <- 0

# we want to use control in any time for both 6h and 24h group, so we calculate them separately in different group table
#24h
samples_24 <- samples[samples$Time==24 | samples$Time==0,]
samples_24$Time[samples_24$Time==0] <-24
grp_table_24 <- data.frame(cbind('group' = paste(samples_24$CellLine,samples_24$Time,sep = '_'),
                                 'condition' = paste(samples_24$CellLine,samples_24$DrugName,samples_24$Conc,samples_24$Time,sep = '_'),
                                 'control' = samples_24$ctrl),stringsAsFactors = F)
rownames(grp_table_24) <- samples_24$well

cnt_24 <- counts[,samples_24$well]
raw_result_24 <- edgeR_wrapper2(cnt_24,grp_table_24)

#same for 6h
samples_6 <- samples[samples$Time==6 | samples$Time==0,]
samples_6$Time[samples_6$Time==0] <-6
grp_table_6 <- data.frame(cbind('group' = paste(samples_6$CellLine,samples_6$Time,sep = '_'),
                                'condition' = paste(samples_6$CellLine,samples_6$DrugName,samples_6$Conc,samples_6$Time,sep = '_'),
                                'control' = samples_6$ctrl),stringsAsFactors = F)
rownames(grp_table_6) <- samples_6$well

cnt_6 <- counts[,samples_6$well]
raw_result_6 <- edgeR_wrapper2(cnt_6,grp_table_6)

logFC_raw <- cbind(raw_result_6$logFC,raw_result_24$logFC[rownames(raw_result_6$logFC),])
pval_raw <- cbind(raw_result_6$pmat,raw_result_24$pmat[rownames(raw_result_6$logFC),])
fdr_raw <- cbind(raw_result_6$fdr_mat,raw_result_24$fdr_mat[rownames(raw_result_6$logFC),])

cnt_sybl <- ens2symbol(rownames(logFC_raw))
cnt_sybl <- cnt_sybl[cnt_sybl$ensembl_gene_id!=''&cnt_sybl$hgnc_symbol!='',]
genes <- cnt_sybl$hgnc_symbol
names(genes) <- cnt_sybl$ensembl_gene_id

logFC <- logFC_raw[rownames(logFC_raw)%in%names(genes),]
rownames(logFC) <- genes[rownames(logFC)]

pval <- pval_raw[rownames(pval_raw)%in%names(genes),]
rownames(pval) <- genes[rownames(pval)]

fdr <- fdr_raw[rownames(fdr_raw)%in%names(genes),]
rownames(fdr) <- genes[rownames(fdr)]

write.csv(logFC,paste(out_folder,'logFC.csv',sep = ''))
write.csv(pval,paste(out_folder,'pval.csv',sep = ''))
write.csv(fdr,paste(out_folder,'fdr.csv',sep = ''))

