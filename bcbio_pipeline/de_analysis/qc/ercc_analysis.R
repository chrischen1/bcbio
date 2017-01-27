source('~/Dropbox/_ChrisProject/workspace/bcbio/src/api_rnaseq.R')

cnt <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/processed_data/CDK46_response_201608.counts", row.names=1, stringsAsFactors=FALSE)
cnt_tgt <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/processed_data/CDK46_response_targeted.counts", row.names=1, stringsAsFactors=FALSE)
ERCC_annotation <- read.delim("~/Dropbox/_ChrisProject/workspace/bcbio/annotations/ercc/ERCC_annotation.tsv", stringsAsFactors=FALSE)
colnames(cnt_tgt) <- gsub('([A-Z])0','\\1',colnames(cnt_tgt))
probe_genes <- unlist(read.delim('~/Dropbox/_ChrisProject/workspace/bcbio/annotations/probe_genes.txt',as.is = T,header = F))
probe_spikes <- unlist(read.delim('~/Dropbox/_ChrisProject/workspace/bcbio/annotations/probe_spikes.txt',as.is = T,header = F))

cnt_ercc <- cnt[ERCC_annotation$ERCC.ID,]
cnt_tgt_ercc <- cnt_tgt[ERCC_annotation$ERCC.ID,]
design_count <- ERCC_annotation$concentration.in.Mix.2..attomoles.ul./ERCC_annotation$length
names(design_count) <- ERCC_annotation$ERCC.ID

cnt_ercc_i <- cnt_ercc[probe_spikes,colnames(cnt_tgt_ercc)]
cnt_tgt_ercc_i <- cnt_tgt_ercc[probe_spikes,]
ercc_ratio <- (cnt_tgt_ercc_i+0.1)/(cnt_ercc_i+0.1)
ercc_ratio <- ercc_ratio[,-1]

boxplot(t(ercc_ratio),las=2,cex.axis=0.7)
par(mfrow=c(1,2))
plot(density(apply(ercc_ratio,1,var)),main='Density for capture efficiency in spikes')
spikes_var_within <- apply(ercc_ratio,1,var)
plot(log(design_count[names(spikes_var_within)]+0.1,2),log(spikes_var_within,2),
     xlab='log2 Concentration in design',ylab= 'log2 Variation for spikes',main='Concentration vs Efficiency Variation')
abline(v=-3,col='red')


spikes_filtered <- names(design_count)[design_count>0.125]
write.csv(spikes_filtered,'spikes_filtered.csv',row.names = F,col.names = F,quote = F)

#filter genes




# par(mfrow=c(1,2))
# plot(log(design_count),log(apply(cnt_ercc,1,sum)+0.1,2),
#      main = 'RNAseq',xlab = 'ERCC conctration/lenth,log2',ylab = 'counts,log2',ylim=c(-1,20))
plot(log(design_count),log(apply(cnt_tgt_ercc,1,sum)+0.1,2),
     main = 'CaptureSeq',xlab = 'ERCC conctration/length,log2',ylab = 'counts,log2',ylim=c(-1,20),cex=0.2)
text(log(design_count),log(apply(cnt_tgt_ercc,1,sum)+0.1,2),gsub('ERCC-','',rownames(cnt_ercc)),cex=0.8)

cnt2 <- cnt[apply(cnt,1,max)>0,]
cnt_tgt2 <- cnt_tgt[apply(cnt_tgt,1,max)>0,]
gene_intersect2 <- intersect(rownames(cnt2),rownames(cnt_tgt2))


gene_cap_all <- target_gene_processing()
gene_cap <- intersect(gene_cap_all[,1],gene_intersect2)
gene_nocap <- gene_intersect2[!gene_intersect2 %in% gene_cap]

gene_cap_test <- gene_cap
gene_nocap_test <- sample(gene_nocap,length(gene_cap))

logc_gene_cap_raw <- log(apply(cnt2[gene_cap_test,],1,sum)+0.1,2)
logc_gene_cap_tgt <- log(apply(cnt_tgt2[gene_cap_test,],1,sum)+0.1,2)
# before filtering 
logc_gene_nocap_raw <- log(apply(cnt2[gene_nocap_test,],1,sum)+0.1,2)
logc_gene_nocap_tgt <- log(apply(cnt_tgt2[gene_nocap_test,],1,sum)+0.1,2)
plot(logc_gene_cap_raw,logc_gene_cap_tgt,
     main = 'Log2 Counts for genes',xlab = 'RNAseq',ylab = 'Captureseq',cex=0.4,xlim=c(-1,22),ylim=c(-1,22))
points(logc_gene_nocap_raw,logc_gene_nocap_tgt,col=2,cex=0.4)
legend('topleft',c('Captured','Not captured'),pch='o',col=c(1,2))



test_set <- sample(colnames(cnt_tgt2),24)
train_set <- colnames(cnt_tgt2)[!colnames(cnt_tgt2) %in% test_set]

logc_gene_cap_raw_train <- log(apply(cnt2[gene_cap_test,train_set],1,sum)+0.1,2)
logc_gene_cap_tgt_train <- log(apply(cnt_tgt2[gene_cap_test,train_set],1,sum)+0.1,2)

gene_dif <- logc_gene_cap_tgt_train - logc_gene_cap_raw_train[names(logc_gene_cap_tgt_train)]

logc_gene_cap_raw_test <- log(apply(cnt2[gene_cap_test,test_set],1,sum)+0.1,2)
logc_gene_cap_tgt_test <- log(apply(cnt_tgt2[gene_cap_test,test_set],1,sum)+0.1,2)

plot(logc_gene_cap_raw_test,logc_gene_cap_tgt_test-gene_dif[names(logc_gene_cap_tgt_test)],
     main = 'Log Counts for spikes',xlab = 'RNAseq',ylab = 'Captureseq',cex=0.4)
points(logc_gene_nocap_raw,logc_gene_nocap_tgt,col=2,cex=0.4)


# after filtering 
f_cap_raw <- names(logc_gene_cap_raw)[logc_gene_cap_raw>12]
f_cap_tgt <- names(logc_gene_cap_tgt)[logc_gene_cap_tgt>12]
logc_gene_cap_rawf <- logc_gene_cap_raw[unique(f_cap_raw,f_cap_tgt)]
logc_gene_cap_tgtf <- logc_gene_cap_tgt[unique(f_cap_raw,f_cap_tgt)]

f_nocap_raw <- names(logc_gene_nocap_raw)[logc_gene_nocap_raw>12]
f_nocap_tgt <- names(logc_gene_nocap_tgt)[logc_gene_nocap_tgt>12]
logc_gene_nocap_rawf <- logc_gene_nocap_raw[unique(f_nocap_raw,f_nocap_tgt)]
logc_gene_nocap_tgtf <- logc_gene_nocap_tgt[unique(f_nocap_raw,f_nocap_tgt)]

plot(logc_gene_cap_rawf,logc_gene_cap_tgtf,
     main = 'Log Counts for genes',xlab = 'RNAseq',ylab = 'Captureseq',cex=0.4)
points(logc_gene_nocap_rawf,logc_gene_nocap_tgtf,col=2,cex=0.4)

cap_model <- data.frame('x'=logc_gene_cap_rawf,'y'=logc_gene_cap_tgtf)
lm(y~x,data = cap_model)

nocap_model <- data.frame('x'=logc_gene_nocap_rawf,'y'=logc_gene_nocap_tgtf)
lm(y~x,data = nocap_model)


logc_ercc_raw <- log(apply(cnt_ercc,1,sum)+0.1,2)
logc_ercc_tgt <- log(apply(cnt_tgt_ercc,1,sum)+0.1,2)

ercc_noprobe <- rownames(cnt_ercc)[!rownames(cnt_ercc) %in% probe_spikes]

logc_ercc_type1_raw <- log(apply(cnt_ercc[probe_spikes,],1,sum)+0.1,2)
logc_ercc_type1_tgt <- log(apply(cnt_tgt_ercc[probe_spikes,],1,sum)+0.1,2)

logc_ercc_type2_raw <- log(apply(cnt_ercc[ercc_noprobe,],1,sum)+0.1,2)
logc_ercc_type2_tgt <- log(apply(cnt_tgt_ercc[ercc_noprobe,],1,sum)+0.1,2)
# ERCC probe vs no probe
plot(logc_ercc_type1_raw,logc_ercc_type1_tgt,main = 'Log Counts for spikes',
     xlab = 'RNAseq',ylab = 'Captureseq',cex=0.4,xlim=c(-1,22),ylim=c(-1,22))
points(logc_ercc_type2_raw,logc_ercc_type2_tgt,col=2,cex=0.4)
legend('topleft',c('Captured','Not captured'),pch='o',col=c(1,2))

cap_model_ercc <- data.frame('x'=logc_ercc_type1_raw,'y'=logc_ercc_type1_tgt)
lm(y~x,data = cap_model_ercc)
nocap_model_ercc <- data.frame('x'=logc_ercc_type2_raw,'y'=logc_ercc_type2_tgt)
lm(y~x,data = nocap_model_ercc)



# ERCC detected in RNAseq but not in CaptureSeq
ercc_cnt2 <- grep('^ERCC',rownames(cnt2),value = T)
ercc_cnt_tgt2 <- grep('^ERCC',rownames(cnt_tgt2),value = T)

ercc_cnt2[!ercc_cnt2 %in%ercc_cnt_tgt2]
