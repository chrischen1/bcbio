# Clustering on parental lines in batch2 and basal cases

mut_path = '/groups/sorger/cchris/vep_files/'
x=list.files(mut_path)
mut_list <- list()

for(i in x){
  vep <- read.table(paste(mut_path,i,sep = ''), quote="\"", stringsAsFactors=FALSE)
  vep2 <- vep[vep$V7 == 'missense_variant',]
  id <- gsub('\\.vep.*','',i)
  id <- gsub('_R1_001-freebayes','',id)
  id <- gsub('_R-freebayes','',id)
  mut_list[[id]] <- unique(paste(vep2$V5,vep2$V10,vep2$V11,sep = '_'))
}
save(list = c('mut_list'),file = '~/mut_list.rda')

# for(i in x){
#   vcf <- read.table(paste(mut_path,i,sep = ''), quote="\"", stringsAsFactors=FALSE)
#   vcf2 <- vcf[vcf$V7 == 'PASS',]
#   mut_list[[gsub('.csv','',i,fixed = T)]] <- unique(paste(vcf2$V1,vcf2$V2,vcf2$V4,vcf2$V5,sep = '_'))
# }
# mut_path = '/groups/sorger/cchris/vcf_gatk_final/'

# for(i in list.files(mut_path)){
#   vcf <- read.table(paste(mut_path,i,sep = ''), quote="\"", stringsAsFactors=FALSE)
#   vcf2 <- vcf[vcf$V7 == 'PASS',]
#   write.csv(unique(paste(vcf2$V1,vcf2$V2,vcf2$V4,vcf2$V5,sep = '_')),paste('~/test_out/',gsub('vcf','csv',i,fixed = T),sep = ''))
# }

# mut_list <- list()
# mut_path <- '~/test_out/'
# for(i in list.files(mut_path)){
#   mut_list[[gsub('.csv','',i,fixed = T)]] <- read.delim(paste(mut_path,i,sep = ''),sep = ',',as.is = T)$x
# }

m <- matrix(0,length(names(mut_list)),length(names(mut_list)))
rownames(m) <- colnames(m) <- names(mut_list)
n <- names(mut_list)
for(a in 1:length(n)){
  for(b in a:length(n)){
    i=n[a]
    j=n[b]
    if(i==j){
      m[i,j] <- 1
    }else{
      m[i,j] <- m[j,i]<- length(intersect(mut_list[[i]],mut_list[[j]]))/length(unique(c(mut_list[[i]],mut_list[[j]])))
    }
  }
}

test <- m
rownames(test)<- colnames(test) <- gsub('-tech-replic','',rownames(test))
rownames(test)<- colnames(test) <- gsub('-bio-replic','',rownames(test))
rownames(test)<- colnames(test) <- gsub('_R$','',rownames(test))
rownames(test)<- colnames(test) <- gsub('20170115_(\\d+)_.+','\\1',rownames(test))
write.csv(test,'~/test.csv')

# 
basal_info <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/RNAseq_BrCaPorfiling_201604/annotation/RNAseq_sample_batch1.tsv", stringsAsFactors=FALSE)
cdk_info <- read.delim("~/Dropbox/_ChrisProject/workspace/RNAseq/20170115/annotation/20161220RNAseqCDKiPR.tsv",as.is = T)

test <- read.csv("~/Dropbox/_ChrisProject/workspace/RNAseq/scratch/test.csv", row.names=1, stringsAsFactors=FALSE)
cells <- c(gsub(' ','_',c(paste(basal_info[,2],basal_info[,3]),cdk_info[25:36,2])),'mcf7_encode')
names(cells) <- c(basal_info[,1],cdk_info[25:36,1],'mcf7_encode')
cells2 <- paste(cells[1:nrow(basal_info)],'GATK',sep = '')
names(cells2) <- paste(names(cells[1:nrow(basal_info)]),'qc',sep = '_')
cells <- c(cells,cells2)
rownames(test)<- colnames(test) <- cells[rownames(test)]
test <- test[-7,-7]

dissimilarity <- 1 - test
distance <- as.dist(dissimilarity)
plot(hclust(distance),main="Dissimilarity = 1 - mutation overlapping", xlab="",cex=0.7)
