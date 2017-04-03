args <- commandArgs(trailingOnly = TRUE)

vep_path = '/groups/sorger/cchris/vep_basal_gatk/'


vep <- read.table(paste(vep_path,args,sep = ''), quote="\"", stringsAsFactors=FALSE)
vep2 <- vep[vep$V7 == 'missense_variant',]
id <- gsub('\\.vep.*','',i)
id <- gsub('_R1_001-freebayes','',id)
id <- gsub('_R-freebayes','',id)
mut_list[[id]] <- unique(paste(vep2$V5,vep2$V10,vep2$V11,sep = '_'))