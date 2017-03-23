keep_types  = c('frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','protein_altering_variant','start_lost','stop_gained','stop_lost')
vep_folder  = "~/Dropbox/_ChrisProject/workspace/RNAseq/20170115/processed_data/paired_calling/"
x <- read.delim('/n/scratch2/cc400/rnaseq_bcbio/rnaseq_batch2/samples/final/2017-01-21_samples/combined.sf',as.is = T)

vep_parser <- function(vep_folder,vep_file,samples){
  sample_id <- gsub('-freebayes.vep','',vep_file)
  vep <- read.table(paste(vep_folder,vep_file,sep = ''), quote="\"",stringsAsFactors = F)
  vep_f <- vep[vep$V7 %in% keep_types,]
  genes <- gsub('.+SYMBOL=(.+);SYMBOL_SOURCE.+','\\1',vep_f$V14)
  mut_table <- cbind(sample_id,genes,vep_f$V4,vep_f$V5,vep_f$V7,vep_f$V2,vep_f$V3,vep_f$V8,vep_f$V9,vep_f$V10,vep_f$V11)
  colnames(mut_table) <- c('SAMPLE','SYMBOL','ENS_GENE','ENS_FEATURE','CONSEQUENCE','LOCATION','ALLELE','cDNA_POS','CDS_POS','PROT_POS','AMINO_ACID_CHANGE')
  return(mut_table)
}

vep_files <- grep('vep$',list.files(vep_folder),value = T)
mut_table <- NULL
for(vep_file in vep_files){
  mut_table <- rbind(mut_table,vep_parser(vep_folder,vep_file,samples))
}

num_read <- x$numreads
names(num_read) <- paste(gsub('.+_S(\\d+)_R1.+','\\1',x$sample),x$id,sep = '_')

pair_table <- rbind(c('pair1_hs578t',25,26),c('pair2_hs578t',25,27),c('pair1_mcf7',28,29),
                    c('pair1_t47d',31,32),c('pair2_t47d',31,33),c('pair1_hcc1806',34,35),c('pair2_hcc1806',34,36))
id_ctr <- pair_table[,2]
id_trt <- pair_table[,3]
names(id_ctr) <- names(id_trt) <- pair_table[,1]
mut_table2 <- cbind(mut_table,'read_ctr'=num_read[paste(id_ctr[mut_table[,1]],mut_table[,4],sep = '_')],'read_trt'=num_read[paste(id_trt[mut_table[,1]],mut_table[,4],sep = '_')])
write.table(mut_table2,'~/variation_results2.tsv',sep = '\t',row.names = F)

mut_table3 <- mut_table2[mut_table2$read_ctr>4 & mut_table2$read_trt>4,c(1,2,4,5,10:13)]
mut <- unique(cbind(mut_table3$SAMPLE,mut_table3$SYMBOL))
a=table(mut[,2])
a[a>1]
mut_table3[mut_table3$SYMBOL=='AHNAK2',]
mut_table3[mut_table3$SYMBOL=='DDX24',]
mut_table3[mut_table3$SYMBOL=='FIP1L1',]
mut_table3[mut_table3$SYMBOL=='MMADHC',]
mut_table3[mut_table3$SYMBOL=='PRR14L',]
mut_table3[mut_table3$SYMBOL=='PSMD1',]
mut_table3[mut_table3$SYMBOL=='SLC29A3',]
mut_table3[mut_table3$SYMBOL=='SNRNP200',]
mut_table3[mut_table3$SYMBOL=='TOP2A',]
mut_table3[mut_table3$SYMBOL=='ZNF594',]

# mut_table2 <- data.frame(mut_table,stringsAsFactors = F)
# mut_table3 <- unique(cbind(mut_table2$SAMPLE,mut_table2$SYMBOL,paste(mut_table2$SYMBOL,mut_table2$PROT_POS,mut_table2$AMINO_ACID_CHANGE,sep = '_'),mut_table2$CONSEQUENCE))
# mut_table3 <- data.frame(mut_table3,stringsAsFactors = F)
# 
# table_pair <- table(mut_table3[,1])
# table_gene <- table(mut_table3[,2])
# table_var <- table(mut_table3[,3])

