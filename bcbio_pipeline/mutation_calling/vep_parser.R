keep_types  = c('frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','protein_altering_variant','start_lost','stop_gained','stop_lost')
vep_folder  = "/home/cc400/vep_files/"

vep_parser <- function(vep_folder,vep_file){
  sample_id <- gsub('.vep','',vep_file,fixed = T)
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
  mut_table <- rbind(mut_table,vep_parser(vep_folder,vep_file))
}
write.table(mut_table,'~/Dropbox/_ChrisProject/workspace/RNAseq/20170115/processed_data/variation_results.tsv',sep = '\t',row.names = F)
