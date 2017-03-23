vcf_path = '~/vcf_filter/'

vcf_all <- NULL
for(i in list.files(vcf_path)){
  vcf_all <- rbind(vcf_all,cbind(gsub('.vcf','',i,fixed = T),read.table(paste(vcf_path,i,sep = ''),as.is = T,header = F)))
}
count_info <- NULL
for(i in 1:nrow(vcf_all)){
  trt_x <- unlist(strsplit(vcf_all$V10[i],':'))
  ctr_x <- unlist(strsplit(vcf_all$V11[i],':'))
  names(trt_x) <- names(ctr_x) <- unlist(strsplit(vcf_all$V9[i],':'))
  count_info <- rbind(count_info,c(paste(vcf_all[i,2],vcf_all[i,3],vcf_all[i,6],vcf_all[i,7],sep = ':'),as.character(vcf_all[i,1]),vcf_all[i,5],vcf_all[i,6],
                                   trt_x['RO'],trt_x['AO'],ctr_x['RO'],ctr_x['AO']))
}
colnames(count_info) <- c('location','sample','ref','alt','trt_count_ref','trt_count_alt','ctr_count_ref','ctr_count_alt')

write.table(count_info,'~/Dropbox/_ChrisProject/workspace/RNAseq/20170115/processed_data/vcf_combined_pass.tsv',row.names = F)