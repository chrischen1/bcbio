vcf_path = '/home/cc400/vcf_files/'
for(i in list.files(vcf_path)){
  var_res <- NULL
  vcf  <- read.table(paste(vcf_path,i,sep = ''), quote="\"", stringsAsFactors=FALSE)
  vcf2 <- vcf[vcf$V7 == 'PASS' & vcf$V10 != '.:.:.:.:.:.:.:.:.' & vcf$V11 != '.:.:.:.:.:.:.:.:.',]
  vcf3 <- vcf2[gsub('\\\\','|',gsub('(.{3}).+','\\1',vcf2$V10)) != gsub('\\\\','|',gsub('(.{3}).+','\\1',vcf2$V11)),]
  var_res <- rbind(var_res,vcf3)
  write.table(var_res,paste('~/vcf_filter/',i,sep = ''),quote = F,col.names = F,row.names = F)
}