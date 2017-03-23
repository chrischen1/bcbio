sample_path = '/n/scratch2/cc400/rnaseq_bcbio/rnaseq_encode_mcf7/'
x=grep('._1\\.fastq',list.files(sample_path),value = T)
y=gsub('(.*)_1\\.fastq','\\1',x)
z=cbind('samplename'=y,'description'=y)
write.csv(z,paste(sample_path,'samples.csv',sep = ''),row.names = F,quote = F)