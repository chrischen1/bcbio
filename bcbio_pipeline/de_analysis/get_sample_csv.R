sample_path = '/n/scratch2/cc400/rnaseq_bcbio/rnaseq_batch2/'
x=grep('\\.fastq',list.files(sample_path),value = T)
y=gsub('\\.fastq','',x)
z=cbind('samplename'=y,'description'=y)
write.csv(z,paste(sample_path,'samples.csv',sep = ''),row.names = F,quote = F)