sample_path = '/n/scratch2/cc400/rnaseq_bcbio/rnaseq_batch1/'
x=grep('\\.R1',list.files(sample_path),value = T)
y=gsub('(.*\\.R)1.+','\\1',x)
z=cbind('samplename'=y,'description'=y)
write.csv(z,paste(sample_path,'samples.csv',sep = ''),row.names = F,quote = F)