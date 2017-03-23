vcfDir  = '/home/cc400/vcf_filter/'
outDir    = '/home/cc400/vep_files/'
bsub_fig  = 'bsub -q short -W 4:00 -o '
vep       = 'perl ~/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl'

dir.create(outDir,showWarnings = F)
vcf_files <- grep('.vcf$',grep('vep',list.files(vcfDir),value = T,invert = T),value = T)

for (i in vcf_files){
  sample_id = gsub('.vcf','',i,fixed = T)
  system(paste(bsub_fig,outDir,sample_id,'.out ',vep,' -i ',vcfDir,i,' --offline -symbol -o ',outDir,sample_id,'.vep',sep = ''))
}