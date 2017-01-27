#Current Pipeline: bcbio-nextgen
http://bcbio-nextgen.readthedocs.io/en/latest/contents/introduction.html

##RPKM/TPM and Differential Expression analysis 

The setup of pipeline is determined by the yaml file. The template yaml file used in this project is based on:
https://github.com/chrischen1/bcbio/blob/master/bcbio_pipeline/de_analysis/SE_PE_align.yaml

Generate yaml file based on template for samples:
bcbio_nextgen.py -w template ./ SE_PE_align.yaml project.csv ./*

cd project/work

Start alignment:
bsub -q long -W 100:0 -R "rusage[mem=10000]" -o output.out bcbio_nextgen.py ../config/project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380

As bcbio doesn’t give RPKM directly(only feature counts and TPM), I have calculated RPKM based on the count of each transcript:

1 Count up the total reads of transcripts in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.

2 Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving reads per million (RPM)

3 Divide the RPM values by the length of the “effective length” given by bcbio, in kilobases. This gives RPKM.

4 Bcbio have a list of all transcripts ID and their corresponding genes ID, sum all RPKMs of transcripts of a gene to get the RPKM of this gene. (Bcbio use this method to get TPM) 

The only difference changed in the format is that I added “Cellline” as the first column name for the raw RPKM before merging(the same format for the file after merging) in order to show correct column names in synapse preview, this should not change the data format in R if both. Also, the result with merged replicates have been updated according to new RPKM results, so future prediction models could be built upon them.

The transformation from TPM(bcbio default output) to PRKM is based on this script:
https://github.com/chrischen1/bcbio/blob/master/bcbio_pipeline/de_analysis/tpm2rpkm.R

##Variant Calling

The setup of pipeline is determined by the yaml file. The template yaml file used in this project is based on:
https://github.com/chrischen1/bcbio/blob/master/bcbio_pipeline/mutation_calling/mutation_calling.yaml
The pipeline uses freebayes approach to call variants in samples:
http://clavius.bc.edu/~erik/CSHL-advanced-sequencing/freebayes-tutorial.html

Generate yaml file based on template for samples:
bcbio_nextgen.py -w template ./mutation_calling.yaml project.csv ./*

cd project/work

Start alignment:
bsub -q long -W 100:0 -R "rusage[mem=10000]" -o output.out bcbio_nextgen.py ../config/project.yaml -n 32 -t ipython -s lsf -q parallel -r mincores=2 -r minconcores=2 '-rW=72:00' --retries 3 --timeout 380

The output of bcbio are .vcf files, use the script below to transform them into readable .vep files by Variant Effect Predictor v85 (http://useast.ensembl.org/info/docs/tools/vep/index.html, required installation)
https://github.com/chrischen1/bcbio/blob/master/bcbio_pipeline/mutation_calling/vcf2vep.R

To keep only certian variant types, use this script below after the last step(you can specify which type to include in the script):
https://github.com/chrischen1/bcbio/blob/master/bcbio_pipeline/mutation_calling/vep_parser.R
