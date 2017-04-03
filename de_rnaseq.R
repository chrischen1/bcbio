## de_rnaseq.R - functions for the differential expression analysis
##
## LSP RNAseq bcbio pipeline 
## https://github.com/chrischen1/bcbio

library(edgeR)
library(biomaRt)

#' transform TPM to RPKM
#'
#' @param combined output file end with .combined from bcbio.
#' @param tx2gene output file which maps ensumble ID to gene from bcbio.
#' @param spikes a vector of string defining the name of spikes.
#' @return p by n matrix for p genes across n samples
tpm2rpkm <- function(combined,tx2gene,spikes = NULL){
  gene_mapping <- cbind('transcript'= c(tx2gene$V1,spikes$GenBank),'gene' = c(tx2gene$V2,spikes$ERCC_ID))
  genes <- gene_mapping[,2]
  names(genes) <- gene_mapping[,1]
  lib_size <- data.frame('numreads'=combined$numreads,'sample'=combined$sample)
  x <- lib_size %>% group_by(sample) %>% summarise_each(funs(sum))
  scale_factor <- x$numreads/1000000
  names(scale_factor) <- x$sample
  
  combined$RPM <- combined$numreads/scale_factor[combined$sample]
  combined$RPKM <- combined$RPM/(combined$effectiveLength/1000)
  combined$gene <- genes[combined$id]
  
  rpkm_combined <- data.frame('sample'=combined$sample,'gene'=combined$gene,'RPKM'=combined$RPKM)
  rpkm_combined_gene <- rpkm_combined %>% group_by(sample,gene)%>% summarise_each(funs(sum))
  
  rpkm_raw <- acast(rpkm_combined_gene,gene~sample)
}

#' get hgnc_symbol from ensembl_gene_id from 
#'
#' @param ens vector of ensembl_gene_ids.
#' @return a dataframe with 2 columns: ensembl_gene_id and hgnc_symbol
ens2symbol <- function(ens){
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  target_gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = ens, mart = ensembl)
  return(target_gene)
}

#' generate .csv file used by bcbio
#' 
#' @param sample_path path of .fastq files
#' @return a csv file contain basic sample meta info required for bcbio
get_sample_csv <- function(sample_path){
  x=grep('\\.fastq',list.files(sample_path),value = T)
  y=gsub('\\.fastq','',x)
  z=cbind('samplename'=y,'description'=y)
  write.csv(z,paste(sample_path,'samples.csv',sep = ''),row.names = F,quote = F)
}

#' wrapper for getting fold change, pvalue and FDR, by per cell line per time point
#' 
#' @param cnt p by n matrix for p genes across n samples
#' @param grp_table dataframe with 3 columns: group, condition and control
#'  group: contains information which treatment samples will be compared against control cases in each group
#'  condition: indicates type of treatment, replicates have same condition
#'  control: TRUE for controls and FALSE for treatments
#'  order of well in samples annotation must be the same as the columns in count table
#'  @param combine_fdr T for combine FDR and p-values with group and F for compute pairwisely
#'  @param w n by p matrix for n samples and p factors for batch effect correction from RUVSeq
#'  @param CommonDisp and TagwiseDisp used internally for passing overal dispersion to comparisons without replicates
#'  @return list of 3 if combine_fdr = F: pmat,fdr_mat and logFC: all are p by m matrix for p genes across m types of treatments
#'          p by m+4 matrix for p genes across m types of treatments and p-value, LR,logCPM and FDR
edgeR_wrapper <- function(cnt,grp_table,combine_fdr = F,w = NULL,CommonDisp = NULL,TagwiseDisp = NULL){
  design <- model.matrix(~condition,data = grp_table)
  if(!is.null(w))  design <- cbind(design,w)
  y <- DGEList(counts=cnt, group=grp_table$condition)
  if(is.null(CommonDisp)){
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    CommonDisp <- y$common.dispersion
    TagwiseDisp <- y$tagwise.dispersion
  }
  if(length(grp_table$condition)==unique(length(grp_table$condition))){
  # When both control and treatment lacking replicates, use overall dispersion instead
    y$common.dispersion <- CommonDisp
    y$tagwise.dispersion <- TagwiseDisp
  }else{
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
  }
  if(combine_fdr){
    y <- calcNormFactors(y)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2:(ncol(design)))
    lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
    colnames(lrt_tab) <- gsub('logFC.condition','',colnames(lrt_tab))
    return(lrt_tab)
  }
  p_mat <- fdr_mat <- logFC <- NULL
  col_names <- c()
  for(i in unique(grp_table$group)){
    grp_table_i <- grp_table[grp_table$group==i,]
    ctr_row <- rownames(grp_table_i)[grp_table_i$control==T]
    for (j in unique(grp_table_i$condition[grp_table_i$control!=T])){
      j_row <- rownames(grp_table_i)[grp_table_i$condition==j]
      grp_new <- rbind(grp_table[c(ctr_row,j_row),])
      cnt_new <- cnt[,rownames(grp_new)]
      result_new <- edgeR_wrapper(cnt_new,grp_new,combine_fdr = T,CommonDisp = CommonDisp,TagwiseDisp = TagwiseDisp)
      if(is.null(p_mat)){
        p_mat <- result_new$PValue
        fdr_mat <- result_new$FDR
        logFC <- result_new$logFC
      }else{
        p_mat <- cbind(p_mat,result_new$PValue)
        fdr_mat <- cbind(fdr_mat,result_new$FDR)
        logFC <- cbind(logFC,result_new$logFC)
      }
      col_names <- c(col_names,j)
    }
  }
  colnames(p_mat) <- colnames(fdr_mat) <- colnames(logFC) <- col_names
  rownames(p_mat) <- rownames(fdr_mat) <- rownames(logFC) <- rownames(cnt)
  return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
}