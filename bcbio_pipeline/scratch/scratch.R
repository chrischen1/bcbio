#' wrapper for getting fold change with RUV normalization implemented, pvalue and FDR, by per cell line per time point
#' 
#' @param cnt_time p by n matrix for p genes across n samples
#' @param grp_table_time dataframe with 3 columns: group, condition and control
#'  group: contains information which treatment samples will be compared against control cases in each group
#'  condition: indicates type of treatment, replicates have same condition
#'  control: TRUE for controls and FALSE for treatments
#'  order of well in samples annotation must be the same as the columns in count table
#'  @param combine_fdr T for combine FDR and p-values with group and F for compute pairwisely
#'  @return list of 3 if combine_fdr = F: pmat,fdr_mat and logFC: all are p by m matrix for p genes across m types of treatments
#'          p by m+4 matrix for p genes across m types of treatments and p-value, LR,logCPM and FDR
# edgeR_wrapper_ruv <- function(cnt,grp_table,w=NULL){
#   if(is.null(w) ){
#     return(edgeR_wrapper(cnt,grp_table))
#   }else{
#     logFC <- NULL
#     p_mat <- NULL
#     fdr_mat <- NULL
#     mat_label <- c()
#     grp_table <- grp_table[colnames(cnt),]
#     
#   }
#   
#  
#   if((!is.null(w)) & (!is.null(com_disp)) & (!is.null(tag_disp))){
#     design <- cbind(model.matrix(~condition,data = grp_table),w)
#     y <- DGEList(counts=cnt, group=grp_table$condition)
#     y <- estimateGLMCommonDisp(y, design)
#     y <- estimateGLMTagwiseDisp(y, design)
#     com_disp <- y$common.dispersion
#     tag_disp <- y$tagwise.dispersion
#   }
#   for (i in unique(grp_table$group)){
#     selected_grp_table <- grp_table[grp_table$group==i,]
#     selected_grp_table <- rbind(selected_grp_table[as.logical(selected_grp_table$control),],selected_grp_table[!as.logical(selected_grp_table$control),])
#     selected_grp <- selected_grp_table$condition
#     selected_cnt <- cnt[,rownames(selected_grp_table)]
#     selected_w <- w[rownames(selected_grp_table),]
#     y_selected <- DGEList(counts=selected_cnt, group=selected_grp)
#     if(is.null(w)){
#       design <- model.matrix(~condition,data = selected_grp_table)
#       if((!is.null(com_disp)) & (!is.null(tag_disp))){
#         y_selected$common.dispersion <- com_disp
#         y_selected$tagwise.dispersion <- tag_disp
#       }else{
#         y_selected <- estimateGLMCommonDisp(y_selected, design)
#         y_selected <- estimateGLMTagwiseDisp(y_selected, design)
#       }
#       y_selected <- calcNormFactors(y_selected)
#       fit <- glmFit(y_selected, design)
#       lrt <- glmLRT(fit, coef=2:(ncol(design)))
#     }else {
#       design <- cbind(model.matrix(~condition,data = selected_grp_table),selected_w)
#       y_selected$common.dispersion <- com_disp
#       y_selected$tagwise.dispersion <- tag_disp
#       y_selected <- calcNormFactors(y_selected)
#       fit <- glmFit(y_selected, design)
#       lrt <- glmLRT(fit, coef=2:(ncol(design)-ncol(w)))
#     }
#     lrt_tab <- topTags(lrt,n = Inf)$table[rownames(cnt),]
#     if(is.null(logFC)){
#       logFC <- lrt_tab[,grep('logFC',colnames(lrt_tab))]
#       p_mat <- lrt_tab$PValue
#       fdr_mat <- lrt_tab$FDR
#     }else{
#       logFC <- cbind(logFC,lrt_tab[,grep('logFC',colnames(lrt_tab))])
#       p_mat <- cbind(p_mat,lrt_tab$PValue)
#       fdr_mat <- cbind(fdr_mat,lrt_tab$FDR)
#     }
#     mat_label <- c(mat_label,i)
#   }
#   if(!is.null(dim(p_mat))){
#     colnames(p_mat) <- colnames(fdr_mat) <- mat_label
#     rownames(p_mat) <- rownames(fdr_mat) <- rownames(cnt)
#     colnames(logFC) <- gsub('logFC.condition','',colnames(logFC))
#   }else{
#     names(logFC) <- names(p_mat) <- names(fdr_mat) <- rownames(cnt)
#   }
#   return(list('pmat'=p_mat,'fdr_mat'=fdr_mat,'logFC'=logFC))
# }

# get_fc_exact <- function(cnt,grp,x='logFC'){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   design <- model.matrix(~grp)
#   DGEx1f <- estimateDisp(DGEx1f,design)
#   logFC <- NULL
#   col_label <- c()
#   for(i in unique(samples$CellLine)){
#     grp_cell <- grep(i,grp,value = T)
#     grp_cell_ctr <- intersect(grep('-',grp,value = T),grp_cell)
#     grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
#     for(j in grp_cell_trt){
#       col_label <- c(col_label,j)
#       fit <-exactTest(DGEx1f, pair=c(grp_cell_ctr,j)) 
#       tags <- topTags(fit, n=Inf)$table
#       val <- which(x==colnames(tags))
#       tags2 <- tags[,val]
#       names(tags2) <- tags$genes
#       tags2 <- tags2[DGEx1$genes$genes]
#       if(is.null(logFC)) logFC <- tags2
#       else logFC <- cbind(logFC,tags2)
#     }
#   }
#   colnames(logFC) <- col_label
#   return(logFC)
# }
# 
# get_fc_glm <- function(cnt,grp){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   logFC <- NULL
#   for(i in unique(samples$CellLine)){
#     grp_cell <- grep(i,grp)
#     grp_cell_ctr <- intersect(grep('-',grp),grp_cell)
#     grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
#     grp_cell_all <- c(grp_cell_ctr,grp_cell_trt)
#     y_d = DGEx1f[,grp_cell_all]
#     grp_c <- y_d$samples$group
#     design <- model.matrix(~grp_c)
#     y_d <- estimateDisp(y_d, design, robust=TRUE)
#     fit <- glmFit(y_d, design)
#     lrt <- glmLRT(fit, coef=2:ncol(fit$design))
#     tags <- topTags(lrt, n=Inf)$table
#     tags2 <- tags[,2:ncol(fit$design)]
#     rownames(tags2) <- tags$genes
#     tags2 <- tags2[DGEx1$genes$genes,]
#     if(is.null(logFC)) logFC <- tags2
#     else logFC <- cbind(logFC,tags2)
#   }
#   colnames(logFC) <- gsub('logFC.grp_c','',colnames(logFC))
#   return(logFC)
# }
# 
# get_fc_qglm <- function(cnt,grp){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   logFC <- NULL
#   for(i in unique(samples$CellLine)){
#     grp_cell <- grep(i,grp)
#     grp_cell_ctr <- intersect(grep('-',grp),grp_cell)
#     grp_cell_trt <- grp_cell[grp_cell!=grp_cell_ctr]
#     grp_cell_all <- c(grp_cell_ctr,grp_cell_trt)
#     y_d = DGEx1f[,grp_cell_all]
#     grp_c <- y_d$samples$group
#     design <- model.matrix(~grp_c)
#     y_d <- estimateDisp(y_d, design, robust=TRUE)
#     fit <- glmQLFit(y_d, design)
#     lrt <- glmQLFTest(fit, coef=2:ncol(fit$design))
#     tags <- topTags(lrt, n=Inf)$table
#     tags2 <- tags[,2:ncol(fit$design)]
#     rownames(tags2) <- tags$genes
#     tags2 <- tags2[DGEx1$genes$genes,]
#     if(is.null(logFC)) logFC <- tags2
#     else logFC <- cbind(logFC,tags2)
#   }
#   colnames(logFC) <- gsub('logFC.grp_c','',colnames(logFC))
#   return(logFC)
# }
# 
# get_raw_fc <- function(cnt,samples){
#   cnt <- cnt+0.01
#   raw_fc <- NULL
#   for (i in unique(samples$CellLine)){
#     well_ctr <- samples$well[samples$CellLine==i & samples$ctrl]
#     well_trt <- samples$well[samples$CellLine==i & (!samples$ctrl)]
#     if(is.null(raw_fc)) raw_fc <- cnt[,well_trt]/cnt[,well_ctr]
#     else raw_fc <- cbind(raw_fc,cnt[,well_trt]/cnt[,well_ctr])
#   }
#   return(log(raw_fc))
# }
# 
# getFC_ruv <- function(cnt,grp,w = NULL){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   if(is.null(w)) design <- model.matrix(~grp)
#   else design <- model.matrix(~grp+w)
#   y_d <- estimateDisp(DGEx1f, design, robust=TRUE)
#   fit <- glmFit(y_d, design)
#   lrt <- glmLRT(fit, coef=2:ncol(fit$design))
#   tags <- topTags(lrt, n=Inf)$table
#   lfc_ind <- grep('logFC',colnames(tags))
#   lfc_ind <- lfc_ind[colnames(tags)[lfc_ind] != 'logFC.w']
#   tags2 <- tags[,lfc_ind]
#   rownames(tags2) <- tags$genes
#   tags2 <- tags2[rownames(cnt),]
#   colnames(tags2) <- gsub('logFC.grp','',colnames(tags2))
#   tags3 <- cbind('_._'=0,tags2)
#   ctr_ind <- grep('_._',colnames(tags3),fixed = T)
#   FC_table <- NULL
#   for(i in ctr_ind){
#     FC_table <- rbind(FC_table,cbind(rep(i,6),(i+1):(i+6)))
#   }
#   logFC <- NULL
#   for (i in 1:nrow(FC_table)){
#     fc_diff <- tags3[,FC_table[i,2]]-tags3[,FC_table[i,1]]
#     if(is.null(logFC)) logFC <- fc_diff
#     else logFC <- cbind(logFC,fc_diff)
#   }
#   colnames(logFC) <- colnames(tags3)[-ctr_ind]
#   rownames(logFC) <- rownames(tags3)
#   return(logFC) 
# }
# 
# get_pval_ruv <- function(cnt,grp,w = NULL){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   if(is.null(w)) {
#     design <- model.matrix(~grp)
#   }else {
#     design <- model.matrix(~grp+w)
#   }
#   y_d <- estimateDisp(DGEx1f, design, robust=TRUE)
#   fit <- glmFit(y_d, design)
#   pvals <- NULL
#   for(i in 2:ncol(fit$design)){
#     lrt <- glmLRT(fit, coef=i)
#     tags <- topTags(lrt, n=Inf)$table
#     tags2 <- tags$PValue
#     names(tags2) <- tags$genes
#     if(is.null(pvals)){
#       pvals <- tags2
#     }else{
#       pvals <- cbind(pvals,tags2)
#     }
#   }
#   colnames(pvals) <- colnames(fit$design)[-1]
#   return(pvals) 
# }
# 
# getFC <- function(cnt,grp,hyper_grp=NULL,w = NULL){
#   DGEx1 = DGEList(counts=cnt, genes = rownames(cnt),group = grp)
#   DGEx1f <- calcNormFactors(DGEx1)
#   if(is.null(hyper_grp)){
#     hyper_grp <- rep(1,length(grp))
#   }else{
#     hyper_grp <- hyper_grp
#   }
#   logFC <- NULL
#   for(i in unique(hyper_grp)){
#     DGEx1f_sub <- DGEx1f[,hyper_grp ==i]
#     grp_sub <- grp[hyper_grp ==i]
#     if(is.null(w)) {
#       design <- model.matrix(~grp_sub)
#       y_d <- estimateDisp(DGEx1f_sub, design, robust=TRUE)
#     }else{
#       w_sub <- w[hyper_grp ==i]
#       design0 <- model.matrix(~grp_sub)
#       design <- model.matrix(~grp_sub+w_sub)
#       y_d <- estimateDisp(DGEx1f_sub, design0, robust=TRUE)
#     }
#     fit <- glmFit(y_d, design)
#     lrt <- glmLRT(fit, coef=2:length(unique(grp_sub)))
#     tags <- topTags(lrt, n=Inf)$table
#     tags2 <- tags[,2:length(unique(grp_sub))]
#     rownames(tags2) <- tags[,1]
#     colnames(tags2) <- gsub('logFC.grp_sub','',colnames(tags2))
#     if(is.null(logFC)){
#       logFC <- tags2[rownames(cnt),]
#     }else{
#       logFC <- cbind(logFC,tags2[rownames(cnt),])
#     }
#   }
#   return(logFC)
# }
