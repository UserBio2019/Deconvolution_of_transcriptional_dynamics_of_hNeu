############################################################################
## This piece of script contains the batch correction and the             ##
## differential gene expression analysis performed on bulk RNA-Seq data   ##
## relative to the paper: "Cellular and transcriptional dynamics of human ##
## neutrophils at steady state and upon stress". It refers to mature      ##
## neutrophils and can be used to analyze the transcriptional changes     ##
## of monocytes without substantial modifications. Please refer to the    ##
## paper (materials and methods) for details on thresholds and filters.   ##
## For differential gene expression analysis only the comparison mature   ## 
## cells (Neu) at steady state vs GCSF treated cells has been reported    ##
## for simplicity but the same template can be used for the other         ##
## comparisons considered in the original paper. For any doubt and        ##
## questions don't hesitate to contact me at:                             ##
## lusito.eleonora83@gmail.com                                            ##
############################################################################


##################### BATCH CORRECTION MATURE Neutrophils Immature Neutrophils and Monocytes #############################

library(gplots)
library(yarn)  
library(gplots)
library(yarn)  
library(edgeR)
library(ggplot2)
library(DESeq2)
library(factoextra) 
library(ComplexHeatmap)
library(data.table)
library(ggrepel)  
library(sva)


      
load("Bulk210S.Robj")
	isexpr = (rowSums(cpm(Bulk210S) > 1) > 32)
    expr = Bulk210S[isexpr, , keep.lib.sizes=FALSE]
    expr$samples$lib.size = colSums(expr$counts)
    expr_N = calcNormFactors(expr)
    expr_Disp = estimateDisp(expr_N)
    expr_CommonDisp = estimateCommonDisp(expr_N)
    Scaling_factors = as.data.frame(expr_N$samples[,3])
    N = expr_N$samples$lib.size 
    qnorm = expr_N$samples$norm.factors
    scale_factor = N*qnorm/exp(mean(log(N*qnorm)))
    scaling_factors = data.frame(1/scale_factor)
    rownames(scaling_factors) = rownames(expr_N$samples)
    names(scaling_factors)="scaling_factors"


names(expr_N$genes) = c("RefSeq", "gene.length", "Gene")	
cpm_expr = cpm(expr_N, normalized.lib.sizes=T)
rpkm_expr = rpkm(expr_N, normalized.lib.sizes=T, log=F)
rpkm_log = rpkm(expr_N, normalized.lib.sizes=T, log=T)


df1=read.delim("Metadata_NDN_LDN_MONO.txt", header= T)
Batch = df1$Index4
Bio_group = df1$Index1
 
adjusted = ComBat_seq(dgel_expr_N$counts, batch=Batch, group=Bio_group, full_mod=TRUE)
dgel_expr_N$counts = adjusted
names(dgel_expr_N$genes) = c("RefSeq", "gene.length", "Gene")
rownames(dgel_expr_N$counts) = dgel_expr_N$genes$Gene
save(dgel_expr_N, file = "Combat_Bulk210S.Robj")

###########################################################################################################################
 
 
###################### BATCH CORRECTION MATURE Neutrophils and Monocytes ##################################################
 
load("Bulk_NDN_MONO_128S.Robj")
    isexpr = (rowSums(cpm(dgel) > 1) > 19)
    dgel_expr = dgel[isexpr, , keep.lib.sizes=FALSE]
    dgel_expr$samples$lib.size = colSums(dgel_expr$counts)
    dgel_expr_N = calcNormFactors(dgel_expr)
    dgel_expr_Disp = estimateDisp(dgel_expr_N)
    dgel_expr_CommonDisp = estimateCommonDisp(dgel_expr_N)
    Scaling_factors = as.data.frame(dgel_expr_N$samples[,3])
    N = dgel_expr_N$samples$lib.size 
    qnorm = dgel_expr_N$samples$norm.factors
    scale_factor = N*qnorm/exp(mean(log(N*qnorm)))
    scaling_factors = data.frame(1/scale_factor)
    rownames(scaling_factors) = rownames(dgel_expr_N$samples)
    names(scaling_factors)="scaling_factors"

names(dgel_expr_N$genes) = c("RefSeq", "gene.length", "Gene")	
cpm_expr = cpm(dgel_expr_N, normalized.lib.sizes=T)
rpkm_expr = rpkm(dgel_expr_N, normalized.lib.sizes=T, log=F)
rpkm_log = rpkm(dgel_expr_N, normalized.lib.sizes=T, log=T)


df1=read.delim("Metadata_NDN_MONO.txt", header= T)
Batch = df1$Index4
Bio_group = df1$Index1
 
adjusted = ComBat_seq(dgel_expr_N$counts, batch=Batch, group=Bio_group, full_mod=TRUE)
dgel_expr_N$counts = adjusted
names(dgel_expr_N$genes) = c("RefSeq", "gene.length", "Gene")
rownames(dgel_expr_N$counts) = dgel_expr_N$genes$Gene
save(dgel_expr_N, file = "Combat_Bulk_NDN_MONO_128S.Robj")



######################### Differential gene expression analysis ####################################################
#DGEA: NDN(HDvsPostG)

load("NDN_70S.Robj")      
dgel_expr_N = NDN_70S
Data_comnparisons = read.delim("Comparisons_NDN.txt", header = T)
Data_comnparisons = cbind(Data_comnparisons,dgel_expr_N$samples$group)
Comparison_selection_1 = Data_comnparisons[Data_comnparisons$Group1 != 0,]

    names.use = colnames(dgel_expr_N$counts)[colnames(dgel_expr_N$counts) %in% Comparison_selection_1[,1]]
    df.subset_1 = dgel_expr_N$counts[, names.use]
    df.subset_1 = setcolorder(as.data.frame(df.subset_1), as.character(Comparison_selection_1[,1]))

group = Comparison_selection_1$Group1
group = factor(group, levels = c("A", "B"))
   batch = Comparison_selection_1[,7]
   batch = factor(batch, levels =c("NextS_NovaS1", "NovaS2"))
   design = model.matrix(~ 0 + group + batch) 
   colnames(design) = gsub("group", "", colnames(design))
       Samples = dgel_expr_N$samples[rownames(dgel_expr_N$samples) %in% colnames(df.subset_1), ]
       dgel = DGEList(counts=df.subset_1, genes=dgel_expr_N$genes)
       dgel$samples$norm.factors = Samples$norm.factors
       isexpr = (rowSums(cpm(dgel) > 1) > 10)
       dgel_expr = dgel[isexpr,]
       rpkm_log_PostG_expressed = rpkm(dgel_expr, normalized.lib.sizes=T, log=T)
       names(dgel_expr$genes) = c("RefSeq", "gene.length", "Gene")
       dgel_expr$samples$lib.size = colSums(dgel_expr$counts)
   dgel_expr$samples$group = c(rep("2", 17), rep("1", 17))
y = estimateDisp(dgel_expr, design, robust=TRUE)
    pdf('Common_dispersion_NDN_PostG.pdf', 12, 16, paper= "USr")  
    plotBCV(y)
    dev.off()
fit = glmQLFit(y, design, robust=TRUE)
    pdf('GroupHD_PostG.pdf', 12, 16, paper= "USr") 
    plotQLDisp(fit)
    dev.off()
con = makeContrasts(A - B, levels = colnames(design))
qlf = glmQLFTest(fit, contrast=con)
Group1 =  topTags(qlf , n=nrow(df.subset_1), sort.by="logFC")
median_gex = t(fit$counts)
cb = cbind(as.data.frame(Comparison_selection_1$Group1), median_gex)
Gex_by_groups= data.frame(do.call("rbind", by(cb[,c(2:6534)],cb[,1],colMeans)))
Gex_by_groups = t(Gex_by_groups)



#Significant DEG selection
Group1_deg = rownames(Group1$table[which(Group1$table$FDR < 0.05 & abs(Group1$table$logFC)  >= 1.5),])
Group1_deg = as.data.frame(Group1_deg)
mtch = match(Group1$table[,3], Group1_deg[,1])
    cb = cbind(mtch, Group1$table)
    cb =  cb[complete.cases(cb[,1]), ] 
    cb = cb[,-1]
    cb = cb[,-c(1,2,5,6)]
    cb= cb[,c("Gene","logFC","PValue", "FDR")]
    cb_ = cb[order(row.names(cb)),] 
names(cb_) = c("Gene", "logFC_NDN_HD_PostG", "Pvalue_NDN_HD_PostG", "FDR_NDN_HD_PostG")




