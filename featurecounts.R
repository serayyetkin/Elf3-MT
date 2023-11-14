```library(Rsubread)
library(BiocParallel)
library(limma)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
setwd("/archive/yetkins/TEZ/bam_files")

#feature_counts####
bams <- dir(path = setwd("/archive/yetkins/TEZ/bam_files/"), pattern = ".bam")

counts <- featureCounts(files = bams,
                        annot.ext = "/archive/yetkins/ESC-MEF_RNA-seq/gencode.vM23.annotation.gtf", 
                        isGTFAnnotationFile= TRUE,
                        GTF.featureType= c("exon"),
                        GTF.attrType="gene_id",
                        useMetaFeatures=T,
                        countMultiMappingReads=F,
                        isPairedEnd=FALSE,
                        nthreads=5
)


counts_df <- as.data.frame(counts$counts)
counts_df[counts_df==0] <- NA
counts_df_2 <- na.omit(counts_df)
counts_df_2$gene_id <- row.names(counts_df_2) 
write.table(counts_df_2, "/archive/yetkins/TEZ/counts_Elf3.tsv", quote = F, col.names = T, sep = "\t")```
