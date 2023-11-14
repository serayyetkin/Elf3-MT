counts <- read_delim("counts_Elf3.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
counts <- counts[,c(1,4,5,2,3)]
colnames(counts) <- c("gene_id","scr_1","scr_2","Elf3kd_1","Elf3kd_2")
gtf <- rtracklayer::import("~/Desktop/Bosphorus/gencode.vM23.annotation.gtf")
gtf_data_frame <- as.data.frame(gtf)

gtf_filtered <- select(gtf_data_frame,
                       c(10,12,11,1:3,5))
gtf_filtered_uniq <- gtf_filtered %>% distinct(gene_id, .keep_all= TRUE)
count_df_annotated <- merge(gtf_filtered_uniq,counts,by="gene_id")

count_df_annotated$gene_length <- count_df_annotated$end-count_df_annotated$start
count_df_annotated <- count_df_annotated[,c(1:6,12,7:11)]
write.table(count_df_annotated, "~/Desktop/TEZ/counts_Elf3_annotated.tsv", quote = F, col.names = T, sep = "\t")

rownames(count_df_annotated) <-count_df_annotated$gene_id

#calculate fpkm values
z <- DGEList(counts=count_df_annotated[,c(9:12)], genes=count_df_annotated[,c("gene_id","gene_length")])
z<- calcNormFactors(z)
#PCA <- plotMDS(z)
RPKM <-rpkm(z)
fp <- as.data.frame(RPKM)
fp$gene_id <- rownames(fp)
FPKM <- fp%>%select(gene_id,everything())
FPKM.annotated <- merge(FPKM, gtf_filtered_uniq, by="gene_id")
FPKM.annotated <- FPKM.annotated[,c(1,6:11,2:5)]
write.table(FPKM.annotated, "~/Desktop/TEZ/fpkm_Elf3_annotated.tsv", quote = F, col.names = T, sep = "\t")
write.csv(FPKM.annotated, "~/Desktop/TEZ/fpkm_Elf3_annotated.csv")
#dge
count_df <- count_df_annotated[,c(9:12)]
row.names(count_df) <- count_df_annotated$gene_id
group <- c(rep("control", 2), rep("treatment", 2)) #CHANGE this according to sample number
dge_1 <- DGEList(counts=count_df, group=group)
dge_1 <- calcNormFactors(dge_1)
dge_1 <- estimateCommonDisp(dge_1, verbose=TRUE)
dge_1 <- estimateTagwiseDisp(dge_1)
et_1 <- exactTest(dge_1, pair=c("control", "treatment"))
etp_1 <- topTags(et_1, n=nrow(et_1$table))
dgeList_1 <- etp_1$table
dgeList_1 <- tibble::rownames_to_column(dgeList_1, "gene_id")
group_elf3_annotated <- merge(dgeList_1, FPKM.annotated, by="gene_id")
group_elf3_annotated <- group_elf3_annotated[,c(1,6:11,2:5,12:15)]
row.names(group_elf3_annotated) <- group_elf3_annotated$gene_id
write.table(group_elf3_annotated, "~/Desktop/TEZ/Elf3kd_DGE_all.tsv", quote = F, col.names = T, sep = "\t")

EnhancedVolcano(group_elf3_annotated,
                lab = group_elf3_annotated$gene_name,
                x = 'logFC',
                y = 'PValue',
                title = 'Scr versus Elf3-kd',
                pointSize = 2.0,
                col=c('gray', 'lightgreen', 'lightblue', 'red'),
                labSize = 6.0,
                pCutoff = 0.05,
                FCcutoff = 1,
                selectLab = c('Cdh1','Elf3','Grhl1','Snai2','Zeb1','Itga5','Elf5','Uhrf1','Gm12925','Acta2','H2afx'))

group_elf3_annotated_Downreg <- group_elf3_annotated[group_elf3_annotated[,c("PValue")]<0.05 & group_elf3_annotated[,c("logFC")]< c(-1),]
group_elf3_annotated_Upreg <- group_elf3_annotated[group_elf3_annotated[,c("PValue")]<0.05 & group_elf3_annotated[,c("logFC")]> c(1),]

group_elf3_annotated_Downreg_protein <- group_elf3_annotated_Downreg[group_elf3_annotated_Downreg$gene_type == 'protein_coding',]
group_elf3_annotated_Downreg_lcnRNA <- group_elf3_annotated_Downreg[group_elf3_annotated_Downreg$gene_type == 'lncRNA',]

group_elf3_annotated_Upreg_protein <- group_elf3_annotated_Upreg[group_elf3_annotated_Upreg$gene_type == 'protein_coding',]
group_elf3_annotated_Upreg_lcnRNA <- group_elf3_annotated_Upreg[group_elf3_annotated_Upreg$gene_type == 'lncRNA',]


countdata <- count_df_annotated[,c(9:12)]
row.names(countdata) <- count_df_annotated$gene_id
countdata_ <- na.omit(countdata)

metadata <- read_delim("metadata.csv", delim = ";", 
                       escape_double = FALSE, trim_ws = TRUE)
dds <- DESeqDataSetFromMatrix(countData=countdata_, 
                              colData=metadata, 
                              design= ~ dex)
dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)
a <- plotPCA(vsdata, intgroup="dex")
aa <- a +theme_minimal()
aa
gene_list <- read_csv("gene_list.txt")
FPKM_gene <- merge(FPKM.annotated[,c(2,8:11)],gene_list, by='gene_name')
FPKM_gene_1 <- log2(FPKM_gene[,c(2:5)]+1)
rownames(FPKM_gene_1) <- FPKM_gene$gene_name
colnames(FPKM_gene_1) <- c('siControl_1','siControl_2','siElf3_1','siElf3_2')
pheatmap(FPKM_gene_1,margin = c(5,5),cluster_cols = FALSE,cluster_rows = TRUE,scale = "row",show_rownames = TRUE)

#FPKM_gene$scr_mean <- (FPKM_gene$scr_1+FPKM_gene$scr_2)/2
#FPKM_gene$siRNA_mean <- (FPKM_gene$Elf3kd_1+FPKM_gene$Elf3kd_2)/2
#FPKM_gene_2 <- log2(FPKM_gene[,c(6,7)]+1)
#rownames(FPKM_gene_2) <- FPKM_gene$gene_name
#pheatmap(FPKM_gene_2,margin = c(2,2),cluster_cols = FALSE,cluster_rows = TRUE,scale = "row",show_rownames = TRUE)

library(msigdbr)
library(fgsea)
msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)



rankData <- group_elf3_annotated$logFC
names(rankData) <- group_elf3_annotated$gene_name

# run fgsea enrichment
fgseaRes <- fgsea(pathways=pathwaysH, stats=rankData)
head(fgseaRes[order(pval), ])
fgseaRes %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, -padj)
library(data.table)
fwrite(fgseaRes, file="fgseaRes_elf3.txt", sep="\t", sep2=c("", " ", ""))


ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05),colour="black") +
  scale_fill_manual(values = c("lightblue", "#F8766D"))+
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +  theme_minimal()

plotEnrichment(pathwaysH[["HALLMARK_TGF_BETA_SIGNALING"]],
               rankData) + labs(title="HALLMARK_TGF_BETA_SIGNALING")

plotEnrichment(pathwaysH[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               rankData) + labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

plotEnrichment(pathwaysH[["HALLMARK_E2F_TARGETS"]],
               rankData) + labs(title="HALLMARK_E2F_TARGETS")

plotEnrichment(pathwaysH[["HALLMARK_APICAL_JUNCTION"]],
               rankData) + labs(title="HALLMARK_APICAL_JUNCTION")

plotEnrichment(pathwaysH[["HALLMARK_APICAL_SURFACE"]],
               rankData) + labs(title="HALLMARK_APICAL_SURFACE")

plotEnrichment(pathwaysH[["HALLMARK_ANGIOGENESIS"]],
               rankData) + labs(title="HALLMARK_ANGIOGENESIS")

plotEnrichment(pathwaysH[["HALLMARK_DNA_REPAIR"]],
               rankData) + labs(title="HALLMARK_DNA_REPAIR")

plotEnrichment(pathwaysH[["HALLMARK_HYPOXIA"]],
               rankData) + labs(title="HALLMARK_HYPOXIA")

library(reshape2)
corr_melted_FPKM.annotated <- round(cor(FPKM.annotated[,c(8:11)]),2)
melted_FPKM.annotated <- melt(corr_melted_FPKM.annotated)
head(melted_FPKM.annotated)
library(corrplot)
corrplot::corrplot(cor(melted_FPKM.annotated))
library(ggplot2)
ggplot(data = melted_FPKM.annotated, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()

organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

group_elf3_annotated_Downreg$gene_id_v <- gsub("\\..*", "", group_elf3_annotated_Downreg$gene_id)

ego_d <- enrichGO(gene     = group_elf3_annotated_Downreg$gene_id_v,
                       OrgDb         = organism,,
                       keyType       = 'ENSEMBL',  #specify the 'keytype' to the input ID type
                       ont           = "BP",
                       pAdjustMethod = "BH")
ego_df <- as.data.frame(ego_d)
barplot(ego_d, drop=TRUE, showCategory=28,color = "pvalue",font.size = 8)


group_elf3_annotated_Upreg$gene_id_v <- gsub("\\..*", "", group_elf3_annotated_Upreg$gene_id)

ego_u <- enrichGO(gene     = group_elf3_annotated_Upreg$gene_id_v,
                  OrgDb         = organism,,
                  keyType       = 'ENSEMBL',  #specify the 'keytype' to the input ID type
                  ont           = "BP",
                  pAdjustMethod = "BH")
ego_df <- as.data.frame(ego)
barplot(ego_u, drop=TRUE, showCategory=28,color = "pvalue",font.size = 8)
