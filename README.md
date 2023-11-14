# Elf3-MT
# ** ELF3 SILENCING IN MESENCHYMAL EPITHELIAL TRANSITION**

## Aim of the study
 This Project mainly aims to investigate the impact of Elf3 silencing on Mesenchymal epithelial transition in NMuMG cell line using an integrated approach, including gene modulation and bioinformatics. This github page show

## Data Collection
 The experiments utilized the Normal Murine Mammary Gland cell line (NMuMG), procured from ATCC (ATCC® CRL-1636) in 2017. This cell line is a widely recognized and frequently employed model for studying EMT and MET (Miettinen et al., 1994). This model is favored due to its distinct gene signature alterations manifested in observable morphological changes. Upon exposure to TGFβ, these cells undergo EMT labeled as TGFβ3, characterized by the downregulation of Cdh1, the rearrangement of cortical actin to stress fibers, and the upregulation of mesenchymal genes such as N-cadherin and Vimentin. Notably, upon removal of TGFβ from the culture medium, these cells transition back to the epithelial state via MET within 72 hours labeled as PT72 (post-transition) (Alotaibi et al., 2015).

# <img src='https://github.com/serayyetkin/Elf3-MT/assets/73422665/9ae07972-bb33-4b11-9c6a-ed2377f53c23' width="400" height="200">

 siRNAs against Elf3 were purchased from Qiagen for siElf3_2; the sequence is: TTGAACCAACTTGTTCGATAA. NMuMG cells were transfected with 50 nM siRNA using the Lipofectamine RNAi Max as recommended by the manufacturer and collected 72 h later. Total RNA was isolated from treated and control cells using the Nucleospin RNA II kit. Concentrations of samples were measured using the instrument NanoDrop 2000. 1 µg RNA was converted to cDNA using the Maxima First Strand cDNA Synthesis Kit. cDNA was then diluted 1:10 and used for qPCR. AMPIGENE qPCR probe mix (Enzo) and Universal Probe Library probes were used for the experiment as recommended by the manufacturer. Analysis and quantification of qPCR were performed using the ΔΔCt method. Validated samples were shipped to sequencing company.

 All data collection and processing operations were conducted on an Apple MacBook Air, 8 GB RAM computer with MacOS Monterey 12.6.3 operating system, and IBG HPC environments have been used for data collection in this study.

## RNA-Seq Analysis 

```mermaid
graph LR;
    FastQC-->Trimmomatic;
    Trimmomatic-->HiSAT2;
    HiSAT2-->Featurecounts;
    Featurecounts-->EdgeR;
    EdgeR-->GSEA;
    GSEA-->Visualisation;
```
 As shown in the figure analysis starts with quality control and FastQC tool is used. Then adaptor and quality trimming is performed with trimmomatic. To align sequences to genome HISAT2 tool is used. Quantification performed with featurecounts and differentially expressed genes defined with EdgeR. Gene set enrichment analysis performed with fgsea R package. Visualisation and statistic analysis performed in R.

