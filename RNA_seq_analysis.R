#Ovchinnikov YB1 RIP analysis GEO accession: GSE130781

#Files for Analysis (Illumina HiSeq 2000)
#YB1-RIP SRR_05/6
#YB1 WT HEK 293 SRR12/13
#YB1 KO HEK293 SRR19/20

library(DESeq2)
library(EnhancedVolcano)
library(dplyr)


#### Read in RNA-seq bam file and form separate matrix for HEK cells and YB1 RIP####
setwd("C:/Users/12298/Desktop/Data_Analytics/Ovchinnikov_2020/Genes_Reads/")
files<-dir(path=".",pattern="tab")
gene.id<-read.table(file=files[1],header = F)[,1] #takes first column of first file
YB1 <-matrix(NA,ncol=6,nrow=length(gene.id)) #Generates blank matrix to insert files into during loop (same length as gene.id.
for(i in 1:6){
  YB1[,i] = read.table(file=files[i],header=F)[,2] #takes column two with +/- strand combination
}

colnames(YB1) = c("YB1_RIP","YB1_RIP","HEK_WT","HEK_WT", "HEK_KO", "HEK_KO")
rownames(YB1) = gene.id
YB1=YB1[5:(length(gene.id)),]

#####Pulling gene ID from biomaRt and merging with YB1 matrix#####
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(YB1)
G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "hgnc_symbol"),
                values = genes, mart = mart)

#Dividing matrix into three (wt HEK and YB1 KO HEK, YB1_RIP and wt HEK, and YB1_RIP alone) 
YB1_DE <- YB1[,3:6] %>% pracma::flipdim(dim=2)
YB1_RIP_DE <- YB1[,1:4] 
YB1_RIP <- YB1[,1:2] #picked back up at line 67


#############################################
######YB1 wt and YB1 KO cell DE Analysis#####
#############################################
sample=colnames(YB1_DE)
group<-c("HEK_KO","HEK_KO", "HEK_WT","HEK_WT")
coldata<-data.frame(sample, group)

#Running DEseq
dds<-DESeqDataSetFromMatrix(countData = YB1_DE, colData = coldata, design= ~ group)
dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds,contrast = c("group","HEK_KO","HEK_WT"))

#Making data frame with Fold change, Gene name, p-value
YB1_DE <- as.data.frame(res)[complete.cases(as.data.frame(res)),] %>% tibble::rownames_to_column( "ensembl_gene_id")
YB1_DE <- merge(YB1_DE, G_list, by="ensembl_gene_id")
saveRDS(YB1_DE, "YB1_DE.rds")

#Volcano plot of DE
EnhancedVolcano(YB1_DE, lab=subset(YB1_DE, hgnc_symbol == "PARP1"), x='log2FoldChange', y='padj', xlim=c(-8, 10), ylim=c(0, 170), pCutoff=0.05)

#############################
#####YB1-RIP DE Analysis#####
#############################
YB1_RIP_NZ <- YB1_RIP[YB1_RIP[,1] != 0 & YB1_RIP[,2] != 0,] #Exclude non-zero values
YB1_RIP_NZ <- cbind(YB1_RIP_NZ, apply(YB1_RIP_NZ, 1, mean)) #take mean of two replicates 
YB1_RIP_NZ <- subset(YB1_RIP_NZ, YB1_RIP_NZ[,3] > 100) #Exclude values < 100
YB1_RIP_DE <- subset(YB1_RIP_DE, rownames(YB1_RIP_DE) %in% rownames(YB1_RIP_NZ)) #Extracting Genes with RIP values > 1000

#DEseq analysis
sample=colnames(YB1_RIP_DE)
group<-c("YB1_RIP","YB1_RIP","HEK_WT","HEK_WT")
coldata<-data.frame(sample, group)

dds<-DESeqDataSetFromMatrix(countData = YB1_RIP_DE, colData = coldata, design= ~ group)
dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds,contrast = c("group","YB1_RIP","HEK_WT"))

#Making data frame with Fold change, Gene name, p-value
YB1_RIP_DE <- as.data.frame(res)[complete.cases(as.data.frame(res)),] %>% tibble::rownames_to_column( "ensembl_gene_id")
YB1_RIP_DE <- merge(YB1_RIP_DE, G_list, by="ensembl_gene_id")
saveRDS(YB1_RIP_DE, "YB1_RIP_DE.rds")

#Excluding fold change < 0
YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$log2FoldChange > 0)

#Volcano Plot for YB1_RIP DE
EnhancedVolcano(YB1_RIP_DE, lab=YB1_RIP_DE$hgnc_symbol, x='log2FoldChange', y='padj', xlim=c(0, 2.5), pCutoff=0.05, FCcutoff = 0.5)


######RIP enrichment and Correlation with DE genes######
setwd("C:/Users/12298/Desktop/Data_Analytics/Ovchinnikov_2020/")
YB1_DE <- readRDS("YB1_DE.rds")
YB1_RIP_DE <- readRDS("YB1_RIP_DE.rds")

YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$log2FoldChange > 0)
YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$padj < 0.05)
YB1_DE <- subset(YB1_DE, YB1_DE$padj < 0.05)
YB1_DE <- subset(YB1_DE, YB1_DE$ensembl_gene_id %in% YB1_RIP_DE$ensembl_gene_id)
YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$ensembl_gene_id %in% YB1_DE$ensembl_gene_id)

#Creating DF for ggplot
RIP_DE_GGplot <- cbind(YB1_DE[, c(3, 7, 8)], YB1_RIP_DE[, c(3, 7, 8)]) %>% `colnames<-`(c("HEK_DE_FC", "HEK_DE_pVal", "HEK_DE_Gene", "RIP_DE_FC", "RIP_DE_pVal", "RIP_DE_Gene"))

#Creating categorical variables to add colors and labels on ggplot 
RIP_DE_GGplot$FC_Cat <- ifelse(abs(RIP_DE_GGplot$HEK_DE_FC) > 1, "-1 < FC > 1", "-1 > FC < 1")
RIP_DE_GGplot$pVal_Cat <- ifelse(RIP_DE_GGplot$HEK_DE_pVal < 0.05, "p < 0.05", "p > 0.05")

#Plotting YB1 RIP vs HEK DE 
ggplot(RIP_DE_GGplot, aes(x=HEK_DE_FC, y=RIP_DE_FC, color = as.factor(FC_Cat), shape = as.factor(pVal_Cat), label=HEK_DE_Gene)) + geom_point(size=1, alpha=0.65) + xlim(c(-9, 3)) + labs(y="RIP Enrichment", x="Differential Seq Genes") + ggtitle("RIP Enrichment vs YB1 HEK Differential") + geom_vline(xintercept = c(-1, 1), linetype = "dotted") + geom_text(aes(label=ifelse(abs(HEK_DE_FC) > 1, as.character(HEK_DE_Gene), '')), vjust=2, hjust=0.5, size = 2, color="black") + theme(legend.position="none")

  
#Correlation
YB1_DE <- subset(YB1_DE, YB1_DE$log2FoldChange <0)
YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$log2FoldChange > 0)
YB1_DE <- subset(YB1_DE, YB1_DE$ensembl_gene_id %in% YB1_RIP_DE$ensembl_gene_id)
YB1_RIP_DE <- subset(YB1_RIP_DE, YB1_RIP_DE$ensembl_gene_id %in% YB1_DE$ensembl_gene_id)

cor.test(x=YB1_DE$log2FoldChange, y=YB1_RIP_DE$log2FoldChange)

RIP_DE_Cor <- cbind(YB1_DE[, c(3, 7, 8)], YB1_RIP_DE[, c(3, 7, 8)]) %>% `colnames<-`(c("HEK_DE_FC", "HEK_DE_pVal", "HEK_DE_Gene", "RIP_DE_FC", "RIP_DE_pVal", "RIP_DE_Gene"))

ggplot(RIP_DE_Cor, aes(x=HEK_DE_FC, y=RIP_DE_FC)) + geom_point(size = 1) + ggpubr::stat_cor(method="pearson") + labs(x="HEK Differential Expression", y="RIP Enrichment") + geom_smooth(method="lm", se=FALSE) + ggtitle("RIP Enrichment : HEK DE Correlation") + ggsave("Correlation.pdf")


###############################
######GSEA of DE analysis######
###############################

#Load libraries and pathways 
library(fgsea)
Hallmark_Pathways <- gmtPathways("C:/Users/12298/Desktop/Data_Analytics/GSEA_Pathways/Standard.gmt")
Hallmark_Pathways %>% head() %>% lapply(head)


#Remove insignificant data points
YB1_DE <- readRDS("YB1_DE.rds")
P_SIG <- subset(YB1_DE, padj<0.05)#Exclude p>0.05


#Create ranks from the Fold Change
res2 <- P_SIG %>% dplyr::select(hgnc_symbol, stat) %>% na.omit() %>% distinct() %>% group_by(hgnc_symbol) %>% summarize(stat=mean(stat))
ranks <- tibble::deframe(res2)


#Perform GSEA
fgseaRes <- fgsea(pathways = Hallmark_Pathways, stats = ranks)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy$pathway <- gsub("^.{0,9}", "", fgseaResTidy$pathway)


#Plotting GSEA
fgseaResTidy <- subset(fgseaResTidy, fgseaResTidy$padj<0.05)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col() + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA") + theme_minimal()


###################################
######GSEA of RIP_DE analysis######
###################################

#Load libraries and pathways 
library(fgsea)
Hallmark_Pathways <- gmtPathways("C:/Users/12298/Desktop/Data_Analytics/GSEA_Pathways/Cancer.gmt")
Hallmark_Pathways %>% head() %>% lapply(head)

#Remove insignificant data points
YB1_RIP_DE <- readRDS("YB1_RIP_DE.rds")
P_SIG <- subset(YB1_RIP_DE, padj<0.05)#Exclude p>0.05
P_SIG <- subset(P_SIG, log2FoldChange > 0)

#Create ranks from the Fold Change
res2 <- P_SIG %>% dplyr::select(hgnc_symbol, stat) %>% na.omit() %>% distinct() %>% group_by(hgnc_symbol) %>% summarize(stat=mean(stat))
ranks <- tibble::deframe(res2)


#Perform GSEA
fgseaRes <- fgsea(pathways = Hallmark_Pathways, stats = ranks)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
#fgseaResTidy$pathway <- gsub("^.{0,9}", "", fgseaResTidy$pathway)

# Show in a nice table:
GSEA_ANALYSIS <- as.data.frame(fgseaResTidy) 


#Plotting GSEA
fgseaResTidy <- subset(fgseaResTidy, fgseaResTidy$pval < 0.05)
fgseaResTidy <- subset(fgseaResTidy, fgseaResTidy$NES >0)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col() + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score", title="Cancer pathways NES from GSEA") + theme_minimal()


######END#######













