library(DESeq2)
library(data.table)
library(dplyr)
library(KEGGprofile)
library(biomaRt)
library(plyr)
library(ggplot2)




sample_table <- function(directory){
  sampleName <- grep("counts",list.files(directory),value=TRUE)
  sampleCondition <- rep(directory, times=length(sampleName))
  sampleFiles = paste0(directory, "/", sampleName)
  return(cbind(sampleName, sampleFiles, sampleCondition))
}


sampleTable <- sub(".tsv", "", list.files(".", pattern = ".tsv")) %>% lapply(., sample_table) %>% do.call(rbind, .) %>% data.frame

sampleTable <- sampleTable[!duplicated(sampleTable$sampleName),]



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design= ~ sampleCondition)

dds <- DESeq(ddsHTSeq)

vst_dds <- vst(dds)

#counts.norm <- assay(vst_dds)
plotPCA(vst_dds, intgroup=c("sampleCondition")) + theme_bw()


c <- counts(dds, normalized=T)
rownames(c) <- gsub("\\..*$", "", rownames(c))

conv <- convertId(c, filters="entrezgene")

ens <- gsub("\\..*$", "", rownames(c))

length(unique(ens)) == length(ens)



mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

ensembl_genes<- rownames(c)

bm <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "entrezgene", "description"),
  values= ensembl_genes,
  mart= mart)

bm$description <- NULL

c <- data.frame(c)

c$ensgene <- rownames(c)

m <- merge(c, bm, by.x = "ensgene", by.y = "ensembl_gene_id")

m <- data.table(m)
m$ensgene <- NULL

a <- aggregate(. ~ entrezgene, data=m, median)

rownames(a) <- a$entrezgene
a$entrezgene <- NULL


plot_pathway(a[,!grepl("13", sampleTable$sampleCondition)], pathway_id = "04012",species='hsa', type="lines", result_name="12.png")


erbb <- c("1956", "2064", "2065", "2066")

erbb_expr <- a[rownames(a) %in% erbb,]


pl <- data.frame(list(erbb=colMeans(as.matrix(erbb_expr)), codon=sampleTable$sampleCondition, file=sampleTable$sampleName))

ggplot(data=pl, aes(x=codon,y=erbb,  fill=codon))+
  geom_violin()

## Pre-Filtering. Remove rows with few counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
