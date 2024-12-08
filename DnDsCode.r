#### Exercise 3: Compare dN/dS ratios among gene-sets

setwd("/Users/martapervil/Desktop/PGB/Project/dn_ds")

## 1- Read the table

d <- read.table("HumanDnDsW.txt", header=T, sep="\t")

## 2- Attach gene names

install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listDatasets(ensembl)
listAttributes(ensembl)
listFilters(ensembl)

geneInfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                  filters ='ensembl_gene_id', values = d$GeneID, mart = ensembl)
#Lo que está haciendo aquí es descargar de ensembl el id, nombre, cromosma, incio y final
#de cada uno de los genes anotados en la tabla HumanDnDsW.txt

d <- left_join(d, geneInfo, by=c("GeneID"="ensembl_gene_id"))
#Aquí tenemos una tabla con los genes de HumanDnDsW.txt y sus respectivos Dn, Ds,
#w, id, nombres...


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instala GO.db desde Bioconductor
BiocManager::install("GO.db")

library("GO.db")
library(annotate)

GOTERM[["GO:0001813"]]
xx <- as.list(GOTERM)
names(xx)

Term(xx[[1]])
terms <- lapply(xx, Term)

immune_terms <- grep("complement", terms, ignore.case = TRUE)
immune_terms <- unlist(terms[immune_terms])
names(immune_terms)


splicing_terms <- grep("splicing", terms, ignore.case = TRUE)
splicing_terms <- unlist(terms[splicing_terms])
names(splicing_terms)

neurogenesis_terms <- grep("neurogenesis", terms, ignore.case = TRUE)
neurogenesis_terms <- unlist(terms[neurogenesis_terms])
names(neurogenesis_terms)

immune_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                   filters = 'go', values = names(immune_terms), mart = ensembl)

splicing_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                     filters = 'go', values = names(splicing_terms), mart = ensembl)               

neurogenesis_gene <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                       filters = 'go', values = names(neurogenesis_terms), mart = ensembl)               

#### Subset dnds table using these two sets of genes


df_immuno <- filter(d, GeneID %in% unique(immune_gene$ensembl_gene_id))[,c(1,2,3,5,6)]
df_immuno$Geneset <- "Immuno"

df_splicing <- filter(d, GeneID %in% unique(splicing_gene$ensembl_gene_id))[,c(1,2,3,5,6)]
df_splicing$Geneset <- "Splicing"

df_neurogenesis <- filter(d, GeneID %in% unique(neurogenesis_gene$ensembl_gene_id))[,c(1,2,3,5,6)]
df_neurogenesis$Geneset <- "Neurogenesis"

df <- rbind(df_immuno, df_splicing, df_neurogenesis)


library(ggplot2)
library(ggstatsplot)
library(hrbrthemes)
library(viridis)
library(dplyr)

df %>% ggplot( aes(x=Geneset, y=Homo_sapiens.w , fill=Geneset)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + ylim(c(0,1))
