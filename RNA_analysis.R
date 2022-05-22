## Code for RNA seq analysis
## @authors: Sneha Rao and Raju Tomer



library("DESeq2")
library("EnhancedVolcano")
library(biomaRt)
library(ViSEAGO)
library(data.table)
library(plotly)


#Define input data dir
setwd('D:/RNAData/Results')

# Read Data
countdata <-read.csv(file.path("D:/RNAData/InputFiles/FullGeneCounts.csv"), header = TRUE, stringsAsFactors=FALSE, row.names = 1)
coldata <- read.table(file.path("D:/RNAData/InputFiles/SampleMetaNet.txt"), header = TRUE, stringsAsFactors=FALSE)
nrow(countdata)
coldata

#### for Setd1 vs. Ctrl (change for different datasets)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~Batch+Genotype)

nrow(ddsMat)
keep <- rowSums(counts(ddsMat)) > 1
nrow(ddsMat)

dds <- DESeq(ddsMat)
resultsNames(dds)
res <- results(dds, contrast=c("Genotype", "setd1a", "wt"))
nrow(res)
res <- na.omit(res)
nrow(res)
summary(res)

## for annotation of genes
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version=99)

### Find up or down regulated genes
m1 = res$padj <= 0.05
sum(m1)
m2 = res$log2FoldChange > 0 ## choose appropriate cutoff
sum(m2)
m3 = res$log2FoldChange < 0
sum(m3)
m4 = m2 | m3
sum(m4)
m = m1 & m4
sum(m)

up_genes <- rownames(res[m1,])
genemap <- getBM( attributes = c("ensembl_gene_id"),
                  filters = "external_gene_name",
                  values = up_genes,
                  mart = ensembl )

selection_up <- genemap$ensembl_gene_id


genemap <- getBM( attributes = c("ensembl_gene_id"),
                  filters = "external_gene_name",
                  values = rownames(ddsMat),
                  mart = ensembl )

background <- genemap$ensembl_gene_id

################ GO analysis
fname = "BP_AllDGE_MoNNet_30vs15D_Batch" # output file name string
go_process = "BP"   ## BP, CC or MF
fname = paste(go_process, "_", fname, sep = "")
fname




###############################################################


## Make ENSEMBL object
biomart = "genes"
host = "https://jan2020.archive.ensembl.org"
Ensembl <- listEnsembl()
match.arg(biomart, Ensembl$biomart)
mart <- useEnsembl(biomart, host = host, version = 99)
Ensembl <- new("genomic_ressource", db = "Ensembl", stamp = paste(host, 
                                                                  Ensembl$version[Ensembl$biomart == biomart]), data = data.table(), 
               mart = list(mart), organisms = data.table(listDatasets(mart)))



#ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "mmusculus_gene_ensembl",
  Ensembl
)


################################ topGO 

##### TopGO analysis
BP_up <- ViSEAGO::create_topGOdata(
  geneSel=selection_up,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont=go_process,
  nodeSize=5
)

BP_down <- ViSEAGO::create_topGOdata(
  geneSel=selection_down,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont=go_process,
  nodeSize=5
)

############################################################################

##### Run tests

##### Setd1 vs wt
classic_up <-topGO::runTest(
  BP_up,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.01
)

classic_down <-topGO::runTest(
  BP_down,
  algorithm ="elim",
  statistic = "fisher",
  cutOff=0.01
)

##############################################################

###### Merge
BP_sResults<-ViSEAGO::merge_enrich_terms(cutoff = 0.01,
                                         Input=list(
                                           down=c("BP_down", "classic_down"),
                                           up=c("BP_up","classic_up")
                                         )
)

BP_sResults

#################


ViSEAGO::show_table(BP_sResults)

# print the merged table in a file
ViSEAGO::show_table(
  BP_sResults,
  paste(fname, ".xls", sep="")
)


# initialize 
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

# compute all available Semantic Similarity (SS) measures
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)

# GOterms heatmap with the default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels=FALSE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)


# Display the clusters-heatmap
a=ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)
a
orca(a,file = paste(fname, ".svg", sep=""))


Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)


# Display the clusters-heatmap
a=ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)
a
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms",
  paste(fname, "_lab", ".pdf", sep="")
)
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms",
  paste(fname, "_lab", ".svg", sep="")
)

orca(a,file = paste(fname, ".png", sep=""))



###############################################################


