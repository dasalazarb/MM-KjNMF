options(connectionObserver = NULL)
library(clusterProfiler)
# require(GenomicRanges)
# library(fgsea)
# library(DOSE)
# library(limma)
# library(edgeR)
library(dplyr)
library(ggplot2)
# library(magrittr)
#library(pathview)
#library(enrichplot)
# library(ReactomePA)
library(org.Hs.eg.db)
library(AnnotationHub)

# _____________________________________________ #
#### ........... Cargar archivos ........... ####
# _____________________________________________ #
nombre <- commandArgs(trailingOnly = TRUE)

pathFeatureLabel <- strsplit(nombre, "__")[[1]][1]
nameDotPlot <- strsplit(nombre, "__")[[1]][2]
iteration_number <- strsplit(nombre, "__")[[1]][3]
# pathFeatureLabel <- "C:/Users/da.salazarb/Downloads/Nueva_carpeta/kmmjnmf/pathFeatureLabel"
# nameDotPlot <- ""; iteration_number <- 0

# _________________________________________ #
#### ........... Co-modulos ........... ####
# ________________________________________ #
## detect the dir of mrna e.g. 2_mrna
dir_mrna <- dir(pathFeatureLabel)[grepl("^[0-9]*\\_mrna$",dir(pathFeatureLabel))]
path <- paste0(pathFeatureLabel, "/", dir_mrna)
archivos_ <- list.files(path = path, pattern = "\\.txt$")
archivos <- c()
for (archivo in archivos_)  {
  if (file.info(paste0(path,"/",archivo))$size > 0) {
    archivos <- c(archivos, archivo)
  }}

geneItems <- function(path, archivo) {
  # > archivos[2]
  # [1] "tcga_mrna_depurado.csv"
  tcga <- read.csv(paste0(path,"/",archivo),header = FALSE)
  tcga$V1 <- as.vector(tcga$V1)
  tcga$logFC <- 0 # as.vector(sapply(tcga, function(x) DEA[x,]$logFC))
  tcga$logFC[is.na(tcga$logFC)] <- 0
  
  ## bitr: Biological Id TranslatoR
  eg <-  bitr(gsub("mrna_","",tcga$V1), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  eg$SYMBOL <- paste0("mrna_",eg$SYMBOL)
  for (j in 1:length(tcga$V1)) {
    if (tcga$V1[j] %in% eg$SYMBOL) {
      tcga$V1[j] <- eg$ENTREZID[grep(tcga$V1[j],eg$SYMBOL)][1]
    } else {
      tcga$V1[j] <- gsub("mrna_","",tcga$V1[j])
    }
  } 
  
  # ## feature 1: numeric vector -> numeric vector: fold change or other type of numerical variable
  geneList <- tcga[,2]
  
  # ## feature 2: named vector -> named vector: every number was named by the corresponding gene ID
  names(geneList) <- as.character(tcga[,1])
  
  # ## feature 3: decreasing order
  geneList <- sort(geneList, decreasing = TRUE)
  gene <- names(geneList) #[abs(geneList) > 1]
  return(list(geneList=geneList, gene=gene))
}

##### Cluster biological comparison #####
gcSample <- list()

for (i in 1:length(archivos)) {
  archivo <- archivos[i]
  itemGenes <- geneItems(path, archivo)
  archivo <- gsub(paste0(dir_mrna, "_co-md_"), "c-", gsub(".txt", "", archivos[[i]]))
  geneList <- itemGenes$geneList
  archivo <- paste0("c-", as.character(as.integer(gsub("c-", "", archivo))+1))
  gcSample[[archivo]] <- itemGenes$gene
}

# _______________________________________________________________ #
#### ........... KEGG comparison for gene clusters ........... ####
# _______________________________________________________________ #
# ## Solucion a enrichGO: http://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
hub <- AnnotationHub()
quer <- query(hub, "Homo Sapiens")
OrgDb <- hub[[quer$ah_id[grep("OrgDb",quer$rdataclass)]]]
 
# ck1 <- compareCluster(geneCluster = gcSample, fun = "enrichGO", OrgDb=OrgDb)
ck2 <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG", pvalueCutoff=0.05)

# terms2remove <- c("Human immunodeficiency virus 1 infection", "Yersinia infection", "Hepatitis B", "Cushing syndrome", "Neurotrophin signaling pathway", 
#                   "Osteoclast differentiation", "Human T-cell leukemia virus 1 infection", "Legionellosis", "Epstein-Barr virus infection", "Prolactin signaling pathway", 
#                   "Measles", "Fluid shear stress and atherosclerosis", "Leishmaniasis", "Longevity regulating pathway", "Viral myocarditis", "Lipid and atherosclerosis", 
#                   "Malaria", "Hepatitis C", "Kaposi sarcoma-associated herpesvirus infection", "Tuberculosis", "Shigellosis", "Coronavirus disease - COVID-19", "Oocyte meiosis", 
#                   "Viral carcinogenesis", "Diabetic cardiomyopathy", "Intestinal immune network for IgA production", "Insulin resistance", "Axon guidance", "Alzheimer disease", 
#                   "Epithelial cell signaling in Helicobacter pylori infection", "Chagas disease","Pertussis", "Pathogenic Escherichia coli infection", 
#                   "Non-alcoholic fatty liver disease", "Human papillomavirus infection", "Longevity regulating pathway - multiple species", "Salmonella infection", 
#                   "Influenza A", "Rheumatoid arthritis", "Herpes simplex virus 1 infection", "Amoebiasis", "African trypanosomiasis", "Type II diabetes mellitus", 
#                   "Viral protein interaction with cytokine and cytokine receptor", "Vascular smooth muscle contraction", "Homologous recombination", 
#                   "Bacterial invasion of epithelial cells", "Allograft rejection", "Parathyroid hormone synthesis, secretion and action", "Apoptosis - multiple species", 
#                   "Insulin signaling pathway", "Human cytomegalovirus infection", "Platelet activation", "Toxoplasmosis")
# ck2@compareClusterResult <- ck2@compareClusterResult[!ck2@compareClusterResult$Description %in% terms2remove,]

pdf(file = paste0(pathFeatureLabel,"/co-mod_R_plots/", nameDotPlot, "_enrichKEGG.pdf"), width = 20, height = 18)
# jpeg(filename = paste0(pathFeatureLabel,"/co-mod_R_plots/", nameDotPlot, "_enrichKEGG.jpg"), width = 20, height = 15, units = "in", res = 350)
dotplot(ck2,font.size=12)
dev.off()

## guardar records_mrna.csv, se usa en 16_clinical_info_TCGABiolinks.R y en 22_Similarity_ccle_tcga_clusters.R
ck2 <- as.data.frame(ck2)

for (i in 1:length(ck2$geneID)) {
  ck2$geneID[i] <- paste(paste0("mrna_", bitr(unlist(strsplit(ck2$geneID[i], "/")), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL), collapse = "/")
}

write.table(ck2, paste0(pathFeatureLabel,"/co-mod_records/record_mrna.csv"), sep=",", row.names = FALSE)

####
# ck2 <- read.csv(paste0(pathFeatureLabel,"/co-mod_records/record_mrna.csv")); dim(ck2)

if (identical(FALSE, file.exists(paste0(pathFeatureLabel,"/co-mod_records/K", length(archivos_),".csv")))) { 
  ck2.temp <- ck2 %>% 
    dplyr::select(ID, Description, GeneRatio, Count) %>% 
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>% 
    dplyr::bind_cols(count.terms=1) %>% 
    dplyr::group_by(ID, Description) %>% 
    dplyr::summarise(dplyr::across(c("GeneRatio", "Count"), list(mean=mean)), count.terms_sum = sum(count.terms)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(iterNum_count.terms = paste0(iteration_number, ".",count.terms_sum, "_")) %>% 
    # dplyr::mutate(iterNum_count.terms = paste0(iteration_number, "_",count.terms_sum,",",round(GeneRatio_mean,2),",",round(Count_mean,2),"__")) %>% 
    dplyr::rename(GeneRatio=GeneRatio_mean, Count=Count_mean, count.terms=count.terms_sum)
  
  write.table(ck2.temp, paste0(pathFeatureLabel,"/co-mod_records/K", length(archivos_),".csv"), sep=",", row.names = FALSE)
  
} else {
  
  ck2.temp1 <- read.csv(paste0(pathFeatureLabel,"/co-mod_records/K", length(archivos_),".csv")); head(ck2.temp1, 10)

  ####### actual ck2 -> generate iterNum_count.terms
  ck2.temp <- ck2 %>% 
    dplyr::select(ID, Description, GeneRatio, Count) %>% 
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>% 
    dplyr::bind_cols(count.terms=1) %>% 
    dplyr::group_by(ID, Description) %>% 
    dplyr::summarise(dplyr::across(c("GeneRatio", "Count"), list(mean=mean)), count.terms_sum = sum(count.terms)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(iterNum_count.terms = paste0(iteration_number, ".",count.terms_sum, "_")) %>% 
    # dplyr::mutate(iterNum_count.terms = paste0(iteration_number, "_",count.terms_sum,",",round(GeneRatio_mean,2),",",round(Count_mean,2),"__")) %>% 
    dplyr::rename(GeneRatio=GeneRatio_mean, Count=Count_mean, count.terms=count.terms_sum)
  head(ck2.temp)
  
  #### update ck2 with previous info
  ck2.ID <- ck2.temp1 %>% 
    dplyr::bind_rows(ck2.temp) %>% 
    dplyr::group_by(ID, Description) %>% 
    dplyr::summarise(dplyr::across(c("GeneRatio", "Count"), list(mean=mean)), count.terms_sum = sum(count.terms), 
                     iterNum_count.terms = paste(iterNum_count.terms, collapse="")) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(GeneRatio=GeneRatio_mean, Count=Count_mean, count.terms=count.terms_sum)
  head(ck2.ID)
  
  write.table(ck2.ID, paste0(pathFeatureLabel,"/co-mod_records/K", length(archivos_),".csv"), sep=",", row.names = FALSE)
}