library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)

pathFeatureLabel <- "C:/Users/da.salazarb/Downloads/Nueva_carpeta/kmmjnmf/pathFeatureLabel"

real_data <- read.csv(paste0(pathFeatureLabel, "/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_['cnv', 'drug', 'mirna', 'mrna'].csv"), header = TRUE)
real_data$X <- NULL
real_data <- unique(real_data)
head(real_data)

real_data <- real_data %>% 
  dplyr::filter(last_stop_control < 1e-3) %>% 
  dplyr::select(K, r1, r2, d1, L1, L2, o1, o2, last_stop_control, contains("rho"), contains("AUC"), contains("WeightedAverage"), BioRelevanceScore) %>% 
  dplyr::mutate(K=as.character(K)) %>% 
  dplyr::arrange(K)

head(real_data)
dim(real_data)

## ----------------------------------------------##
# order_data_kernel <- real_data %>% 
#   dplyr::select(K, contains("rho"), AUC_H_average) %>% 
#   dplyr::mutate(prom_rho = rowMeans(dplyr::select(.,starts_with("rho")),na.rm = TRUE), prom_auc = AUC_H_average) %>% 
#   dplyr::mutate(prom_pond=rowMeans(dplyr::select(., contains("prom")))) %>% 
#   dplyr::select(-starts_with("rho_"), -starts_with("AUC_"), -prom_rho, -prom_auc) %>% 
#   dplyr::group_by(K) %>%
#   dplyr::summarize(mean_K_prom = mean(prom_pond, na.rm = TRUE)) %>% 
#   dplyr::arrange(mean_K_prom) %>% 
#   dplyr::mutate(K=paste0("K: ",K)) %>% 
#   dplyr::arrange(K)
# order_data_kernel

# _______________________________________________ #
#### ........... rho standard jNMF ........... ####
# ______________________________________________ #
### rho
rho_data_ <- real_data %>% 
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, contains("rho"), -contains("WeightedAverage")) %>% 
  dplyr::mutate(sum_rho=rowSums(dplyr::select(.,starts_with("rho")),na.rm = TRUE)) %>%
  dplyr::group_by(K) %>%
  dplyr::filter(sum_rho >= 2) %>%
  dplyr::ungroup(K) %>%
  dplyr::arrange(desc(sum_rho))

rho_data_ <- rho_data_[match(unique(rho_data_$K), rho_data_$K),]

rho_data_ <- rho_data_ %>%
  tidyr::gather(key="rho", value = "rho_value", -K, -r1, -r2, -L1, -L2, -d1, -o1, -o2, -sum_rho) %>%
  dplyr::mutate(K = as.factor(K)) %>%
  dplyr::mutate(rho=gsub("[0-9]*_", "", gsub("rho", "", rho))) %>% 
  dplyr::mutate(K=paste0(K,"_", rho)) %>% 
  dplyr::filter(!is.na(rho_value))

data_profile <- as.data.frame(c("mirna", "cnv", "protein", "mrna", "drug"),
                              c("miRNA", "CNV", "Protein", "mRNA", "Drug"))
colnames(data_profile)[1] <- "profile"

rho_data_ <- rho_data_ %>%
  tidyr::separate(K, c("K", "rho"), "_") %>%
  dplyr::mutate(x=paste0("K: ", K, " - r1: ", r1, " - r2: ", r2, " - L1: ", L1, " - L2: ", L2, " - d1: ", d1, " - o1: ", o1, " - o2: ", o2)) %>%
  dplyr::mutate(x=gsub(" \\- L2", "\n \\- L2", x)) %>%
  dplyr::mutate(x=factor(x, levels = unique(x[order(sum_rho)]))) %>%
  dplyr::mutate(rho_profile=sapply(rho, function(x) row.names(data_profile)[grep(x,data_profile$profile)] ) )

rho_p_jnmf <- ggplot(rho_data_, aes(x=x, y=rho_value, group=rho_profile)) +
  geom_point() + geom_line(aes(linetype=rho_profile)) +
  labs(title="", x="", y="Cophenetic coefficient", linetype="Input matrix") +
  theme(text = element_text(size=18), axis.text.x = element_text(size=10,angle = 90))
rho_p_jnmf

# _______________________________________________ #
#### ........... AUC standard jNMF ........... ####
# ______________________________________________ #
### auc
auc_data_ <- real_data %>%
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, contains("auc"), -contains("WeightedAverage"), -AUC_H_average)
  # dplyr::mutate(sum_auc=rowSums(dplyr::select(.,starts_with("auc")))) %>%
  # dplyr::group_by(K) %>%
  # dplyr::filter(sum_auc > 2) %>%
  # dplyr::ungroup(K) %>%
  # dplyr::arrange(desc(AUC_H_average))

auc_data_ <- auc_data_[match(unique(auc_data_$K), auc_data_$K),]

auc_data_ <- auc_data_ %>%
  # dplyr::filter(sum_rho > quantile(sum_rho,prob=.90)) %>%
  tidyr::gather(key="auc", value = "auc_value", -K, -r1, -r2, -L1, -L2, -d1, -o1, -o2) %>%
  dplyr::mutate(K = as.factor(K)) %>%
  dplyr::mutate(auc=gsub("[0-9]*_", "", gsub("AUC_H_", "", auc))) %>%
  # dplyr::mutate(r1=round(r1,3), L1=round(L1,3), L2=round(L2,3), d1=round(d1,3)) %>%
  dplyr::mutate(K=paste0(K,"_", auc))

data_profile <- as.data.frame(c("mirna", "cnv", "protein", "mrna", "drug"),
                              c("miRNA", "CNV", "Protein", "mRNA", "Drug"))
colnames(data_profile)[1] <- "profile"

auc_data_ <- auc_data_ %>%
  tidyr::separate(K, c("K", "auc"), "_") %>%
  dplyr::mutate(x=paste0("K: ", K, " - r1: ", r1, " - r2: ", r2, " - L1: ", L1, " - L2: ", L2, " - d1: ", d1, " - o1: ", o1, " - o2: ", o2)) %>%
  dplyr::mutate(x=gsub(" \\- L2", "\n \\- L2", x)) %>%
  dplyr::mutate(x=factor(x, levels = unique(x[order(order_data_kernel$K)]))) %>%
  dplyr::mutate(auc_profile=sapply(auc, function(x) row.names(data_profile)[grep(x,data_profile$profile)] ) )

auc_p_jnmf <- ggplot(auc_data_, aes(x=x, y=auc_value, group=auc_profile)) +
  geom_point() + geom_line(aes(linetype=auc_profile)) +
  labs(title="", x="", y="AUC", linetype="Input matrix") +
  theme(text = element_text(size=18), axis.text.x = element_text(size=10,angle = 90))
auc_p_jnmf

# ______________________________________ #
#### ........... BioScore ........... ####
# ______________________________________ #
### auc
bio_data_ <- real_data %>%
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, BioRelevanceScore) %>% 
  dplyr::arrange(desc(BioRelevanceScore))

bio_data_ <- bio_data_[match(unique(bio_data_$K), bio_data_$K),]

bio_data_ <- bio_data_ %>%
  tidyr::gather(key="bio", value = "bio_value", -K, -r1, -r2, -L1, -L2, -d1, -o1, -o2) %>%
  dplyr::mutate(K = as.factor(K))

bio_data_ <- bio_data_ %>%
  # tidyr::separate(K, c("K", "auc"), "_") %>%
  dplyr::mutate(x=paste0("K: ", K, " - r1: ", r1, " - r2: ", r2, " - L1: ", L1, " - L2: ", L2, " - d1: ", d1, " - o1: ", o1, " - o2: ", o2)) %>%
  dplyr::mutate(x=gsub(" \\- L2", "\n \\- L2", x)) %>% 
  dplyr::mutate(x=factor(x, levels = unique(x[order(order_data_kernel$K)]))) # %>%
  # dplyr::mutate(auc_profile=sapply(auc, function(x) row.names(data_profile)[grep(x,data_profile$profile)] ) )

bio_p_jnmf <- ggplot(bio_data_, aes(x=x, y=bio_value, group=1)) +
  geom_line() + geom_point() + 
  labs(title="", x="", y="BioScore", linetype="Input matrix") +
  theme(text = element_text(size=18), axis.text.x = element_text(size=10,angle = 90))
bio_p_jnmf

# __________________________________________________ #
#### ........... individual plots LGG ........... ####
# __________________________________________________ #
pdf(file = paste0(pathFeatureLabel, "/co-mod_plots/SuppFigure5a.pdf"), width = 10, height = 10)
rho_p_jnmf
dev.off()

pdf(file = paste0(pathFeatureLabel, "/co-mod_plots/SuppFigure5b.pdf"), width = 10, height = 10)
auc_p_jnmf
dev.off()