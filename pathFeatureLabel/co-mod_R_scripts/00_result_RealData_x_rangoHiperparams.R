library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(ggalluvial)

pathFeatureLabel <- "C:/Users/da.salazarb/Downloads/Nueva_carpeta/kmmjnmf/pathFeatureLabel"

real_data <- read.csv(paste0(pathFeatureLabel, "/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_['cnv', 'drug', 'mirna', 'mrna'].csv"), header = TRUE)

real_data <- real_data %>% 
  dplyr::select(-X) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(last_stop_control < 1e-3) %>% 
  dplyr::select(K, r1, r2, d1, L1, L2, o1, o2, last_stop_control, contains("rho"), contains("AUC"), contains("WeightedAverage"), BioRelevanceScore) %>% 
  dplyr::mutate(K=as.character(K)) %>% 
  dplyr::arrange(K)

head(real_data)
dim(real_data)

# _______________________________________________ #
#### ........... rho standard jNMF ........... ####
# ______________________________________________ #
### rho
rho_data_ <- real_data %>% 
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, contains("rho"), -contains("WeightedAverage")) %>% 
  dplyr::mutate(sum_rho=rowSums(dplyr::select(.,starts_with("rho")),na.rm = TRUE)) %>%
  # dplyr::group_by(K) %>%
  # dplyr::filter(sum_rho >= 2) %>%
  # dplyr::ungroup(K) %>%
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
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, contains("auc"), -contains("WeightedAverage"), -AUC_H_average) %>% 
  dplyr::mutate(sum_auc=rowSums(dplyr::select(.,starts_with("auc")),na.rm = TRUE)) %>% 
  dplyr::arrange(desc(sum_auc))

auc_data_ <- auc_data_[match(unique(auc_data_$K), auc_data_$K),]

auc_data_ <- auc_data_ %>%
  # dplyr::filter(sum_rho > quantile(sum_rho,prob=.90)) %>%
  tidyr::gather(key="auc", value = "auc_value", -K, -r1, -r2, -L1, -L2, -d1, -o1, -o2, -sum_auc) %>%
  dplyr::mutate(K=as.numeric(K)) %>% 
  dplyr::arrange(K) %>% 
  dplyr::mutate(K = as.factor(K)) %>%
  dplyr::mutate(auc=gsub("[0-9]*_", "", gsub("AUC_H_", "", auc))) %>%
  dplyr::mutate(K=paste0(K,"_", auc))

data_profile <- as.data.frame(c("mirna", "cnv", "protein", "mrna", "drug"),
                              c("miRNA", "CNV", "Protein", "mRNA", "Drug"))
colnames(data_profile)[1] <- "profile"

auc_data_ <- auc_data_ %>%
  tidyr::separate(K, c("K", "auc"), "_") %>%
  dplyr::mutate(x=paste0("K: ", K, " - r1: ", r1, " - r2: ", r2, " - L1: ", L1, " - L2: ", L2, " - d1: ", d1, " - o1: ", o1, " - o2: ", o2)) %>%
  dplyr::mutate(x=gsub(" \\- L2", "\n \\- L2", x)) %>% 
  dplyr::mutate(K=as.numeric(K)) %>% 
  dplyr::arrange(K) %>% 
  dplyr::mutate(x=factor(x, levels = unique(x[order(K)]))) %>% 
  dplyr::mutate(auc_profile=sapply(auc, function(x) row.names(data_profile)[grep(x,data_profile$profile)] ) )
head(auc_data_)

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
  dplyr::select(K, r1, r2, L1, L2, d1, o1, o2, BioRelevanceScore)

bio_data_ <- bio_data_[match(unique(bio_data_$K), bio_data_$K),]

bio_data_ <- bio_data_ %>%
  tidyr::gather(key="bio", value = "bio_value", -K, -r1, -r2, -L1, -L2, -d1, -o1, -o2) %>%
  dplyr::mutate(K = as.factor(K))

bio_data_ <- bio_data_ %>%
  dplyr::mutate(x=paste0("K: ", K, " - r1: ", r1, " - r2: ", r2, " - L1: ", L1, " - L2: ", L2, " - d1: ", d1, " - o1: ", o1, " - o2: ", o2)) %>%
  dplyr::mutate(x=gsub(" \\- L2", "\n \\- L2", x)) %>% 
  dplyr::mutate(x=factor(x, levels = unique(x[order(desc(bio_value))])))

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