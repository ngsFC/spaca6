library("tidyverse")
library("ggpubr")

setwd("/home/francescoc/Desktop/Project/spaca/GSE106487")

# List all .txt files in the directory
file_list <- list.files(pattern = "\\.txt$")

# Read all files into a list of data frames
data_list <- lapply(file_list, function(file) {
  read.table(file, header = TRUE, sep = "\t")  # Adjust 'sep' based on your file format
})

# Join all data frames on the 'Gene' column
combined_data <- Reduce(function(x, y) full_join(x, y, by = "Gene"), data_list)
combined_data_filt <- combined_data %>% 
  filter(grepl("SPACA6", Gene)) %>%
  pivot_longer(!Gene) %>%
  mutate(Diagnosis=case_when(grepl("N4", name) ~ "NOA",
                             !grepl("N4", name) ~ "NormalF")) %>%
  mutate(Diagnosis=factor(Diagnosis, levels=c("NormalF", "NOA")))


a <- ggplot(combined_data_filt %>%
              filter(Gene == "SPACA6P", value>0), aes(x=Diagnosis, y=value, color=Diagnosis)) +
  geom_violin(trim=FALSE) +
  geom_boxplot() + 
  #geom_jitter() + 
  scale_color_manual(values=c("grey", "steelblue")) +
  theme_minimal() +
  ylab("TPM") +
  #theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("SPACA6P gene") +
  #ylim(c(0,600)) +
  #scale_fill_manual(values = c()) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=90, hjust=1),
        text = element_text(size=15)) + 
  stat_compare_means(method = "t.test", paired = F,label = "p.signif", ref.group = "NormalF")
  


a <- ggplot(combined_data_filt %>%
              filter(Gene == "SPACA6P-AS", value >0), aes(x=Diagnosis, y=value, color=Diagnosis)) +
  geom_violin(trim=FALSE) +
  geom_boxplot() +
  #geom_jitter() + 
  scale_color_manual(values=c("grey", "steelblue")) +
  theme_minimal() +
  ylab("TPM") +
  #theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("SPACA6P-AS gene") +
  #ylim(c(0,100)) +
  #scale_fill_manual(values = c()) +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=90, hjust=1),
        text = element_text(size=15)) + 
  stat_compare_means(method = "t.test", paired = F,label = "p.signif", ref.group = "NormalF")


library(plotly)

anno <- combined_data_filt %>%
  filter(Gene == "SPACA6P",
         !name %in% c("GV_2n_4_48", "X4_GV_2N_2_14", "X5_1N_15", "X2JJ_25", "X1JC_15"))


pca <- combined_data %>% 
  filter(grepl("SPACA6", Gene)) %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Samples") %>%
  filter(!Samples %in% c("GV_2n_4_48", "X4_GV_2N_2_14", "X5_1N_15", "X2JJ_25", "X1JC_15")) %>% 
  column_to_rownames("Samples")
pca <- as.data.frame(pca)
pca["Color"] <- as.character(anno$Diagnosis)
ma <- pca[c(1:(length(pca)-1))]
PCma <- prcomp(ma, center=T,scale.=T)
PCma.gr <- data.frame(PCma$x, Color=pca$Color)
PCs <- round(PCma$sdev^2 / sum(PCma$sdev^2) * 100, 2)

PCs <- paste(colnames(PCma.gr), "(", paste(as.character(PCs), "%", ")", sep=""))
a <- ggplot(PCma.gr,aes(x=PC1, y=PC2, col=Color, label=rownames(PCma.gr))) + 
  xlab(PCs[1]) + 
  ylab(PCs[2]) +
  geom_point(size=2,alpha=1,aes(shape=Color)) +
  theme_classic() + 
  #xlim(-60,30) +
  #ylim(-60,60) +
  ggtitle("PCA all probes")
ggplotly(a)






# Load necessary libraries
library(tidyverse) # For data manipulation and visualization
library(readr) # For reading the data file
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(forcats) # For working with factors

# Filter for SPACA6P and SPACA6P-AS genes
genes_of_interest <- combined_data_filt %>% filter(Gene %in% c("SPACA6P", "SPACA6P-AS")) %>% mutate(value=log10(value+1))

# Separate data by diagnosis for easier analysis
data_normal <- genes_of_interest %>% filter(Diagnosis == "NormalF")
data_disease <- genes_of_interest %>% filter(Diagnosis == "NOA") # Replace "Disease" with the actual disease name in your dataset

# Visualization with ggplot2
# Violin plot to show expression distribution
ggplot(genes_of_interest, aes(x=Diagnosis, y=value, fill=Diagnosis)) +
  geom_violin() +
  facet_wrap(~Gene) +
  theme_minimal() +
  labs(title="Gene Expression by Diagnosis", y="TPM", x="") +
  scale_fill_manual(values = c("darkgrey", "steelblue")) + 
  stat_compare_means(method = "wilcox.test", paired = F,label = "p.signif", ref.group = "NormalF")
  

# Box plot for central tendency and variability
ggplot(genes_of_interest, aes(x=Diagnosis, y=value, fill=Diagnosis)) +
  geom_boxplot() +
  facet_wrap(~Gene) +
  theme_minimal() +
  labs(title="Gene Expression by Diagnosis", y="TPM", x="") +
  scale_fill_manual(values = c("darkgrey", "steelblue"))


# Statistical Testing
# Assuming the data is not normally distributed, use Mann-Whitney U test
# Mann-Whitney U test for SPACA6P
spaca6p_values_normal <- data_normal %>% filter(Gene == "SPACA6P") %>% pull(value)
spaca6p_values_disease <- data_disease %>% filter(Gene == "SPACA6P") %>% pull(value)
wilcox.test(spaca6p_values_normal, spaca6p_values_disease)

# Repeat the Mann-Whitney U test for SPACA6P-AS similarly
spaca6p_values_normal <- data_normal %>% filter(Gene == "SPACA6P-AS") %>% pull(value)
spaca6p_values_disease <- data_disease %>% filter(Gene == "SPACA6P-AS") %>% pull(value)
wilcox.test(spaca6p_values_normal, spaca6p_values_disease)

# Note: You might need to adjust the script based on the exact structure of your dataset and the specific names of the diagnoses.





