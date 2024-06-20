library("tidyverse")
library("ggpubr")

setwd("/media/ciccio/dati/scData/spaca6/GSE106487_data")

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


# Load necessary libraries
library(tidyverse) # For data manipulation and visualization
library(readr) # For reading the data file
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(forcats) # For working with factors

# Filter for SPACA6P and SPACA6P-AS genes
genes_of_interest <- combined_data_filt %>% filter(Gene %in% c("SPACA6P", "SPACA6P-AS")) %>% mutate(value=log2(value+1))

# Separate data by diagnosis for easier analysis
data_normal <- genes_of_interest %>% filter(Diagnosis == "NormalF")
data_disease <- genes_of_interest %>% filter(Diagnosis == "NOA") # Replace "Disease" with the actual disease name in your dataset

# Visualization with ggplot2
# Violin plot to show expression distribution
p <- ggplot(genes_of_interest, aes(x=Diagnosis, y=value, fill=Diagnosis)) +
  geom_violin() +
  facet_wrap(~Gene) +
  theme_minimal() +
  labs(title="Gene Expression by Diagnosis", y="TPM", x="") +
  scale_fill_manual(values = c("darkgrey", "steelblue")) + 
  stat_compare_means(method = "wilcox.test", paired = F,label = "p.signif", ref.group = "NormalF")
  
ggsave("./../R/expData.png", device = "png", width=5, height=4, dpi=300, plot = p, bg = "white")

# Box plot for central tendency and variability
p1 <- ggplot(genes_of_interest, aes(x=Diagnosis, y=value, fill=Diagnosis)) +
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





