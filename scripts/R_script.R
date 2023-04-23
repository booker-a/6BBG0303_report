setwd("/Users/anthonybooker")

### IMPORTANT: GIVE THE NMD_METRIC A PROPER NAME IN THE GRAPHS


# PACKAGE CAR PARK ----
library(sp)
library(SeuratObject)
library(BiocManager)
library(broom)
library(magrittr)
library(Seurat)
library(rtracklayer)
library(tidyverse)
library(enrichplot)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(corrplot)
library(corrr)
library(stats)
library(clusterProfiler) ### NOTE: masks other package functions - unload when not needed
library(org.Mm.eg.db) ### NOTE: masks other package functions - unload when not needed
library(car)
library(rstatix)
library(caret)
library(patchwork)
library(rcartocolor)

detach(package:clusterProfiler, unload=TRUE)



# TASK 1 - Pairing up cell and nuclei data points ----

# IMPORTING SEURAT OBJECT
seurat <- readRDS("/Users/anthonybooker/Desktop/fp01_asnmd/data/Bakken_nuclei_cell_seurat.rds")


## Viewing the seurat object
View(seurat)


## Locating the SNN graph within the 'graphs' slot
seurat@graphs$integrated_snn


## Locating the metadata slot, so that we can extract WholeCell/Nuclei data from the 'orig.ident' variable/column
seurat@meta.data


## Extracting data points for "Wholecell" entries, and assigning to a new vector
wholecell_names <- filter(seurat@meta.data, orig.ident == "WholeCell") %>%
  rownames_to_column("cell") %>%
  pull(cell)

View(seurat@meta.data)
## Extracting data points for "Nuclei" entries, and assigning to new a vector
nuclei_names <- filter(seurat@meta.data, orig.ident == "Nuclei") %>%
  rownames_to_column("cell") %>%
  pull(cell)


## Creating a new matrix with data points that match with whole_cell_names and nuclei_names
subsetted_snn <- as.matrix(seurat@graphs$integrated_snn[wholecell_names, nuclei_names])


## Converting from matrix to data frame, and pivoting from a wide to long format
subsetted_snn <- subsetted_snn %>%
  as.data.frame() %>% 
  rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "nuclei", values_to = "correlation") # or pivot_longer(-cell)

View(subsetted_snn)

## Pairing up the cells and nuclei, and keeping unique/distinct entries (repeats or duplicates are removed), and assigning to new variable
cell_nuclei_pair <- subsetted_snn %>%
  group_by(cell) %>%
  arrange(desc(correlation)) %>%
  distinct(cell, .keep_all = T) %>% View()
  ungroup() %>% View()

View(cell_nuclei_pair)

saveRDS(cell_nuclei_pair, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/cell_nuclei_pair.rds")



# TASK 2 - Identifying NMD+/- genes from GTF dataset ----

## Importing the desired GTF file from GENCODE
GTF <- readGFF("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz")


## Extracting gene names and gene IDs from the GTF file that are NMD positive
NMD_plus_GTF <- GTF %>% 
  filter(transcript_type=="nonsense_mediated_decay") %>% 
  group_by(gene_name) %>%
  select(gene_id, gene_name) %>%
  unique()

View(NMD_plus_GTF)

saveRDS(NMD_plus_GTF, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/NMD_plus_GTF.rds")


## Extracting gene names and gene IDs from the GTF file that are NMD negative
NMD_minus_GTF <- GTF %>% 
  filter(!gene_id %in% NMD_plus_GTF$gene_id) %>%
  filter(gene_type=="protein_coding") %>% 
  filter(type == "exon") %>%
  group_by(gene_id, gene_name, transcript_id) %>%
  summarise(num_of_exons = n()) %>% 
  ungroup() %>% 
  filter(num_of_exons >= 3) %>%
  group_by(gene_id, gene_name) %>% 
  filter(n() > 2) %>% View() ## filtering for 2 or more entries
  distinct(gene_id, gene_name)

View(NMD_minus_GTF)

saveRDS(NMD_minus_GTF, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/NMD_minus_GTF.rds")



# TASK 3 - Calculating the log2 fold change difference between cell and nuclei gene counts ----

## Combining the grouped NMD plus and minus data frames and filtering for genes in seurat object
NMD_combined <- NMD_plus_GTF %>%
  mutate(type = "NMD_plus") %>%
  bind_rows(NMD_minus_GTF %>% mutate(type = "NMD_minus")) %>%
  filter(gene_name %in% rownames(seurat)) %>%
  ungroup() %>% View()

saveRDS(NMD_combined, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/NMD_combined.rds")


## Creating a matrix for the counts of genes for nuclei
nuclei_genecounts <- as.matrix(seurat@assays$RNA@data)[NMD_combined$gene_name, cell_nuclei_pair$nuclei]

View(nuclei_genecounts)

## Creating a matrix for the counts of genes for cells
cell_genecounts <- as.matrix(seurat@assays$RNA@data)[NMD_combined$gene_name, cell_nuclei_pair$cell]


## Separating cell and nuclei data points in the same column by using the deliminator comma
colnames(nuclei_genecounts) <- paste0(cell_nuclei_pair$cell, ",", cell_nuclei_pair$nuclei)


## Calculating the log2 fold change difference in gene counts between cell and nuclei (more positive value means more expression in nuclei, vice versa for cells)
logFC_genecounts <- nuclei_genecounts - cell_genecounts

View(logFC_genecounts)



# TASK 4 - Adding the log2 fold change values into the cell and nuclei paired data frame ----

## Pivoting logFC_genecounts from wide format into long format
logFC_genecounts_long <- logFC_genecounts %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  pivot_longer(-gene_name, names_to = "pairs", values_to = "logFC")

View(logFC_genecounts_long)

seurat@meta.data %>% View()

## Separating cell and nuclei data points by the deliminator ","
sep_logFC_genecounts_long <- separate(logFC_genecounts_long, pairs, into = c("cell", "nuclei"), sep = ",", remove = TRUE)

View(sep_logFC_genecounts_long)


## Joining to the cell_nuclei_pair data frame
cell_nuclei_pair <- left_join(sep_logFC_genecounts_long, cell_nuclei_pair)

View(cell_nuclei_pair)



# TASK 5 - Testing the accuracy of nuclei-cell pairing ----

## TASK 5a - Using 'tasic.major.class' ----

## Create df of cell names and tasic.major.class (cell_class was not used as it contained some N/A values which may interfere with the results downstream)
names_major_class <- seurat@meta.data %>% 
  select(tasic.major.class) %>% 
  rownames_to_column("cells") %>% View()


## Creating a new cell-nuclei pair df to work with for testing 
testing_major_cell_nuclei_pair <- subsetted_snn %>%
  group_by(cell) %>%
  arrange(desc(correlation)) %>%
  distinct(cell, .keep_all = T) %>%
  ungroup() %>%
  select(-correlation) %>%
  left_join(names_major_class, by = c("cell"="cells")) %>%
  dplyr::rename(cell_ID = tasic.major.class) %>% ##dplyr package specified as otherwise there are conflicts with the rename() function
  left_join(names_major_class, by = c("nuclei"="cells")) %>%
  dplyr::rename(nuclei_ID = tasic.major.class) 

View(testing_major_cell_nuclei_pair)

(testing_major_cell_nuclei_pair)

## Creating a new column to assess match of tasic.major.class data between cell-nuclei pairs
testing_major_cell_nuclei_pair <- testing_major_cell_nuclei_pair %>% 
  dplyr::mutate(Match = (cell_ID == nuclei_ID))

View(testing_major_cell_nuclei_pair)

seurat_cell@meta.data %>% View()

## Calculating the amount of T/F
summary(testing_major_cell_nuclei_pair$Match)







## ## TASK 5b - Using 'tasic.sub.class' ----

## Creating df of cell names and tasic.major.class
names_sub_class <- seurat@meta.data %>% 
  select(tasic.sub.class) %>% 
  rownames_to_column("cells") 


## Creating a new cell-nuclei pair df to work with for testing tasic.sub.class
testing_sub_cell_nuclei_pair <- subsetted_snn %>%
  group_by(cell) %>%
  arrange(desc(correlation)) %>%
  distinct(cell, .keep_all = T) %>%
  ungroup() %>%
  select(-correlation) %>%
  left_join(names_sub_class, by = c("cell"="cells")) %>%
  dplyr::rename(cell_ID = tasic.sub.class) %>%
  left_join(names_sub_class, by = c("nuclei"="cells")) %>%
  dplyr::rename(nuclei_ID = tasic.sub.class) 

View(testing_sub_cell_nuclei_pair)


## Creating a new column to assess match of tasic.sub.class data between cell-nuclei pairs
testing_sub_cell_nuclei_pair <- testing_sub_cell_nuclei_pair %>% 
  dplyr::mutate(Match = (cell_ID == nuclei_ID))

View(testing_sub_cell_nuclei_pair)


## Calculating the amount of T/F
summary(testing_sub_cell_nuclei_pair$Match)




## TASK 5c - Using 'cell_class' ----

## Creating df of cell_names and cell_class
names_cell_class <- seurat@meta.data %>% 
  select(cell_class) %>% 
  rownames_to_column("cells") 

View(names_cell_class)


## Creating a new cell-nuclei pair df to work with for testing 
testing_cell_class_cell_nuclei_pair <- subsetted_snn %>%
  group_by(cell) %>%
  arrange(desc(correlation)) %>%
  distinct(cell, .keep_all = T) %>%
  ungroup() %>%
  select(-correlation) %>%
  left_join(names_cell_class, by = c("cell"="cells")) %>%
  dplyr::rename(cell_ID = cell_class) %>% ##dplyr package specified as otherwise there are conflicts with the rename() function
  left_join(names_cell_class, by = c("nuclei"="cells")) %>%
  dplyr::rename(nuclei_ID = cell_class) 

View(testing_cell_class_cell_nuclei_pair)


## Creating a new column to assess match of tasic.major.class data between cell-nuclei pairs
testing_cell_class_cell_nuclei_pair <- testing_cell_class_cell_nuclei_pair %>% 
  dplyr::mutate(Match = (cell_ID == nuclei_ID))

View(testing_cell_class_cell_nuclei_pair)


## Calculating the amount of T/F
summary(testing_cell_class_cell_nuclei_pair$Match)

View(seurat_cell@meta.data)



## TASK 5d - making a bar graph for the report ----


Figure6Adata <- data.frame(
  "FALSE" = c(1,14),
  "TRUE" = c(426,413),
  row.names = c("tasic.major.class", "tasic.sub.class")
)

colnames(Figure6Adata) <- c("FALSE", "TRUE")

Figure6Adata <- rownames_to_column(Figure6Adata, "identity")

Figure6Adata <- Figure6Adata %>% 
  pivot_longer(-identity, names_to = "match", values_to = "n") %>% View()

View(Figure6Adata)

ggplot(Figure6Adata, aes(x = identity, y = n), fill = as.factor(match)) +




# TASK 6 - Calculating the expected proportions cell identity classes ----

## TASK 6a - for 'tasic.major.class' ----



## Displaying the proportions of Nuclei/WholeCell entries for Excitatory and Inhibitory cell classes
table(seurat@meta.data$orig.ident, seurat@meta.data$tasic.major.class)


## Making a new df for expected tasic_major_class
expected_tasic_major_class <- seurat@meta.data %>%
  group_by(orig.ident, tasic.major.class) %>%
  tally() %>%
  pivot_wider(names_from = tasic.major.class, values_from = orig.ident)


## Summary of the new df
sum(expected_tasic_major_class$n)

View(expected_tasic_major_class)

View(seurat_cell@meta.data)
## Calculating the percentage of matched expected entries
### p(matched) = p(Exc | WholeCell) * p(Exc | Nuclei) + p(Inh | WholeCell) * p(Inh | Nuclei)
prcnt_expected_tasic_major_class_match <- 
  ( (expected_tasic_major_class[3,1] / sum(expected_tasic_major_class[3:4,1])) * #p(Exc | WholeCell)
      (expected_tasic_major_class[1,1] / sum(expected_tasic_major_class[1:2,1]))) + #p(Exc | Nuclei)
  ( (expected_tasic_major_class[4,1] / sum(expected_tasic_major_class[3:4,1])) * #p(Inh | WholeCell)
      (expected_tasic_major_class[2,1] / sum(expected_tasic_major_class[1:2,1])) ) #p(Inh | Nuclei)

prcnt_expected_tasic_major_class_match




## TASK 6b - for 'tasic.sub.class' ----

## Displaying the proportions of Nuclei/WholeCell entries for all sub classes
expected_tasic_sub_class <- table(seurat@meta.data$orig.ident, seurat@meta.data$tasic.sub.class) %>% 
  as.data.frame() %>%
  #tally(Freq) %>% ## ADDS UP TO 890 CELLS SO THIS IS GOOD
  pivot_wider(names_from = Var2, values_from = Freq) %>% # DO I DO THIS?select(-L5_Chrna6, -L6b) %>% ## These subclasses are removed from the expected calculations as no values were recorded in WholeCell data)
  dplyr::rename("identity" = "Var1")

View(expected_tasic_sub_class)


#### p(matched) = SUM p(subclassi | WholeCell) * p(subclassi | Nuclei) 
prcnt_expected_tasic_sub_class_match <-
  ( (expected_tasic_sub_class[2,2] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L2/3 | WholeCell)
      (expected_tasic_sub_class[1,2] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L2/3 | Nuclei)
  ( (expected_tasic_sub_class[2,3] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L4 | WholeCell)
      (expected_tasic_sub_class[1,3] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L4 | Nuclei)
  ( (expected_tasic_sub_class[2,4] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L5_Chrna6 | WholeCell)
     (expected_tasic_sub_class[1,4] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L5_Chrna6 | Nuclei)
  ( (expected_tasic_sub_class[2,5] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L5a | WholeCell)
      (expected_tasic_sub_class[1,5] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L5a | Nuclei)
  ( (expected_tasic_sub_class[2,6] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L5a2 | WholeCell)
      (expected_tasic_sub_class[1,6] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L5a2 | Nuclei)
  ( (expected_tasic_sub_class[2,7] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L5b | WholeCell)
      (expected_tasic_sub_class[1,7] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L5b | Nuclei)
  ( (expected_tasic_sub_class[2,8] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L6a1 | WholeCell)
      (expected_tasic_sub_class[1,8] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L6a1 | Nuclei)
  ( (expected_tasic_sub_class[2,9] / sum(expected_tasic_sub_class[2, 2:14])) * #p(L6a2 | WholeCell)
      (expected_tasic_sub_class[1,9] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(L6a2 | Nuclei)
  ( (expected_tasic_sub_class[2,10] / sum(expected_tasic_sub_class[2, 2:14])) * #p(Pvalb | WholeCell)
      (expected_tasic_sub_class[1,10] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(Pvalb | Nuclei)
  ( (expected_tasic_sub_class[2,11] / sum(expected_tasic_sub_class[2, 2:14])) * #p(Sst_Cbln4 | WholeCell)
      (expected_tasic_sub_class[1,11] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(Sst_Cbln4 | Nuclei)
  ( (expected_tasic_sub_class[2,12] / sum(expected_tasic_sub_class[2, 2:14])) * #p(Sst_Nrf2f | WholeCell)
      (expected_tasic_sub_class[1,12] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(Sst_Nrf2f | Nuclei)
  ( (expected_tasic_sub_class[2,13] / sum(expected_tasic_sub_class[2, 2:14])) * #p(Sst_Cbln4 | WholeCell)
      (expected_tasic_sub_class[1,13] / sum(expected_tasic_sub_class[1, 2:14]))) + #p(Sst_Cbln4 | Nuclei)
  ( (expected_tasic_sub_class[2,14] / sum(expected_tasic_sub_class[2, 2:14])) * #p(Vip | WholeCell)
      (expected_tasic_sub_class[1,14] / sum(expected_tasic_sub_class[1, 2:14])) ) #p(Vip | Nuclei)
  
prcnt_expected_tasic_sub_class_match
  
# TASK 7 - Calculating Fisher's test p.values for each cell identity classes ----

##TASK 7a - for 'tasic.major.class' ----

## As 78% of 427 cells is 333, the expected number of matched pairs is 333, and 94 unmatched. Comparing against observed matches:

summary(testing_major_cell_nuclei_pair$Match) 

### Discussion:
### as carried out before, 426/427 match, while 1/427 do not match.


## Making a new df to carry out a Fishers test on tasic_major_class
Fishers_major <- data.frame(
  "Matched" = c(426, 333),
  "Unmatched" = c(1, 94),
  row.names = c("Observed", "Expected"),
  stringsAsFactors = FALSE
)
colnames(Fishers_major) <- c("Matched", "Unmatched")

Fishers_major


## Carrying out a Fishers test
fisher.test(Fishers_major, alternative = "two.sided") ## added alternative to double check, two-sided is actually the default


### Discussion:
### As the p.value is < 2.2e-16, this suggests that the results are highly significant



## TASK 7b - for 'tasic.sub.class' ----

## As 21% of 427 cells is 89, the expected number of matched pairs is 89, and 338 unmatched. Comparing against observed matches:

## Making a new df to carry out a Fishers test on tasic_major_class
Fishers_sub <- data.frame(
  "Matched" = c(413, 89),
  "Unmatched" = c(14,338),
  row.names = c("Observed", "Expected"),
  stringsAsFactors = FALSE
)
colnames(Fishers_sub) <- c("Matched", "Unmatched")

Fishers_sub

## Carrying out a Fishers test
fisher.test(Fishers_sub, alternative = "two.sided") ## added alternative to double check, two-sided is actually the default

### Discussion:
### As the p.value is < 2.2e-16, this suggests that the results are highly significant, + odds ratio 111




# Special Task: FIGURE 6 PLOTS ----
Figure6B <- data.frame(
  tasic.major.class = c("Observed", "Expected"),
  Matched = c(426, 333),
  Unmatched = c(1, 94)
)

View(Figure6B)

Figure6B <- gather(Figure6B, key = "type", value = "count", -tasic.major.class)

View(Figure6B)

ggplot(Figure6B, aes(x = tasic.major.class, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  guides(fill=guide_legend(title="")) +
  theme_pubr() +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "'tasic.major.class'") +
  geom_text(aes(label = count), position = position_fill(vjust=1), vjust=-0.25)


Figure6A <- Fishers_sub %>% 
  as.data.frame() %>% 
  rownames_to_column("tasic.sub.class")

View(Figure6A)

Figure6A <- gather(Figure6A, key = "type", value = "count", -tasic.sub.class)

View(Figure6A)

ggplot(Figure6A, aes(x = tasic.sub.class, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  guides(fill=guide_legend(title="")) +
  theme_pubr() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "'tasic.sub.class'") +
  geom_text(aes(label = count), position = position_fill(vjust = 1), vjust=-0.25)


# TASK 8 - Plotting some of the genes to check the skew of whether its more or less expresses in cell/nuclei ----
## (more +ve = more nuclei expression, more -ve = more cell expression)

## Dataframes used to select genes for plotting:
View(NMD_plus_GTF)
View(NMD_minus_GTF)
View(NMD_combined)


## Creating a new df for plotting
cell_nuclei_pair_plot <- cell_nuclei_pair %>%
  as.data.frame() %>% 
  select(gene_name, logFC)

View(cell_nuclei_pair_plot)


## Plotting the logFC distribution for different genes
cell_nuclei_pair_plot %>% 
  filter(gene_name %in% c("Rab11a", "Wsb1")) %>% ## an NMD-plus and NMD-minus gene
  ggplot(aes(x = gene_name, y = logFC, fill = gene_name)) +
  geom_jitter() +
  geom_violin() +
  stat_summary(fun = "mean",
               geom = "crossbar",
               colour = "red", show_guide=FALSE) +
  theme(legend.position = "top") +
  theme_pubr() +
  labs(x = "gene name",
       y = "log2FC") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                    name="gene type",
                    breaks=c("Rab11a", "Wsb1"),
                    labels=c("Control", "NMD"))



# TASK 9 - Calculating non weighted & untransformed NMD metric ---- 

## Importing dataset to be used for calculating the NON WEIGHTED NMD METRIC 
weighting_dataset <- read_tsv("Desktop/fp01_asnmd/data/DEGs_CHX_Day12.tsv")

View(weighting_dataset)


## Subsetting all NMD+ genes (padj < 0.05 & log2FoldChange >= log2(2))
candidate_dataset_plus_NMD <- weighting_dataset %>%
  filter(gene %in% NMD_plus_GTF$gene_name) %>% 
  filter(padj < 0.05 & log2FoldChange >= log2(2)) ## selecting significant genes

View(candidate_dataset_plus_NMD)


## Condensing NMDplus
cond_NMDplus <- cell_nuclei_pair %>%
  filter(gene_name %in% candidate_dataset_plus_NMD$gene) %>%
  mutate(linear_FC = 2^logFC) %>%
  group_by(cell) %>%
  summarise(average_NW_FC_NMD_PLUS=mean(linear_FC))

View(cond_NMDplus)


## Subsetting all NMD- genes - these genes should not be asssociated with NMD

### Discussion:
### making a lapply() loop  to evaluate which random sample of NMD- genes produces the best average_NW_above_base curve
### last attempt at NMD_metric was good - but how can we be sure that was the best random sample that could be achieved?
### hence, let's sample from 5000 seeds to decide which could result in the best NMD outcome


### TASK 9a - Optimising the NMD metric: choosing the best sample of NMD- based on 5000 seed iterations ----


#### Making a new df, filtering for genes with certain parameters
candidate_dataset_minus_NMD <- weighting_dataset %>%
  filter(gene %in% NMD_minus_GTF$gene_name) %>% 
  filter(padj >= 0.05, log2FoldChange < 0.5, log2FoldChange > -0.5)

View(candidate_dataset_minus_NMD)

#### Retrieving average_NW_FC_NMD_MINUS for each cell per seed (5000 iterations)
testing_NMD_metrics <- do.call(bind_rows, lapply(1000:4999, function(x){
  col_name <- paste0("seed",x)
  set.seed(x)
  tmp_sample <- candidate_dataset_minus_NMD %>% 
    slice_sample(n = nrow(candidate_dataset_plus_NMD), replace = F)
  cell_nuclei_pair %>%
    filter(gene_name %in% tmp_sample$gene) %>% 
    mutate(linear_FC = 2^logFC) %>%
    group_by(cell) %>%
    summarise(average_NW_FC_NMD_MINUS=mean(linear_FC)) %>% 
    mutate(seed = x) %>%
    select(cell, average_NW_FC_NMD_MINUS, seed) 
}))

View(testing_NMD_metrics)


#### Calculating the average_NW_above_base
testing_NMD_metrics <- left_join(testing_NMD_metrics, cond_NMDplus, by = "cell")

testing_NMD_metrics <- testing_NMD_metrics %>% 
  mutate(average_NW_above_base = (average_NW_FC_NMD_PLUS / average_NW_FC_NMD_MINUS))

View(testing_NMD_metrics)


## TASK 9b - Selecting the best NMD_metric ----

### Filtering for cells that have average_NW_above_base above 1 + evaluating how many cells per seed
testing_NMD_metrics %>% 
  filter(average_NW_above_base > 1) %>% # a value of less than one means that the level of NMD is below baseline (i.e no NMD minus)
  group_by(seed) %>% 
  tally() %>% View()


### Conducting Shapiro-Wilks test
filtered_testing_NMD_metrics_shapiro <- testing_NMD_metrics %>%
  filter(average_NW_above_base > 1) %>% 
  select(cell, seed, average_NW_above_base) %>% 
  group_by(seed) %>% 
  shapiro_test(average_NW_above_base)

View(filtered_testing_NMD_metrics_shapiro)

shapiro

### Filtering out for the top 10 seeds that have the highest p value (i.e are the most normally distributed)
filtered_testing_NMD_metrics_shapiro %>% 
  top_n(p, n = 10) %>%
  arrange(desc(p)) %>% View()


### Evaluating how many cells there are for our desired seed (2065)
testing_NMD_metrics %>% 
  filter(average_NW_above_base > 1) %>%
  group_by(seed) %>% View()
  filter(seed == "2065") %>%
  tally()

### Discussion:
### This gives 423 cell, and a high degree of normal distribution


### Creating cond_NMDcombined df with our 2065 seed data from testing_NMD_metrics( which was for NMDminus, old condNMDminus)
cond_NMDcombined <- testing_NMD_metrics %>% 
  group_by(seed) %>%
  filter(seed == "2065") %>%
  ungroup() %>% 
  select(-seed)

View(cond_NMDcombined)


### Plotting the untransformed NMD metric
untransformed_NMD_metric_plot <- ggplot(cond_NMDcombined, aes(x = average_NW_above_base)) +
  geom_density() +
  theme_pubr() +
  labs(x = "NMD-es (unscaled)")



# TASK 10 - Transforming the NMD metric ----

## Filtering out cells that have average_NW_above_base < 1 + transforming the data
cond_NMDcombined <- cond_NMDcombined %>% 
  filter(average_NW_above_base > 1) %>% View()## a value of less than one means that the level of NMD is below baseline (i.e no NMD minus)
  mutate(average_NW_above_base = (average_NW_above_base -1))

View(cond_NMDcombined) ## from this, 423 cells identified


## Transforming the data between minimum value and a maximum value of 2
cond_NMDcombined$scaled_average_NW_above_base <- ((cond_NMDcombined$average_NW_above_base) / 
                                                    max(cond_NMDcombined$average_NW_above_base)) * 2

View(cond_NMDcombined)


## Plotting the transformed NMD metric
transformed_NMD_metric_plot <- ggplot(cond_NMDcombined, aes(x = scaled_average_NW_above_base)) +
  geom_histogram(aes(y = after_stat(density)), colour = "white", 
                     fill = "cornflowerblue", bins = 20) + ## 20 equal bins as metric is from lowest value - 2 max
  geom_density(colour = "darkred", size = 1 ) +
  theme_pubr() +
  labs(x ="NMD-es") ## GIVE THE NMD METRIC A NAME

(transformed_NMD_metric_plot | untransformed_NMD_metric_plot)

# TASK 11 - Subsetting the seurat Object for WholeCell and adding our NEW TRANSFORMED chosen NMD metric to the metadata ----

## Making a new df with only cell entries from the seurat object
seurat_cell <- subset(seurat, orig.ident=="WholeCell")


## Confirming that the correct data has been subsetted
View(seurat_cell@meta.data)


## As there are 427 cells previously identified, checking that the number matches:
tally(seurat_cell@meta.data) ### returns 427, so confirms


##Isolating the NMD metric for subsetting
NMD_metric <- cond_NMDcombined %>% 
  select(cell, scaled_average_NW_above_base, average_NW_above_base) %>%
  dplyr::rename("NMD_metric" = scaled_average_NW_above_base,
                "NMD_metric_unscaled" = average_NW_above_base)

View(NMD_metric)

saveRDS(NMD_metric, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/cond_NMDcombined.rds")


### Discussion:
### naming as cells so that left_join can be used later
### need to specify the dplyr package due to conflicts


## Appending the chosen NMD metric to the meta data

View(seurat_cell@meta.data)

seurat_cell@meta.data <- seurat_cell@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::filter(cell %in% NMD_metric$cell) %>%  ## as NMD_metric has 423 cells due to adjustment and transformation, need to filter the meta data
  left_join(NMD_metric, keep = FALSE) %>%
  column_to_rownames("cell")

View(seurat_cell@meta.data)

saveRDS(seurat_cell, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/seurat_cell.rds")



# TASK 12 - Visualsing the distribution of the (transformed) NMD metric within different cell sub-populations ----
## TASK 12a - between 'cell_class' ----

### BOX PLOT - Occurance of NMD within different mouse cell classes
comparisons1 <- 

NMD_metric_cell_class_box <- seurat_cell@meta.data %>% 
  select(cell_class, NMD_metric) %>% 
  ggplot(aes(x = cell_class,
             y = NMD_metric, fill = cell_class)) +
  geom_jitter(alpha = 0.20) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  theme_pubr(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black",
               alpha = 1,
               shape = 15, size = 3) +
  labs(x = "Cell class",y = "NMD-es") +
  scale_fill_brewer(palette = "Dark2") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("GABAergic", "Glutamatergic"))),
                                        c("Glutamatergic", as.character("NA"))), ### HOW TO GET NA TO WORK HERE?
                     label = "p.format") +
  stat_compare_means(label.y = 2.8)
stat_compare
?stat_compare_means

### RIDGE PLOTS
### ADD HERE IF WANT TO

## TASK 12b - between 'dissected_layer' ----

## BOX PLOT - Occurance of NMD within different mouse dissected layers
NMD_metric_dissected_layer_box <- seurat_cell@meta.data %>% 
  select(dissected_layer, NMD_metric) %>% 
  ggplot(aes(x = dissected_layer,
             y = NMD_metric, fill = dissected_layer)) +
  geom_jitter(alpha = 0.20) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  theme_pubr(base_size = 10) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",
               alpha = 1,
               shape = 15, size = 3) +
  stat_compare_means() +
  labs(x = "Dissected layer(s)",
       y = "NMD-es") +
  rcartocolor::scale_fill_carto_d(palette = "Prism")



## TASK 12c - between 'tasic.major.class' ----

## BOX PLOT - Occurance of NMD within different mouse dissected layers

  
NMD_metric_tasic_major_box <- seurat_cell@meta.data %>% 
  select(tasic.major.class, NMD_metric) %>% 
  ggplot(aes(x = tasic.major.class,
             y = NMD_metric, fill = tasic.major.class)) +
  geom_jitter(alpha = 0.20) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  theme_pubr(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",
               alpha = 1,
               shape = 15, size = 3) +
  labs(x = "'tasic.major.class' (cell class)",
       y = "NMD-es") +
  scale_fill_brewer(palette = "Set1") + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Excitatory", "Inhibitory"))) +
  stat_compare_means(label.y = 2.5, label.x = 0.65)


## TASK 12d - between 'tasic.sub.class' ----

## Box plot - Occurance of NMD within different mouse tasic sub classes"


### FOR LAYERS ATTRIUBUTED TO EXCITATORY
NMD_metric_tasic_sub_box_EXCITATORY <- seurat_cell@meta.data %>%
  filter(tasic.major.class == "Excitatory") %>%
  select(tasic.sub.class, NMD_metric) %>% 
  ggplot(aes(x = tasic.sub.class,
             y = NMD_metric, fill = tasic.sub.class)) +
  geom_jitter(alpha = 0.20) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  theme_pubr(base_size = 10) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",
               alpha = 1,
               shape = 15, size = 3) +
  labs(x = "'tasic.sub.class' (excitatory) ",
       y = "NMD-es") +
  rcartocolor::scale_fill_carto_d(palette = "Vivid") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("L4","L2/3"),
                                        c("L4","L5a"),
                                        c("L4","L5a2"),
                                        c("L4","L5b"),
                                        c("L4","L6a1"),
                                        c("L4","L6a2")),
                     label.y = c (2,2.12,2.24,2.36,2.48,2.60,2.72),
                     label = "p.format") +
  stat_compare_means(label.y = 2.8, label.x = 1.0)

ggsave("TASIC_SUB_CLASS_BOX_PLOT_EXC.png", plot = NMD_metric_tasic_sub_box_EXCITATORY, path = "/Users/anthonybooker/Desktop/fp01_asnmd/outputs/for_report")


### FOR LAYERS ATTRIBUTED TO INHIBITORY
NMD_metric_tasic_sub_box_INHIBITORY <- seurat_cell@meta.data %>%
  filter(tasic.major.class == "Inhibitory") %>% 
  select(tasic.sub.class, NMD_metric) %>% 
  ggplot(aes(x = tasic.sub.class,
             y = NMD_metric, fill = tasic.sub.class)) +
  geom_jitter(alpha = 0.20) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  theme_pubr(base_size = 10) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",
               alpha = 1,
               shape = 15, size = 3) +
  labs(x = "'tasic.sub.class' (inhibitory)",
       y = "NMD-es") +
  rcartocolor::scale_fill_carto_d(palette = "Bold") +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Pvalb","Sst_Cbln4"),
                                        c("Pvalb","Sst_Nr2f2"),
                                        c("Pvalb","Vip")),
                     label.y = c (2,2.15,2.3,2.45),
                     label = "p.format") +
  stat_compare_means(label.y = 2.5, label.x = 0.75)
  


# TASK 13 - Testing correlation between new transformed NMD metric and individual gene with plotting method ----

## subsetting a new df for testing
seurat_cell_testing <- seurat_cell@meta.data %>% 
  select(NMD_metric)

View(seurat_cell_testing)


## Filtering for Snrnp70

Upf1_testing <- seurat_cell@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene == "Snrnp70") %>% 
  pivot_longer(-gene, names_to = "cells", values_to = "Snrnp70_expression") %>%
  filter(cells %in% NMD_metric$cell) %>% ## as NMD_metric has 410(old)/423(new) cells
  select(-gene)

View(Upf1_testing)


## Filtering for Upf2
Upf2_testing <- seurat_cell@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>% 
  filter(gene == "Upf2") %>%
  pivot_longer(-gene, names_to = "cells", values_to = "Upf2_expression") %>% 
  filter(cells %in% NMD_metric$cell) %>%  ## as NMD_metric has 410(old)/423(new) cells
  select(-gene)


## Appending Upf1 and Upf2 into the new subsetted df
seurat_cell_testing <- seurat_cell_testing %>%
  rownames_to_column("cell") %>% 
  left_join(Upf1_testing, by = c("cell"="cells"), keep = FALSE) %>%
  left_join(Upf2_testing, by = c("cell"="cells"), keep = FALSE) %>%
  column_to_rownames("cell")

View(seurat_cell_testing)


## Plotting Upf1
Upf1_plot <- ggplot(seurat_cell_testing, 
                    aes(x = Upf1_expression,
                        y = NMD_metric)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  stat_cor(label.y = 2.25) +
  theme_pubr() +
  labs(x = "UPF1 expression",
       y = "NMD-es")


## Plotting Upf2
Upf2_plot <- ggplot(seurat_cell_testing, 
                    aes(x = NMD_metric,
                        y = Snrnp70_expression)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  stat_cor(label.y = 2.12) +
  theme_pubr() +
  labs(x = "NMD-es",
       y = "Snrnp70 expression")



# TASK 14 - Testing correlation for ALL genes ----

## Making a df so that the loop can run faster - filtered out ## filter out columns that do not match cells in meta data, i.e for 410(old)/423(new) cells!
gene_names_v2 <- as.data.frame(seurat_cell@assays$RNA@data) %>%
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>%
  filter(cell %in% NMD_metric$cell) %>%
  pivot_wider(names_from = "cell", values_from = "gene_expression")

View(gene_names_v2)


## Making a df to store the results of the loop
results_gene_names_v2 <- data.frame(gene = character(), estimate = numeric(), p.value = numeric())

View(results_gene_names_v2)


## Running a 'for' loop to calculate the Pearson's estimate and p.value for each gene
## (let the loop run in full (if interrupted run again from beginning of task 17))
for (gene_name in gene_names_v2$gene) {
  result <- gene_names_v2 %>% 
    filter(gene == gene_name) %>%
    pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>%
    column_to_rownames("cell") %>%
    mutate(select(seurat_cell@meta.data, NMD_metric)) %>%
    rownames_to_column("cell") %>%
    select(gene_expression, NMD_metric) %$%
    tidy(cor.test(gene_expression, NMD_metric)) %>%
    mutate(gene = gene_name) %>%
    select(gene, estimate, p.value)
  
  results_gene_names_v2 <- rbind(results_gene_names_v2, result)
}

View(results_gene_names_v2) 

?cor.test
asaveRDS(results_gene_names_v2, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/NEW_results_gene_names_v2.rds")

seurat_cell@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, NMD_metric) %>% View()

?cor.test

## Displaying the top 5 neg regulated genes for estimate (unfiltered)
results_gene_names_v2 %>% 
  arrange(estimate) %>% 
  head(n = 5)


## Displaying the top 5 pos regulated genes for estimate (unfiltered)
results_gene_names_v2 %>% 
  arrange(desc(estimate)) %>% 
  head(n = 25) %>% 
  View()


## Filtering out genes that have p value < 0.001
results_gene_names_v2_pvalue_filtered <- results_gene_names_v2 %>% 
  filter(p.value < 0.001)

View(results_gene_names_v2_pvalue_filtered)


## Total POSITIVE estimate with p.value less than 0.001
sum(results_gene_names_v2_pvalue_filtered$estimate > 0)


## Total NEGATIVE estimate wit p.value less than 0.001
sum(results_gene_names_v2_pvalue_filtered$estimate < 0)


## run p.adjust() on column p.value, create new column
results_gene_names_v2$p.adjust <- p.adjust(results_gene_names_v2$p.value,
                                           method = "fdr")


## Filtering out genes that have p value < 0.01
results_gene_names_v2_padjust_filtered <- results_gene_names_v2 %>% 
  filter(p.adjust < 0.01)

?p.adjust()

## Total POSITIVE and NEGATIVE estimates with p.adjusted values less than 0.001
sum(results_gene_names_v2_padjust_filtered$estimate > 0)
sum(results_gene_names_v2_padjust_filtered$estimate < 0)

View(results_gene_names_v2_padjust_filtered)



# TASK 15 - Testing the GO enrichment of NMD+/- genes ----

## For positive correlation:
positive_corr <- results_gene_names_v2_padjust_filtered %>%
  filter(estimate > 0) %>%
  pull(gene)

View(positive_corr)


## For negative correlation:
negative_corr <- results_gene_names_v2_padjust_filtered %>%
  filter(estimate < 0) %>%
  pull(gene)

View(negative_corr)


## Carrying out GO enrichment for both positive and negatively correlated genes
GO_positive_corr <- enrichGO(gene = positive_corr,
                             org.Mm.eg.db,
                             keyType = "SYMBOL",
                             universe = rownames(seurat_cell))

GO_negative_corr <- enrichGO(gene = negative_corr,
                             org.Mm.eg.db,
                             keyType = "SYMBOL",
                             universe = rownames(seurat_cell))

View(as.data.frame(GO_positive_corr))
View(as.data.frame(GO_negative_corr))

### PLOTTING

saveRDS(GO_negative_corr, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/GO_negative_corr.rds")

dotplot(GO_positive_corr)

dotplot(GO_positive_corr,
        x = "Count",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "Positively correlating genes")

dotplot(GO_negative_corr,
        x = "Count",
        color = "p.adjust",
        showCategory = 10,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "Negatively correlating genes")


# TASK 16 - Annotating the GO enrichment ----

## Filtering genes that are NMD_plus in NMD_combined
NMD_combined_NMDplus <- NMD_combined %>% 
  filter(type == "NMD_plus")

View(NMD_combined_NMDplus)


## Genes that are NMD_plus in NMD_combined (which is filtered for genes that match the seurat object) both NMD_plus and GO_positive - note this is all genes from over-representation test, not with pvalue cutoff
positive_corr_NMDannotation <- GO_positive_corr@gene %>% 
  as.data.frame() %>%
  dplyr::rename("gene_name" = ".") %>%
  filter(gene_name %in% NMD_combined$gene_name) %>%
  mutate(type = ifelse(gene_name %in% NMD_combined_NMDplus$gene_name, # must filter first
                       "NMD_plus", "NMD_minus"))

View(positive_corr_NMDannotation)


## Evaluating the proportion of NMD+/- possitively correlating genes
positive_corr_NMDannotation %>%
  group_by(type) %>%
  tally()


## Genes that are in both NMD_minus and GO_negative-, note this is all genes from over-representation test, not with pvalue cutoff
negative_corr_NMDannotation <- GO_negative_corr@gene %>% 
  as.data.frame() %>%
  dplyr::rename("gene_name" = ".") %>%
  filter(gene_name %in% NMD_combined$gene_name) %>%
  mutate(type = ifelse(gene_name %in% NMD_combined_NMDplus$gene_name, # must filter first
                       "NMD_plus", "NMD_minus"))

View(negative_corr_NMDannotation)


## Evaluating the proportion of NMD+/- negatively correlating genes
negative_corr_NMDannotation %>%
  group_by(type) %>%
  tally()

View(negative_corr_NMDannotation)


## Genes that are in both NMD_minus and GO_negative and NMD_plus and GO_positive
corr_combined_NMDannotation <- positive_corr_NMDannotation %>% 
  mutate(correlation = "positive") %>%
  bind_rows(negative_corr_NMDannotation %>% mutate(correlation = "negative"))

View(corr_combined_NMDannotation)


## Calculating the proportion of NMD+/- for both positively and negatively correlating genes
corr_combined_NMDannotation %>%
  group_by(type, correlation) %>%
  tally() %>%
  group_by(correlation) %>%
  mutate(fraction = n/sum(n))



# TASK 17 - Carrying out Fisher's tests on the GO corrs ----

### hypothesis was that neg genes should be enriched with NMD+ genes, 
### because we expect the expression of NMD+ genes to be low when NMD-es is high.
### compare this proportion to the background proportion of NMD+ genes

## For positive/negative GO corr
Fishers_GOcorr_pos_neg <- data.frame(
  "negative_corr" = c(196, 143),
  "positive_corr" = c(423, 248),
  row.names = c("NMD_minus", "NMD_plus"),
  stringsAsFactors = FALSE
)

View(Fishers_GOcorr_pos_neg)

fisher.test(Fishers_GOcorr_pos_neg)


## For positive GO corr against background
Fishers_GOcorr_pos_backgr <- data.frame(
  "positive_corr" = c(196, 143),
  "background" = c(8878, 3386), ##background is all NMD_plus_GTF genes within all genes in NMD_combined (which has been filtered for genes within the seurat object)
  row.names = c("NMD_minus", "NMD_plus"),
  stringsAsFactors = FALSE
)
lapply
View(Fishers_GOcorr_pos_backgr)

fisher.test(Fishers_GOcorr_pos_backgr)


## For negative GO corr against background
Fishers_GOcorr_neg_backgr <- data.frame(
  "negative_corr" = c(423, 248),
  "background" = c(8878, 3386), ##background is all NMD_plus_GTF genes within all genes in NMD_combined (which has been filtered for genes within the seurat object)
  row.names = c("NMD_minus", "NMD_plus"),
  stringsAsFactors = FALSE
)

View(Fishers_GOcorr_neg_backgr)

fisher.test(Fishers_GOcorr_neg_backgr)


### Discussion:
### The p.value results from positive/negative GO corr against background are highly significant







# Task 18 - Preparing a linear model for NMD: Testing assumptions ----

## TASK 18a - Evaluating that there is no Collinearity ----

### TASK 18ai - Conducting METHOD 1: with unadjusted pos_corr and neg_corr ----


### PREPARATION

## Making combined df of pos_corr and neg_corr to test correlation between
pos_neg_combined_corr <- positive_corr %>%
  as.data.frame() %>%
  bind_rows(negative_corr %>% as.data.frame()) %>%
  dplyr::rename("gene_name" = ".")

View(pos_neg_combined_corr) # genes being used


## Scaling the gene expression and making a df
seurat_cell_corr_filtered_expression <- seurat_cell@assays$RNA@data %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>% #
  filter(cell %in% NMD_metric$cell) %>% #
  pivot_wider(names_from = "cell", values_from = "gene_expression") %>% #
  filter(gene %in% pos_neg_combined_corr$gene_name) %>%
  column_to_rownames("gene") %>%
  t() %>% ## makes a matrix anyway so don't need as.matrix()
  scale()

View(seurat_cell_corr_filtered_expression)


### TESTING COLLINEARITY

## (Using corrr package method)

## Making a new df for correlation between filtered gene expressions
correlation_seurat_cell_corr_filtered_expression <- correlate(seurat_cell_corr_filtered_expression)


## Checking the correlations have been calculated
View(correlation_seurat_cell_corr_filtered_expression)


## Filtering the correlations and plotting
correlation_seurat_cell_corr_filtered_expression %>% 
  shave() %>% 
  stretch(na.rm = TRUE) %>%
  arrange(desc(r)) %>% ## add View if here you want
  ggplot(aes(r)) +
  geom_histogram(bins = 20)



### TASK 18aii - Conducting METHOD 2: with adjusted pos_corr and neg_corr for estimate/r > 0.2 and r < -0.2 ----


### PREPARATION

## Making new df for positive correlation based on new estimate/r filtering:
positive_corr2 <- results_gene_names_v2_padjust_filtered %>%
  filter(estimate > 0.20) %>%
  pull(gene)

View(positive_corr2)


## Making new df for for negative correlation based on new estimate/r filtering:
negative_corr2 <- results_gene_names_v2_padjust_filtered %>%
  filter(estimate < -0.20) %>%
  pull(gene)

View(negative_corr2)


## Making combined df of pos_corr and neg_corr to test correlation between
pos_neg_combined_corr2 <- positive_corr2 %>%
  as.data.frame() %>%
  bind_rows(negative_corr %>% as.data.frame()) %>%
  dplyr::rename("gene_name" = ".")

View(pos_neg_combined_corr2) # genes being used


## Scaling the gene expression and making a df
seurat_cell_corr_filtered_expression2 <- seurat_cell@assays$RNA@data %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>% #
  filter(cell %in% NMD_metric$cell) %>% #
  pivot_wider(names_from = "cell", values_from = "gene_expression") %>% #
  filter(gene %in% pos_neg_combined_corr2$gene_name) %>%
  column_to_rownames("gene") %>%
  t() %>% ## makes a matrix anyway so don't need as.matrix()
  scale()

View(seurat_cell_corr_filtered_expression2)


### TESTING COLLINEARITY

## (Using corrr package method)
correlation_seurat_cell_corr_filtered_expression2 <- correlate(seurat_cell_corr_filtered_expression2)

View(correlation_seurat_cell_corr_filtered_expression2)


## Filtering the correlations and plotting
correlation_seurat_cell_corr_filtered_expression2 %>% 
  shave() %>% 
  stretch(na.rm = TRUE) %>%
  arrange(desc(r)) %>% ## add View if here you want
  ggplot(aes(r)) +
  geom_histogram(bins = 20) 



### TASK 18aiii - Conducting METHOD 3: with adjusted pos_corr and neg_corr for top 50 genes from correlation (loop) testing with lowest p.adjust values ----

### PREPARATION

## Selecting the top 50 genes with the lowest p.adjust values
pos_neg_combined_corr3 <- results_gene_names_v2 %>% 
  slice_min(p.adjust, n = 50, with_ties = F) %>% ## FALSE as some genes had the same p.adjust value, will randomly select one of them
  dplyr::select(gene) %>% ## The library(org.Mm.eg.db) masks select!
  dplyr::rename("gene_name" = "gene")

View(pos_neg_combined_corr3)


## Confirming the above step was carried out correctly
results_gene_names_v2 %>% 
  arrange(p.adjust) %>% View()


## Scaling the gene expression and making a df
seurat_cell_corr_filtered_expression3 <- seurat_cell@assays$RNA@data %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>% #
  filter(cell %in% NMD_metric$cell) %>% #
  pivot_wider(names_from = "cell", values_from = "gene_expression") %>% #
  filter(gene %in% pos_neg_combined_corr3$gene_name) %>%
  column_to_rownames("gene") %>%
  t() %>% ## makes a matrix anyway so don't need as.matrix()
  scale()

?scale
View(seurat_cell_corr_filtered_expression3)


### TESTING COLLINEARITY

### (Using corrr package method)
correlation_seurat_cell_corr_filtered_expression3 <- correlate(seurat_cell_corr_filtered_expression3)

View(correlation_seurat_cell_corr_filtered_expression3)

corrr::corr
## Filtering the correlations and plotting
correlation_seurat_cell_corr_filtered_expression3 %>% 
  shave() %>% 
  stretch(na.rm = TRUE) %>% 
  arrange(desc(r)) %>% ## add View if here you want
  ggplot(aes(r)) +
  geom_histogram(bins = 20, fill = "cornflowerblue") +
  labs(x = "correlation") +
  theme_pubr()

?shave

?correlate()?

### DISCUSSION: Out of all the tests for collinearity, we can pick any as none show exact correlation
### values of 1.0. Hence, it is likely we will choose the METHOD2/expression2 genes, as it is the most stringent method
### (as the r/estimate values for selected genes have been filtered for r > 0.2 and r < -0.2)

### IN ACTUALITY METHOD 3 WAS USED FOR TRAINING THE MODEL


## TASK 18b - Testing homodescasity ----

### TASK 18bi - By using ALL the expression2/method2 collinearlity genes ----

View(pos_neg_combined_corr2$gene_name) ## list of genes being tested


## Making a df with ranked cells per gene based on expression, from 1-423 , NMD_metric, sdev and binned_exp

## as the NMD metric will be used in this step, must make sure to filter out for the 423 cells instead of
## using the 427.

test_homoscedasticity_method2_ALL <- seurat_cell@assays$RNA@data %>% 
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>%
  filter(cell %in% NMD_metric$cell) %>% 
  pivot_wider(names_from = "cell", values_from = "gene_expression") %>% 
  filter(gene %in% pos_neg_combined_corr2$gene_name) %>%
  column_to_rownames("gene") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  dplyr::mutate(NMD_metric = seurat_cell@meta.data$NMD_metric) %>%
  pivot_longer(-c(NMD_metric, cell), names_to = "gene", values_to = "expression") %>% 
  dplyr::select(cell, gene, expression, NMD_metric) %>%
  arrange(gene) %>%
  group_by(gene) %>% arrange(expression) %>%
  mutate(rank  = order(gene)) %>%
  mutate(binned_exp = ceiling(rank / 42.3)) %>% ##splitting into 10 equal bins based on 423 cells
  ungroup() %>% 
  group_by(gene, binned_exp) %>% 
  mutate(sdev = sd(NMD_metric))

View(test_homoscedasticity_method2_ALL)

## Checking if each bin per gene has 42/43 genes
test_homoscedasticity_method2_ALL %>% 
  group_by(gene, binned_exp) %>% 
  tally() %>% View()


## PLOTTING THE SDEV - CAN PUT IN REPORT (DOESNT WORK?)
test_homoscedasticity_method2_ALL %>% 
  ggplot(aes(x = binned_exp,
             y = sdev,
             colour = gene)) +
  geom_point()


## Conducting LeveneTest for all method2 sample genes
test_homoscedasticity_method2_ALL_levenetest <- test_homoscedasticity_method2_ALL %>% 
  group_by(gene) %>%
  levene_test(NMD_metric ~ as.factor(binned_exp)) %>% 
  as.data.frame()## from rstatix package

View(test_homoscedasticity_method2_ALL_levenetest)


## Identifying any genes with p values < 0.05, and removing them
test_homoscedasticity_method2_ALL_levenetest_filtered <- test_homoscedasticity_method2_ALL_levenetest %>% 
  filter(!p < 0.05)

View(test_homoscedasticity_method2_ALL_levenetest_filtered) ## now have 1003! genes (much improved over old 390!!)



### Result: the p values for the filtered variable are greater than p = 0.05, so the null hypothesis is maintained.
### As a result, it is assumed there is no difference between the variances




### TASK 18bii - By using the expression3/method3 collinearlity genes ----

## Making a df with ranked cells per gene based on expression, from 1-423 , NMD_metric, sdev and binned_exp

## as the NMD metric will be used in this step, must make sure to filter out for the 423 cells instead of
## using the 427

View(pos_neg_combined_corr3$gene_name) ## genes being tested

test_homoscedasticity_method3_ALL <- seurat_cell@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "gene_expression") %>%
  filter(cell %in% NMD_metric$cell) %>% 
  pivot_wider(names_from = "cell", values_from = "gene_expression") %>% 
  filter(gene %in% pos_neg_combined_corr3$gene_name) %>%
  column_to_rownames("gene") %>% 
  t() %>%    
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  dplyr::mutate(NMD_metric = seurat_cell@meta.data$NMD_metric) %>%
  pivot_longer(-c(NMD_metric, cell), names_to = "gene", values_to = "expression") %>% 
  dplyr::select(cell, gene, expression, NMD_metric) %>%
  arrange(gene) %>%
  group_by(gene) %>% arrange(expression) %>%
  mutate(rank  = order(gene)) %>%
  mutate(binned_exp = ceiling(rank / 42.3)) %>%
  ungroup() %>% 
  group_by(gene, binned_exp) %>% 
  mutate(sdev = sd(NMD_metric))
ceiling

View(test_homoscedasticity_method3_ALL)

## Checking if each bin per gene has 42/43
test_homoscedasticity_method3_ALL %>% 
  group_by(gene, binned_exp) %>% 
  tally() %>% View()


### PLOTTING, MAYBE INCLUDE MAYBE DONT, ASK FURSHAM?

HOMOSECDASTICITY_SDEVMETHOD3_PLOT <- test_homoscedasticity_method3_ALL %>% 
  group_by(gene, binned_exp) %>% 
  ggplot(aes(x = binned_exp,
             y = sdev,
             colour = gene)) +
  geom_point() +
  them_pubr(legend.position = "none") +
  labs(x = "binned expression") +
  theme_pubr()

ggsave("homodesc_method3_sdev_plot.png", plot = HOMOSECDASTICITY_SDEVMETHOD3_PLOT, path = "/Users/anthonybooker/Desktop/fp01_asnmd/outputs")



## Conducting LeveneTest for all method3 sample genes
test_homoscedasticity_method3_ALL_levenetest <- test_homoscedasticity_method3_ALL %>% 
  group_by(gene) %>%
  levene_test(NMD_metric ~ as.factor(binned_exp)) %>% 
  as.data.frame() ## from rstatix package

View(test_homoscedasticity_method3_ALL_levenetest)

?levene_test()

## Identifying any genes with p values < 0.05, and removing them
test_homoscedasticity_method3_ALL_levenetest_filtered <- test_homoscedasticity_method3_ALL_levenetest %>% 
  filter(!p < 0.05)

View(test_homoscedasticity_method3_ALL_levenetest_filtered) ## now have 49 genes! (improvement over old 48)

### Result: the p values for the filtered variable are greater than p = 0.05, so the null hypothesis is maintained.
### As a result, it is assumed there is no difference between the variances



# TASK 19 - Creating and training a linear model ----

## ATTEMPT 1 - WITH THE 390(old)/1003(NEW) GENES (using a random seed) ----

## Using the (test_homoscedasticity_method2_ALL_levenetest_filtered genes)
## this data includes seurat_cell_corr_filtered_expression2 genes that have been filtered to test
## for homodescasity - so 1003 genes are available for the model (passed assumptions)

## Obtaining the gene expression for the 1003 genes
seurat_cell_corr_filtered_expression2_homodesc_filter <- seurat_cell_corr_filtered_expression2 %>% 
  as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  pivot_longer(-cell, names_to = "gene", values_to = "expression") %>% 
  filter(gene %in% test_homoscedasticity_method2_ALL_levenetest_filtered$gene) %>% 
  pivot_wider(names_from = "gene", values_from = "expression") %>% 
  column_to_rownames("cell")

## Creating df for the model (with NMD_metric and the 1003 genes)
model_df <- seurat_cell_corr_filtered_expression2_homodesc_filter %>% 
  as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  filter(cell %in% (seurat_cell@meta.data %>% 
                      dplyr::select(NMD_metric) %>%
                      rownames_to_column("cell") %>% pull(cell))) %>%
  cbind(seurat_cell@meta.data %>% 
          dplyr::select(NMD_metric) %>%
          rownames_to_column("cell") %>% dplyr::select(NMD_metric)) %>%
  column_to_rownames("cell") %>% 
  relocate(NMD_metric)

View(model_df)


## Using Caret to train the model

# Setting the seed - ensuring the DataPartition code is reproducible
set.seed(950) # (out of all tested, produces the most amount of signif genes)

# Partition data frame into training and testing sets
train_model_df <- createDataPartition(model_df$NMD_metric, 
                                      p = 0.75, list = FALSE)

training_set_df <-model_df[train_model_df, ]

testing_set_df <- model_df[-train_model_df, ]

View(training_set_df)
View(testing_set_df)

saveRDS(training_set_df, "/Users/anthonybooker/Desktop/fp01_asnmd/readme/training_set_df_lm_model.rds")

# Evaluating the number of cells/rows per testing and training sets

nrow(training_set_df) ## 319 (75%) cells for training
nrow(testing_set_df) ## 104 (25%) cells for testing


# Training the linear model with the training_set
lm_model <- caret::train(
  form = NMD_metric ~ . , ##testing NMD against all predictor variables
  data = model_df,
  trControl = trainControl(method = "cv", number = 10, verboseIter = TRUE), ## verbose should track progress
  method = "lm")

summary(lm_model)
View(lm_model$results)

### Evaluating which genes are significantly contributing to the NMD_metric in the lm_model
summary(lm_model)$coefficient %>% 
  as.data.frame() %>%
  filter(`Pr(>|t|)` < 0.05) %>% View()
  rownames_to_column("gene") %>%
  filter(gene != "(Intercept)") %>% View()


# Results from the linear model training
lm_model
lm_model$results

summary(lm_model)

### Discusssion: using seed 950, 18 genes are significantly contributing to the model. However,
### the model may be suffering from the low number of data points (cells); we are only using the first 200 variables out of 1003.
### Hence, we should test the ~49 gene sample as more stringent selection process upstream



 ## ATTEMPT 2 - WITH THE 49 GENES (using optimisation of seed selection) USE THIS ONE ----

## (Using the 49 genes from test_homoscedasticity_method3_ALL_levenetest_filtered)

## Obtaining the gene expression for the 49 genes
seurat_cell_corr_filtered_expression3_homodesc_filter <- seurat_cell_corr_filtered_expression3 %>% 
  as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  pivot_longer(-cell, names_to = "gene", values_to = "expression") %>% 
  filter(gene %in% test_homoscedasticity_method3_ALL_levenetest_filtered$gene) %>% 
  pivot_wider(names_from = "gene", values_from = "expression") %>% 
  column_to_rownames("cell")

View(seurat_cell_corr_filtered_expression3_homodesc_filter)


## Creating df for the model (with NMD_metric and the 49 genes)
model_2_df <- seurat_cell_corr_filtered_expression3_homodesc_filter %>% 
  as.data.frame() %>% 
  rownames_to_column("cell") %>% 
  filter(cell %in% (seurat_cell@meta.data %>% 
                      dplyr::select(NMD_metric) %>%
                      rownames_to_column("cell") %>% pull(cell))) %>%
  cbind(seurat_cell@meta.data %>% 
          dplyr::select(NMD_metric) %>%
          rownames_to_column("cell") %>% dplyr::select(NMD_metric)) %>%
  column_to_rownames("cell") %>% 
  relocate(NMD_metric)

View(model_2_df)


## Using Caret to train the model

### OPTIMISING THE SEED SELECTION
lm_model_2_testing <- do.call(bind_rows, lapply(sample(1:10000, 1000), function(x){
  col_name <- paste0("seed",x)
  set.seed(x)
  train_model_2_df <- createDataPartition(model_2_df$NMD_metric, 
                                          p = 0.75, list = FALSE)
  training_set_2_df <-model_2_df[train_model_2_df, ]
  testing_set_2_df <- model_2_df[-train_model_2_df, ]
  lm_model_2 <- caret::train(
    form = NMD_metric ~ . , ##testing NMD against all predictor variables
    data = training_set_2_df,
    trControl = trainControl(method = "cv", number = 10, verboseIter = TRUE), ## verbose should track progress
    method = "lm")
  temp_df <- data.frame(seed = x, 
                        RMSE = pull(as.data.frame(lm_model_2$results$RMSE)),
                        multiple_r2 = pull(as.data.frame(summary(lm_model_2)$r.squared)), #we select from summary as it accounts for our variables
                        adj_r2 = pull(as.data.frame(summary(lm_model_2)$adj.r.squared)),
                        n_sig_genes = pull(summary(lm_model_2)$coefficients %>% 
                                             as.data.frame() %>%
                                             filter(`Pr(>|t|)` < 0.05) %>%
                                             rownames_to_column("gene") %>%
                                             dplyr::filter(gene != "(Intercept)") %>% 
                                             tally()))
  
}))


summary(lm_model_2)

View(lm_model_2_testing)

#### DISCUSSION:
#### From this, seed 9040 seems to give the best R2 value. 
#### We pick the adjusted R2 as it accounts for....
#### Also, a RMSE value of 0.35 indicates a relatively good fit



### TRAINING THE MODEL

## Setting the seed - ensuring the DataPartition code is reproducible
set.seed(9040) ##9040


## Partition data frame into training and testing sets
train_model_2_df <- createDataPartition(model_2_df$NMD_metric, 
                                        p = 0.75, list = FALSE)
createData
training_set_2_df <-model_2_df[train_model_2_df, ]

testing_set_2_df <- model_2_df[-train_model_2_df, ]

View(training_set_2_df)
View(testing_set_2_df)


## Evaluating the number of cells/rows per testing and training sets

nrow(training_set_2_df) ## 319 (75%) cells for training
nrow(testing_set_2_df) ## 104 (25%) cells for testing


## Training the linear model with the training_set
lm_model_2 <- caret::train(
  form = NMD_metric ~ . , ##testing NMD against all predictor variables
  data = training_set_2_df,
  trControl = trainControl(method = "cv", number = 10, verboseIter = TRUE), ## verbose should track progress
  method = "lm")

summary(lm_model_2)
lm_model_2$results


## Evaluating which genes are significantly contributing to the NMD_metric in the lm_model
summary(lm_model_2)$coefficients %>%
  as.data.frame() %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  rownames_to_column("gene") %>%
  filter(gene != "(Intercept)") %>% View()

### Discussion:
### 5 genes are identified as significantly contirbuting to the model



# TASK 20 - Using the selected model to predict and comparing the prediction to our actual results----

## Running a prediction from the model using the testing data set
predict_lm_model_2 <- as.data.frame(caret::predict.train(lm_model_2, testing_set_2_df))

View(predict_lm_model_2)
View(testing_set_2_df)


## Creating a df with observed vs predicted NMD metrics
## (Below is carried out in this manner to make sure that the cell values are aligned)
obs_exp_lm_model_2_test <- predict_lm_model_2 %>% 
  rename("predict_NMD_metric" = "caret::predict.train(lm_model_2, testing_set_2_df)") %>%
  rownames_to_column("cell") %>% 
  left_join(., testing_set_2_df %>% rownames_to_column("cell") %>% select(NMD_metric, cell)) %>%
  column_to_rownames("cell") %>% 
  rename("obs_NMD_metric" = "NMD_metric")

View(obs_exp_lm_model_2_test)


## Plotting a scatter graph for the observed vs predicted results
ggplot(obs_exp_lm_model_2_test, aes(x = obs_NMD_metric,
                                    y = predict_NMD_metric)) +
  geom_point() +
  geom_abline(colour = "red", alpha = 0.5) +
  stat_cor() +
  labs(x = "observed NMD-es",
       y = "predicted NMD-es") +
  theme_pubr()


### Discussion:
### As the R2 is quite low, there is a pretty poor fit for the model
### This may have been as a result of the low number of data points (not enough cells, 423 not sufficient)
### In future, to make a more reliable model, it is likely that a new dataset for NMD with more cells is required


## Plotting a qqplot for the observed vs predicted results
qqlot(obs_exp_lm_model_2_test$obs_NMD_metric, obs_exp_lm_model_2_test$predict_NMD_metric,
       xlab = "observed NMD-es",
       ylab = "predicted NMD-es")


### Discussion:
### qqplot shows good results
### Therefore, demonstrates that the distribution of the predicted metric is similar to the expected



# TASK 21 - Investigating linkage between cell classification and NMD predictions FOR TESTING PREDICITION ----

## Adding the cell group classifications to the predicition values
View(seurat_cell@meta.data)

obs_exp_lm_model_2_test <- obs_exp_lm_model_2_test %>% 
  rownames_to_column("cell") %>% 
  left_join(., seurat_cell@meta.data %>% 
              rownames_to_column("cell") %>%
              filter(cell %in% .$cell) %>%
              select(cell, cell_class, dissected_layer, tasic.major.class, tasic.sub.class))

View(obs_exp_lm_model_2_test)

## Investgating which sub cell groups account for slight changes in NMD
## Here, pos Range means predicted was higher. neg means predicted was lower than observed

###FOR SUB CLASS
obs_exp_lm_model_2_test %>% 
  filter(!is.na(predict_NMD_metric) | !is.na(obs_NMD_metric)) %>%
  mutate(errorsqrd = (predict_NMD_metric - obs_NMD_metric)^2) %>% ## (error squared so that all data points are positive)
  ggplot(aes(x = tasic.sub.class, ## CHANGE OUT HERE FOR TASIC.SUB.CLASS ETC !!!!!!!!!!!
             y = errorsqrd, fill = tasic.sub.class)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) + 
  geom_jitter(alpha = 0.4, colour = "black") +
  theme_pubr(legend = "none") +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black",
               alpha = 1,
               shape = 15, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  rcartocolor::scale_fill_carto_d(palette = "Bold") +
  labs(x = "'tasic.sub.class'")

###FOR CELL CLASS
obs_exp_lm_model_2_test %>% 
  filter(!is.na(predict_NMD_metric) | !is.na(obs_NMD_metric)) %>%
  mutate(errorsqrd = (predict_NMD_metric - obs_NMD_metric)^2) %>% ## (error squared so that all data points are positive)
  ggplot(aes(x = tasic.major.class, ## CHANGE OUT HERE FOR TASIC.SUB.CLASS ETC !!!!!!!!!!!
             y = errorsqrd, fill = tasic.major.class)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) + 
  geom_jitter(alpha = 0.4, colour = "black") +
  theme_pubr(legend = "none") +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black",
               alpha = 1,
               shape = 15, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  rcartocolor::scale_fill_carto_d(palette = "Bold") +
  labs(x = "'tasic.major.class' (cell class)")


## Evaluating concordence of pearson's correlation and model estimate coefficient
results_gene_names_v2 %>% 
  filter(gene %in% test_homoscedasticity_method3_ALL_levenetest_filtered$gene) %>%
  left_join(summary(lm_model_2)$coefficients %>%
          as.data.frame() %>%
          select(Estimate) %>%
          rownames_to_column("gene") %>% 
          filter(gene != "(Intercept)")) %>%
  rename("model_estimate" = Estimate) %>%
  mutate(concordance = (sign(estimate) == sign(model_estimate))) %>% ##View()
  group_by(concordance) %>% View()
  tally() %>% View()

### Discussion:
### From this, there is a higher proportion of concordance (i.e. more TRUE for same directionality than false)
### This means.....? (maybe makes our model more viable despite low accuracy due to low number of data points)


View(seurat_cell@meta.data)
# TASK 22 - Investigating linkage between cell classification and NMD predictions FOR MODEL_2_DF PREDICITION ----
predict_lm_model_2_NO_PART <- as.data.frame(caret::predict.train(lm_model_2, model_2_df))

View(predict_lm_model_2_NO_PART)
View(model_2_df)

obs_exp_lm_model_2_NO_PART <- model_2_df %>% 
  rownames_to_column("cell") %>% 
  select(cell, NMD_metric) %>%
  left_join( predict_lm_model_2_NO_PART %>% 
               rownames_to_column("cell")) %>%
  rename("predict_NMD_metric" = "caret::predict.train(lm_model_2, model_2_df)") %>% 
  left_join(., seurat_cell@meta.data %>% 
              rownames_to_column("cell") %>%
              select(cell, cell_class, dissected_layer, tasic.major.class, tasic.sub.class)) %>% 
  rename("obs_NMD_metric" = "NMD_metric")
  
View(obs_exp_lm_model_2_NO_PART)



## Plotting a scatter graph for the observed vs predicted results
ggplot(obs_exp_lm_model_2_NO_PART, aes(x = obs_NMD_metric,
                                    y = predict_NMD_metric)) +
  geom_point() +
  geom_abline(colour = "blue", alpha = 0.5) +
  labs(x = "observed NMD-es",
       y = "predicted NMD-es") +
  stat_cor() +
  theme_pubr()

## Plotting a qqplot for the observed vs predicted results

qqplot(obs_exp_lm_model_2_NO_PART$obs_NMD_metric, obs_exp_lm_model_2_NO_PART$predict_NMD_metric,
      xlab = "observed NMD-es",
      ylab = "predicted NMD-es")


## Investgating which sub cell groups account for slight changes in NMD
## Here, pos Range means predicted was higher. neg means predicted was lower than observed
obs_exp_lm_model_2_NO_PART %>% 
  filter(!is.na(predict_NMD_metric) | !is.na(obs_NMD_metric)) %>%
  mutate(errorsqrd = (predict_NMD_metric - obs_NMD_metric)^2) %>% ## (error squared so that all data points are positive)
  ggplot(aes(x = tasic.sub.class, ## CHANGE OUT HERE FOR TASIC.SUB.CLASS ETC !!!!!!!!!
             y = errorsqrd, fill = tasic.sub.class))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 2) +
  geom_jitter(alpha = 0.4, colour = "black") +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black",
               alpha = 1,
               shape = 15, size = 3) +
  labs(x = "'tasic.sub.class'") +
  theme_pubr(legend = "none") +
  rcartocolor::scale_fill_carto_d(palette = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Discussion:
### L6a1 shows the most error, likewise with the previous box plot above.
### This may suggest that this cell sub class is not ideal for the model


## Evaluating concordence of pearson's correlation and model estimate coefficient (SAME AS WITH TESTING DATA ABOVE AS SAME GENES)
results_gene_names_v2 %>% 
  filter(gene %in% test_homoscedasticity_method3_ALL_levenetest_filtered$gene) %>%
  left_join(summary(lm_model_2)$coefficients %>%
              as.data.frame() %>%
              select(Estimate) %>%
              rownames_to_column("gene") %>% 
              filter(gene != "(Intercept)")) %>%
  rename("model_estimate" = Estimate) %>%
  mutate(concordance = (sign(estimate) == sign(model_estimate))) %>% ##View()
  group_by(concordance) %>% 
  tally() %>% View()



# TASK 23 - Finalisng plots for the report

## Scaled vs unscaled NMD_metric plots + Fishers
p1 <- (transformed_NMD_metric_plot | untransformed_NMD_metric_plot)

## NMD metric in different sub-cell populations - OLD
p2 <- (( NMD_metric_cell_class_box / NMD_metric_tasic_major_box)
  |  (NMD_metric_dissected_layer_box / NMD_metric_tasic_sub_box) + theme(plot.margin = unit(c(0,0,0,40), "pt"))) +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(widths = c(1, 4)) ## for tagging with letters

## NMD-es in different cell populations
p3 <- ((NMD_metric_tasic_major_box + plot_layout(widths = c(1, 4))) + NMD_metric_tasic_sub_box_INHIBITORY / NMD_metric_tasic_sub_box_EXCITATORY) +
  plot_annotation(tag_levels = 'A') 
  
?plot_layout
ggsave("NMD-ES IN CELLS FIGURE.png", plot = p3, path = "/Users/anthonybooker/Desktop/fp01_asnmd/outputs/for_report")

## The END :D 
