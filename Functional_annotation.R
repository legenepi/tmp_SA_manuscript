#!/usr/bin/env Rscript


#Rationale: Functional annotation of the 517 variants in credible sets and Supplementary Figure for the manuscript.

library(tidyverse)
library(data.table)
library(dplyr)
library(ggsci) # for scientific publication
library(ggrepel) # for better label placement

annot <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/FAVOR_credset_chrpos38_2023_08_08.csv", sep=",", stringsAsFactors = FALSE, fill=TRUE, header=TRUE)
annot <- annot %>% select("Variant..VCF.","Chromosome","Position","Genecode.Comprehensive.Category") %>%
    rename(id="Variant..VCF.",chromosome="Chromosome",posb38="Position")
annot <- as.data.frame(sapply(annot, function(x) gsub("\"", "", x)))
annot <- annot %>% rename(Functional_annotation=Genecode.Comprehensive.Category)
#remove MHC region credset:
annot <- annot %>% filter(!str_starts(id, '6-3'))

#April 2025:
favor_file <- "/alice-home/3/n/nnp5/PhD/PhD_project/Var_to_Gene/input/Additional_credset_snps_March2025/20250310_FAVOR_output_additional_credset_SNPs_processed.csv.gz"
annot <- fread(favor_file,na.strings = c("",NA))
annot <- annot %>% select("VariantVcf","Chromosome","Position","GenecodeComprehensiveCategory") %>%
    rename(id="VariantVcf",chromosome="Chromosome",posb38="Position")
annot <- as.data.frame(sapply(annot, function(x) gsub("\"", "", x)))
annot <- annot %>% rename(Functional_annotation=GenecodeComprehensiveCategory)
table(annot$Func)

##pie chart:
# Add a count column based on each Functional_annotation category
annot_with_counts <- annot %>%
  group_by(Functional_annotation) %>%
  mutate(count = n()) %>%
  ungroup()

annot_with_counts <- annot_with_counts %>% select(Functional_annotation, count) %>% unique()
annot_with_counts$percentage <- round(annot_with_counts$count / sum(annot_with_counts$count) * 100, 2)
#rename upstream;downstream with upstream as the variant is upstream to the gene (manual check in dbSNP):
annot_with_counts <- annot_with_counts %>%
  mutate(Functional_annotation = ifelse(Functional_annotation == "upstream;downstream", "upstream", Functional_annotation))


# Create pie chart with counts and percentages as labels
# Get the positions
df2 <- annot_with_counts %>%
  mutate(csum = rev(cumsum(rev(percentage))),
         pos = percentage/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percentage/2, pos))

p <- ggplot(annot_with_counts, aes(x = "" , y = percentage, fill = fct_inorder(Functional_annotation))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(percentage, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE,
                   fill = NA) +
  theme_void() +
  # Apply a scientific-friendly color palette
  scale_fill_npg(name = "Functional Annotation")

ggsave("pie_chart_functional_annotation.png",
       plot = p,
       width = 8, height = 6,
       dpi = 300)  # 300 DPI is generally suitable for print quality

# Alternatively, save as a PDF
ggsave("pie_chart_functional_annotation.pdf",
       plot = p,
       width = 8, height = 6)  # PDF does not require DPI settings