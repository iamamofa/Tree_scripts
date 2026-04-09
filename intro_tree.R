# Phylogenetic Tree Annotation with Genomic Metadata
# Converted from R Markdown to R script

# =========================
# Load Required Libraries
# =========================
library(tidyverse)
library(ggtree)
library(ape)
library(treeio)
library(tidytree)
library(readxl)
library(cowplot)
library(RColorBrewer)
library(ggnewscale)
library(tidyr)
library(knitr)
library(kableExtra)

# =========================
# Import Data
# =========================

# Read phylogenetic tree
tree <- read.tree("/home/wsp9/Petra_project/clean.core.aln.treefile")

# Read BAP summary with AMR and MLST data
bap_data <- read_excel("Bap_summary.xlsx")

# Read metadata with source, date of isolation, and group information
metadata <- read_excel("Genome_List__Metadata.xlsx")

# Combine BAP and metadata
bap_data <- bap_data %>%
  rename(Genome_ID = `# s_id`) %>%
  select(Genome_ID, species, mlst, amr_gen, amr_res, amr_cls, dis_gen)

# Create comprehensive metadata file
combined_meta <- metadata %>%
  left_join(bap_data, by = c("Genome ID" = "Genome_ID")) %>%
  rename(
    genome_id = `Genome ID`,
    source = Source,
    date_isolation = `Date of Isolation`,
    group = Group
  ) %>%
  select(genome_id, source, date_isolation, group, species, mlst, amr_gen, amr_res, dis_gen) %>%
  mutate(
    date_isolation = as.Date(date_isolation),
    amr_res = replace_na(amr_res, "None"),
    amr_gen = replace_na(amr_gen, "None"),
    dis_gen = replace_na(dis_gen, "None"),
    species = replace_na(species, "Unknown")
  )

# Display the combined metadata
print(head(combined_meta, 10))

# =========================
# Prepare Data for Tree Annotation
# =========================

# Create dataframe with tree nodes
tree_data <- as_tibble(tree)

# Map metadata to tree tips
tree_annotation <- combined_meta %>%
  rename(label = genome_id) %>%
  mutate(label = as.character(label))

# Check for mismatches between tree tips and metadata
tree_tips <- tree$tip.label
meta_ids <- unique(tree_annotation$label)
missing_in_meta <- setdiff(tree_tips, meta_ids)
missing_in_tree <- setdiff(meta_ids, tree_tips)

if (length(missing_in_meta) > 0) {
  cat("Samples in tree but not in metadata:\n")
  print(missing_in_meta)
}

if (length(missing_in_tree) > 0) {
  cat("Samples in metadata but not in tree:\n")
  print(missing_in_tree)
}

# Create color palettes
n_sources <- length(unique(tree_annotation$source))
source_colors <- setNames(
  brewer.pal(max(n_sources, 3), "Set1")[1:n_sources],
  unique(tree_annotation$source)
)

n_groups <- length(unique(tree_annotation$group))
group_colors <- setNames(
  brewer.pal(max(n_groups, 3), "Set2")[1:n_groups],
  unique(tree_annotation$group)
)

n_species <- length(unique(tree_annotation$species))
species_colors <- setNames(
  brewer.pal(max(n_species, 3), "Dark2")[1:n_species],
  unique(tree_annotation$species)
)

cat("Source categories:\n")
print(names(source_colors))

cat("Group categories:\n")
print(names(group_colors))

cat("Species categories:\n")
print(names(species_colors))

# =========================
# Create Main Annotated Tree
# =========================

p <- ggtree(tree, layout = "rectangular", size = 0.6) +
  geom_tiplab(size = 3.5, hjust = -0.1) +
  xlim(0, max(tree$edge.length) * 1.4)

p <- p %<+% tree_annotation +
  aes(color = source) +
  geom_tippoint(size = 3, stroke = 0.5, aes(color = source)) +
  scale_color_manual(values = source_colors, name = "Source") +
  theme(legend.position = "right")

print(p)

# Optional save
ggsave("main_annotated_tree.png", plot = p, width = 16, height = 12, dpi = 300)

# =========================
# Add Multiple Metadata Layers
# =========================

heatmap_data <- tree_annotation %>%
  select(label, source, date_isolation, group, species, mlst, amr_res, dis_gen) %>%
  mutate(
    amr_res_short = ifelse(
      nchar(amr_res) > 50,
      paste0(substr(amr_res, 1, 47), "..."),
      amr_res
    ),
    dis_gen_short = ifelse(
      nchar(dis_gen) > 50,
      paste0(substr(dis_gen, 1, 47), "..."),
      dis_gen
    )
  )

p_main <- ggtree(tree, layout = "rectangular", size = 0.7) +
  geom_tiplab(size = 3, hjust = -0.1) +
  xlim(0, max(tree$edge.length) * 2.0)

p_full <- p_main %<+% heatmap_data +
  new_scale("fill") +
  geom_tile(aes(x = max(tree$edge.length) * 1.3, fill = source, width = 0.15), height = 0.9) +
  scale_fill_manual(values = source_colors, name = "Source", na.value = "gray90") +
  new_scale("fill") +
  geom_tile(aes(x = max(tree$edge.length) * 1.5, fill = group, width = 0.15), height = 0.9) +
  scale_fill_manual(values = group_colors, name = "Group", na.value = "gray90") +
  new_scale("fill") +
  geom_tile(aes(x = max(tree$edge.length) * 1.7, fill = species, width = 0.15), height = 0.9) +
  scale_fill_manual(values = species_colors, name = "Species", na.value = "gray90") +
  theme(legend.position = "bottom", legend.box = "vertical") +
  theme_tree2()

p_full <- p_full +
  ggtitle("Phylogenetic Tree with Integrated Metadata Annotation") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

print(p_full)

# Optional save
ggsave("faceted_metadata_tree.png", plot = p_full, width = 18, height = 14, dpi = 300)

# =========================
# Create Detailed Metadata Summary
# =========================

summary_table <- tree_annotation %>%
  group_by(source, group, species) %>%
  summarise(
    n_samples = n(),
    mlst_types = n_distinct(mlst),
    with_amr = sum(amr_res != "None"),
    with_dis = sum(dis_gen != "None"),
    date_range = paste(min(date_isolation, na.rm = TRUE), "to", max(date_isolation, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(source, group)

print(summary_table)

write.csv(summary_table, "metadata_summary_statistics.csv", row.names = FALSE)

# =========================
# AMR Profile Annotation
# =========================

tree_annotation_amr <- tree_annotation %>%
  mutate(
    amr_category = case_when(
      amr_res == "None" ~ "No AMR",
      grepl("fluoroquinolone", tolower(amr_res)) ~ "Fluoroquinolone",
      grepl("aminoglycoside", tolower(amr_res)) ~ "Aminoglycoside",
      grepl("beta.lactam|penicillin|cephalosporin", tolower(amr_res)) ~ "Beta-lactam",
      grepl("macrolide", tolower(amr_res)) ~ "Macrolide",
      TRUE ~ "Other/Multi-drug"
    )
  )

amr_colors <- c(
  "No AMR" = "#2ecc71",
  "Fluoroquinolone" = "#e74c3c",
  "Aminoglycoside" = "#f39c12",
  "Beta-lactam" = "#9b59b6",
  "Macrolide" = "#3498db",
  "Other/Multi-drug" = "#e67e22"
)

p_amr <- ggtree(tree, layout = "rectangular", size = 0.7) +
  geom_tiplab(size = 3, hjust = -0.1) +
  xlim(0, max(tree$edge.length) * 1.4)

p_amr <- p_amr %<+% tree_annotation_amr +
  geom_tippoint(size = 4, stroke = 1, aes(color = amr_category, fill = amr_category)) +
  scale_color_manual(values = amr_colors, name = "AMR Profile") +
  scale_fill_manual(values = amr_colors, name = "AMR Profile") +
  ggtitle("Phylogenetic Tree Colored by AMR Profile") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p_amr)

# Optional save
ggsave("amr_profile_tree.png", plot = p_amr, width = 16, height = 12, dpi = 300)

# =========================
# MLST and Disease Gene Annotation
# =========================

p_mlst <- ggtree(tree, layout = "rectangular", size = 0.7) +
  geom_tiplab(size = 2.5, hjust = -0.1, offset = 0.001) +
  xlim(0, max(tree$edge.length) * 2.5)

p_mlst <- p_mlst %<+% tree_annotation +
  geom_tippoint(size = 3.5, stroke = 0.5, aes(color = as.factor(mlst))) +
  scale_color_discrete(name = "MLST Type") +
  ggtitle("Phylogenetic Tree with MLST Type Annotation") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

print(p_mlst)

# Optional save
ggsave("mlst_annotation_tree.png", plot = p_mlst, width = 16, height = 12, dpi = 300)

# =========================
# Temporal Dynamics Visualization
# =========================

tree_annotation_temp <- tree_annotation %>%
  mutate(
    year_month = format(date_isolation, "%Y-%m")
  )

p_temporal <- ggtree(tree, layout = "rectangular", size = 0.7) +
  geom_tiplab(size = 3, hjust = -0.1) +
  xlim(0, max(tree$edge.length) * 1.4)

p_temporal <- p_temporal %<+% tree_annotation_temp +
  geom_tippoint(size = 4, stroke = 0.5, aes(color = date_isolation)) +
  scale_color_gradient(
    name = "Date of Isolation",
    low = "#3498db",
    high = "#e74c3c"
  ) +
  ggtitle("Phylogenetic Tree with Temporal Dynamics") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p_temporal)

# Optional save
ggsave("temporal_tree.png", plot = p_temporal, width = 16, height = 12, dpi = 300)

# =========================
# Export Full Metadata Annotation Table
# =========================

export_table <- tree_annotation %>%
  select(genome_id = label, source, date_isolation, group, species, mlst, amr_res, amr_gen, dis_gen) %>%
  arrange(genome_id) %>%
  rename(
    `Genome ID` = genome_id,
    `Source` = source,
    `Date of Isolation` = date_isolation,
    `Group` = group,
    `Species` = species,
    `MLST` = mlst,
    `AMR Resistance` = amr_res,
    `AMR Genes` = amr_gen,
    `Disease Genes` = dis_gen
  )

print(export_table)

write.csv(export_table, "complete_genome_annotation_table.csv", row.names = FALSE)

# =========================
# Session Information
# =========================

sessionInfo()

# =========================
# Notes
# =========================
cat("
This script integrates the following data sources:
- Phylogenetic Tree: clean_core_aln.treefile
- BAP Summary: species identification, MLST types, AMR resistance classes, AMR genes, disease genes
- Genome Metadata: source information, date of isolation, sample grouping

The annotation includes:
1. Sample source
2. Date of isolation
3. Sample grouping classification
4. Species identification
5. MLST sequence type
6. AMR resistance profiles
7. AMR genes detected
8. Disease/virulence genes detected
")
