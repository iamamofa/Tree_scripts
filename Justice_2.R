# =============================================================================

# PHYLOGENETIC TREE + METADATA (RECTANGULAR + CIRCULAR)

# Clean, publication-ready pipeline

# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(treeio)
  library(readxl)
  library(RColorBrewer)
})

# =============================================================================

# 1. FILE PATHS

# =============================================================================

TREE_FILE <- "/home/wsp9/Petra_project/clean.core.aln.treefile"
BAP_FILE  <- "/home/wsp9/Petra_project/Bap_summary.xlsx"
META_FILE <- "/home/wsp9/Petra_project/Genome_List__Metadata.xlsx"

OUT_DIR <- "/home/wsp9/Petra_project/figures"
dir.create(OUT_DIR, showWarnings = FALSE)

# =============================================================================

# 2. LOAD DATA

# =============================================================================

tree <- read.tree(TREE_FILE)

bap_data <- read_excel(BAP_FILE) %>%
  rename(Genome_ID = `# s_id`) %>%
  select(Genome_ID, species, amr_gen, amr_res, dis_gen)

metadata <- read_excel(META_FILE)

# =============================================================================

# 3. CLEAN + MERGE METADATA

# =============================================================================

combined_meta <- metadata %>%
  left_join(bap_data, by = c("Genome ID" = "Genome_ID")) %>%
  rename(
    genome_id = `Genome ID`,
    source = Source,
    date_isolation = `Date of Isolation`
  ) %>%
  mutate(
    genome_id = as.character(genome_id),
    source = replace_na(as.character(source), "Unknown"),
    species = replace_na(as.character(species), "Unknown"),
    dis_gen = replace_na(as.character(dis_gen), "None"),
    date_isolation = as.Date(date_isolation)
  )

# =============================================================================

# 4. MATCH TREE TIPS

# =============================================================================

tree_annotation <- combined_meta %>%
  rename(label = genome_id) %>%
  filter(label %in% tree$tip.label) %>%
  slice(match(tree$tip.label, label))

# =============================================================================

# 5. BUILD METADATA TABLE

# =============================================================================

meta <- tree_annotation %>%
  mutate(
    date_group = ifelse(is.na(date_isolation),
                        "Unknown",
                        format(date_isolation, "%Y-%m")),
    disease_present = ifelse(dis_gen == "None", "No", "Yes")
  ) %>%
  select(label, source, species, date_group, disease_present)

# Order date factor properly

meta$date_group <- factor(meta$date_group,
                          levels = sort(unique(meta$date_group)))

# =============================================================================

# 6. COLOR PALETTES (CRITICAL)

# =============================================================================

source_cols <- c(
  "Animal"  = "#1b9e77",
  "Human"   = "#d95f02",
  "Unknown" = "#d9d9d9"
)

species_cols <- c(
  "Vibrio cholerae"     = "#7570b3",
  "Vibrio paracholerae" = "#e7298a",
  "Unknown"             = "#d9d9d9"
)

disease_cols <- c(
  "No"  = "#66c2a5",
  "Yes" = "#e31a1c"
)

# Date palette (gradient across time)

date_levels <- levels(meta$date_group)
date_cols <- setNames(
  colorRampPalette(brewer.pal(9, "YlGnBu"))(length(date_levels)),
  date_levels
)

# =============================================================================

# 7. RECTANGULAR TREE + HEATMAP

# =============================================================================

cat("[INFO] Building rectangular tree...\n")

p_rect <- ggtree(tree, layout = "rectangular", size = 0.6) +
  geom_tiplab(size = 2.5, align = TRUE, offset = 0.02)

heatmap_df <- meta %>%
  column_to_rownames("label")

p_rect_hm <- gheatmap(
  p_rect,
  heatmap_df,
  offset = 0.05,
  width = 1.2,
  colnames = TRUE,
  colnames_angle = 90,
  font.size = 3
)

# Apply independent color scales

p_rect_hm <- p_rect_hm +
  new_scale_fill() +
  scale_fill_manual(values = source_cols, name = "Source") +
  
  new_scale_fill() +
  scale_fill_manual(values = species_cols, name = "Species") +
  
  new_scale_fill() +
  scale_fill_manual(values = date_cols, name = "Date") +
  
  new_scale_fill() +
  scale_fill_manual(values = disease_cols, name = "Disease") +
  
  ggtitle("Phylogenetic Tree with Metadata (Rectangular)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "tree_rectangular.png"),
       p_rect_hm, width = 18, height = 12, dpi = 300)

# =============================================================================

# 8. CIRCULAR TREE WITH RINGS (BEST VERSION)

# =============================================================================

cat("[INFO] Building circular tree...\n")

tree <- ladderize(tree)

p_circ <- ggtree(tree, layout = "fan", open.angle = 15, size = 0.4)

# Compute ring positions

xmax <- max(p_circ$data$x)
ring_w <- xmax * 0.03
gap <- xmax * 0.035

x1 <- xmax + gap
x2 <- x1 + ring_w + gap
x3 <- x2 + ring_w + gap
x4 <- x3 + ring_w + gap

# ── Source ring ─────────────────────────────────────

p_circ <- p_circ + new_scale_fill() +
  geom_fruit(
    data = meta,
    geom = geom_tile,
    mapping = aes(y = label, x = x1, fill = source),
    width = ring_w
  ) +
  scale_fill_manual(values = source_cols, name = "Source")

# ── Species ring ────────────────────────────────────

p_circ <- p_circ + new_scale_fill() +
  geom_fruit(
    data = meta,
    geom = geom_tile,
    mapping = aes(y = label, x = x2, fill = species),
    width = ring_w
  ) +
  scale_fill_manual(values = species_cols, name = "Species")

# ── Date ring ───────────────────────────────────────

p_circ <- p_circ + new_scale_fill() +
  geom_fruit(
    data = meta,
    geom = geom_tile,
    mapping = aes(y = label, x = x3, fill = date_group),
    width = ring_w
  ) +
  scale_fill_manual(values = date_cols, name = "Date")

# ── Disease ring ────────────────────────────────────

p_circ <- p_circ + new_scale_fill() +
  geom_fruit(
    data = meta,
    geom = geom_tile,
    mapping = aes(y = label, x = x4, fill = disease_present),
    width = ring_w
  ) +
  scale_fill_manual(values = disease_cols, name = "Disease")

p_circ <- p_circ +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "tree_circular.png"),
       p_circ, width = 12, height = 12, dpi = 300)

# =============================================================================

# DONE

# =============================================================================

cat("\n[ALL DONE] Figures saved to:", OUT_DIR, "\n")

