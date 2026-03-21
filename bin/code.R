### Installations
#options(repos = c(CRAN = "https://cran.rstudio.com"))
#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)

## Ensure ComplexHeatmap is available (Bioconductor)
#if(!requireNamespace("ComplexHeatmap", quietly = TRUE)){
#  if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#  BiocManager::install("ComplexHeatmap", ask = FALSE)
#}

# Load core libraries in original order to avoid masking issues
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(phyloseq)
library(tibble)
library(ggplot2)
library(UpSetR)
library(reshape2)
library(ComplexUpset)
library(vegan)
library(compositions)
# added for diversity stat testing and plot annotations
library(rstatix)
#if(!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)
library(microbiome)
library(scater)
library(readr)
library(tidyverse)
library(cowplot)


### Set working directory 
setwd("D:/#2026#/data-analysis-plan/graphpad/OTUs")

# import sylph output
prokaryote.sylph <- read_tsv("prokaryote_table.tsv")   %>%
  pivot_longer(
    cols = starts_with("../"),
    names_to = "sample",
    values_to = "Abundance",
    values_drop_na = TRUE
  ) %>%
  mutate(
    sample = gsub(".*/.*/.*_([A-Za-z0-9+-]+)_.*", "\\1", sample),
    genome = stringr::str_extract(clade_name, "GCF_[0-9]+\\.[0-9]+")
  )



# OTU table
otu_table <- prokaryote.sylph %>%
  filter(!is.na(genome),
         Abundance > 0) %>%
  select(sample, genome, Abundance) %>%
  pivot_wider(
    names_from = sample, 
    values_from = Abundance, 
    values_fill = 0  # Fill missing values with 0
  )


# Taxonomy table
tax.sylph <- read_tsv("prokaryote_table.tsv") %>%
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), 
           sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(Strain = if_else(Kingdom == "NO_TAXONOMY",Phylum, Strain),
         Phylum = if_else(Kingdom == "NO_TAXONOMY",NA, Phylum),
         genome = stringr::str_extract(clade_name, "GCF_[0-9]+\\.[0-9]+")
  ) %>%
  select("genome","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>%
  filter(!is.na(genome)) %>% distinct()


if (!file.exists("metadata.txt")) {
  stop("metadata.txt not found in the working directory. Please place it in D:/#2026#/data-analysis-plan/graphpad/OTUs before running the analysis.")
}

metadata <- read_tsv("metadata.txt")

# Sanity check sample IDs between the OTU table and metadata before downstream analysis
otu_sample_ids_raw <- colnames(otu_table)[colnames(otu_table) != "genome"]
otu_sample_ids <- otu_sample_ids_raw %>%
  as.character() %>%
  stringr::str_trim() %>%
  .[!is.na(.) & . != ""] %>%
  unique()

metadata <- metadata %>%
  mutate(
    SampleID = as.character(SampleID),
    SampleID = stringr::str_trim(SampleID)
  )

otu_duplicate_ids <- tibble(SampleID = otu_sample_ids) %>%
  add_count(SampleID, name = "n_id") %>%
  filter(n_id > 1) %>%
  distinct(SampleID) %>%
  transmute(
    SampleID,
    reason = "Duplicated sample ID in OTU table"
  )

invalid_otu_ids <- otu_duplicate_ids$SampleID
otu_valid_ids <- setdiff(otu_sample_ids, invalid_otu_ids)

metadata_missing_ids <- metadata %>%
  filter(is.na(SampleID) | SampleID == "") %>%
  transmute(
    SampleID = dplyr::coalesce(SampleID, "NA"),
    reason = "Missing or blank SampleID in metadata"
  )

metadata_duplicate_ids <- metadata %>%
  filter(!is.na(SampleID), SampleID != "") %>%
  add_count(SampleID, name = "n_id") %>%
  filter(n_id > 1) %>%
  distinct(SampleID) %>%
  transmute(
    SampleID,
    reason = "Duplicated SampleID in metadata"
  )

invalid_metadata_ids <- metadata_duplicate_ids$SampleID

metadata_valid <- metadata %>%
  filter(!is.na(SampleID), SampleID != "", !SampleID %in% invalid_metadata_ids)

otu_only_ids <- setdiff(otu_valid_ids, metadata_valid$SampleID)
metadata_only_ids <- setdiff(metadata_valid$SampleID, otu_valid_ids)

excluded_ids <- bind_rows(
  otu_duplicate_ids,
  tibble(SampleID = otu_only_ids, reason = "Present in OTU table but absent from metadata"),
  tibble(SampleID = metadata_only_ids, reason = "Present in metadata but absent from OTU table"),
  metadata_missing_ids,
  metadata_duplicate_ids
) %>%
  filter(!is.na(SampleID), SampleID != "") %>%
  group_by(SampleID) %>%
  summarise(reason = paste(sort(unique(reason)), collapse = "; "), .groups = "drop") %>%
  arrange(reason, SampleID)

readr::write_tsv(excluded_ids, "excluded-IDs.txt")

kept_sample_ids <- intersect(otu_valid_ids, metadata_valid$SampleID)
removed_sample_count <- nrow(excluded_ids)
proceeding_sample_count <- length(kept_sample_ids)

if (length(kept_sample_ids) == 0) {
  stop("No shared sample IDs remain between the OTU table and metadata after exclusions. Please review excluded-IDs.txt.")
}

otu_table <- otu_table %>%
  select(genome, any_of(kept_sample_ids))

readr::write_tsv(otu_table, "filtered-otu-table.tsv")

prokaryote.sylph <- prokaryote.sylph %>%
  filter(sample %in% kept_sample_ids)

metadata <- metadata_valid %>%
  filter(SampleID %in% kept_sample_ids) %>%
  arrange(match(SampleID, kept_sample_ids))

message(
  "Sample ID sanity check complete: removed ",
  removed_sample_count,
  " samples; proceeding with ",
  proceeding_sample_count,
  " matched samples in downstream analysis. Cleaned OTU table written to filtered-otu-table.tsv."
)


############ ############ ############ #######
############ Make Phyloseq Object ############ 
############ ############ ############ #######

# taxonomy table in phyloseq format
TAX = tax.sylph %>%
  tibble::column_to_rownames("genome") %>% as.matrix() %>% tax_table()

# OTU table in phyloseq format (coerce to integer counts)
OTU_mat <- otu_table %>% column_to_rownames("genome") %>% as.matrix()
# round and ensure integer storage (estimate_richness requires integer counts)
OTU_mat <- round(OTU_mat)
storage.mode(OTU_mat) <- "integer"
OTU <- OTU_mat %>% otu_table(taxa_are_rows = T)

# Metadata table in phyloseq format
samples = metadata %>% column_to_rownames("SampleID") %>% sample_data()

carbom <- phyloseq(OTU, TAX, samples)

## Color palettes used by plots
cb_palette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
)

mycolors <- c('grey',
              '#4d9221','#7fbc41','#b8e186','#e6f5d0','#e0f3db','#c7e9b4','#c7eae5','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494',
              "#FCD9D6", "#FBB3BA", "#F99AB3", "#F778A7", "#EE579D", "#DF3797", "#A5017C", "#850178", "#49006A",
              '#bf812d','#ffeda0','#fecc5c','#fd8d3c','#e31a1c',
              "#7570b3")

shared_taxa_palette_curated <- c(
  "#1f77b4","#ff7f0e","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#393b79",
  "#5254a3","#6b6ecf","#9c9ede","#e41a1c","#377eb8",
  "#ff7f00","#984ea3","#a65628","#f781bf","#4d4d4d"
)

#### Phylum-level graphics

phylum_plot_data <- prokaryote.sylph %>% 
  filter(!sample %in% c("CAGE4481") & !grepl("t_GCA", sample)) %>%  
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species", "Strain"), 
           sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(
    Phylum = if_else(Kingdom == "NO_TAXONOMY","NO_TAXONOMY", Phylum),
    Class  = if_else(Kingdom == "NO_TAXONOMY","NO_TAXONOMY", Class),
    Order  = if_else(Kingdom == "NO_TAXONOMY","NO_TAXONOMY", Order)
  ) %>%
  filter(is.na(Class), !is.na(Phylum), is.na(Strain)) %>%
  left_join(metadata, by = c("sample"="SampleID"))

phylum_levels <- sort(unique(phylum_plot_data$Phylum))
if(length(shared_taxa_palette_curated) < length(phylum_levels)){
  phylum_palette_shared <- grDevices::colorRampPalette(shared_taxa_palette_curated)(length(phylum_levels))
} else {
  phylum_palette_shared <- shared_taxa_palette_curated[seq_along(phylum_levels)]
}
phylum_palette_shared <- stats::setNames(phylum_palette_shared, phylum_levels)

### Phylum sample-level plot removed — using only group-summary plot for clarity

## Sample-level phylum plot (compare cases and controls)
phylam <- phylum_plot_data %>%
  ggplot(aes(x = sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = phylum_palette_shared, breaks = phylum_levels, drop = FALSE) +
  facet_wrap(~Status, scales = "free_y", ncol = 2) +
  theme(axis.text.y = element_text(size = 7)) +
  labs(title = "Phylum Taxonomic Composition by TB Status (Samples)",
       x = "Sample",
       y = "Relative Abundance (%)")

ggsave("1-phylum_TB_composition_samples.tiff", plot = phylam, width = 12, height = 8, dpi = 600, device = "tiff", compression = "lzw")

## Create Status-level summary for Phylum to enable easy group comparison
phylum_summary <- phylum_plot_data %>%
  group_by(Status, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Status) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%
  ungroup()

phylum_summary_plot <- ggplot(phylum_summary, aes(x = Status, y = RelAb, fill = Phylum)) +
  geom_col(position = "stack", width = 0.6) +
  theme_classic(base_size = 12) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = phylum_palette_shared, breaks = phylum_levels, drop = FALSE) +
  labs(title = "Mean Phylum Composition by TB Status", x = "Status", y = "Relative Abundance") +
  theme(axis.text.x = element_text(face = "bold"))

# Save only the group-summary phylum plot (single clear figure)
ggsave("2-phylum_TB_composition.tiff", plot = phylum_summary_plot, width = 8, height = 6, dpi = 600, device = "tiff", compression = "lzw")

### Phylum-level statistical testing (per-sample relative abundances)
# For each Phylum, compute per-sample relative abundance and test Status difference with Wilcoxon
phylum_per_sample <- prokaryote.sylph %>%
  filter(!sample %in% c("CAGE4481") & !grepl("t_GCA", sample)) %>%
  separate(clade_name, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"), sep = "\\|", fill = "right", remove = FALSE) %>%
  mutate(Phylum = if_else(Kingdom == "NO_TAXONOMY","NO_TAXONOMY", Phylum)) %>%
  filter(is.na(Class), !is.na(Phylum), is.na(Strain)) %>%
  group_by(sample, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  left_join(metadata, by = c("sample" = "SampleID")) %>%
  filter(!is.na(Status))

# Run Wilcoxon test per Phylum
phylum_stats <- phylum_per_sample %>%
  group_by(Phylum) %>%
  wilcox_test(formula = RelAb ~ Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  ungroup()

# Compute medians per Status to infer direction
medians <- phylum_per_sample %>%
  group_by(Phylum, Status) %>%
  summarise(median_rel = median(RelAb, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = median_rel)

# Merge medians into stats and infer direction (works for two-level Status)
phylum_stats2 <- phylum_stats %>%
  dplyr::left_join(medians, by = "Phylum")

# Infer which Status column corresponds to which group and compute direction robustly
status_cols <- setdiff(colnames(medians), "Phylum")
if(length(status_cols) == 2){
  col1 <- status_cols[1]
  col2 <- status_cols[2]
  phylum_stats2 <- phylum_stats2 %>%
    dplyr::mutate(direction = dplyr::case_when(
      is.na(.data[[col1]]) | is.na(.data[[col2]]) ~ NA_character_,
      .data[[col1]] > .data[[col2]] ~ paste0("higher_in_", col1),
      TRUE ~ paste0("higher_in_", col2)
    ))
} else {
  phylum_stats2 <- phylum_stats2 %>% dplyr::mutate(direction = NA_character_)
}


# Compute per-Phylum sample counts per Status and join to stats
counts_tbl <- phylum_per_sample %>%
  group_by(Phylum, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0, names_prefix = "n_")

# Tidy columns and write to disk (use any_of to avoid select errors if some stats columns are missing)
select_cols <- c("Phylum", setdiff(colnames(counts_tbl), "Phylum"), "statistic", "p", "p.adj", "p.adj.signif", "direction")
phylum_out <- phylum_stats2 %>%
  dplyr::left_join(counts_tbl, by = "Phylum") %>%
  dplyr::select(dplyr::any_of(select_cols))

write.table(phylum_out, file = "phylum_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Simple interpretation text
sig <- phylum_out %>% filter(!is.na(p.adj) & p.adj < 0.05)
if(nrow(sig) == 0){
  interp <- "No phylum shows a significant difference between Status groups after BH correction (alpha = 0.05)."
} else {
  lines <- sig %>% dplyr::arrange(p.adj) %>% dplyr::mutate(line = paste0(Phylum, ': p.adj=', signif(p.adj,3), ', ', direction)) %>% pull(line)
  interp <- paste("Significant phyla after BH correction (alpha=0.05):", paste(lines, collapse = "; "))
}
writeLines(interp, con = "phylum_stats_report.txt")
cat(interp, "\n")

# Combine summary + sample-level phylum plots, with statistical table underneath and a single shared legend
phylum_table_df <- phylum_out %>%
  dplyr::arrange(p.adj) %>%
  dplyr::select(dplyr::any_of(c("Phylum", "p", "p.adj", "p.adj.signif", "direction")))
if(nrow(phylum_table_df) > 15) phylum_table_df <- phylum_table_df %>% head(15)

## Build table grob with white background to avoid dark rendering issues
phylum_table_gg <- gridExtra::tableGrob(
  phylum_table_df,
  rows = NULL,
  theme = gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 7), bg_params = list(fill = "white")),
    colhead = list(fg_params = list(fontsize = 8), bg_params = list(fill = "white"))
  )
)

# Convert table grob to a ggplot object so it composes predictably with cowplot
phylum_table_plot <- ggpubr::as_ggplot(phylum_table_gg)

# Keep taxa legend in the summary plot and remove legend from the sample-level plot
phylum_sum_with_legend <- phylum_summary_plot +
  theme(legend.position = "right", legend.title = element_text(size = 9))
phylum_samples_noleg <- phylam + theme(legend.position = "none")

# Build left area: summary (with legend) and sample plots side-by-side
phylum_left_main <- cowplot::plot_grid(
  phylum_sum_with_legend,
  phylum_samples_noleg,
  ncol = 2,
  rel_widths = c(1, 1.4)
)

# Right column: table only
phylum_right_col <- phylum_table_plot

# Combine into two columns: left plots and right-side table
final_phylum_plot <- cowplot::plot_grid(
  phylum_left_main,
  phylum_right_col,
  ncol = 2,
  rel_widths = c(0.70, 0.30),
  align = "h"
)

# Save the combined phylum figure
cowplot::save_plot(
  "3-phylum_comparison_with_stats.tiff",
  final_phylum_plot,
  base_width = 18,
  base_height = 12,
  dpi = 600
)


##### genus level graphics
plot_data_raw <- prokaryote.sylph %>% 
  filter(sample != "CAGE4481") %>% 
  separate(clade_name, 
           into = c("Kingdom","Phylum","Class","Order",
                    "Family","Genus","Species","Strain"), 
           sep = "\\|", fill = "right", remove = FALSE) %>%
  filter(is.na(Species), is.na(Strain), !is.na(Genus))

top15_genera_main <- plot_data_raw %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

plot_data <- plot_data_raw %>%
  mutate(Genus = if_else(Genus %in% top15_genera_main, Genus, "Other")) %>%
  group_by(sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(metadata, by = c("sample"="SampleID")) %>%
  filter(!is.na(Status))

# Palette for genera shared across the summary plot, sample plot, and legend
n_genera <- dplyr::n_distinct(plot_data$Genus)
palette_genus_curated <- c(
  "#1f77b4","#ff7f0e","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#393b79",
  "#5254a3","#6b6ecf","#9c9ede","#e41a1c","#377eb8",
  "#ff7f00","#984ea3","#a65628","#f781bf","#4d4d4d"
)
genus_levels <- sort(unique(plot_data$Genus))
if(length(palette_genus_curated) < n_genera){
  palette_genus <- grDevices::colorRampPalette(palette_genus_curated)(n_genera)
} else {
  palette_genus <- palette_genus_curated[seq_len(n_genera)]
}
palette_genus <- stats::setNames(palette_genus, genus_levels)

# Horizontal stacked bar plot
p <- ggplot(plot_data, aes(x = sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  theme_classic(base_size = 12) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = palette_genus, breaks = genus_levels, drop = FALSE) +
  facet_wrap(~Status, scales = "free_y", ncol = 2) +
  theme(axis.text.y = element_text(size = 5),
        strip.text = element_text(face = "bold")) +
  labs(title = "Top 15 Genus-level Taxonomic Composition by TB Status",
       x = "Sample",
       y = "Relative Abundance")

# Export larger sample-level genus TIFF for comparison
ggsave("4-genus_TB_composition_samples.tiff",
  plot = p,
  width = 16,
  height = 12,
  dpi = 600,
  device = "tiff",
  compression = "lzw")
# Create Status-level summary for Genus (Top 15) to aid comparison
genus_summary <- plot_data %>%
  group_by(Status, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Status) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%
  ungroup()

genus_summary_plot <- ggplot(genus_summary, aes(x = Status, y = RelAb, fill = Genus)) +
  geom_col(position = "stack", width = 0.6) +
  theme_classic(base_size = 12) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = palette_genus, breaks = genus_levels, drop = FALSE) +
  labs(title = "Mean Genus Composition by TB Status (Top 15)", x = "Status", y = "Relative Abundance") +
  theme(axis.text.x = element_text(face = "bold"))

genus_summary_plot
ggsave("5-genus_TB_composition.tiff",
       plot = genus_summary_plot,
       width = 8,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")
## Genus-level statistical testing (per-sample relative abundances)
genus_per_sample <- plot_data %>%
  group_by(sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(RelAb = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  left_join(metadata, by = c("sample" = "SampleID")) %>%
  filter(!is.na(Status))

# Run Wilcoxon test per Genus
genus_stats <- genus_per_sample %>%
  group_by(Genus) %>%
  wilcox_test(formula = RelAb ~ Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  ungroup()

# Compute medians per Status to infer direction
medians_genus <- genus_per_sample %>%
  group_by(Genus, Status) %>%
  summarise(median_rel = median(RelAb, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = median_rel)

# Merge medians into stats and infer direction robustly
genus_stats2 <- genus_stats %>% dplyr::left_join(medians_genus, by = "Genus")
status_cols_g <- setdiff(colnames(medians_genus), "Genus")
if(length(status_cols_g) == 2){
  gcol1 <- status_cols_g[1]
  gcol2 <- status_cols_g[2]
  genus_stats2 <- genus_stats2 %>%
    dplyr::mutate(direction = dplyr::case_when(
      is.na(.data[[gcol1]]) | is.na(.data[[gcol2]]) ~ NA_character_,
      .data[[gcol1]] > .data[[gcol2]] ~ paste0("higher_in_", gcol1),
      TRUE ~ paste0("higher_in_", gcol2)
    ))
} else {
  genus_stats2 <- genus_stats2 %>% dplyr::mutate(direction = NA_character_)
}

# Compute per-Genus sample counts per Status and join to stats
counts_genus_tbl <- genus_per_sample %>%
  group_by(Genus, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0, names_prefix = "n_")

select_cols_g <- c("Genus", setdiff(colnames(counts_genus_tbl), "Genus"), "statistic", "p", "p.adj", "p.adj.signif", "direction")
genus_out <- genus_stats2 %>%
  dplyr::left_join(counts_genus_tbl, by = "Genus") %>%
  dplyr::select(dplyr::any_of(select_cols_g))

write.table(genus_out, file = "genus_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Simple interpretation text for genera
sig_g <- genus_out %>% filter(!is.na(p.adj) & p.adj < 0.05)
if(nrow(sig_g) == 0){
  interp_g <- "No genus shows a significant difference between Status groups after BH correction (alpha = 0.05)."
} else {
  lines_g <- sig_g %>% dplyr::arrange(p.adj) %>% dplyr::mutate(line = paste0(Genus, ': p.adj=', signif(p.adj,3), ', ', direction)) %>% pull(line)
  interp_g <- paste("Significant genera after BH correction (alpha=0.05):", paste(lines_g, collapse = "; "))
}
writeLines(interp_g, con = "genus_stats_report.txt")
cat(interp_g, "\n")
 

## Combine summary + sample-level genus plots, with statistical table underneath and a single shared legend
table_df <- genus_out %>% dplyr::arrange(p.adj) %>% dplyr::select(dplyr::any_of(c("Genus", "p", "p.adj", "p.adj.signif", "direction")))
if(nrow(table_df) > 15) table_df <- table_df %>% head(15)

## Build table grob with white background to avoid dark rendering issues
table_gg <- gridExtra::tableGrob(
  table_df,
  rows = NULL,
  theme = gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 7), bg_params = list(fill = "white")),
    colhead = list(fg_params = list(fontsize = 8), bg_params = list(fill = "white"))
  )
)

# Convert table grob to a ggplot object so it composes predictably with cowplot
table_plot <- ggpubr::as_ggplot(table_gg)

# Keep taxa legend in the summary plot (original place) and remove shared-legend handling
gen_sum_with_legend <- genus_summary_plot + theme(legend.position = "right", legend.title = element_text(size = 9))
p_noleg <- p + theme(legend.position = "none")

# Build left area: summary (with legend) and sample plots side-by-side (reduced width)
left_main <- cowplot::plot_grid(gen_sum_with_legend, p_noleg, ncol = 2, rel_widths = c(1, 1.4))

# Right column: table only (as ggplot)
right_col <- table_plot

# Combine into two columns: left plots and right (table). tune widths so table aligns alongside plots.
final_genus_plot <- cowplot::plot_grid(
  left_main,
  right_col,
  ncol = 2,
  rel_widths = c(0.70, 0.30),
  align = "h"
)

# Save the combined figure
cowplot::save_plot("6-genus_comparison_with_stats.tiff", final_genus_plot, base_width = 18, base_height = 12, dpi = 600)

# ------------------------------------------------------------------
# Declassify all "Other" genera, plot their relative abundances,
# rerun per-genus tests on the full set, and export results.
# Plots/stats are saved as TSV/TIFF in the working directory.
# ------------------------------------------------------------------

# Determine global top 15 genera (same logic used for main plots)
top15_genera <- prokaryote.sylph %>%
  filter(!sample %in% c("CAGE4481")) %>%
  separate(clade_name,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"),
           sep = "\\|", fill = "right", remove = FALSE) %>%
  filter(is.na(Species), is.na(Strain), !is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

plot_data_other_raw <- prokaryote.sylph %>%
  filter(!sample %in% c("CAGE4481")) %>%
  separate(clade_name,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"),
           sep = "\\|", fill = "right", remove = FALSE) %>%
  filter(is.na(Species), is.na(Strain), !is.na(Genus))

# All other genera (previously lumped as "Other")
others_list <- setdiff(
  plot_data_other_raw %>% pull(Genus) %>% unique(),
  top15_genera
)

# Data for previously-lumped 'Other' genera
others_plot_data <- plot_data_other_raw %>%
  filter(Genus %in% others_list) %>%
  group_by(sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(metadata, by = c("sample" = "SampleID")) %>%
  filter(!is.na(Status))

# Shared palette for 'Other' genera plots
n_others <- dplyr::n_distinct(others_plot_data$Genus)
others_levels <- sort(unique(others_plot_data$Genus))
if(length(shared_taxa_palette_curated) < n_others){
  palette_others <- grDevices::colorRampPalette(shared_taxa_palette_curated)(n_others)
} else {
  palette_others <- shared_taxa_palette_curated[seq_len(n_others)]
}
palette_others <- stats::setNames(palette_others, others_levels)

# Summary composition of previously-lumped 'Other' genera by Status
others_summary <- others_plot_data %>%
  group_by(Status, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Status) %>%
  mutate(RelAb = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

others_summary_plot <- ggplot(others_summary, aes(x = Status, y = RelAb, fill = Genus)) +
  geom_col(position = "stack", width = 0.6) +
  theme_classic(base_size = 12) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = palette_others, breaks = others_levels, drop = FALSE) +
  labs(title = "Composition of Previously-Lumped 'Other' Genera by TB Status", x = "Status", y = "Relative Abundance") +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE, keyheight = grid::unit(0.32, "cm"), keywidth = grid::unit(0.32, "cm"))) +
  theme(axis.text.x = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = grid::unit(0.06, "cm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        legend.key.size = grid::unit(0.32, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

ggsave("7-others_genus_TB_composition.tiff",
       plot = others_summary_plot,
       width = 8,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

# Sample-level composition for 'Other' genera (stacked bars showing relative abundance)
others_sample <- others_plot_data %>%
  group_by(Status, sample) %>%
  mutate(sample_tot = sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sample = fct_reorder(sample, -sample_tot))

# Use Abundance with position = "fill" to ensure relative proportions are shown per sample
p_others_sample <- ggplot(others_sample, aes(x = sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  theme_classic(base_size = 10) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +
  scale_fill_manual(values = palette_others, breaks = others_levels, drop = FALSE) +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE, keyheight = grid::unit(0.28, "cm"), keywidth = grid::unit(0.28, "cm"))) +
  facet_wrap(~Status, scales = "free_y", ncol = 2) +
  theme(axis.text.y = element_text(size = 5), strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0.05, "cm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.key.size = grid::unit(0.28, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)) +
  labs(title = "Sample-level Composition: Previously-Lumped 'Other' Genera", x = "Sample", y = "Relative Abundance")
p_others_sample

ggsave("8-others_genus_TB_composition_samples.tiff",
       plot = p_others_sample,
       width = 16,
       height = 12,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

# ------------------------------------------------------------------
# Run per-genus Wilcoxon tests on the full (declassified) genus set,
# export full results and a version excluding Mycobacterium, and
# subset results to those genera in 'others_list'.
# ------------------------------------------------------------------

plot_data_full <- plot_data_other_raw %>%
  group_by(sample, Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(metadata, by = c("sample" = "SampleID")) %>%
  filter(!is.na(Status))

genus_per_sample_full <- plot_data_full %>%
  group_by(sample) %>%
  mutate(RelAb = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

genus_stats_full <- genus_per_sample_full %>%
  group_by(Genus) %>%
  wilcox_test(formula = RelAb ~ Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  ungroup()

medians_genus_full <- genus_per_sample_full %>%
  group_by(Genus, Status) %>%
  summarise(median_rel = median(RelAb, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = median_rel)

genus_stats_full2 <- genus_stats_full %>% dplyr::left_join(medians_genus_full, by = "Genus")
status_cols_full <- setdiff(colnames(medians_genus_full), "Genus")
if(length(status_cols_full) == 2){
  g1 <- status_cols_full[1]
  g2 <- status_cols_full[2]
  genus_stats_full2 <- genus_stats_full2 %>%
    dplyr::mutate(direction = dplyr::case_when(
      is.na(.data[[g1]]) | is.na(.data[[g2]]) ~ NA_character_,
      .data[[g1]] > .data[[g2]] ~ paste0("higher_in_", g1),
      TRUE ~ paste0("higher_in_", g2)
    ))
} else {
  genus_stats_full2 <- genus_stats_full2 %>% dplyr::mutate(direction = NA_character_)
}

counts_genus_full <- genus_per_sample_full %>%
  group_by(Genus, Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Status, values_from = n, values_fill = 0, names_prefix = "n_")

select_cols_full <- c("Genus", setdiff(colnames(counts_genus_full), "Genus"), "statistic", "p", "p.adj", "p.adj.signif", "direction")
genus_out_full <- genus_stats_full2 %>%
  dplyr::left_join(counts_genus_full, by = "Genus") %>%
  dplyr::select(dplyr::any_of(select_cols_full))

write.table(genus_out_full, file = "genus_stats_full.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Export results for the 'Other' genera specifically
genus_out_others <- genus_out_full %>% dplyr::filter(Genus %in% others_list)
write.table(genus_out_others, file = "genus_stats_others.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Also write a version excluding Mycobacterium
genus_out_no_mtb <- genus_out_full %>% dplyr::filter(!grepl("mycobact", Genus, ignore.case = TRUE))
write.table(genus_out_no_mtb, file = "genus_stats_no_mycobacterium.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Interpretation for others (excluding Mycobacterium already handled above)
sig_others <- genus_out_others %>% filter(!is.na(p.adj) & p.adj < 0.05)
if(nrow(sig_others) == 0){
  interp_others <- "No genus within the previously-lumped 'Other' group shows significant differences after BH correction (alpha = 0.05)."
} else {
  lines_o <- sig_others %>% dplyr::arrange(p.adj) %>% dplyr::mutate(line = paste0(Genus, ': p.adj=', signif(p.adj,3), ', ', direction)) %>% pull(line)
  interp_others <- paste("Significant genera within the previously-lumped 'Other' group:", paste(lines_o, collapse = "; "))
}
writeLines(interp_others, con = "genus_stats_others_report.txt")
cat(interp_others, "\n")

## Combine summary + sample-level 'Other' genus plots with statistical table and shared legend
others_table_df <- genus_out_others %>%
  dplyr::arrange(p.adj) %>%
  dplyr::select(dplyr::any_of(c("Genus", "p", "p.adj", "p.adj.signif", "direction")))
if(nrow(others_table_df) > 15) others_table_df <- others_table_df %>% head(15)

others_table_gg <- gridExtra::tableGrob(
  others_table_df,
  rows = NULL,
  theme = gridExtra::ttheme_default(
    core = list(fg_params = list(fontsize = 7), bg_params = list(fill = "white")),
    colhead = list(fg_params = list(fontsize = 8), bg_params = list(fill = "white"))
  )
)

others_table_plot <- ggpubr::as_ggplot(others_table_gg)
others_sum_with_legend <- others_summary_plot + theme(legend.position = "right", legend.title = element_text(size = 9))
others_samples_noleg <- p_others_sample + theme(legend.position = "none")

others_left_main <- cowplot::plot_grid(
  others_sum_with_legend,
  others_samples_noleg,
  ncol = 2,
  rel_widths = c(1, 1.4)
)

final_others_plot <- cowplot::plot_grid(
  others_left_main,
  others_table_plot,
  ncol = 2,
  rel_widths = c(0.70, 0.30),
  align = "h"
)

cowplot::save_plot(
  "9-others_genus_comparison_with_stats.tiff",
  final_others_plot,
  base_width = 18,
  base_height = 12,
  dpi = 600
)

############# Alpha & Beta diversity (phyloseq) #############
# Compute both genome-level and genus-level alpha diversity metrics from `carbom`
library(phyloseq)
carbom_genus <- tax_glom(carbom, taxrank = "Genus", NArm = FALSE)
carbom_genus <- prune_taxa(taxa_sums(carbom_genus) > 0, carbom_genus)

alpha_df_genome <- phyloseq::estimate_richness(carbom, measures = c("Observed", "Shannon", "Simpson")) %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  left_join(as.data.frame(sample_data(carbom)) %>% tibble::rownames_to_column(var = "SampleID"), by = "SampleID") %>%
  dplyr::rename(observed = Observed, shannon = Shannon, simpson = Simpson) %>%
  dplyr::filter(!is.na(Status) & Status %in% c("Case", "Control"))

alpha_df_genus <- phyloseq::estimate_richness(carbom_genus, measures = c("Observed", "Shannon", "Simpson")) %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  left_join(as.data.frame(sample_data(carbom_genus)) %>% tibble::rownames_to_column(var = "SampleID"), by = "SampleID") %>%
  dplyr::rename(observed = Observed, shannon = Shannon, simpson = Simpson) %>%
  dplyr::filter(!is.na(Status) & Status %in% c("Case", "Control"))

alpha_metrics <- c("shannon", "simpson", "observed")
metric_labels <- c(
  shannon = "Shannon Diversity",
  simpson = "Simpson Diversity",
  observed = "Observed Richness"
)

compute_alpha_stats <- function(alpha_df, level_label) {
  alpha_long <- alpha_df %>%
    dplyr::select(Status, dplyr::all_of(alpha_metrics)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(alpha_metrics), names_to = "metric", values_to = "value")

  alpha_stats <- alpha_long %>%
    group_by(metric) %>%
    wilcox_test(value ~ Status) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    ungroup() %>%
    mutate(level = level_label)

  sig_alpha <- alpha_stats %>%
    filter(!is.na(p.adj) & p.adj < 0.05) %>%
    arrange(p.adj)

  if(nrow(sig_alpha) == 0){
    alpha_interp <- paste0(
      level_label,
      "-level alpha diversity did not differ significantly between Case and Control samples after BH correction ",
      "(Observed richness, Shannon, and Simpson all had adjusted p-values >= 0.05). ",
      "These alpha-diversity statistics describe overall within-sample ",
      tolower(level_label),
      "-level diversity and do not identify a single genus driving the difference."
    )
  } else {
    alpha_lines <- sig_alpha %>%
      mutate(
        metric_label = unname(metric_labels[metric]),
        line = paste0(metric_label, ": p=", signif(p, 3), ", p.adj=", signif(p.adj, 3), ", significance=", p.adj.signif)
      ) %>%
      pull(line)
    alpha_interp <- paste0(
      level_label,
      "-level alpha diversity differed between Case and Control samples for the following metric(s) after BH correction: ",
      paste(alpha_lines, collapse = "; "),
      ". These results indicate differences in overall within-sample ",
      tolower(level_label),
      "-level diversity, not that one specific genus is significantly enriched by alpha-diversity testing alone."
    )
  }

  list(stats = alpha_stats, interp = alpha_interp)
}

make_alpha_plot <- function(alpha_df, metric, title_prefix) {
  ymax_metric <- tryCatch(max(alpha_df[[metric]], na.rm = TRUE), error = function(e) NA_real_)
  label_y <- if(!is.na(ymax_metric)) ymax_metric * 1.08 else NULL

  ggplot(alpha_df, aes(x = Status, y = .data[[metric]], fill = Status)) +
    geom_boxplot() +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 4, label.y = label_y) +
    theme_classic() +
    labs(
      x = "Group",
      y = unname(metric_labels[[metric]]),
      title = paste(title_prefix, "(", unname(metric_labels[[metric]]), ")", sep = "")
    ) +
    theme(legend.position = "none")
}

alpha_genome_res <- compute_alpha_stats(alpha_df_genome, "Genome")
alpha_genus_res <- compute_alpha_stats(alpha_df_genus, "Genus")

write.table(alpha_genus_res$stats, file = "10-alpha_diversity_genus_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(alpha_genome_res$stats, file = "10-supplementary-alpha_diversity_genome_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

alpha_main_interp <- paste(
  "Main analysis:",
  alpha_genus_res$interp,
  "Supplementary analysis:",
  alpha_genome_res$interp
)
writeLines(alpha_main_interp, con = "10-alpha_diversity_interpretation.txt")
cat(alpha_main_interp, "\n")

p_shannon_genus <- make_alpha_plot(alpha_df_genus, "shannon", "Genus-Level Alpha Diversity")
p_simpson_genus <- make_alpha_plot(alpha_df_genus, "simpson", "Genus-Level Alpha Diversity")
p_observed_genus <- make_alpha_plot(alpha_df_genus, "observed", "Genus-Level Alpha Diversity")

p_shannon_genome <- make_alpha_plot(alpha_df_genome, "shannon", "Genome-Level Alpha Diversity")
p_simpson_genome <- make_alpha_plot(alpha_df_genome, "simpson", "Genome-Level Alpha Diversity")
p_observed_genome <- make_alpha_plot(alpha_df_genome, "observed", "Genome-Level Alpha Diversity")

alpha_panel_main <- cowplot::plot_grid(
  p_shannon_genus,
  p_simpson_genus,
  p_observed_genus,
  p_shannon_genome,
  p_simpson_genome,
  p_observed_genome,
  ncol = 3,
  labels = c(
    "Main: Genus Shannon", "Main: Genus Simpson", "Main: Genus Observed",
    "Supplementary: Genome Shannon", "Supplementary: Genome Simpson", "Supplementary: Genome Observed"
  )
)

ggsave("10-alpha_diversity_main_and_supplementary.tiff",
       plot = alpha_panel_main,
       width = 18,
       height = 12,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

alpha_panel_genus <- cowplot::plot_grid(p_shannon_genus, p_simpson_genus, p_observed_genus, ncol = 3)
ggsave("10-alpha_diversity_genus_TB_status.tiff",
       plot = alpha_panel_genus,
       width = 18,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

alpha_panel_genome <- cowplot::plot_grid(p_shannon_genome, p_simpson_genome, p_observed_genome, ncol = 3)
ggsave("10-supplementary-alpha_diversity_genome_TB_status.tiff",
       plot = alpha_panel_genome,
       width = 18,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

# Beta diversity: main ecological view with NMDS + Bray-Curtis,
# supplementary compositional view with Aitchison distance + PCA
carbom_case_control <- subset_samples(carbom, !is.na(Status) & Status %in% c("Case", "Control"))
carbom_case_control <- prune_taxa(taxa_sums(carbom_case_control) > 0, carbom_case_control)

carbom_genus_case_control <- subset_samples(carbom_genus, !is.na(Status) & Status %in% c("Case", "Control"))
carbom_genus_case_control <- prune_taxa(taxa_sums(carbom_genus_case_control) > 0, carbom_genus_case_control)

get_sample_metadata <- function(physeq_obj) {
  meta_df <- as(sample_data(physeq_obj), "data.frame")
  data.frame(
    SampleID = rownames(meta_df),
    Status = meta_df$Status,
    row.names = rownames(meta_df),
    check.names = FALSE
  )
}

get_sample_matrix <- function(physeq_obj) {
  otu_mat <- as(otu_table(physeq_obj), "matrix")
  if (taxa_are_rows(physeq_obj)) {
    otu_mat <- t(otu_mat)
  }
  otu_mat
}

compute_bray_nmds_result <- function(physeq_obj, level_label) {
  sample_mat <- get_sample_matrix(physeq_obj)
  meta_df <- get_sample_metadata(physeq_obj)
  n_samples <- nrow(sample_mat)
  n_unique_profiles <- nrow(unique(as.data.frame(sample_mat)))

  bray_dist <- vegan::vegdist(sample_mat, method = "bray")
  permanova_res <- adonis2(bray_dist ~ meta_df$Status, permutations = 999)
  permanova_p <- tryCatch(as.numeric(permanova_res$`Pr(>F)`[1]), error = function(e) NA_real_)
  permanova_r2 <- tryCatch(as.numeric(permanova_res$R2[1]), error = function(e) NA_real_)

  nmds_warning <- NULL
  nmds_res <- withCallingHandlers(
    vegan::metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE, trace = FALSE),
    warning = function(w) {
      nmds_warning <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  near_zero_stress <- !is.null(nmds_warning) && grepl("stress is \\(nearly\\) zero", nmds_warning)

  nmds_points <- as.data.frame(nmds_res$points) %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
    left_join(meta_df, by = "SampleID")

  nmds_caption_extra <- if (near_zero_stress) {
    " Note: NMDS stress is nearly zero, which often reflects limited or very simple distance structure; interpret the ordination cautiously."
  } else {
    ""
  }

  ord_plot <- ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Status)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_classic() +
    labs(
      title = paste(level_label, "Level Beta Diversity (NMDS Bray-Curtis)"),
      x = "NMDS1",
      y = "NMDS2",
      caption = paste0(
        "n = ", n_samples,
        "; unique profiles = ", n_unique_profiles,
        ". ",
        "Ordination stress = ", signif(nmds_res$stress, 3),
        ". PERMANOVA (Bray-Curtis): p = ", ifelse(is.na(permanova_p), "NA", signif(permanova_p, 3)),
        "; R2 = ", ifelse(is.na(permanova_r2), "NA", signif(permanova_r2, 3)),
        ".",
        nmds_caption_extra
      )
    )

  interp <- if (is.na(permanova_p)) {
    paste0(level_label, "-level Bray-Curtis PERMANOVA could not be evaluated.")
  } else if (permanova_p < 0.05) {
    paste0(
      level_label,
      "-level Bray-Curtis PERMANOVA was significant (p=",
      signif(permanova_p, 3),
      ", R2=",
      signif(permanova_r2, 3),
      "), indicating an overall ecological community-composition difference between Case and Control samples. ",
      "The NMDS plot is an ordination display of that distance structure and should be interpreted visually as supportive rather than as the significance test itself.",
      if (near_zero_stress) paste0(" NMDS stress was nearly zero, suggesting the distance structure is very simple or the dataset is limited (n=", n_samples, ", unique profiles=", n_unique_profiles, "), so the ordination may appear cleaner than the underlying evidence warrants.") else ""
    )
  } else {
    paste0(
      level_label,
      "-level Bray-Curtis PERMANOVA was not significant (p=",
      signif(permanova_p, 3),
      ", R2=",
      signif(permanova_r2, 3),
      "), suggesting no strong overall ecological community-composition separation between Case and Control samples. ",
      "The NMDS plot is still useful for visualization, but apparent visual spacing should not be over-interpreted without PERMANOVA support.",
      if (near_zero_stress) paste0(" NMDS stress was nearly zero, which often happens when there are few samples or very simple distance relationships (n=", n_samples, ", unique profiles=", n_unique_profiles, "), so interpret the ordination cautiously.") else ""
    )
  }

  stats_tbl <- tibble(
    level = level_label,
    analysis = "Main",
    distance = "Bray-Curtis",
    ordination = "NMDS",
    p = permanova_p,
    r2 = permanova_r2,
    stress = nmds_res$stress,
    n_samples = n_samples,
    n_unique_profiles = n_unique_profiles,
    ordination_note = ifelse(near_zero_stress, "NMDS stress nearly zero; ordination may be overfit or based on very simple structure", NA_character_)
  )

  list(plot = ord_plot, stats = stats_tbl, interp = interp)
}

compute_aitchison_pca_result <- function(physeq_obj, level_label, pseudocount = 1) {
  sample_mat <- get_sample_matrix(physeq_obj)
  meta_df <- get_sample_metadata(physeq_obj)

  clr_mat <- log(sample_mat + pseudocount)
  clr_mat <- sweep(clr_mat, 1, rowMeans(clr_mat), FUN = "-")

  aitchison_dist <- dist(clr_mat, method = "euclidean")
  permanova_res <- adonis2(aitchison_dist ~ meta_df$Status, permutations = 999)
  permanova_p <- tryCatch(as.numeric(permanova_res$`Pr(>F)`[1]), error = function(e) NA_real_)
  permanova_r2 <- tryCatch(as.numeric(permanova_res$R2[1]), error = function(e) NA_real_)

  pca_res <- prcomp(clr_mat, center = TRUE, scale. = FALSE)
  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  pca_points <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("SampleID") %>%
    dplyr::rename(PC1 = PC1, PC2 = PC2) %>%
    left_join(meta_df, by = "SampleID")

  ord_plot <- ggplot(pca_points, aes(x = PC1, y = PC2, color = Status)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_classic() +
    labs(
      title = paste(level_label, "Level Beta Diversity (Aitchison PCA)"),
      x = paste0("PC1 (", scales::percent(var_explained[1], accuracy = 0.1), ")"),
      y = paste0("PC2 (", scales::percent(var_explained[2], accuracy = 0.1), ")"),
      caption = paste0(
        "Aitchison PERMANOVA: p = ", ifelse(is.na(permanova_p), "NA", signif(permanova_p, 3)),
        "; R2 = ", ifelse(is.na(permanova_r2), "NA", signif(permanova_r2, 3))
      )
    )

  interp <- if (is.na(permanova_p)) {
    paste0(level_label, "-level Aitchison PERMANOVA could not be evaluated.")
  } else if (permanova_p < 0.05) {
    paste0(
      level_label,
      "-level Aitchison PERMANOVA was significant (p=",
      signif(permanova_p, 3),
      ", R2=",
      signif(permanova_r2, 3),
      "), indicating a compositional difference between Case and Control samples after CLR-based analysis. ",
      "The PCA plot summarizes the dominant compositional gradients, but the PERMANOVA provides the formal group comparison."
    )
  } else {
    paste0(
      level_label,
      "-level Aitchison PERMANOVA was not significant (p=",
      signif(permanova_p, 3),
      ", R2=",
      signif(permanova_r2, 3),
      "), suggesting no strong compositional separation between Case and Control samples under the CLR/Aitchison framework. ",
      "The PCA plot should be interpreted descriptively rather than as proof of group differences."
    )
  }

  stats_tbl <- tibble(
    level = level_label,
    analysis = "Supplementary",
    distance = "Aitchison",
    ordination = "PCA",
    p = permanova_p,
    r2 = permanova_r2,
    stress = NA_real_
  )

  list(plot = ord_plot, stats = stats_tbl, interp = interp)
}

beta_main_res <- compute_bray_nmds_result(carbom_genus_case_control, "Genus")
beta_supp_res <- compute_aitchison_pca_result(carbom_genus_case_control, "Genus")

beta_stats <- bind_rows(beta_main_res$stats, beta_supp_res$stats)
write.table(beta_stats, file = "11-beta_diversity_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

beta_interp <- paste(
  "Main analysis:",
  beta_main_res$interp,
  "Supplementary analysis:",
  beta_supp_res$interp
)
writeLines(beta_interp, con = "11-beta_diversity_interpretation.txt")
cat(beta_interp, "\n")

beta_panel_main <- cowplot::plot_grid(
  beta_main_res$plot,
  beta_supp_res$plot,
  ncol = 2,
  labels = c("Main: Genus-Level NMDS Bray-Curtis", "Supplementary: Genus-Level Aitchison PCA")
)

ggsave("11-beta_diversity_main_and_supplementary.tiff",
       plot = beta_panel_main,
       width = 16,
       height = 7,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

ggsave("11-beta_diversity_genus_TB_status.tiff",
       plot = beta_main_res$plot,
       width = 8,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

ggsave("11-supplementary-beta_diversity_aitchison_genus_TB_status.tiff",
       plot = beta_supp_res$plot,
       width = 8,
       height = 6,
       dpi = 600,
       device = "tiff",
       compression = "lzw")

### Abundance of MTB species between the two groups
prokaryote.sylph %>% 
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), 
           sep = "\\|", fill = "right", remove = FALSE) %>%
  # keep species TB and filter strain to avoid double counting
  filter(Species == "s__Mycobacterium tuberculosis") %>% filter(is.na(Strain)) %>%
  left_join(metadata, by = c("sample"="SampleID")) %>%
  filter(!is.na(Status)) %>%
  ggplot(aes(x = sample, y = Abundance)) + 
  geom_bar(stat="identity", fill="#7570b3") +
  theme_classic() +  # scale_y_continuous(expand = c(0, 0))+  #scale_x_continuous(expand = c(0, 0))+
  scale_fill_manual(values=mycolors)+
  facet_wrap(~Status, scales = "free_x", ncol = 2)+ # 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) + # x axis 45deg 
  #  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + # x axis blank
  labs(title = "Taxonomic composition by TB status") +   ylab("Abundance (%)") 

# #######################################################
# ###########  ANCOMBC - DIFFERENTIAL ABUNDANCE ########
# # ANCOM-BC corrects for compositional bias in microbiome data
# ######################################################
# # #Very important for permission issues when installing
# # unlink("C:/Users/Administrator.DESKTOP-7PUVTGC/AppData/Local/R/win-library/4.5/00LOCK", recursive = TRUE)

# # install.packages("Matrix")
# # install.packages("Rcpp")
# # library(devtools)
# # install_github("biobakery/Maaslin2")
# # remove.packages(c("ANCOMBC", "CVXR"))

# # # 1. Define your library path
# # lib_path <- .libPaths()[1]

# # # 2. Force remove the stuck '00LOCK' folder
# # # The error shows it is likely named '00LOCK-CVXR' or just '00LOCK'
# # lock_dir <- file.path(lib_path, "00LOCK-CVXR")
# # unlink(lock_dir, recursive = TRUE, force = TRUE)

# # # Also check for generic lock
# # unlink(file.path(lib_path, "00LOCK"), recursive = TRUE, force = TRUE)

# # # 3. Install the specific working version (1.0-15) with NO LOCK check
# # # This bypasses the rename error
# # remotes::install_version("CVXR", 
# #                          version = "1.0-15", 
# #                          repos = "http://cran.us.r-project.org",
# #                          INSTALL_opts = '--no-lock')
# # install_github("FrederickHuangLin/ANCOMBC")
library(phyloseq)
library(ANCOMBC)
library(Maaslin2)
library(SIAMCAT)
library(ggplot2)
### Prepare AMCOMBC data - remove samples with missing data
library(ANCOMBC)
carbom <- subset_samples(carbom, !is.na(Status))
carbom <- prune_taxa(taxa_sums(carbom) > 0, carbom)

output = ancombc2(data = carbom, tax_level = "Genus",
                  fix_formula = "Status", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                  group = "Status", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.1, n_cl = 1, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20,
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

install.packages("kableExtra")
library(knitr)


output$res
output$res |>
  dplyr::select(taxon, lfc_StatusControl , q_StatusControl ) |>
  filter(q_StatusControl  < 0.05) |>
  arrange(lfc_StatusControl) |>
  head() |>
  kable()


## Extract significant taxa
res <- output$res

ancom_table <- data.frame(
  taxa = res$taxon,
  logFC = res$lfc_StatusControl,
  qvalue = res$q_StatusControl,
  significant = res$diff_StatusControl
)

ancom_sig <- subset(ancom_table, significant == TRUE)

ancom_sig

# #### interpret direction (TB enriched vs Control enriched)
# #### Because the model is Control vs Case:
# #### positive logFC → higher in Control; negative logFC → higher in TB (Case)

# # Therefore, enrich genera are: g__Dolosigranulum, g__Staphylococcus, g__Eggerthia, g__Desulfobulbus, g__Catonella, g__Leptotrichia_A (+logFC and TRUE for controls)

############ Visuals
library(ggplot2)


## Volcano plot
res <- output$res

volcano_df <- data.frame(
  taxa = res$taxon,
  logFC = res$lfc_StatusControl,
  qvalue = res$q_StatusControl
)

volcano_df$significant <- res$diff_StatusControl

ggplot(volcano_df, aes(x = logFC, y = -log10(qvalue))) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_classic() +
  labs(
    x = "Log Fold Change (Control vs Case)",
    y = "-log10(q-value)",
    title = "ANCOM-BC Differential Abundance"
  ) +
  geom_vline(xintercept = 0, linetype="dashed")

# Save ANCOM-BC results and volcano plot
try({
  write.table(output$res, file = "ancombc_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  ggsave(filename = "12-ancombc_volcano.tiff", plot = last_plot(), device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
}, silent = TRUE)

## Forest plot of significant taxa
library(dplyr)

sig_df <- data.frame(
  taxa = res$taxon,
  logFC = res$lfc_StatusControl,
  se = res$se_StatusControl,
  qvalue = res$q_StatusControl,
  significant = res$diff_StatusControl
)

sig_df <- sig_df %>%
  filter(significant == TRUE)

ggplot(sig_df, aes(x = logFC, y = reorder(taxa, logFC))) +
  geom_point(size=3, color="darkred") +
  geom_errorbarh(aes(xmin = logFC - se, xmax = logFC + se), height=0.2) +
  theme_classic() +
  labs(
    x = "Log Fold Change (Control vs Case)",
    y = "Taxa",
    title = "Significant Taxa from ANCOM-BC"
  ) +
  geom_vline(xintercept = 0, linetype="dashed")

try({
  ggsave(filename = "12-ancombc_forest.tiff", plot = last_plot(), device = "tiff", dpi = 600, width = 8, height = 6, compression = "lzw")
}, silent = TRUE)

#### Transform to genus-level
# Agglomerate to genus level
carbom_genus <- tax_glom(carbom, taxrank = "Genus")
# Remove taxa without genus annotation
carbom_genus <- subset_taxa(carbom_genus, !is.na(Genus))

#### Beta diversity + PERMANOVA
library(vegan)
# Extract OTU table
otu <- as.data.frame(otu_table(carbom_genus))
if(taxa_are_rows(carbom_genus)) otu <- t(otu)
# Bray-Curtis distance
bray_dist <- vegdist(otu, method = "bray")
# PERMANOVA
meta <- data.frame(sample_data(carbom_genus))
adonis_res <- adonis2(bray_dist ~ Status, data = meta)
adonis_res

##### Identify taxa driving differences (SIMPER)
simper_res <- simper(otu, meta$Status)
# Extract Case vs Control comparison
simper_summary <- summary(simper_res)[[1]]
# Convert to dataframe
simper_df <- as.data.frame(simper_summary)
# Add genus names
simper_df$Genus <- rownames(simper_df)
# Top drivers
top_simper <- simper_df %>%
  arrange(desc(average)) %>%
  head(20)
top_simper


#### Differential abundance (DESeq2 validation)
#BiocManager::install("DESeq2")
library(phyloseq)
library(DESeq2)
library(dplyr)
library(tibble)

# Make sure Status is a factor and set reference level
sample_data(carbom_genus)$Status <- factor(sample_data(carbom_genus)$Status)
sample_data(carbom_genus)$Status <- relevel(sample_data(carbom_genus)$Status, ref = "Control")

# Optional but recommended: remove very sparse / low-count taxa first
carbom_deseq <- prune_taxa(taxa_sums(carbom_genus) > 10, carbom_genus)

# Build DESeq2 object
dds <- phyloseq_to_deseq2(carbom_deseq, ~ Status)

# Fix zero-heavy microbiome count data
dds <- estimateSizeFactors(dds, type = "poscounts")

# Run DESeq2
dds <- DESeq(dds, sfType = "poscounts")

# Extract results
# Here, positive log2FoldChange = enriched in Case/TB
# negative log2FoldChange = enriched in Control
res <- results(dds, contrast = c("Status", "Case", "Control"))

# Convert results to data frame and keep taxon IDs
res_df <- as.data.frame(res) %>%
  rownames_to_column("TaxonID")

# Extract taxonomy from the SAME filtered phyloseq object used for DESeq2
tax_df <- as.data.frame(tax_table(carbom_deseq)) %>%
  rownames_to_column("TaxonID")

# Add Genus annotation
res_df <- res_df %>%
  left_join(tax_df %>% select(TaxonID, Genus), by = "TaxonID")

# Optional: classify direction
res_df <- res_df %>%
  mutate(
    direction = case_when(
      log2FoldChange > 0 ~ "Case_enriched",
      log2FoldChange < 0 ~ "Control_enriched",
      TRUE ~ "No_change"
    )
  )

# View results
head(res_df)

# Significant taxa
sig_res <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

sig_res

# Save results
write.table(
  res_df,
  file = "deseq2_results_with_genus.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sig_res,
  file = "deseq2_significant_results.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#######################################################
#### Covariate-adjusted differential abundance (MaAsLin2)
## MaAsLin2 allows adjustment for confounders such as age, sex, HIV status
######################################################

library(Maaslin2)

# OTU table
otu <- as.data.frame(otu_table(carbom))
if(taxa_are_rows(carbom)){
  otu <- t(otu)
}
otu <- as.data.frame(otu)

# Metadata
meta <- data.frame(sample_data(carbom), check.names = FALSE)
meta <- meta[, c("Status","Age","Sex","Study_site","HIV_status")]

meta$Status <- factor(meta$Status)
meta$Sex <- factor(meta$Sex)
meta$Study_site <- factor(meta$Study_site)
meta$HIV_status <- factor(meta$HIV_status)
meta$Age <- as.numeric(meta$Age)
meta <- meta[rownames(otu), ]

# Run MaAsLin2
fit <- Maaslin2(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin2_output",
  fixed_effects = c("Status","Age","Sex","Study_site","HIV_status"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH"
)


############################################
############################################
###########################################
## Machine Learning Biomarker Discovery (SIAMCAT)
## SIAMCAT identifies microbial signatures predicting TB.
##########################################
##########################################
library(SIAMCAT)
library(phyloseq)

# Extract OTU table
otu <- as.data.frame(otu_table(carbom))

# Ensure taxa = rows and samples = columns
if(!taxa_are_rows(carbom)){
  otu <- t(otu)
}

otu <- as.matrix(otu)

# Convert to relative abundance (per sample)
otu_rel <- sweep(otu, 2, colSums(otu), "/")

# Check each sample sums to 1
summary(colSums(otu_rel))

# Prepare metadata
meta <- data.frame(sample_data(carbom))

# Ensure metadata order matches OTU table
meta <- meta[colnames(otu_rel), ]

# Create labels
label <- create.label(
  meta = meta,
  label = "Status",
  case = "Case",
  control = "Control"
)

# Create SIAMCAT object
sc.obj <- siamcat(
  feat = otu_rel,
  meta = meta,
  label = label
)

# Filter low abundance features
sc.obj <- filter.features(
  sc.obj,
  filter.method = "abundance",
  cutoff = 0.001
)

# Normalize features
sc.obj <- normalize.features(
  sc.obj,
  norm.method = "log.std"
)

# Data splitting
sc.obj <- create.data.split(
  sc.obj,
  num.folds = 5,
  num.resample = 10
)

# Train ML model
sc.obj <- train.model(
  sc.obj,
  method = "lasso"
)

# Predictions
sc.obj <- make.predictions(sc.obj)

# Evaluate model
sc.obj <- evaluate.predictions(sc.obj)

# Plot model performance
model.evaluation.plot(sc.obj)

# identify microbial biomarkers
model.interpretation.plot(
  sc.obj,
  fn.plot = "siamcat_biomarkers.pdf"
)

# Save SIAMCAT object and a simple structure summary
try({
  saveRDS(sc.obj, file = "siamcat_object.rds")
  structure_df <- data.frame(component = names(sc.obj))
  write.table(structure_df, file = "siamcat_object_structure.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}, silent = TRUE)

# Generate aggregated report and publication-ready SIAMCAT/ANCOM figures
try({
  source("generate_report.R")
}, silent = TRUE)

