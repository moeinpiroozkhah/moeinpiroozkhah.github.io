# =============================================================================
# GO & KEGG Enrichment Analysis using enrichR
# Author: Moein Piroozkhah, MD
# GitHub: https://github.com/moeinpiroozkhah
# Description: Functional enrichment analysis of a gene list using the Enrichr
#              API. Covers GO (BP/MF/CC) and KEGG pathway analysis with
#              publication-ready visualizations.
# =============================================================================

# ── 0. Install & load packages ───────────────────────────────────────────────

if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats")
if (!requireNamespace("openxlsx",quietly = TRUE)) install.packages("openxlsx")

library(enrichR)
library(ggplot2)
library(dplyr)
library(forcats)
library(openxlsx)

# ── 1. Connect to Enrichr ────────────────────────────────────────────────────
# Options: "Enrichr" (human/mouse), "FlyEnrichr", "WormEnrichr", "YeastEnrichr"

setEnrichrSite("Enrichr")

# Optional: list all available databases
# available_dbs <- listEnrichrDbs()
# head(available_dbs)

# ── 2. Define your gene list ─────────────────────────────────────────────────
# Replace with your own DEG symbols (HGNC gene symbols, e.g. from DESeq2 output)

my_gene_list <- c(
  "TP53", "BRCA1", "MYC", "VEGFA", "EGFR",
  "STAT3", "AKT1", "PTEN", "CDK2", "NOTCH1",
  "HIF1A", "BCL2", "CASPASE3", "TNF", "IL6"
  # Add your full list here ...
)

cat("Gene list loaded:", length(my_gene_list), "genes\n")

# ── 3. Define databases to query ─────────────────────────────────────────────

go_dbs <- c(
  "GO_Biological_Process_2023",
  "GO_Molecular_Function_2023",
  "GO_Cellular_Component_2023"
)

kegg_dbs <- c(
  "KEGG_2021_Human"
)

# ── 4. Run enrichment ────────────────────────────────────────────────────────

cat("Running GO enrichment...\n")
go_results <- enrichr(my_gene_list, go_dbs)

cat("Running KEGG enrichment...\n")
kegg_results <- enrichr(my_gene_list, kegg_dbs)

# ── 5. Filter significant results (adj. p < 0.05) ───────────────────────────

go_bp  <- go_results$GO_Biological_Process_2023 |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

go_mf  <- go_results$GO_Molecular_Function_2023 |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

go_cc  <- go_results$GO_Cellular_Component_2023 |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value)

kegg   <- kegg_results$KEGG_2021_Human |>
  filter(Adjusted.P.value < 0.05) |>
  arrange(Adjusted.P.value) |>
  mutate(GeneRatio = sapply(Overlap, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }))

cat(nrow(go_bp),  "significant GO-BP terms\n")
cat(nrow(go_mf),  "significant GO-MF terms\n")
cat(nrow(go_cc),  "significant GO-CC terms\n")
cat(nrow(kegg),   "significant KEGG pathways\n")

# ── 6. Visualization: GO-BP dotplot ──────────────────────────────────────────

plot_go_bp <- go_bp |>
  slice_min(Adjusted.P.value, n = 15) |>
  mutate(
    Term      = gsub(" \\(GO:.*\\)", "", Term),   # strip GO ID from label
    Term      = fct_reorder(Term, Combined.Score),
    neg_log_p = -log10(Adjusted.P.value)
  ) |>
  ggplot(aes(x = Combined.Score, y = Term,
             size = neg_log_p, color = Odds.Ratio)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#5bbf9a", high = "#1a4a3a",
                       name = "Odds Ratio") +
  scale_size_continuous(name = "-log10(adj. p)", range = c(3, 9)) +
  labs(
    title    = "GO Biological Process Enrichment",
    subtitle = paste0("Top 15 significant terms  |  n = ", length(my_gene_list), " genes"),
    x        = "Combined Score",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    axis.text.y      = element_text(size = 10)
  )

print(plot_go_bp)
ggsave("GO_BP_dotplot.pdf", plot_go_bp, width = 10, height = 7)
ggsave("GO_BP_dotplot.png", plot_go_bp, width = 10, height = 7, dpi = 300)

# ── 7. Visualization: KEGG bar chart ─────────────────────────────────────────

plot_kegg <- kegg |>
  slice_min(Adjusted.P.value, n = 20) |>
  mutate(Term = fct_reorder(Term, GeneRatio)) |>
  ggplot(aes(x = GeneRatio, y = Term,
             fill = -log10(Adjusted.P.value))) +
  geom_col(width = 0.65) +
  scale_fill_gradient(low = "#9fd4c0", high = "#1a4a3a",
                      name = "-log10(adj. p)") +
  labs(
    title    = "KEGG Pathway Enrichment",
    subtitle = paste0("Top 20 significant pathways  |  n = ", length(my_gene_list), " genes"),
    x        = "Gene Ratio (overlap / pathway size)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    axis.text.y      = element_text(size = 10)
  )

print(plot_kegg)
ggsave("KEGG_barchart.pdf", plot_kegg, width = 10, height = 8)
ggsave("KEGG_barchart.png", plot_kegg, width = 10, height = 8, dpi = 300)

# ── 8. Export all results to Excel (one sheet per database) ──────────────────

wb <- createWorkbook()

result_list <- list(
  "GO_BP"  = go_bp,
  "GO_MF"  = go_mf,
  "GO_CC"  = go_cc,
  "KEGG"   = kegg
)

for (sheet_name in names(result_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, result_list[[sheet_name]])

  # Bold header row
  addStyle(wb, sheet_name,
           style     = createStyle(textDecoration = "bold"),
           rows      = 1,
           cols      = 1:ncol(result_list[[sheet_name]]),
           gridExpand = TRUE)
}

saveWorkbook(wb, "enrichment_results.xlsx", overwrite = TRUE)
cat("✓ Results saved to enrichment_results.xlsx\n")
cat("✓ Plots saved as PDF and PNG\n")

# =============================================================================
# Output files:
#   enrichment_results.xlsx  — full results table, 4 sheets
#   GO_BP_dotplot.pdf/.png   — GO Biological Process dotplot
#   KEGG_barchart.pdf/.png   — KEGG pathway bar chart
# =============================================================================
