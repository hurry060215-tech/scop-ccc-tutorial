options(scop_env_init = FALSE, log_message.verbose = FALSE)

suppressPackageStartupMessages({
  library(scop)
  library(ggplot2)
})

data("pancreas_sub", package = "scop")
srt <- pancreas_sub

lr <- data.frame(
  sender = c("Ductal", "Ductal", "Ngn3-low-EP", "Ngn3-low-EP", "Ngn3-high-EP", "Ngn3-high-EP", "Pre-endocrine", "Pre-endocrine", "Endocrine", "Endocrine", "Ductal", "Pre-endocrine", "Endocrine", "Ngn3-low-EP", "Ngn3-high-EP"),
  receiver = c("Ngn3-low-EP", "Ngn3-high-EP", "Ngn3-high-EP", "Pre-endocrine", "Pre-endocrine", "Endocrine", "Endocrine", "Ductal", "Ductal", "Ngn3-high-EP", "Endocrine", "Ngn3-low-EP", "Pre-endocrine", "Ductal", "Ductal"),
  ligand = c("Wnt5a", "Jag1", "Fgf10", "Dll1", "Cxcl12", "Bmp4", "Nrg1", "Tgfb1", "Ins1", "Pdgfa", "Egf", "Vegfa", "Gcg", "Sema3a", "Notch2"),
  receptor = c("Fzd2", "Notch1", "Fgfr2", "Notch2", "Cxcr4", "Bmpr1a", "Erbb3", "Tgfbr2", "Insr", "Pdgfra", "Egfr", "Kdr", "Gcgr", "Nrp1", "Jag1"),
  pathway_name = c("WNT", "NOTCH", "FGF", "NOTCH", "CXCL", "BMP", "NRG", "TGFb", "INSULIN", "PDGF", "EGF", "VEGF", "GLUCAGON", "SEMA", "NOTCH"),
  score = c(.88, .81, .76, .69, .74, .57, .92, .48, .63, .52, .59, .67, .38, .46, .44),
  pvalue = c(.001, .003, .006, .012, .004, .035, .001, .041, .009, .044, .018, .013, .061, .033, .049),
  stringsAsFactors = FALSE
)
lr$interaction_name <- paste(lr$ligand, lr$receptor, sep = "_")
lr$interaction_label <- paste(lr$ligand, lr$receptor, sep = " - ")
lr$pair_lr <- paste(lr$ligand, lr$receptor, sep = "-")
lr$classification <- lr$pathway_name
lr$method <- "LIANA"

srt@tools[["CCC"]] <- list(
  method = "CCC",
  methods = "LIANA",
  long_table = lr,
  pair_table = scop:::aggregate_ccc_long(lr),
  liana_table = scop:::ccc_long_to_liana(lr, sample_col = "method"),
  metadata = list(schema = "scop_ccc_unified_v1")
)
srt@tools[["LIANA"]] <- list(
  method = "LIANA",
  long_table = lr,
  pair_table = scop:::aggregate_ccc_long(lr),
  liana_table = scop:::ccc_long_to_liana(lr, sample_col = "method")
)

try_plot <- function(fun, type, ...) {
  message("TRY ", fun, " ", type)
  res <- tryCatch({
    if (identical(fun, "CCCHeatmap")) {
      do.call(CCCHeatmap, c(list(srt = srt, method = "LIANA", plot_type = type, title = paste("CCCHeatmap:", type), verbose = FALSE), list(...)))
    } else if (identical(fun, "CCCNetworkPlot")) {
      do.call(CCCNetworkPlot, c(list(srt = srt, method = "LIANA", plot_type = type, title = paste("CCCNetworkPlot:", type), verbose = FALSE), list(...)))
    } else {
      do.call(CCCStatPlot, c(list(srt = srt, method = "LIANA", plot_type = type, title = paste("CCCStatPlot:", type), verbose = FALSE), list(...)))
    }
  }, error = function(e) e)
  data.frame(
    fun = fun,
    plot_type = type,
    ok = !inherits(res, "error"),
    class = paste(class(res), collapse = "|"),
    message = if (inherits(res, "error")) conditionMessage(res) else "",
    stringsAsFactors = FALSE
  )
}

heatmap_types <- c("heatmap", "focused_heatmap", "dot", "matrix_dot", "tile", "source_target_dot", "source_target_tile", "sample_dot", "bubble", "bubble_lr", "pathway_bubble", "ligand_target", "role_heatmap", "role_network", "role_network_marsilea", "diff_heatmap")
network_types <- c("circle", "circle_focused", "chord", "lr_chord", "gene_chord", "pathway", "individual_lr", "individual", "individual_outgoing", "individual_incoming", "arrow", "sigmoid", "bipartite", "embedding_network", "diff_network", "diffusion")
stat_types <- c("bar", "sankey", "box", "violin", "role_scatter", "role_network", "role_network_marsilea", "pathway_summary", "comparison", "lr_contribution", "gene", "ranknet", "scatter", "role_change")

results <- do.call(rbind, c(
  lapply(heatmap_types, function(x) try_plot("CCCHeatmap", x, display_by = if (x %in% c("dot", "tile", "source_target_dot", "source_target_tile", "bubble", "bubble_lr")) "interaction" else "aggregation", top_n = 20)),
  lapply(network_types, function(x) try_plot("CCCNetworkPlot", x, display_by = if (x %in% c("arrow", "sigmoid")) "interaction" else "aggregation", top_n = 12, group.by = "CellType", reduction = "umap", signaling = "NOTCH", pairLR.use = "Jag1-Notch1")),
  lapply(stat_types, function(x) try_plot("CCCStatPlot", x, display_by = if (x %in% c("bar", "sankey", "box", "violin")) "interaction" else "aggregation", top_n = 20, signaling = "NOTCH"))
))

write.csv(results, "scripts/scop_plot_type_probe.csv", row.names = FALSE)
print(results)
