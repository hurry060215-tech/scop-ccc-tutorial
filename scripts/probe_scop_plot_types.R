options(scop_env_init = FALSE, log_message.verbose = FALSE)

suppressPackageStartupMessages({
  library(scop)
  library(ggplot2)
})

data("pancreas_sub", package = "scop")
srt <- pancreas_sub

cell_types <- sort(unique(as.character(srt$CellType)))
theta <- seq(0, 2 * pi, length.out = length(cell_types) + 1)[seq_along(cell_types)]
centers <- setNames(
  lapply(seq_along(cell_types), function(i) c(cos(theta[i]), sin(theta[i]))),
  cell_types
)
set.seed(11)
embedding <- t(vapply(as.character(srt$CellType), function(type) {
  centers[[type]] + stats::rnorm(2, sd = 0.13)
}, numeric(2)))
colnames(embedding) <- c("UMAP_1", "UMAP_2")
rownames(embedding) <- colnames(srt)
srt[["umap"]] <- SeuratObject::CreateDimReducObject(
  embeddings = embedding,
  key = "UMAP_",
  assay = SeuratObject::DefaultAssay(srt)
)

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
lr$sample <- rep(c("control", "treated", "day7"), length.out = nrow(lr))
lr$context <- lr$sample

make_bundle <- function(df, method) {
  df$method <- method
  list(
    method = method,
    long_table = df,
    pair_table = scop:::aggregate_ccc_long(df),
    liana_table = scop:::ccc_long_to_liana(df, sample_col = "sample")
  )
}

srt@tools[["LIANA"]] <- make_bundle(lr, "LIANA")
srt@tools[["CellphoneDB"]] <- make_bundle(lr, "CellphoneDB")
srt@tools[["CCC"]] <- list(
  method = "CCC",
  methods = c("LIANA", "CellphoneDB"),
  long_table = rbind(srt@tools[["LIANA"]]$long_table, srt@tools[["CellphoneDB"]]$long_table),
  pair_table = scop:::aggregate_ccc_long(rbind(srt@tools[["LIANA"]]$long_table, srt@tools[["CellphoneDB"]]$long_table)),
  liana_table = scop:::ccc_long_to_liana(rbind(srt@tools[["LIANA"]]$long_table, srt@tools[["CellphoneDB"]]$long_table), sample_col = "method"),
  metadata = list(schema = "scop_ccc_unified_v1")
)

ligand_target_df <- expand.grid(
  ligand = unique(lr$ligand)[1:8],
  target = c("Pdx1", "Mafa", "Neurog3", "Ins1", "Gcg", "Sox9"),
  stringsAsFactors = FALSE
)
ligand_target_df$weight <- round((sin(seq_len(nrow(ligand_target_df)) * 1.7) + 1) / 2, 3)
ligand_target_df$sender <- rep(c("Ductal", "Ngn3-high-EP", "Pre-endocrine"), length.out = nrow(ligand_target_df))
ligand_target_df$receiver <- "Endocrine"
niche_long <- lr
niche_long$method <- "Nichenetr"
srt@tools[["Nichenetr"]] <- list(
  method = "Nichenetr",
  long_table = niche_long,
  pair_table = scop:::aggregate_ccc_long(niche_long),
  liana_table = scop:::ccc_long_to_liana(niche_long, sample_col = "sample"),
  ligand_target_df = ligand_target_df
)

try_plot <- function(fun, type, method = "LIANA", ...) {
  message("TRY ", fun, " ", type)
  res <- tryCatch({
    if (identical(fun, "CCCHeatmap")) {
      do.call(CCCHeatmap, c(list(srt = srt, method = method, plot_type = type, title = paste("CCCHeatmap:", type), verbose = FALSE), list(...)))
    } else if (identical(fun, "CCCNetworkPlot")) {
      do.call(CCCNetworkPlot, c(list(srt = srt, method = method, plot_type = type, title = paste("CCCNetworkPlot:", type), verbose = FALSE), list(...)))
    } else {
      do.call(CCCStatPlot, c(list(srt = srt, method = method, plot_type = type, title = paste("CCCStatPlot:", type), verbose = FALSE), list(...)))
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

heatmap_args <- function(x) {
  args <- list(display_by = if (x %in% c("dot", "tile", "source_target_dot", "source_target_tile", "bubble", "bubble_lr", "pathway_bubble")) "interaction" else "aggregation", top_n = 20)
  if (identical(x, "sample_dot")) args$sample_key <- "sample"
  if (identical(x, "bubble_lr")) args$pairLR.use <- c("Jag1-Notch1", "Dll1-Notch2", "Wnt5a-Fzd2")
  if (identical(x, "ligand_target")) args$method <- "Nichenetr"
  args
}

network_args <- function(x) {
  args <- list(display_by = if (x %in% c("arrow", "sigmoid", "bipartite", "gene_chord", "lr_chord", "pathway", "individual", "individual_lr", "individual_outgoing", "individual_incoming", "embedding_network", "diffusion")) "interaction" else "aggregation", top_n = 12)
  if (identical(x, "lr_chord")) args$pairLR.use <- c("Jag1-Notch1", "Dll1-Notch2", "Wnt5a-Fzd2")
  if (identical(x, "pathway") || identical(x, "individual") || identical(x, "individual_lr")) args$signaling <- "NOTCH"
  if (identical(x, "individual_lr")) args$pairLR.use <- "Jag1-Notch1"
  if (identical(x, "embedding_network")) {
    args$group.by <- "CellType"
    args$reduction <- "umap"
  }
  args
}

stat_args <- function(x) {
  args <- list(display_by = if (x %in% c("bar", "sankey", "box", "violin")) "interaction" else "aggregation", top_n = 20)
  if (identical(x, "comparison")) {
    args$method <- "CCC"
    args$sample_col <- "method"
  }
  if (identical(x, "lr_contribution") || identical(x, "gene")) args$signaling <- "NOTCH"
  args
}

results <- do.call(rbind, c(
  lapply(heatmap_types, function(x) do.call(try_plot, c(list(fun = "CCCHeatmap", type = x), heatmap_args(x)))),
  lapply(network_types, function(x) do.call(try_plot, c(list(fun = "CCCNetworkPlot", type = x), network_args(x)))),
  lapply(stat_types, function(x) do.call(try_plot, c(list(fun = "CCCStatPlot", type = x), stat_args(x))))
))

write.csv(results, "scripts/scop_plot_type_probe.csv", row.names = FALSE)
print(results)
