options(scop_env_init = FALSE, log_message.verbose = FALSE)

suppressPackageStartupMessages({
  library(scop)
  library(ggplot2)
})

out_dir <- file.path("assets", "scop-gallery")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  liana_table = scop:::ccc_long_to_liana(rbind(srt@tools[["LIANA"]]$long_table, srt@tools[["CellphoneDB"]]$long_table), sample_col = "method")
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
srt@tools[["Nichenetr"]] <- c(make_bundle(niche_long, "Nichenetr"), list(ligand_target_df = ligand_target_df))

save_any <- function(name, obj, width = 7.2, height = 5) {
  path <- file.path(out_dir, paste0(name, ".png"))
  if (is.list(obj) && !inherits(obj, "ggplot")) {
    rec_idx <- which(vapply(obj, inherits, logical(1), what = "recordedplot"))
    if (length(rec_idx) > 0L) {
      obj <- obj[[rec_idx[1]]]
    }
  }
  if (inherits(obj, "recordedplot")) {
    png(path, width = width, height = height, units = "in", res = 180, bg = "white")
    replayPlot(obj)
    dev.off()
  } else {
    ggplot2::ggsave(path, obj, width = width, height = height, dpi = 180, bg = "white")
  }
  path
}

run_plot <- function(name, expr, width = 7.2, height = 5) {
  message("Generating ", name)
  out <- tryCatch(force(expr), error = function(e) e)
  if (inherits(out, "error")) {
    message("FAILED ", name, ": ", conditionMessage(out))
    return(data.frame(name = name, ok = FALSE, message = conditionMessage(out), stringsAsFactors = FALSE))
  }
  save_any(name, out, width = width, height = height)
  data.frame(name = name, ok = TRUE, message = "", stringsAsFactors = FALSE)
}

results <- list()
add <- function(name, expr, width = 7.2, height = 5) {
  results[[length(results) + 1L]] <<- run_plot(name, expr, width, height)
}

add("heatmap", CCCHeatmap(srt, method = "LIANA", plot_type = "heatmap", display_by = "aggregation", title = "CCCHeatmap: heatmap", verbose = FALSE), 7, 5)
add("focused_heatmap", CCCHeatmap(srt, method = "LIANA", plot_type = "focused_heatmap", display_by = "aggregation", title = "CCCHeatmap: focused_heatmap", verbose = FALSE), 7, 5)
add("dot", CCCHeatmap(srt, method = "LIANA", plot_type = "dot", display_by = "interaction", top_n = 15, title = "CCCHeatmap: dot", verbose = FALSE), 8.5, 5.2)
add("tile", CCCHeatmap(srt, method = "LIANA", plot_type = "tile", display_by = "interaction", top_n = 15, title = "CCCHeatmap: tile", verbose = FALSE), 8.5, 5.2)
add("source_target_dot", CCCHeatmap(srt, method = "LIANA", plot_type = "source_target_dot", display_by = "interaction", top_n = 15, title = "CCCHeatmap: source_target_dot", verbose = FALSE), 8.5, 5.2)
add("source_target_tile", CCCHeatmap(srt, method = "LIANA", plot_type = "source_target_tile", display_by = "interaction", top_n = 15, title = "CCCHeatmap: source_target_tile", verbose = FALSE), 8.5, 5.2)
add("sample_dot", CCCHeatmap(srt, method = "LIANA", plot_type = "sample_dot", display_by = "interaction", top_n = 15, sample_key = "sample", title = "CCCHeatmap: sample_dot", verbose = FALSE), 8.5, 5.2)
add("bubble", CCCHeatmap(srt, method = "LIANA", plot_type = "bubble", display_by = "interaction", top_n = 12, title = "CCCHeatmap: bubble", verbose = FALSE), 8.5, 5.2)
add("bubble_lr", CCCHeatmap(srt, method = "LIANA", plot_type = "bubble_lr", display_by = "interaction", pairLR.use = c("Jag1-Notch1", "Dll1-Notch2", "Wnt5a-Fzd2"), title = "CCCHeatmap: bubble_lr", verbose = FALSE), 8.5, 5.2)
add("pathway_bubble", CCCHeatmap(srt, method = "LIANA", plot_type = "pathway_bubble", display_by = "interaction", top_n = 12, title = "CCCHeatmap: pathway_bubble", verbose = FALSE), 8.5, 5.2)
add("ligand_target", CCCHeatmap(srt, method = "Nichenetr", plot_type = "ligand_target", top_n = 24, title = "CCCHeatmap: ligand_target", verbose = FALSE), 8, 5)

add("circle", CCCNetworkPlot(srt, method = "LIANA", plot_type = "circle", display_by = "aggregation", title = "CCCNetworkPlot: circle", verbose = FALSE), 7, 5.5)
add("circle_focused", CCCNetworkPlot(srt, method = "LIANA", plot_type = "circle_focused", display_by = "aggregation", min_interaction_threshold = .55, title = "CCCNetworkPlot: circle_focused", verbose = FALSE), 7, 5.5)
add("chord", CCCNetworkPlot(srt, method = "LIANA", plot_type = "chord", display_by = "aggregation", top_n = 12, title = "CCCNetworkPlot: chord", verbose = FALSE), 7, 5.5)
add("gene_chord", CCCNetworkPlot(srt, method = "LIANA", plot_type = "gene_chord", display_by = "interaction", top_n = 10, title = "CCCNetworkPlot: gene_chord", verbose = FALSE), 7, 5.5)
add("pathway", CCCNetworkPlot(srt, method = "LIANA", plot_type = "pathway", display_by = "interaction", signaling = "NOTCH", title = "CCCNetworkPlot: pathway", verbose = FALSE), 7, 5.5)
add("individual_lr", CCCNetworkPlot(srt, method = "LIANA", plot_type = "individual_lr", display_by = "interaction", signaling = "NOTCH", pairLR.use = "Jag1-Notch1", title = "CCCNetworkPlot: individual_lr", verbose = FALSE), 7, 5.5)
add("individual", CCCNetworkPlot(srt, method = "LIANA", plot_type = "individual", display_by = "interaction", signaling = "NOTCH", title = "CCCNetworkPlot: individual", verbose = FALSE), 7, 5.5)
add("individual_outgoing", CCCNetworkPlot(srt, method = "LIANA", plot_type = "individual_outgoing", display_by = "interaction", title = "CCCNetworkPlot: individual_outgoing", verbose = FALSE), 7, 5.5)
add("individual_incoming", CCCNetworkPlot(srt, method = "LIANA", plot_type = "individual_incoming", display_by = "interaction", title = "CCCNetworkPlot: individual_incoming", verbose = FALSE), 7, 5.5)
add("arrow", CCCNetworkPlot(srt, method = "LIANA", plot_type = "arrow", display_by = "interaction", top_n = 12, title = "CCCNetworkPlot: arrow", verbose = FALSE), 7, 5)
add("sigmoid", CCCNetworkPlot(srt, method = "LIANA", plot_type = "sigmoid", display_by = "interaction", top_n = 12, title = "CCCNetworkPlot: sigmoid", verbose = FALSE), 7, 5)
add("bipartite", CCCNetworkPlot(srt, method = "LIANA", plot_type = "bipartite", display_by = "interaction", top_n = 12, title = "CCCNetworkPlot: bipartite", verbose = FALSE), 7, 5)
add("diffusion", CCCNetworkPlot(srt, method = "LIANA", plot_type = "diffusion", display_by = "interaction", top_n = 12, title = "CCCNetworkPlot: diffusion", verbose = FALSE), 7, 5)

add("bar", CCCStatPlot(srt, method = "LIANA", plot_type = "bar", display_by = "interaction", top_n = 12, title = "CCCStatPlot: bar", verbose = FALSE), 7.2, 5)
add("sankey", CCCStatPlot(srt, method = "LIANA", plot_type = "sankey", display_by = "interaction", top_n = 12, title = "CCCStatPlot: sankey", verbose = FALSE), 7.2, 5)
add("box", CCCStatPlot(srt, method = "LIANA", plot_type = "box", display_by = "interaction", top_n = 15, title = "CCCStatPlot: box", verbose = FALSE), 7.2, 5)
add("violin", CCCStatPlot(srt, method = "LIANA", plot_type = "violin", display_by = "interaction", top_n = 15, title = "CCCStatPlot: violin", verbose = FALSE), 7.2, 5)
add("role_scatter", CCCStatPlot(srt, method = "LIANA", plot_type = "role_scatter", display_by = "aggregation", title = "CCCStatPlot: role_scatter", verbose = FALSE), 7, 5)
add("pathway_summary", CCCStatPlot(srt, method = "LIANA", plot_type = "pathway_summary", top_n = 12, title = "CCCStatPlot: pathway_summary", verbose = FALSE), 7, 5)
add("comparison", CCCStatPlot(srt, method = "CCC", plot_type = "comparison", display_by = "interaction", sample_col = "method", title = "CCCStatPlot: comparison", verbose = FALSE), 7, 5)
add("lr_contribution", CCCStatPlot(srt, method = "LIANA", plot_type = "lr_contribution", signaling = "NOTCH", title = "CCCStatPlot: lr_contribution", verbose = FALSE), 7, 5)
add("scatter", CCCStatPlot(srt, method = "LIANA", plot_type = "scatter", display_by = "aggregation", title = "CCCStatPlot: scatter", verbose = FALSE), 7, 5)

status <- do.call(rbind, results)
write.csv(status, file.path("scripts", "scop_gallery_status.csv"), row.names = FALSE)
print(status)
cat("Generated", sum(status$ok), "scop-generated figures in", out_dir, "\n")
