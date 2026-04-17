# Combine conserved footprints across species, generate dimensionality-reduction
# views, build cell-type trees, and reconstruct ancestral TF-marker states.
source("../R/csps_minimal.R")
stop_if_missing_packages(c("ape", "phangorn", "umap", "scales", "adephylo", "fastcluster", "stringr"))
graphics.off()

processed_dir = get_env_string("CSPS_PROCESSED_DIR", "../work/processed")
results_root = get_env_string("CSPS_RESULTS_DIR", "../work/results")
results_dir = get_env_string("CSPS_TREE_RESULTS_DIR", "../work/results_cell_type_trees")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

csps_species = get_env_vector("CSPS_SPECIES", c("Ctel", "Ofus"))
csps_ref = get_env_vector("CSPS_REF_SPECIES", c("Ctel"))[1]
comparison_name = get_env_string("CSPS_COMPARISON_NAME", "annelids")
list_comparisons = list()
list_comparisons[[comparison_name]] = list(csps_species, ref = csps_ref)

fc_thr = get_env_numeric("CSPS_FC_THR", 2)
fraction_genes = get_env_numeric("CSPS_FRACTION_GENES", 0.7)
min_ec = get_env_numeric("CSPS_MIN_EC", 0)
num_bs = get_env_integer("CSPS_NUM_BS", 100)
collapse_bs = get_env_numeric("CSPS_COLLAPSE_BS", num_bs * 0.2)
nthreads = get_env_integer("CSPS_THREADS", 4)
focus = "cell_type"
focid = "cts"
options(mc.cores = nthreads)

build_species_color_map = function(species_id, annot_df) {
	raw_names = paste(species_id, annot_df$cell_type, sep = "|")
	sanitised_names = paste(species_id, make.names(annot_df$cell_type), sep = "|")
	color_map = c(setNames(annot_df$color, raw_names), setNames(annot_df$color, sanitised_names))
	color_map[!duplicated(names(color_map))]
}

resolve_tip_colors = function(labels, color_map, fallback = "black") {
	out = unname(color_map[labels])
	names(out) = labels
	missing = is.na(out)
	if (any(missing)) {
		stripped_labels = sub("\\.[0-9]+$", "", labels[missing])
		out[missing] = unname(color_map[stripped_labels])
	}
	out[is.na(out)] = fallback
	out
}

gene_orthology = lapply(csps_species, function(sp) {
	df = read_tsv(file.path(processed_dir, "orthology", sprintf("%s.gene_orthogroups.tsv", sp)))
	dic_from_vecs(df$gene, df$orthogroup)
})
names(gene_orthology) = csps_species

for (nn in seq_along(list_comparisons)) {
	set_id = names(list_comparisons)[nn]
	sps_list = list_comparisons[[nn]][[1]]
	sps_ref = list_comparisons[[nn]][["ref"]]
	sps_query = sps_list[!sps_list %in% sps_ref]

	# Build a shared colour lookup once so every downstream plot uses the same
	# species-aware palette.
	ann_cts = c()
	for (spi in sps_list) {
		ctt_spi = read_tsv(file.path(processed_dir, spi, sprintf("annot.%s.leiden.tsv", spi)))
		ann_spi_v = build_species_color_map(spi, ctt_spi)
		ann_cts = c(ann_cts, ann_spi_v)
	}

	for (spi in sps_query) {
		message(sprintf("csps %s | load %s-%s ICC data...", set_id, sps_ref, spi))
		icc_genes_i = read_tsv(file.path(results_root, sps_ref, "csps", sprintf("dat.icc.%s-%s.ec_scores.tsv", sps_ref, spi)))
		icc_genes_i = icc_genes_i[icc_genes_i$ec_value >= min_ec, , drop = FALSE]
		if (spi == sps_query[1]) {
			icc_genes_f = icc_genes_i[, c("sp1", "sp2")]
		} else {
			icc_genes_f = merge(icc_genes_f, icc_genes_i[, c("sp1", "sp2")], by = "sp1", all = FALSE)
		}
	}
	colnames(icc_genes_f) = c(sps_ref, sps_query)

	# Merge quantile-normalised footprints on the conserved ortholog pairs that
	# survived ICC filtering.
	for (spi in sps_list) {
		message(sprintf("csps %s | load %s expression data...", set_id, spi))
		fps_i = read_matrix_tsv(file.path(processed_dir, spi, sprintf("dat.%s.expression.cts_fp.tsv.gz", spi)))
		fps_i_f = quantile_normalisation(fps_i)
		fps_i_f = fps_i_f[icc_genes_f[, spi], , drop = FALSE]
		colnames(fps_i_f) = sprintf("%s|%s", spi, colnames(fps_i_f))
		if (spi == sps_list[1]) {
			csps_m = fps_i_f
		} else {
			csps_m = cbind(csps_m, fps_i_f)
		}
	}

	csps_m = csps_m[, !grepl("unknown", colnames(csps_m), ignore.case = TRUE), drop = FALSE]
	csps_m = csps_m[rowSums(!is.na(csps_m)) >= ncol(csps_m) * fraction_genes, , drop = FALSE]
	csps_m[is.na(csps_m)] = 0
	# Convert footprints into a simple present/absent matrix of strong enrichment
	# for tree building and dimensionality reduction.
	csps_m_b = (csps_m >= fc_thr) * 1
	csps_m_b_tree = csps_m_b[apply(csps_m_b, 1, function(vv) length(unique(vv)) > 1), , drop = FALSE]
	message(sprintf("csps %s | informative binary genes for tree/bootstrap = %i of %i", set_id, nrow(csps_m_b_tree), nrow(csps_m_b)))
	if (nrow(csps_m_b_tree) == 0) {
		stop(sprintf("csps %s | no informative binary genes remain for tree construction", set_id))
	}

	icc_ecv = read_tsv(file.path(results_root, sps_ref, "csps", sprintf("dat.icc.%s-%s.ec_scores.tsv", sps_ref, sps_list[length(sps_list)])))
	icc_ecv_v = dic_from_vecs(icc_ecv$sp1, icc_ecv$ec_value)
	com = weighted_cor_columns(csps_m, weights = icc_ecv_v[rownames(csps_m)], method = "pearson")

	ctt = data.frame()
	for (spi in sps_list) {
		ctt_i = read_tsv(file.path(processed_dir, spi, sprintf("annot.%s.leiden.tsv", spi)))
		ctt_i$cell_type_sps = paste(spi, ctt_i$cell_type, sep = "|")
		ctt = rbind(ctt, ctt_i)
	}
	dic_cts_col = c()
	for (spi in sps_list) {
		ctt_i = ctt[grepl(sprintf("^%s\\|", spi), paste(spi, ctt$cell_type, sep = "|")), , drop = FALSE]
	}
	dic_cts_col = c()
	for (spi in sps_list) {
		ctt_i = read_tsv(file.path(processed_dir, spi, sprintf("annot.%s.leiden.tsv", spi)))
		dic_cts_col = c(dic_cts_col, build_species_color_map(spi, ctt_i))
	}

	hm1 = plot_complex_heatmap(
		com,
		name = "wpearson",
		color_min = 0,
		color_max = 0.8,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		do_dotplot = FALSE,
		use_raster = FALSE,
		categories_row = rownames(com),
		categories_col = colnames(com),
		colors_row = dic_cts_col,
		colors_col = dic_cts_col,
		cell_border = grid::gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = grid::gpar(col = "black", lwd = 1, lty = 1)
	)
	hm2 = plot_complex_heatmap(
		com,
		name = "wpearson",
		color_min = 0,
		color_max = 0.8,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = TRUE,
		cluster_col = TRUE,
		do_dotplot = FALSE,
		use_raster = FALSE,
		categories_row = rownames(com),
		categories_col = colnames(com),
		colors_row = dic_cts_col,
		colors_col = dic_cts_col,
		cell_border = grid::gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = grid::gpar(col = "black", lwd = 1, lty = 1)
	)
	grDevices::pdf(file.path(results_dir, sprintf("csps.%s.wpearson.hm.%s.pdf", set_id, focid)), width = 8 + ncol(com) / 20, height = 8 + nrow(com) / 20)
	print(hm1)
	print(hm2)
	grDevices::dev.off()

	# UMAP is kept as a qualitative view; PCA additionally feeds a tree built from
	# the retained cell-type loadings.
	grDevices::pdf(file.path(results_dir, sprintf("csps.%s.dimred.%s.pdf", set_id, focid)))
	csps_u = umap::umap(t(csps_m_b))
	plot(csps_u$layout, col = resolve_tip_colors(rownames(csps_u$layout), ann_cts), pch = c(15, 17, 18, 19)[factor(gsub("\\|.*", "", rownames(csps_u$layout)))])
	title(main = "umap")
	text(csps_u$layout, rownames(csps_u$layout), col = scales::alpha(resolve_tip_colors(rownames(csps_u$layout), ann_cts), 0.6), cex = 0.4)

	if (ncol(csps_m_b) < nrow(csps_m_b)) {
		csps_p = stats::princomp(csps_m_b)
		plot(csps_p$loadings[, c(1, 2)], col = resolve_tip_colors(rownames(csps_p$loadings), ann_cts), pch = c(15, 17, 18, 19)[factor(gsub("\\|.*", "", rownames(csps_p$loadings)))])
		title(main = "pca")
		text(csps_p$loadings[, c(1, 2)], rownames(csps_p$loadings), col = scales::alpha(resolve_tip_colors(rownames(csps_p$loadings), ann_cts), 0.6), cex = 0.4)
		pov = csps_p$sdev^2 / sum(csps_p$sdev^2)
		num_pcs_v = find_pca_elbow(pov)
		n_pcs = num_pcs_v["Perpendicular line"]
	}
	grDevices::dev.off()

	plot_height = ceiling(ncol(csps_m_b) / 6 + 4)
	if (ncol(csps_m_b) < nrow(csps_m_b)) {
		# Cluster cell types in the reduced PC space rather than on the full binary
		# matrix to focus on the dominant conserved axes of variation.
		csps_p_f = csps_p$loadings[, 1:n_pcs, drop = FALSE]
		csps_p_f_h = stats::hclust(stats::dist(csps_p_f, method = "manhattan"), method = "average")
		csps_p_f_t = ape::as.phylo(csps_p_f_h)
		grDevices::pdf(file.path(results_dir, sprintf("csps.%s.dendrogram.%s.PCA.pdf", set_id, focid)), height = plot_height, width = 12)
		ape::plot.phylo(ape::ladderize(csps_p_f_t), tip.color = resolve_tip_colors(csps_p_f_t$tip.label, ann_cts), font = 1, underscore = TRUE, type = "phylo", main = sprintf("upgma from %i PCs", n_pcs))
		ape::add.scale.bar()
		grDevices::dev.off()
		ape::write.tree(ape::ladderize(csps_p_f_t), file.path(results_dir, sprintf("csps.%s.dendrogram.%s.PCA.newick", set_id, focid)))
	}

	ali_f = phangorn::as.phyDat(t(csps_m_b_tree), type = "USER", levels = c("0", "1"), names = rownames(csps_m_b_tree), return.index = TRUE)
	parsimony_fun = function(x) reorder(ape::as.phylo(fastcluster::hclust(phangorn::dist.hamming(x), method = "average")), "postorder")
	ali_f_phy = parsimony_fun(ali_f)
	# Bootstrap the binary-gene tree separately so weak clades can be collapsed in
	# the presentation output without discarding the fully resolved topology.
	boo_list = tryCatch(
		phangorn::bootstrap.phyDat(ali_f, parsimony_fun, multicore = FALSE, bs = num_bs),
		error = function(e) {
			message(sprintf("csps %s | bootstrap failed, continue without support values: %s", set_id, conditionMessage(e)))
			NULL
		}
	)
	if (!is.null(boo_list)) {
		valid_bootstrap = vapply(boo_list, function(x) inherits(x, "phylo"), logical(1))
		if (!all(valid_bootstrap)) {
			message(sprintf("csps %s | dropping %i invalid bootstrap replicates", set_id, sum(!valid_bootstrap)))
			boo_list = boo_list[valid_bootstrap]
		}
		if (length(boo_list) == 0) {
			message(sprintf("csps %s | no valid bootstrap replicates remain, continue without support values", set_id))
			boo_list = NULL
		}
	}
	ali_f_phy_col = ali_f_phy
	if (!is.null(boo_list)) {
		collapse_nodes = which(ape::prop.clades(ali_f_phy, boo_list) < collapse_bs) + length(ali_f_phy$tip.label)
		rownames(ali_f_phy$edge) = as.character(seq_len(nrow(ali_f_phy$edge)))
		collapse_edges = sort(as.numeric(rownames(ali_f_phy$edge)[ali_f_phy$edge[, 2] %in% collapse_nodes]))
		ali_f_phy_col$edge.length[collapse_edges] = 0
		ali_f_phy_col = ape::di2multi(ali_f_phy_col)
	}
	ali_f_phy = ape::makeNodeLabel(ali_f_phy)
	descendants_from_nodes = adephylo::listTips(ali_f_phy)
	species_palette = grDevices::hcl.colors(length(sps_list), palette = "Dark 3")
	names(species_palette) = sps_list
	tip_species = sub("\\|.*", "", ali_f_phy$tip.label)
	tip_cell_type = sub("^[^|]+\\|", "", ali_f_phy$tip.label)
	# Export a supervisor-friendly editable tip table so downstream figure scripts
	# do not have to recover metadata from the tree object itself.
	tip_metadata = data.frame(
		leaf_id = ali_f_phy$tip.label,
		species_abbrev = tip_species,
		cell_state_label = tip_cell_type,
		leaf_display_label = ali_f_phy$tip.label,
		tip_color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts),
		species_color = unname(species_palette[tip_species]),
		clade_assignment = "",
		stringsAsFactors = FALSE
	)
	write_tsv(tip_metadata, file.path(results_dir, sprintf("csps.%s.dendrogram.%s.tip_metadata.tsv", set_id, focid)))
	# Store support values and descendants per node in plain TSV form so the tree
	# can be re-annotated later without re-running the bootstrap step.
	node_support = data.frame(
		node_id = ali_f_phy$node.label,
		support = if (!is.null(boo_list)) ape::prop.clades(ali_f_phy, boo_list) else NA_real_,
		species_represented = unlist(lapply(descendants_from_nodes[ali_f_phy$node.label], function(vv) paste(sort(unique(sub("\\|.*", "", names(vv)))), collapse = ","))),
		tips_in_node = unlist(lapply(descendants_from_nodes[ali_f_phy$node.label], function(vv) paste(sort(names(vv)), collapse = ","))),
		stringsAsFactors = FALSE
	)
	write_tsv(node_support, file.path(results_dir, sprintf("csps.%s.dendrogram.%s.node_support.tsv", set_id, focid)))

	grDevices::pdf(file.path(results_dir, sprintf("csps.%s.dendrogram.%s.UPGMA.pdf", set_id, focid)), height = plot_height, width = 12)
	if (is.null(boo_list)) {
		ape::plot.phylo(ali_f_phy, tip.color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts), font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f", fc_thr), root.edge = TRUE)
		ape::add.scale.bar()
		ape::plot.phylo(ali_f_phy_col, tip.color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts), font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f", fc_thr), root.edge = TRUE, use.edge.length = FALSE)
	} else {
		phangorn::plotBS(ali_f_phy, boo_list, tip.color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts), font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP")
		ape::add.scale.bar()
		phangorn::plotBS(ali_f_phy_col, boo_list, tip.color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts), font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP", use.edge.length = FALSE)
	}
	grDevices::dev.off()
	ape::write.tree(ali_f_phy, file.path(results_dir, sprintf("csps.%s.dendrogram.%s.UPGMA.newick", set_id, focid)))

	# Restrict ancestral-state reconstruction to confident TF markers so the
	# output focuses on interpretable regulatory programs.
	list_markers_tfs = data.frame()
	for (spi in sps_list) {
		mks_i = read_tsv(file.path(processed_dir, spi, sprintf("seu.%s.markers.cts.tsv.gz", spi)))
		mks_i = mks_i[!is.na(mks_i$gene) & mks_i$p_val_adj < 0.01 & mks_i$avg_log2FC > 0 & mks_i$is_tf, , drop = FALSE]
		if (!"annotation" %in% colnames(mks_i)) {
			mks_i$annotation = ""
		}
		focus_nodes = unique(mks_i$focus_node)
		for (foi in focus_nodes) {
			mks_ii = unique(mks_i[mks_i$focus_node == foi, c("gene", "annotation"), drop = FALSE])
			if (nrow(mks_ii) > 0) {
				mks_ii$orthogroup = gene_orthology[[spi]][mks_ii$gene]
				list_markers_tfs = rbind(
					list_markers_tfs,
					data.frame(
						focus_node = sprintf("%s|%s", spi, foi),
						orthogroup = mks_ii$orthogroup,
						gene = mks_ii$gene,
						annotation = mks_ii$annotation,
						stringsAsFactors = FALSE
					)
				)
			}
		}
	}

	if (nrow(list_markers_tfs) > 0) {
		list_markers_tfs = list_markers_tfs[!is.na(list_markers_tfs$orthogroup), , drop = FALSE]
		orthogroup_labels = list_markers_tfs[!is.na(list_markers_tfs$annotation) & nzchar(list_markers_tfs$annotation), c("orthogroup", "annotation"), drop = FALSE]
		orthogroup_labels = orthogroup_labels[!duplicated(orthogroup_labels$orthogroup), , drop = FALSE]
		orthogroup_label_v = if (nrow(orthogroup_labels) > 0) {
			dic_from_vecs(orthogroup_labels$orthogroup, orthogroup_labels$annotation)
		} else {
			c()
		}
		list_markers_tfs$presence = 1
		list_markers_tfs_m = xtabs(presence ~ focus_node + orthogroup, data = list_markers_tfs)
		list_markers_tfs_m_f = list_markers_tfs_m[, apply(list_markers_tfs_m, 2, function(vv) sum(vv > 0) >= 2), drop = FALSE]
		missing_tip_rows = setdiff(ali_f_phy$tip.label, rownames(list_markers_tfs_m_f))
		if (length(missing_tip_rows) > 0) {
			empty_rows = matrix(0, nrow = length(missing_tip_rows), ncol = ncol(list_markers_tfs_m_f), dimnames = list(missing_tip_rows, colnames(list_markers_tfs_m_f)))
			list_markers_tfs_m_f = rbind(list_markers_tfs_m_f, empty_rows)
		}
		list_markers_tfs_m_f = list_markers_tfs_m_f[ali_f_phy$tip.label, , drop = FALSE]
		ali_f_phy_sansnodes = ali_f_phy
		ali_f_phy_sansnodes$node.label = NULL
		grDevices::pdf(file.path(results_dir, sprintf("csps.%s.dendrogram.%s.UPGMA.markers.pdf", set_id, focid)), height = plot_height, width = 24)
		markers_per_node = data.frame()
		for (nnn in seq_len(ncol(list_markers_tfs_m_f))) {
			ggo = colnames(list_markers_tfs_m_f)[nnn]
			orthogroup_name = orthogroup_label_v[ggo]
			if (is.na(orthogroup_name) || orthogroup_name == "") {
				orthogroup_name = ggo
			}
			vv_i = factor(as.character(list_markers_tfs_m_f[ali_f_phy$tip.label, ggo]), levels = c("0", "1"))
			# Fit a simple discrete-character model per orthogroup and record the
			# posterior probability that each internal node carries the marker.
			ace_i = ape::ace(vv_i, ali_f_phy_sansnodes, type = "discrete", method = "ML", marginal = TRUE, model = "ER")
			ace_i_d = ace_i$lik.anc
			rownames(ace_i_d) = ali_f_phy$node.label
				vv_i_m = data.frame(row.names = ali_f_phy$tip.label, "0" = as.numeric(vv_i == 0), "1" = as.numeric(vv_i == 1))
			ape::plot.phylo(ali_f_phy, tip.color = resolve_tip_colors(ali_f_phy$tip.label, ann_cts), font = 1, underscore = TRUE, type = "phylo", root.edge = TRUE, use.edge.length = TRUE, show.node.label = TRUE)
			ape::nodelabels(node = 1:ali_f_phy$Nnode + ape::Ntip(ali_f_phy), pie = ace_i_d, piecol = c("honeydew2", "blue3"), cex = 0.2)
			ape::tiplabels(pie = vv_i_m, piecol = c("honeydew2", "blue3"), cex = 0.1)
			title(main = ggo)
			ace_i_d_i = data.frame(
				orthogroup = ggo,
				orthogroup_name = orthogroup_name,
				# Posterior probability that the internal node carries the marker state.
				probability = ace_i_d[, "1"],
				# Threshold the posterior into a simple yes/no call for quick filtering.
				present = ace_i_d[, "1"] >= 0.7,
				node = rownames(ace_i_d),
				tips_in_node = unlist(lapply(descendants_from_nodes[rownames(ace_i_d)], function(vv) paste(sort(names(vv)), collapse = ",")))
			)
			markers_per_node = rbind(markers_per_node, ace_i_d_i)
		}
		grDevices::dev.off()
		# This table is the machine-readable counterpart of the node-pie PDF above:
		# one row per orthogroup-by-node ancestral-state reconstruction result.
		write_tsv(markers_per_node, file.path(results_dir, sprintf("csps.%s.dendrogram.%s.UPGMA.markers.tsv", set_id, focid)))
		node_markers_split = split(markers_per_node, markers_per_node$node)
		clade_program_template = node_support
		clade_program_template$clade_key = clade_program_template$node_id
		clade_program_template$clade_name = ""
		clade_program_template$label_color = ""
		clade_program_template$show = FALSE
		clade_program_template$anchor_node_id = clade_program_template$node_id
		clade_program_template$tf_text = ""
		clade_program_template$suggested_tfs = vapply(clade_program_template$node_id, function(node_id) {
			node_df = node_markers_split[[node_id]]
			if (is.null(node_df)) {
				return("")
			}
			node_df = node_df[node_df$present, , drop = FALSE]
			if (nrow(node_df) == 0) {
				return("")
			}
			paste(utils::head(unique(node_df[order(node_df$probability, decreasing = TRUE), "orthogroup_name"]), 10), collapse = ", ")
		}, character(1))
		# Pre-populate the rendered TF-program column so the default outputs already
		# show candidate markers without requiring manual backfilling.
		clade_program_template$tf_text = clade_program_template$suggested_tfs
		clade_program_template$clade_name = ifelse(nzchar(clade_program_template$suggested_tfs), clade_program_template$node_id, "")
		clade_program_template$show = nzchar(clade_program_template$suggested_tfs)
		write_tsv(clade_program_template, file.path(results_dir, sprintf("csps.%s.dendrogram.%s.tf_program_template.tsv", set_id, focid)))

		ali_f_phy_is_tip = ali_f_phy$edge[, 2] <= length(ali_f_phy$tip.label)
		ali_f_phy_ordered_tips = ali_f_phy$tip.label[ali_f_phy$edge[ali_f_phy_is_tip, 2]]
		genes_for_plot = intersect(unique(list_markers_tfs$gene), rownames(csps_m))
		csps_m_order_f = csps_m[genes_for_plot, ali_f_phy_ordered_tips, drop = FALSE]
		expr_m = csps_m_order_f
		# Order TFs by where their smoothed peak footprint lies along the tree so the
		# final heatmap reads as a lineage-ordered marker panel.
		top_markers = select_top_markers(expr_m, matrix_thr = 1.5, n_top_markers = 1e6, n_markers_rollmean = 4)
		expr_m = expr_m[top_markers, , drop = FALSE]
		top_cl_per_gene = colnames(expr_m)[apply(expr_m, 1, which.max)]
		rownames(expr_m) = sprintf("%s | %s", rownames(expr_m), gene_orthology[[sps_ref]][rownames(expr_m)])
		rownames(expr_m) = stringr::str_trunc(rownames(expr_m), 100)
		pp1 = plot_complex_heatmap(
			expr_m,
			name = "Footprint",
			cluster_row = FALSE,
			cluster_col = FALSE,
			use_raster = FALSE,
			# Make the lower-but-still-meaningful 1.1-2 band visually legible rather
			# than letting it collapse into the background of the linear 1-4 scale.
			color_breaks = c(0, 1.099, 1.1, 2, 4),
			color_max = 4,
			color_mat = c("gray96", "gray96", "#ffe08a", "orange", "#520c52"),
			cell_border = grid::gpar(col = "white", lwd = 1, lty = 1),
			heatmap_border = grid::gpar(col = "black", lwd = 1, lty = 1),
			categories_row = top_cl_per_gene,
			categories_col = colnames(expr_m),
			colors_row = ann_cts,
			colors_col = ann_cts,
			fontsize = 6
		)
		pdfheight = max(8, (nrow(expr_m) / 8 + 6))
		grDevices::pdf(file.path(results_dir, sprintf("csps.%s.dendrogram.%s.UPGMA.ordered_heatmap_TFs.pdf", set_id, focid)), width = 20, height = pdfheight)
		print(pp1)
		grDevices::dev.off()
	}

	saveRDS(csps_m, file.path(results_dir, sprintf("csps.%s.cspsmatrix.%s.rds", set_id, focid)))
}
