# Shared helpers extracted from the larger coral pipeline for the annelid-only run.
stop_if_missing_packages = function(pkgs) {
	missing = pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
	if (length(missing) > 0) {
		stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")))
	}
}

# Environment helpers keep the stage scripts portable between local and cluster runs.
get_env_string = function(name, default = NULL) {
	value = Sys.getenv(name, unset = "")
	if (nzchar(value)) {
		value
	} else {
		default
	}
}

get_env_integer = function(name, default) {
	value = get_env_string(name, as.character(default))
	as.integer(value)
}

get_env_numeric = function(name, default) {
	value = get_env_string(name, as.character(default))
	as.numeric(value)
}

get_env_flag = function(name, default = FALSE) {
	value = get_env_string(name, if (default) "TRUE" else "FALSE")
	tolower(value) %in% c("1", "true", "t", "yes", "y")
}

get_env_vector = function(name, default, sep = ",") {
	value = get_env_string(name, NULL)
	if (is.null(value) || !nzchar(value)) {
		return(default)
	}
	parts = trimws(strsplit(value, sep, fixed = TRUE)[[1]])
	parts[nzchar(parts)]
}

dic_from_vecs = function(names, terms, make_unique_on_names = TRUE) {
	dic = terms
	names(dic) = names
	if (make_unique_on_names) {
		dic = dic[!duplicated(names(dic))]
	}
	dic
}

read_tsv = function(path) {
	stop_if_missing_packages("data.table")
	data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
}

read_matrix_tsv = function(path) {
	df = read_tsv(path)
	rownames(df) = df[[1]]
	as.matrix(df[, -1, drop = FALSE])
}

# Weighted Pearson is used throughout the pipeline to let ICC-conserved genes
# contribute more strongly than weakly conserved orthologs.
weighted_cor_columns = function(X, Y = NULL, weights = NULL, method = "pearson") {
	X = as.matrix(X)
	if (is.null(Y)) {
		Y = X
	}
	Y = as.matrix(Y)
	if (is.null(weights)) {
		return(stats::cor(X, Y, method = method, use = "pairwise.complete.obs"))
	}
	if (method != "pearson") {
		stop("Weighted correlations are only implemented for method='pearson'.")
	}
	w = as.numeric(weights)
	keep = is.finite(w) & w > 0
	X = X[keep, , drop = FALSE]
	Y = Y[keep, , drop = FALSE]
	w = w[keep]
	if (length(w) == 0) {
		stop("No positive finite weights remain for weighted correlation.")
	}
	w_sum = sum(w)
	mx = as.numeric(crossprod(w, X) / w_sum)
	my = as.numeric(crossprod(w, Y) / w_sum)
	Xc = sweep(X, 2, mx, "-")
	Yc = sweep(Y, 2, my, "-")
	cov_xy = crossprod(Xc * w, Yc) / w_sum
	sx = sqrt(colSums((Xc^2) * w) / w_sum)
	sy = sqrt(colSums((Yc^2) * w) / w_sum)
	out = cov_xy / outer(sx, sy)
	out[!is.finite(out)] = 0
	rownames(out) = colnames(X)
	colnames(out) = colnames(Y)
	out
}

write_tsv = function(x, path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Centralise ComplexHeatmap setup so the analysis stages only need to pass data
# and annotation vectors.
plot_complex_heatmap = function(
	mat,
	name = "heatmap",
	color_mat = c("white", "#d6e72e", "#6fb600", "#003f4d"),
	color_breaks = NULL,
	color_min = 0,
	color_max = 1,
	fontsize = 10,
	hm_body_height = NULL,
	hm_body_width = NULL,
	categories_col = NULL,
	categories_row = NULL,
	colors_col = NULL,
	colors_row = NULL,
	title_row = NULL,
	title_col = NULL,
	name_row_show = TRUE,
	name_col_show = TRUE,
	cluster_row = TRUE,
	cluster_col = TRUE,
	use_raster = TRUE,
	show_legend_row = FALSE,
	show_legend_col = FALSE,
	show_legend_hm = TRUE,
	both_sides_row = TRUE,
	both_sides_col = TRUE,
	cell_border = grid::gpar(col = NA, lwd = 1, lty = 1),
	heatmap_border = grid::gpar(col = NA, lwd = 1, lty = 1),
	do_dotplot = FALSE,
	dot_size_mat = NULL,
	dot_size_min = NULL,
	dot_size_max = NULL,
	cex_dotplot = 0.02,
	col_dotplot_border = NA
) {
	stop_if_missing_packages(c("ComplexHeatmap", "circlize"))
	if (is.null(color_breaks)) {
		color_breaks = seq(color_min, color_max, length.out = length(color_mat))
	}
	if (length(color_breaks) != length(color_mat)) {
		stop("plot_complex_heatmap: color_breaks must have the same length as color_mat")
	}
	col_fun = circlize::colorRamp2(color_breaks, color_mat)
	if (is.null(title_row)) title_row = sprintf("n = %i rows", nrow(mat))
	if (is.null(title_col)) title_col = sprintf("n = %i columns", ncol(mat))
	if (is.null(categories_row)) categories_row = rownames(mat)
	if (is.null(categories_col)) categories_col = colnames(mat)

	right_annotation = NULL
	top_annotation = NULL
	if (!is.null(colors_row)) {
		if (is.null(names(colors_row))) names(colors_row) = unique(categories_row)
		missing_row_levels = setdiff(unique(categories_row), names(colors_row))
		if (length(missing_row_levels) > 0) {
			fallback = rep("gray80", length(missing_row_levels))
			names(fallback) = missing_row_levels
			colors_row = c(colors_row, fallback)
		}
		right_annotation = ComplexHeatmap::HeatmapAnnotation(
			c = categories_row,
			name = title_row,
			col = list(c = colors_row),
			which = "row",
			show_legend = show_legend_row,
			show_annotation_name = FALSE
		)
	}
	if (!is.null(colors_col)) {
		if (is.null(names(colors_col))) names(colors_col) = unique(categories_col)
		missing_col_levels = setdiff(unique(categories_col), names(colors_col))
		if (length(missing_col_levels) > 0) {
			fallback = rep("gray80", length(missing_col_levels))
			names(fallback) = missing_col_levels
			colors_col = c(colors_col, fallback)
		}
		top_annotation = ComplexHeatmap::HeatmapAnnotation(
			c = categories_col,
			name = title_col,
			col = list(c = colors_col),
			which = "column",
			show_legend = show_legend_col,
			show_annotation_name = FALSE
		)
	}

	left_annotation = if (both_sides_row) right_annotation else NULL
	bottom_annotation = if (both_sides_col) top_annotation else NULL

	cell_fun_dotplot = NULL
	if (do_dotplot) {
		if (is.null(dot_size_mat)) dot_size_mat = mat
		rangesize = function(x) {
			if (!is.null(dot_size_max) && !is.null(dot_size_min)) {
				(x - dot_size_min) / (dot_size_max - dot_size_min)
			} else {
				(x - min(x)) / (max(x) - min(x))
			}
		}
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			grid::grid.circle(
				x = x,
				y = y,
				r = sqrt(rangesize(dot_size_mat)[i, j]) * cex_dotplot,
				gp = grid::gpar(col = col_dotplot_border, fill = col_fun(mat[i, j]))
			)
		}
	}

	ComplexHeatmap::Heatmap(
		mat,
		name = name,
		cell_fun = cell_fun_dotplot,
		height = hm_body_height,
		width = hm_body_width,
		use_raster = use_raster,
		cluster_rows = cluster_row,
		cluster_columns = cluster_col,
		row_title = title_row,
		row_title_gp = grid::gpar(fontsize = fontsize),
		column_title = title_col,
		column_title_gp = grid::gpar(fontsize = fontsize),
		show_row_names = name_row_show,
		show_column_names = name_col_show,
		row_names_gp = grid::gpar(fontsize = fontsize),
		column_names_gp = grid::gpar(fontsize = fontsize),
		top_annotation = top_annotation,
		left_annotation = left_annotation,
		right_annotation = right_annotation,
		bottom_annotation = bottom_annotation,
		row_gap = grid::unit(0.5, "mm"),
		column_gap = grid::unit(0.5, "mm"),
		rect_gp = cell_border,
		border_gp = heatmap_border,
		show_heatmap_legend = show_legend_hm,
		col = col_fun
	)
}

# Combine several common elbow heuristics so stage scripts can choose a single
# rule without reimplementing the full diagnostic logic.
find_pca_elbow = function(sdev) {
	x = 1:length(sdev)
	eb = data.frame(x, sdev)
	sse = NULL
	for (k in 1:length(sdev)) {
		D = ifelse(x <= k, 0, 1)
		x2 = (x - k) * D
		sse[k] = sum(lm(sdev ~ x + x2)$residuals^2)
	}
	dim_plm = which.min(sse)
	names(dim_plm) = "Piecewise linear model"
	df1 = diff(sdev, differences = 1)
	den = density(df1)
	rlev = rle(diff(den$y) > 0)$lengths
	if (length(rlev) <= 2) {
		dim_fid = max(which(df1 < mean(df1))) + 1
	} else {
		cutoff = sum(rlev[-((length(rlev) - 1):length(rlev))])
		dim_fid = max(which(df1 < den$x[cutoff])) + 1
	}
	names(dim_fid) = "First derivative"
	df2 = diff(sdev, differences = 2)
	df2p = df2[df2 > 0]
	den = density(df2p)
	rlev = rle(diff(den$y) > 0)$lengths
	if (length(rlev) <= 2) {
		dim_sed = which.max(df2) + 1
	} else {
		cutoff = sum(rlev[1:2])
		dim_sed = max(which(df2 > den$x[cutoff])) + 1
	}
	names(dim_sed) = "Second derivative"
	fit = NULL
	res = NULL
	for (i in 1:(length(sdev) - 2)) {
		fit[[i]] = lm(sdev ~ x, data = eb[(i + 1):length(sdev), ])
		res[i] = sdev[i] - predict(fit[[i]], newdata = data.frame(x = i))
	}
	den = density(res)
	rlev = rle(diff(den$y) > 0)$lengths
	if (length(rlev) <= 2) {
		dim_pr = max(which(res > mean(res))) + 1
	} else {
		cutoff = sum(rlev[1:2])
		dim_pr = max(which(res > den$x[cutoff])) + 1
	}
	names(dim_pr) = "Preceding residual"
	A = c(1, sdev[1])
	B = c(length(sdev), sdev[length(sdev)])
	Dist = NULL
	for (i in 2:(length(sdev) - 1)) {
		C = c(i, sdev[i])
		D = cbind(rbind(A, B, C), rep(1, 3))
		S = 0.5 * abs(det(D))
		Dist[i] = 2 * S / dist(rbind(A, B))
	}
	dim_perl = which.max(Dist)
	names(dim_perl) = "Perpendicular line"
	set.seed(2022)
	dim_clu = min(stats::kmeans(sdev, 2)$size) + 1
	names(dim_clu) = "K-means clustering"
	c(dim_plm, dim_fid, dim_sed, dim_pr, dim_perl, dim_clu)
}

select_top_markers = function(matrix, matrix_thr = 0, n_top_markers = 20, n_markers_rollmean = 2) {
	stop_if_missing_packages("zoo")
	markers = unique(as.vector(unlist(apply(matrix, 2, function(c) names(head(sort(-c[c >= matrix_thr]), n = n_top_markers))))))
	markers_order = order(apply(matrix[markers, , drop = FALSE], 1, function(r) which.max(zoo::rollmean(r, n_markers_rollmean))))
	markers[markers_order]
}

# Quantile normalisation makes the per-species footprint matrices comparable
# before they are merged into one cross-species matrix.
quantile_normalisation = function(x) {
	df_rank = data.frame(apply(x, 2, rank, ties.method = "min"))
	df_sorted = data.frame(apply(x, 2, sort))
	df_mean = apply(df_sorted, 1, mean)
	index_to_mean = function(my_index, my_mean) my_mean[my_index]
	df_final = apply(df_rank, 2, index_to_mean, my_mean = df_mean)
	rownames(df_final) = rownames(x)
	df_final
}

# Keep genes that show strong enrichment in both species so the similarity
# matrix is driven by shared, informative variation.
csps_select_covariable_genes = function(sp1_fp, sp2_fp, merged, method = "min_fc", cross_fc_thrs = 2, cross_n = 1) {
	merged_sp1 = merged[, 1:ncol(sp1_fp), drop = FALSE]
	merged_sp2 = merged[, (ncol(sp1_fp) + 1):ncol(merged), drop = FALSE]
	if (method == "min_fc") {
		top_cross_ix = which(
			apply(merged_sp1, 1, function(x) sort(x, decreasing = TRUE)[cross_n]) > cross_fc_thrs &
			apply(merged_sp2, 1, function(x) sort(x, decreasing = TRUE)[cross_n]) > cross_fc_thrs
		)
		top_cross_sp1 = rownames(sp1_fp)[top_cross_ix]
		top_cross_sp2 = rownames(sp2_fp)[top_cross_ix]
	} else {
		stop("Only method='min_fc' is supported in the annelid adaptation.")
	}
	list(sp1 = top_cross_sp1, sp2 = top_cross_sp2)
}

# Build the cross-species cell-type correlation matrix from a matched gene set.
csps_correlation_matrix_noobj = function(mm, m1, m2, use_var_genes = FALSE, cor_method = "wpearson", gene_weights = NULL, add_sps_prefix = TRUE, prefix_sp1 = NULL, prefix_sp2 = NULL) {
	gene_names = rownames(m1)
	rownames(mm) = gene_names
	rownames(m1) = gene_names
	rownames(m2) = gene_names
	if (add_sps_prefix) {
		if (is.null(prefix_sp1)) prefix_sp1 = "sp1"
		if (is.null(prefix_sp2)) prefix_sp2 = "sp2"
		colnames(m1) = paste(prefix_sp1, colnames(m1), sep = "|")
		colnames(m2) = paste(prefix_sp2, colnames(m2), sep = "|")
		colnames(mm) = c(colnames(m1), colnames(m2))
	}
	if (is.logical(use_var_genes) && any(use_var_genes == FALSE)) {
		var_genes = rownames(mm)
	} else if (is.character(use_var_genes)) {
		var_genes = use_var_genes
	} else {
		stop("use_var_genes has to be FALSE or a character vector")
	}
	m1_f = m1[var_genes, , drop = FALSE]
	m2_f = m2[var_genes, , drop = FALSE]
	if (cor_method == "wpearson") {
		if (is.null(gene_weights) || length(gene_weights) != length(var_genes)) {
			stop("Weighted Pearson requires gene_weights matching var_genes")
		}
		com = weighted_cor_columns(m1_f, m2_f, gene_weights, method = "pearson")
	} else {
		com = stats::cor(m1_f, m2_f, method = cor_method)
	}
	list(cor_matrix = com, overlap_matrix = NULL, overlapping_genes = NULL)
}

# Iteratively reweight coexpression conservation so genes supported by a stable
# conserved neighbourhood keep influence while weak matches are damped down.
csps_calc_icc = function(mat_sp1, mat_sp2, niter = 100, icc_thr = 0.05, method = "pearson", verbose = TRUE, num_cores = 2) {
	if (verbose) message(sprintf("ICC init coexpression correlation matrices, 1st species (%i x %i)", nrow(mat_sp1), ncol(mat_sp1)))
	exc_sp1 = stats::cor(t(mat_sp1), method = method, use = "pairwise.complete.obs")
	if (verbose) message(sprintf("ICC init coexpression correlation matrices, 2nd species (%i x %i)", nrow(mat_sp2), ncol(mat_sp2)))
	exc_sp2 = stats::cor(t(mat_sp2), method = method, use = "pairwise.complete.obs")
	i = 1
	if (verbose) message(sprintf("ICC iteration %i (%i x %i | %i x %i), run", i - 1, nrow(exc_sp1), ncol(exc_sp1), nrow(exc_sp2), ncol(exc_sp2)))
	exc_csp = stats::cor(exc_sp1, exc_sp2, method = method, use = "pairwise.complete.obs")
	if (verbose) message(sprintf("ICC iteration %i (%i x %i | %i x %i), done", i - 1, nrow(exc_sp1), ncol(exc_sp1), nrow(exc_sp2), ncol(exc_sp2)))
	icc_ecv = list()
	icc_ecv[[i]] = diag(exc_csp)
	exp_sp1 = exc_sp1
	exp_sp2 = exc_sp2
	for (i in 2:niter) {
		ixs = which(icc_ecv[[i - 1]] > 0)
		exi_sp1 = matrix(0, nrow = nrow(exp_sp1), ncol = ncol(exp_sp1))
		exi_sp1[ixs, ixs] = exp_sp1[ixs, ixs]
		exi_sp2 = matrix(0, nrow = nrow(exp_sp2), ncol = ncol(exp_sp2))
		exi_sp2[ixs, ixs] = exp_sp2[ixs, ixs]
		exi_weights = rep(0, length.out = length(icc_ecv[[i - 1]]))
		exi_weights[ixs] = icc_ecv[[i - 1]][ixs]
		if (verbose) message(sprintf("ICC iteration %i (%i x %i | %i x %i), run with weights", i - 1, nrow(exi_sp1), ncol(exi_sp1), nrow(exi_sp2), ncol(exi_sp2)))
		icc_eci = weighted_cor_columns(exi_sp1, exi_sp2, exi_weights, method = method)
		icc_ecv[[i]] = diag(icc_eci)
		icc_delta = sum((icc_ecv[[i]][ixs] - icc_ecv[[i - 1]][ixs])^2)
		if (verbose) message(sprintf("ICC iteration %i (%i x %i | %i x %i), done with delta = %.1e", i - 1, nrow(exi_sp1), ncol(exi_sp1), nrow(exi_sp2), ncol(exi_sp2), icc_delta))
		exp_sp1 = exi_sp1
		exp_sp2 = exi_sp2
		if (icc_delta < icc_thr) break
	}
	data.frame(sp1 = rownames(exc_sp1), sp2 = rownames(exc_sp2), ec_value = icc_ecv[[i]])
}

# Start from one-to-one orthologs, then optionally rescue duplicates by scoring
# each ambiguous ortholog block against the conserved reference set.
csps_markers_icc_noobj = function(mat_sp1, mat_sp2, og_pairs, niter = 100, icc_thr = 0.05, method = "pearson", do_duplicates = FALSE, use_variable_o2o_genes = TRUE, variable_o2o_thr = 0.75, nthreads_icc = 2) {
	stop_if_missing_packages(c("igraph", "stringr"))
	mat_sp1 = as.matrix(mat_sp1)
	mat_sp2 = as.matrix(mat_sp2)
	colnames(og_pairs) = c("sp1", "sp2")
	og_pairs = og_pairs[og_pairs[, 1] %in% rownames(mat_sp1) & og_pairs[, 2] %in% rownames(mat_sp2), , drop = FALSE]
	list_o2o_sp1 = names(which(table(og_pairs[, 1]) == 1))
	list_o2o_sp2 = names(which(table(og_pairs[, 2]) == 1))
	bool_o2o = og_pairs[, 1] %in% list_o2o_sp1 & og_pairs[, 2] %in% list_o2o_sp2
	og_pairs_o2o = og_pairs[bool_o2o, , drop = FALSE]
	mar_sp1 = mat_sp1[og_pairs_o2o[, 1], , drop = FALSE]
	mar_sp2 = mat_sp2[og_pairs_o2o[, 2], , drop = FALSE]
	message(sprintf("ICC markers | pairs of one-to-one orthologs = %i", nrow(og_pairs_o2o)))
	ecv_o2o = csps_calc_icc(mat_sp1 = mar_sp1, mat_sp2 = mar_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = TRUE, num_cores = nthreads_icc)
	ecv_o2o$ec_value[ecv_o2o$ec_value < 0 | is.na(ecv_o2o$ec_value)] = 0
	ecv_o2o$is_o2o = 1
	ecv_dup = NULL
	ecv_out = ecv_o2o
	if (do_duplicates) {
		if (use_variable_o2o_genes) {
			max_fcs_sp1 = apply(mat_sp1, 1, max, na.rm = TRUE)
			max_fcs_sp2 = apply(mat_sp2, 1, max, na.rm = TRUE)
			vargenes_sp1 = rownames(mat_sp1)[max_fcs_sp1 > stats::quantile(max_fcs_sp1, variable_o2o_thr)]
			vargenes_sp2 = rownames(mat_sp2)[max_fcs_sp2 > stats::quantile(max_fcs_sp2, variable_o2o_thr)]
			ecv_o2o_ref = ecv_o2o[ecv_o2o[, 1] %in% vargenes_sp1 & ecv_o2o[, 2] %in% vargenes_sp2, , drop = FALSE]
			refgenes_sp1 = ecv_o2o_ref[, 1]
			refgenes_sp2 = ecv_o2o_ref[, 2]
		} else {
			refgenes_sp1 = ecv_o2o[, 1]
			refgenes_sp2 = ecv_o2o[, 2]
		}
		if (length(refgenes_sp1) == 0) stop("No one-to-one orthologs available for duplicate EC scoring.")
		list_dup_sp1 = names(which(table(og_pairs[, 1]) > 1))
		list_dup_sp2 = names(which(table(og_pairs[, 2]) > 1))
		bool_dup = og_pairs[, 1] %in% list_dup_sp1 | og_pairs[, 2] %in% list_dup_sp2
		og_pairs_dup = og_pairs[bool_dup, , drop = FALSE]
		if (nrow(og_pairs_dup) > 0) {
			og_pairs_dup_sps = as.matrix(og_pairs_dup)
			og_pairs_dup_sps[, 1] = paste("sp1", og_pairs_dup_sps[, 1])
			og_pairs_dup_sps[, 2] = paste("sp2", og_pairs_dup_sps[, 2])
			graph_dup = igraph::graph_from_edgelist(og_pairs_dup_sps, directed = TRUE)
			graph_dup_components = igraph::components(graph_dup)
			dict_dups = split(names(graph_dup_components$membership), graph_dup_components$membership)
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes)", length(dict_dups), length(unique(unlist(dict_dups)))))
			dup_scores = list()
			for (idx in seq_along(dict_dups)) {
				d = dict_dups[[idx]]
				dg = stringr::word(d, 2)
				genes1 = dg[stringr::word(d, 1) == "sp1"]
				genes2 = dg[stringr::word(d, 1) == "sp2"]
				gen_dups = og_pairs[og_pairs[, 1] %in% genes1 & og_pairs[, 2] %in% genes2, , drop = FALSE]
				mai_sp1 = rbind(mat_sp1[gen_dups[, 1], , drop = FALSE], mar_sp1[refgenes_sp1, , drop = FALSE])
				rownames(mai_sp1)[seq_len(nrow(gen_dups))] = gen_dups[, 1]
				mai_sp2 = rbind(mat_sp2[gen_dups[, 2], , drop = FALSE], mar_sp2[refgenes_sp2, , drop = FALSE])
				rownames(mai_sp2)[seq_len(nrow(gen_dups))] = gen_dups[, 2]
				ecv_dup_i = csps_calc_icc(mat_sp1 = mai_sp1, mat_sp2 = mai_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = FALSE, num_cores = nthreads_icc)
				ecv_dup_q = ecv_dup_i[ecv_dup_i$sp1 %in% gen_dups[, 1] & ecv_dup_i$sp2 %in% gen_dups[, 2], , drop = FALSE]
				ecv_dup_q$cluster = as.character(idx)
				dup_scores[[idx]] = ecv_dup_q
			}
			ecv_dup = do.call(rbind, dup_scores)
			ecv_dup$ec_value[is.na(ecv_dup$ec_value)] = -1
			best_sp1 = ave(ecv_dup$ec_value, interaction(ecv_dup$sp1, ecv_dup$cluster), FUN = max)
			best_sp2 = ave(ecv_dup$ec_value, interaction(ecv_dup$sp2, ecv_dup$cluster), FUN = max)
			ecv_dup_f = ecv_dup[ecv_dup$ec_value == best_sp1 & ecv_dup$ec_value == best_sp2 & ecv_dup$ec_value >= 0, , drop = FALSE]
			ecv_dup_f = ecv_dup_f[!duplicated(ecv_dup_f$sp1) & !duplicated(ecv_dup_f$sp2), c("sp1", "sp2", "ec_value"), drop = FALSE]
			ecv_dup_f$is_o2o = 0
			message(sprintf("ICC markers | pairs of non-o2o homologs kept = %i", nrow(ecv_dup_f)))
			ecv_out = rbind(ecv_o2o, ecv_dup_f)
		}
	}
	list(ec_markers = ecv_out, ec_duplicates = ecv_dup)
}
