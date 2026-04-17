# Build the weighted cross-species cell-type similarity matrices used for the
# main heatmap outputs.
source("../R/csps_minimal.R")
graphics.off()

sps_list = get_env_vector("CSPS_SPECIES", c("Ctel", "Ofus"))
sps_refl = get_env_vector("CSPS_REF_SPECIES", c("Ctel"))
processed_dir = get_env_string("CSPS_PROCESSED_DIR", "../work/processed")
results_root = get_env_string("CSPS_RESULTS_DIR", "../work/results")
cross_fc_thr = get_env_numeric("CSPS_FC_THR", 2)

for (spi in sps_refl) {
	out_fn = file.path(results_root, spi, "csps")
	dir.create(out_fn, recursive = TRUE, showWarnings = FALSE)
	for (spj in sps_list[sps_list != spi]) {
		message(sprintf("csps | %s | compare to %s, load ICC...", spi, spj))
		icc_ecv = read_tsv(file.path(out_fn, sprintf("dat.icc.%s-%s.ec_scores.tsv", spi, spj)))
		icc_ecv_v = dic_from_vecs(icc_ecv$sp1, icc_ecv$ec_value)
		message(sprintf("csps | %s | compare to %s, load footprints...", spi, spj))
		mc_fp_i = read_matrix_tsv(file.path(processed_dir, spi, sprintf("dat.%s.expression.cts_fp.tsv.gz", spi)))
		mc_fp_j = read_matrix_tsv(file.path(processed_dir, spj, sprintf("dat.%s.expression.cts_fp.tsv.gz", spj)))
		message(sprintf("csps | %s | compare to %s, load cell type annotations...", spi, spj))
		ctt_i = read_tsv(file.path(processed_dir, spi, sprintf("annot.%s.leiden.tsv", spi)))
		ctt_j = read_tsv(file.path(processed_dir, spj, sprintf("annot.%s.leiden.tsv", spj)))
		ctt_i$cell_type_species = sprintf("%s|%s", spi, ctt_i$cell_type)
		ctt_j$cell_type_species = sprintf("%s|%s", spj, ctt_j$cell_type)
		ctt_i_cts_col_v = dic_from_vecs(ctt_i$cell_type_species, ctt_i$color)
		ctt_j_cts_col_v = dic_from_vecs(ctt_j$cell_type_species, ctt_j$color)
		ctt_i_cts_col_vv = dic_from_vecs(ctt_i$cell_type, ctt_i$color)
		ctt_j_cts_col_vv = dic_from_vecs(ctt_j$cell_type, ctt_j$color)
		message(sprintf("csps | %s | compare to %s, conform matrices...", spi, spj))
		mc_fp_i_f = mc_fp_i[icc_ecv$sp1, , drop = FALSE]
		mc_fp_j_f = mc_fp_j[icc_ecv$sp2, , drop = FALSE]
		mc_fp_m_f = cbind(mc_fp_i_f, mc_fp_j_f)
		message(sprintf("csps | %s | compare to %s, get covariable genes...", spi, spj))
		# Restrict the comparison to genes that show strong enrichment in both
		# species, then weight them by their ICC conservation score.
		genes_covariable = csps_select_covariable_genes(
			sp1_fp = mc_fp_i_f,
			sp2_fp = mc_fp_j_f,
			merged = mc_fp_m_f,
			cross_fc_thrs = cross_fc_thr,
			cross_n = 1,
			method = "min_fc"
		)
		message(sprintf("csps | %s | compare to %s, wpearson...", spi, spj))
		icc_out = csps_correlation_matrix_noobj(
			mm = mc_fp_m_f,
			m1 = mc_fp_i_f,
			m2 = mc_fp_j_f,
			prefix_sp1 = spi,
			prefix_sp2 = spj,
			use_var_genes = genes_covariable[[1]],
			gene_weights = icc_ecv_v[genes_covariable[[1]]],
			cor_method = "wpearson"
		)
		pp1 = plot_complex_heatmap(
			icc_out$cor_matrix,
			cluster_row = FALSE,
			cluster_col = FALSE,
			color_min = 0,
			color_max = 0.7,
			fontsize = 8,
			categories_row = rownames(icc_out$cor_matrix),
			categories_col = colnames(icc_out$cor_matrix),
			colors_row = ctt_i_cts_col_v,
			colors_col = ctt_j_cts_col_v,
			color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
			cell_border = grid::gpar(col = "white", lwd = 1, lty = 1),
			heatmap_border = grid::gpar(col = "black", lwd = 1, lty = 1)
		)
		message(sprintf("csps | %s | compare to %s, save...", spi, spj))
		saveRDS(icc_out, file.path(out_fn, sprintf("csps.%s-%s.cts.rds", spi, spj)))
		grDevices::pdf(file.path(out_fn, sprintf("csps.%s-%s.cts.pdf", spi, spj)), width = 8 + ncol(icc_out$cor_matrix) / 20, height = 8 + nrow(icc_out$cor_matrix) / 20)
		print(pp1)
		grDevices::dev.off()
		dup_path = file.path(out_fn, sprintf("dat.icc.%s-%s.ec_scores.duplicates.tsv", spi, spj))
		if (file.exists(dup_path)) {
			message(sprintf("csps | %s | compare to %s, EC duplicate plots...", spi, spj))
			icc_dup_ecv = read_tsv(dup_path)
			cc_xt = table(icc_dup_ecv$cluster)
			# Only plot small ambiguous ortholog groups; large duplicate blocks become
			# unreadable and are better inspected from the TSV output.
			cc_xt_f = names(cc_xt)[cc_xt < 4]
			icc_dup_ecv_f = icc_dup_ecv[icc_dup_ecv$cluster %in% cc_xt_f, , drop = FALSE]
			grDevices::pdf(file.path(out_fn, sprintf("csps.%s-%s.cts.ec_scores.duplicates.pdf", spi, spj)), width = 8, height = 20)
			for (clu in cc_xt_f) {
				icc_dup_ecv_ff = icc_dup_ecv_f[icc_dup_ecv_f$cluster == clu, , drop = FALSE]
				icc_dup_ecv_ff$ec_value[is.na(icc_dup_ecv_ff$ec_value)] = 0
				icc_dup_ecv_ff = icc_dup_ecv_ff[order(icc_dup_ecv_ff$ec_value, decreasing = TRUE), , drop = FALSE]
				layout(matrix(1:20, ncol = 2, byrow = TRUE))
				for (nn in seq_len(nrow(icc_dup_ecv_ff))) {
					ggi = icc_dup_ecv_ff[nn, "sp1"]
					ggj = icc_dup_ecv_ff[nn, "sp2"]
					ecv = icc_dup_ecv_ff[nn, "ec_value"]
					barplot(mc_fp_i[ggi, ], las = 2, ylab = "Footprint", col = ctt_i_cts_col_vv[colnames(mc_fp_i)], main = sprintf("%s\nEC=%.2f", ggi, ecv), cex.names = 0.5, cex.main = 0.8)
					abline(h = 1, lty = 2)
					barplot(mc_fp_j[ggj, ], las = 2, ylab = "Footprint", col = ctt_j_cts_col_vv[colnames(mc_fp_j)], main = ggj, cex.names = 0.5, cex.main = 0.8)
					abline(h = 1, lty = 2)
				}
			}
			grDevices::dev.off()
		}
	}
}

message("All done!")
