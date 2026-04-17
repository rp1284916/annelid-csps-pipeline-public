# Score ortholog pairs by conserved coexpression context before any
# cross-species cell-type comparisons are made.
source("../R/csps_minimal.R")
graphics.off()

sps_list = get_env_vector("CSPS_SPECIES", c("Ctel", "Ofus"))
sps_refl = get_env_vector("CSPS_REF_SPECIES", c("Ctel"))
processed_dir = get_env_string("CSPS_PROCESSED_DIR", "../work/processed")
results_root = get_env_string("CSPS_RESULTS_DIR", "../work/results")
orthology_fn = get_env_string("CSPS_ORTHOLOGY_FILE", file.path(processed_dir, "orthology", "ortholog_pairs.tsv"))
nthreads_icc = get_env_integer("CSPS_THREADS", 8)
gene_total_umi_threshold = get_env_numeric("CSPS_GENE_TOTAL_UMI_THRESHOLD", 20)
icc_do_duplicates = get_env_flag("CSPS_ICC_DO_DUPLICATES", TRUE)
icc_use_variable_o2o_genes = get_env_flag("CSPS_ICC_USE_VARIABLE_O2O_GENES", TRUE)
icc_variable_o2o_thr = get_env_numeric("CSPS_ICC_VARIABLE_O2O_THR", 0.75)
icc_niter = get_env_integer("CSPS_ICC_NITER", 100)
icc_convergence_thr = get_env_numeric("CSPS_ICC_CONVERGENCE_THR", 0.05)

for (spi in sps_refl) {
	out_fn = file.path(results_root, spi, "csps")
	dir.create(out_fn, recursive = TRUE, showWarnings = FALSE)
	for (spj in sps_list[sps_list != spi]) {
		message(sprintf("csps | %s | compare to %s, load footprints...", spi, spj))
		mc_fp_i = read_matrix_tsv(file.path(processed_dir, spi, sprintf("dat.%s.expression.mcs_fp.tsv.gz", spi)))
		mc_fp_j = read_matrix_tsv(file.path(processed_dir, spj, sprintf("dat.%s.expression.mcs_fp.tsv.gz", spj)))
		message(sprintf("csps | %s | compare to %s, load gene totals...", spi, spj))
		gene_totals_i = read_tsv(file.path(processed_dir, spi, sprintf("dat.%s.gene_total_umis.tsv", spi)))
		gene_totals_j = read_tsv(file.path(processed_dir, spj, sprintf("dat.%s.gene_total_umis.tsv", spj)))
		names_i = dic_from_vecs(gene_totals_i$gene, gene_totals_i$total_umis)
		names_j = dic_from_vecs(gene_totals_j$gene, gene_totals_j$total_umis)
		message(sprintf("csps | %s | compare to %s, load ortholog pairs...", spi, spj))
		og_pairs = read_tsv(orthology_fn)[, c(spi, spj)]
		colnames(og_pairs) = c("sp1", "sp2")
		og_pairs = unique(og_pairs[complete.cases(og_pairs), , drop = FALSE])
		message(sprintf("csps | %s | compare to %s, drop genes with too few counts (>= %.1f total UMIs)...", spi, spj, gene_total_umi_threshold))
		# Low-count genes contribute unstable coexpression estimates and are filtered
		# before the ICC step.
		good_genes_i = intersect(names(names_i)[names_i >= gene_total_umi_threshold], rownames(mc_fp_i))
		good_genes_j = intersect(names(names_j)[names_j >= gene_total_umi_threshold], rownames(mc_fp_j))
		message(sprintf("csps | %s | compare to %s, run ICC (duplicates=%s, variable_o2o=%s, variable_o2o_thr=%.2f)...", spi, spj, icc_do_duplicates, icc_use_variable_o2o_genes, icc_variable_o2o_thr))
		icc = csps_markers_icc_noobj(
			mat_sp1 = mc_fp_i[good_genes_i, , drop = FALSE],
			mat_sp2 = mc_fp_j[good_genes_j, , drop = FALSE],
			og_pairs = og_pairs,
			# Keep duplicate rescue enabled so lineage-specific expansions can still
			# retain their best-supported conserved partner.
			niter = icc_niter,
			icc_thr = icc_convergence_thr,
			do_duplicates = icc_do_duplicates,
			use_variable_o2o_genes = icc_use_variable_o2o_genes,
			variable_o2o_thr = icc_variable_o2o_thr,
			nthreads_icc = nthreads_icc
		)
		write_tsv(icc$ec_markers, file.path(out_fn, sprintf("dat.icc.%s-%s.ec_scores.tsv", spi, spj)))
		if (!is.null(icc$ec_duplicates)) {
			write_tsv(icc$ec_duplicates, file.path(out_fn, sprintf("dat.icc.%s-%s.ec_scores.duplicates.tsv", spi, spj)))
		}
		saveRDS(icc, file.path(out_fn, sprintf("dat.icc.%s-%s.obj.rds", spi, spj)))
	}
}

message("All done!")
