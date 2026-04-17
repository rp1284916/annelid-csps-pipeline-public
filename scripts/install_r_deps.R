# Install required packages into a user-writable library so the same script
# works on shared systems such as Appocrita.
user_lib = Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(user_lib)) {
	user_lib = file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."))
	Sys.setenv(R_LIBS_USER = user_lib)
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# Keep the dependency list centralised here so local and cluster runs install
# the same package set.
cran_packages = c(
	"data.table",
	"igraph",
	"stringr",
	"circlize",
	"ape",
	"phangorn",
	"umap",
	"scales",
	"adephylo",
	"fastcluster",
	"zoo",
	"R.utils"
)

bioc_packages = c("ComplexHeatmap")
installed = rownames(installed.packages())
missing_cran = setdiff(cran_packages, installed)
missing_bioc = setdiff(bioc_packages, installed)

if (length(missing_cran) == 0 && length(missing_bioc) == 0) {
	message("All required R packages are already installed.")
} else {
	if (length(missing_cran) > 0) {
		message(sprintf("Installing missing CRAN packages: %s", paste(missing_cran, collapse = ", ")))
		install.packages(missing_cran, repos = "https://cloud.r-project.org", lib = user_lib)
	}
	if (length(missing_bioc) > 0) {
		if (!requireNamespace("BiocManager", quietly = TRUE)) {
			install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = user_lib)
		}
		message(sprintf("Installing missing Bioconductor packages: %s", paste(missing_bioc, collapse = ", ")))
		BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, lib = user_lib)
	}
}
