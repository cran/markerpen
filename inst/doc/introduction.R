## ----data---------------------------------------------------------------------
library(markerpen)
library(scales)

# A subset of the ROSMAP data
# Subsampled from the geneExprRaw.txt file located in
# https://github.com/ellispatrick/CortexCellDeconv/tree/master/CellTypeDeconvAnalysis/Data
load(system.file("examples", "gene_expr.RData", package = "markerpen"))

# Convert gene name to Ensembl name
ind = match(rownames(dat), markerpen::gene_mapping$name)
ind = na.omit(ind)
ensembl = markerpen::gene_mapping$ensembl[ind]

# Get expression matrix - rows are observations, columns are genes
mat_exp = t(dat[markerpen::gene_mapping$name[ind], ])
colnames(mat_exp) = ensembl
print(mat_exp[1:5, 1:5])

## ----marker_list--------------------------------------------------------------
# Read in prior marker genes
load(system.file("examples", "published_markers.RData", package = "markerpen"))
load(system.file("examples", "markers_range.RData", package = "markerpen"))

## ----refine-------------------------------------------------------------------
# Markers for astrocytes
ast_re = refine_markers(mat_exp, markers_range$astrocytes, pub_markers$astrocytes,
                        lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
# Remove selected markers from the expression matrix
mat_rest = mat_exp[, setdiff(colnames(mat_exp), ast_re$markers)]

# Markers for oligodendrocytes
oli_re = refine_markers(mat_rest, markers_range$oligodendrocytes, pub_markers$oligodendrocytes,
                        lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
mat_rest = mat_rest[, setdiff(colnames(mat_rest), oli_re$markers)]

# Markers for microglia
mic_re = refine_markers(mat_rest, markers_range$microglia, pub_markers$microglia,
                        lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
mat_rest = mat_rest[, setdiff(colnames(mat_rest), mic_re$markers)]

# Markers for endothelial
end_re = refine_markers(mat_rest, markers_range$endothelial, pub_markers$endothelial,
                        lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)
mat_rest = mat_rest[, setdiff(colnames(mat_rest), end_re$markers)]

# Markers for neurons
neu_re = refine_markers(mat_rest, markers_range$neurons, pub_markers$neurons,
                        lambda = 0.45, w = 1.5, maxit = 500, eps = 1e-3, verbose = 0)

# Refined markers
markers_re = list(astrocytes       = ast_re$markers,
                  oligodendrocytes = oli_re$markers,
                  microglia        = mic_re$markers,
                  endothelial      = end_re$markers,
                  neurons          = neu_re$markers)

## ----order--------------------------------------------------------------------
# Post-process selected markers
# Pick the first 50 ordered markers
cor_markers = cor(mat_exp[, unlist(markers_re)])
markers_ord = sort_markers(cor_markers, markers_re)
markers_ord = lapply(markers_ord, head, n = 50)

## -----------------------------------------------------------------------------
# Function to visualize the sample correlation matrix
vis_cor = function(mat_exp, markers)
{
    all_genes = colnames(mat_exp)
    markers = intersect(unlist(markers), all_genes)
    cor_markers = cor(mat_exp[, unlist(markers)])
    p = nrow(cor_markers)

    cols = c("#08306b", "#08519c", "#2171b5", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7",
             "#ffffff",
             "#fcf1f1", "#fae1e1", "#facdcd", "#f49c9c", "#f56566", "#f13a3c", "#d00003")

    ncols = length(cols)
    cols = scales::gradient_n_pal(cols, values = (0:ncols) / ncols)((1:100) / 100)
    op = par(mar = c(0, 0, 0, 0))
    image(cor_markers[, p:1], col = cols, breaks = (-50:50) / 50, asp = 1, axes = FALSE)
    par(op)
}

## ----fig.align='center'-------------------------------------------------------
vis_cor(mat_exp, pub_markers)

## ----fig.align='center'-------------------------------------------------------
vis_cor(mat_exp, markers_re)

## ----fig.align='center'-------------------------------------------------------
vis_cor(mat_exp, markers_ord)

