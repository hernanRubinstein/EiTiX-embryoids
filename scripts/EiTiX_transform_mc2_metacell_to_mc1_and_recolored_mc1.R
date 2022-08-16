library("reticulate")
use_python("/net/mraid14/export/data/users/eladch/tools/CO7/python3/3.7.5/bin/python", required=TRUE)
library("anndata")
sc = import("scanpy")

library(metacell)
mc2 = import("metacells")


# +
##### TODO change mat_id and mc2_db_dir to what you need #####
mc2_db_dir = "/net/mraid14/export/tgdata/users/hernan/wd_NEW/scrna_db/"

mat_id = "mat.iETX3.h5ad"
mc_id_2 = "mc.iETX3.h5ad"
######
# -

# 2. Run MC2 using the vignettes, I'm showing here the process where you possibly also generated clusters and 2D projections from MC2, you can ignore the clusters and 2D as needed


# 3. Import MC2 output into a fake mc object and save it (I'm also replacing the egc with the fractions)
mc_2 = anndata::read_h5ad(paste(mc2_db_dir,mc_id_2, sep = ""))
mc_2_cells = anndata::read_h5ad(paste(mc2_db_dir,mat_id, sep = ""))

# + active=""
# ### optional for clusters and colors from mc2
# mc2_cols_tab = fread(paste(mc2_db_dir,"/cluster-colors.csv", sep = ""))
# mc_clusts = as.character(mc_2$obs$cluster)
# mc2_cols = mc2_cols_tab$color
# names(mc2_cols) = mc2_cols_tab$cluster
# mc_colors = mc2_cols[mc_clusts]
# names(mc_colors) = names(mc_clusts)
# mc_colors = mc_colors[order(as.integer(names(mc_colors)))]
# ### end optional
# -

cells_mc = mc_2_cells$obs$metacell
cells_mc = cells_mc + 1
cell_names = mc_2_cells$obs_names
names(cells_mc) = cell_names
gene_names = mc_2$var_names

outlier_cells = names(cells_mc)[cells_mc == 0]
cells_mc_filt = cells_mc[cells_mc > 0]

umis2 = as.matrix(mc_2$X)
fractions = umis2 / rowSums(umis2)
fractions = t(fractions)

scdb_init(mc2_db_dir,force_reinit = T)
mat = scdb_mat("iETX3")
dim(mat@mat)

scdb_add_mc(mc_id_2, tgMCCov(cells_mc_filt, c(outlier_cells, setdiff(mat@cells, c(outlier_cells, names(cells_mc_filt)))), mat))
mc = scdb_mc(mc_id_2)
colnames(fractions) = colnames(mc@e_gc)
mc@e_gc = fractions

scdb_init(mc2_db_dir, force_reinit = T)
mc_ref <- scdb_mc("embexe_recolored")
mc@color_key <- mc_ref@color_key

ct <- mc_2$obs$projected_type
ct_to_col <- mc@color_key$color
names(ct_to_col) <- mc@color_key$group
ct_col <- ct_to_col[match(ct, names(ct_to_col))]
mc@colors <- array(ct_col)

scdb_init(mc2_db_dir,force_reinit = T)
scdb_add_mc("iETX3_recolored", mc)

# You can now use above mc normally in R

# +
setwd("/net/mraid14/export/tgdata/users/hernan/wd_NEW/SynEmbs3/")

# If you need the 2D in R
mc_x = mc2$ut$get_o_numpy(mc_2, 'umap_x')
mc_y = mc2$ut$get_o_numpy(mc_2, 'umap_y')

xrange = 0.03 * (max(mc_x) - min(mc_x))
yrange = 0.03 * (max(mc_y) - min(mc_y))

sc_x = mc_x[mc@mc] + rnorm(length(names(mc@mc)), 0, xrange)
sc_y = mc_y[mc@mc] + rnorm(length(names(mc@mc)), 0, yrange)

xlim = c(min(mc_x), max(mc_x))
ylim = c(min(mc_y), max(mc_y))

sc.cols <- mc@colors[mc@mc]


pdf(file = "figs/SynEms-3.0_umap.pdf",width = 8,height = 8)
plot(sc_x, sc_y, col = sc.cols, 
    pch = 19, cex = 0.8, xaxt = "n", yaxt = "n", xlab = "", 
    ylab = "", axes = T, main = "Synthetic Embryos GD6 + 8 2d-umap", cex.main = 2)
points(mc_x, mc_y, pch = 21, bg = mc@colors, 
    ylim = ylim, xlim = xlim, cex = 2.5, col = "black")
dev.off()

options(repr.plot.width=10, repr.plot.height=10)
plot(sc_x, sc_y, col = sc.cols, 
    pch = 19, cex = 0.8, xaxt = "n", yaxt = "n", xlab = "", 
    ylab = "", axes = T, main = "Synthetic Embryos GD6 + 8 2d-umap", cex.main = 2)
points(mc_x, mc_y, pch = 21, bg = mc@colors, 
    ylim = ylim, xlim = xlim, cex = 2.5, col = "black")
# -

length(unique(mc@mc))
