# +
library("metacell")
library("devtools")
library("Matrix")
library("data.table")
library("tidyverse")

setwd("/net/mraid14/export/tgdata/users/hernan/wd_NEW/SynEmbs3/")
scdb_init("/net/mraid14/export/tgdata/users/hernan/wd_NEW/scrna_db/", force_reinit=T)
fig_dir <- "/net/mraid14/export/tgdata/users/hernan/wd_NEW/figs/SynEmbs3/"
if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

scfigs_init(fig_dir)
mat_id <- "SynEmb3"
mat = scdb_mat(mat_id)
# -

bad_genes = read.table(file = "/net/mraid14/export/tgdata/users/hernan/wd_NEW/SynEmbs3/SynEmbs3_bad_genes.txt",
                       header = T, stringsAsFactors = F)$x

gset = scdb_gset(mat_id)
gset_f = gset_new_restrict_nms(gset=gset, bad_genes, inverse=T, "feat filt") #excluding from constructing the mc. object

scdb_add_gset(id = mat_id, gset = gset_f)

mcell_add_cgraph_from_mat_bknn(mat_id=mat_id,
                               gset_id = mat_id, 
                               graph_id=mat_id, 
                               K=100, 
                               dsamp=T)

mcell_coclust_from_graph_resamp(coc_id=mat_id, 
                                graph_id=mat_id, 
                                min_mc_size=40, 
                                p_resamp=0.75, 
                                n_resamp=500)

mcell_mc_from_coclust_balanced(coc_id=mat_id,
                               mat_id= mat_id, 
                               mc_id= mat_id,
                               K=50, 
                               min_mc_size=30, 
                               alpha=2)

mcell_mc_split_filt(new_mc_id=mat_id, mc_id=mat_id, mat_id=mat_id, T_lfc=3, plot_mats=F)
mcell_gset_from_mc_markers(gset_id = mat_id, mc_id=mat_id)
mcell_mc_plot_marks(mc_id=mat_id, gset_id=mat_id, mat_id=mat_id)

wt_atlas = mcell_gen_atlas(mat_id = "embexe",
                           mc_id = "embexe_recolored",
                           gset_id  = "embexe",
                           mc2d_id= "embexe_recolored_umap")

mcell_proj_on_atlas(mat_id = mat_id, 
                    mc_id = mat_id, 
                    atlas = wt_atlas, 
                    recolor_mc_id = paste(mat_id,"_recolored",sep=""),
                    fig_cmp_dir = fig_dir,max_entropy = 4, burn_cor = 0.6)

tgconfig::override_params(config_file = "../config/sing_emb.yaml",package = "metacell")

mcell_mc2d_force_knn(mc2d_id = mat_id, symmetrize = F,
                     mc_id = paste(mat_id,"_recolored",sep=""),
                     graph_id = mat_id)

mcell_mc2d_plot(mc2d_id=mat_id,plot_edges = T)
