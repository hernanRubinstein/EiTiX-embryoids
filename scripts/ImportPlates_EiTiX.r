library("metacell")
library("devtools")
library("Matrix")
library("data.table")
library("tidyverse")

# +
#load_all("/net/mraid14/export/data/users/atanay/proj/metac/metacell/")
new_dir <- "/net/mraid14/export/tgdata/users/hernan/wd_NEW/SynEmbs3/"
if(!dir.exists(new_dir)) {
    dir.create(new_dir)
  }

setwd(new_dir)
scdb_init("/net/mraid14/export/tgdata/users/hernan/wd_NEW/scrna_db/", force_reinit=T)
fig_dir <- "/net/mraid14/export/tgdata/users/hernan/wd_NEW/SynEmbs3/figs/"
if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

scfigs_init(fig_dir)
# -

plts_dir <- "../plates.csv"
mat_id <- "SynEmb3_new"
metadata = read.csv(file = plts_dir, stringsAsFactors = F)

list_of_experiments <- c("synthetic embryos_MZG")

metadata = metadata[metadata$Experiment %in% list_of_experiments,]

unique(metadata$Plate)
length(unique(metadata$Plate))

excluded_plates = c("MZG_110322_day8_p5", 
                    "MZG_110322_day8_p6", 
                    "MZG_110322_day8_p7", 
                    "MZG_110322_day8_p8", 
                    "MZG_110322_day8_p9", 
                    "MZG_110322_day8_p10", 
                    "MZG_110322_day8_p11", 
                    "MZG_110322_day8_p12", 
                    "MZG_110322_day8_p5", 
                    "MZG_110322_day8_p6", 
                    "MZG_110322_day8_p7", 
                    "MZG_110322_day8_p8", 
                    "MZG_110322_day8_p9", 
                    "MZG_110322_day8_p10", 
                    "MZG_110322_day8_p11", 
                    "MZG_110322_day8_p12")

metadata = metadata[!(metadata$Plate %in% excluded_plates),]

# make competable with name errors ##
metadata$Amp.Batch.ID = metadata$Plate
metadata$Seq.Batch.ID = metadata$Sequencing.Dates
metadata$Batch.Set.ID = metadata$Amp.Batch.ID

write.table(x = metadata,
            file = paste("key_",mat_id,".txt",sep = ""),
            quote = F,
            sep = "\t",
            row.names = F)

##### CREATING THE MAT. OBJECT #####
## multiplexing ##
mcell_import_multi_mars(mat_nm = mat_id,
                        dataset_table_fn = paste("key_",mat_id,".txt",sep = ""),
                        base_dir = "/net/mraid14/export/tgdata/db/tgdb/mars_runs/stelzer_star/plates/umis/",
                        patch_cell_name=T,
                        force=TRUE)

# +
# Filter genes and cells ##
mat = scdb_mat(mat_id)
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = c(grep("^mt\\-", nms, v=T),grep("ERCC", nms,v=T), 
              "Neat1","Atpase6", "Xist", "Malat1", "Cytb",
              "AK018753","AK140265","AK163440","DQ539915")

mcell_mat_ignore_genes(mat_id, mat_id, 
                       bad_genes,reverse=F) 

mcell_mat_ignore_small_cells(mat_id, mat_id, 800)
# -

## Integrating metadata with FCS info ##
mat = scdb_mat(mat_id)

files_with_metadata = list.files(path = "/net/mraid14/export/tgdata/db/tgdb/mars_runs/stelzer_star/plates/metadata/")

files_with_metadata = intersect(files_with_metadata,paste(metadata$Plate,".tsv",sep=""))

extra_cell_metadata = data.frame()

for (plate_id in files_with_metadata) {
  print(plate_id)
  extra_cell_metadata_temp = read.table(paste("/net/mraid14/export/tgdata/db/tgdb/mars_runs/stelzer_star/plates/metadata/",plate_id,sep = ""),
                                        sep = "\t",h = T,stringsAsFactors = F)
  colnames(extra_cell_metadata_temp)[colnames(extra_cell_metadata_temp) == "dtomato_a"] = "tdtomato_a"
  #extra_cell_metadata = rbind(extra_cell_metadata,extra_cell_metadata_temp,fill = TRUE)
  extra_cell_metadata = bind_rows(extra_cell_metadata,extra_cell_metadata_temp)
}

# +
##### Remove the duplicate cells ####
duplicate_cells = names(table(extra_cell_metadata$cell)[table(extra_cell_metadata$cell) >1])
tot_number_of_duplicates = sum(table(extra_cell_metadata$cell)[table(extra_cell_metadata$cell) >1]) - length(duplicate_cells)

for (cell_id in duplicate_cells) {
  duplicates = which(extra_cell_metadata$cell == cell_id)
  # remove the first element from the duplicate list
  duplicates = duplicates[-1]
  extra_cell_metadata = extra_cell_metadata[-duplicates,]
}
extra_cell_metadata$cell = paste(extra_cell_metadata$plate,extra_cell_metadata$cell,sep=".")
# -

temp = mat@cell_metadata
temp$cell = rownames(temp)
temp = left_join(temp,extra_cell_metadata,by="cell")
temp$embryo[temp$embryo == ""] = "empty"
rownames(temp) = temp$cell
mat@cell_metadata = temp

new_ignore_cells = union(mat@ignore_cells,mat@cell_metadata$cell[mat@cell_metadata$embryo == "empty"])
new_ignore_cells = union(new_ignore_cells,duplicate_cells)
mat = scm_ignore_cells(scmat = mat,ig_cells = new_ignore_cells)

scdb_add_mat(id = mat_id,mat = mat)


scdb_ls("mat", regexp = mat_id)
