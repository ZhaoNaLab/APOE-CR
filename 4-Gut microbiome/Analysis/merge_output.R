library(plyr)


input_dir <- "/research/bsi/projects/microbiome/s216158.apoE_mouse_microbiome/processing/kraken_out_cleaned_final/"
output_dir <- "/research/bsi/projects/microbiome/s216158.apoE_mouse_microbiome/processing/kraken_out_cleaned_final/"

table_names <- list.files(input_dir, pattern="*.genus.025.bracken", recursive=TRUE, full.names=TRUE)
tables <- lapply(table_names, function(x) read.table(x, header=TRUE, comment.char="", sep="\t"))
table_prop <- lapply(tables, function(x) x[,c("name", "fraction_total_reads")])
table_ct <- lapply(tables, function(x) x[,c("name", "new_est_reads")])


for(i in 1:length(table_names)){
	names(table_prop[[i]]) <- c("Genus", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " fraction_total_reads"))
	names(table_ct[[i]]) <- c("Genus", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " new_est_reads"))
}
merged_prop <- plyr::join_all(table_prop, type="full")
merged_prop[is.na(merged_prop)] <- 0
merged_ct <- plyr::join_all(table_ct, type="full")
merged_ct[is.na(merged_ct)] <- 0

write.table(merged_prop, file=paste0(output_dir, "allsamples.genus.prop.tab"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged_ct, file=paste0(output_dir, "allsamples.genus.ct.tab"), sep="\t", quote=FALSE, row.names=FALSE)

table_names <- list.files(input_dir, pattern="*.species.025.bracken", recursive=TRUE, full.names=TRUE)
tables <- lapply(table_names, function(x) read.table(x, header=TRUE, comment.char="", sep="\t"))
table_prop <- lapply(tables, function(x) x[,c("name", "fraction_total_reads")])
table_ct <- lapply(tables, function(x) x[,c("name", "new_est_reads")])

for(i in 1:length(table_names)){
  names(table_prop[[i]]) <- c("Species", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " fraction_total_reads"))
  names(table_ct[[i]]) <- c("Species", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " new_est_reads"))
}
merged_prop <- plyr::join_all(table_prop, type="full")
merged_prop[is.na(merged_prop)] <- 0
merged_ct <- plyr::join_all(table_ct, type="full")
merged_ct[is.na(merged_ct)] <- 0


write.table(merged_prop, file=paste0(output_dir, "allsamples.species.prop.tab"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged_ct, file=paste0(output_dir, "allsamples.species.ct.tab"), sep="\t", quote=FALSE, row.names=FALSE)



table_names <- list.files(input_dir, pattern="*.phylum.025.bracken", recursive=TRUE, full.names=TRUE)
tables <- lapply(table_names, function(x) read.table(x, header=TRUE, comment.char="", sep="\t"))
table_prop <- lapply(tables, function(x) x[,c("name", "fraction_total_reads")])
table_ct <- lapply(tables, function(x) x[,c("name", "new_est_reads")])

for(i in 1:length(table_names)){
  names(table_prop[[i]]) <- c("Phylum", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " fraction_total_reads"))
  names(table_ct[[i]]) <- c("Phylum", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " new_est_reads"))
}
merged_prop <- plyr::join_all(table_prop, type="full")
merged_prop[is.na(merged_prop)] <- 0
merged_ct <- plyr::join_all(table_ct, type="full")
merged_ct[is.na(merged_ct)] <- 0


write.table(merged_prop, file=paste0(output_dir, "allsamples.phylum.prop.tab"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged_ct, file=paste0(output_dir, "allsamples.phylum.ct.tab"), sep="\t", quote=FALSE, row.names=FALSE)




table_names <- list.files(input_dir, pattern="*.family.025.bracken", recursive=TRUE, full.names=TRUE)
tables <- lapply(table_names, function(x) read.table(x, header=TRUE, comment.char="", sep="\t"))
table_prop <- lapply(tables, function(x) x[,c("name", "fraction_total_reads")])
table_ct <- lapply(tables, function(x) x[,c("name", "new_est_reads")])

for(i in 1:length(table_names)){
  names(table_prop[[i]]) <- c("Family", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " fraction_total_reads"))
  names(table_ct[[i]]) <- c("Family", paste0(unlist(strsplit(table_names[[i]], "/"))[10], " new_est_reads"))
}
merged_prop <- plyr::join_all(table_prop, type="full")
merged_prop[is.na(merged_prop)] <- 0
merged_ct <- plyr::join_all(table_ct, type="full")
merged_ct[is.na(merged_ct)] <- 0


write.table(merged_prop, file=paste0(output_dir, "allsamples.family.prop.tab"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(merged_ct, file=paste0(output_dir, "allsamples.family.ct.tab"), sep="\t", quote=FALSE, row.names=FALSE)
