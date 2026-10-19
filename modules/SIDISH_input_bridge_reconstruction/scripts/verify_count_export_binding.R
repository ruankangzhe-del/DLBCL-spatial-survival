args<-commandArgs(TRUE);suppressPackageStartupMessages(library(Matrix))
o<-readRDS(args[1]);a<-attr(o,'assays')[['RNA']];layers<-attr(a,'layers');out<-args[3]
capture.output(lapply(layers,dim),file=file.path(out,'RDS_count_layer_dimensions.txt'))
stopifnot('counts'%in%names(layers))
cc<-unclass(attr(a,'cells'));ff<-unclass(attr(a,'features'))
cells<-rownames(cc)[cc[,'counts']];genes<-rownames(ff)[ff[,'counts']]
ec<-read.delim(file.path(args[2],'sc_cells.tsv'),stringsAsFactors=FALSE)[[1]]
eg<-read.delim(file.path(args[2],'sc_genes.tsv'),stringsAsFactors=FALSE)[[1]]
stopifnot(identical(cells,ec),all(eg%in%genes),!anyDuplicated(eg))
saved<-layers[['counts']];m<-readMM(file.path(args[2],'sc_counts_genes_by_cells.mtx'))
sub<-saved[match(eg,genes),,drop=FALSE]
delta<-drop0(sub-m)
omitted<-setdiff(genes,eg);omitted_nnz<-if(length(omitted))nnzero(saved[match(omitted,genes),,drop=FALSE]) else 0
r<-data.frame(check=c('same_cell_order','export_genes','RDS_genes','different_count_entries','max_count_difference','omitted_genes','omitted_gene_nonzeros'),value=c(TRUE,length(eg),length(genes),nnzero(delta),if(nnzero(delta))max(abs(delta@x)) else 0,paste(omitted,collapse=','),omitted_nnz))
write.table(r,file.path(out,'RDS_export_count_identity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(r);stopifnot(nnzero(delta)==0,omitted_nnz==0)
