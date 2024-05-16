diffChIP_to_bed <- function(path_to_file){
  df = read.csv(path_to_file)
  df_bed = df[2:20]
  write.csv(df_bed, file=glue::glue('{path_to_file}.bed'), row.names = FALSE)
}

annotate_bed <- function(in_path, out_dir, filename){
  
  bed <- read.csv(in_path)
  
  peak <- makeGRangesFromDataFrame(bed)
  
  ann <- annotatePeak(peak=peak, TxDb=txdb, tssRegion=c(-1000, 1000))
  
  annot <- data.frame(ann@anno)
  
  entrez <- annot$geneId
  
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                           keys = entrez,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  
  annot %>% 
    left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
    write.table(file=glue::glue("{out_dir}/{filename}"), sep="\t", quote=F, row.names=F)
}



