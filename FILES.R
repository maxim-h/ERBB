library(GenomicDataCommons)
library(magrittr)

download_reads <- function(f){
  l <- read.table(paste0(f, ".tsv"), header = T, sep="\t")

  ge_manifest = files() %>% 
    filter( ~ cases.project.project_id == 'TCGA-LUAD' &
              type == 'gene_expression' &
              analysis.workflow_type == 'HTSeq - Counts' &
              cases.case_id %in% rownames(l)) %>%
    manifest()
  
  lapply(ge_manifest$id,gdcdata, destination_dir=f, overwrite=FALSE, progress=TRUE)
}

lapply(sub(".tsv", "", list.files(pattern = "*.tsv")), download_reads)



