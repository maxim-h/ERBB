library(GenomicDataCommons)
library(magrittr)

ge_manifest = files() %>% 
  filter( ~ cases.project.project_id == 'TCGA-OV' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()


files()$cases.case_id
?manifest
cases()$fields

g12c <- read.table("g12c.tsv", header = T)


ge_manifest = files() %>% 
  filter( ~ cases.case_id == "TCGA-05-441") %>%
  manifest()


ge_manifest = files() %>% 
  filter( ~ cases.project.project_id == 'TCGA-LUAD' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts' &
            cases.case_id %in% g12c$Case)

files() %>% available_fields()



