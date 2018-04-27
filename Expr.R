library(DESeq2)
library(data.table)
library(dplyr)


directory <- "g12c"

sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*counts).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)


c <- read.table(gzfile("g12C/06db23e3-2dc3-4da9-915e-fe3d4a091de5.htseq.counts.gz"))
d <- read.table(gzfile("g12C/09572b33-5517-462d-9a46-ec3d9b11b59b.htseq.counts.gz"))
e <- read.table(gzfile("g12C/12c4510a-fbad-4bff-97b5-846d4da395bd.htseq.counts.gz"))


l <- lapply(paste0("g12c/",list.files("g12c")), function(x) read.table(gzfile(x)))
r <- bind_cols(l)
r <- r[, c(1, seq(2, ncol(r), 2))]
