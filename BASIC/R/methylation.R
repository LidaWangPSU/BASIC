if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k["cg07730886", ]
annEPIC["cg07730886", ]



