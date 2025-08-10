library(dplyr)
tair_pfma <- read.table("data/blast/tair/all.domain.txt",sep = "\t") %>%
  select(-c(2,3,8,9))
tair_pep <- Biostrings::readAAStringSet("data/blast/tair/Arabidopsis_thaliana.TAIR10.pep.all.fa")

rice_pfam <- read.table("data/blast/rice/all.pfam",header = T)
rice_pep <- Biostrings::readAAStringSet("data/blast/rice/Osativa_323_v7.0.protein.fa")

usethis::use_data(tair_pfma, tair_pep, rice_pfam, rice_pep, compress = "xz", overwrite = T)
