library(dplyr)
# hmm
hmm_search <- function(
    mode = "raw", threshold = 1e-30,
    pfam = NULL, query_fa = NULL,
    prefix = NULL, outdir= getwd()){
  library(dplyr)
  domtblout_file <- paste0(prefix, ".domtblout")
  hmmout_file <- paste0(prefix, ".hmmout")

  args <- c("--cut_tc","--domtblout", domtblout_file,
            "-o", hmmout_file, pfam, query_fa)

  system2("hmmsearch", args = args,stdout = T)

  hmm <- read.table(file = domtblout_file, skip = "#") %>%
    dplyr::filter(V7 < threshold)

  seq <- Biostrings::readAAStringSet(query_fa)
  names(seq) <- sub("\\s.*", "", names(seq))
  hmm_seq <- seq[hmm$V1]

  if(mode == "raw"){

    return(hmm_seq)
  }

  if(mode == "trim"){

    in_maff_file <- file.path(outdir, paste0(prefix, ".fa"))
    out_maff_file <- file.path(outdir, paste0(prefix, ".aln"))
    new_hmm <- file.path(outdir, paste0("new",basename(hmmout_file)))

    two_hmm_seq <- paste0()
    Biostrings::writeXStringSet(x = hmm_seq, filepath = in_maff_file)

    system2("mafft", args = c("--auto", in_maff_file, ">", out_maff_file))
    system2("hmmbuild",args = c(new_hmm, out_maff_file))

  }
}

