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

  domtblout <- read.table(file = domtblout_file, skip = "#") %>%
    dplyr::filter(V7 < threshold)

  query_seq <- Biostrings::readAAStringSet(query_fa)
  names(query_seq) <- sub("\\s.*", "", names(query_seq))
  hmm_seq <- query_seq[domtblout$V1]

  if(mode == "raw"){

    return(hmm_seq)
  }

  if(mode == "trim"){

    in_maff_file <- file.path(outdir, paste0(prefix, ".fa"))
    out_maff_file <- file.path(outdir, paste0(prefix, ".aln"))
    new_hmm <- file.path(outdir, paste("new",basename(pfam), sep = "_"))
    new_domtblout_file <- file.path(outdir, paste("new",basename(domtblout_file), sep = "_"))
    new_hmmout_file <- file.path(outdir, paste("new",basename(hmmout_file), sep = "_"))

    Biostrings::writeXStringSet(x = hmm_seq, filepath = in_maff_file)

    system2("mafft", args = c("--auto", in_maff_file, ">", out_maff_file))
    system2("hmmbuild",args = c(new_hmm, out_maff_file))

    args <- c("--domtblout", new_domtblout_file,
              "-o", new_hmmout_file, new_hmm, query_fa)
    system2("hmmsearch", args = args,stdout = T)

    new_domtblout <- read.table(file = new_domtblout_file, skip = "#") %>%
      dplyr::filter(V7 < threshold)
    new_hmm_seq <- query_seq[new_domtblout$V1]

    return(new_hmm_seq)
  }
}

new_test <- hmm_search(mode = "trim",pfam = "~/genefam/hmm/PF00759.hmm",query_fa = "~/genefam/genome/wm82_a4/Gmax_508_Wm82.a4.v1.protein_primaryTranscriptOnly.fa",prefix = "test",outdir = "~/genefam/test/")
test <- hmm_search(mode = "raw",pfam = "~/genefam/hmm/PF00759.hmm",query_fa = "~/genefam/genome/wm82_a4/Gmax_508_Wm82.a4.v1.protein_primaryTranscriptOnly.fa",prefix = "test",outdir = "~/genefam/test/")

diff_seq_count <- function(a = NULL, b = NULL){
  a_set <- names(a)
  b_set <- names(b)
  a_uniq <- a_set[!a_set %in% b_set]
  b_uniq <- b_set[!b_set %in% a_set]
  diff_summary <- list(a_uniq = a_uniq, b_uniq = b_uniq)
  return(diff_summary)
}

sumy <- diff_seq_count(test, new_test)






