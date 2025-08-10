extract_seq_from_pfam <- function(pfamID = NULL, save = TRUE, outdir = "./"){
  data("tair_pep")
  data("tair_pfma")
  data("rice_pep")
  data("rice_pfam")

  pattern <- paste0("^", pfamID)

  tair_blast_id <- tair_pfma$V1[grepl(pattern = pattern, x = tair_pfma$V4)]
  rice_blast_id <- rice_pfam$model[grepl(pattern = pattern, x = rice_pfam$hmm_acc)]

  tair_blast_seq <- tair_pep[tair_blast_id]
  rice_blast_seq <- rice_pep[rice_blast_id]

  merge_seq <- c(tair_blast_seq, rice_blast_seq)

  if(save){

    Biostrings::writeXStringSet(merge_seq, filepath = paste0(outdir, "/", "merge_blast_seq.fa"))

  }
  return(seq_list)
}

test <- extract_seq_from_pfam(pfamID = "PF02365", outdir = "../genefam/test/")



make_blastdb <- function(in_seq = NULL,
                         dbtype = "prot", db_name = NULL,
                         verbose = TRUE){

  args <- c(
    "-in", in_seq,
    "-dbtype", dbtype,
    ifelse(!is.null(db_name),
           paste("-out", db_name))
    )

  system2(
    "makeblastdb",
    args = args,
    stdout = ifelse(verbose, "", FALSE)
    )
  }

make_blastdb(in_seq = "~/genefam/test/merge_blast_seq.fa", db_name = "merge")

blast <- function(type = "blastp", db_name = NULL,
                  query_seq = NULL, evalue = 1e-10,
                  outfmt = 6, threads = 1, num_alignments = 10, outdir = NULL,
                  other_args = NULL){
  args = c("-query", query_seq, "-db", db_name,
           "-evalue", evalue, "-num_threads",
           threads, "-num_alignments", num_alignments,
           "-out", file.path(outdir, paste0(db_name, ".outfmt", outfmt)),
           ifelse(!is.null(other_args), other_args)
           )
  if(type == blastp){
    system2(command = "blastp", args = args)
  }
  if(type == blastn){
    system2(command = "blastn", args = args)
  }

}

extract_blastout <- function()
