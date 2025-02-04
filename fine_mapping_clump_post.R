  library(stringr)
  library(data.table)
  library(dplyr)

  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if (!require("plyranges", quietly = TRUE)) {
    BiocManager::install("plyranges")
  }

  library(plyranges) # reduce_ranges

  # Get command-line arguments
  args <- commandArgs(trailingOnly = TRUE)

  # The first argument is the nested_clumps string
  nested_clumps <- args[1]

  nested_clumps <- "[[EUR, fixed, N06A, /exports/eddie/scratch/aedmond3/ad/work/ca/f5a554f9730a7e45262a0815630bbf/fixed-N06A-EUR.ma, /exports/eddie/scratch/aedmond3/ad/work/ca/f5a554f9730a7e45262a0815630bbf/fixed-N06A-EUR.3.clumps, /exports/eddie/scratch/aedmond3/ad/work/ca/f5a554f9730a7e45262a0815630bbf/fixed-N06A-EUR.3.log, 3], [AFR, fixed, N06A, /exports/eddie/scratch/aedmond3/ad/work/69/2731963a60ec1c8d7c77acddd5b7b4/fixed-N06A-AFR.ma, /exports/eddie/scratch/aedmond3/ad/work/69/2731963a60ec1c8d7c77acddd5b7b4/fixed-N06A-AFR.3.clumps, /exports/eddie/scratch/aedmond3/ad/work/69/2731963a60ec1c8d7c77acddd5b7b4/fixed-N06A-AFR.3.log, 3], [SAS, fixed, N06A, /exports/eddie/scratch/aedmond3/ad/work/a1/56c1329b83ed8e0edf449a317ca361/fixed-N06A-SAS.ma, /exports/eddie/scratch/aedmond3/ad/work/a1/56c1329b83ed8e0edf449a317ca361/fixed-N06A-SAS.3.clumps, /exports/eddie/scratch/aedmond3/ad/work/a1/56c1329b83ed8e0edf449a317ca361/fixed-N06A-SAS.3.log, 3]]"

  cat(nested_clumps)

  nested_clumps_trim <- str_remove_all(nested_clumps, "\\[\\[|\\]\\]")

  nested_clumps_list <- as.list(unlist(str_split(nested_clumps_trim, "\\], \\[")))

  nested_clumps_clean <- lapply(nested_clumps_list, function(l){
                        str_split(l, ",") %>%
                        unlist() %>%
                        str_trim() %>%
                        str_remove_all(., "\\'")
                      })

  # Get details on CHR
  chr <- unique ( sapply(nested_clumps_clean, function(l){ l[7] }) )

  # start lapply here over nested_clumps_clean
  granges_list <- lapply(nested_clumps_clean, function(nested_clump_clean){

    # Read in clumps results
    clumped_data <- fread(nested_clump_clean[5]) # index 5 is the .clumps file (see order in output from CLUMP process)

    if(nrow(clumped_data) == 0){
      # skip this and go to next element in lapply
      return() # returns NULL
    } else {

      loci <- clumped_data %>%
        dplyr::select(CHR= `#CHROM`, POS, SNP=ID)

      hg19 <- genome_info('hg19')
      # pull out lengths manually since seqnames uses "chrN" instead of "N"
      hg19_chr_lengths <- as_tibble(hg19) |> slice(1:23) |> pull(width)

      # add genome info for autosomes and X
      grng <- loci %>%
                arrange(CHR) %>%
                as_granges(seqnames = CHR,
                              start = POS,
                              end = POS) 

      seqlevels(grng) <- as.character(1:23)
      seqlengths(grng) <- hg19_chr_lengths
      grng <- set_genome_info(grng, 'hg19', is_circular = rep(FALSE, 23))

      # stretch each region, by 100 kb upstream and downstream, then trim back to position boundaries
      grng_stretched <- stretch(anchor_center(grng), 200000) %>% trim()

      return(grng_stretched)
    }

  })

  granges_list

  if (all(sapply(granges_list, is.null))) {
    # If the list is all NULL, create a file with NULL content
    writeLines("NULL", paste0(chr, ".finemapRegions"))
  } else {


# # Test code to get overlapping regions for multiple ancestries
# gr1 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 100, end = 200), 
#   strand = "*", 
#   SNP = "rs111111"
# )
# gr2 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 100, end = 200), 
#   strand = "*", 
#   SNP = "rs111112"
# )
# # Example 2: Non-Overlapping Ranges
# gr3 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 300, end = 400), 
#   strand = "*", 
#   SNP = "rs222222"
# )
# gr4 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 500, end = 600), 
#   strand = "*", 
#   SNP = "rs222223"
# )
# # Example 3: Partly Overlapping Ranges
# gr5 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 150, end = 250), 
#   strand = "*", 
#   SNP = "rs333333"
# )
# gr6 <- GRanges(
#   seqnames = "1", 
#   ranges = IRanges(start = 180, end = 280), 
#   strand = "*", 
#   SNP = "rs333334"
# )
# # Combine the GRanges objects into a list
# granges_list_test <- list(gr1, gr2, gr3, gr4, gr5, gr6, NULL)


  # Reduce list of granges into one granges object, then reduce any overlapping regions (from different ancestries)
  if(inherits(do.call(c, granges_list),  "GRanges")){
    grng_streched <- do.call(c, granges_list) %>%
                    reduce()
  }else{
    grng_streched <- do.call(c, granges_list)[[1]] %>%
                    reduce()
  }

  # If the above throws an error look into this further, as with test data the following works:
    # do.call(c, granges_list_test) %>%
    #                 reduce()


  # Reduce ranges to collapse overlapping or nearby regions
  grng_reduced <- grng_streched %>%
    reduce_ranges(min.gapwidth = 5000)  %>% # What should this be set to?
    as_tibble() %>%
    mutate(WIDTH = end - start + 1) %>%
    dplyr::select(CHR = seqnames, BP_START = start, BP_END = end, WIDTH)
  

  # Save a table of min and max BP positions for SuSiEx, per chr
  write.table(grng_reduced, paste0("chr", chr, ".finemapRegions"), row.names = F, quote = F, sep = "\t")
  }


# Output chr
writeLines(chr, "chr.txt")
  
  