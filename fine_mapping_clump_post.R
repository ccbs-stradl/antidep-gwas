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

  nested_clumps <- '${nested_clumps}'

  nested_clumps_trim <- str_remove_all(nested_clumps, "\\[\\[|\\]\\]")

  nested_clumps_list <- as.list(unlist(str_split(nested_clumps_trim, "\\], \\[")))

  nested_clumps_clean <- lapply(nested_clumps_list, function(l){
                        str_split(l, ",") %>%
                        unlist() %>%
                        str_trim() %>%
                        str_remove_all(., "\\'")
                      })

  # Get details on CHR
  chr <- unique ( sapply(nested_clumps_clean, function(l){ l[[7]] }) )

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

  if (all(sapply(granges_list, is.null))) {
    # If the list is all NULL, create a file with NULL content
    writeLines("NULL", paste0(, ".finemapRegions"))
  } else {

  # Reduce list of granges into one granges object, then reduce any overlapping regions (from different ancestries)
  grng_streched <- do.call(c, granges_list) %>%
                    reduce()


  # Reduce ranges to collapse overlapping or nearby regions
  grng_reduced <- grng_streched %>%
    reduce_ranges(min.gapwidth = 5000)  %>% # What should this be set to?
    as_tibble() %>%
    mutate(WIDTH = end - start + 1) %>%
    dplyr::select(CHR = seqnames, BP_START = start, BP_END = end, WIDTH)
  

  # Save a table of min and max BP positions for SuSiEx, per chr
  write.table(grng_reduced, paste0(chr, ".finemapRegions"), row.names = F, quote = F, sep = "\t")
  }

  # capture the chr
  cat(chr)
  