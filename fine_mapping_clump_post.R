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

  # Remove the inner and outer "[[" and "]]"
  nested_clumps_trim <- str_remove_all(nested_clumps, "\\[\\[|\\]\\]")

  # Create a list of the collected channels, so each item in the list is a channel before they were collected here:  CLUMP_CLUSTER_CH = CLUMP_CH.collect(flat : false)
  # each item in the list corresponds to a different ancestry
  nested_clumps_list <- as.list(unlist(str_split(nested_clumps_trim, "\\], \\[")))

  # clean up white space and spit each item in the list by "," into a vector
  nested_clumps_clean <- lapply(nested_clumps_list, function(l){
                        str_split(l, ",") %>%
                        unlist() %>%
                        str_trim() %>%
                        str_remove_all(., "\\'")
                      })

  # Get details on CHR (7th item in channel "CLUMP_CH" - this code will break if the output order of CLUMP_CH changes)
  chr <- unique ( sapply(nested_clumps_clean, function(l){ l[7] }) )

  # start lapply here over nested_clumps_clean, get a list of GRanges for each ancestry
  # return NULL if there is no GRanges object
  granges_list <- lapply(nested_clumps_clean, function(nested_clump_clean){

    # Read in clumps results
    clumped_data <- fread(nested_clump_clean[5]) # index 5 is the .clumps file (see order in output from CLUMP process)

    if(nrow(clumped_data) == 0){
      # skip this and go to next element in lapply if there is no clumped data
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

  # Reduce list of granges into one granges object, then reduce any overlapping regions (from different ancestries)
  # First deal with all NULLs (ie. no clumps data, and no regions to fine map for all ancestries)
  if (all(sapply(granges_list, is.null))) {
    # If the list is all NULL, create a file with NULL content
    writeLines("NULL", paste0(chr, ".finemapRegions"))
  } else {

  # If the first item in granges_list is a GRanges object then the following code runs:
  # do.call(c, granges_list) %>%
  #                 reduce() 
  # However, if the first item in the list is NULL then that code will result in an error
  # granges_list_null_first <- list(granges_list[[2]], granges_list[[1]], granges_list[[3]])
  #   do.call(c, granges_list_null_first) %>%
  #                   reduce()
  # Error in (function (classes, fdef, mtable)  : 
  # unable to find an inherited method for function ‘reduce’ for signature ‘"list"’

  # Therefore remove NULL items from list to avoid this error
  granges_list <- granges_list[granges_list != "NULL"]

  # Now there are no NULL items, and at least one GRanges object (because of previous "if (all(sapply(granges_list, is.null)))") we
  grng_streched <- do.call(c, granges_list) %>%
                    reduce()

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
  
  