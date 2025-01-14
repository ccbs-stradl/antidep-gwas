# Loop over each fine mapped region
for(i in 1:length(results$cs)){
  tryCatch({
  cs <- results$cs[[i]] %>%
              rename(CHROM = CHR, POS = BP) %>%
              separate(`-LOG10P`, into = paste0("logP_", ancestries), sep = ",") %>%
              mutate(P = as.numeric(logP_EUR))

  snp <- results$snp[[i]]

  # Get CHR:BP:BP for the fine mapped region
  CHR <- unique(cs$CHR)
  BP_START <- unique(cs$BP_START)
  BP_END <- unique(cs$BP_END)
  region <- paste0(CHR, ":", BP_START, ":", BP_END)

  # Region plot -log10(p) vs snp requires R^2 between SNPs
  # Get a list of SNPs that were included in the fine mapping
  # arrange by PIP so the first row is the SNP with highest PIP
  snp_list <- snp %>%
    arrange(desc(`PIP(CS1)`)) %>%
    pull(SNP)

  writeLines(snp_list, paste0("tmp/snp_list_", region, ".txt"))

  # Per ancestry:
  region_plots <- lapply(ancestries, function(ancestry){
    # Add r2 column to results, conduct in PLINK
    system2("plink2", # edit to make this a variable for where plink2 is stored locally
            args = c(
              "--bfile", paste0("reference/ukb_imp_v3.qc.geno02.mind02_", ancestry ,"_", CHR),
              "--r2-phased",
              "--ld-snp", snp_list[1],
              "--ld-window-kb", "1000",
              "--ld-window-r2", "0.2", # r^2 < 0.2 is filtered out
              "--extract", paste0("tmp/snp_list_", region, ".txt"),
              "--out", paste0("tmp/ld_results",ancestry, region, ".txt")
            )
    )

    # Define the path to the LD results file
    ld_results_file <- paste0("tmp/ld_results", ancestry, region, ".txt.vcor")

    # Check if the LD results file exists
    if (file.exists(ld_results_file)) {
      # If the file exists, load the LD results
      ld_results <- fread(ld_results_file)
      
      # Perform the usual joining logic
      joined_results <- snp %>%
        left_join(., ld_results, by = c("SNP" = "ID_B")) %>%
        left_join(., sumstats[[ancestry]], by = "SNP") %>%
        dplyr::select(-ends_with(".y")) %>%
        rename_with(~ str_remove(., "\\.x$"), ends_with(".x")) %>%
        left_join(cs, by = "SNP") %>%
        dplyr::select(-ends_with(".y")) %>%
        rename_with(~ str_remove(., "\\.x$"), ends_with(".x")) %>%
        dplyr::select(SNP, 
                      CHROM = CHR,
                      POS = BP,
                      P,
                      R2 = PHASED_R2,
                      PIP_CS1 = `PIP(CS1)`,
                      OVRL_PIP,
                      CS_PIP) %>%
        mutate(R2 = ifelse(SNP == snps[1], 1, R2)) %>%
        mutate(R2 = ifelse(is.na(R2), 0, R2)) %>%
        mutate(PIP_color = ifelse(PIP_CS1 == OVRL_PIP, "CS_SNP", NA))
    } else {
      # If the file does not exist, create a default joined_results with R2 set to NA for all rows
      joined_results <- snp %>%
        left_join(., sumstats[[ancestry]], by = "SNP") %>%
        dplyr::select(-ends_with(".y")) %>%
        rename_with(~ str_remove(., "\\.x$"), ends_with(".x")) %>%
        left_join(cs, by = "SNP") %>%
        dplyr::select(-ends_with(".y")) %>%
        rename_with(~ str_remove(., "\\.x$"), ends_with(".x")) %>%
        mutate(R2 = NA) %>%
        mutate(R2 = ifelse(SNP == snps[1], 1, R2)) %>%
        mutate(R2 = ifelse(is.na(R2), 0, R2)) %>%
        dplyr::select(SNP, 
                      CHROM = CHR,
                      POS = BP,
                      P,
                      R2,  # R2 column with NA for all rows
                      PIP_CS1 = `PIP(CS1)`,
                      OVRL_PIP,
                      CS_PIP) %>%
        mutate(PIP_color = ifelse(PIP_CS1 == OVRL_PIP, "CS_SNP", NA))
    }

    # Region plot (one per ancestry):
    region_plot <- topr::locuszoom(
      df = joined_results,
      chr = CHR,
      xmin = BP_START,
      xmax = BP_END,
      log_trans_p = TRUE,
      build = 37,
      alpha = 0.5,
      extract_plots = TRUE,
      title = ancestry
    )$main_plot

    return(region_plot)
  })

  # PIP plot (one for all ancestries):
  bounding_boxes <- joined_results %>%
    filter(PIP_color == "CS_SNP") %>%
    summarize(
      xmin = min(POS)-5000,
      xmax = max(POS)+5000,
      ymin = min(PIP_CS1)-0.01,
      ymax = max(PIP_CS1)+0.01
    )

  pip_plot <- ggplot(joined_results) +
    geom_point(aes(x = POS, y = PIP_CS1, color = PIP_color))+
    geom_rect(data = bounding_boxes, 
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                alpha = 0, color = "black", inherit.aes = FALSE) +
    labs(x = paste0("Position (BP) on CHR: ", CHR),
    y = "PIP",
    color = "SNP in credible set") +
    theme_minimal()


  png(paste0("fineMapping/plots/region_plot", region, ".png"), width = 2300, height = 2300, res = 300)
  cowplot::plot_grid(plotlist = c(region_plots, list(pip_plot)), ncol = 1)
  dev.off()

  }, error = function(e) {
    cat("Error in fine map region ID:", i, "., For region: ", region, "., and ancestry: ", ancestry, "\n")
    cat("Error message:", e$message, "\n")
  })

}