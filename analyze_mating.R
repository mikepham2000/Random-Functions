`analyze_mating` <-
function(file_path, farm_code, breed_code, 
                                      inbreeding_threshold = 0.0625,
                                      output_dir = NULL) {
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- dirname(file_path)
  }
  
  # Read data
  cat("Reading Excel file...\n")
  dt <- read_excel(file_path, sheet = 1)
  
  # Remove duplicates
  dt <- dt[!duplicated(dt$`Mã số heo`), ]
  
  # Rename columns to English for easier handling
  dt <- dt %>%
    rename(
      Code = Code,
      Animal_ID = `Mã số heo`,
      Sex = `Giới tính`,
      Breed = Giống,
      Age = `Tuần tuổi`,
      Sire = `Mã số heo cha`,
      Dam = `Mã số heo mẹ`,
      Status = `Hiện trạng`,
      ebv_nba = `NBA (ebv)`,
      ebv_fi = `FI (ebv)`,
      ebv_adg = `ADG (ebv)`,
      ebv_bf = `BF (ebv)`,
      Index = Index,
      group = Grp,
      Marker = Marker
    )
  
  # Show available farms and breeds
  cat("\nAvailable farms:\n")
  print(table(dt$Code))
  
  cat("\nAvailable breeds:\n")
  print(table(dt$Breed))
  
  # Filter for specific farm and breed
  # For pedigree construction, use all animals; for mating, apply status filters
  dt_pedigree <- dt %>%
    filter(Code == farm_code, 
           Breed == breed_code)
  
  dt_sires <- dt %>%
    filter(Code == farm_code, 
           Breed == breed_code,
           Sex == "1",
           Status == "Onfarm")  # Only Onfarm sires
  
  dt_dams <- dt %>%
    filter(Code == farm_code, 
           Breed == breed_code,
           Sex == "2",
           Status %in% c("Onfarm", "Piglet"))  # Onfarm and Piglet dams
  
  if (nrow(dt_sires) == 0 || nrow(dt_dams) == 0) {
    stop("No sires or dams found for analysis")
  }
  
  cat(sprintf("\nAnalyzing farm %s, breed %s\n", farm_code, breed_code))
  cat(sprintf("Total animals for pedigree: %d\n", nrow(dt_pedigree)))
  cat(sprintf("Sires (Onfarm only): %d\n", nrow(dt_sires)))
  cat(sprintf("Dams (Onfarm + Piglet): %d\n", nrow(dt_dams)))
  
  # Create and edit pedigree using all animals
  cat("\nCreating pedigree structure...\n")
  ped <- editPed(sire = dt_pedigree$Sire,
                 dam = dt_pedigree$Dam,
                 label = dt_pedigree$Animal_ID)
  
  cat(sprintf("Pedigree includes %d animals (with ancestors)\n", nrow(ped)))
  cat(sprintf("Maximum generation depth: %d\n", max(ped$gene)))
  
  # Create pedigree object
  ped.complete <- pedigree(sire = ped$sire, 
                          dam = ped$dam, 
                          label = ped$label)
  
  # Calculate A matrix
  cat("\nCalculating relationship matrix (A)...\n")
  A <- getA(ped.complete)
  
  # Get sires and dams IDs
  sire_id <- dt_sires %>% pull(Animal_ID)
  dam_id <- dt_dams %>% pull(Animal_ID)
  
  cat(sprintf("\n%d sires (Onfarm) | %d dams (Onfarm+Piglet) on farm\n", 
              length(sire_id), length(dam_id)))
  
  if (length(sire_id) == 0 || length(dam_id) == 0) {
    stop("No sires or dams found for analysis")
  }
  
  # Get marker information
  sire_markers <- dt_sires %>% select(Animal_ID, Marker)
  dam_markers <- dt_dams %>% select(Animal_ID, Marker)
  
  # Extract relationship submatrix
  relationship <- as.matrix(t(A[sire_id, dam_id]))
  
  # Convert to data frame and add marker information
  rel_table <- as.data.frame(as.table(relationship)) %>%
    rename(Dam = Var1, Sire = Var2, Rel = Freq) %>%
    left_join(dam_markers, by = c("Dam" = "Animal_ID")) %>%
    rename(Dam_Marker = Marker) %>%
    left_join(sire_markers, by = c("Sire" = "Animal_ID")) %>%
    rename(Sire_Marker = Marker)
  
  # Calculate expected offspring inbreeding
  rel_table <- rel_table %>%
    mutate(Offspring_Inbreeding = Rel / 2,
           # Check for AA x AA exception at 0.125 relationship
           Marker_Exception = (Dam_Marker == "AA" & Sire_Marker == "AA" & Rel == 0.125))
  
  # Create summary statistics
  cat("\n=== INBREEDING STATISTICS ===\n")
  cat(sprintf("Mean expected offspring inbreeding: %.3f%%\n", 
              mean(rel_table$Offspring_Inbreeding) * 100))
  cat(sprintf("Max expected offspring inbreeding: %.3f%%\n", 
              max(rel_table$Offspring_Inbreeding) * 100))
  cat(sprintf("Min expected offspring inbreeding: %.3f%%\n", 
              min(rel_table$Offspring_Inbreeding) * 100))
  
  # Identify restricted matings (including marker exceptions for display)
  mating_restricted <- rel_table %>%
    filter(Offspring_Inbreeding >= inbreeding_threshold | Marker_Exception) %>%
    arrange(Dam, desc(Rel))  # Sort by relationship, not offspring inbreeding
  
  # Separate truly restricted from exceptions
  truly_restricted <- mating_restricted %>%
    filter(!Marker_Exception)
  
  aa_exceptions <- mating_restricted %>%
    filter(Marker_Exception)
  
  cat(sprintf("\nTotal possible matings: %d\n", nrow(rel_table)))
  cat(sprintf("Matings exceeding threshold (%.2f%%): %d (%.1f%%)\n", 
              inbreeding_threshold * 100,
              nrow(truly_restricted),
              nrow(truly_restricted) / nrow(rel_table) * 100))
  
  if (nrow(aa_exceptions) > 0) {
    cat(sprintf("AA × AA matings with relationship 0.125 (allowed): %d\n", nrow(aa_exceptions)))
  }
  
  # Create wide format table for restricted matings
  if (nrow(mating_restricted) > 0) {
    # Group by dam and create wide format
    # Add asterisk to sire IDs that are AA exceptions
    mating_tab <- mating_restricted %>%
      mutate(Sire_Display = ifelse(Marker_Exception, paste0(Sire, "*"), Sire)) %>%
      group_by(Dam) %>%
      mutate(rank = row_number()) %>%
      ungroup() %>%
      select(Dam, rank, Sire_Display) %>%
      pivot_wider(
        id_cols = Dam,
        names_from = rank,
        values_from = Sire_Display,
        names_prefix = "Sire"
      )
    
    # Extract last 4 digits for sorting
    mating_tab <- mating_tab %>% 
      mutate(Dam_id = as.integer(str_extract(Dam, "\\d{4}$")) %||% 0) %>% 
      arrange(Dam_id) %>%
      select(-Dam_id)
    
    # Save results in the specified format
    output_file <- file.path(output_dir, 
                           sprintf("Cam_phoi_%s_breed%s.csv", 
                                   farm_code, breed_code))
    write.csv(mating_tab, output_file, row.names = FALSE, na = "")
    cat(sprintf("\nRestricted matings saved to: %s\n", output_file))
    
    # Also save AA exceptions separately if they exist
    if (nrow(aa_exceptions) > 0) {
      aa_file <- file.path(output_dir, 
                          sprintf("AA_exceptions_%s_breed%s.csv", 
                                  farm_code, breed_code))
      write.csv(aa_exceptions %>% 
                select(Dam, Dam_Marker, Sire, Sire_Marker, Rel, Offspring_Inbreeding), 
                aa_file, row.names = FALSE)
      cat(sprintf("AA × AA exceptions saved to: %s\n", aa_file))
    }
    
    # Create visualization
    plot_file <- file.path(output_dir, 
                         sprintf("inbreeding_distribution_%s_breed%s.png", 
                                 farm_code, breed_code))
    
    p <- ggplot(rel_table, aes(x = Offspring_Inbreeding * 100)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = inbreeding_threshold * 100, 
                 color = "red", linetype = "dashed", size = 1) +
      labs(title = sprintf("Distribution of Expected Offspring Inbreeding\nFarm: %s, Breed: %s", 
                          farm_code, breed_code),
           x = "Expected Offspring Inbreeding (%)",
           y = "Frequency") +
      theme_minimal() +
      annotate("text", x = inbreeding_threshold * 100, 
               y = Inf, label = "Threshold", 
               vjust = 2, hjust = -0.1, color = "red")
    
    ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
    cat(sprintf("Distribution plot saved to: %s\n", plot_file))
    
  } else {
    cat("\nNo matings exceed the inbreeding threshold.\n")
    mating_tab <- NULL
  }
  
  # Return results
  results <- list(
    farm = farm_code,
    breed = breed_code,
    threshold = inbreeding_threshold,
    n_sires = length(sire_id),
    n_dams = length(dam_id),
    n_total_matings = nrow(rel_table),
    n_restricted = nrow(mating_restricted),
    relationship_table = rel_table,
    restricted_matings = mating_restricted,
    restricted_matings_wide = mating_tab,
    A_matrix = A
  )
  
  invisible(results)
}
