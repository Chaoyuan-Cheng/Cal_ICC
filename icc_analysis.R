# =============================================================================
# ICC Analysis for Species-level Relative Abundance Data
# Model: log_abund ~ s(Period_in_Diet, re) + s(Host, re) + factor(Diet)
# Random effects: Host (individual), Period_in_Diet (time position within diet phase)
# Fixed effect: Diet (HF vs LF)
#
# Two independent analyses:
#   - D group (DD): HF -> LF (Period 1-3 = HF, Period 4-6 = LF)
#   - R group (DR): LF -> HF (Period 1-3 = LF, Period 4-6 = HF)
#
# Period_in_Diet maps: Period 1,4->1; Period 2,5->2; Period 3,6->3
# =============================================================================

library(tidyverse)
library(mgcv)
library(ggthemes)
library(wesanderson)

setwd("/home/ps/Desktop/work/work_260129/R0202/data/")

# =============================================================================
# Parameters (adjust here)
# =============================================================================

# Data files
ABUNDANCE_FILE <- "relative abundance.csv"

# Two independent analyses with opposite diet transitions
ANALYSIS_CONFIG <- list(
  D = list(
    metadata_file = "1_metadata-D.txt",
    output_prefix = "ICC_species_D",
    description = "D group: HF -> LF"
  ),
  R = list(
    metadata_file = "1_metadata-R.txt",
    output_prefix = "ICC_species_R",
    description = "R group: LF -> HF"
  )
)

# Filtering
MIN_MEAN_ABUNDANCE <- 0.001    # Minimum mean relative abundance (0.001 = 0.1%)
TOP_N_SPECIES <- 25            # Number of top abundant species to plot

# Plot settings
PLOT_WIDTH <- 8                # Figure width (inches)
PLOT_HEIGHT_PER_SPECIES <- 0.25  # Height per species (inches)
PLOT_MIN_HEIGHT <- 8           # Minimum figure height (inches)
X_AXIS_LIMIT <- 1              # X-axis upper limit for ICC
SPECIES_FONT_SIZE <- 7         # Font size for species names
BASE_FONT_SIZE <- 11           # Base font size
LINE_WIDTH <- 1.6              # Error bar line width
POINT_SIZE <- 3.0              # Point size

# =============================================================================
# 1. Load abundance data (shared across analyses)
# =============================================================================

cat("Loading abundance data...\n")

abundance_raw <- read.csv(ABUNDANCE_FILE, row.names = 1, check.names = FALSE)
cat(sprintf("Raw data: %d species, %d samples\n", nrow(abundance_raw), ncol(abundance_raw)))

# =============================================================================
# 2. Define analysis function
# =============================================================================

run_icc_analysis <- function(abundance_raw, metadata_file, output_prefix, description) {

  cat(sprintf("\n=== Running analysis: %s ===\n", description))

  # Load metadata
  metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
  rownames(metadata) <- metadata$ID
  cat(sprintf("Metadata: %d samples\n", nrow(metadata)))

  # Check sample overlap
  common_samples <- intersect(colnames(abundance_raw), metadata$ID)
  cat(sprintf("Common samples: %d\n", length(common_samples)))

  # Subset to common samples
  abundance <- abundance_raw[, common_samples]
  metadata <- metadata[common_samples, ]

  # Convert factors
  metadata$Host <- as.factor(metadata$Host)
  metadata$Period <- as.factor(metadata$Period)
  metadata$Diet <- as.factor(metadata$Diet)
  metadata$Group <- as.factor(metadata$Group)
  metadata$Day <- as.factor(metadata$Day)

  # Create Period_in_Diet variable (time position within diet phase)
  # Period 1,4 -> 1; Period 2,5 -> 2; Period 3,6 -> 3
  # This separates time-within-diet variation from Diet effect
  metadata$Period_in_Diet <- as.factor(((as.numeric(as.character(metadata$Period)) - 1) %% 3) + 1)

  # Print data structure for verification
  cat("\nData structure verification:\n")
  print(table(metadata$Period, metadata$Diet))
  cat("\nPeriod_in_Diet mapping:\n")
  print(table(metadata$Period, metadata$Period_in_Diet))

  # Filter species by abundance
  cat("\nFiltering species...\n")
  mean_abund <- rowMeans(abundance)
  species_keep <- names(mean_abund[mean_abund > MIN_MEAN_ABUNDANCE])
  abundance_filtered <- abundance[species_keep, ]
  cat(sprintf("Filtered: %d species (mean abundance > %.1f%%)\n",
              nrow(abundance_filtered), MIN_MEAN_ABUNDANCE * 100))

  # Calculate ICC for each species
  ICC_list <- list()
  species_names <- rownames(abundance_filtered)

  cat(sprintf("\nCalculating ICC for %d species...\n", length(species_names)))
  pb <- txtProgressBar(min = 0, max = length(species_names), style = 3)

  for (i in seq_along(species_names)) {

    tryCatch({

      setTxtProgressBar(pb, i)
      sp_name <- species_names[i]

      # Prepare data for this species
      sp_data <- data.frame(
        abundance = as.numeric(abundance_filtered[sp_name, ]),
        metadata
      )

      # Log transform (add small value to handle zeros)
      sp_data$log_abund <- log10(sp_data$abundance + 1e-6)

      # Fit GAM model with random effects
      # Period_in_Diet has 3 levels (time position within each diet phase)
      # This is orthogonal to Diet effect (no longer confounded)
      gam_model <- mgcv::bam(
        log_abund ~
          s(Period_in_Diet, bs = "re") +
          s(Host, bs = "re") +
          factor(Diet),
        data = sp_data,
        family = gaussian,
        discrete = TRUE
      )

      # Get variance components from gam.vcomp (profile likelihood CI)
      vc_mat <- gam.vcomp(gam_model)
      sd_time <- vc_mat["s(Period_in_Diet)", "std.dev"]
      sd_time_lower <- vc_mat["s(Period_in_Diet)", "lower"]
      sd_time_upper <- vc_mat["s(Period_in_Diet)", "upper"]
      sd_host <- vc_mat["s(Host)", "std.dev"]
      sd_host_lower <- vc_mat["s(Host)", "lower"]
      sd_host_upper <- vc_mat["s(Host)", "upper"]
      sd_scale <- vc_mat["scale", "std.dev"]

      # Calculate variances
      var_time <- sd_time^2
      var_host <- sd_host^2
      var_resid <- sd_scale^2
      var_time_lower <- sd_time_lower^2
      var_time_upper <- sd_time_upper^2
      var_host_lower <- sd_host_lower^2
      var_host_upper <- sd_host_upper^2

      # Calculate total variance
      var_total <- var_time + var_host + var_resid

      # Calculate ICC point estimates directly from variance components
      icc_host <- var_host / var_total
      icc_time <- var_time / var_total

      # Calculate CI for ICC
      icc_time_lower <- var_time_lower / (var_time_lower + var_host + var_resid)
      icc_time_upper <- var_time_upper / (var_time_upper + var_host + var_resid)
      icc_host_lower <- var_host_lower / (var_host_lower + var_time + var_resid)
      icc_host_upper <- var_host_upper / (var_host_upper + var_time + var_resid)

      # Calculate p-values using likelihood ratio test (compare full vs reduced models)
      # Reduced model without Host random effect
      gam_no_host <- mgcv::bam(
        log_abund ~ s(Period_in_Diet, bs = "re") + factor(Diet),
        data = sp_data, family = gaussian, discrete = TRUE
      )
      # Reduced model without Period_in_Diet random effect
      gam_no_time <- mgcv::bam(
        log_abund ~ s(Host, bs = "re") + factor(Diet),
        data = sp_data, family = gaussian, discrete = TRUE
      )

      # LRT: 2 * (logLik_full - logLik_reduced), df = 1, chi-squared test
      lrt_host <- 2 * (logLik(gam_model) - logLik(gam_no_host))
      lrt_time <- 2 * (logLik(gam_model) - logLik(gam_no_time))
      p_host <- pchisq(as.numeric(lrt_host), df = 1, lower.tail = FALSE) / 2  # one-sided
      p_time <- pchisq(as.numeric(lrt_time), df = 1, lower.tail = FALSE) / 2  # one-sided

      # Extract Diet fixed effect
      summ <- summary(gam_model)
      diet_coef <- summ$p.table["factor(Diet)LF", ]
      diet_estimate <- diet_coef["Estimate"]
      diet_se <- diet_coef["Std. Error"]
      diet_p <- diet_coef["Pr(>|t|)"]

      # Calculate Diet effect as proportion of total variance (pseudo-ICC for visualization)
      # Use absolute standardized effect: |estimate| / total_sd
      total_sd <- sqrt(var_time + var_host + var_resid)
      diet_effect_std <- abs(diet_estimate) / total_sd
      # Cap at 1 for visualization
      diet_effect_std <- min(diet_effect_std, 1)

      # CI for Diet effect (using delta method approximation)
      diet_effect_lower <- max(0, (abs(diet_estimate) - 1.96 * diet_se) / total_sd)
      diet_effect_upper <- min(1, (abs(diet_estimate) + 1.96 * diet_se) / total_sd)

      # Get mean abundance for this species
      sp_mean_abund <- mean_abund[sp_name]

      # Create results dataframe
      result <- data.frame(
        ICC = c(icc_host, icc_time, diet_effect_std),
        lower = c(icc_host_lower, icc_time_lower, diet_effect_lower),
        upper = c(icc_host_upper, icc_time_upper, diet_effect_upper),
        effect = c("Host", "Time", "Diet"),
        p_val = c(p_host, p_time, diet_p),
        Species = sp_name,
        mean_abundance = sp_mean_abund
      )
      result$Significant <- result$p_val < 0.05

      # Ensure lower bound is not negative
      result$lower <- pmax(0, result$lower)

      ICC_list[[i]] <- result

    }, error = function(e) {
      message(paste("Error in species", sp_name, ":", e$message))
    })
  }

  close(pb)

  # Combine results
  icc_df <- do.call(rbind, ICC_list)
  cat(sprintf("\nSuccessfully calculated ICC for %d species\n", length(unique(icc_df$Species))))

  # Save full results CSV first (before filtering for plot)
  output_csv <- paste0(output_prefix, "_results.csv")
  write.csv(icc_df, output_csv, row.names = FALSE)
  cat(sprintf("Full results saved to: %s\n", output_csv))

  # Select top N species by mean abundance for plotting
  top_species <- icc_df %>%
    filter(effect == "Host") %>%
    arrange(-mean_abundance) %>%
    head(TOP_N_SPECIES) %>%
    pull(Species)

  icc_df_plot <- icc_df %>%
    filter(Species %in% top_species)

  cat(sprintf("Plotting top %d species by mean abundance\n", length(top_species)))

  # Prepare data for plotting (order by mean abundance, highest at top)
  abundance_order <- icc_df_plot %>%
    filter(effect == "Host") %>%
    arrange(-mean_abundance) %>%
    pull(Species)

  icc_df_plot$Species <- factor(icc_df_plot$Species, levels = rev(abundance_order))
  icc_df_plot$effect <- factor(icc_df_plot$effect, levels = c("Host", "Time", "Diet"))

  # Generate ICC figure (top N species only)
  cat("\nGenerating ICC figure...\n")

  pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

  # Split data by significance for cleaner plotting
  icc_sig <- icc_df_plot %>% filter(Significant == TRUE)
  icc_nonsig <- icc_df_plot %>% filter(Significant == FALSE)

  fig_icc <- ggplot(icc_df_plot, aes(y = Species, x = ICC)) +
    # Non-significant: grey, dotted, circle
    geom_errorbarh(data = icc_nonsig,
                   aes(xmin = lower, xmax = upper),
                   height = 0, linewidth = LINE_WIDTH, col = "grey60", linetype = "dotted") +
    geom_point(data = icc_nonsig,
               size = POINT_SIZE, col = "grey40", fill = "grey70", shape = 21) +
    # Significant: colored, solid, triangle
    geom_errorbarh(data = icc_sig,
                   aes(xmin = lower, xmax = upper, col = ICC),
                   height = 0, linewidth = LINE_WIDTH, linetype = "solid") +
    geom_point(data = icc_sig,
               aes(fill = ICC), size = POINT_SIZE, col = "grey40", shape = 24) +
    # Scales
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_gradientn(colours = pal) +
    scale_fill_gradientn(colours = pal) +
    # Facet by effect type
    facet_wrap(~effect, ncol = 3) +
    # Theme
    theme_clean(base_size = BASE_FONT_SIZE) +
    xlab("ICC / Standardized Effect") +
    ylab("") +
    ggtitle(paste0(description, " (Top ", TOP_N_SPECIES, " species)")) +
    labs(fill = "Effect Size", colour = "Effect Size") +
    theme(axis.text.y = element_text(size = SPECIES_FONT_SIZE)) +
    coord_cartesian(xlim = c(0, X_AXIS_LIMIT)) +
    theme(legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "grey80")) +
    theme(plot.background = element_rect(color = "white", fill = "white"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5, size = 14))

  # Save - adjust height based on number of species in plot
  n_species_plot <- length(unique(icc_df_plot$Species))
  fig_height <- max(PLOT_MIN_HEIGHT, n_species_plot * PLOT_HEIGHT_PER_SPECIES)

  # Save PDF
  ggsave(paste0(output_prefix, ".pdf"), fig_icc, width = PLOT_WIDTH, height = fig_height)

  # Save PNG
  png(paste0(output_prefix, ".png"), width = PLOT_WIDTH, height = fig_height, units = "in", res = 300)
  print(fig_icc)
  dev.off()

  cat(sprintf("\nFigure saved to:\n"))
  cat(sprintf("  - %s.pdf\n", output_prefix))
  cat(sprintf("  - %s.png\n", output_prefix))

  # Summary statistics (for all species)
  cat("\n========== Summary Statistics (All Species) ==========\n\n")

  summary_stats <- icc_df %>%
    group_by(effect) %>%
    summarise(
      mean_ICC = round(mean(ICC, na.rm = TRUE), 3),
      sd_ICC = round(sd(ICC, na.rm = TRUE), 3),
      median_ICC = round(median(ICC, na.rm = TRUE), 3),
      n_significant = sum(Significant, na.rm = TRUE),
      n_total = n(),
      pct_significant = round(sum(Significant, na.rm = TRUE) / n() * 100, 1)
    )

  print(summary_stats)
  cat("\n=======================================================\n")

  return(list(icc_df = icc_df, summary_stats = summary_stats))
}

# =============================================================================
# 3. Run analyses for both groups
# =============================================================================

results <- list()

for (group_name in names(ANALYSIS_CONFIG)) {
  config <- ANALYSIS_CONFIG[[group_name]]
  results[[group_name]] <- run_icc_analysis(
    abundance_raw = abundance_raw,
    metadata_file = config$metadata_file,
    output_prefix = config$output_prefix,
    description = config$description
  )
}

cat("\n\n========================================\n")
cat("All analyses completed!\n")
cat("========================================\n")
