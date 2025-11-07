# -----------------------------------------------------------------
# Simulation for occuN model
#
# -----------------------------------------------------------------

# --- Installation ---
# Set a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

packages <- c("devtools", "ggplot2", "patchwork")
for(p in packages){
    if (!requireNamespace(p, quietly = FALSE)) install.packages(p)
}

# Quietly install forked 'unmarked' package
suppressMessages(
    devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE, quiet = FALSE)
)

# --- Load required libraries ---
library(unmarked)
library(ggplot2)
library(patchwork)

##########
# 1. Set Simulation Parameters
##########

set.seed(123) # For reproducibility

# --- Simulation repetitions ---
n_sims <- 100 # Number of full datasets to generate

# --- Model fitting repetitions ---
n_reps <- 30 # Number of random-start repetitions for each model fit

# --- Full Landscape parameters (200x200) ---
full_grid_dim <- 200 # Landscape is 200x200 cells
full_n_cells <- full_grid_dim * full_grid_dim # Total number of cells (40000)

# --- Site parameters (reference clustering) ---
site_dim <- 5 # Sites are 5x5 cell blocks (5-cellSq)
full_n_sites_x <- full_grid_dim / site_dim # 40
full_n_sites_y <- full_grid_dim / site_dim # 40
full_M <- full_n_sites_x * full_n_sites_y # Total number of sites (1600)

# --- Observation parameters ---
J_obs <- 3 # Number of surveys per site

# --- True parameter values ---
# Detection (alphas) for formula ~obs_cov1
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)

# State (betas) for formula ~cell_cov1
true_betas <- c(beta_int = -4.0, beta_cov = 2.0) 

# --- Model settings ---
# optimizers = "BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "SANN", "nlminb" 
selected_optimizer <- "nlminb"

# --- Ablation Study Parameters ---
M_values_to_test <- c(100, 225, 400, 900, 1600)


cat("--- Simulation Starting ---\n")
cat(sprintf("Running %d full simulations.\n", n_sims))
cat(sprintf("Each simulation tests %d M-values.\n", length(M_values_to_test)))
cat(sprintf("Each model fit uses %d random-start repetitions.\n", n_reps))
cat(sprintf("TOTAL MODEL FITS: %d\n\n", n_sims * length(M_values_to_test) * n_reps))
cat("True Alphas (Detection):", true_alphas, "\n")
cat("True Betas (State/Abundance):", true_betas, "\n")


##########
# 2. Define FULL Sites & Create Weight Matrix (w)
# (This is static geometry, so it lives OUTSIDE the loop)
##########

cat("Generating static full landscape geometry (1600 sites from 40000 cells)...\n")

# Assign each cell (1 to 40000) to a site (1 to 1600)
full_cell_row <- (0:(full_n_cells - 1) %/% full_grid_dim) + 1
full_cell_col <- (0:(full_n_cells - 1) %% full_grid_dim) + 1
full_site_row <- (full_cell_row - 1) %/% site_dim + 1
full_site_col <- (full_cell_col - 1) %/% site_dim + 1

# This vector has length 40000, with values from 1 to 1600
full_site_id_for_cell <- (full_site_row - 1) * full_n_sites_x + full_site_col

# Create the full weight matrix w (1600 x 40000)
full_w <- matrix(0, full_M, full_n_cells)
for (i in 1:full_M) {
    full_w[i, full_site_id_for_cell == i] = 1
}

##########
# 3. Helper Function for Subsetting
# (This is a static helper function, lives OUTSIDE the loop)
##########

#' Get Subset Landscape
#'
#' Slices the top-left quadrant of a full landscape to match a target
#' number of sites (M_target).
#'
#' @param M_target The desired number of sites (e.g., 100, 225)
#' @param site_dim Dimension of a single site (e.g., 5 for 5x5)
#' @param full_cell_row Vector of row numbers for all cells
#' @param full_cell_col Vector of col numbers for all cells
#' @param full_cellCovs Data frame of covariates for all cells
#' @param full_site_id_for_cell Vector mapping all cells to their site ID
#' @param full_w The full weight matrix (full_M x full_n_cells)
#' @param full_lambda_j Vector of latent abundance for all cells
#' @return A list containing all necessary subsetted data structures for occuN.
get_subset_landscape <- function(M_target, site_dim, full_cell_row, full_cell_col, full_cellCovs, full_site_id_for_cell, full_w, full_lambda_j) {
    
    # 1. Calculate dimensions of the target subset
    n_sites_per_dim_target <- sqrt(M_target)
    grid_dim_target <- n_sites_per_dim_target * site_dim
    
    cat(sprintf("\n--- Subsetting for M = %d ---\n", M_target))
    
    # 2. Find which *full* cell indices are in the top-left quadrant
    cell_indices_to_keep <- which(full_cell_row <= grid_dim_target & full_cell_col <= grid_dim_target)
    
    # 3. Find which *full* site indices are in the top-left quadrant
    site_ids_in_subset <- unique(full_site_id_for_cell[cell_indices_to_keep])
    site_ids_to_keep <- sort(site_ids_in_subset) # Should be 1:M_target
    
    if (length(site_ids_to_keep) != M_target) {
        warning(sprintf("Expected %d sites, but found %d unique sites in subset.", M_target, length(site_ids_to_keep)))
    }
    
    # --- 4. Subset all the data ---
    
    # Subset cellCovs
    sub_cellCovs <- full_cellCovs[cell_indices_to_keep, , drop = FALSE]
    rownames(sub_cellCovs) <- NULL
    
    # Subset weight matrix 'w'
    # Rows = sites to keep, Cols = cells to keep
    sub_w <- full_w[site_ids_to_keep, cell_indices_to_keep]
    
    # Create a re-mapped site ID vector for the subsetted cells
    # This maps global site IDs (e.g., 1, 2, 41, 42) to local IDs (1, 2, 3, 4)
    global_site_ids_for_subset_cells <- full_site_id_for_cell[cell_indices_to_keep]
    local_site_id_map <- 1:M_target
    names(local_site_id_map) <- as.character(site_ids_to_keep)
    sub_site_id_for_cell <- as.numeric(local_site_id_map[as.character(global_site_ids_for_subset_cells)])
    
    # Subset cell coordinates (relative to the new grid)
    sub_cell_row <- full_cell_row[cell_indices_to_keep]
    sub_cell_col <- full_cell_col[cell_indices_to_keep]
    
    # Subset lambda_j
    sub_lambda_j <- full_lambda_j[cell_indices_to_keep]
    
    # Return everything in a list
    return(list(
        M = M_target,
        n_cells = length(cell_indices_to_keep),
        grid_dim = grid_dim_target,
        n_sites_x = n_sites_per_dim_target,
        n_sites_y = n_sites_per_dim_target,
        cellCovs = sub_cellCovs,
        w = sub_w,
        site_id_for_cell = sub_site_id_for_cell,
        cell_row = sub_cell_row,
        cell_col = sub_cell_col,
        lambda_j = sub_lambda_j
    ))
}

##########
# 4. Initialize Loop & Storage
##########

all_results_df <- data.frame()
all_plots_list <- list() # Only for sim == 1

cat("\n--- Starting Main Simulation Loop ---\n")

# --- Main loop over SIMULATIONS ---
for (sim in 1:n_sims) {
    
    cat(sprintf("\n===========================================\n"))
    cat(sprintf("=== STARTING SIMULATION %d of %d ===\n", sim, n_sims))
    cat(sprintf("===========================================\n"))

    ##########
    # 5. Create FULL Landscape (Cells & Covariates)
    # (This is random, so it's INSIDE the sim loop)
    ##########
    
    cat("Generating new random cell covariates for full landscape...\n")
    
    # Create cell-level covariates for all 40000 cells
    full_cellCovs <- data.frame(
        cell_cov1 = rnorm(full_n_cells)
    )
    
    # Create design matrix for state (X_design) for all cells
    full_X_cell <- model.matrix(~cell_cov1, data = full_cellCovs)
    
    # Calculate true latent abundance (lambda_j) for each cell
    full_log_lambda_j <- full_X_cell %*% true_betas
    full_lambda_j <- exp(full_log_lambda_j)


    # --- Main loop over M values ---
    for (M_i in M_values_to_test) {
        
        # --- 6.1. Get subset data for M = M_i ---
        subset_data <- get_subset_landscape(
            M_target = M_i,
            site_dim = site_dim,
            full_cell_row = full_cell_row,
            full_cell_col = full_cell_col,
            full_cellCovs = full_cellCovs,
            full_site_id_for_cell = full_site_id_for_cell,
            full_w = full_w,
            full_lambda_j = full_lambda_j
        )
        
        # --- 6.2. Unpack the list into current simulation variables ---
        M <- subset_data$M
        n_cells <- subset_data$n_cells
        grid_dim <- subset_data$grid_dim
        n_sites_x <- subset_data$n_sites_x
        n_sites_y <- subset_data$n_sites_y
        cellCovs <- subset_data$cellCovs
        w <- subset_data$w
        site_id_for_cell <- subset_data$site_id_for_cell
        cell_row <- subset_data$cell_row
        cell_col <- subset_data$cell_col
        lambda_j <- subset_data$lambda_j
        
        cat(sprintf("Running simulation for M=%d (%d sites from %d cells).\n", M, M, n_cells))
        
        ##########
        # 7. Simulate True State (Occupancy Z)
        # (Random, so inside M loop)
        ##########
        
        # Calculate site-level latent abundance (lambda_tilde_i)
        # \tilde{\lambda}_i = \sum_j w_{ij} \lambda_j
        lambda_tilde_i <- w %*% lambda_j
        
        # Calculate site-level occupancy probability (psi_i)
        # \psi_i = 1 - e^{-\tilde{\lambda}_i}
        psi_i <- 1 - exp(-lambda_tilde_i)
        
        # Simulate true occupancy state (Z_i)
        Z_i <- rbinom(M, 1, psi_i)
        
        # cat("Simulated true occupancy states:\n")
        # print(table(Z_i))
        
        ##########
        # 8. Simulate Observation Data (y)
        # (Random, so inside M loop)
        ##########
        
        # Create observation-level covariates
        obs_cov1 <- matrix(rnorm(M * J_obs), M, J_obs)
        obsCovs <- list(obs_cov1 = obs_cov1)
        
        # Create the detection history matrix y
        y <- matrix(NA, M, J_obs)
        
        for (i in 1:M) {
            if (Z_i[i] == 0) {
                y[i, ] <- 0
                next
            }
            for (k in 1:J_obs) {
                logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * obsCovs$obs_cov1[i, k]
                p_ik <- plogis(logit_p_ik)
                y[i, k] <- rbinom(1, 1, p_ik)
            }
        }
        
        ##########
        # 9. Bundle Data into unmarkedFrameOccuN
        ##########
        
        umf <- unmarkedFrameOccuN(
            y = y,
            obsCovs = obsCovs,
            cellCovs = cellCovs,
            w = w
        )
        
        ##########
        # 10. Fit the occuN Model (with repetitions)
        ##########
        
        cat(sprintf("Fitting occuN model for M=%d (running %d reps)...\n", M, n_reps))
        
        best_fm <- NULL
        min_nll <- Inf
        n_params <- length(true_alphas) + length(true_betas)
    
        for (rep in 1:n_reps) {
            
            # Generate random starts from a uniform distribution
            rand_starts <- runif(n_params, -5, 5) 
            
            # Try to fit the model, suppress warnings/errors
            fm_rep <- try(occuN(
                formula = ~obs_cov1 ~ cell_cov1,
                data = umf,
                starts = rand_starts, # Use random starts
                se = TRUE,
                method = selected_optimizer
            ), silent = TRUE)
            
            # Check if the fit was successful
            if (inherits(fm_rep, "try-error")) {
                # if(rep %% 10 == 0) cat(sprintf("    (Rep %d failed to converge)\n", rep))
                next
            }
            
            # Check if this fit is better than the current best
            current_nll <- fm_rep@negLogLike
            if (current_nll < min_nll) {
                min_nll <- current_nll
                best_fm <- fm_rep
            }
        } # --- End of rep loop ---
        
        # Rename the best model to 'fm' to match existing code
        fm <- best_fm
        
        # Check if any model successfully converged
        if (is.null(fm)) {
                cat(sprintf("\n!!! ALL %d REPS FAILED for M=%d, sim=%d. Skipping. !!!\n", n_reps, M, sim))
                
                # Add NA results so this run is still recorded
                loop_results <- data.frame(
                    Parameter = c(
                        "alpha (det_int)", "alpha (det_cov1)",
                        "beta (state_int)", "beta (state_cov1)"
                    ),
                    True_Value = c(true_alphas, true_betas),
                    Estimated_Value = c(NA, NA, NA, NA),
                    M = M,
                    sim_id = sim
                )
                all_results_df <- rbind(all_results_df, loop_results)

                # If it's sim 1, add empty plots to keep the grid aligned
                if(sim == 1) {
                        all_plots_list <- c(all_plots_list, list(ggplot(), ggplot(), ggplot()))
                }
                next # Skip to the next M_i
        }
        
        cat(sprintf("    Best model for M=%d, sim=%d found with NLL: %.2f\n", M, sim, min_nll))
    
        ##########
        # 11. Store Results
        ##########
        
        loop_results <- data.frame(
            Parameter = c(
                "alpha (det_int)", "alpha (det_cov1)",
                "beta (state_int)", "beta (state_cov1)"
            ),
            True_Value = c(true_alphas, true_betas),
            Estimated_Value = c(coef(fm, 'det'), coef(fm, 'state')),
            M = M, # Add the M column
            sim_id = sim # Add the sim_id column
        )
        
        all_results_df <- rbind(all_results_df, loop_results)
        
        ##########
        # 12. Generate and Store Plots (ONLY FOR SIM 1)
        ##########
        
        if (sim == 1) {
            cat("    (sim=1, generating landscape plots...)\n")
            
            # Create data frame for cell covariates
            cell_df <- data.frame(
                x = cell_col,
                y = cell_row,
                covariate = cellCovs$cell_cov1
            )
            
            # Map site-level results back to the cell_df for plotting
            cell_df$site_latent_abundance <- as.numeric(lambda_tilde_i[site_id_for_cell])
            cell_df$site_occupancy_prob <- as.numeric(psi_i[site_id_for_cell])
            cell_df$site_true_occupancy <- as.factor(Z_i[site_id_for_cell])
            
            # Create data frame for site boundary boxes
            xmin_vals <- seq(from = 0.5, to = grid_dim - site_dim + 0.5, by = site_dim)
            xmax_vals <- seq(from = site_dim + 0.5, to = grid_dim + 0.5, by = site_dim)
            ymin_vals <- xmin_vals
            ymax_vals <- xmax_vals
            
            site_boxes <- expand.grid(xmin = xmin_vals, ymin = ymin_vals)
            site_boxes$xmax <- site_boxes$xmin + site_dim
            site_boxes$ymax <- site_boxes$ymin + site_dim
            
            
            # Plot 1: Original Cell Covariate
            p_covariate <- ggplot(cell_df, aes(x = x, y = y, fill = covariate)) +
                geom_raster() +
                scale_fill_viridis_c() +
                coord_fixed(expand = FALSE) +
                geom_rect(data = site_boxes, 
                                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                    color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                labs(title = sprintf(" Sites and Covariate (M = %d)", M), fill = "Covariate") +
                theme_minimal()
            
            # Plot 2: Site-Level Latent Abundance
            p_abundance <- ggplot(cell_df, aes(x = x, y = y, fill = site_latent_abundance)) +
                geom_raster() +
                scale_fill_viridis_c(option = "magma") +
                coord_fixed(expand = FALSE) +
                geom_rect(data = site_boxes, 
                                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                    color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                labs(title = sprintf("Site Abundance (M = %d)", M), fill = "Latent Abund.") +
                theme_minimal()
            
            # Plot 3: True Simulated Occupancy State
            p_occupancy <- ggplot(cell_df, aes(x = x, y = y, fill = site_true_occupancy)) +
                geom_raster() +
                scale_fill_manual(values = c("0" = "navyblue", "1" = "yellow")) +
                coord_fixed(expand = FALSE) +
                geom_rect(data = site_boxes, 
                                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                    color = "red", fill = NA, linewidth = 0.5, inherit.aes = FALSE) +
                labs(title = sprintf("Site Occupancy (M = %d)", M), fill = "Occupied") +
                theme_minimal()
            
            # --- Store plots for this iteration ---
            all_plots_list <- c(all_plots_list, list(p_covariate, p_abundance, p_occupancy))
        } # --- End of if(sim == 1) for plots ---
        
    } # --- End of for (M_i) loop ---
    
} # --- End of for (sim) loop ---

##########
# 13. Save Aggregate Results
##########

cat("\n--- Simulation Study Complete ---\n")

output_dir <- "output"
if (!dir.exists(output_dir)) {
        dir.create(output_dir)
}


# Save the full results data frame
write.csv(all_results_df, file.path(output_dir, "params.csv"), row.names = FALSE)

cat(sprintf("\n--- All parameters saved to %s/params.csv ---\n", output_dir))


##########
# 14. Combine and Save Landscape Plots (from Sim 1)
##########

cat("\nGenerating combined 5x3 landscape plot grid (from Sim 1)...\n")

# Combine all plots into a 5-row (for M) by 3-column (for plot type) grid
combined_plot <- patchwork::wrap_plots(all_plots_list, 
                                                                             nrow = length(M_values_to_test), 
                                                                             ncol = 3)

# Save the combined plot (needs to be tall)
ggsave(file.path(output_dir, "plot.png"), 
             plot = combined_plot, 
             dpi = 300, 
             width = 18, 
             height = 26) 

cat(sprintf("\n--- Combined landscape plot saved to %s/plot.png ---\n", output_dir))

##########
# 15. Create and Save Error Boxplots
##########

cat("\nGenerating 2x2 error boxplot grid (all sims)...\n")

# Calculate error
all_results_df$Error <- all_results_df$True_Value - all_results_df$Estimated_Value

# Convert M to factor for plotting
all_results_df$M_factor <- as.factor(all_results_df$M)

# Plot 1: State Intercept
p_err_beta_int <- ggplot(all_results_df[all_results_df$Parameter == "beta (state_int)", ], 
                                                 aes(x = M_factor, y = Error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "State Intercept", 
             x = "M (Sites)", 
             y = "Error (True - Est.)") +
    theme_bw()

# Plot 2: State Slope
p_err_beta_cov <- ggplot(all_results_df[all_results_df$Parameter == "beta (state_cov1)", ],
                                                 aes(x = M_factor, y = Error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "State Slope", 
             x = "M (Sites)", 
             y = "Error (True - Est.)") +
    theme_bw()

# Plot 3: Detection Intercept
p_err_alpha_int <- ggplot(all_results_df[all_results_df$Parameter == "alpha (det_int)", ],
                                                    aes(x = M_factor, y = Error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Detection Intercept", 
             x = "M (Sites)", 
             y = "Error (True - Est.)") +
    theme_bw()

# Plot 4: Detection Slope
p_err_alpha_cov <- ggplot(all_results_df[all_results_df$Parameter == "alpha (det_cov1)", ],
                                                    aes(x = M_factor, y = Error)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Detection Slope", 
             x = "M (Sites)", 
             y = "Error (True - Est.)") +
    theme_bw()

# Combine the 4 error plots into a 2x2 grid
# (State Intercept | State Slope) / (Detection Intercept | Detection Slope)
combined_error_plot <- (p_err_beta_int | p_err_beta_cov) / (p_err_alpha_int | p_err_alpha_cov)

# Save the combined error plot
ggsave(file.path(output_dir, "error_boxplots.png"), 
             plot = combined_error_plot, 
             dpi = 300, 
             width = 12, 
             height = 10)

cat(sprintf("\n--- Error boxplots saved to %s/error_boxplots.png ---\n", output_dir))
cat("--- Script Finished ---\n")