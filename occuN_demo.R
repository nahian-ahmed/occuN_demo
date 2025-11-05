# -----------------------------------------------------------------
# Simulation for occuN model
#
# This script simulates data according to the occuN process
# and then fits the occuN model to check for parameter recovery.
# -----------------------------------------------------------------



# Quietly install forked 'unmarked' package
suppressMessages(
    devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE, quiet = FALSE)
)

# Quietly install plotting packages
suppressMessages(
  if (!require("ggplot2")) {
    install.packages("ggplot2", quiet = TRUE)
    library(ggplot2)
  }
)

# Quietly install patchwork
suppressMessages(
  if (!require("patchwork")) {
    install.packages("patchwork", quiet = TRUE)
    library(patchwork)
  }
)


# Load required libraries
library(unmarked)

##########
# 1. Set Simulation Parameters
##########

set.seed(123) # For reproducibility

# Landscape parameters
grid_dim <- 50 # Landscape is 50x50 cells
n_cells <- grid_dim * grid_dim # Total number of cells

# Site parameters (reference clustering)
site_dim <- 5 # Sites are 5x5 cell blocks (5-cellSq)
n_sites_x <- grid_dim / site_dim
n_sites_y <- grid_dim / site_dim
M <- n_sites_x * n_sites_y # Total number of sites (100)

# Observation parameters
J_obs <- 3 # Number of surveys per site

# True parameter values
# Detection (alphas) for formula ~obs_cov1
true_alphas <- c(alpha_int = 0.5, alpha_cov = -1.0)

# State (betas) for formula ~cell_cov1
true_betas <- c(beta_int = -4.0, beta_cov = 2.0) 


# optimizers = "BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "SANN", "nlminb" 
selected_optimizer <- "nlminb"



cat("--- Simulation Starting ---\n")
cat(sprintf("Simulating %d sites (from %d cells) with %d surveys.\n", M, n_cells, J_obs))
cat("True Alphas (Detection):", true_alphas, "\n")
cat("True Betas (State/Abundance):", true_betas, "\n")

##########
# 2. Create Landscape (Cells & Covariates)
##########

# Create cell-level covariates
cellCovs <- data.frame(
  cell_cov1 = rnorm(n_cells)
)



# Create design matrix for state (X_design)
X_cell <- model.matrix(~cell_cov1, data = cellCovs)

# Calculate true latent abundance (lambda_j) for each cell
log_lambda_j <- X_cell %*% true_betas
lambda_j <- exp(log_lambda_j)


##########
# 3. Define Sites & Create Weight Matrix (w)
##########

# This creates the non-overlapping 5x5 sites

# Assign each cell to a site
cell_row <- (0:(n_cells - 1) %/% grid_dim) + 1
cell_col <- (0:(n_cells - 1) %% grid_dim) + 1
site_row <- (cell_row - 1) %/% site_dim + 1
site_col <- (cell_col - 1) %/% site_dim + 1

# This vector has length n_cells, with values from 1 to M
site_id_for_cell <- (site_row - 1) * n_sites_x + site_col

# Create the weight matrix w (M x n_cells)
w <- matrix(0, M, n_cells)
for (i in 1:M) {
  w[i, site_id_for_cell == i] = 1
}
# In this simple non-overlapping case, w_ij is 1 if cell j is in site i, 0 otherwise.

##########
# 4. Simulate True State (Occupancy Z)
##########


# Calculate site-level latent abundance (lambda_tilde_i)
# \tilde{\lambda}_i = \sum_j w_{ij} \lambda_j
lambda_tilde_i <- w %*% lambda_j

# Calculate site-level occupancy probability (psi_i)
# \psi_i = 1 - e^{-\tilde{\lambda}_i}
psi_i <- 1 - exp(-lambda_tilde_i)

# Simulate true occupancy state (Z_i)
Z_i <- rbinom(M, 1, psi_i)

cat("\nSimulated true occupancy states:\n")
print(table(Z_i)) # <-- This should now show both 0s and 1s!

##########
# 5. Simulate Observation Data (y)
##########

# Create observation-level covariates
obs_cov1 <- matrix(rnorm(M * J_obs), M, J_obs)
obsCovs <- list(obs_cov1 = obs_cov1)

# Create the detection history matrix y
y <- matrix(NA, M, J_obs)

for (i in 1:M) {
  # If site is not occupied, detection is 0
  if (Z_i[i] == 0) {
    y[i, ] <- 0
    next
  }
  
  # If site is occupied, simulate detections
  for (k in 1:J_obs) {
    # Build linear predictor for detection
    logit_p_ik <- true_alphas[1] * 1 + true_alphas[2] * obsCovs$obs_cov1[i, k]
    
    # Convert to probability
    p_ik <- plogis(logit_p_ik)
    
    # Simulate detection
    y[i, k] <- rbinom(1, 1, p_ik)
  }
}

cat("\nSimulated detection history (first 6 sites):\n")
print(head(y))

##########
# 6. Bundle Data into unmarkedFrameOccuN
##########

# We need to provide all the pieces:
# y = detection data (M x J_obs)
# obsCovs = list of obs covariates (M x J_obs matrices)
# cellCovs = data.frame of cell covariates (n_cells x n_covs)
# w = weight matrix (M x n_cells)

umf <- unmarkedFrameOccuN(
  y = y,
  obsCovs = obsCovs,
  cellCovs = cellCovs,
  w = w
)

cat("\nCreated unmarkedFrameOccuN:\n")
print(umf)


##########
# 7. Fit the occuN Model
##########

cat("\nFitting occuN model...\n")

# The formula matches our simulation design:
# ~obs_cov1 for detection (alphas)
# ~cell_cov1 for state (betas)

# Provide starting values to help the optimizer
starts <- c(true_alphas, true_betas)

# Fit occuN model
fm <- occuN(
  formula = ~obs_cov1 ~ cell_cov1,
  data = umf,
  starts = starts,
  se = TRUE,
  method = selected_optimizer
)

##########
# 8. Compare Results
##########

cat("\n--- Simulation Complete ---")
cat("\n\n--- Parameter Recovery Check ---\n\n")

# Create a simplified comparison table (no SEs)
results <- data.frame(
  Parameter = c(
    "alpha (det_int)", "alpha (det_cov1)",
    "beta (state_int)", "beta (state_cov1)"
  ),
  True_Value = c(true_alphas, true_betas),
  Estimated_Value = c(coef(fm, 'det'), coef(fm, 'state'))
)

print(results)

write.csv(results, "params.csv", row.names = FALSE)

##########
# 9. Plot Landscape and Sites
##########

# Create data frame for 50x50 cell covariates
cell_df <- data.frame(
  x = cell_col,
  y = cell_row,
  covariate = cellCovs$cell_cov1
)

cat("\n\nGenerating plots...\n")

# Map site-level results back to the cell_df for plotting
cell_df$site_latent_abundance <- as.numeric(lambda_tilde_i[site_id_for_cell])
cell_df$site_occupancy_prob <- as.numeric(psi_i[site_id_for_cell])
cell_df$site_true_occupancy <- as.factor(Z_i[site_id_for_cell]) # Factor for discrete colors

# Create data frame for site boundary boxes
# We have n_sites_x (10) and n_sites_y (10)
# And site_dim (5)
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
  
  # Add red boxes for sites
  geom_rect(data = site_boxes, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "red", 
            fill = NA, # NA fill = transparent
            linewidth = 0.5,
            inherit.aes = FALSE) + # IMPORTANT: stops geom_rect from using the 'fill' aesthetic
  
  labs(title = "Cell Covariate and 5-cellSq Sites",
       fill = "Covariate" ) +
  theme_minimal()


# Plot 2: Site-Level Latent Abundance
p_abundance <- ggplot(cell_df, aes(x = x, y = y, fill = site_latent_abundance)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") + # Use a different palette
  coord_fixed(expand = FALSE) +
  # Add red boxes for sites
  geom_rect(data = site_boxes, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "red", 
            fill = NA, 
            linewidth = 0.5,
            inherit.aes = FALSE) +
  labs(title = "Site Abundance",
       fill = "Latent Abund.") +
  theme_minimal()

# Plot 3: True Simulated Occupancy State
p_occupancy <- ggplot(cell_df, aes(x = x, y = y, fill = site_true_occupancy)) +
  geom_raster() +
  scale_fill_manual(values = c("0" = "navyblue", "1" = "yellow")) + # Clear discrete colors
  coord_fixed(expand = FALSE) +
  # Add red boxes for sites
  geom_rect(data = site_boxes, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "red", 
            fill = NA, 
            linewidth = 0.5,
            inherit.aes = FALSE) +
  labs(title = "Site Occupancy",
       fill = "Occupied") +
  theme_minimal()

# Combine and Save Plots

# Use patchwork to combine plots in a single row
combined_plot <- p_covariate + p_abundance + p_occupancy + 
  plot_layout(nrow = 1)


# Save the combined plot
ggsave("plot.png", 
       plot = combined_plot, 
       dpi = 300, 
       width = 18, 
       height = 6)

cat("\n--- Plot saved to plot.png ---\n")
cat("--- Script Finished ---\n")