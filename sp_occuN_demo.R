#########################################################
# Script for occuN Clustering Comparison Plot
#
# (UPDATED BASED ON RECOMMENDATIONS)
#
#########################################################


######################################################
# 1. SETUP: Install and load all necessary packages
######################################################

options(repos = c(CRAN = "https://cloud.r-project.org/"))
# Added 'aricode' for ARI and 'PRROC' for AUROC calculation
packages <- c("devtools", "sf", "terra", "dplyr", "scales", "viridis", "ClustGeo", "aricode", "PRROC")
for(p in packages){
    if (!requireNamespace(p, quietly = FALSE)) install.packages(p)
}

# Quietly install forked 'unmarked' package
suppressMessages(
    devtools::install_github("nahian-ahmed/unmarked", ref = "occuN", force = TRUE, quiet = FALSE)
)

library(unmarked)
library(sf)
library(terra)
library(dplyr)
library(scales)
library(viridis)
library(ClustGeo)
library(aricode)
library(PRROC)


set.seed(123) # for reproducibility

#######################
# 2. HELPER FUNCTIONS
#######################


# Fits the occuN model and calculates metrics.
fit_and_evaluate_occuN <- function(y_data, obs_covs, w_matrix, landscape_raster, val_data, n_restarts = 30, optimizer_method = "BFGS") {

    # a. Create the unmarkedFrame
    umf <- unmarked::unmarkedFrameOccuN(y = y_data, cellCovs = as.data.frame(landscape_raster),
                                obsCovs = obs_covs, w = w_matrix)

    

    # Multiple restarts logic
    best_fit <- NULL
    min_nll <- Inf

    num_params <- 4 

    cat("  Starting optimization with", n_restarts, "restarts...\n")

    for (i in 1:n_restarts) {
        # random_starts <- runif(num_params, min = -5, max = 5)
        # random_startsstarts <- rnorm(num_params, mean = 0, sd = 0.1)
        random_starts <- runif(num_params, min = -0.1, max = 0.1)

        current_fit <- try(
            unmarked::occuN(~obs_cov ~state_cov, data = umf, starts = random_starts, se = FALSE, method = optimizer_method),
            silent = TRUE
        )

        if (!inherits(current_fit, "try-error")) {
            current_nll <- current_fit@negLogLike
            if (is.finite(current_nll) && current_nll < min_nll) {
                min_nll <- current_nll
                best_fit <- current_fit
                cat("    New best fit found at restart", i, "with NLL:", format(min_nll, digits = 15), "\n")
            }
        }
    }

    if (is.null(best_fit)) {
        warning("All model fitting attempts failed or resulted in non-finite NLL.")
        return(list(fit = NULL, metrics = data.frame(beta_intercept=NA, beta_1=NA, alpha_intercept=NA, alpha_1=NA,
                                                     lambda_pt_mse=NA, psi_pt_mse=NA, lambda_ls_mse=NA, psi_ls_mse=NA, auroc=NA),
                    pred_lambda_map = setValues(landscape_raster, NA),
                    pred_psi_map = setValues(landscape_raster, NA)))
    }

    final_fit <- unmarked::occuN(~obs_cov ~state_cov, data = umf, starts = coef(best_fit), se = TRUE, method = optimizer_method)

    beta_est <- coef(final_fit, type = 'state')
    alpha_est <- coef(final_fit, type = 'det')

    # Predictions on validation points
    pred_lambda_val <- exp(val_data$X_val %*% beta_est)
    pred_psi_val <- 1 - exp(-pred_lambda_val)
    
    val_obs_cov_df <- data.frame(obs_cov = val_data$obs_cov)
    X_p_val <- model.matrix(~obs_cov, data=val_obs_cov_df)
    
    pred_p_val_logit <- (X_p_val %*% alpha_est)
    pred_p_val <- plogis(as.vector(pred_p_val_logit))

    # Marginal probability of detection at validation points (for AUC)
    pred_marginal_det <- as.vector(pred_psi_val * pred_p_val)
    
    # Calculate AUROC using PRROC
    scores.class1 <- pred_marginal_det[val_data$y_val_obs == 1]
    scores.class0 <- pred_marginal_det[val_data$y_val_obs == 0]
    
    if(length(scores.class1) > 0 && length(scores.class0) > 0){
        roc_result <- PRROC::roc.curve(scores.class1 = scores.class1, scores.class0 = scores.class0, curve = FALSE)
        auroc_val <- roc_result$auc
    } else {
        auroc_val <- NA
    }
    
    # Get X matrix for all landscape cells
    cell_covs_df_for_pred <- as.data.frame(landscape_raster, cells = TRUE)
    X_pred_all_cells <- model.matrix(~ state_cov, data = cell_covs_df_for_pred)

    # Predictions on the entire landscape
    pred_lambda_landscape <- exp(X_pred_all_cells %*% beta_est)
    pred_psi_landscape <- 1 - exp(-pred_lambda_landscape)

    return(list(
        fit = final_fit,
        metrics = data.frame(
            beta_intercept = beta_est[1], beta_1 = beta_est[2],
            alpha_intercept = alpha_est[1], alpha_1 = alpha_est[2],
            lambda_pt_mse = mean((val_data$true_lambda_val - pred_lambda_val)^2),
            psi_pt_mse = mean((val_data$true_psi_val - pred_psi_val)^2),
            lambda_ls_mse = mean((val_data$true_lambda_j - pred_lambda_landscape)^2),
            psi_ls_mse = mean((val_data$true_psi_j - pred_psi_landscape)^2),
            auroc = auroc_val
        ),
        pred_lambda_map = setValues(landscape_raster, pred_lambda_landscape),
        pred_psi_map = setValues(landscape_raster, pred_psi_landscape)
    ))
}

# Helper function to generate cellSq clusterings.
create_cellSq_clusters <- function(points_df, landscape_raster, cell_side_length) {
    ext <- ext(landscape_raster)
    res_x <- xres(landscape_raster)
    res_y <- yres(landscape_raster)
    n_cols_landscape <- ncol(landscape_raster)
    grid_col <- floor((points_df$x - ext$xmin) / (res_x * cell_side_length))
    grid_row <- floor((ext$ymax - points_df$y) / (res_y * cell_side_length))
    n_cols_new_grid <- ceiling(n_cols_landscape / cell_side_length)
    site_ids <- grid_row * n_cols_new_grid + grid_col + 1
    return(as.numeric(factor(site_ids)))
}

###############################################################
# 3. SIMULATE LANDSCAPE AND SHARED DATA
###############################################################

n_cells_side = 50
landscape <- rast(nrows = n_cells_side, ncols = n_cells_side, xmin = 0, xmax = n_cells_side, ymin = 0, ymax = n_cells_side, crs="EPSG:32610")
n_cells <- ncell(landscape)

# Create a "blob-like" covariate by smoothing random noise
state_cov_noise <- rast(landscape)
values(state_cov_noise) <- runif(n_cells)
# Use a 3-cell circular moving window to smooth the noise into blobs
focal_window <- focalMat(state_cov_noise, 3, "circle") 
state_cov_blobby <- terra::focal(state_cov_noise, w = focal_window, fun = "mean", na.rm = TRUE)

# Rescale and name the final covariate
state_cov <- rast(landscape)
values(state_cov) <- scales::rescale(values(state_cov_blobby))
names(state_cov) <- "state_cov"

# Prepare shared data (constant across species)
cell_covs_df <- as.data.frame(state_cov, cells = TRUE)
X_truth_all_cells <- model.matrix(~ state_cov, data = cell_covs_df) # Used for both species


# Weighted sampling for initial points (to create co-located points) ---
n_train_points <- 2000
n_val_points <- 1000
n_points_total <- n_train_points + n_val_points
n_unique_locations <- 1500 # Define # of unique locations to sample from

cat("Sampling", n_unique_locations, "unique locations with weights based on state_cov...\n")

# Get covariate values to use as weights
# prob_weights <- values(state_cov)
# prob_weights[is.na(prob_weights)] <- 0 
# prob_weights <- prob_weights^1 # Power to increase bias
# prob_weights <- prob_weights + 1e-6 

# Sample cell indices based on weights for UNIQUE locations
# sampled_cell_indices <- sample(1:ncell(state_cov), 
#                                size = n_unique_locations, 
#                                replace = TRUE, # Allow multiple unique locations per cell
#                                prob = prob_weights)

sampled_cell_indices <- sample(1:ncell(state_cov), 
                               size = n_unique_locations, 
                               replace = TRUE) # Allow multiple unique locations per cell
                               


# Get center coordinates of sampled cells
cell_coords <- xyFromCell(state_cov, sampled_cell_indices)
res_x <- xres(state_cov)
res_y <- yres(state_cov)

# Jitter points randomly within their sampled cell
rand_x <- runif(n_unique_locations, min = cell_coords[, "x"] - res_x/2, max = cell_coords[, "x"] + res_x/2)
rand_y <- runif(n_unique_locations, min = cell_coords[, "y"] - res_y/2, max = cell_coords[, "y"] + res_y/2)

unique_locations_df <- data.frame(x = rand_x, y = rand_y)

# Create total points by sampling UNIQUE locations WITH REPLACEMENT ---
cat("Generating", n_points_total, "total points by sampling unique locations with replacement...\n")
sampled_rows_idx <- sample(1:n_unique_locations, n_points_total, replace = TRUE)
points_df <- unique_locations_df[sampled_rows_idx, ]
points_df$id <- 1:n_points_total # Add ID
cat("Done.\n")

# --- CHANGE 1: Calculate n_unique_pts ONCE ---
n_unique_pts <- nrow(unique(points_df[, c("x", "y")]))

# Split into train/validation (constant across experiments)
train_pts_df_base <- points_df %>% sample_n(n_train_points) # Base 2000 training points
val_pts_df <- points_df %>% anti_join(train_pts_df_base, by = "id") # Remaining 1000 validation points

# Prepare validation data covariates (constant across experiments)
val_pts_covs <- terra::extract(state_cov, val_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
X_val <- model.matrix(~state_cov, data=val_pts_covs)
val_cell_indices <- cellFromXY(state_cov, val_pts_df[,c("x","y")])
val_pts_df$obs_cov <- runif(nrow(val_pts_df)) # Constant obs_cov for validation

###############################################################
# 4. DEFINE SPECIES SCENARIOS
###############################################################

species_scenarios <- list(
    
    "Species_A_Generalist" = list(
        # Generalist: No spatial variation
        # --- CHANGE 2: BETA PARAMETERS SCALED FOR AREA ---
        # Original: beta_true = c("(Intercept)" = -1.0, "state_cov" = 1.0),
        beta_true = c("(Intercept)" = -3.7, "state_cov" = 1.0),
        alpha_true = c("(Intercept)" = 2.0, "obs_cov" = -1.0),
        
        # The "true" or "reference" clustering for a generalist
        # is grids (cellSq)
        ref_clustering_name = "5-cellSq"
    ),
    
    "Species_B_Specialist" = list(
        # Specialist: Strong, patchy spatial variation
        # --- CHANGE 3: BETA PARAMETERS SCALED FOR AREA ---
        # Original: beta_true = c("(Intercept)" = -3.0, "state_cov" = 5.0),
        beta_true = c("(Intercept)" = -3.9, "state_cov" = 5.0),
        alpha_true = c("(Intercept)" = 1.5, "obs_cov" = -1.0),
        
        # The "true" or "reference" clustering for a specialist
        # follows the patches (clustGeo)
        ref_clustering_name = "clustGeo-50-5"
    )
)

###############################################################
# 5. MAIN LOOP: Iterate through species
###############################################################

output_dir <- "output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# optimizers = "BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "SANN", "nlminb" 

selected_optimizer <- "Nelder-Mead"
# --- REMOVED: training_sizes loop and training_sets_list ---

# Initialize lists to store summary statistics for efficient CSV export
descriptive_stats_list <- list()
performance_metrics_list <- list()
estimated_params_list <- list()
perf_counter <- 1

# OUTER LOOP FOR SPECIES
for (species_name in names(species_scenarios)) {
    
    cat(paste("\n=============================================\n"))
    cat(paste("--- STARTING SIMULATION FOR:", species_name, "---\n"))
    cat(paste("=============================================\n"))
    
    # a. Get species-specific parameters
    current_species_params <- species_scenarios[[species_name]]
    beta_true <- current_species_params$beta_true
    alpha_true <- current_species_params$alpha_true
    ref_scenario_name <- current_species_params$ref_clustering_name # This is the "true" clustering

    # b. Generate SPECIES-SPECIFIC ground truth
    true_lambda_j <- as.vector(exp(X_truth_all_cells %*% beta_true))
    true_psi_j <- 1 - exp(-true_lambda_j)
    true_N_j <- rpois(n_cells, true_lambda_j) # Realized abundance

    true_lambda_map <- setValues(landscape, true_lambda_j)
    true_psi_map <- setValues(landscape, true_psi_j)
    true_N_map <- setValues(landscape, true_N_j)

    # Create sf object for all individuals for spatial intersection
    if(any(true_N_j > 0)) {
        cell_coords <- xyFromCell(landscape, 1:ncell(landscape))
        res_x <- xres(landscape); res_y <- yres(landscape)
        individual_points_list <- list()
        cells_with_individuals <- which(true_N_j > 0)
        for(cell_idx in cells_with_individuals) {
            n_inds <- true_N_j[cell_idx]
            center_x <- cell_coords[cell_idx, "x"]
            center_y <- cell_coords[cell_idx, "y"]
            rand_x <- runif(n_inds, min = center_x - res_x/2, max = center_x + res_x/2)
            rand_y <- runif(n_inds, min = center_y - res_y/2, max = center_y + res_y/2)
            individual_points_list[[length(individual_points_list) + 1]] <- data.frame(x = rand_x, y = rand_y)
        }
        all_individuals_df <- do.call(rbind, individual_points_list)
        individuals_sf <- st_as_sf(all_individuals_df, coords=c("x","y"), crs=st_crs(landscape))
    } else {
        all_individuals_df <- data.frame(x=numeric(0), y=numeric(0))
        individuals_sf <- st_as_sf(data.frame(x=numeric(0), y=numeric(0)), coords=c("x","y"), crs=st_crs(landscape))
    }

    # c. Generate SPECIES-SPECIFIC validation data (detections)
    true_lambda_val <- true_lambda_j[val_cell_indices]
    true_psi_val <- true_psi_j[val_cell_indices]
    
    p_det_val_logit <- alpha_true[1] + alpha_true[2] * val_pts_df$obs_cov
    p_det_val <- plogis(p_det_val_logit)
    
    true_z_val <- true_N_j[val_cell_indices] > 0
    y_val_obs <- rbinom(length(true_z_val), 1, p_det_val * true_z_val)
    
    # This list is SPECIES-SPECIFIC and passed to the fit function
    val_data_list <- list(
        y_val_obs = y_val_obs,
        X_val = X_val,
        obs_cov = val_pts_df$obs_cov,
        true_lambda_val = true_lambda_val,
        true_psi_val = true_psi_val,
        true_lambda_j = true_lambda_j, # For landscape-level MSE
        true_psi_j = true_psi_j      # For landscape-level MSE
    )

    # --- REMOVED: INNER LOOP FOR TRAINING SIZE ---
    
    cat(paste("\n--- Running simulation for Species:", species_name, "---\n"))
        
    # Get the pre-sampled training set
    train_pts_df <- train_pts_df_base # --- MODIFIED: Use base training set directly

    # d. Generate SPECIES-SPECIFIC reference clustering & training data
    
    cat(paste("\n  --- REFERENCE (TRUTH):", ref_scenario_name, "---\n"))
    
    # Create the "true" reference clustering for this species
    if (startsWith(ref_scenario_name, "clustGeo")) {
        train_pts_covs <- terra::extract(state_cov, train_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
        D0 <- dist(train_pts_covs)
        D1 <- dist(train_pts_df[, c("x", "y")])
        D0 <- (D0 - min(D0)) / (max(D0) - min(D0))
        D1 <- (D1 - min(D1)) / (max(D1) - min(D1))
        ref_tree <- hclustgeo(D0, D1, alpha = 0.5)

        ref_k_proportion <- as.numeric(sub(".*-", "", ref_scenario_name)) / 100
        # --- CHANGE 4: Use n_unique_pts for k calculation ---
        ref_k <- floor(ref_k_proportion * n_unique_pts) 
        cat(paste("Using k =", ref_k, "for clustGeo reference clustering.\n"))
        ref_clusters <- cutree(ref_tree, k = ref_k)
    } else {
        cell_side <- as.numeric(gsub("-cellSq", "", ref_scenario_name))
        cat(paste("Using", cell_side, "x", cell_side, "grid for reference clustering.\n"))
        ref_clusters <- create_cellSq_clusters(train_pts_df, landscape, cell_side)
    }

    train_pts_df$ref_site_id <- ref_clusters

    # --- START: CHANGE 5: MODIFIED POLYGON CREATION ---
    # Conditionally create polygons based on ref_scenario_name
    if (startsWith(ref_scenario_name, "clustGeo")) {
        # Determine the true occupancy of reference sites using CONVEX HULLS
        ref_polygons_sf <- train_pts_df %>%
            st_as_sf(coords = c("x", "y"), crs = st_crs(landscape)) %>%
            group_by(ref_site_id) %>%
            summarise(geometry = st_combine(geometry), .groups = 'drop') %>%  # 1. Combine points
            mutate(geometry = st_convex_hull(geometry)) %>%                  # 2. Get hull (might be POINT, LINE, or POLYGON)
            mutate(geometry = st_buffer(geometry, dist = 0.5))               # 3. Buffer to *guarantee* it's a POLYGON
    
    } else {
        # Use the *actual* grid squares for reference polygons
        cell_side <- as.numeric(gsub("-cellSq", "", ref_scenario_name))
        ext <- ext(landscape); res_x <- xres(landscape); res_y <- yres(landscape)
        site_representatives <- train_pts_df %>% group_by(ref_site_id) %>% slice(1) %>% ungroup()
        
        ref_polygons_sf <- site_representatives %>%
            mutate(
                grid_col = floor((x - ext$xmin) / (res_x * cell_side)),
                grid_row = floor((ext$ymax - y) / (res_y * cell_side)),
                xmin_sq = ext$xmin + grid_col * res_x * cell_side,
                xmax_sq = xmin_sq + res_x * cell_side,
                ymax_sq = ext$ymax - grid_row * res_y * cell_side,
                ymin_sq = ymax_sq - res_y * cell_side
            ) %>%
            rowwise() %>%
            mutate(geometry = list(st_polygon(list(matrix(
                    c(xmin_sq, ymin_sq, xmax_sq, ymin_sq, xmax_sq, ymax_sq, xmin_sq, ymax_sq, xmin_sq, ymin_sq),
                    ncol = 2, byrow = TRUE
            ))))) %>%
            ungroup() %>% st_as_sf(crs = st_crs(landscape)) %>% select(ref_site_id, geometry)
    }
    
    # Check for empty/invalid geometries and remove them
    ref_polygons_sf <- ref_polygons_sf[!st_is_empty(ref_polygons_sf), ]
    ref_polygons_sf <- ref_polygons_sf[st_is_valid(ref_polygons_sf), ]
    # --- END: MODIFIED POLYGON CREATION ---
    
    ref_site_true_occupancy <- apply(st_intersects(ref_polygons_sf, individuals_sf, sparse = FALSE), 1, any)

    # Calculate true site-level lambda and psi by extracting weighted means from true maps
    true_site_lambda <- terra::extract(true_lambda_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]
    true_site_psi <- terra::extract(true_psi_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]

    # Generate SPECIES-SPECIFIC training detections
    site_lambda_map <- data.frame(ref_site_id = ref_polygons_sf$ref_site_id, true_site_lambda = true_site_lambda)
    
    # Ensure ref_site_id is present for the join
    ref_polygons_sf <- ref_polygons_sf %>% mutate(ref_site_id = seq_len(nrow(.)))
    site_lambda_map <- data.frame(ref_site_id = ref_polygons_sf$ref_site_id, true_site_lambda = true_site_lambda)
    
    train_pts_df <- train_pts_df %>%
        left_join(site_lambda_map, by = "ref_site_id")
        
    train_pts_df <- train_pts_df %>%
        mutate(
            site_is_occupied = ref_site_true_occupancy[ref_site_id],
            obs_cov = runif(n())
        )
        
    p_detection_logit <- alpha_true[1] + alpha_true[2] * train_pts_df$obs_cov
    p_detection <- plogis(p_detection_logit)
        
    train_pts_df$detection <- rbinom(nrow(train_pts_df), 1, p_detection)
    train_pts_df$detection <- train_pts_df$detection * train_pts_df$site_is_occupied
        
    train_pts_df <- train_pts_df %>% select(-true_site_lambda)

    points_per_cluster <- train_pts_df %>% count(ref_site_id) %>% pull(n)

    # Add descriptive stats for this Species
    new_desc_row <- data.frame(
        "Species" = species_name,
        "Reference.Clustering" = ref_scenario_name,
        "Number.of.Clusters" = length(unique(train_pts_df$ref_site_id)),
        "Mean.Number.of.Points.in.Cluster" = round(mean(points_per_cluster), 4),
        "SD.Number.of.points.in.Clusters" = round(sd(points_per_cluster), 4),
        "Site.Lambda.Mean" = round(mean(true_site_lambda, na.rm = TRUE), 4),
        "Site.Lambda.SD" = round(sd(true_site_lambda, na.rm = TRUE), 4),
        "Site.Psi.Mean" = round(mean(true_site_psi, na.rm = TRUE), 4),
        "Site.Psi.SD" = round(sd(true_site_psi, na.rm = TRUE), 4),
        "Prevalence" = round(sum(train_pts_df$detection) / nrow(train_pts_df), 4)
    )

    descriptive_stats_list[[length(descriptive_stats_list) + 1]] <- new_desc_row

    ########################################################################
    # 6. ANALYSIS CLUSTERING
    #    Apply all 5 analysis methods to this species' data
    ########################################################################

    # a. Define the different clustGeo scenarios to test
    #    (We must re-calculate the tree for each training set)
    train_pts_covs <- terra::extract(state_cov, train_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
    D0 <- dist(train_pts_covs); D1 <- dist(train_pts_df[, c("x", "y")])
    D0 <- (D0 - min(D0)) / (max(D0) - min(D0)); D1 <- (D1 - min(D1)) / (max(D1) - min(D1))
    analysis_tree <- hclustgeo(D0, D1, alpha = 0.5) # Use this tree for all analysis clusters
    
    # NOTE: n_unique_pts is now defined globally (Change 1)
    k_values <- c("clustGeo-50-25" = floor(0.25 * n_unique_pts), 
                  "clustGeo-50-90" = floor(0.90 * n_unique_pts), 
                  "clustGeo-50-5" = floor(0.05 * n_unique_pts))
    k_values <- k_values[c("clustGeo-50-25", "clustGeo-50-90", "clustGeo-50-5")] # Re-order
    
    cluster_scenarios <- lapply(k_values, function(k) cutree(analysis_tree, k))
    
    # b. Add cellSq clustering scenarios
    cellSq_1_clusters <- create_cellSq_clusters(train_pts_df, landscape, 1)
    cellSq_5_clusters <- create_cellSq_clusters(train_pts_df, landscape, 5)

    # c. Combine all 5 ANALYSIS scenarios into a single named list
    all_scenarios <- c(cluster_scenarios, 
                       "1-cellSq" = list(cellSq_1_clusters), 
                       "5-cellSq" = list(cellSq_5_clusters))

    # d. Run the pipeline for each analysis scenario
    results <- list()

    for (current_scenario_name in names(all_scenarios)) {
        cat(paste("    - Analysis Method:", current_scenario_name, "\n"))
        cl <- all_scenarios[[current_scenario_name]]
        train_pts_df$current_site_id <- cl

        site_data <- train_pts_df %>%
            group_by(current_site_id) %>%
            summarise(detections = list(detection), obs_covs_list = list(obs_cov), n_pts = n(), .groups = 'drop')

        max_surveys <- max(site_data$n_pts)
        n_sites <- nrow(site_data)
        y_data <- matrix(NA, nrow = n_sites, ncol = max_surveys)
        obs_cov_matrix <- matrix(NA, nrow = n_sites, ncol = max_surveys)
        for(i in 1:n_sites){
            y_data[i, 1:site_data$n_pts[i]] <- unlist(site_data$detections[i])
            obs_cov_matrix[i, 1:site_data$n_pts[i]] <- unlist(site_data$obs_covs_list[i])
        }
        obs_covs_for_umf <- list(obs_cov = obs_cov_matrix)

        if (grepl("-cellSq", current_scenario_name, fixed = TRUE)) {
            cell_side <- as.numeric(gsub("-cellSq", "", current_scenario_name))
            ext <- ext(landscape); res_x <- xres(landscape); res_y <- yres(landscape)
            site_representatives <- train_pts_df %>% group_by(current_site_id) %>% slice(1) %>% ungroup()
            polygons_df <- site_representatives %>%
                mutate(
                    grid_col = floor((x - ext$xmin) / (res_x * cell_side)),
                    grid_row = floor((ext$ymax - y) / (res_y * cell_side)),
                    xmin_sq = ext$xmin + grid_col * res_x * cell_side,
                    xmax_sq = xmin_sq + res_x * cell_side,
                    ymax_sq = ext$ymax - grid_row * res_y * cell_side,
                    ymin_sq = ymax_sq - res_y * cell_side
                ) %>%
                rowwise() %>%
                mutate(geometry = list(st_polygon(list(matrix(
                        c(xmin_sq, ymin_sq, xmax_sq, ymin_sq, xmax_sq, ymax_sq, xmin_sq, ymax_sq, xmin_sq, ymin_sq),
                        ncol = 2, byrow = TRUE
                ))))) %>%
                ungroup() %>% st_as_sf(crs = st_crs(landscape)) %>% select(current_site_id, geometry)
        } else {
            # --- Use the same robust logic for analysis polygons ---
            polygons_df <- train_pts_df %>%
                st_as_sf(coords = c("x", "y"), crs = st_crs(landscape)) %>%
                group_by(current_site_id) %>%
                summarise(geometry = st_combine(geometry), .groups = 'drop') %>%
                mutate(geometry = st_convex_hull(geometry)) %>%
                mutate(geometry = st_buffer(geometry, dist = 0.5))
        }
        
        # Remove empty/invalid geometries
        polygons_df <- polygons_df[!st_is_empty(polygons_df), ]
        polygons_df <- polygons_df[st_is_valid(polygons_df), ]
        
        # Ensure site_data and polygons_df align
        site_data <- site_data %>% filter(current_site_id %in% polygons_df$current_site_id)
        n_sites <- nrow(site_data)
        
        # Re-filter y_data and obs_covs_for_umf if sites were dropped
        y_data <- y_data[site_data$current_site_id %in% polygons_df$current_site_id, , drop = FALSE]
        obs_cov_matrix <- obs_cov_matrix[site_data$current_site_id %in% polygons_df$current_site_id, , drop = FALSE]
        obs_covs_for_umf <- list(obs_cov = obs_cov_matrix)
        
        
        weights_extract <- terra::extract(state_cov, vect(polygons_df), weights = TRUE, cells = TRUE)
        w_matrix <- matrix(0, nrow = n_sites, ncol = ncell(landscape))
        
        if(nrow(weights_extract) > 0) {
          for(i in seq_len(nrow(weights_extract))){
              site_idx <- weights_extract[i, "ID"]
              cell_idx <- weights_extract[i, "cell"]
              actual_site_id <- polygons_df$current_site_id[site_idx]
              w_matrix_row_idx <- which(site_data$current_site_id == actual_site_id)
              w_matrix[w_matrix_row_idx, cell_idx] <- weights_extract[i, "weight"]
          }
        }
        
        # Pass the SPECIES-SPECIFIC validation data
        analysis_results <- fit_and_evaluate_occuN(y_data, obs_covs_for_umf, w_matrix, state_cov, val_data=val_data_list, optimizer_method = selected_optimizer)

        # Calculate ARI against the "true" clustering for this species
        current_ari <- ARI(ref_clusters, cl)
        analysis_results$ari <- current_ari

        # Store estimated parameters for CSV export
        new_params_row <- data.frame(
            "Species" = species_name, # --- NEW COLUMN ---
            "Reference.Clustering" = ref_scenario_name,
            # "Training.Size" = n_train_points, # --- REMOVED ---
            "Clustering" = current_scenario_name,
            "beta_intercept" = analysis_results$metrics$beta_intercept,
            "beta_1" = analysis_results$metrics$beta_1,
            "alpha_intercept" = analysis_results$metrics$alpha_intercept,
            "alpha_1" = analysis_results$metrics$alpha_1
        )
        estimated_params_list[[length(estimated_params_list) + 1]] <- new_params_row

        # Calculate and store performance metrics
        pred_cell_stats <- global(c(analysis_results$pred_lambda_map, analysis_results$pred_psi_map), c("mean", "sd"), na.rm=TRUE)
        
        # Use the "true" polygons for this species to extract site stats
        pred_site_lambda <- terra::extract(analysis_results$pred_lambda_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]
        pred_site_psi <- terra::extract(analysis_results$pred_psi_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]

        new_perf_row <- data.frame(
            "No." = perf_counter,
            "Species" = species_name,
            "Reference.Clustering" = ref_scenario_name,
            "Clustering" = current_scenario_name,
            "Optimizer" = selected_optimizer,
            "Cell.Lambda.Mean" = round(pred_cell_stats[1, "mean"], 4),
            "Cell.Lambda.SD" = round(pred_cell_stats[1, "sd"], 4),
            "Cell.Psi.Mean" = round(pred_cell_stats[2, "mean"], 4),
            "Cell.Psi.SD" = round(pred_cell_stats[2, "sd"], 4),
            "Site.Lambda.Mean" = round(mean(pred_site_lambda, na.rm = TRUE), 4),
            "Site.Lambda.SD" = round(sd(pred_site_lambda, na.rm = TRUE), 4),
            "Site.Psi.Mean" = round(mean(pred_site_psi, na.rm = TRUE), 4),
            "Site.Psi.SD" = round(sd(pred_site_psi, na.rm = TRUE), 4),
            "Lambda.LS.MSE" = analysis_results$metrics$lambda_ls_mse,
            "Psi.LS.MSE" = analysis_results$metrics$psi_ls_mse,
            "ARI" = round(current_ari, 4),
            "AUROC" = round(analysis_results$metrics$auroc, 4)
        )
        performance_metrics_list[[length(performance_metrics_list) + 1]] <- new_perf_row
        perf_counter <- perf_counter + 1

        results[[current_scenario_name]] <- c(list(polygons = polygons_df), analysis_results)
    }

    #################################
    # 7. PLOT GENERATION (7x4 GRID)
    #################################

    # Re-order plots so the "true" analysis method is first, if it exists
    if(ref_scenario_name %in% names(results)) {
        other_scenarios <- setdiff(names(results), ref_scenario_name)
        plot_order <- c(ref_scenario_name, other_scenarios)
    } else {
        # This happens if, e.g., ref is clustGeo-50-25 but analysis is only 10, 90
        plot_order <- names(results) 
    }
    # Ensure we only have 5
    plot_order <- plot_order[1:5]
    results <- results[plot_order]

    # Save PNG to output folder
    # Filename includes SPECIES and ref (training size removed)
    png_filename <- file.path(output_dir, paste0("species=", species_name, "_ref=", ref_scenario_name, ".png"))
    png(png_filename, width = 14, height = 20, units = "in", res = 300)
    cat(paste("Generating plot:", png_filename, "\n"))

    # Calculate the maximum abundance from ONLY the reference map
    ref_max_val <- global(true_lambda_map, "max", na.rm = TRUE)$max
    lambda_max <- ref_max_val * 1.00
    if (!is.finite(lambda_max) || lambda_max == 0) {
        warning("Could not determine a finite reference max; falling back to 10.")
        lambda_max <- 10 
    }

    # Robust legend axis logic
    abun_breaks_internal <- pretty(c(0, lambda_max * 0.9))
    abun_breaks_internal <- abun_breaks_internal[abun_breaks_internal >= 0 & abun_breaks_internal < lambda_max]
    abun_breaks <- c(abun_breaks_internal, lambda_max)
    abun_labels <- as.character(abun_breaks_internal)
    lambda_max_rounded <- round(lambda_max, 1)
    abun_labels <- c(abun_labels, paste0(">", lambda_max_rounded))
    abun_axis_args <- list(at = abun_breaks, labels = abun_labels)

    # Define color palettes
    covariate_col_palette <- viridis(100, option = "turbo")
    abundance_col_palette <- viridis(100, option = "viridis")
    occupancy_col_palette <- viridis(100, option = "magma")
   
    par(cex = 2.0)
    layout(matrix(1:28, nrow = 7, byrow = TRUE), widths = c(1, 1, 1, 0.8), heights = c(1,1, 1, 1, 1, 1,1))
    margin_default <- c(0,0,0,0)
    text_offset = -0.75

    # ROW 1: Ground Truth Point Pattern, Abundance, and Occupancy
    par(mar = margin_default)
    plot(state_cov, main = "Point Pattern", legend = FALSE, col = "white")
    if(nrow(all_individuals_df) > 0) {
      points(all_individuals_df$x, all_individuals_df$y, pch = 16, col = "red", cex = 0.5)
    }

    par(mar = margin_default)
    plot(true_N_map, main = "Abundance Pattern", col = viridis(100, option="cividis"), axes = TRUE)
    text(x = xyFromCell(true_N_map, 1:ncell(true_N_map))[,1],
        y = xyFromCell(true_N_map, 1:ncell(true_N_map))[,2],
        labels = values(true_N_map),
        cex = 0.6, col = "white")

    true_occ_map <- ifel(true_N_map > 0, 1, 0)
    plot(true_occ_map, main = "Occurrence Pattern", col = c("black", "white"), legend = FALSE)

    par(mar=c(0,0,0,0))
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
    text(text_offset, 0.7, "Reference Species", pos=4, font=2, cex=1.2)
    text(text_offset, 0.55, species_name, pos=4, cex=1.2, col="blue") 
    text(text_offset, 0.35, paste("Total Individuals:", sum(true_N_j)), pos=4, cex=1.2)
    text(text_offset, 0.2, paste("Realized Mean Abundance:", round(mean(true_N_j), 2)), pos=4, cex=1.2)
    text(text_offset, 0.05, paste("Realized Occupancy:", round(mean(true_N_j > 0), 2)), pos=4, cex=1.2)

    # ROW 2: Reference
    par(mar = margin_default)
    plot(state_cov, main = "Covariate & Observations", col = covariate_col_palette, axes = TRUE)
    points(train_pts_df$x, train_pts_df$y, pch = 16, col = "blue", cex = 0.5)
    points(val_pts_df$x, val_pts_df$y, pch = 17, col = "magenta1", cex = 0.5)
    legend("bottom", legend = c("Train", "Val"), col = c("blue", "magenta1"), pch = c(16, 17), horiz = TRUE, bg="white", cex=1.0)

    par(mar = margin_default)
    plot(true_lambda_map, main = "Reference Lambda", col = abundance_col_palette, range = c(0, lambda_max), extend = "max", pax = abun_axis_args, axes = TRUE)
    
    plot(true_psi_map, main = "Reference Psi", col = occupancy_col_palette, range=c(0,1), axes = TRUE)

    # --- MODIFICATION 1: ADD TRUE PARAMETERS TO ROW 2, COL 4 ---
    par(mar=c(0,0,0,0))
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
    
    beta_0_true <- round(beta_true[1], 2)
    beta_1_true <- round(beta_true[2], 2)
    alpha_0_true <- round(alpha_true[1], 2)
    alpha_1_true <- round(alpha_true[2], 2)
    
    b0_text <- bquote(paste(beta[0], ": ", .(beta_0_true)))
    b1_text <- bquote(paste(beta[1], ": ", .(beta_1_true)))
    a0_text <- bquote(paste(alpha[0], ": ", .(alpha_0_true)))
    a1_text <- bquote(paste(alpha[1], ": ", .(alpha_1_true)))
    
    text(text_offset, 0.9, "Simulation Stats", pos=4, font=2, cex=1.2)
    text(text_offset, 0.7, "True Parameters", pos=4, font=2, cex=1.1)
    text(text_offset, 0.55, b0_text, pos=4, cex=1.1)
    text(text_offset, 0.4, b1_text, pos=4, cex=1.1)
    text(text_offset, 0.25, a0_text, pos=4, cex=1.1)
    text(text_offset, 0.1, a1_text, pos=4, cex=1.1)
    
    # nrow(train_pts_df) will correctly show 2000
    text(text_offset, -0.1, paste("Training Pts.:", nrow(train_pts_df)), pos=4, cex=1.1) 
    text(text_offset, -0.25, paste("Validation Pts.:", nrow(val_pts_df)), pos=4, cex=1.1)
    text(text_offset, -0.4, paste("Cell Mean Exp. Abun.:", round(mean(.GlobalEnv$true_lambda_j), 2)), pos=4, cex=1.1)
    text(text_offset, -0.55, paste("Cell Mean Occu. Prob.:", round(mean(.GlobalEnv$true_psi_j), 2)), pos=4, cex=1.1)
    # --- END MODIFICATION 1 ---


    # ROWS 3-7: Clustering Scenarios
    scenario_names <- names(results)
    for(i in 1:5){
        res <- results[[i]]
        metrics <- res$metrics
        current_scenario_name <- scenario_names[i]

        trailing_ref_str <- paste(current_scenario_name, "Sites")
        if (current_scenario_name == ref_scenario_name) {
            # Highlight the "true" analysis method
            trailing_ref_str <- paste(trailing_ref_str, "(SPECIES TRUTH)")
        }

        par(mar = margin_default)
        plot(state_cov, main = trailing_ref_str, legend = FALSE, col = "white")
        plot(res$polygons$geometry, add = TRUE, border = "gray30", col = alpha("gray", 0.3))
        points(train_pts_df$x, train_pts_df$y, pch = 16, col = "blue", cex = 0.5)

        if(is.null(res$fit) || all(is.na(values(res$pred_lambda_map)))) {
            par(mar = margin_default)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", main = expression(paste(lambda, "Estimated from", current_scenario_name)))
            text(0, 0, "Fit Failed", cex=1.5)
            par(mar = margin_default)
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", main = expression(paste(psi, "Estimated from", current_scenario_name)))
            text(0, 0, "Fit Failed", cex=1.5)
        } else {
            par(mar = margin_default)
            lambda_map_to_plot <- res$pred_lambda_map
            v_cleaned <- values(lambda_map_to_plot)
            v_cleaned[!is.finite(v_cleaned)] <- lambda_max
            v_cleaned[v_cleaned > lambda_max] <- lambda_max
            v_cleaned[v_cleaned < 0] <- 0
            values(lambda_map_to_plot) <- v_cleaned
            
            plot(lambda_map_to_plot, main = paste("Abundance Estimated from", current_scenario_name), 
                 col = abundance_col_palette, range = c(0, lambda_max), 
                 extend = "max", pax = abun_axis_args, axes = TRUE)
            
            plot(res$pred_psi_map, main = paste("Occupancy Estimated from", current_scenario_name), col = occupancy_col_palette, range=c(0,1), axes = TRUE)
        }

        par(mar=c(0,0,0,0))
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
        
        text(text_offset, 0.95, "Metrics", pos=4, font=2, cex=1.2)
        # ARI is vs. the "true" ref_clusters for this SPECIES
        text(text_offset, 0.8, paste("ARI (vs. truth):", round(res$ari, 4)), pos=4, cex=1.1) 
        
        b0_est_val <- if(is.na(metrics$beta_intercept)) "NA" else round(metrics$beta_intercept, 2)
        b1_est_val <- if(is.na(metrics$beta_1)) "NA" else round(metrics$beta_1, 2)
        a0_est_val <- if(is.na(metrics$alpha_intercept)) "NA" else round(metrics$alpha_intercept, 2)
        a1_est_val <- if(is.na(metrics$alpha_1)) "NA" else round(metrics$alpha_1, 2)
        
        b0_est_text <- bquote(paste(hat(beta)[0], ": ", .(b0_est_val)))
        b1_est_text <- bquote(paste(hat(beta)[1], ": ", .(b1_est_val)))
        a0_est_text <- bquote(paste(hat(alpha)[0], ": ", .(a0_est_val)))
        a1_est_text <- bquote(paste(hat(alpha)[1], ": ", .(a1_est_val)))
        
        text(text_offset, 0.65, b0_est_text, pos=4, cex=1.1)
        text(text_offset, 0.5, b1_est_text, pos=4, cex=1.1)
        text(text_offset, 0.35, a0_est_text, pos=4, cex=1.1)
        text(text_offset, 0.2, a1_est_text, pos=4, cex=1.1)
        
        lambda_ls_mse_text <- if (is.na(metrics$lambda_ls_mse)) "NA" else formatC(metrics$lambda_ls_mse, format = "e", digits = 2)
        psi_ls_mse_text <- if (is.na(metrics$psi_ls_mse)) "NA" else formatC(metrics$psi_ls_mse, format = "e", digits = 2)
        lambda_pt_mse_text <- if (is.na(metrics$lambda_pt_mse)) "NA" else formatC(metrics$lambda_pt_mse, format = "e", digits = 2)
        psi_pt_mse_text <- if (is.na(metrics$psi_pt_mse)) "NA" else formatC(metrics$psi_pt_mse, format = "e", digits = 2)
        
        text(text_offset, 0.0, paste("Cell Abun. MSE:", lambda_ls_mse_text), pos=4, cex=1.1)
        text(text_offset, -0.15, paste("Cell Occu. MSE:", psi_ls_mse_text), pos=4, cex=1.1)
        text(text_offset, -0.3, paste("Val. Pts. Abun. MSE:", lambda_pt_mse_text), pos=4, cex=1.1)
        text(text_offset, -0.45, paste("Val. Pts. Occ. MSE:", psi_pt_mse_text), pos=4, cex=1.1)
        
        auroc_text <- if (is.na(metrics$auroc)) "NA" else round(metrics$auroc, 4)
        text(text_offset, -0.6, paste("AUC:", auroc_text), pos=4, cex=1.1)
  
        par(mar = margin_default)
    }
    dev.off()
    
    
} 
#############################
# 8. WRITE SUMMARY CSV FILES
#############################

# Convert lists to data frames
descriptive_stats_df <- do.call(rbind, descriptive_stats_list)
performance_metrics_df <- do.call(rbind, performance_metrics_list)
estimated_params_df <- do.call(rbind, estimated_params_list)

# Save CSVs to output folder
write.csv(descriptive_stats_df, file.path(output_dir, "descriptive_stats.csv"), row.names = FALSE)
write.csv(performance_metrics_df, file.path(output_dir, "performance_metrics.csv"), row.names = FALSE)
write.csv(estimated_params_df, file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)


cat("\n--- All simulations complete. ---\n")







