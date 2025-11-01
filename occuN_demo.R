#########################################################
# Script for occuN Clustering Comparison Plot
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
# --- MODIFIED: Removed lambda_reg_alpha and lambda_reg_beta arguments ---
fit_and_evaluate_occuN <- function(y_data, obs_covs, w_matrix, landscape_raster, val_data, n_restarts = 30, optimizer_method = "BFGS") {

    # a. Create the unmarkedFrame
    umf <- unmarked::unmarkedFrameOccuN(y = y_data, cellCovs = as.data.frame(landscape_raster),
                                obsCovs = obs_covs, w = w_matrix)

    # Multiple restarts logic
    best_fit <- NULL
    min_nll <- Inf

    # --- MODIFIED: num_params is 4 (2 state, 2 det) ---
    num_params <- 4 

    cat("  Starting optimization with", n_restarts, "restarts...\n")

    for (i in 1:n_restarts) {
        random_starts <- runif(num_params, min = -5, max = 5)
        # --- MODIFIED: Removed regularization parameters from occuN call ---
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
        # --- MODIFIED: Removed alpha_lambda from failure dataframe ---
        return(list(fit = NULL, metrics = data.frame(beta_intercept=NA, beta_1=NA, alpha_intercept=NA, alpha_1=NA,
                                                     lambda_pt_mse=NA, psi_pt_mse=NA, lambda_ls_mse=NA, psi_ls_mse=NA, auroc=NA),
                    pred_lambda_map = setValues(landscape_raster, NA),
                    pred_psi_map = setValues(landscape_raster, NA)))
    }

    # --- MODIFIED: Removed regularization parameters from occuN call ---
    final_fit <- unmarked::occuN(~obs_cov ~state_cov, data = umf, starts = coef(best_fit), se = TRUE, method = optimizer_method)

    # --- UPDATED PARAMETER EXTRACTION ---
    beta_est <- coef(final_fit, type = 'state')
    alpha_est <- coef(final_fit, type = 'det')
    # --- MODIFIED: Removed alpha_lambda_est ---
    # --- END UPDATE ---

    # Predictions on validation points
    pred_lambda_val <- exp(val_data$X_val %*% beta_est)
    pred_psi_val <- 1 - exp(-pred_lambda_val)
    
    # --- UPDATED PREDICTION LOGIC FOR AUROC ---
    # Predicted detection probability for validation points
    val_obs_cov_df <- data.frame(obs_cov = val_data$obs_cov)
    X_p_val <- model.matrix(~obs_cov, data=val_obs_cov_df)

    # --- MODIFIED: Removed log-abundance from detection probability calculation ---
    # Get predicted log-abundance for validation points (NO LONGER NEEDED FOR P)
    # pred_log_lambda_val <- log(pred_lambda_val)
    # Handle -Inf values from log(0) (NO LONGER NEEDED)
    # pred_log_lambda_val[is.infinite(pred_log_lambda_val)] <- -100 
        
    # Calculate logit-p using the
    pred_p_val_logit <- (X_p_val %*% alpha_est)
    pred_p_val <- plogis(as.vector(pred_p_val_logit))
    # --- END UPDATE ---

    # Marginal probability of detection at validation points (for AUC)
    pred_marginal_det <- as.vector(pred_psi_val * pred_p_val)
    
    # Calculate AUROC using PRROC
    scores.class1 <- pred_marginal_det[val_data$y_val_obs == 1]
    scores.class0 <- pred_marginal_det[val_data$y_val_obs == 0]
    
    # Handle cases where one class is missing in the validation data
    if(length(scores.class1) > 0 && length(scores.class0) > 0){
        roc_result <- PRROC::roc.curve(scores.class1 = scores.class1, scores.class0 = scores.class0, curve = FALSE)
        auroc_val <- roc_result$auc
    } else {
        auroc_val <- NA
    }

    # Predictions on the entire landscape
    pred_lambda_landscape <- exp(X_truth_all_cells %*% beta_est)
    pred_psi_landscape <- 1 - exp(-pred_lambda_landscape)

    return(list(
        fit = final_fit,
        # --- UPDATED METRICS DATAFRAME ---
        # --- MODIFIED: Removed alpha_lambda ---
        metrics = data.frame(
            beta_intercept = beta_est[1], beta_1 = beta_est[2],
            alpha_intercept = alpha_est[1], alpha_1 = alpha_est[2],
            lambda_pt_mse = mean((true_lambda_val - pred_lambda_val)^2),
            psi_pt_mse = mean((true_psi_val - pred_psi_val)^2),
            lambda_ls_mse = mean((true_lambda_j - pred_lambda_landscape)^2),
            psi_ls_mse = mean((true_psi_j - pred_psi_landscape)^2),
            auroc = auroc_val
        ),
        # --- END UPDATE ---
        pred_lambda_map = setValues(landscape_raster, pred_lambda_landscape),
        pred_psi_map = setValues(landscape_raster, pred_psi_landscape)
    ))
}

# Helper function to generate cellSq clusterings.
create_cellSq_clusters <- function(points_df, landscape_raster, cell_side_length) {

    # Get landscape metadata
    ext <- ext(landscape_raster)
    res_x <- xres(landscape_raster)
    res_y <- yres(landscape_raster)
    n_cols_landscape <- ncol(landscape_raster)

    # Calculate the column and row index for each point in the new, larger grid
    grid_col <- floor((points_df$x - ext$xmin) / (res_x * cell_side_length))
    grid_row <- floor((ext$ymax - points_df$y) / (res_y * cell_side_length))

    # Calculate the number of columns in the new grid
    n_cols_new_grid <- ceiling(n_cols_landscape / cell_side_length)

    # Create a unique ID for each grid cell (site)
    site_ids <- grid_row * n_cols_new_grid + grid_col + 1

    # Return a dense, sequential vector of cluster IDs (like cutree does) to
    # ensure they can be used as row indices later.
    return(as.numeric(factor(site_ids)))
}

###############################################################
# 3. SIMULATE GROUND-TRUTH LANDSCAPE AND SPECIES DISTRIBUTION
###############################################################


n_cells_side = 50
landscape <- rast(nrows = n_cells_side, ncols = n_cells_side, xmin = 0, xmax = n_cells_side, ymin = 0, ymax = n_cells_side, crs="EPSG:32610")
n_cells <- ncell(landscape)

state_cov <- rast(landscape)
xy <- xyFromCell(state_cov, 1:ncell(state_cov))
values(state_cov) <- scales::rescale(xy[,1] + xy[,2] + rnorm(n_cells, 0, 10))
names(state_cov) <- "state_cov"

# True parameters
beta_true <- c("(Intercept)" = 1.0, "state_cov" = -3.0)
alpha_true <- c("(Intercept)" = 1.0, "obs_cov" = -1.0)
# --- MODIFIED: Removed true abundance effect parameter ---
# --- END ADDITION ---


# True landscape values
cell_covs_df <- as.data.frame(state_cov, cells = TRUE)
X_truth_all_cells <- model.matrix(~ state_cov, data = cell_covs_df)
true_lambda_j <- as.vector(exp(X_truth_all_cells %*% beta_true))
true_psi_j <- 1 - exp(-true_lambda_j)
true_N_j <- rpois(n_cells, true_lambda_j)

# Create rasters for true lambda and psi
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

        # Fetch no. of individuals for this cell
        n_inds <- true_N_j[cell_idx]

        # Calculate cell centers
        center_x <- cell_coords[cell_idx, "x"]
        center_y <- cell_coords[cell_idx, "y"]

        # Randomly distribute individuals in each cell based on N_j (n_inds).
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

# Simulate total points and split into train/validation (constant across experiments)
n_points_total <- 5000
points_df <- data.frame(id = 1:n_points_total, x = runif(n_points_total, 0, n_cells_side), y = runif(n_points_total, 0, n_cells_side))
n_train_total <- 4000
train_pts_df_base <- points_df %>% sample_n(n_train_total) # Base 4000 training points
val_pts_df <- points_df %>% anti_join(train_pts_df_base, by = "id") # Remaining 1000 validation points

# Prepare validation data (constant across experiments)
val_pts_covs <- terra::extract(state_cov, val_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
X_val <- model.matrix(~state_cov, data=val_pts_covs)
val_cell_indices <- cellFromXY(state_cov, val_pts_df[,c("x","y")])
true_lambda_val <- true_lambda_j[val_cell_indices]
true_psi_val <- true_psi_j[val_cell_indices]

# --- UPDATED SIMULATION OF VALIDATION DETECTIONS (y) ---
val_pts_df$obs_cov <- runif(nrow(val_pts_df))

# --- MODIFIED: Removed log-abundance from detection probability calculation ---
# Get true log-abundance for validation points (NO LONGER NEEDED)
# true_log_lambda_val <- log(true_lambda_val)
# Handle -Inf (NO LONGER NEEDED)
# true_log_lambda_val[is.infinite(true_log_lambda_val)] <- -100

# Calculate logit-p using the
p_det_val_logit <- alpha_true[1] + alpha_true[2] * val_pts_df$obs_cov
p_det_val <- plogis(p_det_val_logit)

true_z_val <- true_N_j[val_cell_indices] > 0
y_val_obs <- rbinom(length(true_z_val), 1, p_det_val * true_z_val)
# --- END UPDATE ---

val_data_list <- list(
    y_val_obs = y_val_obs,
    X_val = X_val,
    obs_cov = val_pts_df$obs_cov
)

###############################################################
# 4. MAIN LOOP: Iterate through each reference scenario
###############################################################


output_dir <- "output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}


# Select optimizer
# Possible values: "BFGS", "nlminb", "Nelder-Mead", "CG", "L-BFGS-B", "SANN"
selected_optimizer <- "Nelder-Mead"

# Define training set sizes for the learning curve
training_sizes <- c(500, 1000, 2000, 4000)

# Pre-sample all training sets to ensure consistency
cat("Pre-sampling training sets...\n")
training_sets_list <- list()
for (size in training_sizes) {
    training_sets_list[[as.character(size)]] <- train_pts_df_base %>% sample_n(size)
}
cat("Done.\n")


reference_scenarios <- c("clustGeo-50-25", "clustGeo-50-10", "clustGeo-50-90", "1-cellSq", "5-cellSq")

# --- MODIFIED: Removed regularization parameter definitions ---
# lambda_reg_beta = 0.1
# lambda_reg_alpha = 0.1


# Initialize lists to store summary statistics for efficient CSV export
descriptive_stats_list <- list()
performance_metrics_list <- list()
estimated_params_list <- list() # New dataframe for estimated parameters
perf_counter <- 1

# Loop order is inverted. Outer loop is now training_sizes.
for (current_train_size in training_sizes) {
    
    cat(paste("\n--- Running simulation for Training Size:", current_train_size, "---\n"))
    
    # Get the pre-sampled training set
    train_pts_df <- training_sets_list[[as.character(current_train_size)]]

    # Inner loop is reference_scenarios
    for (ref_scenario_name in reference_scenarios) {

        cat(paste("\n  --- REFERENCE:", ref_scenario_name, "---\n"))
        
        # NOTE: We now use the 'train_pts_df' defined in the outer loop,
        # ensuring all reference scenarios for a given size use the same points.

        # a. Simulate ground-truth point observation dataset
        if (startsWith(ref_scenario_name, "clustGeo")) {
            train_pts_covs <- terra::extract(state_cov, train_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
            D0 <- dist(train_pts_covs)
            D1 <- dist(train_pts_df[, c("x", "y")])
            D0 <- (D0 - min(D0)) / (max(D0) - min(D0))
            D1 <- (D1 - min(D1)) / (max(D1) - min(D1))
            ref_tree <- hclustgeo(D0, D1, alpha = 0.5)

            ref_k_proportion <- as.numeric(sub(".*-", "", ref_scenario_name)) / 100
            ref_k <- floor(ref_k_proportion * nrow(train_pts_df))
            cat(paste("Using k =", ref_k, "for clustGeo reference clustering.\n"))
            ref_clusters <- cutree(ref_tree, k = ref_k)
        } else {
            cell_side <- as.numeric(gsub("-cellSq", "", ref_scenario_name))
            cat(paste("Using", cell_side, "x", cell_side, "grid for reference clustering.\n"))
            ref_clusters <- create_cellSq_clusters(train_pts_df, landscape, cell_side)
        }

        train_pts_df$ref_site_id <- ref_clusters

        # Determine the true occupancy of reference sites
        ref_polygons_sf <- train_pts_df %>%
            st_as_sf(coords = c("x", "y"), crs = st_crs(landscape)) %>%
            group_by(ref_site_id) %>%
            summarise(n_pts = n(), geometry = st_combine(geometry), .groups = 'drop') %>%
            mutate(
                geometry = case_when(
                    n_pts >= 3 ~ st_convex_hull(geometry),
                    n_pts > 0    ~ st_buffer(geometry, dist = 0.5)
                )
            )
        ref_site_true_occupancy <- apply(st_intersects(ref_polygons_sf, individuals_sf, sparse = FALSE), 1, any)

        # Calculate true site-level lambda and psi by extracting weighted means from true maps
        true_site_lambda <- terra::extract(true_lambda_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]
        true_site_psi <- terra::extract(true_psi_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]

        # --- UPDATED SIMULATION OF TRAINING DETECTIONS ---
        # Get true lambda for each point's reference site
        site_lambda_map <- data.frame(ref_site_id = ref_polygons_sf$ref_site_id, true_site_lambda = true_site_lambda)
        # Ensure ref_site_id is present for the join
        ref_polygons_sf <- ref_polygons_sf %>% mutate(ref_site_id = 1:nrow(.))
        site_lambda_map <- data.frame(ref_site_id = ref_polygons_sf$ref_site_id, true_site_lambda = true_site_lambda)
        train_pts_df <- train_pts_df %>%
            left_join(site_lambda_map, by = "ref_site_id")
            
        # --- MODIFIED: Removed true_log_lambda ---
        train_pts_df <- train_pts_df %>%
            mutate(
                site_is_occupied = ref_site_true_occupancy[ref_site_id],
                obs_cov = runif(n())
            )
            
        # Handle -Inf (NO LONGER NEEDED)
        # train_pts_df$true_log_lambda[is.infinite(train_pts_df$true_log_lambda)] <- -100
            
        # --- MODIFIED: Removed abundance effect from detection logit ---
        p_detection_logit <- alpha_true[1] + alpha_true[2] * train_pts_df$obs_cov
        p_detection <- plogis(p_detection_logit)
            
        train_pts_df$detection <- rbinom(nrow(train_pts_df), 1, p_detection)
        train_pts_df$detection <- train_pts_df$detection * train_pts_df$site_is_occupied
            
        # --- MODIFIED: Removed true_log_lambda from cleanup ---
        # Clean up extra columns
        train_pts_df <- train_pts_df %>% select(-true_site_lambda)
        # --- END UPDATE ---


        points_per_cluster <- train_pts_df %>% count(ref_site_id) %>% pull(n)

        new_desc_row <- data.frame(
            "Reference.Clustering" = ref_scenario_name,
            "Training.Size" = current_train_size,
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
        # 5. CLUSTERING SCENARIOS: Cluster Points, Aggregate Data, and Analyze
        ########################################################################

        # a. Define the different clustGeo scenarios to test
        if (exists("ref_tree")) {
                n_unique_pts <- nrow(unique(train_pts_df[, c("x", "y")]))
                k_values <- c("clustGeo-50-25" = floor(0.25 * n_unique_pts), "clustGeo-50-90" = floor(0.90 * n_unique_pts), "clustGeo-50-10" = floor(0.10 * n_unique_pts))
                k_values <- k_values[c("clustGeo-50-25", "clustGeo-50-90", "clustGeo-50-10")]
                cluster_scenarios <- lapply(k_values, function(k) cutree(ref_tree, k))
        } else {
                train_pts_covs <- terra::extract(state_cov, train_pts_df[,c("x","y")])[, "state_cov", drop=FALSE]
                D0 <- dist(train_pts_covs); D1 <- dist(train_pts_df[, c("x", "y")])
                D0 <- (D0 - min(D0)) / (max(D0) - min(D0)); D1 <- (D1 - min(D1)) / (max(D1) - min(D1))
                ref_tree <- hclustgeo(D0, D1, alpha = 0.5)
                n_unique_pts <- nrow(unique(train_pts_df[, c("x", "y")]))
                k_values <- c("clustGeo-50-25" = floor(0.25 * n_unique_pts), "clustGeo-50-90" = floor(0.90 * n_unique_pts), "clustGeo-50-10" = floor(0.10 * n_unique_pts))
                k_values <- k_values[c("clustGeo-50-25", "clustGeo-50-90", "clustGeo-50-10")]
                cluster_scenarios <- lapply(k_values, function(k) cutree(ref_tree, k))
        }

        # b. Add cellSq clustering scenarios
        cellSq_1_clusters <- create_cellSq_clusters(train_pts_df, landscape, 1)
        cellSq_5_clusters <- create_cellSq_clusters(train_pts_df, landscape, 5)

        # c. Combine all scenarios into a single named list for processing
        all_scenarios <- c(cluster_scenarios, "1-cellSq" = list(cellSq_1_clusters), "5-cellSq" = list(cellSq_5_clusters))

        # d. Run the pipeline for each scenario using a for loop to access scenario names
        results <- list()

        for (current_scenario_name in names(all_scenarios)) {
            cat(paste("    - Method:", current_scenario_name, "\n"))
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
                polygons_df <- train_pts_df %>%
                    st_as_sf(coords = c("x", "y"), crs = st_crs(landscape)) %>%
                    group_by(current_site_id) %>%
                    summarise(geometry = st_combine(geometry), .groups = 'drop') %>%
                    mutate(
                        n_pts = site_data$n_pts[match(current_site_id, site_data$current_site_id)],
                        geometry = case_when(n_pts >= 3 ~ st_convex_hull(geometry), n_pts > 0 ~ st_buffer(geometry, dist = 0.5))
                    )
            }

            weights_extract <- terra::extract(state_cov, vect(polygons_df), weights = TRUE, cells = TRUE)
            w_matrix <- matrix(0, nrow = n_sites, ncol = ncell(landscape))
            for(i in seq_len(nrow(weights_extract))){
                site_idx <- weights_extract[i, "ID"]
                cell_idx <- weights_extract[i, "cell"]
                actual_site_id <- polygons_df$current_site_id[site_idx]
                w_matrix_row_idx <- which(site_data$current_site_id == actual_site_id)
                w_matrix[w_matrix_row_idx, cell_idx] <- weights_extract[i, "weight"]
            }
            
            # --- MODIFIED: Removed regularization parameters from function call ---
            analysis_results <- fit_and_evaluate_occuN(y_data, obs_covs_for_umf, w_matrix, state_cov, val_data=val_data_list, optimizer_method = selected_optimizer)

            # Calculate ARI and add it to the results list for plotting
            current_ari <- ARI(ref_clusters, cl)
            analysis_results$ari <- current_ari

            # --- UPDATED PARAMETER ROW ---
            # Store estimated parameters for CSV export
            # --- MODIFIED: Removed alpha_lambda ---
            new_params_row <- data.frame(
                "Reference.Clustering" = ref_scenario_name,
                "Training.Size" = current_train_size,
                "Clustering" = current_scenario_name,
                "beta_intercept" = analysis_results$metrics$beta_intercept,
                "beta_1" = analysis_results$metrics$beta_1,
                "alpha_intercept" = analysis_results$metrics$alpha_intercept,
                "alpha_1" = analysis_results$metrics$alpha_1
            )
            estimated_params_list[[length(estimated_params_list) + 1]] <- new_params_row
            # --- END UPDATE ---

            # Calculate and store performance metrics
            pred_cell_stats <- global(c(analysis_results$pred_lambda_map, analysis_results$pred_psi_map), c("mean", "sd"), na.rm=TRUE)
            pred_site_lambda <- terra::extract(analysis_results$pred_lambda_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]
            pred_site_psi <- terra::extract(analysis_results$pred_psi_map, vect(ref_polygons_sf), fun = "mean", weights = TRUE)[,2]

            new_perf_row <- data.frame(
                "No." = perf_counter,
                "Reference.Clustering" = ref_scenario_name,
                "Training.Size" = current_train_size,
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
        # 6. PLOT GENERATION (7x4 GRID)
        #################################

        original_scenario_order <- names(results)
        other_scenarios <- setdiff(original_scenario_order, ref_scenario_name)
        plot_order <- c(ref_scenario_name, other_scenarios)
        results <- results[plot_order]

        # Save PNG to output folder
        # Filename includes training size
        png_filename <- file.path(output_dir, paste0("ref=", ref_scenario_name, "_train=", current_train_size, ".png"))
        png(png_filename, width = 14, height = 20, units = "in", res = 300)
        cat(paste("Generating plot:", png_filename, "\n"))

        # # Make color scaling robust to NA/Inf values
        # valid_pred_maps <- Filter(function(res) !all(is.na(values(res$pred_lambda_map))), results)
        # max_vals <- c(global(true_lambda_map, "max", na.rm = TRUE)$max,
        #               sapply(valid_pred_maps, function(res) global(res$pred_lambda_map, "max", na.rm = TRUE)$max))
        # lambda_max <- max(unlist(max_vals), na.rm = TRUE)
        # if (!is.finite(lambda_max)) lambda_max <- 10 # Fallback value

        # Calculate the maximum abundance from ONLY the reference map
        ref_max_val <- global(true_lambda_map, "max", na.rm = TRUE)$max

        # Set the new lambda_max to be % of the reference max
        lambda_max <- ref_max_val * 1.00
        
        # Fallback in case ref_max_val was NA, Inf, or 0
        if (!is.finite(lambda_max) || lambda_max == 0) {
            warning("Could not determine a finite reference max; falling back to 10.")
            lambda_max <- 10 
        }

        # --- *** FIX 2: ROBUST LEGEND AXIS LOGIC *** ---
        # Generate a "pretty" set of breaks, but *exclude* the max
        # This gets us the "internal" ticks
        abun_breaks_internal <- pretty(c(0, lambda_max * 0.9))
        
        # Ensure breaks don't go below 0 or at/past the max
        abun_breaks_internal <- abun_breaks_internal[abun_breaks_internal >= 0 & abun_breaks_internal < lambda_max]
        
        # Add the lambda_max value as the *final* break
        abun_breaks <- c(abun_breaks_internal, lambda_max)
        
        # Create labels for all internal breaks
        abun_labels <- as.character(abun_breaks_internal)
        
        # Get the rounded max value for the last label
        lambda_max_rounded <- round(lambda_max, 1)
        
        # Add the final label with a ">" sign
        abun_labels <- c(abun_labels, paste0(">", lambda_max_rounded))
        
        # This list will be passed to the pax argument
        abun_axis_args <- list(
            at = abun_breaks,
            labels = abun_labels
        )
        # --- *** END FIX 2 *** ---


        # Define color palettes to be used in plots
        covariate_col_palette <- viridis(100, option = "turbo")
        abundance_col_palette <- viridis(100, option = "viridis")
        occupancy_col_palette <- viridis(100, option = "magma")
       

        par(cex = 2.0)
        layout(matrix(1:28, nrow = 7, byrow = TRUE), widths = c(1, 1, 1, 0.8), heights = c(1,1, 1, 1, 1, 1,1))

        # c(bottom, left, top, right) vector.
        margin_default <- c(0,0,0,0)
        text_offset = -0.75

        # ROW 1: Ground Truth Point Pattern, Abundance, and Occupancy
        par(mar = margin_default)

        # Plot state_cov as a background with no legend to ensure alignment
        # This plot (col 1) HAS axes by default because legend=FALSE
        plot(state_cov, main = "Point Pattern", legend = FALSE, col = "white")

        # Add the individual points to the plot, only if there are any
        if(nrow(all_individuals_df) > 0) {
          points(all_individuals_df$x, all_individuals_df$y, pch = 16, col = "red", cex = 0.5)
        }

        par(mar = margin_default)
        # Abundance Pattern
        # --- *** AXIS FIX: Added axes = TRUE *** ---
        plot(true_N_map, main = "Abundance Pattern", col = viridis(100, option="cividis"), axes = TRUE)
        text(x = xyFromCell(true_N_map, 1:ncell(true_N_map))[,1],
            y = xyFromCell(true_N_map, 1:ncell(true_N_map))[,2],
            labels = values(true_N_map),
            cex = 0.6, col = "white")

        # Occurrence Pattern
        true_occ_map <- ifel(true_N_map > 0, 1, 0)
        # This plot (col 3) HAS axes by default because legend=FALSE
        plot(true_occ_map, main = "Occurrence Pattern", col = c("black", "white"), legend = FALSE)


        # Text column for ground truth
        par(mar=c(0,0,0,0))
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
        text(text_offset, 0.7, "Reference Species", pos=4, font=2, cex=1.2)
        text(text_offset, 0.5, paste("Total Individuals:", sum(true_N_j)), pos=4, cex=1.2)
        text(text_offset, 0.35, paste("Realized Mean Abundance:", round(mean(true_N_j), 2)), pos=4, cex=1.2)
        text(text_offset, 0.2, paste("Realized Occupancy:", round(mean(true_N_j > 0), 2)), pos=4, cex=1.2)

        # ROW 2: Reference
        par(mar = margin_default)
        # --- *** AXIS FIX: Added axes = TRUE *** ---
        plot(state_cov, main = "Covariate & Observations", col = covariate_col_palette, axes = TRUE)
        points(train_pts_df$x, train_pts_df$y, pch = 16, col = "blue", cex = 0.5)
        points(val_pts_df$x, val_pts_df$y, pch = 17, col = "magenta1", cex = 0.5)
        legend("bottom", legend = c("Train", "Val"), col = c("blue", "magenta1"), pch = c(16, 17), horiz = TRUE, bg="white", cex=1.0)

        par(mar = margin_default)
        # --- *** AXIS FIX: Added axes = TRUE *** ---
        # Use pax (plotaxt) to control the legend's axis
        plot(true_lambda_map, main = "Reference Lambda", col = abundance_col_palette, range = c(0, lambda_max), extend = "max", pax = abun_axis_args, axes = TRUE)
        
        # --- *** AXIS FIX: Added axes = TRUE *** ---
        plot(true_psi_map, main = "Reference Psi", col = occupancy_col_palette, range=c(0,1), axes = TRUE)

        # MODIFIED: Text block at (2,4) shows simulation stats instead of parameters
        par(mar=c(0,0,0,0))
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
        text(text_offset, 0.7, "Simulation Stats", pos=4, font=2, cex=1.2)
        text(text_offset, 0.5, paste("Training Pts.:", nrow(train_pts_df)), pos=4, cex=1.2)
        text(text_offset, 0.35, paste("Validation Pts.:", nrow(val_pts_df)), pos=4, cex=1.2)
        text(text_offset, 0.2, paste("Cell Mean Exp. Abun.:", round(mean(true_lambda_j), 2)), pos=4, cex=1.2)
        text(text_offset, 0.05, paste("Cell Mean Occu. Prob.:", round(mean(true_psi_j), 2)), pos=4, cex=1.2)

        # ROWS 3-7: Clustering Scenarios
        scenario_names <- names(results)
        for(i in 1:5){
            res <- results[[i]]
            metrics <- res$metrics
            current_scenario_name <- scenario_names[i]

            trailing_ref_str <- paste(current_scenario_name, "Sites")
            if (current_scenario_name == ref_scenario_name) {
                trailing_ref_str <- paste(trailing_ref_str, "(reference clustering)")
            }

            par(mar = margin_default)
            # This plot (col 1) HAS axes by default because legend=FALSE
            plot(state_cov, main = trailing_ref_str, legend = FALSE, col = "white")
            plot(res$polygons$geometry, add = TRUE, border = "gray30", col = alpha("gray", 0.3))
            points(train_pts_df$x, train_pts_df$y, pch = 16, col = "blue", cex = 0.5)

            # Check for failed fits before plotting rasters
            if(is.null(res$fit) || all(is.na(values(res$pred_lambda_map)))) {
                
                # This is a failed fit. Plot "Fit Failed" for both abundance and occupancy.
                par(mar = margin_default)
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", main = expression(paste(lambda, "Estimated from", current_scenario_name)))
                text(0, 0, "Fit Failed", cex=1.5)
                
                par(mar = margin_default)
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", main = expression(paste(psi, "Estimated from", current_scenario_name)))
                text(0, 0, "Fit Failed", cex=1.5)
                
            } else {
                
                # This is a successful fit. Plot the rasters.
                par(mar = margin_default)

                # Create a copy of the map for plotting
                lambda_map_to_plot <- res$pred_lambda_map
                
                # Get the original values
                v_cleaned <- values(lambda_map_to_plot)
                
                # 1. Replace all non-finite values (NA, Inf, NaN) with lambda_max.
                #    This ensures the vector is now fully finite.
                v_cleaned[!is.finite(v_cleaned)] <- lambda_max
                
                # 2. Now that it's finite, clamp any values > lambda_max
                #    (This catches the huge MSE numbers, e.g., 5.86e+29)
                v_cleaned[v_cleaned > lambda_max] <- lambda_max
                
                # 3. Clamp any values < 0 (just in case)
                v_cleaned[v_cleaned < 0] <- 0
                
                # 4. Set the fully cleaned and clamped values back into the raster
                values(lambda_map_to_plot) <- v_cleaned
                
                # 5. Plot the clamped map. 
                # --- *** AXIS FIX: Added axes = TRUE *** ---
                plot(lambda_map_to_plot, main = paste("Abundance Estimated from", current_scenario_name), 
                     col = abundance_col_palette, range = c(0, lambda_max), 
                     extend = "max", pax = abun_axis_args, axes = TRUE)
                
                # Plot the occupancy map
                # --- *** AXIS FIX: Added axes = TRUE *** ---
                plot(res$pred_psi_map, main = paste("Occupancy Estimated from", current_scenario_name), col = occupancy_col_palette, range=c(0,1), axes = TRUE)
            }

            # MODIFIED: Text blocks for scenarios no longer show estimated parameters
            par(mar=c(0,0,0,0))
            plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab="", ylab="")
            
            text(text_offset, 0.8, "Metrics", pos=4, font=2, cex=1.2)
            text(text_offset, 0.65, paste("ARI:", round(res$ari, 4)), pos=4, cex=1.2)
            
            # Check for NA before formatting each MSE value.
            lambda_ls_mse_text <- if (is.na(metrics$lambda_ls_mse)) "NA" else formatC(metrics$lambda_ls_mse, format = "e", digits = 2)
            psi_ls_mse_text <- if (is.na(metrics$psi_ls_mse)) "NA" else formatC(metrics$psi_ls_mse, format = "e", digits = 2)
            lambda_pt_mse_text <- if (is.na(metrics$lambda_pt_mse)) "NA" else formatC(metrics$lambda_pt_mse, format = "e", digits = 2)
            psi_pt_mse_text <- if (is.na(metrics$psi_pt_mse)) "NA" else formatC(metrics$psi_pt_mse, format = "e", digits = 2)
            
            text(text_offset, 0.45, paste("Cell Abun. MSE:", lambda_ls_mse_text), pos=4, cex=1.2)
            text(text_offset, 0.3, paste("Cell Occu. MSE:", psi_ls_mse_text), pos=4, cex=1.2)
            text(text_offset, 0.15, paste("Val. Pts. Abun. MSE:", lambda_pt_mse_text), pos=4, cex=1.2)
            text(text_offset, 0.00, paste("Val. Pts. Occ. MSE:", psi_pt_mse_text), pos=4, cex=1.2)

            # Check for NA before formatting AUROC
            auroc_text <- if (is.na(metrics$auroc)) "NA" else round(metrics$auroc, 4)

            text(text_offset, -0.2, paste("AUC:", auroc_text), pos=4, cex=1.2)


            par(mar = margin_default)
        }
        dev.off()
    
    } # END of reference_scenarios loop
} # END of training_sizes loop

#############################
# 7. WRITE SUMMARY CSV FILES
#############################

# Convert lists to data frames
descriptive_stats_df <- do.call(rbind, descriptive_stats_list)
performance_metrics_df <- do.call(rbind, performance_metrics_list)
estimated_params_df <- do.call(rbind, estimated_params_list)

# Save CSVs to output folder
write.csv(descriptive_stats_df, file.path(output_dir, "descriptive_stats.csv"), row.names = FALSE)
write.csv(performance_metrics_df, file.path(output_dir, "performance_metrics.csv"), row.names = FALSE)
write.csv(estimated_params_df, file.path(output_dir, "estimated_parameters.csv"), row.names = FALSE)


#############################
# 8. GENERATE SUMMARY PLOT
#############################

cat("\n--- Generating summary plot... ---\n")

# Filter for the "diagonal" cases where reference clustering matches analysis clustering
summary_plot_df <- performance_metrics_df %>%
  filter(as.character(Reference.Clustering) == as.character(Clustering)) %>%
  mutate(Training.Size = as.numeric(Training.Size)) %>%
  arrange(Training.Size)

png_filename_summary <- file.path(output_dir, "summary_mse_vs_training_size.png")
png(png_filename_summary, width = 14, height = 7, units = "in", res = 300)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)

methods <- unique(summary_plot_df$Clustering)
colors <- viridis(length(methods))
method_names <- as.character(methods)

# --- *** FIX 1: HANDLE NON-FINITE YLIM VALUES *** ---
# Plot 1: Landscape Lambda MSE
finite_lambda_mse <- summary_plot_df$Lambda.LS.MSE[is.finite(summary_plot_df$Lambda.LS.MSE)]
y_range_lambda <- if(length(finite_lambda_mse) > 0) range(finite_lambda_mse, na.rm = TRUE) else c(0, 1)
# Handle case where range might still be Inf
if (any(!is.finite(y_range_lambda))) y_range_lambda <- c(0, 1)


plot(NULL, xlim = range(training_sizes), ylim = y_range_lambda, 
     xlab = "Training Set Size", ylab = "Landscape Lambda MSE", 
     main = "Lambda MSE vs. Training Size")

for (i in 1:length(methods)) {
  method_data <- summary_plot_df %>% 
    filter(Clustering == methods[i]) %>% 
    arrange(Training.Size)
  lines(method_data$Training.Size, method_data$Lambda.LS.MSE, col = colors[i], type = "b", pch = 19, lwd = 2)
}
legend("topright", legend = method_names, col = colors, lty = 1, pch = 19, bty = "n", cex = 1.0)


# Plot 2: Landscape Psi MSE
finite_psi_mse <- summary_plot_df$Psi.LS.MSE[is.finite(summary_plot_df$Psi.LS.MSE)]
y_range_psi <- if(length(finite_psi_mse) > 0) range(finite_psi_mse, na.rm = TRUE) else c(0, 1)
# Handle case where range might still be Inf
if (any(!is.finite(y_range_psi))) y_range_psi <- c(0, 1)

plot(NULL, xlim = range(training_sizes), ylim = y_range_psi, 
     xlab = "Training Set Size", ylab = "Landscape Psi MSE", 
     main = "Psi MSE vs. Training Size")
# --- *** END FIX 1 *** ---

for (i in 1:length(methods)) {
  method_data <- summary_plot_df %>% 
    filter(Clustering == methods[i]) %>% 
    arrange(Training.Size)
  lines(method_data$Training.Size, method_data$Psi.LS.MSE, col = colors[i], type = "b", pch = 19, lwd = 2)
}
legend("topright", legend = method_names, col = colors, lty = 1, pch = 19, bty = "n", cex = 1.0)

dev.off()


cat("\n--- All simulations complete. ---\n")