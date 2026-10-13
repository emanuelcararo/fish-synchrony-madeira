################################################################################
# Manuscript: Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir 
# Author's: Emanuel R. Cararo
# Last update: 06/08/2026
# It's necessary to run the script "1_FUNCTIONS" before starting this data manipulation!!
# Working directory must be the repository root (folder containing Data/ and Scripts/).
################################################################################

# Load data ----
cpue = readxl::read_excel("Data/1_CPUE.xlsx")

# Create a list of matrices
quarters_species_matrices = create_species_matrices(
  wide_data = cpue, 
  local_col = "Local", 
  year_col = "Year",
  month_col = "Month",
  first_species_col = 4,
  impute_method = "none")

annual_species_matrices = lapply(quarters_species_matrices, function(mat) {
  row_names = rownames(mat)
  matrix_years  = as.numeric(substring(row_names, 1, 4))
  matrix_months = as.numeric(substring(row_names, 6, 7))
  hydro_groups = ifelse(matrix_months == 3, matrix_years - 1, matrix_years)
  rowsum(mat, group = hydro_groups, na.rm = TRUE)
})

quarters_species_matrices = quarters_species_matrices[
  sapply(quarters_species_matrices, function(mat) {
    occurs_all_locals = all(colSums(mat > 0, na.rm = TRUE) > 0)
    occurrence_prop = mean(mat > 0, na.rm = TRUE)
    occurs_all_locals && occurrence_prop >= 0.5
  })
]

annual_species_matrices = annual_species_matrices[
  sapply(annual_species_matrices, function(mat) {
    occurs_all_locals = all(colSums(mat > 0, na.rm = TRUE) > 0)
    occurrence_prop = mean(mat > 0, na.rm = TRUE)
    occurs_all_locals && occurrence_prop >= 0.5
  })
]


# Load the sampling site data 
locs = read.csv("Data/2_LOCALS.csv", header = TRUE, sep = ";")

coord = data.frame(lat = locs$lat, lon = locs$long, Locals = locs$Local)

# Extract habitat type and dam position (upstream/downstream)
habitat_type = locs$Habitat
dam_position = locs$Dam

# Compute geographic distance matrix (great-circle distance, km; ncf::gcdist) 
 geographic_distance = ncf::gcdist(locs$long, locs$lat)
 rownames(geographic_distance) = locs$Local
 colnames(geographic_distance) = locs$Local


# Create connectivity matrix: 1 if sites are on the same side of the dam, 0 otherwise 
n_sites = length(dam_position)
connectivity_matrix = matrix(0, n_sites, n_sites)

for (i in 1:n_sites) {
  for (j in 1:n_sites) {
    if (dam_position[i] == dam_position[j]) {
      connectivity_matrix[i, j] = 1
    }
  }
}

# Create dispersal difficulty matrix: geographic distance penalized by lack of connectivity
# If not connected, assign a high value (twice the maximum observed distance)
dispersal_difficulty = geographic_distance * connectivity_matrix
dispersal_difficulty[dispersal_difficulty == 0] = 2 * max(geographic_distance)



# 7.Load environmental data (water level time series) ----
hydro_data = read.csv("Data/3_WATER_LEVEL.csv", header = T,
                      sep = ";")

# Convert 'HydroYear' column to row names and convert to matrix
hydro_data = hydro_data |> 
  column_to_rownames(var = "HydroYear") |>   
  as.matrix() 

# Apply log transformation to stabilize variance
hydro_data = log(hydro_data)

# Detrend each site's time series to remove long-term trends
hydro_data = pracma::detrend(hydro_data, tt = "linear")

# 8.Compute environmental synchrony)----
env_synchrony = cor(hydro_data, method = "spearman")

#mean(ecodist::lower(env_synchrony))



# 9.Load trait and phylogeny data ----
sp_trait = read_xlsx("Data/4_TRAITS.xlsx", sheet = 1)
sp_trait$Species = gsub(" ", "_", sp_trait$Species)

tree = phytools::read.newick("Data/fish_tree.tre")

# 10.Calculate density dependence matrices ---- 
DD.list = list() 
for(s in 1:length(annual_species_matrices)){
  species_name = names(annual_species_matrices)[s]
  matriz.r = t(annual_species_matrices[[s]])
  corr.dd = numeric(length = nrow(matriz.r))
  names(corr.dd) = rownames(matriz.r)
  beta.dd = numeric(length = nrow(matriz.r))
  p_value = numeric(length = nrow(matriz.r))
  names(beta.dd) = rownames(matriz.r)
  carrying.capacity = rowMeans(log(matriz.r+1)) 
  for(j in 1:nrow(matriz.r)){
    ts_t = c(matriz.r[j,], NA)
    ts_tmais1 = c(matriz.r[j,-1], NA, NA)
    corr.dd[j] = cor(log(ts_t+1), log(ts_tmais1+1), method = "sp", use = "pairwise.complete.obs")
    if(is.na(corr.dd[j])){
      next
    }
    modeldd = summary(lm(log(ts_tmais1+1) ~ log(ts_t+1)))
    beta.dd[j] = modeldd$coefficients[2,1]
    p_value[j] = modeldd$coefficients[2,4]
  }
  DD.list[[species_name]] = cbind(corr.dd, beta.dd, p_value, carrying.capacity)
}

DD.betacc = list()
for (i in 1:length(DD.list)) {
  metrics_matrix = DD.list[[i]]
  data_for_dist = na.omit(metrics_matrix[, c(2, 3)])
  DD.betacc[[i]] = dist(scale(data_for_dist), method = 'euclidean')
  names(DD.betacc)[i] = names(DD.list)[i] 
}

density_dependence = DD.betacc


