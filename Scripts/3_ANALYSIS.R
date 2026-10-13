##############################################################################
# Manuscript: Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir
# Author's: Emanuel R. Cararo
# Last update: 13/10/2026
# It's necessary to run the script "1_FUNCTIONS and 2_DATA" before starting this analysis!!
# Results tables are written to Outputs/ (created if missing).
##############################################################################

dir.create("Outputs", showWarnings = FALSE)

# 1. Cross-correlation Synchrony ----
annual_spatial_synchrony_summary = data.frame(
  Species = names(annual_species_matrices),
  Mean_sync = numeric(length(annual_species_matrices)),
  CI_lower = numeric(length(annual_species_matrices)),
  CI_upper = numeric(length(annual_species_matrices)))

for(i in 1:length(annual_species_matrices)){
  sp_matrix = annual_species_matrices[[i]]
  times = seq(1:nrow(sp_matrix))
  sp_matrix_cleaned = cleandat(t(sp_matrix), times, clev = 3)$cdat
  res = mSynch2(sp_matrix_cleaned, resamp = 1000)
  annual_spatial_synchrony_summary$Mean_sync[i] = res$real
  annual_spatial_synchrony_summary$CI_lower[i] = quantile(res$boot, probs = 0.025, na.rm = TRUE)
  annual_spatial_synchrony_summary$CI_upper[i] = quantile(res$boot, probs = 0.975, na.rm = TRUE)
}


# 2. Asymmetry of synchrony ----
# Define bounds for partial Spearman correlation
lower_bounds = c(0, 0.5)  # Low abundance tail
upper_bounds = c(0.5, 1)  # High abundance tail

# Store correlation matrices by species
tail_sync_lower = list()
tail_sync_upper = list()
sync_total = list()
tail_strength = list()

# Summary of synchrony measures 

annual_tail_synchrony_summary = data.frame( Species = names(annual_species_matrices), Mean_sync = numeric(length(annual_species_matrices)), Lower_tail_sync = numeric(length(annual_species_matrices)), Upper_tail_sync = numeric(length(annual_species_matrices)), Tail_association = numeric(length(annual_species_matrices)) )

for (k in seq_along(annual_species_matrices)) {
  sp_name = names(annual_species_matrices)[k]
  sp = annual_species_matrices[[sp_name]]
  sp = cleandat(t(sp), seq_len(nrow(sp)), clev = 3)$cdat
    # Total correlation
  corr_total = cor(t(sp), method = "spearman")
  # Initialize tail correlation matrices
  corr_lb = matrix(NA, nrow(sp), nrow(sp))
  corr_ub = matrix(NA, nrow(sp), nrow(sp))
  # Calculate lower- and upper-tail correlations
  for (i in 2:nrow(sp)) {
    for (j in 1:(i - 1)) {
      corr_lb[i, j] = partialSpearman(
        sp[i, ], 
        sp[j, ], 
        lower_bounds
      )
      corr_ub[i, j] = partialSpearman(
        sp[i, ], 
        sp[j, ], 
        upper_bounds
      )
    }
  }
 # Tail association = upper-tail - lower-tail
  corr_strength = corr_ub - corr_lb
  # Store matrices
  sync_total[[sp_name]] = corr_total
  tail_sync_lower[[sp_name]] = corr_lb
  tail_sync_upper[[sp_name]] = corr_ub
  tail_strength[[sp_name]] = corr_strength
  # Extract lower triangle
  lower_tri = lower.tri(corr_total)
  # Calculate mean synchrony measures
  annual_tail_synchrony_summary$Mean_sync[k] =
    mean(corr_total[lower_tri], na.rm = TRUE)
  annual_tail_synchrony_summary$Lower_tail_sync[k] =
    mean(corr_lb[lower_tri], na.rm = TRUE)
  annual_tail_synchrony_summary$Upper_tail_sync[k] =
    mean(corr_ub[lower_tri], na.rm = TRUE)
  annual_tail_synchrony_summary$Tail_association[k] =
    mean(corr_strength[lower_tri], na.rm = TRUE)
}

writexl::write_xlsx(annual_spatial_synchrony_summary, "Outputs/mean_spatial_synchrony.xlsx")
writexl::write_xlsx(annual_tail_synchrony_summary, "Outputs/tail_synchrony_summary.xlsx")

# Manual check of tail dependence 

#sp1 = annual_species_matrices$Plagioscion.squamosissimus
# sp1 = cleandat(t(sp1), times, clev = 3)$cdat
# sp1 = t(sp1)
# head(sp1)
# lb = c(0, 0.5)
# ub = c(0.5, 1)
# n=nrow(sp1)
# x=rank(sp1[,1])/(n+1)
# y=rank(sp1[,2])/(n+1)
# x_b= x[x+y > (2*min(ub)) & x+y < (2*max(ub))]
# y_b=y[x+y > (2*min(ub)) & x+y < (2*max(ub))]
# r=sum((x_b-mean(x))*(y_b-mean(y)))/((n-1)*sqrt(var(x)*var(y)))
# tail_sync_upper$Plagioscion.squamosissimus

# Surrogate-based statistical test
test_results = data.frame(Species = character(),
                          Res = numeric(),
                          lowerq = numeric(),
                          upperq = numeric(),
                          empDiff = numeric(),
                          stringsAsFactors = FALSE)

for (i in 1:length(annual_species_matrices)) {
  set.seed(263)
  sp = annual_species_matrices[[i]]
  sp = cleandat(t(sp), seq_len(nrow(sp)), clev = 3)$cdat
  sp = t(sp)
  tail_diff_emp = tail_strength[[i]]
  lb_mat = tail_sync_lower[[i]]
  ub_mat = tail_sync_upper[[i]]
  surrogates = ncsurrog(m = sp, numsurrog = 1000)
  sum_diff_emp = sum(tail_diff_emp, na.rm = TRUE)
  cor_emp = cor(lower(lb_mat), lower(ub_mat), use = "pairwise.complete.obs")
  surrogate_diffs = rep(NA, dim(surrogates)[3])
  surrogate_corrs = rep(NA, dim(surrogates)[3])
  for(r in 1:dim(surrogates)[3]) {
    N = dim(surrogates)[2]
    dat = t(surrogates[,,r])
    corr_lb = matrix(0, N, N)
    corr_ub = matrix(0, N, N)
    filt = ifelse(cor(t(dat), method = "spearman") >= 0, 1, NA)
    for (ii in 2:N) {
      for (jj in 1:(ii - 1)) {
        corr_lb[ii, jj] = partialSpearman(dat[ii,], dat[jj,], lower_bounds)
        corr_ub[ii, jj] = partialSpearman(dat[ii,], dat[jj,], upper_bounds)
      }
    }
    corr_lb = corr_lb + t(corr_lb) * filt
    corr_ub = corr_ub + t(corr_ub) * filt
    surrogate_diffs[r] = sum(corr_ub - corr_lb, na.rm = TRUE)
    surrogate_corrs[r] = 1 - cor(corr_lb[lower.tri(corr_lb)],
                                 corr_ub[lower.tri(corr_ub)],
                                 use = "pairwise.complete.obs")
    }
  p_tail_dependence = rank(c(sum_diff_emp, surrogate_diffs))[1]/1001
  percentile_limits = quantile(surrogate_diffs, probs = c(0.025, 0.975))
  test_results = rbind(test_results,
                       data.frame(Species = names(annual_species_matrices)[i],
                                  Res = p_tail_dependence,
                                  lowerq = percentile_limits[1],
                                  upperq = percentile_limits[2],
                                  empDiff = sum_diff_emp))
  }

writexl::write_xlsx(test_results, "Outputs/tail_dependence_surrogate_test.xlsx")

# Is there a significant difference between lower and upper tail synchrony?
signif_tail_dependence = subset(test_results, empDiff <= lowerq | empDiff >= upperq)
signif_tail_dependence

# 3. Multiple Regression on distance Matrices (MRM) ----
# Prepare a data frame to store results from regression models
regression_results = data.frame(
  Specie = as.character(),
  dd_m2 = as.numeric(), se_dd_m2 = as.numeric(),
  disp_diff_m2 = as.numeric(), se_disp_diff_m2 = as.numeric(),
  env_m2 = as.numeric(), se_env_m2 = as.numeric(),
  f_stat_m2 = as.numeric(), r_2_m2 = as.numeric())

# Run the loop for each species
for(i in 1:length(sync_total)){
  mod_1 = MRM(lower(sync_total[[i]]) ~lower(density_dependence[[i]])+lower(dispersal_difficulty)  + lower(env_synchrony))
  mod_2 = lm(lower(sync_total[[i]]) ~ lower(density_dependence[[i]])+lower(dispersal_difficulty) + lower(env_synchrony))
  std.coefs = lm.beta(mod_2)
  summ.std = summary(std.coefs)
  std.se = summ.std$coefficients[,3]*(summ.std$coefficients[,2]/mod_2$coefficients)
  f_stat = mod_1$F.test[2] # p-value of the MRM F test (permutation)
  r_2 = mod_1$r.squared[1]
  # Coefficient rows follow the formula order: 2 = density dependence, 3 = dispersal difficulty, 4 = water level
  dd = summ.std$coefficients[2,2] 
  disp_diff = summ.std$coefficients[3,2] 
  env  = summ.std$coefficients[4,2]
  regression_results = rbind(regression_results,
                             data.frame(
                               Specie = names(sync_total)[i],
                               dd = dd, se_dd = std.se[2],
                               disp_diff = disp_diff, se_disp_diff = std.se[3],
                               env = env, se_env = std.se[4],
                               f_stat = f_stat, r_2 = r_2))
}

regression_results$Specie = stringr::str_replace_all(regression_results$Specie, "\\.", " ")
writexl::write_xlsx(regression_results, "Outputs/mrm_results.xlsx")

# density_dependence[[i]] from the last loop iteration instead of the matching species.
# m1 = MRM(lower(sync_total$Pterygoplichthys.disjunctivus) ~lower(density_dependence[[i]])+lower(dispersal_difficulty)  + lower(env_synchrony))
# m2 = MRM(lower(sync_total$Pellona.flavipinnis) ~lower(density_dependence[[i]])+lower(dispersal_difficulty)  + lower(env_synchrony))


# 4. Phylogentic last square models ----
# Creating distance matrix of species traits 

trait1= as.data.frame(sp_trait$Diet)
rownames(trait1) = sp_trait$Species

trait2 = as.data.frame(sp_trait$Habitat)
rownames(trait2) = sp_trait$Species

trait3 = as.data.frame(sp_trait$Migratory)
rownames(trait3) = sp_trait$Species

trait4 = as.data.frame(sp_trait$`Parental care`)
rownames(trait4) = sp_trait$Species

trait5 = as.data.frame(sp_trait$Fertilization)
rownames(trait5) = sp_trait$Species


ktab_list = ktab.list.df(list(trait1, trait2,trait3,trait4,trait5))

mat_dist_traits = dist.ktab(ktab_list, type = c("N", "N", "N", "N", "N"))

# Performing the PCoA to synthesize the data set on axes 

pcoa.traits = pcoa(mat_dist_traits)
pcoa.traits$values

sp_trait$References = NULL
sp_trait$PC1 = pcoa.traits$vectors[,1]
sp_trait$PC2 = pcoa.traits$vectors[,2]
sp_trait$Traits = paste(sp_trait$Migratory, sp_trait$Habitat)

# Models 
# Create dataframe to store the data
annual_spatial_synchrony_summary$Species = gsub("[. ]", "_", annual_spatial_synchrony_summary$Species)
#tree$tip.label

merged_data = merge(
  annual_spatial_synchrony_summary, 
  sp_trait, 
  by = "Species"
)


setdiff(annual_spatial_synchrony_summary$Species, sp_trait$Species)

# 2. Select and rename the specific columns you want
data = data.frame(
  Species     = merged_data$Species,
  PC1         = merged_data$PC1,
  PC2         = merged_data$PC2,
  Synchrony   = merged_data$Mean_sync,
  Body_length = merged_data$`Body length`,
  Troph_lvl   = merged_data$`Trophic level`
)

writexl::write_xlsx(data, "Outputs/data_pgls.xlsx")

# Linear model 
linear_model = lm(Synchrony ~ log(Body_length) + log(Troph_lvl) + PC1 + PC2, data = data)
anova(linear_model)
summary(linear_model)

# pgls models 

# Brownian Motion (BM)

pglsModel1 = gls(
  Synchrony ~ log(Body_length) + log(Troph_lvl) + PC1 + PC2,
  correlation = corBrownian(
    phy = tree,
    form = ~ Species
  ),
  data = data,
  method = "ML")

anova(pglsModel1)
coef(pglsModel1)

# Pagel's λ
pglsModel2 = gls(
  Synchrony ~ log(Body_length) + log(Troph_lvl) + PC1 + PC2,
  correlation = corPagel(
    value = 0.5,
    phy = tree,
    form = ~ Species,
    fixed = FALSE
  ),
  data = data,
  method = "ML",
  control = glsControl(
    opt = "optim"
  )
)

anova(pglsModel2)
coef(pglsModel2)
summary(pglsModel2)

# Ornstein-Uhlenbeck (OU)
pglsModel3 = gls(
  Synchrony ~ log(Body_length) + log(Troph_lvl) + PC1 + PC2,
  correlation = corMartins(
    value = 0.5,
    phy = tree,
    form = ~ Species,
    fixed = FALSE),
  data = data,
  method = "ML")

anova(pglsModel3)
coef(pglsModel3)

AIC(pglsModel1, pglsModel2, pglsModel3, linear_model)

# Save model summaries
sink("Outputs/trait_models_summary.txt")
cat("## PCoA of traits (relative eigenvalues)\n"); print(head(pcoa.traits$values))
cat("\n## Linear model\n"); print(anova(linear_model)); print(summary(linear_model))
cat("\n## PGLS Brownian motion\n"); print(anova(pglsModel1)); print(summary(pglsModel1))
cat("\n## PGLS Pagel's lambda\n"); print(anova(pglsModel2)); print(summary(pglsModel2))
cat("\n## PGLS Ornstein-Uhlenbeck\n"); print(anova(pglsModel3)); print(summary(pglsModel3))
cat("\n## AIC\n"); print(AIC(pglsModel1, pglsModel2, pglsModel3, linear_model))
sink()

# Wavelets ----

res_list = list()

  for (i in seq_along(quarters_species_matrices)) {
    current_sp_matrix = quarters_species_matrices[[i]]
    current_sp_name = names(quarters_species_matrices)[i]
    times = 1:nrow(current_sp_matrix)
    sp_cleaned = cleandat(t(current_sp_matrix), times, clev = 3)$cdat
    wpmf_result = wpmf(sp_cleaned, times, sigmethod = "fft")
    res_list[[current_sp_name]] = wpmf_result
  }


saveRDS(res_list, "Outputs/wpmf_results.rds")


### Synchrony networks  ----

cluster_results = vector("list", length(annual_species_matrices))
modularity_results = vector("list", length(annual_species_matrices))
adj_matrices_results = vector("list", length(annual_species_matrices))

names(cluster_results) = names(annual_species_matrices)
names(modularity_results) = names(annual_species_matrices)

for(i in 1:length(annual_species_matrices)){
  
  sp = annual_species_matrices[[i]]
  sp_d = cleandat(t(sp), seq_len(nrow(sp)), clev = 3)$cdat 
  sp_d = t(sp_d)
  sp_cor = cor(sp_d, method = "spearman")
  adj = psych::corr.test(sp_d,
                  method = "spearman",
                  adjust = "none")
  adj_signif = ifelse(adj$p <= 0.05, 1, 0)
  diag(adj_signif) = 0
  adj_final = ifelse(adj_signif == 1, sp_cor, 0)
  adj_matrices_results[[i]] = adj_final
  m = cluseigen(adj_final)
  n_clust = length(m)
  cluster_results[[i]] = m[[length(m)]]
  
  modularity_results[[i]] = wsyn::modularity(adj_final,
                                             membership = m[[length(m)]],
                                             decomp = TRUE)
}


# 2. Null model of total Modularity (Only for >1 cluster) ---

nperm = 1000
n_matrices = length(annual_species_matrices)

results_df = data.frame(
  Matrix = names(annual_species_matrices),
  Num_Clusters = NA_integer_,
  Observed_Q = NA_real_,
  Mean_Null_Q = NA_real_,
  SD_Null_Q = NA_real_,
  Null_Q_2.5 = NA_real_,
  Null_Q_97.5 = NA_real_,
  Z_score = NA_real_,
  P_value = NA_real_,
  stringsAsFactors = FALSE
)

for(i in seq_len(n_matrices)) {
  
  # Extract observed membership and count unique clusters
  # (If cluster_results[[i]] is the membership vector):
  obs_membership = cluster_results[[i]]
  n_clusters = length(unique(obs_membership[!is.na(obs_membership)]))
  
  results_df$Num_Clusters[i] = n_clusters
  
  # Condition: Run permutation ONLY if more than one cluster exists
  if (n_clusters > 1) {
    
    # Extract observed modularity
    modularity_OBS = modularity_results[[i]]$totQ
    
    # Detrend / clean time series
    sp = annual_species_matrices[[i]]
    sp_d = t(wsyn::cleandat(t(sp), seq_len(nrow(sp)), clev = 3)$cdat)
    
    null_Q = numeric(nperm)
    
    for(p in seq_len(nperm)) {
      # Reshuffle time series independently across time points
      sp_null = apply(sp_d, 2, sample)
      
      # Spearman correlation & significance test
      test_res = psych::corr.test(sp_null, method = "spearman", adjust = "none")
      # Threshold adjacency matrix
      adj_null_final = test_res$r * (test_res$p <= 0.05)
      diag(adj_null_final) = 0
      
      # Modularity & clustering
      m_null = wsyn::cluseigen(adj_null_final)
      membership_null = m_null[[length(m_null)]]
      
      if (length(unique(membership_null)) > 1) {
        null_Q[p] = wsyn::modularity(
          adj_null_final,
          membership = membership_null,
          decomp = TRUE
        )$totQ
      } else {
        null_Q[p] <- 0
      }
    }
    
    # Summary statistics
    mean_null = mean(null_Q, na.rm = TRUE)
    sd_null = sd(null_Q, na.rm = TRUE)
    
    results_df$Observed_Q[i]  = modularity_OBS
    results_df$Mean_Null_Q[i] = mean_null
    results_df$SD_Null_Q[i]  = sd_null
    results_df$Null_Q_2.5[i]  = quantile(null_Q, 0.025, na.rm = TRUE)
    results_df$Null_Q_97.5[i] = quantile(null_Q, 0.975, na.rm = TRUE)
    
    results_df$Z_score[i] = if (!is.na(sd_null) && sd_null > 0) {
      (modularity_OBS - mean_null) / sd_null
    } else {
      NA_real_
    }
    
    results_df$P_value[i] <- (sum(null_Q >= modularity_OBS, na.rm = TRUE) + 1) / (sum(!is.na(null_Q)) + 1)
  }
}

# View only the evaluated networks
evaluated_df = subset(results_df, Num_Clusters > 1)
writexl::write_xlsx(results_df, "Outputs/modularity_null_model.xlsx")

