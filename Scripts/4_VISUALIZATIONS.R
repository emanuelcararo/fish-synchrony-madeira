################################################################################
# Manuscript: Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir
# Author's: Emanuel R. Cararo
# Last update: 13/10/2026
# It's necessary to run the script "1_FUNCTIONS, 2_DATA and 3_ANALYSIS" before starting this graphic visualizations!!
################################################################################

dir.create("Plots", showWarnings = FALSE)

# 1. Plot of Mean Spatial Synchrony ----
annual_spatial_synchrony_summary$Crosses_Zero = ifelse( annual_spatial_synchrony_summary$CI_lower < 0 & annual_spatial_synchrony_summary$CI_upper > 0, "crosses_0", "does_not_cross_0")

annual_spatial_synchrony_summary$Species = gsub("[._]", " ", annual_spatial_synchrony_summary$Species)

g1 = ggplot(annual_spatial_synchrony_summary, aes(x = Mean_sync, y = reorder(Species, Mean_sync))) +
  geom_point(aes(shape = Crosses_Zero), color = "black", size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed") + 
  coord_cartesian()+
  ggstatsplot::theme_ggstatsplot() +
  labs(x = "Mean Spatial Synchrony", y = "") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "italic")) +
  ggeasy::easy_all_text_size(12) +
  scale_shape_manual(values = c("does_not_cross_0" = 16,
                                "crosses_0" = 1)) +  
  ggeasy::easy_remove_legend()

ggsave("Plots/Forest_plot_mean_synchrony.tiff", plot = g1, width = 5, height = 7, dpi = 300)

# 2. Plot MRM coefficients ----
regression_results$significant = ifelse(regression_results$f_stat <= 0.001, 
                                        "signif", 
                                        "does_not_signif")

# Columns dd/disp_diff were relabelled in 3_ANALYSIS.R; each panel shows the same coefficient as before.
g2 = ggplot(regression_results, aes(x = dd, y = reorder(Specie, annual_spatial_synchrony_summary$Mean_sync))) +
geom_point(shape=21, size = 3, color = "black", fill = "white",show.legend = TRUE) + geom_errorbarh(aes(xmin = dd - se_dd, xmax = dd + se_dd), color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_ggstatsplot() +
  easy_all_text_size(15) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "italic"))+
  ylab(NULL)+
  easy_text_size(which = "axis.title.x", size = 20)+
  easy_text_size(which = "axis.text.x", size = 20)+
  xlab("DD dissimilarity")+
  easy_remove_legend()+
  easy_y_axis_labels_size(20)

g3 = ggplot(regression_results, aes(x = disp_diff, y = reorder(Specie, annual_spatial_synchrony_summary$Mean_sync))) +
  geom_point(shape=21, size = 3, color = "black", fill = "white",show.legend = TRUE) + geom_errorbarh(aes(xmin = disp_diff - se_disp_diff, xmax = disp_diff + se_disp_diff), color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_ggstatsplot() +
  easy_all_text_size(15) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "italic"))+easy_remove_y_axis()+
  easy_text_size(which = "axis.title.x", size = 20)+
  easy_text_size(which = "axis.text.x", size = 20)+
  xlab("Dispersal difficulty")+
  easy_remove_legend()

g4 = ggplot(regression_results, aes(x = env ,y=reorder(Specie, annual_spatial_synchrony_summary$Mean_sync))) +
  geom_point(shape=21, size = 3, color = "black", fill = "white",show.legend = TRUE) + geom_errorbarh(aes(xmin = env - se_env, xmax = env + se_env), color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_ggstatsplot() +
  easy_all_text_size(15) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "italic"))+easy_remove_y_axis()+
  easy_text_size(which = "axis.title.x", size = 20)+
  easy_text_size(which = "axis.text.x", size = 20)+
  xlab("Water level synchrony")+
  easy_remove_legend()

g5 = ggplot(regression_results, aes(x = r_2 ,y=reorder(Specie, annual_spatial_synchrony_summary$Mean_sync))) +
  geom_point(aes(shape= significant), size = 3, color = "black", show.legend = TRUE)+
  theme_ggstatsplot() +
  easy_all_text_size(15) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = "italic"))+easy_remove_y_axis()+
  easy_text_size(which = "axis.title.x", size = 20)+
  easy_text_size(which = "axis.text.x", size = 20)+
  xlab("R²")+
  scale_shape_manual(values = c("signif" = 16,  "does_not_signif" = 1))+
  easy_remove_legend()+
  scale_x_continuous(
    name = "R²",
    breaks = seq(0, 0.3, by = 0.1),      
    limits = c(0, 0.3),                  
    labels = scales::label_number()    
  )

fig_mrm = g2 | g3 | g4 | g5 + 
  plot_layout(widths = c(2, 1, 1, -1.5))

ggsave("Plots/Forest_plot_coefficients.tiff", plot = fig_mrm, width = 17, height = 10, dpi = 300)

# 3. Plot PCOA and model residuals ----

g6 = sp_trait |>
  ggplot(aes(x = PC1, y = PC2, shape= Traits)) +
  geom_point(size = 3.8) +
  scale_shape_manual(values = c(15, 17, 16, 18)) +
  xlab(paste("Axis 1 (41%)")) +
  ylab(paste("Axis 2: (21%)")) +
  theme_bw()+
  easy_move_legend(to="left")+
  theme(
    legend.justification = c("top")
  )+
  easy_remove_legend_title()+
  easy_all_text_size(20)

g7 = sp_trait |>
  ggplot(aes(x = PC1, y = PC2, shape= Diet)) +
  geom_point(size = 3.8) +
  scale_shape_manual(values = c(15, 17, 16, 18, 0, 1, 2)) +
  xlab(paste("Axis 1: (41%)"))+
  ylab(paste("Axis 2: (21%)")) +
  theme_bw()+
  easy_move_legend(to="left")+
  theme(
    legend.justification = c("top")
  )+
  easy_remove_legend_title()+
  easy_all_text_size(20)

g8 = (g6 + g7) 
ggsave("Plots/PCoA_plots.tiff", plot = g8, width = 15, height = 5,dpi = 600)


png("Plots/Supplementary_LinearModel_Diagnostics.png",
    width = 2400,
    height = 2400,
    res = 300)
par(
  mfrow = c(2, 2),
  mar = c(5, 5, 4, 2) + 0.1, 
  cex.lab = 1.4,             # Axis titles (X and Y labels)
  cex.axis = 1.2,            # Tick mark labels
  cex.main = 1.5,            # Main plot titles
  cex.sub = 1.2              # Subtitles / captions (e.g., leverage notes)
)
plot(linear_model)
dev.off()
par(mfrow = c(1, 1))

# 4. Timescale decomposition plots ----

# Define scale parameters
ano_inicial = 2013 + (5 / 12)

# Transform time and timescale
res_list_corrigido = lapply(res_list, function(res_wpmf) {
  res_wpmf$times = ((res_wpmf$times - 1) / 4) + ano_inicial
  res_wpmf$timescales = res_wpmf$timescales / 4
  return(res_wpmf)
})

pdf("Plots/all_wavelets_grouped.pdf", width = 12, height = 11)
par(
  mfrow = c(4, 3), 
  mar = c(4.2, 3.8, 2.5, 4.5), 
  oma = c(1, 1, 1, 1), 
  mgp = c(2.5, 1, 0),
  cex.axis = 2,  
  cex.lab = 1.9    
)
grouped_indices = c(9, 11, 7, 1, 4, 5, 6, 8, 10, 3, 2)
panel_labels = paste0("(", tolower(LETTERS[1:length(grouped_indices)]), ")")
for(i in seq_along(grouped_indices)) {
  idx = grouped_indices[i]
  is_left_col = (i %% 3 == 1)          
  is_bottom   = (i %in% c(9, 10, 11)) 
  has_legend  = (i %in% c(3, 6, 9, 11))
  
  par(ann = FALSE) 
  if (is_left_col) {
    plotmag(res_list_corrigido[[idx]], zlim = c(0, 1), colorbar = has_legend)
  } else {
    plotmag(res_list_corrigido[[idx]], zlim = c(0, 1), yaxt = "n", colorbar = has_legend)
  }
  par(ann = TRUE)  
  if (is_left_col) {
    title(ylab = "Timescale (Years)") 
  }
  if (is_bottom) {
    title(xlab = "Time (Years)")
  }
  mtext(panel_labels[i], side = 3, line = 0.5, adj = 0, font = 2, cex = 1.5)
}
dev.off()

# # Species identification 
# # Extract species names from the list and clean formatting
# raw_names = names(res_list_corrigido)[grouped_indices]
# species_names = gsub("\\.", " ", raw_names)
# # Create a clear mapping dataframe
# species_legend = data.frame(
#   Panel = panel_labels,
#   Index = grouped_indices,
#   Species = species_names,
#   stringsAsFactors = FALSE
# )
#write.csv(species_legend, "Plots/wavelet_panels_legend.csv", row.names = FALSE)

# Synchrony Networks ----

# 1. Load Shapefile & Prep Base Map ---
madeira_river = readRDS("Data/mad_shape.rds")

sampling_map = ggplot() +
  geom_sf(data = madeira_river, fill = NA, color = "#4A90E2", alpha = 0.6) +
  geom_point(
    data = coord,
    aes(x = lon, y = lat),
    size = 3, shape = 21, color = "black"
  ) +
  geom_text_repel(
    data = coord,
    aes(x = lon, y = lat, label = Locals),
    size = 3.5, fontface = "bold",
    box.padding = 0.5, point.padding = 0.5,
    min.segment.length = 0, segment.color = "grey30",
    max.overlaps = Inf
  ) +
  coord_sf(xlim = c(-65.5, -64.5), ylim = c(-9.8, -9), expand = FALSE) +
  theme_classic() +
  labs(title = "Sampling Locations", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 14))

# Initialize plot container with the base map
clust_plots = list(sampling_map)

# 2. Generate Network Plots for Species with >1 Module ---
sp_names = names(annual_species_matrices)

for (i in seq_along(modularity_results)) {
  # Check number of unique clusters/modules
  curr_clusters = cluster_results[[i]]
  n_modules = length(unique(curr_clusters[!is.na(curr_clusters)]))
  # Filter: Only plot if more than 1 module was identified
  if (n_modules > 1) {
    sp_name = if (!is.null(sp_names)) gsub("[._]", " ", sp_names[i]) else paste("Species", i)
    # Copy coordinates and apply reproducible jitter
    set.seed(123)
    coord_curr = coord %>%
      mutate(
        jit_lon = jitter(lon, amount = 0.01),
        jit_lat = jitter(lat, amount = 0.01),
        module  = as.factor(curr_clusters),
        weight  = modularity_results[[i]]$nodeQ
      )
    
    g = ggplot() +
      geom_sf(data = madeira_river, fill = NA, color = "#4A90E2", alpha = 0.6) +
      geom_point(
        data = coord_curr,
        aes(x = jit_lon, y = jit_lat, fill = module, size = weight),
        shape = 21
      ) +
      ggeasy::easy_all_text_size(20)+
      scale_fill_brewer(palette = "Set2", name = "Module") +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      coord_sf(xlim = c(-65.5, -64.5), ylim = c(-9.8, -9), expand = FALSE) +
      theme_classic() +
      labs(title = sp_name, x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(face = "italic", size = 14))
    clust_plots[[length(clust_plots) + 1]] = g
  }
}

# 3. Export PDF Pages ---
plots_per_page = 6
num_plots = length(clust_plots)

if (num_plots > 0) {
  plot_groups = split(seq_len(num_plots), ceiling(seq_len(num_plots) / plots_per_page))
  
  pdf("Plots/Synchrony_networks.pdf", width = 16, height = 12)
  for (i in seq_along(plot_groups)) {
    indices = plot_groups[[i]]
    sub_list = clust_plots[indices]
    current_tags = as.character(indices)
    page_plot = wrap_plots(sub_list, ncol = 3, nrow = 2) +
      plot_annotation(tag_levels = list(current_tags)) & 
      theme(legend.position = "right")
    print(page_plot)
  }
  dev.off()
}


