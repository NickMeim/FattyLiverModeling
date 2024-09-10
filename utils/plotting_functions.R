# Script with wrapper functions for ggplot with our specific style
# This is a nice tutorial to do similar things, probably better than what I have
# https://rpubs.com/hadley/97970

### Load libraries
library(RColorBrewer)
library(ggpubr)
library(tidyverse, quietly = T)
library(ggfortify, quietly = T)
library(ggrepel)
library(extrafont)
#font_import()
loadfonts(device = "win")

### Palettes and styles for plot elements
font_family = "sans"
size_dot = 1.1 # dot size for plots
size_error = 0.4 # size (thickness of line) for error bars
size_text = 8
size_col = 0.1 # Thickness of line for barplots
size_legend = 0.3 # size of legend boxes in cm
size_line = 0.6 # thickness of lines
size_annotation = 3 # size for asterisks and other annotations
size_stroke = 0.1 # Thickness of stroke for dot plots
color_palette = "Set2" # Palette to be used. Choose from the Color Brewer palettes


### Wrapper functions - these functions return a plot without any stats or annotations, since that part of the code is more customized
### depending on the analysis. These functions keep a consistent style and theme between plots; only the arguments (data) changes

# ==========================================================================
# Function 0 - Simply add style to an existing plot. Useful when
# interfacing with plots produced by other packages
add_theme <- function(plt){
  plt <- plt +
        theme_bw() +
        theme(text = element_text(size = size_text, family = font_family),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black"),
          legend.key.size = unit(size_legend, "cm"))
}


# ==========================================================================
# Function 1 - Plotting columns with error bars and individual reps as dots

plot_bars <- function(data_bar, mapping_bar, data_points = NULL, mapping_points = NULL, data_error = NULL,
                      mapping_error = NULL, plt_labs = NULL, col_pos = "dodge2",
                      palette = color_palette, extra_args = NULL, pos_points = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)){


  # Create plot - column, individual points, and error bar
  plt <- ggplot(data = data_bar, mapping = mapping_bar) +
          geom_col(position = col_pos, size = size_col, color = "black")
          
  
  if (!is.null(data_points)){
    plt <- plt + geom_point(data = data_points, mapping = mapping_points,
                            size = size_dot, show.legend = extra_args$show.legend.pts,
                            color = "#666666",
                            position = pos_points)
  }
  
  if (!is.null(data_error)){
    plt <- plt + geom_errorbar(data = data_error, mapping = mapping_error,
                               position = position_dodge(width=0.9), width = 0.25,
                               size = size_error, colour = "black")
  }

  # Choose between brewer pal or manual depending on the given palettes
  if (length(palette) > 1){
    fill_aes <- scale_fill_manual(values = palette, labels = plt_labs$fill_labs)
  }
  else {
    fill_aes <- scale_fill_brewer(palette = palette, labels = plt_labs$fill_labs)
  }
  
  # If there is an extra argument called fill, we override the fill of the bars
  if (!is.null(extra_args$fill)){
    plt$layers[[1]]$aes_params$fill <- extra_args$fill
  }

  # Add tick marks
  if (is.null(plt_labs$x_labs)){
    scale_x <- NULL
  }
  else {
    scale_x <- scale_x_discrete(labels = plt_labs$x_labs)
  }

  # Customize plot - axis labels are part of the plt_labs argument (a list)
  plt <- plt +
          fill_aes +
          scale_x +
          labs(x = plt_labs$xaxis, y = plt_labs$yaxis, fill = plt_labs$fill)

  # Return
  plt <- add_theme(plt)
  return(plt)
}

# Function 2: Dummy function for making simple XY scatter plots (ggplot wrapper)

plot_scatter <- function(x, y, xlab = NULL, ylab = NULL, args_extra = NULL, aes_extra = NULL, xjitter = F, yjitter = F){
  # Defaults
  size_dot <- 1
  def_color <- "steelblue"
  # If no x, use indices
  if (is.null(x)){
    x <- 1:length(y)
  }
  # Define data frame - TO DO: check if what was passed is a dataframe or single values, which wouldn't map to aesthetics
  df <- data.frame(x = x, y = y)
  if (!is.null(args_extra)){
    df <- cbind(df, as.data.frame(args_extra))
  }
  
  # Define aesthetics
  aes_plt <- aes(x = x, y = y)
  if (!is.null(aes_extra)){
    aes_plt <- modifyList(aes_plt, aes_extra)
    # Plot
    plt <- ggplot(df, mapping = aes_plt) +
      geom_point() +
      labs(x = xlab, y = ylab)
  } else {
    # Plot defaults
    plt <- ggplot(df, mapping = aes_plt) +
      geom_point(size = size_dot, color = def_color, 
                 position = position_jitter(width = ifelse(xjitter, 0.1, 0),
                                            height = ifelse(yjitter, 0.1, 0))) +
      labs(x = xlab, y = ylab)
  }
  
  # Add theme
  plt <- add_theme(plt)
  
  return(plt)
}

# Function 3: Dummy function for running PCA and making plots

# Samples in columns, features in rows. Metadata has sample info to add aesthetics
# aesthetics MUST correspond to a column in metadata
plot_pca <- function(data, metadata = NULL, color_by = NULL, shape_by = NULL, 
                     scale_data = T, return_data = F, label_features = F){
  # Do PCA and extract variances
  pca <- prcomp(t(data), scale. = scale_data)
  pca_var <- round(100*pca$sdev^2/sum(pca$sdev^2),2)
  # Dummy dataframe with metadata
  df <- cbind(metadata, pca$x)
  # Make plot - scores
  plt <- ggplot(df, aes_string(x = "PC1", y = "PC2", 
                               color = color_by, shape = shape_by)) +
    geom_point(size = 1.5) +
    labs(x = paste0("PC1 (",pca_var[1],"%)"),
         y = paste0("PC2 (",pca_var[2],"%)"))
  
  # Make plot - loadings
  df <- data.frame(pca$rotation, feature = rownames(pca$rotation))
  sd_PC1 <- sd(df$PC1)
  sd_PC2 <- sd(df$PC2)
  df <- df %>% mutate(is.big = (abs(PC1) > 2*sd_PC1) | (abs(PC2) > 2*sd_PC2))
  
  # Override if we want tolabel all loadings
  if (label_features == T){
    plt_load <- ggplot(df, aes(x = PC1, y = PC2, label = feature)) +
      geom_point(size = 1.5) +
      geom_text_repel() 
  }
  else {
    plt_load <- ggplot(df, aes(x = PC1, y = PC2, alpha = is.big, label = feature)) +
      geom_point(size = 1.5) +
      geom_vline(xintercept = 2*sd_PC1, color = "red", linetype = 2) +
      geom_vline(xintercept = -2*sd_PC1, color = "red", linetype = 2) +
      geom_hline(yintercept = 2*sd_PC2, color = "red", linetype = 2) +
      geom_hline(yintercept = -2*sd_PC2, color = "red", linetype = 2) +
      geom_text_repel(data = df %>% filter(is.big == T))
    
  }
  plt_load <- plt_load + 
    labs(x = paste0("PC1 (",pca_var[1],"%)"),
         y = paste0("PC2 (",pca_var[2],"%)")) +
    theme(legend.position = "none")
  
  
  if (return_data){
    return(list(data = pca, plt_score = plt, plt_load = plt_load))
  }
  else {
    return(list(plt_score = plt, plt_load = plt_load))
  }
}

# Function 4: Make scatter plots with side boxplots
scatter_box_plot <- function(df,legend_title,
                             font_size =20,font_family = 'Arial',
                             point_shape = 21,point_size=2.8,point_stroke=1.2,
                             x_axis='LV extra 1',y_axis='LV extra 2',
                             box_y_width = 0.2,jitter_y_width=0.1,
                             jitter_x_height = 0.2,
                             theme_use = 'minimal',
                             plotting=TRUE){
  if (theme_use=='bw'){
    th <- ggplot2::theme_bw(base_size = font_size, base_family = font_family)
  }else{
    th <- ggplot2::theme_minimal(base_size = font_size, base_family = font_family)
  }
  scatter_plot <- ggplot(df, aes(x = V1, y = V2, fill = pheno)) +
    geom_point(size = point_size, shape = point_shape, stroke = point_stroke) +
    scale_fill_viridis_c() +
    xlab(x_axis) +
    ylab(y_axis) +
    labs(fill=legend_title)+
    th +
    theme(
      text = element_text(size = font_size, family = font_family),
      legend.position = 'top',
    )
  
  # Boxplot for V1 (x-axis)
  boxplot_x <- ggplot(df, aes(x = V1, y = "", color = pheno)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(height = jitter_x_height, size = point_size*0.5)+
    scale_color_viridis_c() +
    theme_minimal(base_size = font_size, base_family = font_family) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    )
  
  # Boxplot for V2 (y-axis)
  boxplot_y <- ggplot(df, aes(x = "", y = V2, color = pheno)) +
    geom_boxplot(outliers = FALSE,width = box_y_width) +
    geom_jitter(width=jitter_y_width, size = point_size*0.5)+
    scale_color_viridis_c() +
    theme_minimal(base_size = font_size, base_family = font_family) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
    )
  # Combine plots for NAS
  combo1 <- scatter_plot + boxplot_y + plot_layout(widths = c(3, 1))
  combo2 <- boxplot_x + plot_spacer() + plot_layout(widths = c(3, 1))
  combined_plot <- (combo1) / (combo2) +
    plot_layout(heights = c(4, 1))
  if (plotting==TRUE){
    print(combined_plot)
  }
  return(combined_plot)
}

# Function 5: Bar plots for pathway activity. Assumes only one condition is there, but can
# facet externally

plot_pwy_activity <- function(pwy_act_table, plt_lim = NULL, show_fill_legend = F, offset_annot = 2, n.breaks = 8){
  # Set plot limits
  if (is.null(plt_lim)){
    plt_lim = max(abs(pwy_act_table$score)) + 1
    }
  
  # Annotate pvals with asterisks
  pwy_act_table <- pwy_act_table %>%
                    mutate(annot = ifelse(p_value <= 0.0001, "****",
                                          ifelse(p_value <= 0.001,"***",
                                                 ifelse(p_value<=0.01,"**",
                                                        ifelse(p_value<=0.05,'*',
                                                               ifelse(p_value<=0.1,'\u2022', #\u2219
                                                                        'ns')))))) %>%
                   mutate(offset = ifelse(annot == 'ns', offset_annot, 0.5*offset_annot))
  # Make plot
  plt <- pwy_act_table %>%
          ggplot(aes(x = score,y = reorder(Pathway, score),fill = score, label = annot)) + 
            geom_bar(stat = 'identity', show.legend = show_fill_legend, color = "black", size = size_col) +
            scale_fill_gradient2(low='darkblue',high = 'indianred',mid = 'whitesmoke',
                               midpoint = 0,limits = c(-plt_lim,plt_lim))+
            scale_x_continuous(n.breaks = n.breaks,limits = c(-plt_lim,plt_lim))+
            labs(y = "Pathway", x = "Activity score", fill = "Score") +
            geom_text(aes(x = ifelse(score < 0, score - offset, score + offset)),
                      color = 'black', size = size_annotation*0.7,
                      angle = 90, show.legend = FALSE, hjust = "center", vjust = "center") 
  
  return(plt)
}


# Function 6: Make scatter plots with side boxplots - cowplot version
scatter_box_plot_V2 <- function(df,legend_title,
                             font_size =20,font_family = 'Arial',
                             point_shape = 21,point_size=2.8,point_stroke=1.2,
                             x_axis='LV extra 1',y_axis='LV extra 2',
                             box_y_width = 0.2,jitter_y_width=0.1,
                             jitter_x_height = 0.2,
                             theme_use = 'minimal',
                             plotting=TRUE){

  scatter_plot <- ggplot(df, aes(x = V1, y = V2, fill = pheno)) +
    geom_point(size = point_size, shape = point_shape, stroke = point_stroke) +
    scale_fill_viridis_c() +
    xlab(x_axis) +
    ylab(y_axis) +
    labs(fill=legend_title)
  
  scatter_plot <- add_theme(scatter_plot) + theme(legend.position = "top")

  # Boxplot for V1 (x-axis)
  boxplot_x <- ggplot(df, aes(x = V1, y = "", color = pheno)) +
    geom_boxplot(outliers = FALSE, size = size_line, width = 1) +
    geom_jitter(height = jitter_x_height, size = point_size*0.3, show.legend = F)+
    scale_color_viridis_c() 
  
  boxplot_x <- add_theme(boxplot_x) +
                theme_void() + theme(legend.position = "none")
                
  # Boxplot for V2 (y-axis)
  boxplot_y <- ggplot(df, aes(x = "", y = V2, color = pheno)) +
    geom_boxplot(outliers = FALSE,width = box_y_width, size = size_line, width = 1) +
    geom_jitter(width=jitter_y_width, size = point_size*0.3, show.legend = F)+
    scale_color_viridis_c() 
  
  boxplot_y <- add_theme(boxplot_y) +
                theme_void() + theme(legend.position = "none")

  # Combine plots for NAS
  combo1 <- scatter_plot + boxplot_y + plot_layout(widths = c(3, 1))
  combo2 <- boxplot_x + plot_spacer() + plot_layout(widths = c(3, 1))
  combined_plot <- (combo1) / (combo2) +
    plot_layout(heights = c(4, 1))
  if (plotting==TRUE){
    print(combined_plot)
  }
  return(combined_plot)
}