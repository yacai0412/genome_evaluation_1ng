setwd("D:/onedrive1/OneDrive - mails.ucas.ac.cn/Zhanglab/wmx/Revio/20250724_PTA/data_process")
library("ggplot2")
library("patchwork")



### 2025.7.24
XT1_0_relative_depth_df <- read.table("XT1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT1_1_relative_depth_df <- read.table("XT1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_0_relative_depth_df <- read.table("XT2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_1_relative_depth_df <- read.table("XT2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_0_relative_depth_df <- read.table("E1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_1_relative_depth_df <- read.table("E1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_0_relative_depth_df <- read.table("E2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_1_relative_depth_df <- read.table("E2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")


plot_relative_depth <- function(df, title = "Relative Depth Plot") {
  color_A <- c("chr2L", "chr3L", "chr4", "chrY")
  color_B <- c("chr2R", "chr3R", "chrX")
  chrom_levels <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
  
  # Ensure proper column names and add color column
  df$color <- "NA"
  names(df) <- c("chromosome", "position", "depth", "pos", "color")
  df$color[df$chromosome %in% color_A] <- "A"
  df$color[df$chromosome %in% color_B] <- "B"
  
  # Set factor levels
  df$chromosome <- factor(df$chromosome, levels = chrom_levels)
  
  # Create the plot
  p <- ggplot(data = df, aes(x = pos, y = depth, color = color)) +
    geom_point(shape = 1, size = 1.5) +
    facet_grid(. ~ chromosome, space = "free_x", scales = "free_x") +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(limits = c(0, 4)) +
    labs(x = "", y = "relative depth", title = title) + 
    theme_bw() + 
    theme(panel.border = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10, angle = 90))
  
  return(p)
}


p_XT1_0 <- plot_relative_depth(XT1_0_relative_depth_df, title = "XT1_0 relative depth (bin size: 100 kb)")
print(p_XT1_0)
p_XT1_1 <- plot_relative_depth(XT1_1_relative_depth_df, title = "XT1_1 relative depth (bin size: 100 kb)")
print(p_XT1_1)
p_XT2_0 <- plot_relative_depth(XT2_0_relative_depth_df, title = "XT2_0 relative depth (bin size: 100 kb)")
print(p_XT2_0)
p_XT2_1 <- plot_relative_depth(XT2_1_relative_depth_df, title = "XT2_1 relative depth (bin size: 100 kb)")
print(p_XT2_1)
p_E1_0 <- plot_relative_depth(E1_0_relative_depth_df, title = "E1_0 relative depth (bin size: 100 kb)")
print(p_E1_0)
p_E1_1 <- plot_relative_depth(E1_1_relative_depth_df, title = "E1_1 relative depth (bin size: 100 kb)")
print(p_E1_1)
p_E2_0 <- plot_relative_depth(E2_0_relative_depth_df, title = "E2_0 relative depth (bin size: 100 kb)")
print(p_E2_0)
p_E2_1 <- plot_relative_depth(E2_1_relative_depth_df, title = "E2_1 relative depth (bin size: 100 kb)")
print(p_E2_1)


p_XT1_XT2_1_0_relative_depth <- (p_XT1_0 / p_XT1_1 / p_XT2_0 / p_XT2_1)
p_XT1_XT2_1_0_relative_depth
ggsave(plot = p_XT1_XT2_1_0_relative_depth, filename = "XT1_XT2_1_0.relative_depth.pdf", width = 24, height = 16)
ggsave(plot = p_XT1_XT2_1_0_relative_depth, filename = "XT1_XT2_1_0.relative_depth.png", width = 24, height = 16)

p_E1_E2_1_0_relative_depth <- (p_E1_0 / p_E1_1 / p_E2_0 / p_E2_1)
p_E1_E2_1_0_relative_depth
ggsave(plot = p_E1_E2_1_0_relative_depth, filename = "E1_E2_1_0.relative_depth.pdf", width = 24, height = 16)
ggsave(plot = p_E1_E2_1_0_relative_depth, filename = "E1_E2_1_0.relative_depth.png", width = 24, height = 16)




#### cv
cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

calculate_cv_autosomes <- function(data_list, autosomes = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4")) {
  result_list <- list()
  
  for (sample_name in names(data_list)) {
    df <- data_list[[sample_name]]
    
    # Combine all autosomal bins
    autosome_depths <- df[df$V1 %in% autosomes, ]$V3
    
    # Calculate CV over combined autosomal depth
    cv_value <- cv(autosome_depths)
    
    result_list[[length(result_list) + 1]] <- data.frame(
      sample = sample_name,
      chromosome = "autosomes",
      cv = cv_value
    )
  }
  
  result_df <- do.call(rbind, result_list)
  return(result_df)
}


XT1_0_relative_depth_df <- read.table("XT1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT1_1_relative_depth_df <- read.table("XT1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_0_relative_depth_df <- read.table("XT2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_1_relative_depth_df <- read.table("XT2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_0_relative_depth_df <- read.table("E1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_1_relative_depth_df <- read.table("E1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_0_relative_depth_df <- read.table("E2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_1_relative_depth_df <- read.table("E2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")


XT12_E12_01_data_list <- list(
  XT1_0 = XT1_0_relative_depth_df,
  XT1_1 = XT1_1_relative_depth_df,
  XT2_0 = XT2_0_relative_depth_df,
  XT2_1 = XT2_1_relative_depth_df,
  E1_0 = E1_0_relative_depth_df,
  E1_1 = E1_1_relative_depth_df,
  E2_0 = E2_0_relative_depth_df,
  E2_1 = E2_1_relative_depth_df
)

XT12_E12_01_cv_df <- calculate_cv_autosomes(XT12_E12_01_data_list)
XT12_E12_01_cv_df$platform <- "PTA"
XT12_E12_01_cv_df[grepl("XT1_0", XT12_E12_01_cv_df$sample), ]$platform <- "MDA" 
XT12_E12_01_cv_df[grepl("XT2_0", XT12_E12_01_cv_df$sample), ]$platform <- "MDA" 
XT12_E12_01_cv_df[grepl("E1_0", XT12_E12_01_cv_df$sample), ]$platform <- "MDA" 
XT12_E12_01_cv_df[grepl("E2_0", XT12_E12_01_cv_df$sample), ]$platform <- "MDA" 
XT12_E12_01_cv_df$sample <- factor(XT12_E12_01_cv_df$sample, levels = c("XT1_0", "XT1_1", "XT2_0", "XT2_1", "E1_0", "E1_1", "E2_0", "E2_1"))


p_XT12_E12_01_autosomes_cv_results <- ggplot(data = XT12_E12_01_cv_df, aes(x = sample, y = cv, fill = platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(mapping = aes(label = round(cv, 3)), 
            position = position_dodge(width = 0.9)) +
  # facet_grid( ~ Platform, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Coefficient of Variation (CV)", title = "XT_12 E_12 1/0 CV of Autosomes Relative Depth") +
  theme_bw() +
  scale_fill_brewer(palette = "YlGnBu")
# theme(legend.position = "none")
p_XT12_E12_01_autosomes_cv_results
ggsave(plot = p_XT12_E12_01_autosomes_cv_results, filename = "XT12_E12_01.autosomes_cv_results.pdf", width = 10, height = 6)
ggsave(plot = p_XT12_E12_01_autosomes_cv_results, filename = "XT12_E12_01.autosomes_cv_results.png", width = 10, height = 6)




#### gini index
gini_index <- function(x) {
  if (all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  G <- sum((2 * (1:n) - n - 1) * x)
  G / (n * sum(x))
}

calculate_gini_autosomes <- function(data_list, autosomes = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4")) {
  result_list <- list()
  
  for (sample_name in names(data_list)) {
    df <- data_list[[sample_name]]
    
    # Combine all autosomal bins
    autosome_depths <- df[df$V1 %in% autosomes, ]$V3
    
    # Calculate Gini index over combined autosomal depth
    gini_value <- gini_index(autosome_depths)
    
    result_list[[length(result_list) + 1]] <- data.frame(
      sample = sample_name,
      chromosome = "autosomes",
      gini = gini_value
    )
  }
  
  result_df <- do.call(rbind, result_list)
  return(result_df)
}



XT1_0_relative_depth_df <- read.table("XT1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT1_1_relative_depth_df <- read.table("XT1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_0_relative_depth_df <- read.table("XT2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
XT2_1_relative_depth_df <- read.table("XT2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_0_relative_depth_df <- read.table("E1_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E1_1_relative_depth_df <- read.table("E1_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_0_relative_depth_df <- read.table("E2_0.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")
E2_1_relative_depth_df <- read.table("E2_1.ccs.fq.dm6.F4.s.bam.100000.relative.depth", header = F, sep = "\t")


XT12_E12_01_data_list <- list(
  XT1_0 = XT1_0_relative_depth_df,
  XT1_1 = XT1_1_relative_depth_df,
  XT2_0 = XT2_0_relative_depth_df,
  XT2_1 = XT2_1_relative_depth_df,
  E1_0 = E1_0_relative_depth_df,
  E1_1 = E1_1_relative_depth_df,
  E2_0 = E2_0_relative_depth_df,
  E2_1 = E2_1_relative_depth_df
)

XT12_E12_01_gini_df <- calculate_gini_autosomes(XT12_E12_01_data_list)
XT12_E12_01_gini_df$platform <- "PTA"
XT12_E12_01_gini_df[grepl("XT1_0", XT12_E12_01_gini_df$sample), ]$platform <- "MDA" 
XT12_E12_01_gini_df[grepl("XT2_0", XT12_E12_01_gini_df$sample), ]$platform <- "MDA" 
XT12_E12_01_gini_df[grepl("E1_0", XT12_E12_01_gini_df$sample), ]$platform <- "MDA" 
XT12_E12_01_gini_df[grepl("E2_0", XT12_E12_01_gini_df$sample), ]$platform <- "MDA" 
XT12_E12_01_gini_df$sample <- factor(XT12_E12_01_gini_df$sample, levels = c("XT1_0", "XT1_1", "XT2_0", "XT2_1", "E1_0", "E1_1", "E2_0", "E2_1"))


p_XT12_E12_01_autosomes_gini_results <- ggplot(data = XT12_E12_01_gini_df, aes(x = sample, y = gini, fill = platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(mapping = aes(label = round(gini, 3)), 
            position = position_dodge(width = 0.9)) +
  # facet_grid( ~ Platform, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Coefficient of Variation (gini)", title = "XT_12 E_12 1/0 Gini index of Autosomes Relative Depth") +
  theme_bw() +
  scale_fill_brewer(palette = "YlGnBu")
# theme(legend.position = "none")
p_XT12_E12_01_autosomes_gini_results
ggsave(plot = p_XT12_E12_01_autosomes_gini_results, filename = "XT12_E12_01.autosomes_gini_results.pdf", width = 10, height = 6)
ggsave(plot = p_XT12_E12_01_autosomes_gini_results, filename = "XT12_E12_01.autosomes_gini_results.png", width = 10, height = 6)





















