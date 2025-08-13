setwd("D:/onedrive1/OneDrive - mails.ucas.ac.cn/Zhanglab/wmx/Revio/20250724_PTA/data_process")
library("ggplot2")
library("tidyr")
library("dplyr")



### 2025.6.24
#### all chimeric rate
XT12_E12_01_PTA_chimera_sup_df <- read.table("XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.summary.chimericrate.sup.tsv", header = F, sep = "\t")
XT12_E12_01_PTA_chimera_sup_df$V6 <- gsub("\\.", "\n", sub("\\.[^.]*$", "", XT12_E12_01_PTA_chimera_sup_df$V1))
XT12_E12_01_PTA_chimera_sup_df$V6 <- factor(XT12_E12_01_PTA_chimera_sup_df$V6, levels = c("XT1_0", "XT1_1", "XT2_0", "XT2_1", "E1_0", "E1_1", "E2_0", "E2_1"))
XT12_E12_01_PTA_chimera_sup_df$V7 <- XT12_E12_01_PTA_chimera_sup_df$V6
XT12_E12_01_PTA_chimera_sup_df$V8 <- "PTA"
XT12_E12_01_PTA_chimera_sup_df[grepl("E1_0", XT12_E12_01_PTA_chimera_sup_df$V1), ]$V8 <- "MDA" 
XT12_E12_01_PTA_chimera_sup_df[grepl("E2_0", XT12_E12_01_PTA_chimera_sup_df$V1), ]$V8 <- "MDA" 
XT12_E12_01_PTA_chimera_sup_df[grepl("XT1_0", XT12_E12_01_PTA_chimera_sup_df$V1), ]$V8 <- "MDA" 
XT12_E12_01_PTA_chimera_sup_df[grepl("XT2_0", XT12_E12_01_PTA_chimera_sup_df$V1), ]$V8 <- "MDA" 
XT12_E12_01_PTA_chimera_sup_df


p_XT12_E12_01_chimera <- ggplot(data = XT12_E12_01_PTA_chimera_sup_df, mapping = aes(x = V6, y = V5, fill = V8, label = paste(round(V5 * 100, 2), "%", sep = ""))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(size = 8) +
  # facet_grid( ~ V8, scales = "free_x", space = "free_x") +
  labs(x = "samples", y = "chimeric rate", title = "XT_12 E_12 1/0 chimeric rate (supplementary mapping)") +
  theme_bw() +
  theme(legend.title = element_blank())
  # theme(legend.position = "none")
# theme(axis.text.x = element_text(vjust = jitter(x=0.1)))
p_XT12_E12_01_chimera
ggsave(plot = p_XT12_E12_01_chimera, filename = "XT12_E12_01.chimeric_rate.sup.pdf", width = 12, height = 8)
ggsave(plot = p_XT12_E12_01_chimera, filename = "XT12_E12_01.chimeric_rate.sup.png", width = 12, height = 8)



##### inverted repeat
XT12_E12_01_PTA_inverted_repeat_df <- read.table("XT1_XT2_1_0.E1_E2_1_0.PTA.ccs.inverted_repeat.count", header = F, sep = "\t")
XT12_E12_01_PTA_inverted_repeat_df$V5 <- XT12_E12_01_PTA_inverted_repeat_df$V2 / XT12_E12_01_PTA_inverted_repeat_df$V3
XT12_E12_01_PTA_inverted_repeat_df$V6 <- XT12_E12_01_PTA_inverted_repeat_df$V2 / XT12_E12_01_PTA_inverted_repeat_df$V4
XT12_E12_01_PTA_inverted_repeat_df$V7 <- XT12_E12_01_PTA_inverted_repeat_df$V3 / XT12_E12_01_PTA_inverted_repeat_df$V4

XT12_E12_01_PTA_inverted_repeat_df
colnames(XT12_E12_01_PTA_inverted_repeat_df) <- c("name", "inverted_repeat_count", "all_chimera", "total_reads", 
                                                   "inverted_repeat_count/all_chimera", "inverted_repeat_count/total_reads", "all_chimera/total_reads")
XT12_E12_01_PTA_inverted_repeat_df$name1 <- XT12_E12_01_PTA_inverted_repeat_df$name
XT12_E12_01_PTA_inverted_repeat_df$Platform <- "PTA"
XT12_E12_01_PTA_inverted_repeat_df[grepl("E1_0", XT12_E12_01_PTA_inverted_repeat_df$name), ]$Platform <- "MDA" 
XT12_E12_01_PTA_inverted_repeat_df[grepl("E2_0", XT12_E12_01_PTA_inverted_repeat_df$name), ]$Platform <- "MDA" 
XT12_E12_01_PTA_inverted_repeat_df[grepl("XT1_0", XT12_E12_01_PTA_inverted_repeat_df$name), ]$Platform <- "MDA" 
XT12_E12_01_PTA_inverted_repeat_df[grepl("XT2_0", XT12_E12_01_PTA_inverted_repeat_df$name), ]$Platform <- "MDA" 
XT12_E12_01_PTA_inverted_repeat_df


XT12_E12_01_PTA_inverted_repeat_plot_df <- XT12_E12_01_PTA_inverted_repeat_df %>%
  mutate(
    all_chimera_ratio = all_chimera / total_reads,
    inverted_repeat_ratio = inverted_repeat_count / total_reads
  ) %>%
  select(name1, Platform, all_chimera_ratio, inverted_repeat_ratio) %>%
  pivot_longer(cols = c(all_chimera_ratio, inverted_repeat_ratio),
               names_to = "Metric",
               values_to = "Value")
XT12_E12_01_PTA_inverted_repeat_plot_df$Percentage <- paste(round(XT12_E12_01_PTA_inverted_repeat_plot_df$Value * 100, 2), "%", sep = "")
XT12_E12_01_PTA_inverted_repeat_plot_df$name1 <- factor(XT12_E12_01_PTA_inverted_repeat_plot_df$name1, levels = c("XT1_0", "XT1_1", "XT2_0", "XT2_1", "E1_0", "E1_1", "E2_0", "E2_1"))
XT12_E12_01_PTA_inverted_repeat_plot_df

p_XT12_E12_01_PTA_inverted_repeat_plot <- ggplot(XT12_E12_01_PTA_inverted_repeat_plot_df, aes(x = name1, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(size = 8, mapping = aes(label = Percentage), stat = "identity", position = position_dodge(width = 0.8)) +
  # facet_grid( .~ Platform, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Ratio", title = "XT_12 E_12 1/0 inverted repeat ratio") +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(1, 1),             # top-right inside
    legend.justification = c(1, 1), 
    legend.background = element_blank()
  )
p_XT12_E12_01_PTA_inverted_repeat_plot
ggsave(plot = p_XT12_E12_01_PTA_inverted_repeat_plot, filename = "XT12_E12_01.inverted_repeat_ratio.pdf", width = 16, height = 8)
ggsave(plot = p_XT12_E12_01_PTA_inverted_repeat_plot, filename = "XT12_E12_01.inverted_repeat_ratio.png", width = 16, height = 8)


XT12_E12_01_PTA_inverted_repeat_df$Percentage <- paste(round(XT12_E12_01_PTA_inverted_repeat_df$'inverted_repeat_count/all_chimera' * 100, 2), "%", sep = "")
XT12_E12_01_PTA_inverted_repeat_df$name1 <- factor(XT12_E12_01_PTA_inverted_repeat_df$name1, levels = c("XT1_0", "XT1_1", "XT2_0", "XT2_1", "E1_0", "E1_1", "E2_0", "E2_1"))
p_XT12_E12_01_PTA_inverted_repeat_plot1 <- ggplot(XT12_E12_01_PTA_inverted_repeat_df, aes(x = name1, y = inverted_repeat_count/all_chimera, fill = Platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(size = 8, mapping = aes(label = Percentage), stat = "identity", position = position_dodge(width = 0.8)) +
  # facet_grid( .~ Platform, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Ratio", title = "XT_12 E_12 1/0 inverted_repeat/all_chimeria ratio") +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  # theme(legend.title = element_blank())
  theme(legend.title = element_blank())
p_XT12_E12_01_PTA_inverted_repeat_plot1
ggsave(plot = p_XT12_E12_01_PTA_inverted_repeat_plot1, filename = "XT12_E12_01.inverted_repeat_ratio1.pdf", width = 16, height = 8)
ggsave(plot = p_XT12_E12_01_PTA_inverted_repeat_plot1, filename = "XT12_E12_01.inverted_repeat_ratio1.png", width = 16, height = 8)



