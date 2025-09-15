###########Fig2#################
##Fig2A alpha diversity
library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 读取细菌丰度表
abund <- read.table("Bac_exp_Fig2A.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
abund_t <- t(abund)

# 读取 metadata
metadata <- read.table("metadata_Fig2A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- metadata[match(rownames(abund_t), metadata$Sample), ]
rownames(metadata) <- metadata$Sample

# 确保顺序一致
stopifnot(all(rownames(abund_t) == rownames(metadata)))

# 计算 Shannon 指数
shannon <- diversity(abund_t, index = "shannon")
alpha_df <- data.frame(Sample = names(shannon), shannon = shannon)
alpha_df <- merge(alpha_df, metadata, by = "Sample")

# 设置 Als 分组顺序
alpha_df$Als <- factor(alpha_df$Als, 
                       levels = c("Rhizosphere", "Rhizome", "Bulk soil"),
                       ordered = TRUE)

# 两两比较组
my_comparisons <- list(c("Bulk soil","Rhizome"), c("Rhizosphere", "Rhizome"))

# 绘制 Shannon 指数小提琴图（细菌）
p1 <- ggviolin(alpha_df, x = "Als", y = "shannon", fill = "Als",
               add = "boxplot",
               add.params = list(fill = "white", alpha = 0.7),
               palette = "jco", 
               trim = FALSE, 
               size = 1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 读取真菌丰度表
abundf <- read.table("Fun_exp_Fig2A.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
abundf_t <- t(abundf)

# 读取 metadata
metadata <- read.table("metadata_Fig2A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- metadata[match(rownames(abundf_t), metadata$Sample), ]
rownames(metadata) <- metadata$Sample

# 确保顺序一致
stopifnot(all(rownames(abundf_t) == rownames(metadata)))

# 计算 Shannon 指数
shannon <- diversity(abundf_t, index = "shannon")
alpha_df <- data.frame(Sample = names(shannon), shannon = shannon)
alpha_df <- merge(alpha_df, metadata, by = "Sample")

# 设置 Als 分组顺序
alpha_df$Als <- factor(alpha_df$Als, 
                       levels = c("Rhizosphere", "Rhizome", "Bulk soil"),
                       ordered = TRUE)

# 绘制 Shannon 指数小提琴图（真菌）
p2 <- ggviolin(alpha_df, x = "Als", y = "shannon", fill = "Als",
               add = "boxplot",
               add.params = list(fill = "white", alpha = 0.7),
               palette = "jco", 
               trim = FALSE, 
               size = 1) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

## Fig2B
library(vegan)
library(ggplot2)
library(phyloseq)
library(reshape2)

# 读取微生物丰度表
merged_df <- read.table("micro_exp_Fig2B.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data_t <- t(merged_df)

# 匹配 metadata
metadata <- metadata[match(rownames(data_t), metadata$Sample), ]

# 组合分组信息
metadata$Als_Chemotype <- paste(metadata$Als, metadata$Chemotype, sep = "_")
metadata$chandi_wildcul <- paste(metadata$chandi, metadata$Wild.Cultivation, sep = "_")

# Bray-Curtis 距离
bc_dist <- vegdist(data_t, method = "bray")

# PCoA 分析
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
points <- as.data.frame(pcoa_res$points)
colnames(points) <- c("PCoA1", "PCoA2")
points$Sample <- rownames(points)

# 主成分解释度
eigen_values <- pcoa_res$eig
var_explained <- eigen_values[1:2] / sum(eigen_values[eigen_values > 0]) * 100

# 合并 metadata
plot_df <- merge(points, metadata, by = "Sample")

library(ggplot2)
library(RColorBrewer)

# 绘制 PCoA 图
ggplot(plot_df, aes(x = PCoA1, y = PCoA2)) +
  stat_ellipse(
    aes(group = Als_Chemotype, color = Als_Chemotype),
    type = "t",
    linewidth = 0.6,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      shape = Wild.Cultivation,
      fill = chandi,
      color = Als_Chemotype
    ),
    size = 5,
    stroke = 0.5,
    alpha = 0.8
  ) +
  scale_shape_manual(values = c("W" = 21, "P" = 24)) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Set2") +
  guides(
    color = guide_legend(override.aes = list(fill = NA)),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Diversity (Bray-Curtis, PCoA)",
    x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_explained[2]),
    shape = "Wild / Cultivated",
    fill = "chandi",
    color = "Als_Chemotype"
  ) +
  theme(
    panel.grid = element_line(color = "grey90"),
    legend.position = "right",
    legend.key = element_rect(fill = "white")
  )
##Fig2C
library(tidyverse)   # dplyr + tidyr + ggplot2
library(RColorBrewer)
library(scales)      # 百分比坐标轴

# 转置数据并添加 Sample 列
merged_df_t <- as.data.frame(t(merged_df))
merged_df_t$Sample <- rownames(merged_df_t)

# 合并 metadata
data_merged <- merged_df_t %>% 
  left_join(metadata[, c("Sample", "Als", "Chemotype")], by = "Sample")

# 提取物种列并转为数值型
taxa_cols <- setdiff(colnames(data_merged), c("Sample", "Als", "Chemotype"))
data_merged[taxa_cols] <- lapply(data_merged[taxa_cols], function(x) as.numeric(as.vector(x)))

# 计算每组 (Als × Chemotype) 的均值
df_summary <- data_merged %>% 
  group_by(Als, Chemotype) %>% 
  summarise(across(all_of(taxa_cols), mean, na.rm = TRUE), .groups = "drop")

# 计算总体平均丰度，选取前 11 个属
top11_taxa <- df_summary %>% 
  summarise(across(all_of(taxa_cols), mean, na.rm = TRUE)) %>% 
  pivot_longer(cols = everything(),
               names_to  = "Genus",
               values_to = "MeanAbun") %>% 
  arrange(desc(MeanAbun)) %>% 
  slice_head(n = 11) %>% 
  pull(Genus)

# 宽转长，非前 11 属归类为 Other
df_summary_long <- df_summary %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  mutate(
    Genus = if_else(Genus %in% top11_taxa, Genus, "Other")
  )

# 构建复合分组变量 Als+Chemotype
df_summary_long_alt <- df_summary_long %>%
  mutate(
    als_chemotype = paste(Als, Chemotype, sep = "_"),
    als_chemotype = factor(
      als_chemotype,
      levels = df_summary_long %>%
        distinct(Als, Chemotype) %>%
        arrange(Als, Chemotype) %>%
        mutate(combo = paste(Als, Chemotype, sep = "_")) %>%
        pull(combo)
    ),
    Genus = fct_relevel(Genus, "Other", after = Inf)
  )

# 绘制堆积柱状图
p <- ggplot(df_summary_long_alt,
            aes(x = als_chemotype, y = Abundance, fill = Genus)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    x = "Als + Chemotype",
    y = "Relative Abundance (%)",
    fill = "Genus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )
##Fig2D was drawn online http://bioinformatics.psb.ugent.be/webtools/Venn/
##Fig2E
# 读取数据
df <- read_table("lefse_Fig2E.txt")
#head(df)
# 筛选 group_genji == group_soil 的行
df_filtered <- df |>
  filter(group_genji == group_soil)

# 整理数据
df_long <- df_filtered |>
  pivot_longer(
    cols = c(LDA_genji, LDA_soil),
    names_to = "source",
    values_to = "value"
  ) |>
  mutate(
    group = if_else(source == "LDA_genji", group_genji, group_soil),
    source_group = paste0(if_else(source == "LDA_genji", "genji", "soil"), "_", group),
    value = if_else(group == "MSA", value, -value)
  )

# 自定义颜色
color_map <- c(
  "genji_HBA" = "#4e8872",
  "soil_HBA"  = "#79a019",
  "genji_MSA" = "#884c45",
  "soil_MSA"  = "#f4a638"
)

# 绘图
ggplot(df_long, aes(
  y = fct_reorder(genus, value),
  x = value,
  fill = source_group
)) +
  geom_col(
    position = position_dodge2(width = 0.95, padding = 0.3),
    width = 0.85
  ) +
  # 添加 x = 0 的虚线
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  scale_fill_manual(values = color_map) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Group Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    # 去掉 y 轴线、刻度、标签
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    # 保留 x 轴刻度和线
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    axis.line.x = element_line(),
    
    # 其他外观
    legend.position = "top",
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

write.table(df_filtered, file = "common_lefse.tsv", sep = "\t", quote = FALSE, col.names = F)
##Fig2F
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
#绘制Rhizome时只需要将Rhizosphere关键字替换即可
# 1. 筛选 Rhizosphere 样本及 Chemotype
meta_rhizo <- metadata %>% filter(Als == "Rhizosphere")
samples_rhizo <- meta_rhizo$Sample
chemotype_info <- meta_rhizo %>% select(Sample, Chemotype)

# 2. 属丰度矩阵 + presence/absence
abund_df <- merged_df[, samples_rhizo]
pa_df <- abund_df > 0

# 3. 核心属筛选
core_taxa_list <- chemotype_info %>%
  group_split(Chemotype) %>%
  setNames(unique(chemotype_info$Chemotype)) %>%
  map(function(group_df) {
    group_samples <- group_df$Sample
    pa_sub <- pa_df[, group_samples, drop = FALSE]
    rowSums(pa_sub) == length(group_samples)
  })

core_taxa_names <- map(core_taxa_list, ~ rownames(pa_df)[.x])

# 保存所有核心属
core_taxa_all <- map2(core_taxa_names, names(core_taxa_names), ~ tibble(Genus = .x, Chemotype = .y)) %>% bind_rows()
write.table(core_taxa_all, "core_taxa_all_Rhizosphere.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. 每组核心属平均丰度
core_taxa_abund <- map2(core_taxa_names, names(core_taxa_names), function(taxa, chemotype) {
  samples <- chemotype_info %>% filter(Chemotype == chemotype) %>% pull(Sample)
  sub_df <- abund_df[taxa, samples, drop = FALSE]
  mean_abund <- rowMeans(sub_df)
  tibble(Genus = names(mean_abund), Abundance = mean_abund, Chemotype = chemotype)
})

# 5. 合并并构建 Top6 + Others，同时保留所有属
all_abund <- bind_rows(core_taxa_abund) %>%
  group_by(Chemotype) %>%
  mutate(Rank = rank(-Abundance)) %>%
  mutate(Genus_label = if_else(Rank <= 6, Genus, "Others")) %>%
  group_by(Chemotype, Genus_label) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 找出所有组中出现的 Genus_label 项
genus_levels <- all_abund$Genus_label %>% unique() %>% sort()

# 6. 构建完整矩阵（缺失属填 0）用于统一图例
plot_df <- expand.grid(Chemotype = unique(all_abund$Chemotype),
                       Genus_label = genus_levels) %>%
  left_join(all_abund, by = c("Chemotype", "Genus_label")) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  group_by(Chemotype) %>%
  mutate(Percent = Abundance / sum(Abundance) * 100,
         ymax = cumsum(Percent),
         ymin = lag(ymax, default = 0),
         label_pos = (ymax + ymin) / 2)
# 7. 设置统一颜色
colors <- brewer.pal(n = max(8, length(genus_levels)), name = "Paired")
colors <- rep(colors, length.out = length(genus_levels))
names(colors) <- genus_levels

# 8. 绘图函数
make_donut <- function(data, chemotype_label) {
  ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Genus_label)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    scale_fill_manual(values = colors, breaks = genus_levels) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(chemotype_label) +
    geom_text(aes(x = 3.5, y = label_pos,
                  label = ifelse(Percent >= 5, paste0(round(Percent, 1), "%"), "")),
              color = "black", size = 3)
}
# 9. 分组绘图
p_list <- split(plot_df, plot_df$Chemotype) %>%
  imap(~ make_donut(.x, .y))
# 10. 合并图形，共用统一 legend
final_plot <- ggarrange(plotlist = p_list, ncol = 2, common.legend = TRUE, legend = "bottom")

##Fig2G
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

picrust_data <- read.table("Picrust2_EC_Fig2G.tsv",
                           sep = "\t", header = TRUE, check.names = FALSE, quote = "")
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(picrust_data) <- picrust_data$"function"
picrust_data2 <- picrust_data[, -c(1,2)]
count_data <- t(picrust_data2)
count_data <- round(count_data)
metadata <- metadata[match(rownames(count_data), metadata$Sample), ]

dds <- DESeqDataSetFromMatrix(countData = t(count_data),
                              colData = metadata,
                              design = ~ Chemotype)
gene_counts <- rowSums(counts(dds))
threshold <- quantile(gene_counts, 0.90)
dds <- dds[gene_counts > threshold, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("Chemotype", "HBA", "MSA"))

desc_df <- picrust_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("function_id") %>%
  dplyr::select(function_id, description)

res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("function_id") %>%
  left_join(desc_df, by = "function_id") %>%
  arrange(padj)

normalized_counts <- counts(dds, normalized = TRUE)
log2_normalized <- log2(normalized_counts + 1)

top20_ids <- head(res_df$function_id, 15)
top20_desc <- res_df$description[match(top20_ids,res_df$function_id)]
mat <- log2_normalized[top20_ids, ]
rownames(mat) <- top20_desc

annotation_col <- data.frame(
  Chemotype = metadata$Chemotype,
  chandi = metadata$chandi,
  row.names = metadata$Sample
)

ordered_samples <- annotation_col %>%
  rownames_to_column("SampleID") %>%
  arrange(factor(Chemotype, levels = c("HBA", "MSA"))) %>%
  pull(SampleID)

mat <- mat[, ordered_samples]
annotation_col <- annotation_col[ordered_samples, , drop = FALSE]

n_chemotype <- length(unique(annotation_col$Chemotype))
n_als <- length(unique(annotation_col$chandi))
ann_colors <- list(
  Chemotype = setNames(brewer.pal(n = max(3, n_chemotype), name = "Set3")[1:n_chemotype],
                       unique(annotation_col$Chemotype)),
  chandi = setNames(brewer.pal(n = max(3, 8), name = "Set3")[(n_chemotype + 1):(n_chemotype + n_als)],
                    unique(annotation_col$chandi))
)

pheatmap(
  mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0,
  color = colorRampPalette(c("#6e8fb2","#7da494","white","#de7833","#912c2c"))(100),
  main = "Top 15 Significant EC"
)


