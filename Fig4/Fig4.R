##Fig4
library(dplyr)
library(ggplot2)

# 读取数据
df <- read.table("root_Fig4.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 按输入顺序设置 x 轴顺序
df$group <- factor(df$group, levels = unique(df$group))

# 配对列表
pairs <- data.frame(
  group1 = c("TSB_CK1", "TSB_CK1","TSB_CK1", "PDB_CK", "TSA_CK", "TSA_CK", "TSA_CK", "TSA_CK", "TSA_CK",
             "TTC_CK", "TTC_CK", "YMA_CK", "YMA_CK"),
  group2 = c("S.chromofuscus", "S. cyaneochromogenes","S. mayteni", "Paenibacillus alvei",
             "TSA_Re359", "TSA_Re528", "TSA_Re1030", "TSA_Sa564", "TSA_Sa789",
             "TTC_Rs3041", "TTC_Rs7061", "YMA_By037", "YMA_By261")
)

# 计算均值和标准差（用于显著性标注高度）
df_stat <- df %>%
  group_by(group, color) %>%
  summarise(
    mean_val = mean(Root, na.rm = TRUE),
    sd_val = sd(Root, na.rm = TRUE),
    .groups = "drop"
  )

# 显著性检验
test_results <- pairs %>%
  rowwise() %>%
  mutate(
    p_value = t.test(
      Root ~ group,
      data = df %>% filter(group %in% c(group1, group2))
    )$p.value
  ) %>%
  mutate(
    signif_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 合并显著性信息到统计数据
df_stat <- df_stat %>%
  left_join(test_results %>% select(group2, signif_label),
            by = c("group" = "group2"))

# 绘图 —— 使用原始数据 df 画箱线 + 散点
p2 <- ggplot(df, aes(x = group, y = Root, fill = factor(color))) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +   # 箱线图
  geom_jitter(aes(color = factor(color)), width = 0.2, size = 2, alpha = 0.7) +  # 散点
  # 添加显著性标注
  geom_text(data = df_stat,
            aes(x = group, y = mean_val + sd_val + 0.05, label = signif_label),
            inherit.aes = FALSE, vjust = 0) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

library(dplyr)
library(ggplot2)

# 读取数据
df <- read.table("plant_Fig4.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(df)
# 按输入顺序设置 x 轴顺序
df$group <- factor(df$group, levels = unique(df$group))

# 配对列表
pairs <- data.frame(
  group1 = c("TSB_CK1", "TSB_CK1","TSB_CK1", "PDB_CK", "TSA_CK", "TSA_CK", "TSA_CK", "TSA_CK", "TSA_CK",
             "TTC_CK", "TTC_CK", "YMA_CK", "YMA_CK"),
  group2 = c("S.chromofuscus", "S. cyaneochromogenes","S. mayteni", "Paenibacillus alvei",
             "TSA_Re359", "TSA_Re528", "TSA_Re1030", "TSA_Sa564", "TSA_Sa789",
             "TTC_Rs3041", "TTC_Rs7061", "YMA_By037", "YMA_By261")
)

# 显著性检验
test_results <- pairs %>%
  rowwise() %>%
  mutate(
    p_value = t.test(
      plant ~ group,
      data = df %>% filter(group %in% c(group1, group2))
    )$p.value
  ) %>%
  mutate(
    signif_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 创建箱线图 + 散点图
p3 <- ggplot(df, aes(x = group, y = plant, fill = factor(color))) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +  # 箱线图
  geom_jitter(aes(color = factor(color)), width = 0.2, size = 2, alpha = 0.7) +  # 散点
  # 添加显著性标注
  geom_text(data = df_stat,
            aes(x = group, y = mean_val + sd_val + 0.05, label = signif_label),
            inherit.aes = FALSE, vjust = 0) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(
    panel.grid = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
library(dplyr)
library(tidyr)
df <- read.table("corebactest_Fig4.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df_wide <- df %>%
  select(sample, value, chemo, group, color) %>%
  pivot_wider(names_from = chemo, values_from = value)
# 计算比值
df_wide <- df_wide %>%
  mutate(ratio = (Atractylon + Atractylodin) / (`β_Eudesmol` + Hineson))
df_wide$group <- factor(df_wide$group, levels = unique(df_wide$group))
df_stat <- df_wide %>%
  group_by(group, color) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    .groups = "drop"
  )
library(ggplot2)
library(RColorBrewer)

# 先把 color 列转换为因子，保证 Paired 调色板映射
df_stat$color <- factor(df_stat$color, levels = unique(df_stat$color))
df_wide$color <- factor(df_wide$color, levels = levels(df_stat$color))

p1<-ggplot(df_stat, aes(x = group, y = mean_ratio, fill = color)) +
  # 每个样本的小点（透明）
  geom_point(data = df_wide, aes(x = group, y = ratio, fill = color),
             size = 2, shape = 21, alpha = 0.8, color = "black",
             inherit.aes = FALSE) +
  # 大点（均值）
  geom_point(shape = 21, size = 4, color = "black") +
  # 误差条
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), width = 0.2) +
  # 数值标签（在误差条上方，旋转60度）
  geom_text(aes(
    y = mean_ratio + sd_ratio + 5,# 在误差条顶端再往上加5个单位
    label = round(mean_ratio, 2)
  ),
  angle = 60, vjust = 0.5,hjust=0,size = 3, color = "black"
  ) +
  # 使用 Paired 调色板
  scale_fill_brewer(palette = "Paired") +
  labs(
    x = NULL,
    y = NULL,
    title = "Ratio (Atractylon+Atractylodin)/(β_Eudesmol+Hineson)"
  ) +
  coord_cartesian(ylim = c(-2, 200)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
library(patchwork)
(p3 / p2 / p1) + plot_layout(heights = c(1, 1, 2))