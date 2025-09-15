#############################Fig1##################
##Fig1A was generated using Adobe Illustrator.
##Fig1B NJ tree visualization was performed using the Interactive Tree Of Life (iTOL) tool.
##Fig1B
library(tidyverse)
library(patchwork)
# 1. 读取数据
df <- read_csv("feng_Fig1B.csv")
# 2. 转换为长格式
df_long <- df %>%
  pivot_longer(cols = c(Hinesol, `β-Eudesmol`, Atractylon, Atractylodin),
               names_to = "Compound", values_to = "Value")
# 3. 设置因子顺序
spe_order <- c("HBYX_W", "SXXY_W", "HNSX_W", "JSTS_W", "AHSZ_W", "HBLT_W")
df_long <- df_long %>%
  mutate(
    SPECI = factor(SPECI, levels = spe_order),
    Compound = factor(Compound, levels = c("Hinesol", "β-Eudesmol", "Atractylon", "Atractylodin"))
  )
# 4. 计算均值和标准差
df_stat <- df_long %>%
  group_by(SPECI, Compound) %>%
  summarise(
    mean = mean(Value),
    sd = sd(Value),
    .groups = "drop"
  ) %>%
  mutate(center = as.numeric(Compound))
# 5. 定义正态峰函数
generate_peak <- function(mean, sd, center, compound, SPECI) {
  x_seq <- seq(-0.5, 0.5, length.out = 100)
  x_full <- x_seq + center
  y <- dnorm(x_seq, mean = 0, sd = 0.2)
  y_scaled <- y / max(y) * mean
  tibble(x = x_full, y = y_scaled, Compound = compound, center = center, mean = mean, sd = sd, SPECI = SPECI)
}
# 6. 生成所有峰数据
df_peaks <- df_stat %>%
  pmap_dfr(~generate_peak(..3, ..4, ..5, ..2, ..1))
df_peaks$SPECI <- factor(df_peaks$SPECI, levels = spe_order)
df_stat$SPECI <- factor(df_stat$SPECI, levels = spe_order)
# 7. 设置颜色
compound_colors <- c(
  "Hinesol" = "#cd76ab",
  "β-Eudesmol" = "#ffcc00",
  "Atractylon" = "#52B6e6",
  "Atractylodin" = "#3eb0a2"
)
# 8. 分样本绘制
plot_list <- df_peaks %>%
  split(.$SPECI) %>%
  imap(~{
    df_stat_sub <- df_stat %>% filter(SPECI == .y)
    ggplot(.x, aes(x = x, y = y, fill = Compound)) +
      geom_area(alpha = 0.7, color = "black") +
      geom_errorbar(data = df_stat_sub,
                    aes(x = center, ymin = mean - sd, ymax = mean + sd),
                    width = 0.1, color = "black", inherit.aes = FALSE) +
      geom_point(data = df_stat_sub,
                 aes(x = center, y = mean),
                 shape = 21, fill = "white", color = "black", size = 2,
                 inherit.aes = FALSE) +
      scale_x_continuous(breaks = 1:4,
                         labels = levels(df_stat$Compound)) +
      scale_fill_manual(values = compound_colors) +
      labs(title = NULL, x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )
  })
# 9. 拼接子图
wrap_plots(plot_list[spe_order], ncol = 1)

##Fig1C
library(ggplot2)
library(ggprism)
library(patchwork)
library(ggrain)

# 读取数据
df <- read.csv("huanjingFig1C.csv")

# 绘制单个雨滴图函数
plot_rain <- function(df, yvar, ylab) {
  ggplot(df, aes(1, .data[[yvar]], fill = Genotype, color = Genotype)) +
    geom_rain(alpha = 1) +
    scale_fill_manual(values = c("MA" = "#D9EDE4", "SA" = "#246333")) +
    scale_color_manual(values = c("MA" = "#D9EDE4", "SA" = "#246333")) +
    guides(fill = 'none', color = 'none') +
    labs(x = ylab, y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.border = element_blank()
    ) +
    coord_flip()
}

# 三个性状的雨滴图
p1 <- plot_rain(df, "Altitude", "Altitude")
p2 <- plot_rain(df, "Precipitation", "Precipitation")
p3 <- plot_rain(df, "Temperature", "Temperature")

# 拼接为一列三行
p1 / p2 / p3

##Fig1D
library(tidyverse)

# 读取数据
df <- read_csv("AA_hbFig1D.csv")

# 计算中位数和四分位距
stat_df <- df %>%
  group_by(group, Genotype) %>%
  summarise(
    median = median(AA_HB),
    Q1 = quantile(AA_HB, 0.25),
    Q3 = quantile(AA_HB, 0.75),
    .groups = "drop"
  )

# 自定义颜色
color_map <- c("MA" = "#BCE8D3", "SA" = "#246333")

# 绘图
ggplot() +
  # 样本点
  geom_point(
    data = df,
    aes(y = fct_reorder(group, AA_HB), x = AA_HB, color = Genotype),
    size = 4, alpha = 0.7
  ) +
  # 四分位距
  geom_errorbarh(
    data = stat_df,
    aes(y = fct_reorder(group, median), xmin = Q1, xmax = Q3, group = Genotype, color = Genotype),
    height = 0.2, linewidth = 0.7
  ) +
  # 中位数点
  geom_point(
    data = stat_df,
    aes(y = fct_reorder(group, median), x = median, fill = Genotype),
    shape = 21, color = "black", size = 6, stroke = 0.3
  ) +
  # 配色与样式
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw() +
  labs(x = "AA_HB", y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )


