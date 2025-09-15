########Fig3#########
#Fig3A/C,将读取数据从raw_Fig3A.txt换成raw_Fig3C.txt即可
rawgen <- read.table("raw_Fig3A.txt", header = TRUE, sep = "\t", check.names = FALSE)

library(tidyverse)
library(ggpubr)
library(patchwork)
group_colors <- c(
  a = "#febb8b",
  b = "#ffe19d",
  c = "#8c96c7",
  d = "#b4cde5"
)
# 将数据变长
rawgen_long <- rawgen %>%
  pivot_longer(cols = c(a, b, c, d), names_to = "group", values_to = "value") %>%
  mutate(group = factor(group, levels = c("a", "b", "c", "d")))
# 创建每个 type 的图
make_type_plot <- function(type_val) {
  df_sub <- rawgen_long %>% filter(type == type_val)
  
  # 计算每个组的平均值（用于柱状图高度）
  mean_df <- df_sub %>%
    group_by(group) %>%
    summarise(mean_val = mean(value, na.rm = TRUE))
  
  ggplot(df_sub, aes(x = group, y = value)) +
    geom_col(
      data = mean_df,
      aes(y = mean_val, color = group),
      fill = NA,
      width = 0.5,
      linewidth = 1.2
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar", 
                 width = 0.3, size = 1.2, color = "black") +
    # 均值点
    stat_summary(fun = mean, geom = "point", 
                 size = 4, shape = 21, color = "black", fill = "white") +
    # 原始散点
    geom_jitter(width = 0.15, size = 2, aes(color = group)) +
    scale_color_manual(values = group_colors) +
    scale_fill_identity() +
    stat_compare_means(
      comparisons = list(c("a", "b"), c("c", "d")),
      method = "t.test",
      label = "p.signif"
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title.x = element_blank()
    )
}

# 创建每个 type 的图
type_list <- unique(rawgen$type)
type_plots <- map(type_list, make_type_plot)

# 合并图，统一 legend
final_plot <- wrap_plots(type_plots, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")
#Fig3BD
library(tidyverse)
library(patchwork)

#读取数据
df <- read_tsv("Feng_Fig3BDE.txt")

#转为长格式
df_long <- df %>%
  pivot_longer(cols = c(a, b, c, d), names_to = "SPECI", values_to = "Value")

#设置因子顺序
compound_levels <- c("Hinesol", "β-Eudesmol", "Atractylon", "Atractylodin")
spe_order <- c("a", "b", "c", "d")
df_long <- df_long %>%
  mutate(
    type = factor(type, levels = compound_levels),
    SPECI = factor(SPECI, levels = spe_order)
  )

#计算均值、标准差
df_stat <- df_long %>%
  group_by(SPECI, type) %>%
  summarise(
    mean = mean(Value),
    sd = sd(Value),
    .groups = "drop"
  ) %>%
  mutate(center = as.numeric(type))

#构建正态峰函数
generate_peak <- function(mean, sd, center, compound, SPECI) {
  x_seq <- seq(-0.5, 0.5, length.out = 100)
  x_full <- x_seq + center
  y <- dnorm(x_seq, mean = 0, sd = 0.2)
  y_scaled <- y / max(y) * mean
  tibble(x = x_full, y = y_scaled, Compound = compound,
         center = center, mean = mean, sd = sd, SPECI = SPECI)
}

#生成所有峰形数据
df_peaks <- df_stat %>%
  pmap_dfr(~generate_peak(..3, ..4, ..5, ..2, ..1))
df_peaks$SPECI <- factor(df_peaks$SPECI, levels = spe_order)
df_stat$SPECI <- factor(df_stat$SPECI, levels = spe_order)

#自定义颜色
compound_colors <- c(
  "Hinesol" = "#cd76ab",
  "β-Eudesmol" = "#ffcc00",
  "Atractylon" = "#52B6e6",
  "Atractylodin" = "#3eb0a2"
)

#分样本绘图
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
                         labels = levels(df_stat$type)) +
      scale_fill_manual(values = compound_colors) +
      labs(title = paste("Sample", .y), x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )
  })

#拼图
wrap_plots(plot_list[spe_order], nrow = 1)
#Fig3E
library(tidyverse)
library(patchwork)

#读取数据
df <- read.table("Feng_Fig3BDE.txt", header = TRUE, sep = "\t", check.names = FALSE)

#计算比值
calculate_ratio <- function(colname) {
  df_ratio <- df %>%
    select(sample, type, group, !!sym(colname)) %>%
    pivot_wider(names_from = type, values_from = !!sym(colname)) %>%
    mutate(ratio = (Atractylodin + Atractylon) / (`Hinesol` + `β-Eudesmol`)) %>%
    select(sample, group, ratio) %>%
    mutate(Type = ifelse(group == "n", "bn", "bs"))
  return(df_ratio)
}

#绘图函数
make_ratio_plot <- function(ratio_df, label) {
  df_stat <- ratio_df %>%
    group_by(Type) %>%
    summarise(mean = mean(ratio), sd = sd(ratio), .groups = "drop")
  
  earth_colors <- c("bn" = "#8B5A2B", "bs" = "#CDAA7D")
  
  ggplot(ratio_df, aes(x = ratio, y = Type)) +
    geom_point(aes(fill = Type), shape = 21, size = 4, alpha = 0.8, stroke = 0) +
    geom_errorbarh(data = df_stat,
                   aes(xmin = mean - sd, xmax = mean + sd, y = Type),
                   height = 0.1, color = "black", linewidth = 0.6,
                   inherit.aes = FALSE) +
    geom_point(data = df_stat,
               aes(x = mean, y = Type, fill = Type),
               shape = 21, color = "black", size = 5, stroke = 1.2,
               inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    scale_fill_manual(values = earth_colors) +
    labs(
      x = paste0("Ratio (", label, ")"),
      y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "bold"),
      legend.position = "none"
    )
}

#生成图
plot_a <- make_ratio_plot(calculate_ratio("a"), "a")
plot_b <- make_ratio_plot(calculate_ratio("b"), "b")
plot_c <- make_ratio_plot(calculate_ratio("c"), "c")
plot_d <- make_ratio_plot(calculate_ratio("d"), "d")

#拼图
final_plot <- (plot_a | plot_b) / (plot_c | plot_d)

