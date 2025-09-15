##Fig5
#Fig5A
df <- read.table("huifayou_Fig5A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df[, -1] <- lapply(df[, -1], as.numeric)
# 按组计算平均值
library(dplyr)
group_mean <- df %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
# 如果需要矩阵形式
group_mean_mat <- as.data.frame(group_mean)
rownames(group_mean_mat) <- group_mean_mat$Group
group_mean_mat <- group_mean_mat[, -1]
# 输出结果
huifay<-group_mean_mat

df <- read.table("turang_Fig5A.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(ggplot2)
# 去掉含 NA 的行
df <- na.omit(merge)
# 指定 x 和 y 列
x_vars <- c("Hinesol", "β.Eudesmol", "Atractylon", "Atractylodin", "aa_hb")
y_vars <- setdiff(colnames(df), c("Row.names", x_vars))
# 初始化结果表
cor_res <- expand.grid(Xvar = x_vars, Yvar = y_vars, stringsAsFactors = FALSE)
cor_res$rho <- NA
cor_res$pvalue <- NA
# 循环做 Spearman 相关
for (i in seq_len(nrow(cor_res))) {
  x <- df[[cor_res$Xvar[i]]]
  y <- df[[cor_res$Yvar[i]]]
  if (var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) next
  test <- suppressWarnings(cor.test(x, y, method = "spearman"))
  cor_res$rho[i] <- test$estimate
  cor_res$pvalue[i] <- test$p.value
}

# FDR 校正
cor_res$padj <- p.adjust(cor_res$pvalue, method = "fdr")

# 标注显著性，只保留 padj < 0.05 的
cor_res$sig <- ifelse(cor_res$pvalue < 0.001, "***",
                      ifelse(cor_res$pvalue < 0.01, "**",
                             ifelse(cor_res$pvalue < 0.05, "*", "")))

# 作图
ggplot(cor_res, aes(x = Xvar, y = Yvar, fill = rho, size = abs(rho))) +
  geom_point(shape = 21, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Spearman\nrho", size = "|rho|", x = "", y = "")

#Fig5B
library(vegan)
library(ggplot2)
library(patchwork)
library(ggrepel)

# 修改后的CCA绘图函数（仅显示数值型变量）
create_cca_plot <- function(abundance, metadata, title = "") {
  # 转置数据并匹配样本
  abundance_t <- t(abundance)
  common_samples <- intersect(rownames(abundance_t), rownames(metadata))
  abundance_t <- abundance_t[common_samples, ]
  metadata <- metadata[common_samples, ]
  
  # 排除非环境变量列（保留所有数值型变量）
  env_vars <- setdiff(colnames(metadata), c("Sample", "Genotype", "Chemotype")) 
  
  # 执行CCA分析
  formula <- as.formula(paste("abundance_t ~", paste(env_vars, collapse = "+")))
  ord <- cca(formula, data = metadata)
  
  # 提取样本点坐标
  sites <- scores(ord, display = "sites", choices = 1:2)
  sites_df <- as.data.frame(sites)
  colnames(sites_df) <- c("CCA1", "CCA2")
  sites_df$Chemotype <- metadata$Chemotype
  sites_df$Sample <- rownames(sites_df)
  
  # 提取环境变量箭头（仅数值型变量）
  env_arrows <- scores(ord, display = "bp", choices = 1:2)
  env_df <- if(!is.null(env_arrows)) {
    data.frame(
      CCA1 = env_arrows[,1] * 5,  # 调整箭头缩放比例
      CCA2 = env_arrows[,2] *5,
      Var = rownames(env_arrows)
    )
  } else {
    data.frame(CCA1 = numeric(0), CCA2 = numeric(0), Var = character(0))
  }
  
  # 计算解释的方差
  eig <- ord$CCA$eig
  variance1 <- ifelse(length(eig)>0, round(eig[1]/sum(eig)*100, 1), 0)
  variance2 <- ifelse(length(eig)>1, round(eig[2]/sum(eig)*100, 1), 0)
  
  # 创建基础绘图（使用theme_bw，隐藏网格线）
  p <- ggplot() +
    geom_vline(xintercept = 0,linetype = "dashed", color = "gray60", alpha = 0.6) +
    geom_vline(yintercept = 0,linetype = "dashed", color = "gray60", alpha = 0.6) +
    geom_point(data = sites_df, 
               aes(x = CCA1, y = CCA2, color = Chemotype), 
               size = 2, alpha = 0.8) +
    scale_color_manual(values = c("HBA" = "#1B9E77", "MSA" = "#D95F02")) +
    labs(title = paste0(title, " (CCA)"),
         x = paste0("CCA1 (", variance1, "%)"),
         y = paste0("CCA2 (", variance2, "%)")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold")
    )
  
  # 添加环境变量箭头（数值型）
  if(nrow(env_df) > 0){
    p <- p + 
      geom_segment(
        data = env_df,
        aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
        arrow = arrow(length = unit(0.15, "cm")),
        color = "black", alpha = 0.8, linewidth = 0.6
      ) +
      geom_text_repel(
        data = env_df,
        aes(x = CCA1, y = CCA2, label = Var),
        color = "black", size = 2,
        min.segment.length = 0.3,
        box.padding = 0.4
      )
  }
  
  # 添加总解释率
  p <- p + 
    annotate("text", x = Inf, y = Inf,
             label = paste("Total variance explained:", 
                           round(sum(ord$CCA$eig)/ord$tot.chi*100, 1), "%"),
             hjust = 1.05, vjust = 1.2, size = 3, color = "gray30")
  
  return(p)
}

# 数据读取
bac_gennei <- read.table("bac_gennei.tsv", header = TRUE, row.names = 1, sep = "\t")
bac_genjie <- read.table("bac_genjie.tsv", header = TRUE, row.names = 1, sep = "\t")
fun_gennei <- read.table("fun_gennei.tsv", header = TRUE, row.names = 1, sep = "\t")
fun_genjie <- read.table("fun_genjie.tsv", header = TRUE, row.names = 1, sep = "\t")

gennei_meta_bac <- read.table("gennei_meta_bac.tsv", header = TRUE, row.names = 1, sep = "\t")
gennei_meta_fun <- read.table("gennei_meta_fun.tsv", header = TRUE, row.names = 1, sep = "\t")
genjie_meta_bac <- read.table("genjie_meta_bac.tsv", header = TRUE, row.names = 1, sep = "\t")
genjie_meta_fun <- read.table("genjie_meta_fun.tsv", header = TRUE, row.names = 1, sep = "\t")

# 创建CCA图（仅显示数值型变量）
p1 <- create_cca_plot(bac_gennei, gennei_meta_bac, "Bacteria (Gennei)")
p2 <- create_cca_plot(bac_genjie, genjie_meta_bac, "Bacteria (Genjie)")
p3 <- create_cca_plot(fun_gennei, gennei_meta_fun, "Fungi (Gennei)")
p4 <- create_cca_plot(fun_genjie, genjie_meta_fun, "Fungi (Genjie)")

# 组合图形
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# 保存图形（高分辨率）
ggsave("cca_numeric_vars.png", combined_plot, 
       width = 12, height = 10, dpi = 300)