#调用R包
library(piecewiseSEM) 
library(nlme)
library(lme4)
library(semPlot)
library(psych)
library(BBmisc)
library(dplyr)
library(ggcorrplot)#加载包
library(ggpubr)
#读取数据
setwd("C:/Users/JN/Desktop")
mydata1 <- read.csv("rhizosphere.csv",header = T,row.names = 1)#读取数据
#mydata2 <- read.csv("rhizosphere.csv",header = T,row.names = 1)#读取数据

#mydata2<-normalize(mydata1[,3:12], method="standardize", margin=2)##标准化数据
#mydata3 <- mydata1[, (ncol(mydata1)-2):ncol(mydata1)]
#mydata<- cbind(mydata1[,1:2], mydata2, mydata3)

library(corrplot)#加载包
# 计算相关性
M <- cor(mydata1,method = "spearman")
res1 <- cor.mtest(mydata1, conf.level = .95)
corrplot(M,   type = 'upper', tl.pos = 'tp',
         method="circle",
         p.mat = res1$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05),
         pch.cex = 1.5,
         pch.col = "black",
         tl.cex=1)

PCA<-mydata1[,8:11]

head(PCA)

PCA1 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA1
PCA1<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
#head(PCA1)
summary(PCA1)

mydata1$env <- PCA1$x[,"PC1"]

##ID1
##根据模型建议
model1 <- psem(
  lm(env ~ Al , mydata1 ),
  lm(Bacteria_diversity ~ Al  + env+ Genotype , mydata1 ),
  lm(Fungi_diversity ~  Al  + env + Genotype , mydata1 ),
  lm(Bacteria_composition ~  Al  + env + Bacteria_diversity + Fungi_diversity + Fungi_composition +Genotype , mydata1 ),
  lm(Fungi_composition ~ Al  + env + Fungi_diversity + Genotype , mydata1 ),
  lm(Genotype ~  Al  + env , mydata1 ),
  lm(Chemotype ~  Al  + env + Bacteria_diversity + Fungi_diversity + Bacteria_composition + Fungi_composition + Genotype, mydata1 )
)
#查看模型结果
summary(model1, .progressBar = F)

#根际
mydata2 <- read.csv("rhizome.csv",header = T,row.names = 1)#读取数据

# 计算相关性
M <- cor(mydata2,method = "spearman")
res1 <- cor.mtest(mydata2, conf.level = .95)
corrplot(M,   type = 'upper', tl.pos = 'tp',
         method="circle",
         p.mat = res1$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05),
         pch.cex = 1.5,
         pch.col = "black",
         tl.cex=1)

PCA<-mydata2[,8:11]

head(PCA)

PCA1 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA1
PCA1<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA1)
summary(PCA1)

mydata2$env <- PCA1$x[,"PC1"]
##ID1
##根据模型建议
model2 <- psem(
  lm(env ~ Al , mydata2 ),
  lm(Bacteria_diversity_Shannon ~ Al  + env+ Genotype , mydata2 ),
  lm(Fungi_diversity_Shannon ~  Al  + env + Genotype , mydata2 ),
  lm(Bacteria_composition ~  Al  + env + Bacteria_diversity_Shannon + Fungi_diversity_Shannon+ Fungi_composition + Genotype , mydata2 ),
  lm(Fungi_composition ~ Al  + env + Fungi_diversity_Shannon + Genotype , mydata2 ),
  lm(Genotype ~  Al  + env , mydata2 ),
  lm(Chemotype ~  Al  + env + Bacteria_diversity_Shannon + Fungi_diversity_Shannon + Bacteria_composition + Fungi_composition + Genotype, mydata2 )
)
#查看模型结果
summary(model2, .progressBar = F)


##根茎土壤
mydata1 <- read.csv("rhizome_soil.csv",header = T,row.names = 1)#读取数据

#mydata2<-normalize(mydata1[,3:12], method="standardize", margin=2)##标准化数据
#mydata3 <- mydata1[, (ncol(mydata1)-2):ncol(mydata1)]
#mydata<- cbind(mydata1[,1:2], mydata2, mydata3)

library(corrplot)#加载包
# 计算相关性
M <- cor(mydata1,method = "spearman")
res1 <- cor.mtest(mydata1, conf.level = .95)
corrplot(M,   type = 'upper', tl.pos = 'tp',
         method="circle",
         p.mat = res1$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05),
         pch.cex = 1.5,
         pch.col = "black",
         tl.cex=1)

PCA<-mydata1[,1:5]

head(PCA)

PCA1 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA1
PCA1<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA1)
summary(PCA1)

PCA<-mydata1[,7:8]

head(PCA)

PCA2 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA2
PCA2<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA2)
summary(PCA2)


PCA<-mydata1[,9:10]

head(PCA)

PCA3 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA3
PCA3<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA3)
summary(PCA3)


PCA<-mydata1[,12:15]

head(PCA)

PCA4 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA4
PCA4<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA4)
summary(PCA4)

mydata1$soil <- PCA1$x[,"PC1"]
mydata1$Bacteria <- PCA2$x[,"PC1"]
mydata1$Fungi <- PCA3$x[,"PC1"]
mydata1$Climate <- PCA4$x[,"PC1"]

##根据模型建议
model1 <- psem(
  #lm(Climate ~ Al , mydata1 ),
  lm(soil ~ Climate , mydata1 ),
  lm(Bacteria ~ Climate+ soil , mydata1 ),
  lm(Fungi ~  Climate  + soil , mydata1 ),
  #lm(Bacteria_composition ~  Al  + Climate + Bacteria_diversity_shannon + Fungi_diversity_shannon+ Fungi_composition + soil , mydata1 ),
  #lm(Fungi_composition ~ Al  + Climate + soil + Fungi_diversity_shannon  , mydata1 ),
  lm(Chemotype ~  Climate + Bacteria + Fungi+ soil, mydata1 )
)

#查看模型结果
summary(model1, .progressBar = F)


###根际土壤
mydata2 <- read.csv("rhizosphere_soil.csv",header = T,row.names = 1)#读取数据

#mydata2<-normalize(mydata1[,3:12], method="standardize", margin=2)##标准化数据
#mydata3 <- mydata1[, (ncol(mydata1)-2):ncol(mydata1)]
#mydata<- cbind(mydata1[,1:2], mydata2, mydata3)

library(corrplot)#加载包
# 计算相关性
M <- cor(mydata2,method = "spearman")
res1 <- cor.mtest(mydata2, conf.level = .95)
corrplot(M,   type = 'upper', tl.pos = 'tp',
         method="circle",
         p.mat = res1$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05),
         pch.cex = 1.5,
         pch.col = "black",
         tl.cex=1)

PCA<-mydata2[,1:5]

head(PCA)

PCA1 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA1
PCA1<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA1)
summary(PCA1)

PCA<-mydata2[,7:8]

head(PCA)

PCA2 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA2
PCA2<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA2)
summary(PCA2)


PCA<-mydata2[,9:10]

head(PCA)

PCA3 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA3
PCA3<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA3)
summary(PCA3)


PCA<-mydata2[,12:15]

head(PCA)

PCA4 <- principal(PCA, nfactors = 1, rotate ="varimax")
PCA4
PCA4<-prcomp(PCA, scale= TRUE)
#write.table (PCA1$x, file ="pc1.csv",sep =",", quote =FALSE)
head(PCA4)
summary(PCA4)

mydata2$soil <- PCA1$x[,"PC1"]
mydata2$Bacteria <- PCA2$x[,"PC1"]
mydata2$Fungi <- PCA3$x[,"PC1"]
mydata2$Climate <- PCA4$x[,"PC1"]

##ID1
##根据模型建议
model2 <- psem(
  #lm(Climate ~ Al , mydata1 ),
  lm(soil ~ Climate , mydata2 ),
  lm(Bacteria ~ Climate+ soil , mydata2 ),
  lm(Fungi ~  Climate  + soil , mydata2 ),
  #lm(Bacteria_composition ~  Al  + Climate + Bacteria_diversity_shannon + Fungi_diversity_shannon+ Fungi_composition + soil , mydata1 ),
  #lm(Fungi_composition ~ Al  + Climate + soil + Fungi_diversity_shannon  , mydata1 ),
  lm(Chemotype ~  Climate + Bacteria + Fungi+ soil, mydata2 )
)

#查看模型结果
summary(model2, .progressBar = F)


#提取计算数据
kk1 <- summary(model2, .progressBar = F)
kk2 <- kk1$coefficients[,-(3:6)]
kk3 <- as.data.frame(kk2)
kk4 <- kk3[!grepl("~~", kk3$Predictor), ]
kk4

# 结果整理为便于计算的数据框
result1 <- kk4 %>%
  dplyr::relocate(
    from   = Predictor,
    to     = Response,
    weight = Std.Estimate,
    p      = P.Value
  )
#查看整理后结果
result1  
#开始计算
#===========================直接效应============================
# 示例数据
data <-result1 

# 定义计算直接效应的函数
calculateDirectEffects <- function(data, factors) {
  # 初始化直接效应向量
  direct_effects <- numeric(nrow(data))
  
  # 遍历数据集中的每一行
  for (i in 1:nrow(data)) {
    # 检查路径是否在因子列表中，如果是则直接效应即为weight值，否则为0
    if (paste(data$from[i], data$to[i], sep = "_") %in% factors) {
      direct_effects[i] <- data$weight[i]
    } else {
      direct_effects[i] <- 0
    }
  }
  
  # 过滤 Direct_Effect 等于 0 的行
  non_zero_rows <- direct_effects != 0
  direct_effects <- direct_effects[non_zero_rows]
  from_to <- paste(data$from[non_zero_rows], "→", data$to[non_zero_rows])
  
  # 创建包含直接效应的数据框
  direct_effects_df <- data.frame(
    from_to = from_to,
    Direct_Effect = direct_effects
  )
  
  return(direct_effects_df)
}


# 定义因子列表
# 自动化生成直接效应的因子列表的函数
generateDirectEffectFactors <- function(data) {
  # 从数据集中获取唯一的 from 和 to 组合，符合直接效应的条件
  unique_direct_combinations <- unique(paste(data$from, data$to, sep = "_"))
  
  # 提取起点和终点
  from_factors <- unique(data$from)
  to_factors <- unique(data$to)
  
  # 拆分因子并转换为向量
  direct_factors <- paste(rep(from_factors, each = length(to_factors)), rep(to_factors, length(from_factors)), sep = "_")
  
  return(direct_factors)
}

# 调用函数生成直接效应的因子列表
direct_factors <- generateDirectEffectFactors(data)

# 调用函数计算直接效应
direct_effects_result <- calculateDirectEffects(data, direct_factors)

print(direct_effects_result)
#==========================================间接效应==========================================
# 示例数据
data <- result1  

# 从数据中生成路径列表
unique_from <- unique(data$from)
unique_to <- unique(data$to)
paths <- expand.grid(from = unique_from, to = unique_to)
paths <- split(paths, seq(nrow(paths)))
paths <- lapply(paths, function(x) c(as.character(x$from), as.character(x$to), "DM"))

# 定义计算间接效应的函数
calculateIndirectEffects <- function(data, from_factor, through_factor, to_factor) {
  through_weight <- data$weight[data$from == from_factor & data$to == through_factor]
  to_weight <- data$weight[data$from == through_factor & data$to == to_factor]
  
  if (length(through_weight) == 0 | length(to_weight) == 0) {
    message(paste("Skipping invalid path:", from_factor, "→", through_factor, "→", to_factor))
    return(NULL)
  }
  
  indirect_effect <- through_weight * to_weight
  
  return(indirect_effect)
}

# 创建空的数据框来存储间接效应结果
indirect_effects_df <- data.frame(from_to = character(), Indirect_Effect = numeric())

# 循环计算间接效应并存储结果
for (path in paths) {
  from_factor <- path[1]
  through_factor <- path[2]
  to_factor <- path[3]
  
  indirect_effect <- calculateIndirectEffects(data, from_factor, through_factor, to_factor)
  if (!is.null(indirect_effect)) {
    from_to <- paste(from_factor, through_factor, to_factor, sep = " → ")
    indirect_effects_df <- rbind(indirect_effects_df, data.frame(from_to = from_to, Indirect_Effect = indirect_effect))
  }
}

# 打印结果
print(indirect_effects_df)
#然后将间接效应的和求出来
# 创建示例数据框
indirect_effects_df <- indirect_effects_df

# 定义自动化函数
calculate_total_indirect_effect <- function(data) {
  # 提取开头和结尾相同的值，并对其"Indirect_Effect"列求和
  total_indirect_effect <- data %>%
    mutate(
      start_pattern = sub(" → .*", "", from_to),
      end_pattern = sub(".* → ", "", from_to)
    ) %>%
    group_by(start_pattern, end_pattern) %>%
    summarise(total_indirect_effect = sum(Indirect_Effect), .groups = "drop") %>%
    ungroup() %>%
    arrange(start_pattern, end_pattern) %>%
    mutate(from_to = paste0(start_pattern, " → ", end_pattern)) %>%
    select(from_to, total_indirect_effect)
  
  return(total_indirect_effect)
}

# 调用自动化函数并打印结果
result_total_indirect_effect <- calculate_total_indirect_effect(indirect_effects_df)
print(result_total_indirect_effect)

#合并三个结果
# 直接效应数据框
direct_effects_result <- direct_effects_result

# 间接效应数据框
indirect_effects_df <- indirect_effects_df
# 总间接效应数据框
total_indirect_effect_df <- as.data.frame(result_total_indirect_effect)

# 合并直接效应和间接效应数据框
total_effects_df <- bind_rows(direct_effects_result, indirect_effects_df, total_indirect_effect_df)

# 打印结果
print(total_effects_df)

# 创建数据框
data <- total_effects_df

# 使用dplyr包中的group_by和summarise函数，将相同因子的数据转移到相同行并求和
data_processed <- data %>%
  group_by(from_to) %>%
  summarise(Direct_Effect = sum(Direct_Effect, na.rm = TRUE),
            Indirect_Effect = sum(Indirect_Effect, na.rm = TRUE),
            total_indirect_effect = sum(total_indirect_effect, na.rm = TRUE))

# 打印处理后的数据框
print(data_processed)

# 使用dplyr包中的mutate函数，将Direct_Effect和total_indirect_effect列求和
data_processed2 <- data_processed %>%
  mutate(Total_Effect = Direct_Effect + total_indirect_effect)

# 打印处理后的数据框
print(data_processed2)

# 过滤不需要的路径，我们的目标因子是TRAD，所以我们不需要那些中间路径的值
df_filtered1 <- data_processed2 %>%
  filter(grepl("Chemotype$", from_to))

# 输出结果
df_filtered1
write.table (df_filtered1, file ="total.csv",sep =",", quote =FALSE)
