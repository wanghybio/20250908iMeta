library(randomForest)
library(A3)
library(rfPermute)
library(BBmisc)
library(ggplot2)
library(extrafont)
library(viridis)
library(RColorBrewer)
#BiocManager::install("tibble")
# Create the forest.
setwd("C:/Users/JN/Desktop")
tolerance<-read.csv(file="rhizome_soil.csv",row.names = 1)
#tolerance<-normalize(tolerance, method="range", margin=2)
set.seed(1234)
#tolerance=adjusted_data
output.forest <-randomForest(tolerance$Chemotype ~ .,
                             ntree = 1000,nrep=1000,data=tolerance,importance=T,na.action=na.roughfix)
print(output.forest) # View the forest results.
round(importance(output.forest), 2)
#print(importance(output.forest, type = 2))# Importance of each predictor.
#R^2 and p-value  library(A3)

a3(tolerance$Chemotype ~ .,
   data=tolerance, randomForest, p.acc = 0.01)

ozone.rfP <- rfPermute(tolerance$Chemotype ~ .,
                       data=tolerance, ntree = 1000, na.action = na.omit, nrep = 1000,num.cores=1)
layout(matrix(1:6, nrow = 2))
#plotNull(ozone.rfP) 
layout(matrix(1))
#RFimportance=importance(output.forest, scale = T)
RFimportance=importance(ozone.rfP, scale = T)
plot(importance(output.forest, scale =T))#贡献度和显著性散点图
write.csv(RFimportance, "rhizome_soil随机森林Chemotype.csv")######random forest
##excel中处理下数据名称删除%符号再导入，手动调整下分组_soil
H<-read.csv(file="rhizome_soil随机森林Chemotype.csv")
#提取预测变量的显著
H$stars <- ifelse(H$IncMSE.pval < 0.001, "***",
                  ifelse(H$IncMSE.pval < 0.01, "**",
                         ifelse(H$IncMSE.pval < 0.05, "*", "")))

#font_import()
#loadfonts(device = "win")
# 创建柱状图
rainbow_colors <- colorRampPalette(rainbow(7))(22)
colors <- colorRampPalette(brewer.pal(27, "Spectral"))(length(H$IncMSE))
colors1 <- rev(c("#FDFFC0", "#FDCF91","#E05069", "#340E61", "#0B0411"))
p1<-ggplot(H, aes(x=IncMSE, y=reorder(ID, IncMSE), fill=IncMSE)) +#factor(ID,levels=rev(ID))
  geom_bar(stat="identity",width=0.7) +
  #scale_color_manual(values = rainbow_colors) +
  #scale_fill_viridis(discrete = TRUE) +
  scale_fill_gradientn(colors = "#FDCF91") +
  geom_text(aes(label=stars), position=position_dodge(width=0.7),hjust=-1, vjust=0.7, size=6) +
  theme_bw() +
  labs( x="Increase in MSE (%)", size=14) +
  theme(legend.position="none", 
        panel.grid = element_blank(),
        axis.text.y=element_text(color="black", size=14), 
        axis.title.x=element_text(size=14,vjust=-.5),
        axis.text.x=element_text(hjust=1, color="black", size=14), 
        axis.title.y=element_blank()
        )+ # 移除y轴标题
  xlim(0, 24) + # 设置x轴的范围为0到40
  annotate("text", x=15, y=3, label="R² = 99.6%, p < 0.001", size=6) 
  #annotate("text", x=20, y=3, label="", size=6) # 在图的右下角添加额外的信息

p1
ggsave(p1,file="rhizome_soil随机森林Chemotype.pdf",width =6,height = 5,unit="in",dpi=300) #储存


##
tolerance<-read.csv(file="rhizosphere_soil.csv",row.names = 1)
#tolerance<-normalize(tolerance, method="range", margin=2)
set.seed(1234)
#tolerance=adjusted_data
output.forest <-randomForest(tolerance$Chemotype ~ .,
                             ntree = 1000,nrep=1000,data=tolerance,importance=T,na.action=na.roughfix)
print(output.forest) # View the forest results.
round(importance(output.forest), 2)
#print(importance(output.forest, type = 2))# Importance of each predictor.
#R^2 and p-value  library(A3)

a3(tolerance$Chemotype ~ .,
   data=tolerance, randomForest, p.acc = 0.01)

ozone.rfP <- rfPermute(tolerance$Chemotype ~ .,
                       data=tolerance, ntree = 1000, na.action = na.omit, nrep = 1000,num.cores=1)
layout(matrix(1:6, nrow = 2))
#plotNull(ozone.rfP) 
layout(matrix(1))
#RFimportance=importance(output.forest, scale = T)
RFimportance=importance(ozone.rfP, scale = T)
plot(importance(output.forest, scale =T))#贡献度和显著性散点图
write.csv(RFimportance, "rhizosphere_soil随机森林Chemotype.csv")######random forest
##excel中处理下数据名称删除%符号再导入，手动调整下分组
H<-read.csv(file="rhizosphere_soil随机森林Chemotype.csv")
#提取预测变量的显著
H$stars <- ifelse(H$IncMSE.pval < 0.001, "***",
                  ifelse(H$IncMSE.pval < 0.01, "**",
                         ifelse(H$IncMSE.pval < 0.05, "*", "")))

#font_import()
#loadfonts(device = "win")
# 创建柱状图
rainbow_colors <- colorRampPalette(rainbow(7))(22)
colors <- colorRampPalette(brewer.pal(27, "Spectral"))(length(H$IncMSE))
colors1 <- rev(c("#FDFFC0", "#FDCF91","#E05069", "#340E61", "#0B0411"))
p2<-ggplot(H, aes(x=IncMSE, y=reorder(ID, IncMSE), fill=IncMSE)) +#factor(ID,levels=rev(ID))
  geom_bar(stat="identity",width=0.7) +
  #scale_color_manual(values = rainbow_colors) +
  #scale_fill_viridis(discrete = TRUE) +
  scale_fill_gradientn(colors = "#FDCF91") +
  geom_text(aes(label=stars), position=position_dodge(width=0.7),hjust=-1, vjust=0.7, size=6) +
  theme_bw() +
  labs( x="Increase in MSE (%)", size=14) +
  theme(legend.position="none", 
        panel.grid = element_blank(),
        axis.text.y=element_text(color="black", size=14), 
        axis.title.x=element_text(size=14,vjust=-.5),
        axis.text.x=element_text(hjust=1, color="black", size=14), 
        axis.title.y=element_blank()
  )+ # 移除y轴标题
  xlim(0, 25) + # 设置x轴的范围为0到40
  annotate("text", x=15, y=3, label="R² = 99.6%, p < 0.001", size=6) 
#annotate("text", x=20, y=3, label="", size=6) # 在图的右下角添加额外的信息

p2
ggsave(p2,file="rhizosphere_soil随机森林Chemotype.pdf",width =6,height = 5,unit="in",dpi=300) #储存


####
tolerance<-read.csv(file="rhizome.csv",row.names = 1)
#tolerance<-normalize(tolerance, method="range", margin=2)
set.seed(1234)
#tolerance=adjusted_data
output.forest <-randomForest(tolerance$Chemotype ~ .,
                             ntree = 1000,nrep=1000,data=tolerance,importance=T,na.action=na.roughfix)
print(output.forest) # View the forest results.
round(importance(output.forest), 2)
#print(importance(output.forest, type = 2))# Importance of each predictor.
#R^2 and p-value  library(A3)

a3(tolerance$Chemotype ~ .,
   data=tolerance, randomForest, p.acc = 0.01)

ozone.rfP <- rfPermute(tolerance$Chemotype ~ .,
                       data=tolerance, ntree = 1000, na.action = na.omit, nrep = 1000,num.cores=1)
layout(matrix(1:6, nrow = 2))
#plotNull(ozone.rfP) 
layout(matrix(1))
#RFimportance=importance(output.forest, scale = T)
RFimportance=importance(ozone.rfP, scale = T)
plot(importance(output.forest, scale =T))#贡献度和显著性散点图
write.csv(RFimportance, "rhizome随机森林Chemotype.csv")######random forest
##excel中处理下数据名称删除%符号再导入，手动调整下分组
H<-read.csv(file="rhizome随机森林Chemotype.csv")
#提取预测变量的显著
H$stars <- ifelse(H$IncMSE.pval < 0.001, "***",
                  ifelse(H$IncMSE.pval < 0.01, "**",
                         ifelse(H$IncMSE.pval < 0.05, "*", "")))

#font_import()
#loadfonts(device = "win")
# 创建柱状图
rainbow_colors <- colorRampPalette(rainbow(7))(22)
colors <- colorRampPalette(brewer.pal(27, "Spectral"))(length(H$IncMSE))
colors1 <- rev(c("#FDFFC0", "#FDCF91","#E05069", "#340E61", "#0B0411"))
p3<-ggplot(H, aes(x=IncMSE, y=reorder(ID, IncMSE), fill=IncMSE)) +#factor(ID,levels=rev(ID))
  geom_bar(stat="identity",width=0.7) +
  #scale_color_manual(values = rainbow_colors) +
  #scale_fill_viridis(discrete = TRUE) +
  scale_fill_gradientn(colors = "#FDCF91") +
  geom_text(aes(label=stars), position=position_dodge(width=0.7),hjust=-1, vjust=0.7, size=6) +
  theme_bw() +
  labs( x="Increase in MSE (%)", size=14) +
  theme(legend.position="none", 
        panel.grid = element_blank(),
        axis.text.y=element_text(color="black", size=14), 
        axis.title.x=element_text(size=14,vjust=-.5),
        axis.text.x=element_text(hjust=1, color="black", size=14), 
        axis.title.y=element_blank()
  )+ # 移除y轴标题
  xlim(-4, 28) + # 设置x轴的范围为0到40
  annotate("text", x=15, y=3, label="R² = 84.01%, p < 0.001", size=6) 
#annotate("text", x=20, y=3, label="", size=6) # 在图的右下角添加额外的信息

p3
ggsave(p3,file="rhizome随机森林Chemotype.pdf",width =6,height = 5,unit="in",dpi=300) #储存


####
tolerance<-read.csv(file="rhizosphere.csv",row.names = 1)
#tolerance<-normalize(tolerance, method="range", margin=2)
set.seed(1234)
#tolerance=adjusted_data
output.forest <-randomForest(tolerance$Chemotype ~ .,
                             ntree = 1000,nrep=1000,data=tolerance,importance=T,na.action=na.roughfix)
print(output.forest) # View the forest results.
round(importance(output.forest), 2)
#print(importance(output.forest, type = 2))# Importance of each predictor.
#R^2 and p-value  library(A3)

a3(tolerance$Chemotype ~ .,
   data=tolerance, randomForest, p.acc = 0.01)

ozone.rfP <- rfPermute(tolerance$Chemotype ~ .,
                       data=tolerance, ntree = 1000, na.action = na.omit, nrep = 1000,num.cores=1)
layout(matrix(1:6, nrow = 2))
#plotNull(ozone.rfP) 
layout(matrix(1))
#RFimportance=importance(output.forest, scale = T)
RFimportance=importance(ozone.rfP, scale = T)
plot(importance(output.forest, scale =T))#贡献度和显著性散点图
write.csv(RFimportance, "rhizosphere随机森林Chemotype.csv")######random forest
##excel中处理下数据名称删除%符号再导入，手动调整下分组
H<-read.csv(file="rhizosphere随机森林Chemotype.csv")
#提取预测变量的显著
H$stars <- ifelse(H$IncMSE.pval < 0.001, "***",
                  ifelse(H$IncMSE.pval < 0.01, "**",
                         ifelse(H$IncMSE.pval < 0.05, "*", "")))

#font_import()
#loadfonts(device = "win")
# 创建柱状图
rainbow_colors <- colorRampPalette(rainbow(7))(22)
colors <- colorRampPalette(brewer.pal(27, "Spectral"))(length(H$IncMSE))
colors1 <- rev(c("#FDFFC0", "#FDCF91","#E05069", "#340E61", "#0B0411"))
p4<-ggplot(H, aes(x=IncMSE, y=reorder(ID, IncMSE), fill=IncMSE)) +#factor(ID,levels=rev(ID))
  geom_bar(stat="identity",width=0.7) +
  #scale_color_manual(values = rainbow_colors) +
  #scale_fill_viridis(discrete = TRUE) +
  scale_fill_gradientn(colors = "#FDCF91") +
  geom_text(aes(label=stars), position=position_dodge(width=0.7),hjust=-1, vjust=0.7, size=6) +
  theme_bw() +
  labs( x="Increase in MSE (%)", size=14) +
  theme(legend.position="none", 
        panel.grid = element_blank(),
        axis.text.y=element_text(color="black", size=14), 
        axis.title.x=element_text(size=14,vjust=-.5),
        axis.text.x=element_text(hjust=1, color="black", size=14), 
        axis.title.y=element_blank()
  )+ # 移除y轴标题
  xlim(-4, 28) + # 设置x轴的范围为0到40
  annotate("text", x=15, y=3, label="R² = 84.01%, p < 0.001", size=6) 
#annotate("text", x=20, y=3, label="", size=6) # 在图的右下角添加额外的信息

p4
ggsave(p,file="rhizosphere随机森林Chemotype.pdf",width =6,height = 5,unit="in",dpi=300) #储存
