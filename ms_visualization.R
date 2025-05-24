library(SingleR)
library(celldex)
library(uwot)
library(ggplot2)
library(patchwork)
library(caret)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(scales)
library(biomaRt)

# p测试数据集基因表达矩阵
ms_expression <- read.csv('ms_scgpt_expression.csv', header = TRUE,row.names = 1)
dim(ms_expression)

# sc-GPT结果
ms_metadata <- read.csv('ms_scgpt_metadata.csv')
head(scgpt_metadata)

# 用于singleR的训练和测试数据集
train_expression <- read.csv('ms_train_expression.csv', header = TRUE,row.names = 1)
test_expression <- read.csv('ms_test_expression.csv', header = TRUE,row.names = 1)
dim(train_expression)
dim(test_expression)
train_metadata <- read.csv('ms_train_metadata.csv')
test_metadata <- read.csv('ms_test_metadata.csv')


# TOSICA结果
tosica_expression<-read.csv("ms_tosica_expression.csv",header=TRUE)
tosica_metadata <- read.csv('ms_tosica_metadata.csv')
dim(tosica_expression) #13468 300

# 转换为矩阵
train_expression_matrix <- as.matrix(train_expression)
test_expression_matrix <- as.matrix(ms_expression)
dim(train_expression_matrix)
tosica_expression_matrix <- as.matrix(tosica_expression)
rownames(test_expression_matrix)[1:5]
rownames(train_expression_matrix)[1:5]

train_expression_matrix <- t(train_expression_matrix)
test_expression_matrix <- t(test_expression_matrix)

# 提取训练集的细胞类型标签
train_labels <- train_metadata$celltype
length(unique(train_labels))
# # 连接 Ensembl 数据库
# mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")
# ensembl_ids <- rownames(test_expression_matrix)
# conversion <- getBM(
#   attributes = c("ensembl_gene_id", "hgnc_symbol"),
#   filters = "ensembl_gene_id",
#   values = ensembl_ids,
#   mart = mart
# )
# # 去除空 symbol，并去重
# conversion <- conversion[conversion$hgnc_symbol != "", ]
# conversion <- conversion[!duplicated(conversion$ensembl_gene_id), ]
# # 替换 test 的行名为 gene symbol
# rownames(test_expression_matrix) <- conversion$hgnc_symbol[
#   match(rownames(test_expression_matrix), conversion$ensembl_gene_id)
# ]
# # 删除未匹配的行（NA symbol）
# test_expression_matrix <- test_expression_matrix[!is.na(rownames(test_expression_matrix)), ]

# 使用 SingleR 进行注释
singleR_results <- SingleR(test = test_expression_matrix, ref = train_expression_matrix, labels = train_labels)

# 将注释结果添加到 test_metadata 中
test_metadata$SingleR_annotation <- singleR_results$labels

# 转换为矩阵
ms_expression_matrix <- as.matrix(ms_expression)
tosica_expression_matrix <- as.matrix(tosica_expression)

# 提取训练集的细胞类型标签
test_labels <- test_metadata$celltype
length(unique(test_labels))

# 提取所有可能的标签（真实 + 预测）
cell_types <- unique(train_labels)
# 自动生成颜色向量
cell_colors <- c(
  "#4E79A7",  # 柔和蓝
  "#F28E2B",  # 橙金色
  "#E15759",  # 珊瑚红
  "#76B7B2",  # 薄荷绿
  "#59A14F",  # 草绿色
  "#EDC948",  # 柠檬黄
  "#B07AA1",  # 淡紫
  "#FF9DA7",  # 粉玫瑰
  "#9C755F",  # 摩卡棕
  "#BAB0AC",  # 淡灰褐
  "#8CD17D",  # 青草绿
  "#FAA43A",  # 亮橘
  "#D4A6C8",  # 淡粉紫
  "#6B4C9A",  # 深丁香紫
  "#F7C6C7",  # 樱花粉
  "#A1D99B",  # 柔绿
  "#FDD0A2",  # 杏桃色
  "#C6DBEF"   # 雾霭蓝
)


# PCA
expr_mat <- ms_expression_matrix # 细胞 × 基因
pca <- prcomp(expr_mat, center = FALSE, scale. = FALSE)
pca_mat <- pca$x[, 1:40]  # 取前 40 个主成分
# Step 5: UMAP (same n_neighbors and min_dist as scanpy)
set.seed(35)
umap_coords <- umap(pca_mat, n_neighbors = 10, min_dist = 0.5, metric = "euclidean")


# 构建数据框用于可视化
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = test_metadata$celltype,
  Predicted = test_metadata$SingleR_annotation
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Real)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "ms Real Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

# Predicted 标签 UMAP
p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Predicted)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "singleR Predicted Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )


p1 + p2 + plot_layout(ncol = 2)
ggsave("R_ms_2.png", width = 12, height = 5, dpi = 300)

# SingleR confusion matrix
# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(test_metadata$celltype, levels = label_level)
pred_labels <- factor(test_metadata$SingleR_annotation, levels = label_level)

# 保证 levels 一致
common_levels <- union(levels(true_labels), levels(pred_labels))
real <- factor(true_labels, levels = common_levels)
predictions <- factor(pred_labels, levels = common_levels)

# 计算混淆矩阵
cm <- confusionMatrix(predictions, real)
cm_table <- as.data.frame(cm$table)

# 基础 R 实现按 Reference（真实标签）归一化
cm_table$Proportion <- unlist(
  by(cm_table, cm_table$Reference, function(sub) {
    sub$Freq / sum(sub$Freq)
  })
)
# 已有 confusionMatrix 对象 cm
accuracy <- cm$overall['Accuracy']
print(accuracy)
# 取出各类 recall（也称 sensitivity）
recall <- cm$byClass[ , "Sensitivity"]
# 取出各类 precision（也称 Positive Predictive Value）
precision <- cm$byClass[ , "Pos Pred Value"]
# F1-score = 2 * (precision * recall) / (precision + recall)
f1 <- 2 * (precision * recall) / (precision + recall)
# 打印结果
data.frame(Label = rownames(cm$byClass), Recall = recall, Precision = precision, F1 = f1)
macro_recall <- mean(recall, na.rm = TRUE)
macro_precision <- mean(precision, na.rm = TRUE)
macro_f1 <- mean(f1, na.rm = TRUE)
cat("Macro Recall:", macro_recall, "\n")
cat("Macro F1-score:", macro_f1, "\n")


# 绘制归一化混淆矩阵热图
ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Proportion)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )+
  labs(
    title = "ms SingleR Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("R_ms_confusion_matrix.png", width = 8, height = 8, dpi = 300)


# scGPT
label_level = unique(train_labels)
true_labels <- factor(ms_metadata$celltype, levels = label_level)
pred_labels <- factor(ms_metadata$predictions, levels = label_level)
all_labels <- union(levels(factor(true_labels)), levels(factor(pred_labels)))
true_labels <- factor(true_labels, levels = all_labels)
pred_labels <- factor(pred_labels, levels = all_labels)

cell_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#8CD17D", "#FAA43A", "#D4A6C8", "#6B4C9A", "#F7C6C7",
  "#A1D99B", "#FDD0A2", "#C6DBEF"
)
names(cell_colors) <- all_labels

plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = ms_metadata$celltype,
  Predicted = ms_metadata$predictions
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = true_labels)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "Real ms Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

# Predicted 标签 UMAP
p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = pred_labels)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "scGPT Predicted Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

p1 + p2 + plot_layout(ncol = 2)
ggsave("scGPT_ms_2.png", width = 12, height = 5, dpi = 300)


# scgpt confusion matrix
# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(ms_metadata$celltype, levels = label_level)
pred_labels <- factor(ms_metadata$predictions, levels = label_level)
length(unique(true_labels))
# 保证 levels 一致
common_levels <- union(levels(true_labels), levels(pred_labels))
real <- factor(true_labels, levels = common_levels)
predictions <- factor(pred_labels, levels = common_levels)

# 计算混淆矩阵
cm <- confusionMatrix(predictions, real)
cm_table <- as.data.frame(cm$table)

# 基础 R 实现按 Reference（真实标签）归一化
cm_table$Proportion <- unlist(
  by(cm_table, cm_table$Reference, function(sub) {
    sub$Freq / sum(sub$Freq)
  })
)
# 已有 confusionMatrix 对象 cm
accuracy <- cm$overall['Accuracy']
print(accuracy)
# 取出各类 recall（也称 sensitivity）
recall <- cm$byClass[ , "Sensitivity"]
# 取出各类 precision（也称 Positive Predictive Value）
precision <- cm$byClass[ , "Pos Pred Value"]
# F1-score = 2 * (precision * recall) / (precision + recall)
f1 <- 2 * (precision * recall) / (precision + recall)
# 打印结果
data.frame(Label = rownames(cm$byClass), Recall = recall, Precision = precision, F1 = f1)
macro_recall <- mean(recall, na.rm = TRUE)
macro_precision <- mean(precision, na.rm = TRUE)
macro_f1 <- mean(f1, na.rm = TRUE)
cat("Macro Recall:", macro_recall, "\n")
cat("Macro F1-score:", macro_f1, "\n")


# 绘制归一化混淆矩阵热图
ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Proportion)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )+
  labs(
    title = "ms scGPT Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("scgpt_ms_confusion_matrix.png", width = 8, height = 8, dpi = 300)




# TOSICA

plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = tosica_metadata$celltype,
  Predicted = tosica_metadata$Prediction
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Real)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "ms Real Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

# Predicted 标签 UMAP
p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Predicted)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "TOSICA Predicted Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

p1 + p2 + plot_layout(ncol = 2)
ggsave("tosica_ms_2.png", width = 12, height = 5, dpi = 300)

# tosica confusion matrix
# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(tosica_metadata$celltype, levels = label_level)
pred_labels <- factor(tosica_metadata$Prediction, levels = label_level)

# 保证 levels 一致
common_levels <- union(levels(true_labels), levels(pred_labels))
real <- factor(true_labels, levels = common_levels)
predictions <- factor(pred_labels, levels = common_levels)

# 计算混淆矩阵
cm <- confusionMatrix(predictions, real)
cm_table <- as.data.frame(cm$table)

# 基础 R 实现按 Reference（真实标签）归一化
cm_table$Proportion <- unlist(
  by(cm_table, cm_table$Reference, function(sub) {
    sub$Freq / sum(sub$Freq)
  })
)
# 已有 confusionMatrix 对象 cm
accuracy <- cm$overall['Accuracy']
print(accuracy)
# 取出各类 recall（也称 sensitivity）
recall <- cm$byClass[ , "Sensitivity"]
# 取出各类 precision（也称 Positive Predictive Value）
precision <- cm$byClass[ , "Pos Pred Value"]
# F1-score = 2 * (precision * recall) / (precision + recall)
f1 <- 2 * (precision * recall) / (precision + recall)
# 打印结果
data.frame(Label = rownames(cm$byClass), Recall = recall, Precision = precision, F1 = f1)
macro_recall <- mean(recall, na.rm = TRUE)
macro_precision <- mean(precision, na.rm = TRUE)
macro_f1 <- mean(f1, na.rm = TRUE)
cat("Macro Recall:", macro_recall, "\n")
cat("Macro F1-score:", macro_f1, "\n")

# 绘制归一化混淆矩阵热图
ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Proportion)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )+
  labs(
    title = "ms TOSICA Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("tosica_ms_confusion_matrix.png", width = 8, height = 8, dpi = 300)


