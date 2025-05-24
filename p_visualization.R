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

# p测试数据集基因表达矩阵
p_expression <- read.csv('p_scgpt_expression.csv', header = TRUE,row.names = 1)
dim(p_expression)

# sc-GPT结果
scgpt_metadata <- read.csv('p_scgpt_metadata.csv')
head(scgpt_metadata)

# 用于singleR的训练和测试数据集
train_expression <- read.csv('demo_train_expression.csv', header = TRUE)
test_expression <- read.csv('demo_test_expression.csv', header = TRUE)
dim(train_expression)
dim(test_expression)
train_metadata <- read.csv('demo_train_metadata.csv')
test_metadata <- read.csv('demo_test_metadata.csv')


# TOSICA结果
tosica_expression<-read.csv("tosica_expression.csv",header=TRUE)
tosica_metadata <- read.csv('tosica_metadata.csv')
dim(tosica_expression) #4218  299

# 转换为矩阵
train_expression_matrix <- as.matrix(train_expression)
test_expression_matrix <- as.matrix(test_expression)
tosica_expression_matrix <- as.matrix(tosica_expression)

# 提取训练集的细胞类型标签
train_labels <- train_metadata$Celltype
length(unique(train_labels))

train_expression_matrix <- t(train_expression_matrix)
dim(train_expression_matrix)
dim(train_metadata)
test_expression_matrix <- t(test_expression_matrix)

# 使用 SingleR 进行注释
singleR_results <- SingleR(test = test_expression_matrix, ref = train_expression_matrix, labels = train_labels)

# 将注释结果添加到 test_metadata 中
test_metadata$SingleR_annotation <- singleR_results$labels

# 转换为矩阵
p_expression_matrix <- as.matrix(p_expression)
tosica_expression_matrix <- as.matrix(tosica_expression)

# 提取训练集的细胞类型标签
test_labels <- test_metadata$Celltype
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
  "#D4A6C8"   # 淡粉紫
)

# 如果 test_expression_matrix 是 genes × cells 的 matrix，先转置
expr_mat <- tosica_expression_matrix # 细胞 × 基因

# Step 1: Normalize total counts to 1e4 per cell
norm_mat <- sweep(expr_mat, 1, rowSums(expr_mat), FUN = "/") * 1e4

# Step 2: Log1p transformation
log_mat <- log1p(norm_mat)

# Step 3: Z-score scaling per gene
scaled_mat <- scale(log_mat, center = TRUE, scale = TRUE)
# 去掉只有一个唯一值的列（包括全 0）
clean_mat <- scaled_mat[, apply(scaled_mat, 2, function(x) length(unique(x)) > 1)]
pca <- prcomp(expr_mat, center = FALSE, scale. = FALSE)
pca_mat <- pca$x[, 1:40]  # 取前 40 个主成分
# Step 5: UMAP (same n_neighbors and min_dist as scanpy)
set.seed(36)
umap_coords <- umap(pca_mat, n_neighbors = 10, min_dist = 0.5, metric = "euclidean")


# SingleR
# 构建数据框用于可视化
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = test_metadata$Celltype,
  Predicted = test_metadata$SingleR_annotation
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Real)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "p Real Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

# Predicted 标签 UMAP
p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Predicted)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "SingleR Predicted Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )


p1 + p2 + plot_layout(ncol = 2)
ggsave("R_p.png", width = 8, height = 3, dpi = 300)


# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(test_metadata$Celltype, levels = label_level)
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
    title = "p SingleR Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("R_p_confusion_matrix.png", width = 6, height = 6, dpi = 300)


# scGPT

plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = scgpt_metadata$celltype,
  Predicted = scgpt_metadata$predictions
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Real)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "p Real Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

# Predicted 标签 UMAP
p2 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Predicted)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "scGPT Predicted Cell Types") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)  # <-- 使标题居中
  )

p1 + p2 + plot_layout(ncol = 2)
ggsave("scGPT_p.png", width = 8, height = 3, dpi = 300)

# scgpt confusion matrix
# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(scgpt_metadata$celltype, levels = label_level)
pred_labels <- factor(scgpt_metadata$predictions, levels = label_level)

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
    title = "p scGPT Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("scgpt_p_confusion_matrix.png", width = 6, height = 6, dpi = 300)



# TOSICA

plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Real = tosica_metadata$Celltype,
  Predicted = tosica_metadata$Prediction
)

# Real 标签 UMAP
p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Real)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_colors) +
  labs(title = "p Real Cell Types") +
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
ggsave("tosica_p.png", width = 8, height = 3, dpi = 300)

# tosica confusion matrix
# 显式指定 levels
label_level = unique(train_labels)
true_labels <- factor(tosica_metadata$Celltype, levels = label_level)
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
    title = "p TOSICA Confusion Matrix",
    x = "Predicted Label",
    y = "True Label",
    fill = "Proportion"
  )
# 保存图片（可选）
ggsave("tosica_p_confusion_matrix.png", width = 6, height = 6, dpi = 300)


