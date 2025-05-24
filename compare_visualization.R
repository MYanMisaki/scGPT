# 加载必要包
library(ggplot2)
library(dplyr)
library(ggthemes)
library(RColorBrewer)

p_metadata <- read.csv('p_scgpt_metadata.csv')
ms_metadata <- read.csv('ms_scgpt_metadata.csv')

# 统计每种 cell type 数量
celltype_counts <- p_metadata %>%
  count(celltype) %>%
  arrange(desc(n))

# 绘图
ggplot(celltype_counts, aes(x = reorder(celltype, -n), y = n, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal(base_size = 13) +
  labs(title = "Distribution of hPancreas Cell Types",
       x = "Cell Type",
       y = "Number of Cells") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave("p_celltype.png", width =7, height = 5, dpi = 300)

# 加载必要包
library(ggplot2)
library(dplyr)

# 自定义颜色向量
custom_colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                   "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
                   "#8CD17D", "#FAA43A", "#D4A6C8", "#6B4C9A", "#F7C6C7",
                   "#A1D99B", "#FDD0A2", "#C6DBEF")

# 读取数据
p_metadata <- read.csv('p_scgpt_metadata.csv')

# 统计 celltype
celltype_counts <- ms_metadata %>%
  count(celltype) %>%
  arrange(desc(n))

# 绘制纵向柱状图
ggplot(celltype_counts, aes(x = reorder(celltype, -n), y = n, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "MS Cell Type Distribution",
       x = "Cell Type",
       y = "Number of Cells") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave("ms_celltype.png", width =12, height = 6, dpi = 300)


library(ggplot2)
library(tidyr)

# 构造原始数据
df <- data.frame(
  dataset = rep(c("hPancreas", "MS"), each = 3),
  metric = rep(c("Accuracy", "macro recall", "macro F1"), times = 2),
  scGPT = c(0.9644, 0.9103, 0.8531, 0.8305, 0.6532, 0.7963),
  TOSICA = c(0.9556, 0.9136, 0.8465, 0.6411, 0.5724, 0.5422),
  SingleR = c(0.9877, 0.9833, 0.9415, 0.6884, 0.6153, 0.6211)
)

library(ggplot2)
library(tidyr)

# 转换为长格式
df_long <- pivot_longer(df, cols = c("scGPT", "TOSICA", "SingleR"),
                        names_to = "Method", values_to = "Score")

# 折线图
ggplot(df_long, aes(x = metric, y = Score, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ dataset) +
  theme_minimal() +
  labs(
    title = "Performance Comparison of Methods",
    x = "Metric",
    y = "Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
ggsave("metrix.png", width =8, height = 4, dpi = 300)

