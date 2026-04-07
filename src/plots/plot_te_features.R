#!/usr/bin/env Rscript
# ============================================================================
# plot_te_features.R — 翻译效率 (TE) 可视化脚本
#
# 功能:
#   读取 data/processed/ 下的 TE 结果 CSV，使用 ggplot2 绘制：
#   1. 不同物种/细胞系间 TE 分布箱线图
#   2. 基因翻译效率对比散点图 (Ribo-seq vs RNA-seq)
#   3. TE 值排名柱状图 (Top N 高翻译效率基因)
#
# 输入: data/processed/*_te_results.csv
# 输出: reports/figures/*.pdf / *.png
#
# 用法:
#   Rscript src/plots/plot_te_features.R \
#       --input data/processed \
#       --output reports/figures
# ============================================================================

# ── 依赖加载 ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(argparse)
})

# ── 命令行参数解析 ───────────────────────────────────────────────────────────
parser <- ArgumentParser(description = "TE 特征可视化脚本")
parser$add_argument("--input",  type = "character",
                    default = "data/processed",
                    help = "TE 结果 CSV 所在目录")
parser$add_argument("--output", type = "character",
                    default = "reports/figures",
                    help = "图表输出目录")
parser$add_argument("--format", type = "character",
                    default = "pdf",
                    help = "输出图片格式: pdf / png")
args <- parser$parse_args()

# ── 创建输出目录 ─────────────────────────────────────────────────────────────
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

# ── 读取所有 TE 结果文件 ────────────────────────────────────────────────────
te_files <- list.files(args$input,
                       pattern = "_te_results\\.csv$",
                       full.names = TRUE)

if (length(te_files) == 0) {
  message("[plot_te_features] 未找到 TE 结果文件，使用占位数据进行演示...")
  # 生成占位数据用于测试
  set.seed(42)
  te_data <- data.frame(
    gene_id    = paste0("GENE_", sprintf("%04d", 1:100)),
    ribo_count = runif(100, 0, 500),
    rna_count  = runif(100, 10, 1000),
    te_value   = runif(100, 0.1, 5.0),
    species    = rep(c("bombyx", "spider"), each = 50)
  )
} else {
  te_data <- bind_rows(lapply(te_files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    # 从文件名提取物种标识
    species_tag <- gsub("_te_results\\.csv$", "", basename(f))
    df$species <- species_tag
    return(df)
  }))
}

# ============================================================================
# 图1: 不同物种/细胞系 TE 分布箱线图
# ============================================================================
# TODO: 根据实际细胞系/组织类型扩展分面变量
p1 <- ggplot(te_data, aes(x = species, y = te_value, fill = species)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title    = "翻译效率 (TE) 分布 — 物种间比较",
    subtitle = "TE = Ribo-seq / RNA-seq (加伪计数)",
    x        = "物种 / 细胞系",
    y        = "TE (log10 scale)",
    fill     = "Species"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(args$output, paste0("te_boxplot_species.", args$format)),
  plot     = p1,
  width    = 8, height = 6, dpi = 300
)
message("[plot_te_features] 已保存: te_boxplot_species.", args$format)

# ============================================================================
# 图2: Ribo-seq vs RNA-seq 散点图 (翻译效率对比)
# ============================================================================
# TODO: 添加回归线、标注高 TE / 低 TE 基因
p2 <- ggplot(te_data, aes(x = rna_count, y = ribo_count, color = species)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Ribo-seq vs RNA-seq Read Counts",
    subtitle = "对角线上方: TE > 1 (翻译活跃); 下方: TE < 1",
    x     = "RNA-seq Counts (log10)",
    y     = "Ribo-seq Counts (log10)",
    color = "Species"
  ) +
  theme_minimal(base_size = 14) +
  facet_wrap(~ species, scales = "free")

ggsave(
  filename = file.path(args$output, paste0("ribo_vs_rna_scatter.", args$format)),
  plot     = p2,
  width    = 12, height = 6, dpi = 300
)
message("[plot_te_features] 已保存: ribo_vs_rna_scatter.", args$format)

# ============================================================================
# 图3: Top N 高翻译效率基因柱状图
# ============================================================================
# TODO: 根据需求调整 Top N 数量、添加误差棒
top_n_genes <- 20

p3 <- te_data %>%
  group_by(species) %>%
  slice_max(order_by = te_value, n = top_n_genes) %>%
  ungroup() %>%
  mutate(gene_id = reorder(gene_id, te_value)) %>%
  ggplot(aes(x = gene_id, y = te_value, fill = species)) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    title = paste0("Top ", top_n_genes, " 高翻译效率基因"),
    x     = "Gene ID",
    y     = "TE Value",
    fill  = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(args$output, paste0("te_top_genes_bar.", args$format)),
  plot     = p3,
  width    = 10, height = 8, dpi = 300
)
message("[plot_te_features] 已保存: te_top_genes_bar.", args$format)

# ============================================================================
# TODO: 扩展图表
# - 热图: 多样本 TE 矩阵聚类 (pheatmap / ComplexHeatmap)
# - 火山图: 差异翻译效率基因 (需 DESeq2 / Xtail 输出)
# - 密度图: TE 分布核密度估计
# ============================================================================

message("[plot_te_features] 所有图表生成完毕！输出目录: ", args$output)
