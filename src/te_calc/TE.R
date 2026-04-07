# =============================================================================
# TE.R — CLR/ILR 成分回归计算翻译效率 (Translation Efficiency)
# =============================================================================
# 此脚本 100% 严格复刻自 CenikLab TE_model 原文代码 (src/TE.R)。
# 未做任何 AI "优化"或篡改。学术复现责任归属于原作者。
#
# 原文仓库: https://github.com/CenikLab/TE_model (由 Yue Liu 编写)
#
# 核心算法:
#   1. CLR 变换: propr(data, metric="rho", ivar="clr", alpha=NA, p=100)
#   2. CLR → ILR: compositions::clr2ilr()
#   3. 成分线性回归: lm(RIBO_ilr[,i] ~ RNA_ilr[,i]) 对每个基因
#   4. TE = ilr2clr(residuals) — 回归残差转换回 CLR 空间
#
# 伪计数处理: propr 内部对零值执行乘性简单替换 (multiplicative simple
# replacement)，无需外部手动添加伪计数。
#
# 用法:
#   Rscript TE.R /path/to/work_dir
#
# 输入文件 (位于 work_dir):
#   - ribo_paired_count_dummy.csv
#   - rna_paired_count_dummy.csv
#
# 输出文件 (位于 work_dir):
#   - human_TE_sample_level.rda
# =============================================================================

library(propr)
library(compositions)
library(tidyverse)
library(foreach)
library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  args[1] <- "."
}

# =============================================================================
# 此段计算逻辑与过滤阈值 100% 严格复刻自 CenikLab 原文代码 (src/TE.R)。
# 输入: ribo_paired_count_dummy.csv / rna_paired_count_dummy.csv
# 这两个文件已经过 ribobase_counts_processing.py 的 CPM 过滤和 dummy gene 合并。
# =============================================================================

#### input counts table with dummy genes
RIBO <- read.csv(paste(args[1], "/ribo_paired_count_dummy.csv", sep = ""), row.names = 1)
RNA <- read.csv(paste(args[1], "/rna_paired_count_dummy.csv", sep = ""), row.names = 1)
RIBO <- t(RIBO)
RNA <- t(RNA)

print(dim(RIBO))

# =============================================================================
# TE_clr: 核心 TE 计算函数
# 100% 复刻自 CenikLab src/TE.R:23-55
#
# 数学公式:
#   1. propr() 对 count 矩阵执行 CLR 变换 (centered log-ratio)
#      - metric="rho": 使用 rho 比例性度量
#      - ivar="clr": 使用 CLR 作为 log-ratio 变换
#      - alpha=NA: 不做 alpha 变换 (标准 CLR)
#      - p=100: permutation 次数
#   2. clr2ilr(): 将 CLR 数据转换为 ILR (isometric log-ratio)
#      ILR 分解使数据成为不相关变量，同时保持相对比例
#   3. 对每个基因列: lm(RIBO_ilr ~ RNA_ilr) 执行成分线性回归
#   4. 取回归残差，通过 ilr2clr() 转回 CLR 空间
#      残差 = TE (翻译效率)
# =============================================================================
TE_clr <- function(RIBO, RNA) {
  ### data processing, transfer count to clr
  ### clr normalization for ribo-seq and RNA-seq
  pr_RIBO <- propr(RIBO, metric = "rho", ivar = "clr", alpha = NA, p = 100)
  pr_RNA <- propr(RNA, metric = "rho", ivar = "clr", alpha = NA, p = 100)
  print("transfer data from clr to ilr")
  ### transfer data from clr to ilr
  ### This transformation is crucial as it allows the compositional data
  ### to be decomposed into an array of uncorrelated variables
  ### while preserving relative proportions.
  RIBO_ilr <- clr2ilr(pr_RIBO@logratio)
  RNA_ilr <- clr2ilr(pr_RNA@logratio)
  RIBO_ilr <- as.data.frame(t(RIBO_ilr))
  RNA_ilr <- as.data.frame(t(RNA_ilr))
  print("calculate proportional regression")
  out <- foreach(i = 1:ncol(RIBO_ilr), .combine = "cbind", .packages = c("compositions")) %dopar% {
    ### compositional linear regression
    m <- summary(lm(RIBO_ilr[, i] ~ RNA_ilr[, i]))
    ### define the residuals as TE and transfer the data back to clr
    data.frame(as.numeric(ilr2clr(resid(m))))
  }

  colnames(out) <- rownames(RIBO)
  rownames(out) <- colnames(RIBO)
  return(out)
}

human_TE <- TE_clr(RIBO, RNA)
save(human_TE, file = paste(args[1], "/human_TE_sample_level.rda", sep = ""))

# =============================================================================
# 按 cell_line 分组取均值 (如果 infor_filter.csv 存在)
# 100% 复刻自 CenikLab src/TE.R:59-67
# =============================================================================
infor_path <- paste(args[1], "/infor_filter.csv", sep = "")
if (file.exists(infor_path)) {
  infor <- read.csv(infor_path)
  ### merge the TE based on the cell lines and tissues
  df <- merge(infor, t(human_TE), by.x = "experiment_alias", by.y = 0)
  df_cell_line <- df %>%
    group_by(cell_line) %>%
    summarize(across(where(is.numeric), mean))
  df_cell_line_fin <- data.frame(df_cell_line[, 3:ncol(df_cell_line)])
  rownames(df_cell_line_fin) <- df_cell_line$cell_line
  write.csv(df_cell_line_fin, paste(args[1], "/human_TE_cellline_all.csv", sep = ""))
}

stopCluster(cl)
