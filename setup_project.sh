#!/usr/bin/env bash
# ============================================================================
# setup_project.sh — 一键初始化 CCDS 标准目录骨架
# 适用于 Snakescale + TE_model 联合分析项目
# 用法: bash setup_project.sh
# ============================================================================
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo ">>> 正在初始化 CCDS 项目骨架: ${PROJECT_ROOT}"

# ── 1. data/ 四个标准子目录 ──────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/data/raw/fastq"
mkdir -p "${PROJECT_ROOT}/data/raw/genome"
mkdir -p "${PROJECT_ROOT}/data/interim"
mkdir -p "${PROJECT_ROOT}/data/processed"
mkdir -p "${PROJECT_ROOT}/data/external"

# ── 2. src/ 源代码目录 ──────────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/src/te_calc"
mkdir -p "${PROJECT_ROOT}/src/plots"
mkdir -p "${PROJECT_ROOT}/src/data"
mkdir -p "${PROJECT_ROOT}/src/features"
mkdir -p "${PROJECT_ROOT}/src/models"

# ── 3. workflow/ 流水线 ─────────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/workflow/snakescale"

# ── 4. reports/ 报告与图表 ──────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/reports/figures"

# ── 5. notebooks/ 探索性分析 ────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/notebooks"

# ── 6. models/ 序列化模型 ──────────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/models"

# ── 7. references/ 数据字典/手册 ────────────────────────────────────────────
mkdir -p "${PROJECT_ROOT}/references"

# ── 8. 创建占位 Dummy FASTQ / 基因组文件 ────────────────────────────────────
# 家蚕 (Bombyx mori) Ribo-seq / RNA-seq 占位
touch "${PROJECT_ROOT}/data/raw/fastq/bombyx_riboseq_R1.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/bombyx_riboseq_R2.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/bombyx_rnaseq_R1.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/bombyx_rnaseq_R2.fastq.gz"

# 蜘蛛 (Spider) Ribo-seq / RNA-seq 占位
touch "${PROJECT_ROOT}/data/raw/fastq/spider_riboseq_R1.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/spider_riboseq_R2.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/spider_rnaseq_R1.fastq.gz"
touch "${PROJECT_ROOT}/data/raw/fastq/spider_rnaseq_R2.fastq.gz"

# 参考基因组占位
touch "${PROJECT_ROOT}/data/raw/genome/bombyx_mori_genome.fa"
touch "${PROJECT_ROOT}/data/raw/genome/bombyx_mori_genome.fa.fai"
touch "${PROJECT_ROOT}/data/raw/genome/bombyx_mori_annotation.gtf"
touch "${PROJECT_ROOT}/data/raw/genome/spider_genome.fa"
touch "${PROJECT_ROOT}/data/raw/genome/spider_genome.fa.fai"
touch "${PROJECT_ROOT}/data/raw/genome/spider_annotation.gtf"

# ── 9. 在每个重要目录放一个 .gitkeep 保证 Git 追踪空目录 ────────────────────
for d in data/interim data/processed data/external \
         models notebooks references reports/figures; do
    touch "${PROJECT_ROOT}/${d}/.gitkeep"
done

# ── 10. 创建 Python __init__.py ─────────────────────────────────────────────
touch "${PROJECT_ROOT}/src/__init__.py"
touch "${PROJECT_ROOT}/src/te_calc/__init__.py"
touch "${PROJECT_ROOT}/src/plots/__init__.py"
touch "${PROJECT_ROOT}/src/data/__init__.py"
touch "${PROJECT_ROOT}/src/features/__init__.py"
touch "${PROJECT_ROOT}/src/models/__init__.py"

echo ">>> CCDS 项目骨架初始化完成！"
echo ""
echo "目录结构概览:"
# 尝试使用 tree; 若不可用则回退至 find
if command -v tree &>/dev/null; then
    tree -L 3 "${PROJECT_ROOT}"
else
    find "${PROJECT_ROOT}" -maxdepth 3 | sort | sed "s|${PROJECT_ROOT}|.|"
fi
