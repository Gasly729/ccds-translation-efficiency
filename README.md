# Multi-Species Translation Efficiency Analysis Pipeline
# 多物种翻译效率分析流水线

[![Python](https://img.shields.io/badge/python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/Gasly729/ccds-translation-efficiency)](https://github.com/Gasly729/ccds-translation-efficiency)

---

## 一句话简介

基于 **Ribo-seq + RNA-seq**，使用 CenikLab CLR/ILR 成分回归算法，定量分析 human / mouse / arabidopsis / yeast 四物种翻译效率（TE），探究跨物种翻译调控保守模式。

> Quantifies **translation efficiency (TE)** across human, mouse, arabidopsis and yeast using CenikLab's CLR/ILR compositional regression on Ribo-seq + RNA-seq data, enabling cross-species analysis of conserved translational regulation.

---

## 核心结果预览 / Key Results

> 所有结果基于 2026-04-09 修复后的流水线重算（各物种独立计算 + outer-join 合并）。
> All results recomputed with the bug-fixed pipeline (2026-04-09): per-species independent TE + outer-join merge.

| Species | Samples | Genes | Mean TE Variance |
|---------|---------|-------|-----------------|
| Human | 18 | 10,611 | 0.465 |
| Mouse | 11 | 10,864 | 0.652 |
| Arabidopsis | 4 | 15,947 | 0.943 |
| Yeast | 3 | 5,436 | 0.538 |
| **Total (outer join)** | **36** | **42,852** | **0.582** |

---

## 项目结构 / Project Structure

```
project/
│
├── Makefile                          ← 顶层统一调度器
├── requirements.txt                  ← Python 依赖
├── README.md
│
├── data/
│   ├── raw/                          ← 原始数据（只读）
│   ├── interim/                      ← 中间产物
│   ├── processed/                    ← 最终输出
│   │   ├── te_human.csv              ←   Human TE 结果
│   │   ├── te_mouse.csv              ←   Mouse TE 结果
│   │   ├── te_arabidopsis.csv        ←   Arabidopsis TE 结果
│   │   ├── te_yeast.csv              ←   Yeast TE 结果
│   │   ├── te_results_final.csv      ←   跨物种 outer join 合并
│   │   └── sample_annotation.csv     ←   样本元数据（已追踪）
│   └── external/                     ← 外部辅助数据
│
├── src/
│   ├── te_calc/
│   │   ├── te_calculator.py          ← ★ TE 计算主入口（Python）
│   │   └── TE.R                      ← ★ CLR/ILR 核心算法（R，复刻 CenikLab）
│   ├── plots/
│   │   └── generate_pub_style_figures.py  ← ★ 出版级图形生成
│   └── data/
│       └── download_sra.py           ← SRA 数据下载
│
├── workflow/
│   └── snakescale/                   ← 上游 Ribo-seq/RNA-seq 流水线
│       ├── Snakefile                 ←   规则定义（QC→Trim→Align→Count）
│       ├── config.yaml               ←   全局参数与 SLURM 配置
│       └── scripts/
│           ├── generate_yaml.py      ←   CSV 驱动的元数据→YAML 转换
│           └── references.yaml       ←   21+ 物种 Bowtie2 索引路由表
│
└── reports/
    └── figures/pub_style/            ← 出版级图形输出
        ├── Fig1_PCA_*.pdf
        ├── Fig2_TE_median_*.pdf
        ├── Fig3_violin_*.pdf
        └── Fig_UMAP_TE_v3.pdf        ← 跨物种 UMAP（新增）
```

---

## 快速开始 / Quick Start

### 环境依赖 / Dependencies

- **Python**：conda base 环境，含 `ribopy`、`pandas`、`numpy`、`scikit-learn`、`umap-learn`
- **R 环境**：`snakemake-ribo` conda 环境，含 `propr`、`compositions`、`foreach`、`doParallel`
- **Rscript 路径**：`/home/xrx/miniconda3/envs/snakemake-ribo/bin/Rscript`

#### 验证 R 依赖 / Verify R dependencies

```bash
/home/xrx/miniconda3/envs/snakemake-ribo/bin/Rscript -e \
  "library(propr); library(compositions); library(foreach); library(doParallel); cat('OK\n')"
```

### 运行全流程 TE 计算 / Run full TE pipeline

```bash
cd /home/xrx/my_project/project
RSCRIPT_BIN=/home/xrx/miniconda3/envs/snakemake-ribo/bin/Rscript \
python3 -m src.te_calc.te_calculator \
    --ribo_dir data/processed/ribo_files \
    --output_dir data/processed \
    --cpm_cutoff 1 \
    --overall_cutoff 70
```

### 跳过 Stage 0（已有 count 矩阵时）/ Skip extraction

```bash
RSCRIPT_BIN=/home/xrx/miniconda3/envs/snakemake-ribo/bin/Rscript \
python3 -m src.te_calc.te_calculator \
    --skip_extract \
    --output_dir data/processed \
    --cpm_cutoff 1 \
    --overall_cutoff 70
```

### 重新生成图形 / Regenerate figures

```bash
python3 src/plots/generate_pub_style_figures.py \
    --output-dir reports/figures/pub_style
```

---

## 流水线架构 / Pipeline Architecture

```
data/raw/fastq/*.fastq.gz
        │
        ▼  Snakemake: cutadapt → Bowtie2 (rRNA filter) → STAR → featureCounts
data/processed/ribo_files/*.ribo
        │
        ▼  Stage 0: te_calculator.py — 按物种分组提取 raw counts
data/processed/ribo_raw_{species}.csv
data/processed/rnaseq_raw_{species}.csv
        │
        ▼  Stage 1: CPM 标准化 → 低表达过滤 → Dummy gene 合并
data/processed/ribo_paired_count_dummy_{species}.csv
data/processed/rna_paired_count_dummy_{species}.csv
        │
        ▼  Stage 2: TE.R — CLR → ILR → lm(RIBO_ilr ~ RNA_ilr) → ilr2clr(residuals)
data/processed/human_TE_sample_level_{species}.rda
        │
        ▼  Stage 3: postprocess → te_{species}.csv（各物种独立输出）
data/processed/te_human.csv  /  te_mouse.csv  /  te_arabidopsis.csv  /  te_yeast.csv
        │
        ▼  merge_cross_species_te(): outer join（缺失基因填 NaN，不填 0）
data/processed/te_results_final.csv
        │
        ▼  generate_pub_style_figures.py → PCA / Median TE / Violin / UMAP
reports/figures/pub_style/*.pdf
```

**关键设计**：各物种独立运行 Stage 0–3，最后 outer join 合并。跨物种缺失基因位置为 NaN（不填 0，填 0 会引入虚假低 TE 信号破坏成分数据结构）。

---

## 重要配置 / Configuration Notes

| 配置项 | 值 | 说明 |
|--------|-----|------|
| TE.R 并行核数 | `min(64, detectCores()-2)` | R socket 上限 128，64 核为安全上限 |
| CPM 阈值 | `1` | `--cpm_cutoff` |
| overall_cutoff | `70%` | 样本覆盖率阈值 `--overall_cutoff` |
| 跨物种合并策略 | outer join，缺失填 NaN | 不填 0，避免破坏成分数据结构 |
| RSCRIPT_BIN | 环境变量 | 指定 R 可执行路径，默认 `Rscript` |

---

## 已知问题 / Known Issues

- **6 个 GSE 上游处理未完成**：GSE116221、GSE119615、GSE48140、GSE48933、GSE65778、GSE79664 的 Ribo-seq 处理尚未完成，暂未纳入 TE 计算。
- ***C. elegans***：GSM1169554 仅有 Ribo-seq，无配对 RNA-seq，暂未纳入 TE 计算。
- **Yeast 数据冗余**：`all.ribo`、`GSM3561545.ribo`、`GSE125038_all.ribo` 三个文件内容相同，当前作为独立文件处理，后续需去重。

---

## 版本历史 / Changelog

### v2.1（2026-04-09）

- **Bug fix**：修复跨物种 TE 计算 bug——原 `pd.concat().fillna(0)` 导致 arabidopsis/yeast 所有样本 TE 方差塌缩至 ~1e-16（结构零问题）
- **架构重构**：实现各物种独立 TE 计算（Stage 0–3）+ outer-join 合并，每物种输出单独 `te_{species}.csv`
- **重算全物种 TE**：human、mouse、arabidopsis、yeast 全部重新计算并验证
- **新增跨物种 UMAP**：基于 8,452 个共享 gene symbol 的 TE-UMAP（human + mouse，n=29）
- **规范核心脚本注释**：`te_calculator.py`、`TE.R`、`generate_pub_style_figures.py` 注释清理与完善
- **清理中间数据**：删除 33 个已完成 GSE 的中间数据（释放约 3.25 TB）

### v2.0（之前）

- 原始多物种 TE 流水线上线
- 整合 snakescale 上游流水线（无 SQLite 版）
- Makefile 统一调度

---

## 参考文献 / References

- Liu Y *et al.* "Translation efficiency covariation across cell types is a conserved organizing principle of mammalian transcriptomes." CenikLab TE_model. [GitHub](https://github.com/CenikLab/TE_model)
- Cenik *et al.*, NBT 2025（方法论基础）
- Cookiecutter Data Science. [https://cookiecutter-data-science.drivendata.org/](https://cookiecutter-data-science.drivendata.org/)

---

## 许可与声明 / License

MIT License

TE 计算核心算法（`src/te_calc/TE.R` 及 `te_calculator.py` 中 Stage 1–2 的数学逻辑）100% 严格复刻自 CenikLab TE_model 原文代码，未做任何 AI 优化或篡改。学术引用与知识产权归属于原作者。

> The core TE algorithm in `src/te_calc/TE.R` and `te_calculator.py` (Stages 1–2) is a strict reproduction of CenikLab TE_model source code with no AI modifications. All academic attribution belongs to the original authors.

