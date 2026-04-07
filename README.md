# Silk Protein Translation Efficiency Analysis Platform

家蚕（*Bombyx mori*）及蜘蛛丝腺细胞丝蛋白基因翻译效率（TE）分析的端到端自动化工程平台。

## V2.0 更新说明 (Release Notes)

- **并行稳定性增强**: 修复 Nextflow `trace.txt`/运行日志相关的并发写入竞争问题，提升多 Study 并行执行稳定性。
- **TE 结果路由增强**: 升级 TE 计算模块的列名映射逻辑，支持 `all` 容器与 GSM/SRR ID 的智能路由，不破坏原有 CoDA 数学逻辑。
- **线虫参考自动化**: 新增并跑通 *C. elegans* 自动化参考基因组构建与转录本坐标映射脚本。
- **Bowtie2 日志健壮性修复**: 修复 Bowtie2 日志解析在遇到 `[WARNING]` 行时的崩溃问题。
- **Adapter 容错增强**: 完善 Adapter 预检机制的容错策略，在低 adapter presence 或缺少显式 adapter 注释时仍可稳定生成参数。

## 项目背景

本项目旨在构建一条从原始测序数据到翻译效率定量结果的全自动化分析流水线，用于探究丝蛋白基因在丝腺细胞中实现高效翻译的分子机制。项目以 [Cookiecutter Data Science (CCDS)](https://cookiecutter-data-science.drivendata.org/) 为目录规范，整合三个核心模块：

1. **数据获取** — 基于 `sradownloader` 的 SRA 自动化下载与 SRX→SRR 映射归档
2. **上游比对定量** — 基于 `snakescale`（无 SQLite 版）的 Ribo-seq / RNA-seq QC、比对与计数
3. **下游 TE 计算** — 严格复刻 CenikLab TE_model 的成分数据分析（CoDA）方法，通过 CLR/ILR 变换与回归残差定义翻译效率

三段式架构通过顶层 `Makefile` 统一调度，支持 SLURM 集群投递，实现全链路可追溯的学术级可复现性。

## 目录结构

```
project/
│
├── Makefile                       ← 顶层统一调度器 (12 个目标)
├── setup_project.sh               ← CCDS 目录骨架初始化脚本
├── link_real_data.sh              ← 服务器归档数据软链接脚本
├── requirements.txt               ← Python 依赖声明
├── README.md
│
├── data/                          ← 数据层 (不可变原则: raw/ 只读)
│   ├── raw/                       ←   原始数据 (FASTQ, 基因组索引)
│   │   ├── fastq/                 ←     SRA 下载的 FASTQ 文件
│   │   └── genome/                ←     参考基因组
│   ├── interim/                   ←   中间产物 (BAM, adapter check)
│   ├── processed/                 ←   最终输出 (Count 矩阵, TE 结果)
│   └── external/                  ←   外部辅助数据 (SRR 映射表等)
│
├── workflow/                      ← 流水线层
│   └── snakescale/                ←   Snakemake 上游流水线
│       ├── Snakefile              ←     规则定义 (QC→Trim→Align→Count)
│       ├── config.yaml            ←     全局参数与 SLURM 配置
│       ├── scripts/               ←     辅助脚本
│       │   ├── generate_yaml.py   ←       CSV 驱动的元数据→YAML 转换
│       │   ├── references.yaml    ←       21+物种 Bowtie2 索引路由表
│       │   ├── guess_adapters.py  ←       接头序列自动检测
│       │   └── download_reference.py ←    参考序列下载
│       └── schemas/               ←     配置校验 schema
│
├── src/                           ← 源代码层
│   ├── te_calc/                   ←   TE 计算模块 (核心)
│   │   ├── te_calculator.py       ←     Python 预处理 + 调度 (661 行)
│   │   └── TE.R                   ←     R CLR/ILR 成分回归 (原样复刻)
│   ├── data/                      ←   数据获取工具
│   │   ├── download_sra.py        ←     SRA 下载与 SRX→SRR 映射
│   │   └── sradownloader/         ←     sradownloader 集成
│   ├── plots/                     ←   可视化脚本
│   │   └── plot_te_features.R     ←     TE 特征分布图 (R)
│   ├── features/                  ←   特征工程 (预留)
│   └── models/                    ←   模型代码 (预留)
│
├── reports/
│   └── figures/                   ← 生成的图表输出
│
├── notebooks/                     ← 探索性分析 Jupyter Notebooks
├── models/                        ← 序列化/训练好的模型
└── references/                    ← 参考文献、数据字典与操作手册
```

**不可变原则**：`data/raw/` 目录为只读区，一旦写入不可修改；所有中间结果和最终输出分别存放于 `data/interim/` 和 `data/processed/`，可通过 `make clean` 安全清除并重新生成。

## 核心模块工作流

### 三段式架构概览

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                    Stage 1: 数据注入                             │
 │  make download_data                                             │
 │  SRA → sradownloader → data/raw/fastq/*.fastq.gz               │
 └────────────────────────────┬────────────────────────────────────┘
                              │
                              ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │                    Stage 2: 上游比对定量                         │
 │  make run_snakescale                                            │
 │  FASTQ → snakescale (QC → Trim → Bowtie2 rRNA Filter           │
 │         → STAR Align → featureCounts) → Count 矩阵              │
 │  物种路由: references.yaml (21 种, 无 fallback 硬中断)           │
 │  元数据: CSV 驱动 (已脱离 SQLite 依赖)                           │
 └────────────────────────────┬────────────────────────────────────┘
                              │
                              ▼
 ┌─────────────────────────────────────────────────────────────────┐
 │                    Stage 3: 翻译效率计算                         │
 │  make calc_te                                                   │
 │  Count 矩阵 → CPM 标准化 → 低表达基因过滤 (CPM<1, >70%样本)     │
 │  → Dummy gene 合并 (成分闭合) → 非 polyA 基因移除                │
 │  → CLR 变换 → ILR 投影 → 成分线性回归                           │
 │  → TE = ilr2clr(residuals) → te_results_final.csv              │
 │  核心算法: 100% 复刻 CenikLab TE_model (CoDA 方法)              │
 └─────────────────────────────────────────────────────────────────┘
```

### 模块 1: sradownloader 数据注入

- 读取 `data/external/metadata.csv` 中的 SRX 样本编号
- 通过 NCBI E-utilities 自动查询并生成 `srx_to_srr_mapping.csv`
- 调用 `sradownloader` 批量下载 FASTQ 至 `data/raw/fastq/`
- 触发命令: `make download_data`

### 模块 2: snakescale 上游流水线（无 SQLite 版）

- 基于 [RiboBase/snakescale](https://github.com/CenikLab/snakescale) 重构
- **关键改造**: 剥离 SQLite 数据库依赖，改为纯 CSV 驱动的元数据解析
- **物种索引路由**: `references.yaml` 覆盖 21 个已建索引物种 + 2 个业务目标占位（家蚕、蜘蛛）
- **安全机制**: 未知物种触发硬异常中断，杜绝静默 fallback 至人类参考基因组的误路由
- 支持 SLURM 集群投递（`make run_snakescale`）与本地调试（`make run_local`）

### 模块 3: CenikLab TE 计算（CoDA 成分回归）

核心算法严格复刻自参考文献:

> Liu *et al.*, *Translation efficiency covariation across cell types is a conserved organizing principle of mammalian transcriptomes*, CenikLab TE_model

**数学流程**:

| 步骤 | 操作 | 实现 |
|------|------|------|
| 1 | CPM 标准化 | `bioinfokit.analys.norm().cpm()` |
| 2 | 低表达基因过滤 | CPM < 1 的样本数 > 70% × n_samples → dummy |
| 3 | Dummy gene 合并 | Ribo ∪ RNA dummy 基因 count 求和，保持成分闭合 |
| 4 | CLR 变换 | `propr(data, metric="rho", ivar="clr", alpha=NA, p=100)` |
| 5 | ILR 投影 | `compositions::clr2ilr()` |
| 6 | 成分线性回归 | `lm(RIBO_ilr[,i] ~ RNA_ilr[,i])` 对每个基因 |
| 7 | TE 定义 | `ilr2clr(residuals)` — 残差即为翻译效率 |

零值处理由 `propr` R 包内部执行乘性简单替换（multiplicative simple replacement），无需外部伪计数。

## 快速启动 (Quick Start)

### 环境准备

```bash
# Python 依赖
pip install pandas numpy bioinfokit ribopy

# R 依赖 (在 R 终端中执行)
install.packages(c("tidyverse", "foreach", "doParallel"))
# propr 和 compositions 需从 Bioconductor/CRAN 安装
BiocManager::install("propr")
install.packages("compositions")

# Snakemake (推荐 conda 环境)
conda activate snakescale
```

### 初始化

```bash
make init              # 创建 CCDS 目录骨架
make setup_slurm       # 初始化 SLURM profile (集群用户)
```

### 数据准备

```bash
make link_data              # 软链接服务器已有数据到 data/raw/
# 或
make download_data          # 从 SRA 自动下载 FASTQ
make download_prepare       # (仅生成 SRX→SRR 映射，不下载)
make download_fetch         # (仅执行下载)
```

### 上游比对定量

```bash
make dry_run               # 干运行: 验证 Snakemake DAG (不执行)
make run_snakescale        # ★ 生产模式: SLURM 集群投递
make run_local FORCE=1     # 本地调试 (交互节点, 需 FORCE=1 绕过安全检查)
```

### 翻译效率计算

```bash
# 模式 A: 从 .ribo 文件执行全流程 (Stage 0-3)
make calc_te RIBO_DIR=data/processed/ribo_files

# 模式 B: 已有 count 矩阵 (自动检测 data/processed/ribo_raw.csv)
make calc_te

# 模式 C: 无上游数据时自动生成假数据 (管线测试)
make calc_te

# 可选参数
make calc_te CPM_CUTOFF=1 OVERALL_CUTOFF=70     # 调整过滤阈值
make calc_te NONPOLYA_CSV=data/external/nonpolyA_gene.csv  # 移除非 polyA 基因
make calc_te SKIP_R=1                             # 仅预处理, 跳过 R 脚本
```

### 可视化

```bash
make plot              # 生成 TE 特征分布图 → reports/figures/
```

### 一键全流程 & 维护

```bash
make all               # 按顺序执行: init → link_data → download → snakescale → calc_te → plot
make clean             # 清理 interim/, processed/, figures/ (raw/ 不受影响)
make help              # 显示所有可用目标及参数说明
```

## 数据流转

```
data/raw/fastq/*.fastq.gz
        │
        ▼  Snakemake: cutadapt (adapter trimming)
data/interim/adapter_check/
        │
        ▼  Snakemake: Bowtie2 (rRNA filtering)
data/interim/filtered/
        │
        ▼  Snakemake: STAR (genome alignment)
data/interim/*_sorted.bam
        │
        ▼  Snakemake: featureCounts / .ribo generation
data/processed/ribo_raw.csv + rnaseq_raw.csv
        │
        ▼  Python: te_calculator.py (CPM → Filter → Dummy gene)
data/processed/ribo_paired_count_dummy.csv
        │
        ▼  R: TE.R (CLR → ILR → lm() → residuals)
data/processed/human_TE_sample_level.rda
        │
        ▼  Python: postprocess (transpose, group by cell_line)
data/processed/te_results_final.csv
        │
        ▼  R: plot_te_features.R
reports/figures/*.pdf
```

## 测试验证记录

| 测试项 | 命令 | 结果 | 备注 |
|--------|------|------|------|
| Snakemake DAG 构建 | `make dry_run` | ✅ DAG 解析成功 | 缺少物理 FASTQ 时报 MissingInputException（预期行为） |
| TE Dummy 模式 | `make calc_te` | ✅ 200 genes × 6 samples | 输出 `data/processed/te_results_final.csv` |
| Makefile 语法 | `make -n calc_te` | ✅ 三模式路由正确 | 自动识别 A/B/C 模式 |
| 物种路由安全 | 未注册物种测试 | ✅ 硬异常中断 | 不再静默 fallback 至 human |

## 依赖清单

### Python

| 包 | 版本 | 用途 |
|----|------|------|
| pandas | ≥ 2.0 | 数据处理 |
| numpy | ≥ 1.24 | 数值计算 |
| bioinfokit | latest | CPM 标准化 |
| ribopy | latest | .ribo 文件解析 (Stage 0) |
| snakemake | ≥ 7.0 | 上游流水线 |

### R

| 包 | 用途 |
|----|------|
| propr | CLR 变换与 rho 比例性度量 |
| compositions | CLR ↔ ILR 变换 |
| tidyverse | 数据整理与分组聚合 |
| foreach + doParallel | 基因级并行回归 |

### 系统工具

| 工具 | 用途 |
|------|------|
| Bowtie2 | rRNA 过滤比对 |
| STAR | 基因组比对 |
| cutadapt | 接头修剪 |
| samtools | BAM 处理 |
| sradownloader | SRA 数据下载 |

## 参考文献

- Liu Y *et al.* Translation efficiency covariation across cell types is a conserved organizing principle of mammalian transcriptomes. *CenikLab TE_model*. [GitHub](https://github.com/CenikLab/TE_model)
- Cookiecutter Data Science. [https://cookiecutter-data-science.drivendata.org/](https://cookiecutter-data-science.drivendata.org/)
- RiboBase/snakescale. [https://github.com/CenikLab/snakescale](https://github.com/CenikLab/snakescale)

## 许可与声明

TE 计算核心算法（`src/te_calc/TE.R` 及 `te_calculator.py` 中 Stage 1–2 的数学逻辑）100% 严格复刻自 CenikLab TE_model 原文代码，未做任何 AI 优化或篡改。学术引用与知识产权归属于原作者。
