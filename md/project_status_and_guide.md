# 项目现状与自动化管线使用指南

> **文档版本**: 2.0  
> **更新日期**: 2026-04-07  
> **适用项目**: Silk Protein Translation Efficiency Analysis Platform  
> **状态**: 四项关键系统级 Bug 已修复，管线进入可运行状态

---

## 目录

1. [项目整体架构与数据流转概览](#1-项目整体架构与数据流转概览)
2. [关键 Bug 修复总结](#2-关键-bug-修复总结)
3. [标准操作流程 (SOP)](#3-标准操作流程-sop)
4. [常见问题与故障排除](#4-常见问题与故障排除)
5. [依赖环境速查](#5-依赖环境速查)

---

## 1. 项目整体架构与数据流转概览

### 1.1 三段式架构

本项目是家蚕（*Bombyx mori*）及多物种丝蛋白基因翻译效率（TE）分析的端到端自动化平台，采用三段式架构：

```
┌─────────────────────────────────────────────────────────────────────┐
│  Stage 1: 数据注入                                                   │
│  make download_data                                                  │
│  SRA → sradownloader → data/raw/fastq/*.fastq.gz                    │
└────────────────────────────┬────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Stage 2: 上游比对定量 (Snakescale + RiboFlow)                       │
│  make run_snakescale                                                 │
│  Makefile → Snakemake DAG → Nextflow RiboFlow 引擎                   │
│    ├─ Ribo-seq: cutadapt → Bowtie2 rRNA过滤 → Bowtie2 转录组比对    │
│    │           → MAPQ过滤 → BAM→BED → 合并去重 → ribopy create      │
│    ├─ RNA-seq:  cutadapt → Bowtie2 rRNA过滤 → Bowtie2 转录组比对    │
│    │           → MAPQ过滤 → BAM→BED → 合并去重 → ribopy rnaseq set  │
│    └─ 合并: ribopy merge → all.ribo                                  │
└────────────────────────────┬────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────┐
│  Stage 3: 翻译效率计算 (TE Calculator)                               │
│  make calc_te                                                        │
│  Stage 0: .ribo → CSV 提取                                          │
│  Stage 1: CPM 标准化 → 低表达过滤 → Dummy Gene 合并 → 非polyA移除   │
│  Stage 2: CLR 变换 → ILR 投影 → 成分线性回归 (R)                    │
│  Stage 3: TE = ilr2clr(residuals) → te_results_final.csv            │
└─────────────────────────────────────────────────────────────────────┘
```

### 1.2 编排层级

| 层级 | 工具 | 职责 |
|------|------|------|
| **顶层调度** | `Makefile` (12 个目标) | 统一入口，参数传递，模式路由 |
| **DAG 编排** | Snakemake (`Snakefile`) | adapter 检测 → 读长统计 → YAML 生成 → RiboFlow 调用 |
| **核心引擎** | Nextflow (`RiboFlow.groovy`, 2053 行) | Ribo-seq / RNA-seq 双轨并行处理 + .ribo 封装 |
| **TE 计算** | Python (`te_calculator.py`) + R (`TE.R`) | CoDA 成分回归，100% 复刻 CenikLab TE_model |

### 1.3 关键数据流转路径

```
data/raw/fastq/*.fastq.gz           ← 原始 FASTQ (只读)
    ↓ cutadapt
data/interim/adapter_check/          ← 接头检测结果
    ↓ Bowtie2 (rRNA filter)
data/interim/filtered/               ← rRNA 过滤后 clean reads
    ↓ Bowtie2 (transcriptome)
data/interim/*_sorted.bam            ← 转录组比对 BAM
    ↓ bamToBed → merge → dedup
*.merged.post_dedup.bed              ← 去重后 BED
    ↓ ribopy create + rnaseq set + merge
data/processed/ribo_files/all.ribo   ← HDF5 容器 (Ribo+RNA 配对封装)
    ↓ te_calculator.py Stage 0
data/processed/ribo_raw.csv          ← Ribo-seq CDS 计数矩阵
data/processed/rnaseq_raw.csv        ← RNA-seq CDS 计数矩阵
    ↓ te_calculator.py Stage 1
data/processed/*_paired_count_dummy.csv  ← CPM 过滤 + Dummy Gene 处理
    ↓ TE.R (Stage 2)
data/processed/human_TE_sample_level.rda ← CLR/ILR 回归结果
    ↓ te_calculator.py Stage 3
data/processed/te_results_final.csv  ← ★ 最终 TE 结果
    ↓ plot_te_features.R
reports/figures/*.pdf                ← 可视化图表
```

### 1.4 Ribo↔RNA 配对保证机制

配对通过三层机制确保不会错配：

| 层级 | 位置 | 机制 |
|------|------|------|
| L1: 元数据驱动 | `generate_yaml.py` | `matched_RNA-seq_experiment_alias` 列建立 GSM↔GSM 映射 |
| L2: YAML 键对齐 | `project.yaml` | `input.fastq[GSM]` 和 `rnaseq.fastq[GSM]` 使用相同 key |
| L3: Nextflow Join | `RiboFlow.groovy` | Channel `.join()` 按 sample name 物理绑定 |

---

## 2. 关键 Bug 修复总结

以下四项系统级 Bug 已完成修复，管线稳定性和多物种支持能力显著提升。

### 2.1 Nextflow Trace 并发隔离

| 项目 | 详情 |
|------|------|
| **问题** | 当 Snakemake 并行调度多个 study 的 RiboFlow 实例时，所有实例共用同一个 `trace.txt` 文件，导致文件锁竞争、数据覆盖和运行崩溃 |
| **修复位置** | `workflow/snakescale/Snakefile` (第 714-716 行) |
| **修复方案** | 为每个 study 分配独立的 trace/report/timeline 文件路径 |

**修复后代码**:
```python
shell(
    "nextflow riboflow/RiboFlow.groovy -params-file {input.study_yaml} "
    "-profile {params.profile} "
    "-with-trace nextflow_logs/{wildcards.study}_trace.txt "
    "-with-report nextflow_logs/{wildcards.study}_report.html "
    "-with-timeline nextflow_logs/{wildcards.study}_timeline.html"
)
```

**效果**: 多 study 并行运行不再产生文件冲突，每个 study 的执行日志完整独立可追溯。

---

### 2.2 TE 计算列名智能映射

| 项目 | 详情 |
|------|------|
| **问题** | `.ribo` 文件提取的 count 矩阵列名可能是 SRR、GSM、SRX 或特殊值 `"all"`（合并后的 .ribo 文件），与下游配对逻辑期望的列名不匹配，导致 KeyError 或静默丢失样本 |
| **修复位置** | `src/te_calc/te_calculator.py` (第 376-487 行，`validate_and_align_columns()` 函数) |
| **修复方案** | 实现多候选智能匹配 + `"all"` 列名自动重映射 |

**核心逻辑**:
```
匹配优先级: SRR accession → GSM alias → SRX accession
特殊处理:   合并后 .ribo 的 "all" 列 → 自动映射到配对 GSM
```

**效果**: 无论上游输出列名为何种格式，TE 计算都能正确识别并配对 Ribo-seq / RNA-seq 样本。

---

### 2.3 线虫 (C. elegans) 参考基因组构建脚本

| 项目 | 详情 |
|------|------|
| **问题** | `references.yaml` 中线虫的 `regions` 和 `transcript_lengths` 路径为空，导致 Nextflow 文件暂存 (file staging) 时产生目录碰撞错误，管线无法处理线虫样本 |
| **修复位置** | `workflow/snakescale/scripts/build_celegans_ref.sh` (新建)；`workflow/snakescale/scripts/generate_yaml.py` (第 252-259 行) |
| **修复方案** | (1) 创建完整的线虫参考文件构建脚本；(2) 在 YAML 生成阶段增加参考路径空值校验 |

**构建脚本功能**:
1. 从 Ensembl Release 112 下载基因组 FASTA 和 GFF3 注释
2. 使用 gffread 提取 mRNA 序列
3. 构建 Bowtie2 转录组索引 (8 线程)
4. 生成 `regions.bed` (CDS/UTR 坐标) 和 `transcript_lengths.tsv`

**空值校验**:
```python
for ref_key in ["filter", "regions", "transcript_lengths", "transcriptome"]:
    ref_path = ref_contents.get(ref_key, "")
    if not ref_path or str(ref_path).strip() == "":
        raise RuntimeError(
            "物种 '{}' 的参考路径 '{}' 为空!\n"
            "请先构建该物种的参考文件并更新 scripts/references.yaml。\n"
            "对于 C. elegans, 可运行 scripts/build_celegans_ref.sh")
```

**效果**: 线虫参考路径缺失时提前报错并给出清晰修复指引，不再让管线运行到 Nextflow 阶段才崩溃。

> **注意**: 当前 `references.yaml` 中线虫的 `regions` 和 `transcript_lengths` 仍为空占位，需在首次处理线虫样本前运行 `build_celegans_ref.sh` 构建参考文件。

---

### 2.4 人类样本 Adapter 强制执行

| 项目 | 详情 |
|------|------|
| **问题** | 特定人类样本（如 GSM1211429/SRR953773）adapter 含量极低，`guess_adapter()` 自动检测失败，导致 Snakemake 的 `check_adapter` 规则阻断整条管线 |
| **修复位置** | `workflow/snakescale/config.yaml` (第 35 行)；`data/external/test_metadata.csv` (相关行) |
| **修复方案** | (1) 在 config 中启用 `override: True` 允许绕过检查失败；(2) 在元数据中显式注入已知 adapter 序列 |

**Config 修复**:
```yaml
override: True    # 允许 adapter/length 检查失败时继续执行
```

**Metadata 修复**: 为问题样本显式指定 Illumina 3' adapter:
```
GSM1211429: threep_adapter = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

**效果**: adapter 含量低的样本不再阻断管线，同时通过显式注入确保修剪参数正确。

---

## 3. 标准操作流程 (SOP)

### 3.1 环境准备

#### Python 环境
```bash
# 推荐使用 conda 隔离环境
conda create -n riboseq python=3.10
conda activate riboseq

# 安装 Python 依赖
pip install pandas numpy bioinfokit ribopy snakemake
```

#### R 环境
```r
# 在 R 终端中执行
install.packages(c("tidyverse", "foreach", "doParallel", "compositions"))
BiocManager::install("propr")
```

#### 系统工具
确保以下工具已安装并在 `$PATH` 中：
- `bowtie2`, `samtools`, `cutadapt`, `bamToBed` (bedtools)
- `nextflow`, `sradownloader`
- STAR (如使用基因组比对模式)

### 3.2 项目初始化

```bash
cd /home/xrx/my_project/project

# 1. 创建 CCDS 标准目录骨架
make init

# 2. (集群用户) 初始化 SLURM 投递配置
make setup_slurm
```

### 3.3 数据准备

#### 方式 A: 软链接服务器已有数据
```bash
make link_data
# 该命令将服务器归档数据软链接到 data/raw/ 目录
```

#### 方式 B: 从 SRA 自动下载
```bash
# 一键下载 (包含 SRX→SRR 映射 + FASTQ 下载)
make download_data

# 或分步执行:
make download_prepare   # 仅生成 SRX→SRR 映射表
make download_fetch     # 仅执行 FASTQ 下载
```

#### 元数据要求

确保 `data/external/metadata.csv` 包含以下必要列：

| 列名 | 说明 | 示例 |
|------|------|------|
| `experiment_alias` | GSM 样本 ID (主键) | GSM1211426 |
| `experiment_accession` | SRX 实验 ID | SRX327597 |
| `corrected_type` | 实验类型 | `Ribo-Seq` / `RNA-Seq` |
| `matched_RNA-seq_experiment_alias` | 配对 RNA-seq 的 GSM | GSM1211429 |
| `organism` | 物种名 | `Homo sapiens` |
| `threep_adapter` | 3' adapter 序列 | AGATCGGAAG... |

### 3.4 上游比对定量 — 全核心运行

#### 步骤 1: 干运行验证

```bash
make dry_run
# 检查 Snakemake DAG 是否解析成功
# 预期: 列出所有将执行的规则，无报错
# 注意: 缺少物理 FASTQ 文件时报 MissingInputException 为正常行为
```

#### 步骤 2: 生产模式运行 (SLURM 集群)

```bash
# ★ 推荐: 集群投递模式
make run_snakescale
```

该命令将：
1. `check_adapter` — 采样 25000 reads 检测 adapter 存在性
2. `check_lengths` — 统计读长分布
3. `generate_yaml` — 从 CSV 元数据生成 per-study 的 `project.yaml`
4. `run_riboflow` — 启动 Nextflow RiboFlow 引擎，每个 study 独立并发执行

#### 步骤 2 (备选): 本地调试模式

```bash
make run_local FORCE=1
# FORCE=1 用于绕过交互节点安全检查
# 适用于小数据集测试和调试
```

#### 步骤 3: 检查运行结果

运行完成后，检查以下输出：
```bash
# 检查 Nextflow 执行日志 (每个 study 独立)
ls project/workflow/snakescale/nextflow_logs/

# 检查 .ribo 输出文件
ls data/processed/ribo_files/

# 检查 count 矩阵 (如已提取)
ls data/processed/ribo_raw.csv data/processed/rnaseq_raw.csv
```

### 3.5 下游 TE 计算

TE 计算模块支持三种输入模式，自动检测并路由：

#### 模式 A: 从 .ribo 文件全流程计算 (推荐)

```bash
make calc_te RIBO_DIR=data/processed/ribo_files
```

执行流程：Stage 0 (.ribo→CSV) → Stage 1 (CPM/过滤/Dummy) → Stage 2 (CLR/ILR 回归) → Stage 3 (后处理)

#### 模式 B: 从已有 count 矩阵计算

```bash
# 确保 data/processed/ribo_raw.csv 和 rnaseq_raw.csv 已存在
make calc_te
```

跳过 Stage 0，直接进入 Stage 1 预处理。

#### 模式 C: 无上游数据时的管线测试

```bash
make calc_te
# 当无 .ribo 文件也无 count 矩阵时，自动生成 200 genes × 6 samples 假数据
# 用于验证 TE 计算管线完整性
```

#### 可选参数

```bash
# 调整低表达过滤阈值 (默认: CPM<1, 超过70%样本)
make calc_te CPM_CUTOFF=1 OVERALL_CUTOFF=70

# 移除非 polyA 基因
make calc_te NONPOLYA_CSV=data/external/nonpolyA_gene.csv

# 仅执行 Python 预处理，跳过 R 回归 (调试用)
make calc_te SKIP_R=1
```

#### 检查 TE 结果

```bash
# 最终输出
head data/processed/te_results_final.csv

# R 中间结果
ls data/processed/human_TE_sample_level.rda
```

### 3.6 可视化

```bash
make plot
# 输出: reports/figures/*.pdf
```

### 3.7 一键全流程

```bash
# 按顺序执行: init → link_data → download → snakescale → calc_te → plot
make all
```

### 3.8 清理与重置

```bash
make clean
# 清理 data/interim/, data/processed/, reports/figures/
# data/raw/ (原始数据) 不受影响，符合不可变原则
```

---

## 4. 常见问题与故障排除

### Q1: 线虫 (C. elegans) 样本报参考路径为空？

```bash
# 运行参考基因组构建脚本
cd workflow/snakescale/scripts
bash build_celegans_ref.sh

# 构建完成后，更新 references.yaml 中线虫的 regions 和 transcript_lengths 路径
```

### Q2: Adapter 检测失败导致管线中断？

两种解决方案：
1. 在 `config.yaml` 中设置 `override: True`（已默认启用）
2. 在 `metadata.csv` 的 `threep_adapter` 列中显式填写已知 adapter 序列

### Q3: TE 计算报列名不匹配 (KeyError)？

`validate_and_align_columns()` 已实现智能映射（SRR → GSM → SRX），正常情况下不会再出现此问题。如仍遇到，检查：
- `metadata.csv` 和 `srx_to_srr_mapping.csv` 中的样本 ID 是否一致
- `.ribo` 文件中的 experiment name 是否与元数据匹配

### Q4: 多 study 并行运行时 Nextflow 报文件冲突？

已修复。每个 study 的 trace/report/timeline 文件现在独立存放在 `nextflow_logs/{study}_trace.txt`。

### Q5: 如何添加新物种支持？

1. 在 `references.yaml` 中添加物种条目（需提供 filter, transcriptome, regions, transcript_lengths 四项路径）
2. 在 `generate_yaml.py` 的 `ORGANISM_ALIAS` 字典中添加物种名归一化映射
3. 构建对应的 Bowtie2 索引和区域注释文件

---

## 5. 依赖环境速查

### Python 包

| 包 | 版本要求 | 用途 |
|----|----------|------|
| pandas | ≥ 2.0 | 数据处理 |
| numpy | ≥ 1.24 | 数值计算 |
| bioinfokit | latest | CPM 标准化 |
| ribopy | latest | .ribo 文件解析 |
| snakemake | ≥ 7.0 | 上游流水线编排 |

### R 包

| 包 | 用途 |
|----|------|
| propr | CLR 变换与 rho 比例性度量 |
| compositions | CLR ↔ ILR 变换 |
| tidyverse | 数据整理与分组聚合 |
| foreach + doParallel | 基因级并行回归 |

### 系统工具

| 工具 | 用途 | 所在阶段 |
|------|------|----------|
| bowtie2 | rRNA 过滤 + 转录组比对 | Stage 2 |
| samtools | BAM 处理全流程 | Stage 2 |
| cutadapt | Adapter 修剪 | Stage 2 |
| bamToBed (bedtools) | BAM → BED 转换 | Stage 2 |
| Nextflow | RiboFlow 工作流引擎 | Stage 2 |
| sradownloader | SRA 数据下载 | Stage 1 |
| ribopy | .ribo 创建/合并/提取 | Stage 2-3 |

---

## Makefile 目标速查表

| 目标 | 功能 | 备注 |
|------|------|------|
| `make init` | 创建目录骨架 | 首次使用时运行 |
| `make setup_slurm` | 初始化 SLURM profile | 集群用户 |
| `make link_data` | 软链接已有数据 | — |
| `make download_data` | SRA 全流程下载 | 包含映射+下载 |
| `make download_prepare` | 仅生成 SRX→SRR 映射 | — |
| `make download_fetch` | 仅执行 FASTQ 下载 | — |
| `make dry_run` | Snakemake 干运行 | 验证 DAG |
| `make run_snakescale` | **生产模式: SLURM 投递** | ★ 主要入口 |
| `make run_local FORCE=1` | 本地调试模式 | 小数据测试 |
| `make calc_te` | **TE 计算** | 自动三模式路由 |
| `make plot` | 生成可视化图表 | → reports/figures/ |
| `make all` | 一键全流程 | 串联所有阶段 |
| `make clean` | 清理中间/最终产物 | raw/ 不受影响 |
| `make help` | 显示帮助信息 | — |

---

> **文档生成日期**: 2026-04-07  
> **基于源码审计**: `Snakefile` (722行), `RiboFlow.groovy` (2053行), `generate_yaml.py` (407行), `te_calculator.py` (1035行), `TE.R` (121行), `build_celegans_ref.sh` (85行)
