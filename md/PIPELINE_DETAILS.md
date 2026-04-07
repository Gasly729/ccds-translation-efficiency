# 端到端数据流转与运算逻辑白皮书

## Pipeline End-to-End Data Flow & Computational Logic

> **文档版本**: 1.0  
> **生成日期**: 2026-03-25  
> **适用项目**: CCDS Ribo-seq / RNA-seq 翻译效率分析平台  
> **源码依据**: `workflow/snakescale/` (Snakemake + Nextflow) + `src/te_calc/` (Python + R)

---

## 目录

1. [架构总览](#1-架构总览)
2. [阶段 A — 数据获取与元数据准备](#2-阶段-a--数据获取与元数据准备)
3. [阶段 B — Snakescale 编排层](#3-阶段-b--snakescale-编排层)
4. [阶段 C — RiboFlow 核心处理引擎 (Ribo-seq)](#4-阶段-c--riboflow-核心处理引擎-ribo-seq)
5. [阶段 D — RiboFlow RNA-seq 并行处理](#5-阶段-d--riboflow-rna-seq-并行处理)
6. [阶段 E — `.ribo` 容器封装 (核心)](#6-阶段-e--ribo-容器封装-核心)
7. [阶段 F — TE 计算](#7-阶段-f--te-计算)
8. [关键发现：Ribo↔RNA 配对机制](#8-关键发现riborn-配对机制)
9. [完整数据流图](#9-完整数据流图)
10. [工具与版本清单](#10-工具与版本清单)

---

## 1. 架构总览

整条流水线分为三个编排层级和六个执行阶段：

```
┌─────────────────────────────────────────────────────────────────────┐
│  Makefile (顶层入口)                                                │
│    make download_data → make run_snakescale → make calc_te          │
├─────────────────────────────────────────────────────────────────────┤
│  Snakescale (Snakemake DAG)                                         │
│    rule check_adapter → rule check_lengths → rule generate_yaml     │
│    → rule run_riboflow                                              │
├─────────────────────────────────────────────────────────────────────┤
│  RiboFlow (Nextflow Engine: RiboFlow.groovy)                        │
│    clip → filter → transcriptome_alignment → quality_filter         │
│    → bam_to_bed → deduplicate → create_ribo → put_rnaseq_into_ribo │
│    → merge_ribos                                                    │
├─────────────────────────────────────────────────────────────────────┤
│  TE Calculator (Python + R)                                         │
│    Stage 0: .ribo → CSV  |  Stage 1: CPM/Filter/Dummy              │
│    Stage 2: CLR/ILR 回归 |  Stage 3: 后处理                        │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 2. 阶段 A — 数据获取与元数据准备

### 2.1 SRA 数据下载

| 项目 | 详情 |
|------|------|
| **脚本** | `src/data/download_sra.py` |
| **入口命令** | `make download_data` |
| **输入** | `data/external/metadata/metadata.csv` (只读存档) |
| **输出** | `data/raw/fastq/*.fastq.gz`, `data/external/srx_to_srr_mapping.csv` |

**执行逻辑**：

1. **解析元数据**: 从 `metadata.csv` (header=1) 提取所有 `experiment_accession` (SRX/ERX) 列
2. **Entrez API 批量映射**: 调用 NCBI Entrez E-utilities，将 SRX → SRR 批量解析
   - 批次大小: 200 IDs/请求
   - 速率限制: 无 API key 3次/秒, 有 API key 10次/秒
3. **生成映射表**: 输出 `srx_to_srr_mapping.csv`，包含列:
   - `Run` (SRR accession)
   - `srx` (SRX accession)
   - `gsm` (GEO Sample accession)
   - `organism`, `corrected_type` (Ribo-Seq / RNA-Seq), `study_name`
4. **调用 sradownloader**: 生成 `sradownloader_input.txt`，执行:
   ```bash
   sradownloader --outdir data/raw/fastq < sradownloader_input.txt
   ```

### 2.2 元数据结构

**`metadata.csv`** (header 在第 2 行, pandas `header=1`):

| 关键列 | 用途 |
|--------|------|
| `experiment_alias` | GEO Sample ID (GSM), 作为样本主键 |
| `experiment_accession` | SRA Experiment ID (SRX) |
| `corrected_type` | 实验类型: `Ribo-Seq` 或 `RNA-Seq` |
| `matched_RNA-seq_experiment_alias` | **配对锚点**: Ribo-Seq 样本对应的 RNA-Seq GSM |
| `organism` | 物种名 (需归一化) |
| `threep_adapter` / `fivep_adapter` | 3'/5' adapter 序列 |

**`srx_to_srr_mapping.csv`**:

| 关键列 | 用途 |
|--------|------|
| `Run` | SRR accession (FASTQ 文件名前缀) |
| `srx` | SRX accession |
| `gsm` | GSM accession (与 metadata 的 join key) |
| `study_name` | GEO Series (GSE) |
| `corrected_type` | Ribo-Seq / RNA-Seq |

---

## 3. 阶段 B — Snakescale 编排层

| 项目 | 详情 |
|------|------|
| **主文件** | `workflow/snakescale/Snakefile` |
| **配置** | `workflow/snakescale/config.yaml` |
| **入口命令** | `make run_snakescale` |

### 3.1 Snakemake DAG 规则链

```
rule check_adapter     检测 adapter 存在性 (采样 25000 reads)
       ↓
rule check_lengths     统计读长分布
       ↓
rule generate_yaml     ★ 从 CSV 生成 per-study project.yaml
       ↓
rule run_riboflow      执行 Nextflow RiboFlow 引擎
```

### 3.2 generate_yaml.py — YAML 生成器 (关键节点)

**脚本**: `workflow/snakescale/scripts/generate_yaml.py`

此脚本是 **Ribo↔RNA 配对** 的第一个锚定点。它从 CSV 元数据构建 RiboFlow 所需的 per-study YAML：

**执行逻辑**：

1. **按 study 过滤**: 从 `srx_to_srr_mapping.csv` 筛选 `study_name == GSE_ID`
2. **区分实验类型**: `corrected_type == 'Ribo-Seq'` vs `'RNA-Seq'`
3. **合并 metadata**: 通过 `gsm ↔ experiment_alias` join 获取 adapter 和配对信息
4. **物种归一化**: 通过 `ORGANISM_ALIAS` 字典将不一致物种名映射到 `references.yaml` 键
5. **构建 FASTQ 路径字典**:
   ```python
   # Ribo-seq: GSM → [SRR_1.fastq.gz, ...]
   ribo_fastq_dict[gsm] = [download_path/SRR_1.fastq.gz, ...]
   
   # RNA-seq: 通过 matched_RNA-seq_experiment_alias 映射
   # ribo_gsm → matched_rna_gsm → rna_fastq_paths
   rna_fastq_dict[ribo_gsm] = rna_fastq_by_gsm[matched_rna_gsm]
   ```
6. **生成 adapter 参数**:
   - Ribo-Seq: `-u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28 -a <ADAPTER> --overlap=4 --trimmed-only`
   - RNA-Seq: `-u 5 -l 40 --quality-cutoff=28 -a <ADAPTER> --overlap=4`
7. **Reference 路径注入**: 从 `references.yaml` 读取物种对应的 filter / transcriptome / regions / transcript_lengths 路径
8. **输出 YAML**: `data/interim/project/{GSE}/{study}.yaml`

**★ 配对关键点**: YAML 中 `input.fastq` 和 `rnaseq.fastq` 使用 **相同的 GSM 作为 key**。这意味着在下游 RiboFlow 中，同一个 sample name 下既有 Ribo-seq FASTQ 又有 RNA-seq FASTQ，从而实现了 **天然的样本名配对**。

### 3.3 RiboFlow 调用

```bash
nextflow riboflow/RiboFlow.groovy \
    -params-file {study}.yaml \
    -profile stampede_local
```

Docker 容器: `ceniklab/riboflow`

---

## 4. 阶段 C — RiboFlow 核心处理引擎 (Ribo-seq)

**主脚本**: `workflow/snakescale/riboflow/RiboFlow.groovy` (Nextflow DSL1, 2053 行)

### 4.1 Step C1: Adapter 修剪 (`process clip`)

| 项目 | 详情 |
|------|------|
| **工具** | `cutadapt` |
| **输入** | `{sample}.{index}.fastq.gz` |
| **输出** | `{sample}.{index}.clipped.fastq.gz` |

**命令**:
```bash
cutadapt --cores=${task.cpus} ${params.clip_arguments} ${fastq} \
    | gzip -c > ${sample}.${index}.clipped.fastq.gz
```

**典型参数** (来自 `generate_yaml.py`):
```
-u 1                    # 去除 5' 端 1 个碱基
--maximum-length=40     # 丢弃 >40bp 的读段
--minimum-length=15     # 丢弃 <15bp 的读段
--quality-cutoff=28     # 3' 端质量修剪 (Q≥28)
-a AGATCGGAAG...        # 3' adapter 序列
--overlap=4             # adapter 最小匹配长度
--trimmed-only          # (仅 Ribo-Seq) 只保留成功修剪的读段
```

**生物信息学意义**: Ribo-seq 读段来自核糖体足迹 (~28-32nt)，修剪后的读长范围 15-40bp 覆盖了典型 RPF 长度。`--trimmed-only` 确保只保留确实包含 adapter 的读段（即完整穿过 insert 的读段），这对于短片段文库至关重要。

### 4.2 Step C2: rRNA/tRNA 过滤 (`process filter`)

| 项目 | 详情 |
|------|------|
| **工具** | `bowtie2` + `samtools` |
| **输入** | `{sample}.{index}.clipped.fastq.gz` |
| **关键输出** | `{sample}.{index}.unaligned.filter.fastq.gz` (下游使用) |
| **参考索引** | `reference/filter/{species}/{species}_rtRNA*` (Bowtie2 索引) |

**命令**:
```bash
bowtie2 ${params.alignment_arguments.filter} \
        -x ${filter_index} -q ${fastq} \
        --threads ${task.cpus} \
        --al-gz ${sample}.aligned.filter.fastq.gz \
        --un-gz ${sample}.unaligned.filter.fastq.gz \
        2> ${sample}.filter.log \
    | samtools view -bS - \
    | samtools sort -o ${sample}.filter.bam \
    && samtools index ${sample}.filter.bam \
    && samtools idxstats ${sample}.filter.bam > ${sample}.filter.stats
```

**运算逻辑**:
- 将 clipped reads 比对到 rRNA/tRNA 参考序列
- **核心操作**: 取 `--un-gz` (未比对的读段) 作为下游输入
- 比对上 rRNA/tRNA 的读段被视为污染，丢弃
- 输出 BAM + idxstats 用于统计 rRNA 占比

### 4.3 Step C3: 转录组比对 (`process transcriptome_alignment`)

| 项目 | 详情 |
|------|------|
| **工具** | `bowtie2` + `samtools` |
| **输入** | `{sample}.{index}.unaligned.filter.fastq.gz` (过滤后的 clean reads) |
| **输出** | `{sample}.{index}.transcriptome_alignment.bam` |
| **参考索引** | `reference/transcriptome/{species}/appris_{species}_v2_selected*` |

**命令**:
```bash
bowtie2 ${params.alignment_arguments.transcriptome} \
        -x ${transcriptome_reference} -q ${fastq} \
        --threads ${task.cpus} \
        --al-gz ${sample}.aligned.transcriptome_alignment.fastq.gz \
        --un-gz ${sample}.unaligned.transcriptome_alignment.fastq.gz \
        2> ${sample}.transcriptome_alignment.log \
    | samtools view -bS - \
    | samtools sort -o ${sample}.transcriptome_alignment.bam \
    && samtools index ${sample}.transcriptome_alignment.bam \
    && samtools idxstats ${sample}.transcriptome_alignment.bam > ${sample}.transcriptome_alignment.stats
```

**关键点**:
- 比对到 **转录组** (非基因组)，参考序列为 APPRIS 筛选的代表性转录本
- 这意味着比对坐标是 **转录本坐标**，直接对应 CDS/UTR 区域
- 后续通过 `regions.bed` 文件标注 5'UTR / CDS / 3'UTR 边界

### 4.4 Step C4: MAPQ 过滤 (`process quality_filter`)

| 项目 | 详情 |
|------|------|
| **工具** | `samtools view -q` |
| **输入** | 转录组比对 BAM |
| **输出** | `{sample}.{index}.transcriptome_alignment.qpass.bam` |

**命令**:
```bash
samtools view -b -q ${params.mapping_quality_cutoff} ${bam} \
    | samtools sort -o ${sample}.transcriptome_alignment.qpass.bam
samtools view -b -c ${sample}.transcriptome_alignment.qpass.bam > ${sample}.qpass.count
samtools index ${sample}.transcriptome_alignment.qpass.bam
```

**参数**: `mapping_quality_cutoff` (典型值: 2, 来自 `project.yaml`)

### 4.5 Step C5: BAM → BED 转换 (`process bam_to_bed`)

| 项目 | 详情 |
|------|------|
| **工具** | `bamToBed` (bedtools) |
| **输入** | quality-filtered BAM |
| **输出** | `{sample}.{index}.bed` |

**命令**:
```bash
if [ $(samtools view -c ${bam}) -eq 0 ]; then
    touch ${sample}.${index}.bed
else
    bamToBed -i ${bam} > ${sample}.${index}.bed
fi
wc -l ${sample}.${index}.bed > ${sample}.${index}_nodedup_count.txt
```

**输出格式**: BED6 (chrom, start, end, name, score, strand)
- 每一行代表一个比对到转录组的读段位置

### 4.6 Step C6: 多 Lane 合并与去重

**合并** (`process merge_bed`):
```bash
# 先给每行添加 lane index 列 (用于去重后分离)
awk -v newcol=${sample}.${index} '{print($0"\t"newcol)}' ${bed}
# 合并排序
cat ${bed_files} | sort -k1,1 -k2,2n -k3,3n > ${sample}.merged.pre_dedup.bed
```

**去重** (`process deduplicate`, 可选, 由 `params.deduplicate` 控制):
```bash
rfc dedup -i ${bed} -o ${sample}.merged.post_dedup.bed
```

`rfc` 是 RiboFlow 自带的命令行工具集 (riboflow-commands)。

---

## 5. 阶段 D — RiboFlow RNA-seq 并行处理

RNA-seq 数据在 RiboFlow 内部经过 **完全独立但结构相同** 的处理管线：

```
rnaseq_clip → rnaseq_filter → rnaseq_transcriptome_alignment
    → rnaseq_quality_filter → rnaseq_bam_to_bed
    → rnaseq_merge_bed → (rnaseq_deduplicate)
```

### 关键差异

| 步骤 | Ribo-seq | RNA-seq |
|------|----------|---------|
| **cutadapt 参数** | `-u 1 --maximum-length=40 --minimum-length=15 --trimmed-only` | `-u 5 -l 40` (无 `--trimmed-only`) |
| **Adapter 策略** | 必须有 adapter (短片段) | 可选 adapter |
| **去重** | `params.deduplicate` | `params.rnaseq.deduplicate` (独立控制) |
| **Bowtie2 参数** | `params.alignment_arguments.transcriptome` | `params.rnaseq.bt2_arguments` |

### 触发条件

RNA-seq 处理仅在以下条件 **同时满足** 时执行:
```groovy
do_rnaseq = params.get("do_rnaseq", false) && params.get("rnaseq", false)
```
即 `project.yaml` 中必须同时设置 `do_rnaseq: true` 和 `rnaseq:` 节点。

---

## 6. 阶段 E — `.ribo` 容器封装 (核心)

这是本文档的 **核心发现**。`.ribo` 文件是基于 HDF5 的二进制容器格式，由 `ribopy` 库创建和操作。

### 6.1 Step E1: 创建 .ribo 文件 (`process create_ribo`)

| 项目 | 详情 |
|------|------|
| **工具** | `ribopy create` |
| **输入** | Ribo-seq BED 文件, 转录本长度文件, 区域注释文件, (可选) 元数据 YAML |
| **输出** | `{sample}.ribo` |

**完整命令** (来自 `RiboFlow.groovy:1093-1106`):
```bash
ribopy create \
    -n ${sample} \
    --reference ${params.ribo.ref_name} \
    --lengths ${transcript_length_file} \
    --annotation ${annotation_file} \
    --radius ${params.ribo.metagene_radius} \
    -l ${params.ribo.left_span} \
    -r ${params.ribo.right_span} \
    --lengthmin ${params.ribo.read_length.min} \
    --lengthmax ${params.ribo.read_length.max} \
    ${sample_meta_argument} \
    ${root_meta_argument} \
    ${coverage_argument} \
    -n ${task.cpus} \
    --alignmentfile ${bed_file} \
    ${sample}.ribo
```

**参数说明**:

| 参数 | 含义 | 典型值 |
|------|------|--------|
| `--reference` | 参考基因组/转录组名称 | `appris_human_v2` |
| `--lengths` | 转录本长度 TSV | `appris_human_v2_transcript_lengths.tsv` |
| `--annotation` | CDS/UTR 区域 BED | `appris_human_v2_actual_regions.bed` |
| `--radius` | metagene 分析半径 | 通常 50 |
| `-l` / `-r` | metagene 左/右跨度 | 通常 35 / 15 |
| `--lengthmin` / `--lengthmax` | 有效读长范围 | 通常 28 / 35 |
| `--alignmentfile` | Ribo-seq BED 比对文件 | `{sample}.merged.post_dedup.bed` |

**内部运算**:
1. 解析 BED 文件中每条比对记录的 **5' 端位置** (作为 A-site 偏移的基准)
2. 按读长 (lengthmin ~ lengthmax) 分层统计每个转录本上的覆盖度
3. 根据 `annotation.bed` 将转录本坐标分配到 **5'UTR / CDS / 3'UTR** 区域
4. 计算 metagene profile (起始密码子和终止密码子周围的平均覆盖度)
5. 将所有数据写入 HDF5 格式的 `.ribo` 文件

**HDF5 内部结构**:
```
{sample}.ribo (HDF5)
├── experiments/
│   └── {sample}/
│       ├── region_counts/     # 每个区域 (5UTR/CDS/3UTR) 的原始计数
│       ├── metagene/          # 起始/终止密码子周围的 metagene profile
│       ├── coverage/          # (可选) 逐碱基覆盖度
│       └── rnaseq/            # (后续注入) RNA-seq 计数
├── reference/
│   ├── transcript_lengths     # 转录本长度
│   └── annotation             # 区域注释
└── metadata/                  # 样本/实验元数据
```

### 6.2 Step E2: 注入 RNA-seq 数据 (`process put_rnaseq_into_ribo`)

> **★ 这是 RNA-seq 数据进入 .ribo 文件的唯一入口 ★**

| 项目 | 详情 |
|------|------|
| **工具** | `ribopy rnaseq set` |
| **输入** | 已创建的 `{sample}.ribo` + RNA-seq BED 文件 |
| **输出** | 同一个 `{sample}.ribo` (原地修改) |

**命令** (来自 `RiboFlow.groovy:1998`):
```bash
ribopy rnaseq set -n ${sample} -a ${rnaseq} -f bed --force ${ribo}
```

**参数说明**:
| 参数 | 含义 |
|------|------|
| `-n ${sample}` | 实验名称 (必须与 create 时一致) |
| `-a ${rnaseq}` | RNA-seq BED 比对文件路径 |
| `-f bed` | 输入格式为 BED |
| `--force` | 覆盖已有 RNA-seq 数据 |
| `${ribo}` | 目标 .ribo 文件 |

**内部运算**:
1. 解析 RNA-seq BED 文件，统计每个转录本的读段计数
2. 根据 `.ribo` 文件内已有的区域注释，将计数分配到 5'UTR / CDS / 3'UTR
3. 将 RNA-seq 计数矩阵写入 HDF5 的 `experiments/{sample}/rnaseq/` 节点

**配对逻辑** (来自 `RiboFlow.groovy:1974-2006`):
```groovy
// 通过 Nextflow channel join 操作，按 sample name 匹配
RIBO_FOR_RNASEQ                    // [sample, ribo_file]
    .join(RNASEQ_BED_FOR_RIBO_FINAL, remainder: true)  // [sample, bed_file]
    .choice(RNASEQ_FOR_RIBOPY, RIBO_FOR_RNASEQ_EXCLUDED)
    { it[2] != null ? 0 : 1 }     // 有 RNA-seq → 注入; 无 → 跳过
```

### 6.3 Step E3: 合并所有样本 (`process merge_ribos`)

| 项目 | 详情 |
|------|------|
| **工具** | `ribopy merge` |
| **输入** | 所有 per-sample `.ribo` 文件 |
| **输出** | `all.ribo` |

**命令**:
```bash
# 多样本: 合并
ribopy merge all.ribo ${sample_ribo_files}

# 单样本: 符号链接
ln -s ${sample_ribo} all.ribo
```

---

## 7. 阶段 F — TE 计算

| 项目 | 详情 |
|------|------|
| **脚本** | `src/te_calc/te_calculator.py` + `src/te_calc/TE.R` |
| **入口命令** | `make calc_te RIBO_DIR=data/processed/ribo_files` |
| **输入** | `.ribo` 文件 或 预提取的 `ribo_raw.csv` / `rnaseq_raw.csv` |
| **最终输出** | `data/processed/te_results_final.csv` |

### 7.1 Stage 0: .ribo → CSV 提取

**函数**: `stage0_extract()`

```python
ribo = Ribo(ribo_path)
# Ribo-seq CDS 计数
ribo_counts = ribo.get_region_counts(region_name="CDS",
                                      sum_lengths=True,
                                      sum_references=False)
# RNA-seq CDS 计数 (如果存在)
if ribo.has_rnaseq(experiment_name):
    rnaseq_counts = ribo.get_rnaseq()  # 提取 CDS RNA-seq 计数
```

**运算**:
- `sum_lengths=True`: 将所有读长的计数求和 (不区分 28nt / 29nt / ... / 35nt)
- `sum_references=False`: 保留每个转录本的独立计数
- 基因名提取: APPRIS 格式转录本 ID 的第 5 个 `|` 分隔字段
- 同一基因多个转录本: 取均值后取整

**输出**: `ribo_raw.csv`, `rnaseq_raw.csv` (genes × samples 矩阵)

### 7.2 Stage 0.5: 样本配对验证 (可选)

当提供 `--metadata_csv` 和 `--mapping_csv` 时，执行强约束配对:

```
metadata.csv: Ribo-Seq GSM → matched_RNA-seq GSM
srx_to_srr_mapping.csv: GSM → SRX → SRR
最终映射: { bio_sample: { ribo_srr: SRR_x, rna_srr: SRR_y } }
```

- 无法追踪到 SRR 的配对 → DROP + 日志
- 输出终端配对报告: `[配对成功] Sample X: Ribo=SRR123 ↔ RNA=SRR124`

### 7.3 Stage 1: 预处理 (CPM + 过滤 + Dummy Gene)

**函数**: `stage1_preprocess()`

100% 复刻 CenikLab `ribobase_counts_processing.py`:

1. **CPM 标准化**:
   $$CPM_{i,j} = \frac{count_{i,j}}{\sum_k count_{k,j}} \times 10^6$$

2. **低表达基因过滤** (`dummy_gene_df()`):
   - 对每个基因: 统计在多少个样本中 CPM ≥ `cpm_cutoff` (默认 1)
   - 若该比例 < `overall_cutoff`% (默认 70%), 标记为 dummy gene

3. **Dummy Gene 合并**:
   - Ribo 和 RNA 的 dummy gene 取并集 (`outer merge`)
   - 所有 dummy gene 的原始计数求和，合并为一行 `dummy_gene`
   - 保留 dummy_gene 行确保成分数据分析 (CoDA) 的闭合性

4. **非 polyA 基因移除** (可选):
   - 若提供 `--nonpolya_csv`，移除列表中的基因
   - 这些基因在 polyA 富集的 RNA-seq 文库中可能系统性低估

**输出**: `ribo_paired_count_dummy.csv`, `rna_paired_count_dummy.csv`

### 7.4 Stage 2: CLR/ILR 成分回归 (R)

**脚本**: `src/te_calc/TE.R`

100% 复刻 CenikLab `src/TE.R`:

**数学公式**:

1. **CLR 变换** (Centered Log-Ratio):
   ```r
   pr_RIBO <- propr(RIBO, metric="rho", ivar="clr", alpha=NA, p=100)
   ```
   $$CLR(x_i) = \ln\frac{x_i}{g(\mathbf{x})}$$
   其中 $g(\mathbf{x})$ 是所有组分的几何均值。
   
   `propr` 内部对零值执行 **乘性简单替换** (multiplicative simple replacement)。

2. **CLR → ILR 变换** (Isometric Log-Ratio):
   ```r
   RIBO_ilr <- clr2ilr(pr_RIBO@logratio)
   RNA_ilr  <- clr2ilr(pr_RNA@logratio)
   ```
   ILR 将 D 维成分数据映射到 (D-1) 维欧氏空间，消除闭合约束。

3. **成分线性回归** (逐基因/逐 ILR 分量):
   ```r
   m <- lm(RIBO_ilr[,i] ~ RNA_ilr[,i])
   ```
   对每个 ILR 分量，以 RNA 为自变量、Ribo 为因变量进行线性回归。

4. **TE = 回归残差的 CLR 表示**:
   ```r
   TE <- ilr2clr(resid(m))
   ```
   残差代表 Ribo-seq 中 **不能被 RNA 丰度解释** 的部分，即翻译效率。

**输出**: `human_TE_sample_level.rda` (R 数据对象)

### 7.5 Stage 3: 后处理

**函数**: `stage3_postprocess()`

1. 将 `.rda` 转换为 CSV
2. (可选) 按 cell_line 分组取均值
3. 转置输出: `te_results_final.csv` (genes × cell_lines)

---

## 8. 关键发现：Ribo↔RNA 配对机制

### 8.1 配对实现方式

原作者通过 **三层机制** 确保 Ribo-seq 和 RNA-seq 的天然配对：

| 层级 | 位置 | 机制 |
|------|------|------|
| **L1: 元数据驱动** | `generate_yaml.py` | `matched_RNA-seq_experiment_alias` 列建立 GSM↔GSM 映射 |
| **L2: YAML 键对齐** | `project.yaml` | `input.fastq[GSM]` 和 `rnaseq.fastq[GSM]` 使用相同 key |
| **L3: Nextflow Channel Join** | `RiboFlow.groovy:1981` | `RIBO_FOR_RNASEQ.join(RNASEQ_BED_FOR_RIBO_FINAL)` 按 sample name 匹配 |

### 8.2 打包命令总结

| 步骤 | 命令 | 作用 |
|------|------|------|
| **创建 .ribo** | `ribopy create -n SAMPLE --alignmentfile RIBO.bed ...` | 从 Ribo-seq BED 创建 HDF5 容器 |
| **注入 RNA-seq** | `ribopy rnaseq set -n SAMPLE -a RNA.bed -f bed --force SAMPLE.ribo` | 将配对的 RNA-seq BED 写入同一个 .ribo 文件 |
| **合并** | `ribopy merge all.ribo *.ribo` | 将所有 per-sample .ribo 合并为一个文件 |

### 8.3 为什么 .ribo 能规避配对错误

`.ribo` 文件的 HDF5 结构将同一个 experiment name 下的 Ribo-seq 和 RNA-seq 数据 **物理绑定** 在一起。当下游使用 `ribopy` API 提取数据时：

```python
ribo = Ribo("sample.ribo")
ribo_counts = ribo.get_region_counts("CDS")     # Ribo-seq
rna_counts  = ribo.get_rnaseq()                  # RNA-seq (同一 experiment)
```

**两种数据共享相同的 experiment name 和转录本索引**，不可能发生错配。这是比任何元数据表匹配都更可靠的配对机制。

---

## 9. 完整数据流图

```
NCBI SRA
   │
   ▼
┌──────────────────────────────────┐
│  download_sra.py                 │
│  Entrez API: SRX → SRR           │
│  sradownloader: SRR → FASTQ      │
└──────────────┬───────────────────┘
               │  data/raw/fastq/*.fastq.gz
               ▼
┌──────────────────────────────────┐
│  Snakescale (Snakefile)          │
│  check_adapter → check_lengths   │
│  generate_yaml.py:               │
│    CSV → per-study project.yaml  │
│    ★ Ribo↔RNA GSM 配对写入 YAML  │
└──────────────┬───────────────────┘
               │  data/interim/project/{GSE}/{study}.yaml
               ▼
┌──────────────────────────────────────────────────────────────┐
│  RiboFlow (Nextflow: RiboFlow.groovy)                        │
│                                                              │
│  ┌─── Ribo-seq Track ────┐   ┌─── RNA-seq Track ────┐       │
│  │ cutadapt              │   │ cutadapt             │       │
│  │   ↓                   │   │   ↓                  │       │
│  │ bowtie2 filter (rRNA) │   │ bowtie2 filter       │       │
│  │   ↓                   │   │   ↓                  │       │
│  │ bowtie2 transcriptome │   │ bowtie2 transcriptome│       │
│  │   ↓                   │   │   ↓                  │       │
│  │ samtools -q (MAPQ)    │   │ samtools -q (MAPQ)   │       │
│  │   ↓                   │   │   ↓                  │       │
│  │ bamToBed              │   │ bamToBed             │       │
│  │   ↓                   │   │   ↓                  │       │
│  │ merge + dedup         │   │ merge + dedup        │       │
│  │   ↓                   │   │   ↓                  │       │
│  │ RIBO BED              │   │ RNA BED              │       │
│  └───────┬───────────────┘   └───────┬──────────────┘       │
│          │                           │                       │
│          ▼                           │                       │
│  ┌─────────────────────┐             │                       │
│  │ ribopy create       │             │                       │
│  │ → {sample}.ribo     │◄────────────┘                       │
│  │                     │  ribopy rnaseq set                  │
│  │ (HDF5 容器)         │  ★ RNA-seq 注入同一文件              │
│  └─────────┬───────────┘                                     │
│            │                                                 │
│            ▼                                                 │
│  ┌─────────────────────┐                                     │
│  │ ribopy merge        │                                     │
│  │ → all.ribo          │                                     │
│  └─────────┬───────────┘                                     │
└────────────┼─────────────────────────────────────────────────┘
             │  data/processed/ribo_files/*.ribo
             ▼
┌──────────────────────────────────────────────────────────────┐
│  te_calculator.py                                            │
│                                                              │
│  Stage 0:   ribopy API → ribo_raw.csv + rnaseq_raw.csv      │
│  Stage 0.5: (可选) 元数据强约束配对验证                       │
│  Stage 1:   CPM → 低表达过滤 → dummy gene → 非polyA移除      │
│  Stage 2:   TE.R: CLR/ILR 成分回归 → TE = 残差               │
│  Stage 3:   RDA → CSV → 分组均值 → 转置                      │
│                                                              │
│  Output: data/processed/te_results_final.csv                 │
└──────────────────────────────────────────────────────────────┘
```

---

## 10. 工具与版本清单

| 工具 | 用途 | 所在阶段 |
|------|------|----------|
| **sradownloader** | SRA FASTQ 批量下载 | A |
| **cutadapt** | Adapter 修剪 + 长度/质量过滤 | C1, D |
| **bowtie2** | 短读段比对 (rRNA 过滤 + 转录组) | C2, C3, D |
| **samtools** | BAM 排序/索引/MAPQ过滤/统计 | C2–C4, D |
| **bamToBed** (bedtools) | BAM → BED 格式转换 | C5, D |
| **rfc** (riboflow-commands) | 去重, 日志合并, 统计编译 | C6, 统计 |
| **ribopy create** | 从 BED 创建 .ribo (HDF5) | E1 |
| **ribopy rnaseq set** | 注入 RNA-seq 数据到 .ribo | E2 |
| **ribopy merge** | 合并多个 .ribo 文件 | E3 |
| **ribopy** (Python API) | 从 .ribo 提取计数矩阵 | F (Stage 0) |
| **R: propr** | CLR 变换 + rho 比例性分析 | F (Stage 2) |
| **R: compositions** | clr2ilr / ilr2clr 变换 | F (Stage 2) |
| **R: tidyverse** | 数据整理 (dplyr, tidyr) | F (Stage 2) |
| **R: foreach + doParallel** | 并行逐基因回归 | F (Stage 2) |
| **Nextflow** | RiboFlow 工作流引擎 | C–E |
| **Snakemake** | Snakescale DAG 编排 | B |
| **Docker** | 容器化执行环境 (`ceniklab/riboflow`) | C–E |

---

> **文档结束**  
> 本白皮书基于对 `workflow/snakescale/riboflow/RiboFlow.groovy` (2053行)、`workflow/snakescale/scripts/generate_yaml.py` (407行)、`workflow/snakescale/Snakefile` (722行)、`src/te_calc/te_calculator.py` (1035行) 和 `src/te_calc/TE.R` (121行) 的逐行源码审计生成。
