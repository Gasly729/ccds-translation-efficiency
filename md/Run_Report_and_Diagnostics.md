# TE Pipeline 运行状态总结报告

> 生成时间：2026-04-02 18:30  
> 最后一次运行：`snakemake_run8.log` (2026-04-02 17:49 – 18:16)  
> 服务器配置：208 核 / ~1.2 TB RAM

---

## 一、总体进度概览

| 维度 | 数值 |
|------|------|
| 原始 FASTQ 文件数 | 4504 |
| 配置的研究 (studies) | 4 个 |
| Snakemake 总步骤 | 101 |
| 已完成步骤 | 97（含 1 个 `run_riboflow` 成功写入"失败"状态） |
| 失败步骤 | 4（均为 `run_riboflow` 规则） |
| Snakemake 完成率 | **96%**（前处理阶段 100%，Nextflow 阶段受阻） |

---

## 二、四个研究的逐一状态

### 2.1 GSE50597 — 拟南芥 (Arabidopsis)

| 项目 | 详情 |
|------|------|
| 物种 | Arabidopsis thaliana |
| Ribo-seq 样本 | 6 个 (GSM1224473–GSM1224478) |
| RNA-seq 样本 | 6 个 (同上，配对) |
| 前处理 (check_adapter) | ✅ 10/10 TSV 完成 |
| adapter 警告 | 4 个 SRR 低 adapter 存在，但通过验证 |
| Nextflow 状态 | ❌ 在 `individual_alignment_stats` 阶段因 `rfc` 命令缺失退出 (exit 127) |
| 已落盘中间文件 | 110 个文件，41 GB |

**Ribo-seq 阶段推进：**

| 阶段 | 文件数 | 状态 |
|------|--------|------|
| clip | 10 | ✅ 完成 |
| filter | 30 | ✅ 完成 |
| transcriptome_alignment | 18 | ✅ 完成 |
| quality_filter | 4 | ⚠️ 部分完成 (1/6 样本) |

**RNA-seq 阶段推进：**

| 阶段 | 文件数 | 状态 |
|------|--------|------|
| rnaseq/clip | 12 | ✅ 完成 |
| rnaseq/filter | 18 | ✅ 完成 |
| rnaseq/transcriptome_alignment | 12 | ✅ 完成 |
| rnaseq/quality_filter | 4 | ⚠️ 部分完成 (1/6 样本) |
| rnaseq/bam_to_bed | 2 | ⚠️ 部分完成 (1/6 样本) |

---

### 2.2 GSE115647 — 人类 (Human)

| 项目 | 详情 |
|------|------|
| 物种 | Homo sapiens |
| Ribo-seq 样本 | 2 个 (GSM3186648, GSM3186650) |
| RNA-seq 样本 | 2 个 (同上，配对) |
| 前处理 (check_adapter) | ✅ 3/3 TSV 完成 |
| adapter 警告 | 无 |
| Nextflow 状态 | ❌ 在 `individual_alignment_stats` 阶段因 `rfc` 命令缺失退出 (exit 127) |
| 已落盘中间文件 | 29 个文件，1.2 GB |

**Ribo-seq 阶段推进（仅 GSM3186650 完整推进）：**

| 阶段 | 文件数 | 状态 |
|------|--------|------|
| clip | 4 | ✅ 2/2 样本完成 |
| filter | 12 | ✅ 2/2 样本完成 |
| transcriptome_alignment | 6 | ✅ 1/2 样本完成 |
| quality_filter | 4 | ✅ 1/2 样本完成 |
| bam_to_bed | 3 | ✅ 1/2 样本完成 |

> RNA-seq 阶段：尚未启动（Nextflow 在 Ribo-seq 侧的 `individual_alignment_stats` 就已崩溃）

---

### 2.3 GSE48140 — 线虫 (Caenorhabditis)

| 项目 | 详情 |
|------|------|
| 物种 | C. elegans / C. briggsae / C. brenneri / C. remanei（多物种混合） |
| Ribo-seq 样本 | 23 个 |
| RNA-seq 样本 | 23 个 (配对) |
| 前处理 (check_adapter) | ✅ 46/46 TSV 完成 |
| adapter 警告 | 无 |
| 已排除样本 | SRR914374（CRC 校验损坏，已从 mapping CSV 移除） |
| Nextflow 状态 | ❌ 仅推进到 clip + filter 阶段即因级联失败中止 |
| 已落盘中间文件 | 126 个文件，20 GB |

**Ribo-seq 阶段推进：**

| 阶段 | 文件数 | 状态 |
|------|--------|------|
| clip | 36 | ✅ 完成 |
| filter | 90 | ✅ 完成 |
| transcriptome_alignment | — | ❌ 未启动 |

> 注：此 study 的 Nextflow 运行属于第一个被调度的（最先开始），但因 `rfc` 缺失导致后续 study 级联中断，Nextflow 在 `individual_alignment_stats` 处终止后，剩余样本被 ABORT。

---

### 2.4 GSE102659 — 小鼠 (Mouse)

| 项目 | 详情 |
|------|------|
| 物种 | Mus musculus |
| Ribo-seq 样本 | 10 个 (GSM2742552–GSM2742572) |
| RNA-seq 样本 | 10 个 (配对) |
| 前处理 (check_adapter) | ✅ 20/20 TSV 完成 |
| adapter 问题 | ⛔ 10 个 SRR 低 adapter 存在；无法猜测 adapter（检测到多种冲突 adapter） |
| 分类结果 | **Failed pre-run checks** — 未进入 Nextflow |
| Nextflow 状态 | ⛔ 未运行 |
| 已落盘中间文件 | 0（仅有 adapter_check TSV） |

**失败原因日志：**
```
Unable to guess adapter. Study invalid.
Guessed adapters: ['AAAAAAAAAAAAAAA', 'CTGTAGGCACCATCAAT', None]
Multiple adapters detected. Study invalid.
```

---

## 三、已落盘文件汇总

### 3.1 `data/interim/` 目录（101 个文件）

| 子目录 | 文件数 | 说明 |
|--------|--------|------|
| `adapter_check/` | 79 个 TSV | 4 个 study 的 cutadapt adapter 检测报告 |
| `project/` | 8 个 YAML | 4 个 study × 2（原始 + 模板） |
| `modified_project/` | 4 个 YAML | 修改后的 Nextflow 参数文件 |
| `modifications/` | 4 个 YAML | adapter/length 修正记录 |
| `log/` | 9 个文件 | status.txt, valid_studies.txt, yaml_status.txt 等 |

### 3.2 `intermediates/` 目录（storeDir 缓存，61 GB）

| Study | 大小 | 文件数 | 最远阶段 |
|-------|------|--------|----------|
| GSE50597 | 41 GB | 110 | bam_to_bed (部分) |
| GSE48140 | 20 GB | 126 | filter (全部完成) |
| GSE115647 | 1.2 GB | 29 | bam_to_bed (部分) |
| GSE102659 | — | 0 | 未进入 Nextflow |

### 3.3 `data/processed/` 目录

| 文件 | 说明 |
|------|------|
| `riboflow_status/GSE102659/riboflow_status.txt` | "Study GSE102659 failed." |
| `bombyx_te_results.csv` | 历史遗留文件 |
| `te_results_final.csv` | 历史遗留文件 |

> **未生成任何 `.ribo` 文件。**

### 3.4 Nextflow `work/` 目录（16 GB）

- 143 个哈希子目录
- 包含临时 BAM、FASTQ、日志等过程文件

### 3.5 Nextflow trace (最后一次运行 — GSE50597)

| 状态 | 进程数 |
|------|--------|
| COMPLETED | 25 |
| ABORTED | 22 |
| FAILED | 4 |

**FAILED 进程均为：**
- `individual_alignment_stats` (exit 127) × 2
- `rnaseq_individual_alignment_stats` (exit 127) × 2

---

## 四、已完成的 Snakemake 规则统计

| 规则 | 完成数 | 说明 |
|------|--------|------|
| `check_adapter` | 79 | 全部 4 study 的所有样本 |
| `check_adapter_stats` | 4 | 每个 study 一个汇总 |
| `check_lengths` | 4 | 每个 study 一个 |
| `guess_adapter` | 4 | 每个 study 一个 |
| `modify_yaml` | 4 | 每个 study 一个 |
| `classify_studies` | 1 | 生成 valid_studies.txt |
| `run_riboflow` | 1（成功写入失败状态） | GSE102659 写入 "failed" |
| **合计** | **97 / 101** | |

---

## 五、阻塞原因分析

### 5.1 唯一阻塞点：`rfc` 命令缺失 (exit 127)

所有 3 个通过预检的 study（GSE50597、GSE115647、GSE48140）在 Nextflow 的 `individual_alignment_stats` 进程处统一失败，原因是 **`rfc`（riboflow-commands）工具未安装**。

`rfc` 是 RiboFlow 管线的配套 CLI 工具，用于：
- `rfc compile-step-stats` — 汇编各步骤对齐统计
- `rfc bt2-log-to-csv` — 将 bowtie2 日志转换为 CSV

该命令在 genome_alignment 及 post_genome_alignment 阶段也被调用，即使安装后仍需确认这些阶段能否正常运行。

### 5.2 GSE102659 被分类为无效 study

原因：adapter 自动检测失败（检测到冲突的多种 adapter 序列），pipeline 将其标记为 invalid 并跳过 Nextflow 执行。需人工指定 adapter 或调整检测参数。

---

## 六、已实施的修复（本次会话）

| # | 修复内容 | 状态 |
|---|----------|------|
| 1 | 从 `srx_to_srr_mapping.csv` 移除损坏的 SRR914374 | ✅ |
| 2 | 降级 Nextflow 至 22.10.6（DSL1 兼容） | ✅ |
| 3 | 安装 fastqc、hisat2 | ✅ |
| 4 | 修复 `RiboFlow.groovy` 中 `{task.cpus}` → `${task.cpus}` | ✅ |
| 5 | 移除 `samtools index/idxstats` 的 `-@` 参数（samtools 1.6 不支持） | ✅ |
| 6 | 添加 `nextflow.config` 中 trace/timeline/report 的 `overwrite = true` | ✅ |
| 7 | 重装 cutadapt/bowtie2（conda 依赖冲突修复） | ✅ |
| 8 | 安装 biopython（Snakefile 依赖） | ✅ |

---

## 七、下一步行动建议

| 优先级 | 行动 | 预计影响 |
|--------|------|----------|
| **P0** | 安装 `rfc`（riboflow-commands CLI） | 解除 3 个 study 的统一阻塞 |
| **P1** | 重跑 `run_riboflow`（可 `-resume` 复用已有 storeDir 缓存） | 直接跳过已完成的 clip/filter/alignment |
| **P2** | 人工处理 GSE102659 的 adapter 配置 | 使第 4 个 study 可运行 |
| **P3** | 验证 `.ribo` 文件生成 → 进入下游 TE 计算 | 最终目标 |

---

## 八、关键文件路径索引

```
data/interim/
├── adapter_check/          # 79 个 cutadapt TSV 报告
│   ├── GSE102659/          # 20 个
│   ├── GSE115647/          # 3 个
│   ├── GSE48140/           # 46 个
│   └── GSE50597/           # 10 个
├── log/
│   ├── status.txt          # study 分类结果
│   ├── valid_studies.txt   # 通过预检的 study 列表
│   └── yaml_status.txt
├── modified_project/       # 4 个修改后的 Nextflow YAML
├── modifications/          # 4 个修正记录 YAML
└── project/                # 4 个原始 project YAML

workflow/snakescale/
├── intermediates/          # 61 GB Nextflow storeDir 缓存
│   ├── GSE115647/          # 1.2 GB (到 bam_to_bed)
│   ├── GSE48140/           # 20 GB  (到 filter)
│   └── GSE50597/           # 41 GB  (到 bam_to_bed)
├── work/                   # 16 GB Nextflow 临时工作目录
├── nextflow_logs/          # trace.txt, timeline.html, report.html
└── .nextflow.log           # 最后一次 Nextflow 运行日志

logs/
├── snakemake_run2–8.log    # 历次 Snakemake 运行日志
```
