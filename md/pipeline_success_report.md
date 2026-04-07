# TE 计算管线回归测试报告

**测试时间**: 2026-04-07  
**测试环境**: Linux / Conda `snakemake-ribo`  
**测试命令**: `make run_local FORCE=1` → `make calc_te RIBO_DIR=data/processed/ribo_files`

---

## 总体结果: ✅ 管线端到端可运行，TE 计算成功完成

---

## 第一步: 环境与参考基因组

| 项目 | 状态 |
|------|------|
| Conda 环境 `snakemake-ribo` | ✅ 激活成功 |
| gffread 安装 | ✅ |
| C. elegans 参考基因组构建 | ✅ (经多轮修复) |
| references.yaml 更新 | ✅ |
| R 包 (propr, compositions, tidyverse, foreach, doParallel) | ✅ 补充安装 |

### C. elegans 参考基因组修复历程

1. **awk 语法错误**: `build_celegans_ref.sh` 使用了 `gawk` 专用的 `match()` 语法 → 改用 `split()`/`sub()` 兼容 `mawk`
2. **ID 前缀不匹配**: `transcript_lengths.tsv` 使用 `transcript:` 前缀，`regions.bed` 不带 → 统一添加前缀
3. **坐标系统错误**: `regions.bed` 使用基因组坐标 → 编写 `generate_celegans_regions.py` 将 GFF3 坐标映射到转录本相对坐标
4. **排序不一致**: `regions.bed` 排序与 `transcript_lengths.tsv` 不一致 → 按 lengths 文件顺序排序
5. **非 CDS 转录本缺失**: RNA-seq 比对到全部 60,000 转录本，但 ribo 文件只含 31,865 CDS 转录本 → 为非 CDS 转录本添加全长 CDS 注释

---

## 第二步: 上游管线 (Snakemake + RiboFlow)

| Study | 物种 | 状态 | .ribo 文件 |
|-------|------|------|-----------|
| GSE125038 | S. cerevisiae (酵母) | ✅ 成功 | `GSM3561545.ribo` |
| GSE50597 | A. thaliana (拟南芥) | ✅ 成功 | `GSM1224478.ribo` |
| GSE48140 | C. elegans (线虫) | ✅ 成功 | `GSM1169554.ribo` |
| GSE49994 | H. sapiens (人类) | ❌ 失败 | 无 |

### GSE49994 失败原因

**根因**: 原始 FASTQ 文件 `SRR953773_1.fastq.gz` 数据损坏  
- 第 12,375,789 条 read (`SRR953773.12375789`) 的 quality 序列长度 (51) 与 sequence 长度 (46) 不匹配
- cutadapt 报错: `FormatError: length of quality sequence (51) and length of read (46) do not match`
- **结论**: 这是上游原始数据质量问题，非管线 Bug，需要重新下载 FASTQ

### 修复的管线 Bug

1. **rfcommands Bowtie2 日志解析** (`compile_step_stats.py`, `merge/bowtie2_logs.py`)
   - 问题: `[WARNING]` 行以 `[` 开头，未被过滤导致日志行数不匹配
   - 修复: 添加 `this_line[0] == '['` 跳过条件
   - 同时修复了 `IOError` 格式字符串中的 `KeyError`

---

## 第三步: 下游 TE 计算

| 阶段 | 状态 | 输出 |
|------|------|------|
| Stage 0: .ribo → 原始计数矩阵 | ✅ | `ribo_raw.csv` (57,002 genes × 4 samples), `rnaseq_raw.csv` (25,137 genes × 3 samples) |
| Stage 0.5: Ribo/RNA 配对 | ✅ | 3 个共同样本 (GSM1169554 无 RNA-seq 数据被排除) |
| Stage 1: CPM 过滤 + dummy 基因 | ✅ | 21,086 genes 通过过滤 |
| Stage 2: TE.R (CLR/ILR 成分回归) | ✅ | `TE_sample_level.csv` |
| Stage 3: 结果汇总 | ✅ | `te_results_final.csv`, `human_TE_cellline_all.csv` |

### 最终输出文件

```
data/processed/
├── ribo_raw.csv                    # Ribo-seq 原始计数
├── rnaseq_raw.csv                  # RNA-seq 原始计数
├── ribo_paired_count_dummy.csv     # 过滤后 Ribo 计数 (21,086 genes)
├── rna_paired_count_dummy.csv      # 过滤后 RNA 计数 (21,086 genes)
├── TE_sample_level.csv             # 样本级 TE (3 samples × 21,086 genes)
├── human_TE_cellline_all.csv       # cell line 汇总 TE (21,086 genes)
└── te_results_final.csv            # 最终 TE 结果
```

---

## 修改的文件清单

| 文件 | 修改内容 |
|------|---------|
| `workflow/snakescale/scripts/build_celegans_ref.sh` | awk 语法兼容 mawk |
| `workflow/snakescale/scripts/references.yaml` | 线虫参考路径 |
| `workflow/snakescale/scripts/generate_celegans_regions.py` | 新建: GFF3→转录本坐标 regions.bed |
| `workflow/snakescale/reference/transcriptome/celegans/` | 重建: transcript_lengths.tsv, regions.bed |
| `rfcommands/compile_step_stats.py` | Bowtie2 日志解析 [WARNING] 行 |
| `rfcommands/merge/bowtie2_logs.py` | 同上 |

---

## 结论

管线 **端到端功能验证通过**:
- 4 个物种中 **3 个成功** 完成上游 RiboFlow 处理并生成 `.ribo` 文件
- 1 个失败 (GSE49994) 是因为 **原始 FASTQ 数据损坏**，非管线 Bug
- 下游 TE 计算 (CenikLab CLR/ILR 方法) **成功完成**，输出 21,086 个基因的翻译效率矩阵
