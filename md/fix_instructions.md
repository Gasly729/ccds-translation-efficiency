请阅读当前目录下的 `pipeline_failure_report.md`、`PIPELINE_DETAILS.md` 和 `README.md`。前序测试（使用 snakemake-ribo 环境）已暴露出流水线中的 4 个阻断性 Bug。请你作为资深生信研发工程师，直接修改项目中对应的代码和配置文件来修复这些问题。

请按以下优先级逐一修复并验证代码逻辑：

### 1. 修复 Nextflow Trace 文件竞争条件 (P0)
**问题**：Snakemake 并发执行多个 `run_riboflow` 规则时，所有的 Nextflow 实例都试图写入同一个硬编码的 `./nextflow_logs/trace.txt`（由 `riboflow/nextflow.config` 定义），导致 `overwrite = true` 发生竞态崩溃。
**修复方案**：
* 修改 `workflow/snakescale/Snakefile` 中的 `run_riboflow` 规则（约第 712 行），在启动命令中增加 CLI 参数覆盖，例如使用 `-with-trace nextflow_logs/{wildcards.study}_trace.txt`。
* 确保同时处理 `-with-report` 和 `-with-timeline`（如果有的话），使其输出路径均包含 `{wildcards.study}` 以实现物理隔离。

### 2. 修复 TE 计算列名映射断裂 (P0)
**问题**：`src/te_calc/te_calculator.py`（约 444 行 `validate_and_align_columns()`）强制期望 `.ribo` 提取出的 Count 矩阵列名为 SRR ID。但根据上游 `generate_yaml.py` 的设计，`.ribo` 容器内默认使用的是实验的 GSM ID（天然配对锚点），或者当使用 `all.ribo` 时列名可能是 `"all"`。
**修复方案**：
* 审查并修改 `te_calculator.py` 中的列名对齐逻辑。
* 让验证函数能够正确识别和处理基于 GSM ID 的列名，并通过 `metadata.csv` 或 `mapping_csv` 成功将其与配对字典（强约束模式）关联，不要让其因为找不到 SRR ID 而抛出 RuntimeError。

### 3. 修复 H. sapiens Adapter 检测失败 (P1)
**问题**：人类样本 SRR953773 的 Adapter 含量极低，触发了 `guess_adapter` 失败，导致整个 Study 被 Snakemake `classify_studies` 规则判定为 Invalid 并跳过。
**修复方案**：
* 在不破坏全局动态预测逻辑的前提下，修改 `workflow/snakescale/config.yaml`，设置（或提供选项）允许 `override: True` 以跳过预检强行执行，**或者**提供一个脚本动态向 `data/external/test_metadata.csv` 中为该样本注入标准的 Illumina 3' Adapter 序列（如 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`）。

### 4. 修复 C. elegans 参考基因组不完整导致的目录冲突 (P1)
**问题**：`workflow/snakescale/scripts/references.yaml` 中 C. elegans 的 `regions` 和 `transcript_lengths` 为空字符串。YAML 生成器将其解析为裸目录 `reference/`，导致 Nextflow 在 `create_ribo` 阶段 stage 文件时发生同名文件（目录）碰撞。同时缺少线虫的转录组 Bowtie2 索引。
**修复方案**：
* 修改 `workflow/snakescale/scripts/generate_yaml.py`，增加异常捕获：如果 `references.yaml` 中的路径解析为空或不存在，抛出明确的提示错误，而不是将根目录 `reference/` 写入 project.yaml。
* 编写一个独立的辅助脚本（如 `scripts/build_celegans_ref.sh`），包含使用 `gffread` 等工具从 Ensembl/WormBase 下载并构建线虫 Transcriptome Bowtie2 索引、生成 `regions.bed` 和 `transcript_lengths.tsv` 的标准命令，以备后续运行。

完成上述修改后，请给出你修改的差异总结（git diff 摘要）。注意：**绝对不能**破坏 `PIPELINE_DETAILS.md` 中描述的 CoDA (CLR/ILR) 数学逻辑以及通过 YAML 字典键值实现 Ribo-RNA 天然配对的机制。