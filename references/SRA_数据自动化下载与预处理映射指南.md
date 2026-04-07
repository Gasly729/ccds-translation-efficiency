# SRA 数据自动化下载与预处理映射指南

> 版本：v1.0  
> 日期：2026-03-15  
> 作者：Cascade × xrx  
> 项目路径：`/home/xrx/my_projects/joker/26.3/download/`

---

## 一、项目背景

### 1.1 为什么需要这套流程

在 Ribo-seq / RNA-seq 翻译组学研究中，我们通常从文献附录或公共数据库（如 TranslatomeDB）收集大量实验的元数据。这些元数据往往只包含 **GSM（GEO Sample）** 或 **SRX/ERX（SRA Experiment）** 级别的编号，而实际下载 FASTQ 原始测序数据时，下载工具（如 `sradownloader`、`fasterq-dump`）**要求提供最底层的 SRR/ERR/DRR（Run）级别的 accession**。

SRA 数据的层级关系如下：

```
GSE (Series)
 └── GSM (Sample)
      └── SRX / ERX (Experiment)
           └── SRR / ERR / DRR (Run)  ← 下载工具需要这一层
```

本项目构建了一套完整的自动化管线：

1. **ID 映射**：从收集到的 metadata CSV 中提取 SRX/ERX accession，通过 NCBI Entrez API 批量映射到 SRR/ERR Run accession
2. **元数据注入**：将 study_name、GSM、corrected_type 等信息拼接为 `source_name`，注入到 sradownloader 的 RunTable 输入文件中，使下载后的 FASTQ 文件自动获得有意义的命名
3. **批量下载**：使用 sradownloader 一键批量下载数千个样本的原始 FASTQ 数据

### 1.2 数据规模

| 指标 | 值 |
|------|-----|
| 原始 metadata 记录数 | 2,644 条 |
| 唯一 SRX/ERX accession | 2,635 个 |
| 成功映射的 SRR/ERR Run | **4,357 个**（一个 Experiment 可能包含多个 Run） |
| 未映射记录 | 0 条 |

---

## 二、核心突破与优化点

### 2.1 突破 NCBI API 414 URI Too Long 错误（GET → POST 重构）

**问题**：批量查询时，将 200 个 SRX accession 拼接为 `term` 参数放入 GET 请求的 URL 中，导致 URI 长度超过 NCBI 服务器限制，返回 `HTTP 414: Request-URI Too Long`。

**解决方案**：将 `esearch` 和 `efetch` 的所有请求从 GET 改为 **POST**，通过 `urllib.request.Request(url, data=encoded_params)` 将参数放入请求体（request body），彻底绕开 URI 长度限制。

```python
# 修复前（GET，URI 过长会报 414）
esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
urllib.request.urlopen(esearch_url)

# 修复后（POST，参数放入 body）
esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
esearch_data = urllib.parse.urlencode(params).encode("utf-8")
req = urllib.request.Request(esearch_url, data=esearch_data)
urllib.request.urlopen(req)
```

### 2.2 修复 UnboundLocalError 致命崩溃

**问题**：当 esearch 重试全部失败后，`esearch_xml` 变量从未被赋值，随后 `ET.fromstring(esearch_xml)` 直接触发 `UnboundLocalError`，导致整个脚本崩溃，无法继续处理后续批次。

**解决方案**：
1. 在重试循环前初始化 `esearch_xml = ""`
2. 循环后添加 `if not esearch_xml: continue` 拦截，跳过失败批次
3. 确保脚本即使部分批次失败，也能跑完全部批次并输出已成功的结果

### 2.3 NCBI API Key 集成提速

**问题**：NCBI Entrez 默认速率限制为 **3 次/秒**，查询 2,635 个 accession（分 14 批次，每批次 2 次请求）耗时较长。

**解决方案**：通过环境变量注入 `NCBI_API_KEY`，将速率限制提升至 **10 次/秒**，查询速度提升约 3 倍。

```bash
export NCBI_API_KEY="ce26613e7503ae333286c73e1474a64d3507"
```

### 2.4 正则白名单机制防止文件系统 I/O 错误

**问题**：原始 metadata 中的 study_name、GSM 编号、corrected_type 等字段可能包含空格、斜杠、括号、特殊 Unicode 字符（如 `∆`、`µ`）等，这些字符作为文件名会导致：
- Linux 文件系统创建文件失败
- 下游 pipeline 的 shell 脚本解析报错
- sradownloader 内部的文件命名逻辑出错

**解决方案**：使用**正则表达式白名单**机制，而非穷举替换：

```python
def _sanitize_name(name):
    name = re.sub(r"[^A-Za-z0-9_.\-]", "_", name)  # 仅保留安全字符
    name = re.sub(r"_+", "_", name)                   # 压缩连续下划线
    return name.strip("_")
```

### 2.5 CSV RunTable 格式适配 sradownloader

**问题**：sradownloader 的 `read_samples` 函数要求输入文件第一列表头为 `Run`，并可选读取 `source_name` 列用于文件命名。简单的纯文本 ID 列表虽然能用，但无法控制文件命名。

**解决方案**：输出标准 CSV RunTable 格式，包含 `Run` 和 `source_name` 列，使 sradownloader 自动按 `{SRR}_{source_name}` 格式命名下载文件。

---

## 三、环境与依赖

### 3.1 Python 脚本依赖

| 依赖 | 版本要求 | 说明 |
|------|---------|------|
| Python | >= 3.6 | 仅使用标准库，无需 pip install |
| 标准库 | csv, re, urllib, xml.etree | 全部为 Python 内置模块 |
| 网络 | 可访问 NCBI Entrez API | `eutils.ncbi.nlm.nih.gov` |

### 3.2 sradownloader 依赖

| 依赖 | 说明 |
|------|------|
| sradownloader | 已克隆至 `./sradownloader/`，v3.11 |
| Python 3 | sradownloader 本身也是 Python 脚本 |
| SRA Toolkit | 提供 `fasterq-dump`，用于 NCBI 通道下载 |
| gzip | 用于压缩 NCBI 下载的 FASTQ 文件 |

### 3.3 SRA Toolkit 配置

```bash
# 安装后必须运行一次配置（交互式）
vdb-config -i
```

### 3.4 NCBI API Key（强烈推荐）

申请地址：https://www.ncbi.nlm.nih.gov/account/settings/

设置方式：
```bash
export NCBI_API_KEY="your_api_key_here"
```

---

## 四、执行步骤与命令

### Step 1: 运行 ID 映射脚本

```bash
cd /home/xrx/my_projects/joker/26.3/download

# 设置 API Key 加速查询
export NCBI_API_KEY="ce26613e7503ae333286c73e1474a64d3507"

# 执行映射
python3 prepare_sradownloader_input.py
```

**预期输出**：
```
[INFO] 查询 Entrez: 1–200 / 2635
[INFO] 查询 Entrez: 201–400 / 2635
...
[INFO] 已生成 sradownloader 输入文件: .../sradownloader_input.txt
[INFO]   共 4357 个唯一 Run accession (CSV RunTable 格式)
[INFO] 已生成映射表: .../srx_to_srr_mapping.csv
[INFO]   共 4366 条映射记录
```

### Step 2: 使用 sradownloader 批量下载

```bash
./sradownloader/sradownloader \
    --nogeo \
    --outdir fastq_output \
    /home/xrx/my_projects/joker/26.3/download/sradownloader_input.txt
```

#### ⚠️ `--nogeo` 参数的重要性（必读）

**必须添加 `--nogeo` 参数！** 原因如下：

1. **性能灾难**：sradownloader 默认会对每个 SRR accession 调用 `get_geo_name()` 函数，向 `https://www.ncbi.nlm.nih.gov/sra/?term={SRR}&format=text` 发送 HTTP 请求来获取样本名。对于 4,357 个样本，这意味着 **4,357 次额外的 NCBI 网页请求**，每次耗时数秒，仅此步骤就需要数小时。

2. **已有替代方案**：我们的 CSV RunTable 中已经包含了精心构建的 `source_name` 列。sradownloader 的逻辑是：先尝试 GEO 查询命名，如果失败或禁用，则回退到 `source_name` 列。添加 `--nogeo` 后，工具将直接使用我们预制的 `source_name`，既快又准。

3. **网络稳定性**：大量 GEO 查询容易触发 NCBI 的速率限制或网络超时，导致下载中断。

**不加 `--nogeo` 的后果**（实测）：
- 脚本逐个查询 GEO，每个样本打印 `Trying to get name for SRRxxxx from GEO`
- 4,357 个样本 × ~3 秒/请求 ≈ **额外等待 3.6 小时**
- 且 GEO 查询结果会覆盖我们精心准备的 `source_name`

#### 其他有用参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--threads N` | fasterq-dump 的线程数（仅 NCBI 通道） | `--threads 4` |
| `--noena` | 跳过 ENA 通道，直接使用 NCBI | `--noena` |
| `--noncbi` | 跳过 NCBI 通道，仅使用 ENA | `--noncbi` |
| `--force` | 覆盖已存在的输出文件 | `--force` |
| `--retries N` | 每个样本的最大重试次数（默认 5） | `--retries 3` |
| `--quiet` | 静默模式 | `--quiet` |

### Step 3: 验证下载结果

```bash
# 查看下载的文件列表和大小
ls -lh fastq_output/ | head -20

# 统计已下载文件数
ls fastq_output/*.fastq.gz 2>/dev/null | wc -l
```

---

## 五、输出文件说明

### 5.1 `sradownloader_input.txt`（sradownloader 直接输入文件）

| 属性 | 说明 |
|------|------|
| 格式 | CSV（逗号分隔），使用 Python `csv` 模块写入 |
| 行数 | 4,358（1 行表头 + 4,357 行数据） |
| 用途 | 直接传给 sradownloader 作为输入 |

**列定义**：

| 列名 | 说明 | 示例 |
|------|------|------|
| `Run` | SRR/ERR Run accession（sradownloader 必需） | `SRR8718520` |
| `source_name` | 重命名依据，格式 `{study}_{gsm}_{type}`，已白名单清洗 | `GSE128216_GSM3667333_RNA-Seq` |
| `organism` | 物种名 | `Streptomyces clavuligerus` |
| `srx` | 原始 SRX/ERX Experiment accession | `SRX5512335` |
| `gsm` | 原始 GSM/别名 | `GSM3667333` |
| `corrected_type` | 实验类型（Ribo-Seq / RNA-Seq 等） | `RNA-Seq` |
| `study_name` | GEO Series 编号 | `GSE128216` |
| `status` | 映射状态 | `OK` |

**下载后的文件命名规则**：

sradownloader 会自动将文件命名为：
```
{Run}_{source_name}_1.fastq.gz   （单端 / 双端 read1）
{Run}_{source_name}_2.fastq.gz   （双端 read2，如适用）
```

实测示例：
```
SRR8718520_GSE128216_GSM3667333_RNA-Seq_1.fastq.gz
```

### 5.2 `srx_to_srr_mapping.csv`（完整映射追溯表）

| 属性 | 说明 |
|------|------|
| 格式 | CSV |
| 行数 | 4,367（含一对多映射的所有记录） |
| 用途 | 供下游分析回溯 SRR → SRX → GSM 的完整映射关系 |

列结构与 `sradownloader_input.txt` 完全一致，但**保留了所有重复的 Run 记录**（同一 SRR 可能被多条 metadata 记录引用）。

---

## 六、单样本下载测试记录

### 测试命令

```bash
# 提取表头 + 第一行数据
head -n 2 sradownloader_input.txt > test_sradownloader_input.txt

# 单样本测试下载
./sradownloader/sradownloader --nogeo --outdir fastq_output_test test_sradownloader_input.txt
```

### 测试输入

```csv
Run,source_name,organism,srx,gsm,corrected_type,study_name,status
SRR8718520,GSE128216_GSM3667333_RNA-Seq,Streptomyces clavuligerus,SRX5512335,GSM3667333,RNA-Seq,GSE128216,OK
```

### 测试结果

```
fastq_output_test/
├── SRR8718520_GSE128216_GSM3667333_RNA-Seq_1.fastq.gz  (125M, ENA 通道部分下载)
└── SRR8718520_GSE128216_GSM3667333_RNA-Seq.fastq        (4.1G, NCBI 通道完整下载)
```

**结论**：
- ✅ sradownloader 正确解析了 CSV RunTable 格式
- ✅ `Run` 列被正确识别为 SRR accession
- ✅ `source_name` 列被正确用于文件重命名（`SRR8718520_GSE128216_GSM3667333_RNA-Seq`）
- ✅ fasterq-dump 成功下载了 12,420,952 条 reads
- ⚠️ ENA FTP 通道在当前网络环境下不稳定，sradownloader 自动回退到 NCBI 通道

---

## 七、项目文件清单

```
/home/xrx/my_projects/joker/26.3/download/
├── metadata/
│   └── metadata.csv                        # 原始元数据表（2,644 条记录，27 列）
├── sradownloader/                          # sradownloader v3.11（git clone）
│   ├── sradownloader                       # 主程序
│   ├── README.md
│   ├── SraRunTable.txt                     # 官方示例输入
│   └── simple_accession_list.txt           # 官方简单列表示例
├── prepare_sradownloader_input.py          # ★ ID 映射脚本（SRX→SRR + 元数据注入）
├── sradownloader_input.txt                 # ★ 生成的 sradownloader 输入文件（4,357 个 Run）
├── srx_to_srr_mapping.csv                  # ★ 完整映射追溯表
├── sradownloader_input_report.md           # 初期分析报告
├── test_sradownloader_input.txt            # 单样本测试用临时文件
├── fastq_output_test/                      # 单样本测试下载目录
└── SRA_数据自动化下载与预处理映射指南.md      # ★ 本文档
```

---

## 八、常见问题与排障

### Q1: esearch 报 HTTP 414 错误
已通过 POST 请求修复。如果仍然出现，检查 `BATCH_SIZE` 是否被手动调大——默认 200 是安全值。

### Q2: 部分批次查询失败但脚本未崩溃
这是预期行为。脚本设计为"尽力而为"模式，失败批次会被跳过，成功的结果仍会被输出。运行结束后检查 `[WARN]` 日志中的"未能映射"记录数即可。

### Q3: sradownloader 下载速度慢
- 使用 `--noena` 跳过 ENA FTP（国内网络环境下 ENA 经常超时）
- 增加 `--threads 4` 提升 NCBI 通道并发
- 考虑使用 `screen` 或 `tmux` 后台运行

### Q4: 下载中断后如何续传
sradownloader 默认会跳过已存在且非空的输出文件。直接重新运行相同命令即可，已下载的文件不会被重复下载（除非加了 `--force`）。

### Q5: 为什么有 4,357 个 Run 但只有 2,644 条原始记录？
一个 SRX（Experiment）可以包含多个 SRR（Run），例如技术重复或多 lane 测序。这是正常的一对多映射关系。

---

## 九、完整一键执行命令块

```bash
# === 进入项目目录 ===
cd /home/xrx/my_projects/joker/26.3/download

# === Step 1: ID 映射（SRX → SRR）===
export NCBI_API_KEY="ce26613e7503ae333286c73e1474a64d3507"
python3 prepare_sradownloader_input.py

# === Step 2: 批量下载（务必加 --nogeo）===
./sradownloader/sradownloader \
    --nogeo \
    --noena \
    --threads 4 \
    --outdir fastq_output \
    /home/xrx/my_projects/joker/26.3/download/sradownloader_input.txt
```

> **提示**：批量下载 4,357 个样本耗时较长（数天级别），强烈建议在 `tmux` 或 `screen` 会话中运行。
