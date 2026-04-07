# sradownloader 输入准备报告

> 生成时间：2026-03-13  
> 数据来源：`metadata/metadata.csv`（2644 条记录）  
> 目标工具：[sradownloader v3.11](https://github.com/s-andrews/sradownloader)

---

## 一、sradownloader 官方输入格式要求总结

通过阅读 sradownloader 源码（`read_samples` 函数，第 265–353 行），该工具支持 **三种输入方式**：

### 格式 1：SRA RunTable 文件（推荐）
- 来源：NCBI SRA Run Selector 导出的 `SraRunTable.txt`
- 格式：**CSV**（逗号分隔）
- **第一列表头必须为 `Run`**，值必须是 **SRR/ERR/DRR** 开头的 Run accession
- 可选列 `source_name`：用于为下载文件命名
- 示例：
  ```
  Run,Assay Type,source_name,...
  SRR5413015,ChIP-Seq,Spt3 ChEC-seq,...
  SRR5413016,ChIP-Seq,Spt3 ChEC-seq,...
  ```

### 格式 2：简单 accession 列表
- 纯文本文件，**每行一个 SRR/ERR/DRR accession**
- 无表头
- 示例：
  ```
  SRR5413015
  SRR5413016
  ```

### 格式 3：单个 accession（命令行参数）
- 直接在命令行传入一个 SRR/ERR accession
- 示例：`sradownloader SRR5413015`

### ⚠️ 关键约束
| 约束项 | 说明 |
|--------|------|
| **ID 类型** | **仅支持 SRR / ERR / DRR (Run accession)**，不支持 GSE、GSM、SRX、SRS、PRJNA 等其他编号 |
| **第一列** | 工具始终读取每行第一列作为 accession |
| **表头识别** | 仅通过第一列是否为 `"Run"` 来判断表头行 |
| **空行** | 自动跳过 |

---

## 二、当前数据表格状态评估

### 2.1 表格基本信息
| 项目 | 值 |
|------|-----|
| 文件 | `metadata/metadata.csv` |
| 总行数 | 2646（含 1 行标题注释 + 1 行表头 + 2644 行数据） |
| 总列数 | 27 |
| 分隔符 | 逗号（CSV） |

### 2.2 关键列分析

| 列索引 | 列名 | 内容 | 与下载的关系 |
|--------|------|------|-------------|
| col0 | `experiment_alias` | GSM 编号（2391）/ SRX 编号（199）/ ERX 编号（45）/ "SRR_ERS*" 格式（9） | 非 Run accession，不可直接使用 |
| col9 | `experiment_accession` | SRX/ERX 编号（2635 个有效值）| **Experiment accession，需转换为 Run accession** |
| col12 | `sample_accession` | SRS 编号 | Sample accession，非 Run |
| col10 | `study_name` | GSE 编号（如 GSE128216）| Series 编号，非 Run |

### 2.3 ⛔ 核心问题：缺少 SRR/ERR/DRR Run Accession

**全表 2644 条记录中，未发现任何有效的 SRR/ERR/DRR Run accession。**

现有 ID 层级关系：
```
GSE (Series) → GSM (Sample) → SRX/ERX (Experiment) → SRR/ERR/DRR (Run) ← sradownloader 需要这个！
```

你的表格停留在 **SRX/ERX (Experiment)** 层级，还需要向下映射一层才能获得 Run accession。

### 2.4 异常数据

| 异常类型 | 数量 | 详情 |
|----------|------|------|
| col9 (experiment_accession) 为空 | 9 行 | 行 1557–1565，alias 为 `SRR_ERS*` 格式（注意：这**不是**有效 SRR accession，只是命名格式） |
| `library_construction_protocol` 含逗号 | 大量 | CSV 引号已正确包裹，Python csv 模块可正常解析；但 awk/cut 直接切分会出错 |
| 特殊字符（如 `∆`, `µ`） | 少量 | 存在于描述性字段中，不影响 ID 提取 |

---

## 三、解决方案：SRX/ERX → SRR/ERR 映射 + 生成 sradownloader 输入文件

### 方案概述

1. 从 `metadata.csv` 提取所有唯一的 SRX/ERX accession（col9）
2. 通过 NCBI Entrez E-utilities API 批量映射 SRX/ERX → SRR/ERR/DRR
3. 处理 9 条缺失 SRX 的记录（通过 ERS accession 查询）
4. 生成 sradownloader 可直接读取的输入文件

### 执行脚本

详见同目录下的 **`prepare_sradownloader_input.py`** 脚本。

运行方式：
```bash
cd /home/xrx/my_projects/joker/26.3/download
python3 prepare_sradownloader_input.py
```

脚本将生成两个文件：
- **`sradownloader_input.txt`** — 简单 SRR accession 列表（每行一个），可直接被 sradownloader 读取
- **`srx_to_srr_mapping.csv`** — SRX → SRR 完整映射表，便于回溯

然后即可使用 sradownloader 下载：
```bash
./sradownloader/sradownloader --outdir fastq_output sradownloader_input.txt
```

---

## 四、备选方案

如果 Entrez API 速度太慢或网络不稳定，可以考虑：

### 方案 B：使用 pysradb（本地数据库查询，更快）
```bash
pip install pysradb
pysradb srx-to-srr --saveto srx_to_srr.tsv SRX5512335 SRX5512334 ...
```

### 方案 C：手动从 SRA Run Selector 批量导出
1. 访问 https://www.ncbi.nlm.nih.gov/Traces/study/
2. 搜索对应的 BioProject / GSE 编号
3. 导出 RunInfo Table（自带 `Run` 列）
4. 直接将导出文件传给 sradownloader
