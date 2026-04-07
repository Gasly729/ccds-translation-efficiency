#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_sra.py — SRA 数据下载与 SRX→SRR 映射模块 (CCDS 重构版)

从 metadata.csv 提取 SRX/ERX accession，批量映射到 SRR/ERR Run accession，
生成 sradownloader 可直接读取的输入文件，并调用 sradownloader 下载 FASTQ。

CCDS 路径约定:
  - 元数据输入:        data/external/metadata/metadata.csv    (只读存档)
  - 映射表输出:        data/external/srx_to_srr_mapping.csv   (只读存档)
  - sradownloader 输入: data/external/sradownloader_input.txt  (只读存档)
  - 下载报告/日志:     references/sradownloader_input_report.md
  - FASTQ 输出:        data/raw/fastq/
  - sradownloader 工具: src/data/sradownloader/sradownloader

用法:
  # 仅生成映射 (不下载)
  python -m src.data.download_sra --prepare

  # 生成映射 + 下载 FASTQ
  python -m src.data.download_sra --prepare --download

  # 仅下载 (已有映射文件)
  python -m src.data.download_sra --download

  # 指定自定义路径
  python -m src.data.download_sra --prepare \
      --metadata data/external/metadata/metadata.csv \
      --outdir data/raw/fastq

注意:
  - 需要联网访问 NCBI Entrez E-utilities API
  - 大批量查询建议设置 NCBI_API_KEY 环境变量以提升速率限制（无 key: 3次/秒, 有 key: 10次/秒）
  - 如有 API key，设置方式: export NCBI_API_KEY="your_key_here"
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


# ============================================================================
# CCDS 路径配置 — 所有路径相对于项目根目录 (PROJECT_ROOT)
# ============================================================================
def _find_project_root() -> Path:
    """向上查找包含 Makefile 的项目根目录。"""
    current = Path(__file__).resolve().parent
    for _ in range(10):
        if (current / "Makefile").exists() and (current / "data").exists():
            return current
        current = current.parent
    # 回退: 假设 src/data/download_sra.py 在项目根下两层
    return Path(__file__).resolve().parent.parent.parent


PROJECT_ROOT = _find_project_root()

# ── CCDS 标准路径 ───────────────────────────────────────────────────────────
METADATA_CSV        = PROJECT_ROOT / "data" / "external" / "metadata" / "metadata.csv"
SRR_LIST_FILE       = PROJECT_ROOT / "data" / "external" / "sradownloader_input.txt"
MAPPING_FILE        = PROJECT_ROOT / "data" / "external" / "srx_to_srr_mapping.csv"
DOWNLOAD_REPORT     = PROJECT_ROOT / "references" / "sradownloader_input_report.md"
FASTQ_OUTPUT_DIR    = PROJECT_ROOT / "data" / "raw" / "fastq"
SRADOWNLOADER_BIN   = PROJECT_ROOT / "src" / "data" / "sradownloader" / "sradownloader"
LOG_DIR             = PROJECT_ROOT / "logs" / "download"

# ── Entrez API 配置 ─────────────────────────────────────────────────────────
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")
BATCH_SIZE = 200  # 每次 Entrez 查询的 ID 数量
RETRY_MAX = 3
RETRY_DELAY = 5  # 秒


# ============================================================================
# 日志
# ============================================================================
def log(msg):
    print(f"[INFO] {msg}", flush=True)


def warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr, flush=True)


# ============================================================================
# Step 1: 从 metadata.csv 提取 accession
# ============================================================================
def extract_accessions(csv_path):
    """
    返回:
      experiments: list of dict, 每条记录包含 {gsm, srx, organism, corrected_type, study_name}
      missing_srx: list of dict, col9 为空的异常记录
    """
    experiments = []
    missing_srx = []

    with open(csv_path, encoding="utf-8") as f:
        reader = csv.reader(f)
        # 跳过第一行（标题注释行: "Curated Data,..."）
        next(reader)
        # 第二行是真正的表头
        headers = next(reader)

        for row_num, row in enumerate(reader, start=3):
            if len(row) == 0:
                continue

            alias = row[0].strip()  # experiment_alias (GSM/SRX/ERX/SRR_ERS*)
            srx = row[9].strip() if len(row) > 9 else ""
            organism = row[3].strip() if len(row) > 3 else ""
            corrected_type = row[5].strip() if len(row) > 5 else ""
            study_name = row[10].strip() if len(row) > 10 else ""

            record = {
                "row": row_num,
                "gsm": alias,
                "srx": srx,
                "organism": organism,
                "corrected_type": corrected_type,
                "study_name": study_name,
            }

            if srx and re.match(r"^(SRX|ERX|DRX)\d+$", srx):
                experiments.append(record)
            elif srx == "":
                missing_srx.append(record)
                # 尝试从 alias 提取 ERS accession 用于后续查询
                ers_match = re.search(r"(ERS\d+)", alias)
                if ers_match:
                    record["ers"] = ers_match.group(1)
            else:
                # srx 列有值但不是标准 SRX/ERX/DRX 格式
                warn(f"行{row_num}: experiment_accession 格式异常: '{srx}' (alias={alias})")
                experiments.append(record)

    return experiments, missing_srx


# ============================================================================
# Step 2: 通过 Entrez 批量查询 SRX → SRR
# ============================================================================
def entrez_srx_to_srr(srx_list):
    """
    使用 NCBI Entrez E-utilities 将 SRX/ERX accession 批量映射到 SRR/ERR Run accession。
    返回: dict {srx: [srr1, srr2, ...]}
    """
    mapping = defaultdict(list)
    total = len(srx_list)

    for batch_start in range(0, total, BATCH_SIZE):
        batch = srx_list[batch_start: batch_start + BATCH_SIZE]
        batch_end = min(batch_start + BATCH_SIZE, total)
        log(f"查询 Entrez: {batch_start + 1}–{batch_end} / {total}")

        # Step 2a: esearch — 用 SRX accession 搜索 SRA 数据库，获取 UID
        query = " OR ".join(f"{acc}[Accession]" for acc in batch)
        esearch_params = {
            "db": "sra",
            "term": query,
            "retmax": len(batch) * 5,  # 一个 SRX 可能对应多个 SRR
            "usehistory": "y",
        }
        if NCBI_API_KEY:
            esearch_params["api_key"] = NCBI_API_KEY

        esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        esearch_data = urllib.parse.urlencode(esearch_params).encode("utf-8")

        esearch_xml = ""
        for attempt in range(1, RETRY_MAX + 1):
            try:
                req = urllib.request.Request(esearch_url, data=esearch_data)
                with urllib.request.urlopen(req, timeout=60) as resp:
                    esearch_xml = resp.read().decode("utf-8")
                break
            except Exception as e:
                warn(f"esearch 尝试 {attempt}/{RETRY_MAX} 失败: {e}")
                if attempt < RETRY_MAX:
                    time.sleep(RETRY_DELAY * attempt)
                else:
                    warn(f"esearch 最终失败，跳过本批次")

        if not esearch_xml:
            _rate_limit_sleep()
            continue

        root = ET.fromstring(esearch_xml)
        webenv = root.findtext("WebEnv", "")
        query_key = root.findtext("QueryKey", "")
        count = int(root.findtext("Count", "0"))

        if count == 0:
            warn(f"本批次未找到任何结果 (batch {batch_start + 1}–{batch_end})")
            _rate_limit_sleep()
            continue

        # Step 2b: efetch — 获取详细的 XML，从中提取 SRX → SRR 映射
        efetch_params = {
            "db": "sra",
            "query_key": query_key,
            "WebEnv": webenv,
            "rettype": "xml",
            "retmax": count,
        }
        if NCBI_API_KEY:
            efetch_params["api_key"] = NCBI_API_KEY

        efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        efetch_data = urllib.parse.urlencode(efetch_params).encode("utf-8")

        efetch_xml = ""
        for attempt in range(1, RETRY_MAX + 1):
            try:
                req = urllib.request.Request(efetch_url, data=efetch_data)
                with urllib.request.urlopen(req, timeout=120) as resp:
                    efetch_xml = resp.read().decode("utf-8")
                break
            except Exception as e:
                warn(f"efetch 尝试 {attempt}/{RETRY_MAX} 失败: {e}")
                if attempt < RETRY_MAX:
                    time.sleep(RETRY_DELAY * attempt)
                else:
                    warn(f"efetch 最终失败，跳过本批次")

        if not efetch_xml:
            _rate_limit_sleep()
            continue

        # 解析 efetch XML
        try:
            exp_root = ET.fromstring(efetch_xml)
        except ET.ParseError as e:
            warn(f"XML 解析失败: {e}")
            _rate_limit_sleep()
            continue

        for exp_pkg in exp_root.iter("EXPERIMENT_PACKAGE"):
            exp_elem = exp_pkg.find("EXPERIMENT")
            if exp_elem is not None:
                srx_acc = exp_elem.get("accession", "")
            else:
                continue

            for run in exp_pkg.iter("RUN"):
                srr_acc = run.get("accession", "")
                if srr_acc:
                    mapping[srx_acc].append(srr_acc)

        _rate_limit_sleep()

    return mapping


def entrez_ers_to_srr(ers_list):
    """
    使用 NCBI Entrez 将 ERS accession 映射到 SRR Run accession。
    返回: dict {ers: [srr1, srr2, ...]}
    """
    mapping = defaultdict(list)
    if not ers_list:
        return mapping

    log(f"查询 {len(ers_list)} 个 ERS accession ...")

    query = " OR ".join(f"{acc}[Accession]" for acc in ers_list)
    esearch_params = {
        "db": "sra",
        "term": query,
        "retmax": len(ers_list) * 5,
        "usehistory": "y",
    }
    if NCBI_API_KEY:
        esearch_params["api_key"] = NCBI_API_KEY

    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_data = urllib.parse.urlencode(esearch_params).encode("utf-8")

    try:
        req = urllib.request.Request(esearch_url, data=esearch_data)
        with urllib.request.urlopen(req, timeout=60) as resp:
            esearch_xml = resp.read().decode("utf-8")
    except Exception as e:
        warn(f"ERS esearch 失败: {e}")
        return mapping

    root = ET.fromstring(esearch_xml)
    webenv = root.findtext("WebEnv", "")
    query_key = root.findtext("QueryKey", "")
    count = int(root.findtext("Count", "0"))

    if count == 0:
        warn("ERS 查询未返回结果")
        return mapping

    efetch_params = {
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "xml",
        "retmax": count,
    }
    if NCBI_API_KEY:
        efetch_params["api_key"] = NCBI_API_KEY

    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    efetch_data = urllib.parse.urlencode(efetch_params).encode("utf-8")

    try:
        req = urllib.request.Request(efetch_url, data=efetch_data)
        with urllib.request.urlopen(req, timeout=120) as resp:
            efetch_xml = resp.read().decode("utf-8")
    except Exception as e:
        warn(f"ERS efetch 失败: {e}")
        return mapping

    try:
        exp_root = ET.fromstring(efetch_xml)
    except ET.ParseError as e:
        warn(f"ERS XML 解析失败: {e}")
        return mapping

    for exp_pkg in exp_root.iter("EXPERIMENT_PACKAGE"):
        # 查找 sample accession (ERS)
        sample_elem = exp_pkg.find(".//SAMPLE")
        ers_acc = sample_elem.get("accession", "") if sample_elem is not None else ""

        for run in exp_pkg.iter("RUN"):
            srr_acc = run.get("accession", "")
            if srr_acc and ers_acc:
                mapping[ers_acc].append(srr_acc)

    return mapping


def _rate_limit_sleep():
    """遵守 NCBI 速率限制"""
    if NCBI_API_KEY:
        time.sleep(0.1)  # 有 API key: 10 次/秒
    else:
        time.sleep(0.35)  # 无 API key: 3 次/秒


# ============================================================================
# Step 3: 生成输出文件
# ============================================================================
def _sanitize_name(name):
    """正则白名单清洗：仅保留字母、数字、下划线、点、中划线，其余替换为下划线，并压缩连续下划线。"""
    name = re.sub(r"[^A-Za-z0-9_.\-]", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def write_outputs(experiments, missing_srx, srx_mapping, ers_mapping,
                  srr_list_file, mapping_file):
    """
    生成:
      1. sradownloader_input.txt — CSV 格式，表头 Run,source_name,...，适配 sradownloader RunTable 模式
      2. srx_to_srr_mapping.csv  — 完整映射表，便于回溯

    输出路径由参数传入 (CCDS 标准路径)。
    """
    all_rows = []

    # 处理正常记录 (有 SRX/ERX)
    unmapped_srx = []
    for rec in experiments:
        srx = rec["srx"]
        srr_list = srx_mapping.get(srx, [])
        if not srr_list:
            unmapped_srx.append(rec)
            continue
        for srr in srr_list:
            source_name = _sanitize_name(
                f"{rec['study_name']}_{rec['gsm']}_{rec['corrected_type']}"
            )
            all_rows.append({
                "Run": srr,
                "source_name": source_name,
                "organism": rec["organism"],
                "srx": srx,
                "gsm": rec["gsm"],
                "corrected_type": rec["corrected_type"],
                "study_name": rec["study_name"],
                "status": "OK",
            })

    # 处理缺失 SRX 的记录 (通过 ERS 查询)
    unmapped_ers = []
    for rec in missing_srx:
        ers = rec.get("ers", "")
        if ers:
            srr_list = ers_mapping.get(ers, [])
            if srr_list:
                for srr in srr_list:
                    source_name = _sanitize_name(
                        f"{rec['study_name']}_{rec['gsm']}_{rec['corrected_type']}"
                    )
                    all_rows.append({
                        "Run": srr,
                        "source_name": source_name,
                        "organism": rec["organism"],
                        "srx": f"(via {ers})",
                        "gsm": rec["gsm"],
                        "corrected_type": rec["corrected_type"],
                        "study_name": rec["study_name"],
                        "status": "OK (via ERS)",
                    })
            else:
                unmapped_ers.append(rec)
        else:
            unmapped_ers.append(rec)

    # 去重（基于 Run accession），保留首次出现的行
    seen_runs = set()
    unique_rows = []
    for row in all_rows:
        if row["Run"] not in seen_runs:
            seen_runs.add(row["Run"])
            unique_rows.append(row)

    # 确保输出目录存在
    srr_list_file.parent.mkdir(parents=True, exist_ok=True)
    mapping_file.parent.mkdir(parents=True, exist_ok=True)

    # 写入 sradownloader_input.txt（CSV 格式，适配 sradownloader RunTable 模式）
    fieldnames = ["Run", "source_name", "organism", "srx", "gsm", "corrected_type", "study_name", "status"]
    with open(srr_list_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()
        for row in unique_rows:
            writer.writerow(row)
    log(f"已生成 sradownloader 输入文件: {srr_list_file}")
    log(f"  共 {len(unique_rows)} 个唯一 Run accession (CSV RunTable 格式)")

    # 写入 srx_to_srr_mapping.csv（完整映射表，含重复 Run）
    with open(mapping_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)
    log(f"已生成映射表: {mapping_file}")
    log(f"  共 {len(all_rows)} 条映射记录")

    # 报告未映射的记录
    if unmapped_srx:
        warn(f"有 {len(unmapped_srx)} 条 SRX/ERX 记录未能映射到 SRR:")
        for rec in unmapped_srx[:10]:
            warn(f"  行{rec['row']}: {rec['gsm']} → {rec['srx']}")
        if len(unmapped_srx) > 10:
            warn(f"  ... 以及另外 {len(unmapped_srx) - 10} 条")

    if unmapped_ers:
        warn(f"有 {len(unmapped_ers)} 条缺失 SRX 且 ERS 也未能映射的记录:")
        for rec in unmapped_ers:
            warn(f"  行{rec['row']}: {rec['gsm']}")

    return unique_rows, unmapped_srx, unmapped_ers


# ============================================================================
# Step 4: 调用 sradownloader 下载 FASTQ
# ============================================================================
def run_sradownloader(srr_list_file, fastq_outdir, sradownloader_bin):
    """调用 sradownloader 工具下载 FASTQ 文件到 CCDS data/raw/fastq/。"""
    fastq_outdir.mkdir(parents=True, exist_ok=True)

    if not sradownloader_bin.exists():
        warn(f"sradownloader 工具未找到: {sradownloader_bin}")
        warn(f"请将 sradownloader 放置于: {sradownloader_bin}")
        warn(f"或手动运行: sradownloader --outdir {fastq_outdir} {srr_list_file}")
        return False

    if not srr_list_file.exists():
        warn(f"输入文件不存在: {srr_list_file}")
        warn("请先运行 --prepare 模式生成输入文件")
        return False

    cmd = [
        str(sradownloader_bin),
        "--outdir", str(fastq_outdir),
        str(srr_list_file),
    ]
    log(f"执行下载命令: {' '.join(cmd)}")
    log(f"FASTQ 输出目录: {fastq_outdir}")

    # 创建日志目录
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / "sradownloader.log"

    try:
        with open(log_file, "w") as lf:
            proc = subprocess.run(
                cmd,
                stdout=lf,
                stderr=subprocess.STDOUT,
                cwd=str(PROJECT_ROOT),
            )
        if proc.returncode == 0:
            log(f"下载完成！日志: {log_file}")
            return True
        else:
            warn(f"sradownloader 退出码: {proc.returncode}，请查看日志: {log_file}")
            return False
    except Exception as e:
        warn(f"下载执行失败: {e}")
        return False


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="SRA 数据下载与 SRX→SRR 映射工具 (CCDS 版)"
    )
    parser.add_argument(
        "--prepare", action="store_true",
        help="从 metadata.csv 生成 SRX→SRR 映射与 sradownloader 输入文件",
    )
    parser.add_argument(
        "--download", action="store_true",
        help="调用 sradownloader 下载 FASTQ 到 data/raw/fastq/",
    )
    parser.add_argument(
        "--metadata", type=str, default=str(METADATA_CSV),
        help=f"元数据 CSV 路径 (默认: {METADATA_CSV})",
    )
    parser.add_argument(
        "--outdir", type=str, default=str(FASTQ_OUTPUT_DIR),
        help=f"FASTQ 输出目录 (默认: {FASTQ_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--sradownloader", type=str, default=str(SRADOWNLOADER_BIN),
        help=f"sradownloader 可执行文件路径 (默认: {SRADOWNLOADER_BIN})",
    )

    args = parser.parse_args()

    if not args.prepare and not args.download:
        parser.print_help()
        print("\n错误: 请至少指定 --prepare 或 --download 之一", file=sys.stderr)
        sys.exit(1)

    metadata_path = Path(args.metadata)
    fastq_outdir = Path(args.outdir)
    sradownloader_bin = Path(args.sradownloader)

    if args.prepare:
        log("=" * 60)
        log("sradownloader 输入文件准备工具 (CCDS 版)")
        log("=" * 60)
        log(f"  项目根目录: {PROJECT_ROOT}")
        log(f"  元数据文件: {metadata_path}")
        log(f"  映射表输出: {MAPPING_FILE}")
        log(f"  输入文件:   {SRR_LIST_FILE}")
        log(f"  FASTQ 目标: {fastq_outdir}")

        # Step 1: 提取 accession
        if not metadata_path.exists():
            warn(f"元数据文件不存在: {metadata_path}")
            warn(f"请将 metadata.csv 放置于: {METADATA_CSV}")
            sys.exit(1)

        log(f"\n[Step 1] 从 {metadata_path} 提取 accession ...")
        experiments, missing_srx = extract_accessions(str(metadata_path))
        log(f"  有效 SRX/ERX 记录: {len(experiments)}")
        log(f"  缺失 SRX 的记录: {len(missing_srx)}")

        # 收集所有唯一的 SRX/ERX
        unique_srx = list(dict.fromkeys(rec["srx"] for rec in experiments if rec["srx"]))
        log(f"  唯一 SRX/ERX accession: {len(unique_srx)}")

        # Step 2: Entrez 批量查询
        log(f"\n[Step 2] 通过 NCBI Entrez 批量查询 SRX/ERX → SRR/ERR ...")
        if NCBI_API_KEY:
            log("  检测到 NCBI_API_KEY，使用增强速率限制 (10次/秒)")
        else:
            log("  未设置 NCBI_API_KEY，使用默认速率限制 (3次/秒)")
            log("  提示: export NCBI_API_KEY='your_key' 可加速查询")

        srx_mapping = entrez_srx_to_srr(unique_srx)
        log(f"  已映射 {len(srx_mapping)} 个 SRX/ERX")

        # Step 2b: 处理缺失 SRX 的 ERS 记录
        ers_list = [rec["ers"] for rec in missing_srx if "ers" in rec]
        ers_mapping = entrez_ers_to_srr(ers_list)

        # Step 3: 生成输出
        log(f"\n[Step 3] 生成输出文件 ...")
        unique_rows, unmapped_srx, unmapped_ers = write_outputs(
            experiments, missing_srx, srx_mapping, ers_mapping,
            SRR_LIST_FILE, MAPPING_FILE,
        )

        # 最终摘要
        log("\n" + "=" * 60)
        log("完成摘要")
        log("=" * 60)
        log(f"  输入记录总数:       {len(experiments) + len(missing_srx)}")
        log(f"  成功映射到 SRR:     {len(unique_rows)} 个唯一 Run accession")
        log(f"  未能映射:           {len(unmapped_srx) + len(unmapped_ers)} 条")
        log(f"\n  输入文件 (存档): {SRR_LIST_FILE}")
        log(f"  映射表 (存档):   {MAPPING_FILE}")

    if args.download:
        log("\n" + "=" * 60)
        log("启动 FASTQ 下载")
        log("=" * 60)
        success = run_sradownloader(SRR_LIST_FILE, fastq_outdir, sradownloader_bin)
        if success:
            log(f"\n下载完成！FASTQ 文件位于: {fastq_outdir}")
        else:
            warn("下载过程中出现问题，请检查日志")
            sys.exit(1)


if __name__ == "__main__":
    main()
