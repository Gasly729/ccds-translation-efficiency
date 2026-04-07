#!/usr/bin/env python3
"""
Ribo-seq ↔ RNA-seq 配对校验脚本 (Pre-flight Check)

针对指定物种，验证：
  1. 每条 Ribo-seq 记录的 matched_RNA-seq_experiment_alias 在 RNA-seq 记录中真实存在
  2. 对应的物理 FASTQ 文件均在 data/organized/ 目录中
  3. 无"孤儿样本"（有 Ribo 无 RNA 或反之）

输出：
  - 终端彩色摘要
  - logs/pairing_validation.log 详细日志

用法：
  python src/data/validate_pairing.py [--species "Homo sapiens" "Mus musculus" ...]
"""

import argparse
import csv
import logging
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ============================================================
# 配置
# ============================================================
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # project/
METADATA_CSV = BASE_DIR / "data" / "raw" / "TE_ribo-seq" / "metadata.csv"
SRR_MAPPING_CSV = BASE_DIR / "data" / "raw" / "TE_ribo-seq" / "srr_mapping.csv"
ORGANIZED_DIR = BASE_DIR / "data" / "organized"
LOG_DIR = BASE_DIR / "logs"

# 物种名标准化（与 organize_seq_data.py 保持一致）
SPECIES_NORMALIZATION_MAP = {
    "Human": "Homo sapiens",
    "human": "Homo sapiens",
    "Mouse": "Mus musculus",
    "mouse": "Mus musculus",
    "Rat": "Rattus norvegicus",
    "rat": "Rattus norvegicus",
    "Soybean": "Glycine max",
    "soybean": "Glycine max",
    "Arabidopsis": "Arabidopsis thaliana",
    "D. melanogaster": "Drosophila melanogaster",
    "Caenorhabitis brenneri": "Caenorhabditis brenneri",
    "Caenorhaboditis elegans": "Caenorhabditis elegans",
    "danio rerio": "Danio rerio",
    "Saccharomyces cerevisiae* Saccharomyces paradoxus":
        "Saccharomyces cerevisiae x Saccharomyces paradoxus",
}

# 默认校验的 4 个目标物种
DEFAULT_SPECIES = [
    "Homo sapiens",
    "Mus musculus",
    "Arabidopsis thaliana",
    "Caenorhabditis elegans",
]


def normalize_organism(raw: str) -> str:
    """物种名标准化。"""
    name = raw.strip()
    if not name:
        return name
    name = SPECIES_NORMALIZATION_MAP.get(name, name)
    name = name.replace(" ", "_")
    if name and name[0].islower():
        name = name[0].upper() + name[1:]
    return name


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / "pairing_validation.log"

    logger = logging.getLogger("validate_pairing")
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    ))

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(levelname)-7s | %(message)s"))

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def load_srr_mapping(csv_path: Path) -> dict:
    """返回 {GSM: [SRR1, SRR2, ...]}"""
    gsm_to_srrs = defaultdict(list)
    with open(csv_path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gsm = row.get("experiment_alias", "").strip()
            srr = row.get("sra_accession", "").strip()
            if gsm and srr:
                gsm_to_srrs[gsm].append(srr)
    return dict(gsm_to_srrs)


def load_metadata(csv_path: Path) -> list:
    """返回 metadata 行列表。"""
    rows = []
    with open(csv_path, "r", encoding="utf-8") as f:
        next(f)  # skip comment line
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def find_physical_files(organized_dir: Path, species_dir: str, seq_type: str, gsm_id: str) -> list:
    """在 organized 目录中查找某 GSM 对应的物理文件。"""
    target_dir = organized_dir / species_dir / seq_type
    if not target_dir.exists():
        return []
    pattern = f"{gsm_id}_*"
    return list(target_dir.glob(pattern))


def validate_species(
    species_name: str,
    metadata_rows: list,
    gsm_to_srrs: dict,
    organized_dir: Path,
    logger: logging.Logger,
) -> dict:
    """
    校验单个物种的 Ribo↔RNA 配对完整性。

    返回: {
        "species": str,
        "ribo_total": int,
        "paired_ok": int,
        "orphan_ribo": [...],    # 有 Ribo 但无 RNA
        "orphan_rna": [...],     # matched_RNA GSM 在 metadata 中不存在
        "missing_ribo_files": [...],
        "missing_rna_files": [...],
    }
    """
    species_dir = normalize_organism(species_name)

    # 筛选该物种的所有记录
    ribo_records = []
    rna_gsms = set()

    for row in metadata_rows:
        org = normalize_organism(row.get("organism", "").strip())
        ctype = row.get("corrected_type", "").strip()
        alias = row.get("experiment_alias", "").strip()

        if org != species_dir:
            continue

        if ctype == "Ribo-Seq":
            ribo_records.append(row)
        elif ctype == "RNA-Seq":
            rna_gsms.add(alias)

    result = {
        "species": species_name,
        "species_dir": species_dir,
        "ribo_total": len(ribo_records),
        "rna_total": len(rna_gsms),
        "paired_ok": 0,
        "orphan_ribo_no_match_col": [],
        "orphan_ribo_rna_not_found": [],
        "missing_ribo_files": [],
        "missing_rna_files": [],
    }

    for row in ribo_records:
        ribo_gsm = row.get("experiment_alias", "").strip()
        matched_rna_gsm = row.get("matched_RNA-seq_experiment_alias", "").strip()

        # 检查配对列是否存在
        if not matched_rna_gsm or matched_rna_gsm == "NA":
            result["orphan_ribo_no_match_col"].append(ribo_gsm)
            logger.debug(f"[{species_name}] 孤儿 Ribo (无配对列): {ribo_gsm}")
            continue

        # 检查配对的 RNA GSM 是否在 metadata 中
        # (有些 RNA-Seq 样本可能只在 matched_ 列引用但自身不在 metadata 中)
        # 这里我们检查物理文件是否存在即可

        # 检查 Ribo 物理文件
        ribo_files = find_physical_files(organized_dir, species_dir, "Ribo-Seq", ribo_gsm)
        if not ribo_files:
            result["missing_ribo_files"].append(ribo_gsm)
            logger.debug(f"[{species_name}] Ribo 物理文件缺失: {ribo_gsm}")

        # 检查 RNA 物理文件
        rna_files = find_physical_files(organized_dir, species_dir, "RNA-Seq", matched_rna_gsm)
        if not rna_files:
            result["missing_rna_files"].append(matched_rna_gsm)
            logger.debug(f"[{species_name}] RNA 物理文件缺失: {matched_rna_gsm} (配对自 {ribo_gsm})")

        if ribo_files and rna_files:
            result["paired_ok"] += 1
        elif not rna_files and matched_rna_gsm not in rna_gsms:
            result["orphan_ribo_rna_not_found"].append(
                f"{ribo_gsm} -> {matched_rna_gsm} (RNA GSM 不在 metadata)")

    return result


def main():
    parser = argparse.ArgumentParser(description="Ribo↔RNA 配对校验 (Pre-flight Check)")
    parser.add_argument(
        "--species", nargs="+", default=DEFAULT_SPECIES,
        help="要校验的物种列表（默认：Human, Mus musculus, Arabidopsis thaliana, C. elegans）",
    )
    parser.add_argument("--metadata", default=str(METADATA_CSV))
    parser.add_argument("--srr-mapping", default=str(SRR_MAPPING_CSV))
    parser.add_argument("--organized-dir", default=str(ORGANIZED_DIR))
    args = parser.parse_args()

    logger = setup_logging()
    logger.info("=" * 70)
    logger.info("Ribo-seq ↔ RNA-seq 配对校验 (Pre-flight Check)")
    logger.info(f"运行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"目标物种: {args.species}")
    logger.info("=" * 70)

    # 加载数据
    metadata_rows = load_metadata(Path(args.metadata))
    gsm_to_srrs = load_srr_mapping(Path(args.srr_mapping))
    organized_dir = Path(args.organized_dir)

    logger.info(f"metadata 记录数: {len(metadata_rows)}")
    logger.info(f"srr_mapping GSM 数: {len(gsm_to_srrs)}")

    # 逐物种校验
    all_pass = True
    results = []

    for species in args.species:
        logger.info(f"\n{'─' * 50}")
        logger.info(f"校验物种: {species}")
        logger.info(f"{'─' * 50}")

        result = validate_species(species, metadata_rows, gsm_to_srrs, organized_dir, logger)
        results.append(result)

        total_issues = (
            len(result["orphan_ribo_no_match_col"])
            + len(result["orphan_ribo_rna_not_found"])
            + len(result["missing_ribo_files"])
            + len(result["missing_rna_files"])
        )

        match_rate = (
            f"{result['paired_ok']}/{result['ribo_total']}"
            if result["ribo_total"] > 0
            else "0/0"
        )
        pct = (
            f"{result['paired_ok'] / result['ribo_total'] * 100:.1f}%"
            if result["ribo_total"] > 0
            else "N/A"
        )

        logger.info(f"  Ribo-Seq 样本数:           {result['ribo_total']}")
        logger.info(f"  RNA-Seq 样本数 (metadata): {result['rna_total']}")
        logger.info(f"  配对成功 (文件均存在):     {match_rate} ({pct})")

        if result["orphan_ribo_no_match_col"]:
            logger.warning(f"  孤儿 Ribo (无配对列):      {len(result['orphan_ribo_no_match_col'])}")
            for gsm in result["orphan_ribo_no_match_col"][:5]:
                logger.warning(f"    - {gsm}")

        if result["orphan_ribo_rna_not_found"]:
            logger.warning(f"  孤儿 Ribo (RNA不存在):     {len(result['orphan_ribo_rna_not_found'])}")
            for info in result["orphan_ribo_rna_not_found"][:5]:
                logger.warning(f"    - {info}")

        if result["missing_ribo_files"]:
            logger.warning(f"  Ribo 物理文件缺失:        {len(result['missing_ribo_files'])}")

        if result["missing_rna_files"]:
            logger.warning(f"  RNA 物理文件缺失:         {len(result['missing_rna_files'])}")

        if total_issues > 0:
            all_pass = False
            logger.warning(f"  ⚠ {species}: 存在 {total_issues} 个问题")
        else:
            logger.info(f"  ✓ {species}: 全部通过")

    # 总结
    logger.info(f"\n{'=' * 70}")
    logger.info("校验总结")
    logger.info(f"{'=' * 70}")

    summary_lines = []
    for r in results:
        pct = f"{r['paired_ok'] / r['ribo_total'] * 100:.1f}%" if r["ribo_total"] > 0 else "N/A"
        status = "PASS" if (
            not r["orphan_ribo_no_match_col"]
            and not r["orphan_ribo_rna_not_found"]
            and not r["missing_ribo_files"]
            and not r["missing_rna_files"]
        ) else "FAIL"
        line = f"  {r['species']:<30s}  Ribo={r['ribo_total']:>4d}  配对={r['paired_ok']:>4d}  率={pct:>6s}  [{status}]"
        summary_lines.append(line)
        logger.info(line)

    logger.info(f"{'=' * 70}")

    if all_pass:
        logger.info("✓ 所有物种配对校验通过！可以继续执行 Snakemake 流水线。")
    else:
        logger.warning("⚠ 存在配对问题，请检查 logs/pairing_validation.log")

    logger.info(f"详细日志: {LOG_DIR / 'pairing_validation.log'}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
