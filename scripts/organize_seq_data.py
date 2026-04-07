#!/usr/bin/env python3
"""
测序数据自动整理脚本

功能：根据 metadata.csv 和 srr_mapping.csv，将离线下载的测序数据
     按物种（organism）和实验类型（corrected_type）分类整理到标准目录结构中。

映射链路：
  文件名 -> SRR/ERR ID -> (srr_mapping.csv) -> GSM ID -> (metadata.csv) -> organism + corrected_type

输入：
  - metadata.csv：全局元数据表（第1行注释，第2行列头），关键列：
      experiment_alias (GSM ID), organism, corrected_type
  - srr_mapping.csv：SRR/ERR 到 GSM 的映射表，列：
      experiment_alias, sra_accession
  - 离线数据包目录（sradownloader_output）：包含多个 sra_part_*_output 子目录

输出：
  - 按 [Organism]/[corrected_type]/ 分类的目录树（空格替换为下划线）
  - logs/data_organization.log：记录每个文件的源路径和目标路径
  - logs/unmatched_log.txt：无法匹配的异常数据清单

用法：
  python organize_seq_data.py [--dry-run] [--symlink|--hardlink|--copy|--move]
"""

import argparse
import csv
import logging
import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path


# ============================================================
# 配置区
# ============================================================
BASE_DIR = Path(__file__).resolve().parent.parent  # project/
METADATA_CSV = BASE_DIR / "data" / "raw" / "TE_ribo-seq" / "metadata.csv"
SRR_MAPPING_CSV = BASE_DIR / "data" / "raw" / "TE_ribo-seq" / "srr_mapping.csv"
RAW_DATA_DIR = BASE_DIR / "data" / "raw" / "TE_ribo-seq" / "sradownloader_output"
OUTPUT_DIR = BASE_DIR / "data" / "organized"
LOG_DIR = BASE_DIR / "logs"

# 需要处理的数据文件扩展名
DATA_EXTENSIONS = {".fastq.gz", ".fq.gz", ".bam", ".sra"}

# 物种名标准化映射：俗名/缩写/拼写错误 -> 标准拉丁学名
SPECIES_NORMALIZATION_MAP = {
    # 俗名
    "Human":                        "Homo sapiens",
    "human":                        "Homo sapiens",
    "Mouse":                        "Mus musculus",
    "mouse":                        "Mus musculus",
    "Rat":                          "Rattus norvegicus",
    "rat":                          "Rattus norvegicus",
    "Soybean":                      "Glycine max",
    "soybean":                      "Glycine max",
    # 缩写
    "Arabidopsis":                  "Arabidopsis thaliana",
    "D. melanogaster":              "Drosophila melanogaster",
    # 拼写错误
    "Caenorhabitis brenneri":       "Caenorhabditis brenneri",
    "Caenorhaboditis elegans":      "Caenorhabditis elegans",
    # 大小写不规范
    "danio rerio":                  "Danio rerio",
    # 特殊混合名（保留为联合标注）
    "Saccharomyces cerevisiae* Saccharomyces paradoxus":
        "Saccharomyces cerevisiae x Saccharomyces paradoxus",
}


def normalize_organism(raw: str) -> str:
    """
    物种名标准化清洗链：
      1. 去除首尾空白
      2. 查映射字典替换俗名/缩写/拼写错误
      3. 空格替换为下划线
      4. 首字母大写（修复纯小写），其余保持原样
    """
    name = raw.strip()
    if not name:
        return name
    # 查映射字典
    name = SPECIES_NORMALIZATION_MAP.get(name, name)
    # 空格 -> 下划线
    name = name.replace(" ", "_")
    # 首字母大写（仅修正第一个字符，保留学名中种加词的小写）
    if name[0].islower():
        name = name[0].upper() + name[1:]
    return name


# ============================================================
# 工具函数
# ============================================================
def setup_logging(log_dir: Path) -> logging.Logger:
    """配置日志系统，同时输出到文件和控制台。"""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "data_organization.log"

    logger = logging.getLogger("organize_seq_data")
    logger.setLevel(logging.DEBUG)

    # 文件 handler
    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    ))

    # 控制台 handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(levelname)-7s | %(message)s"))

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def load_srr_mapping(csv_path: Path, logger: logging.Logger) -> dict:
    """
    读取 srr_mapping.csv，构建 sra_accession -> experiment_alias 的映射。
    返回: {SRR_ID: GSM_ID, ...}
    """
    mapping = {}
    if not csv_path.exists():
        logger.error(f"SRR 映射文件不存在: {csv_path}")
        sys.exit(1)

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = row.get("sra_accession", "").strip()
            gsm = row.get("experiment_alias", "").strip()
            if srr and gsm:
                mapping[srr] = gsm

    logger.info(f"从 srr_mapping.csv 加载了 {len(mapping)} 条 SRR->GSM 映射")
    return mapping


def load_metadata(csv_path: Path, logger: logging.Logger) -> dict:
    """
    读取 metadata.csv，构建 experiment_alias -> {organism, corrected_type} 的字典。

    注意：metadata.csv 第 1 行是注释行（"Curated Data,..."），第 2 行是真正的列头。
    同时将 matched_RNA-seq_experiment_alias 反向映射，使 RNA-Seq 对应的 GSM 也能查到。
    """
    metadata = {}
    if not csv_path.exists():
        logger.error(f"元数据文件不存在: {csv_path}")
        sys.exit(1)

    with open(csv_path, "r", encoding="utf-8") as f:
        # 跳过第 1 行注释
        next(f)
        reader = csv.DictReader(f)

        if "experiment_alias" not in reader.fieldnames:
            logger.error(f"metadata.csv 缺少 'experiment_alias' 列。现有列: {reader.fieldnames}")
            sys.exit(1)

        rna_reverse_count = 0
        for row in reader:
            alias = row.get("experiment_alias", "").strip()
            organism = row.get("organism", "").strip()
            corrected_type = row.get("corrected_type", "").strip()
            matched_rna = row.get("matched_RNA-seq_experiment_alias", "").strip()

            # 物种名标准化
            organism = normalize_organism(organism)

            if alias:
                metadata[alias] = {
                    "organism": organism,
                    "corrected_type": corrected_type,
                }

            # 反向映射：RNA-Seq 对应的 GSM -> 同物种 + RNA-Seq
            if matched_rna and matched_rna != "NA" and matched_rna not in metadata:
                metadata[matched_rna] = {
                    "organism": organism,
                    "corrected_type": "RNA-Seq",
                }
                rna_reverse_count += 1

    logger.info(f"从 metadata.csv 加载了 {len(metadata)} 条样本记录"
                f"（直接 {len(metadata) - rna_reverse_count} 条 + RNA-Seq反向 {rna_reverse_count} 条）")
    return metadata


def get_file_extension(filename: str) -> str:
    """获取文件的完整扩展名（如 .fastq.gz）。"""
    name = filename
    exts = []
    while True:
        name, ext = os.path.splitext(name)
        if ext:
            exts.append(ext)
        else:
            break
    exts.reverse()
    return "".join(exts)


def is_data_file(filename: str) -> bool:
    """判断文件是否为测序数据文件。"""
    full_ext = get_file_extension(filename)
    return full_ext in DATA_EXTENSIONS


def extract_run_id(filename: str) -> str:
    """从文件名中提取 SRR/ERR/DRR accession ID。"""
    match = re.match(r"((?:SRR|ERR|DRR)\d+)", filename)
    return match.group(1) if match else ""


def sanitize_dirname(name: str) -> str:
    """将物种名等字符串转为合法目录名（空格替换为下划线）。"""
    return re.sub(r"\s+", "_", name.strip())


def build_target_filename(experiment_alias: str, corrected_type: str, original_filename: str) -> str:
    """
    构建标准化目标文件名：{experiment_alias}_{corrected_type}_{SRR_ID}{_1/_2}.[原扩展名]

    保留原始 SRR/ERR ID（同一 GSM 可能有多个 run）和 paired-end 后缀。
    """
    full_ext = get_file_extension(original_filename)
    basename = original_filename
    for _ in range(full_ext.count(".")):
        basename = os.path.splitext(basename)[0]

    run_id = extract_run_id(basename)

    # 提取 paired-end 后缀 (_1 or _2)
    pe_match = re.search(r"_([12])$", basename)
    pe_suffix = f"_{pe_match.group(1)}" if pe_match else ""

    safe_type = sanitize_dirname(corrected_type)

    if run_id:
        return f"{experiment_alias}_{safe_type}_{run_id}{pe_suffix}{full_ext}"
    else:
        return f"{experiment_alias}_{safe_type}{pe_suffix}{full_ext}"


def parse_sradownloader_results(txt_path: Path) -> dict:
    """解析 sradownloader_results.txt，返回 {accession: status} 字典。"""
    results = {}
    if not txt_path.exists():
        return results
    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("ACCESSION"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                results[parts[0].strip()] = parts[1].strip()
    return results


def transfer_file(src: Path, dst: Path, method: str, logger: logging.Logger) -> bool:
    """根据指定方法将文件从 src 转移到 dst。"""
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)

        if dst.exists():
            logger.warning(f"目标文件已存在，跳过: {dst}")
            return False

        if method == "symlink":
            dst.symlink_to(src.resolve())
        elif method == "hardlink":
            os.link(src, dst)
        elif method == "copy":
            shutil.copy2(src, dst)
        elif method == "move":
            shutil.move(str(src), str(dst))
        else:
            logger.error(f"未知的传输方法: {method}")
            return False

        return True
    except Exception as e:
        logger.error(f"文件传输失败 [{method}] {src} -> {dst}: {e}")
        return False


# ============================================================
# 主逻辑
# ============================================================
def organize_data(
    srr_mapping: dict,
    metadata: dict,
    raw_data_dir: Path,
    output_dir: Path,
    transfer_method: str,
    dry_run: bool,
    logger: logging.Logger,
) -> dict:
    """
    遍历原始数据目录，按 metadata 信息分类整理文件。

    映射链路：文件名 -> SRR ID -> GSM ID -> {organism, corrected_type}
    """
    stats = {
        "total_files": 0,
        "matched": 0,
        "unmatched": 0,
        "skipped_non_data": 0,
        "transfer_errors": 0,
        "already_exists": 0,
    }
    unmatched_records = []

    if not raw_data_dir.exists():
        logger.error(f"原始数据目录不存在: {raw_data_dir}")
        sys.exit(1)

    # 遍历所有 sra_part_*_output 子目录
    part_dirs = sorted([d for d in raw_data_dir.iterdir() if d.is_dir()])
    if not part_dirs:
        logger.warning(f"在 {raw_data_dir} 下未找到任何子目录")
        return stats

    logger.info(f"发现 {len(part_dirs)} 个数据包目录")

    for part_dir in part_dirs:
        logger.info(f"--- 处理数据包: {part_dir.name} ---")

        # 检查 sradownloader_results.txt
        results_txt = part_dir / "sradownloader_results.txt"
        if not results_txt.exists():
            unmatched_records.append({
                "type": "缺失说明文件",
                "path": str(part_dir),
                "detail": "该目录下未找到 sradownloader_results.txt",
            })
            logger.warning(f"目录 {part_dir.name} 缺少 sradownloader_results.txt")

        # 遍历目录中的所有文件
        for file_path in sorted(part_dir.iterdir()):
            if not file_path.is_file():
                continue

            filename = file_path.name

            # 跳过非数据文件
            if not is_data_file(filename):
                stats["skipped_non_data"] += 1
                continue

            stats["total_files"] += 1

            # 第一步：从文件名提取 SRR/ERR ID
            run_id = extract_run_id(filename)
            if not run_id:
                stats["unmatched"] += 1
                unmatched_records.append({
                    "type": "无法提取Run_ID",
                    "path": str(file_path),
                    "detail": f"无法从文件名 '{filename}' 中提取 SRR/ERR/DRR ID",
                })
                continue

            # 第二步：SRR -> GSM（通过 srr_mapping）
            gsm_id = srr_mapping.get(run_id)
            if not gsm_id:
                stats["unmatched"] += 1
                unmatched_records.append({
                    "type": "srr_mapping无匹配",
                    "path": str(file_path),
                    "detail": f"Run ID '{run_id}' 在 srr_mapping.csv 中未找到",
                })
                continue

            # 第三步：GSM -> organism + corrected_type（通过 metadata）
            info = metadata.get(gsm_id)
            if not info:
                stats["unmatched"] += 1
                unmatched_records.append({
                    "type": "metadata无匹配",
                    "path": str(file_path),
                    "detail": f"GSM ID '{gsm_id}'（来自 {run_id}）在 metadata.csv 中未找到",
                })
                continue

            organism = info["organism"]
            corrected_type = info["corrected_type"]

            if not organism or not corrected_type:
                stats["unmatched"] += 1
                unmatched_records.append({
                    "type": "metadata信息不完整",
                    "path": str(file_path),
                    "detail": f"GSM '{gsm_id}': organism='{organism}', corrected_type='{corrected_type}'",
                })
                continue

            # 构建目标路径
            organism_dir = sanitize_dirname(organism)
            type_dir = sanitize_dirname(corrected_type)
            target_filename = build_target_filename(gsm_id, corrected_type, filename)
            target_path = output_dir / organism_dir / type_dir / target_filename

            if dry_run:
                logger.info(f"[DRY-RUN] {file_path} -> {target_path}")
                stats["matched"] += 1
            else:
                success = transfer_file(file_path, target_path, transfer_method, logger)
                if success:
                    stats["matched"] += 1
                    logger.info(f"[{transfer_method.upper()}] {file_path} -> {target_path}")
                elif target_path.exists():
                    stats["already_exists"] += 1
                else:
                    stats["transfer_errors"] += 1

    # 写入 unmatched_log.txt
    unmatched_log_path = LOG_DIR / "unmatched_log.txt"
    unmatched_log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(unmatched_log_path, "w", encoding="utf-8") as f:
        f.write(f"# 未匹配数据清单 - 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# 共 {len(unmatched_records)} 条记录\n")
        f.write("=" * 100 + "\n\n")
        for rec in unmatched_records:
            f.write(f"类型: {rec['type']}\n")
            f.write(f"路径: {rec['path']}\n")
            f.write(f"详情: {rec['detail']}\n")
            f.write("-" * 80 + "\n")

    logger.info(f"未匹配记录已写入: {unmatched_log_path}")
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="测序数据自动整理脚本：按物种和实验类型分类",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    transfer_group = parser.add_mutually_exclusive_group()
    transfer_group.add_argument(
        "--symlink", action="store_const", const="symlink", dest="method",
        help="创建符号链接（默认）",
    )
    transfer_group.add_argument(
        "--hardlink", action="store_const", const="hardlink", dest="method",
        help="创建硬链接",
    )
    transfer_group.add_argument(
        "--copy", action="store_const", const="copy", dest="method",
        help="复制文件",
    )
    transfer_group.add_argument(
        "--move", action="store_const", const="move", dest="method",
        help="移动文件（谨慎使用！）",
    )

    parser.add_argument(
        "--dry-run", action="store_true",
        help="试运行模式，只打印操作不实际执行",
    )
    parser.add_argument(
        "--metadata", type=str, default=str(METADATA_CSV),
        help=f"metadata.csv 路径（默认: {METADATA_CSV}）",
    )
    parser.add_argument(
        "--srr-mapping", type=str, default=str(SRR_MAPPING_CSV),
        help=f"srr_mapping.csv 路径（默认: {SRR_MAPPING_CSV}）",
    )
    parser.add_argument(
        "--input-dir", type=str, default=str(RAW_DATA_DIR),
        help=f"原始数据根目录（默认: {RAW_DATA_DIR}）",
    )
    parser.add_argument(
        "--output-dir", type=str, default=str(OUTPUT_DIR),
        help=f"输出目录（默认: {OUTPUT_DIR}）",
    )

    parser.set_defaults(method="symlink")
    args = parser.parse_args()

    # 初始化日志
    logger = setup_logging(LOG_DIR)
    logger.info("=" * 60)
    logger.info("测序数据自动整理脚本启动")
    logger.info(f"运行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"元数据文件:   {args.metadata}")
    logger.info(f"SRR映射文件:  {args.srr_mapping}")
    logger.info(f"输入目录:     {args.input_dir}")
    logger.info(f"输出目录:     {args.output_dir}")
    logger.info(f"传输方式:     {args.method}")
    logger.info(f"试运行模式:   {args.dry_run}")
    logger.info("=" * 60)

    # 加载映射表
    srr_mapping = load_srr_mapping(Path(args.srr_mapping), logger)
    metadata = load_metadata(Path(args.metadata), logger)

    # 执行整理
    stats = organize_data(
        srr_mapping=srr_mapping,
        metadata=metadata,
        raw_data_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        transfer_method=args.method,
        dry_run=args.dry_run,
        logger=logger,
    )

    # 输出统计摘要
    logger.info("=" * 60)
    logger.info("整理完成！统计摘要：")
    logger.info(f"  数据文件总数:       {stats['total_files']}")
    logger.info(f"  成功匹配并处理:     {stats['matched']}")
    logger.info(f"  未匹配:             {stats['unmatched']}")
    logger.info(f"  目标已存在(跳过):   {stats['already_exists']}")
    logger.info(f"  传输错误:           {stats['transfer_errors']}")
    logger.info(f"  跳过的非数据文件:   {stats['skipped_non_data']}")
    logger.info("=" * 60)

    if stats["unmatched"] > 0:
        logger.warning(f"有 {stats['unmatched']} 个文件未能匹配，请查看 logs/unmatched_log.txt")


if __name__ == "__main__":
    main()
