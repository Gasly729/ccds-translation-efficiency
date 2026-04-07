#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
te_calculator.py — 翻译效率 (Translational Efficiency, TE) 计算模块
=====================================================================
整合 CenikLab TE_model 原版核心算法与 CCDS 工程外壳。

核心流程 (4 阶段):
  Stage 0: 从 .ribo 文件提取 CDS 原始计数 → ribo_raw.csv, rnaseq_raw.csv
  Stage 1: CPM 标准化 → 低表达基因过滤 → dummy gene 合并 → 非 polyA 基因移除
  Stage 2: 调用 TE.R 执行 CLR/ILR 成分回归，TE = 回归残差的 CLR 表示
  Stage 3: 后处理，转置汇总，输出 te_results_final.csv

核心算法声明:
  Stage 1–2 的全部数学计算、过滤阈值与标准化策略 100% 严格复刻自
  CenikLab TE_model 原文代码 (ribobase_counts_processing.py + TE.R)。
  未做任何 AI "优化"或篡改。学术复现责任归属于原作者。

用法:
    python -m src.te_calc.te_calculator \\
        --ribo_dir data/processed/ribo_files \\
        --output_dir data/processed \\
        [--cpm_cutoff 1] [--overall_cutoff 70] \\
        [--nonpolya_csv data/external/nonpolyA_gene.csv] \\
        [--skip_r] [--skip_extract]
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from bioinfokit.analys import norm as bioinfokit_norm
    HAS_BIOINFOKIT = True
except ImportError:
    HAS_BIOINFOKIT = False

try:
    from ribopy import Ribo
    HAS_RIBOPY = True
except ImportError:
    HAS_RIBOPY = False


# ── 常量 ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent          # src/te_calc/
TE_R_PATH = SCRIPT_DIR / "TE.R"                       # 同目录下嵌入的 R 脚本
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_CPM_CUTOFF = 1
DEFAULT_OVERALL_CUTOFF = 70


# =============================================================================
# Stage 0: 从 .ribo 文件提取 CDS 原始计数
# =============================================================================
# I/O 外壳采用 Reference B (用户 RiboTE_Pipeline) 的文件遍历逻辑。
# 计数提取参数 (region="CDS", sum_lengths=True, sum_references=False)
# 与 Reference A 一致。
# =============================================================================

def extract_gene_name(transcript_id: str) -> str:
    """从完整转录本ID中提取基因名 (APPRIS 格式: 第5字段)."""
    parts = transcript_id.split("|")
    if len(parts) >= 5:
        return parts[4]
    return transcript_id


def extract_counts_from_ribo(ribo_path: str) -> tuple:
    """从单个 .ribo 文件提取 Ribo-seq 和 RNA-seq 的 CDS 计数."""
    sample_name = os.path.basename(ribo_path).replace(".ribo", "")
    ribo = Ribo(ribo_path)
    experiments = ribo.experiments
    if len(experiments) == 0:
        raise ValueError(f"文件 {ribo_path} 中没有实验数据")
    experiment_name = experiments[0]

    # Ribo-seq CDS counts
    ribo_counts_df = ribo.get_region_counts(
        region_name="CDS",
        sum_lengths=True,
        sum_references=False,
    )
    ribo_counts_series = ribo_counts_df[experiment_name].rename(sample_name)

    # RNA-seq CDS counts (if available)
    rna_counts_series = None
    has_rna = False
    if ribo.has_rnaseq(experiment_name):
        has_rna = True
        rnaseq_df = ribo.get_rnaseq()
        rna_data = rnaseq_df.loc[experiment_name, "CDS"]
        rna_counts_series = rna_data.rename(sample_name)
    else:
        print(f"  [警告] 样本 '{sample_name}' 没有内置 RNA-seq 数据")

    return sample_name, ribo_counts_series, rna_counts_series, has_rna


def stage0_extract(ribo_dir: str, output_dir: str,
                   extract_gene_name_flag: bool = True) -> tuple:
    """
    Stage 0: 遍历 ribo_dir 下所有 .ribo 文件，提取 CDS 计数矩阵。

    输出:
      - {output_dir}/ribo_raw.csv
      - {output_dir}/rnaseq_raw.csv
      - {output_dir}/infor_filter.csv
    """
    import re

    if not HAS_RIBOPY:
        raise ImportError("ribopy 未安装，无法提取 .ribo 文件。请 pip install ribopy")

    ribo_dir = os.path.abspath(ribo_dir)
    os.makedirs(output_dir, exist_ok=True)

    ribo_files = sorted(glob.glob(os.path.join(ribo_dir, "*.ribo")))
    if not ribo_files:
        raise FileNotFoundError(f"在 {ribo_dir} 中未找到 .ribo 文件")
    print(f"[Stage 0] 找到 {len(ribo_files)} 个 .ribo 文件")

    ribo_counts_list = []
    rna_counts_list = []
    sample_names = []

    for ribo_path in ribo_files:
        sname = os.path.basename(ribo_path).replace(".ribo", "")
        print(f"  处理: {sname}")
        try:
            sname, ribo_s, rna_s, has_rna = extract_counts_from_ribo(ribo_path)
            ribo_counts_list.append(ribo_s)
            if has_rna and rna_s is not None:
                rna_counts_list.append(rna_s)
            sample_names.append(sname)
        except Exception as e:
            print(f"  [错误] {sname}: {e}")
            continue

    if not ribo_counts_list:
        raise RuntimeError("没有成功提取任何数据")

    df_ribo = pd.concat(ribo_counts_list, axis=1).fillna(0).astype(int)
    if extract_gene_name_flag:
        df_ribo.index = df_ribo.index.map(extract_gene_name)
        df_ribo = df_ribo.groupby(df_ribo.index).mean().astype(int)
    df_ribo.index.name = ""

    df_rna = None
    if rna_counts_list:
        df_rna = pd.concat(rna_counts_list, axis=1).fillna(0)
        if extract_gene_name_flag:
            df_rna.index = df_rna.index.map(extract_gene_name)
            df_rna = df_rna.groupby(df_rna.index).mean()
        df_rna.index.name = ""

    # 保存
    ribo_out = os.path.join(output_dir, "ribo_raw.csv")
    df_ribo.to_csv(ribo_out)
    print(f"  [OK] {ribo_out} ({df_ribo.shape[0]} genes × {df_ribo.shape[1]} samples)")

    if df_rna is not None:
        rna_out = os.path.join(output_dir, "rnaseq_raw.csv")
        df_rna.to_csv(rna_out)
        print(f"  [OK] {rna_out} ({df_rna.shape[0]} genes × {df_rna.shape[1]} samples)")

    # infor_filter.csv
    info_list = []
    for i, s in enumerate(sample_names, 1):
        cell_line = re.sub(r'_rep\d+$', '', s)
        info_list.append({"": i, "experiment_alias": s, "cell_line": cell_line})
    df_info = pd.DataFrame(info_list)
    info_out = os.path.join(output_dir, "infor_filter.csv")
    df_info.to_csv(info_out, index=False)
    print(f"  [OK] {info_out}")

    return df_ribo, df_rna


# =============================================================================
# Stage 0.5: 样本配对验证 — 强约束 Ribo↔RNA 映射
# =============================================================================
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  配对锚点链路:                                                         │
# │    metadata.csv: experiment_alias (GSM) → experiment_accession (SRX)   │
# │    metadata.csv: Ribo-Seq GSM → matched_RNA-seq_experiment_alias (GSM) │
# │    srx_to_srr_mapping.csv: SRX → SRR (Run accession)                  │
# │                                                                         │
# │  最终映射: { bio_sample: {ribo_srr: ..., rna_srr: ...} }              │
# │                                                                         │
# │  规则:                                                                  │
# │    - 仅 corrected_type == 'Ribo-Seq' 且 matched_RNA != 'NA' 的行有效  │
# │    - 无法追踪到 SRR 的配对将被 DROP 并记录日志                         │
# │    - 最终 ribo_raw.csv 和 rnaseq_raw.csv 的列必须按此映射严格对齐     │
# └─────────────────────────────────────────────────────────────────────────┘
# =============================================================================

def build_sample_pairing(metadata_csv: str, mapping_csv: str) -> dict:
    """
    从 metadata.csv 和 srx_to_srr_mapping.csv 构建确定性 Ribo↔RNA 配对字典。

    返回:
        {
            "GSM1234_ribo_label": {
                "ribo_gsm": "GSM1234",
                "rna_gsm":  "GSM5678",
                "ribo_srx": "SRX111",
                "rna_srx":  "SRX222",
                "ribo_srr": ["SRR111"],
                "rna_srr":  ["SRR222"],
                "cell_line": "HeLa",
                "organism":  "Homo sapiens",
            },
            ...
        }
    """
    import csv

    # ── Step A: 解析 metadata.csv ────────────────────────────────────────
    gsm_to_srx = {}      # GSM → SRX
    gsm_to_info = {}     # GSM → {type, matched, cell_line, organism}
    ribo_entries = []     # (gsm, srx, matched_gsm, cell_line, organism)

    with open(metadata_csv, "r", newline="") as f:
        reader = csv.reader(f)
        # 跳过描述行 (第1行是 "Curated Data,...")
        first_row = next(reader)
        # 第2行是真实表头
        header = next(reader)

        col_map = {h.strip(): i for i, h in enumerate(header)}
        idx_alias = col_map.get("experiment_alias")
        idx_matched = col_map.get("matched_RNA-seq_experiment_alias")
        idx_type = col_map.get("corrected_type")
        idx_srx = col_map.get("experiment_accession")
        idx_cell = col_map.get("cell_line")
        idx_org = col_map.get("organism")

        required = {"experiment_alias": idx_alias, "corrected_type": idx_type,
                    "experiment_accession": idx_srx,
                    "matched_RNA-seq_experiment_alias": idx_matched}
        missing = [k for k, v in required.items() if v is None]
        if missing:
            raise ValueError(
                f"metadata.csv 缺少必需列: {missing}\n"
                f"实际列名: {header}")

        for row in reader:
            if len(row) <= max(idx_alias, idx_srx, idx_type, idx_matched):
                continue
            alias = row[idx_alias].strip()
            srx = row[idx_srx].strip()
            ctype = row[idx_type].strip()
            matched = row[idx_matched].strip()
            cell = row[idx_cell].strip() if idx_cell is not None else ""
            org = row[idx_org].strip() if idx_org is not None else ""

            if alias and srx:
                gsm_to_srx[alias] = srx
                gsm_to_info[alias] = {
                    "type": ctype, "matched": matched,
                    "cell_line": cell, "organism": org,
                }

            if "Ribo" in ctype and matched and matched != "NA":
                ribo_entries.append((alias, srx, matched, cell, org))

    print(f"[配对] metadata.csv: {len(gsm_to_srx)} 个实验, "
          f"{len(ribo_entries)} 个 Ribo-Seq 有配对 RNA")

    # ── Step B: 解析 srx_to_srr_mapping.csv ──────────────────────────────
    srx_to_srr = {}   # SRX → [SRR, ...]
    gsm_to_srr = {}   # GSM → [SRR, ...]  (backup lookup)

    with open(mapping_csv, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        col_map_m = {h.strip(): i for i, h in enumerate(header)}
        idx_run = col_map_m.get("Run")
        idx_msrx = col_map_m.get("srx")
        idx_mgsm = col_map_m.get("gsm")

        if idx_run is None or idx_msrx is None:
            raise ValueError(
                f"srx_to_srr_mapping.csv 缺少 'Run' 或 'srx' 列\n"
                f"实际列名: {header}")

        for row in reader:
            srr = row[idx_run].strip()
            srx = row[idx_msrx].strip()
            gsm = row[idx_mgsm].strip() if idx_mgsm is not None else ""
            if srx and srr:
                srx_to_srr.setdefault(srx, []).append(srr)
            if gsm and srr:
                gsm_to_srr.setdefault(gsm, []).append(srr)

    print(f"[配对] srx_to_srr_mapping.csv: {len(srx_to_srr)} 个 SRX, "
          f"{sum(len(v) for v in srx_to_srr.values())} 个 SRR")

    # ── Step C: 构建配对字典 ─────────────────────────────────────────────
    pairing = {}
    dropped = []

    for ribo_gsm, ribo_srx, rna_gsm, cell_line, organism in ribo_entries:
        rna_srx = gsm_to_srx.get(rna_gsm)
        ribo_srrs = srx_to_srr.get(ribo_srx, gsm_to_srr.get(ribo_gsm, []))
        rna_srrs = srx_to_srr.get(rna_srx, gsm_to_srr.get(rna_gsm, [])) if rna_srx else []

        if not ribo_srrs or not rna_srrs:
            reason = ""
            if not ribo_srrs:
                reason += f"Ribo SRR 缺失 (SRX={ribo_srx}); "
            if not rna_srx:
                reason += f"RNA SRX 缺失 (matched_GSM={rna_gsm}); "
            elif not rna_srrs:
                reason += f"RNA SRR 缺失 (SRX={rna_srx}); "
            dropped.append((ribo_gsm, reason))
            continue

        pairing[ribo_gsm] = {
            "ribo_gsm": ribo_gsm,
            "rna_gsm": rna_gsm,
            "ribo_srx": ribo_srx,
            "rna_srx": rna_srx,
            "ribo_srr": sorted(ribo_srrs),
            "rna_srr": sorted(rna_srrs),
            "cell_line": cell_line,
            "organism": organism,
        }

    if dropped:
        print(f"[配对] 警告: {len(dropped)} 个 Ribo 样本因缺失 SRR 被 DROP:")
        for gsm, reason in dropped[:10]:
            print(f"  [DROP] {gsm}: {reason}")
        if len(dropped) > 10:
            print(f"  ... 及其余 {len(dropped) - 10} 个")

    print(f"[配对] 成功构建 {len(pairing)} 个 Ribo↔RNA 配对")
    return pairing


def print_pairing_report(pairing: dict) -> None:
    """
    防呆日志: 在终端输出完整的配对报告。
    在任何 TE 计算开始前强制调用。
    """
    print()
    print("=" * 74)
    print("  样本配对报告 (Sample Pairing Report)")
    print("=" * 74)
    print(f"  {'生物学样本':<16s} {'Ribo SRR':<16s} {'RNA SRR':<16s} "
          f"{'Cell Line':<16s} {'Organism'}")
    print("-" * 74)

    for key in sorted(pairing.keys()):
        p = pairing[key]
        ribo_str = ",".join(p["ribo_srr"])
        rna_str = ",".join(p["rna_srr"])
        print(f"  {p['ribo_gsm']:<16s} {ribo_str:<16s} {rna_str:<16s} "
              f"{p['cell_line']:<16s} {p['organism']}")
        print(f"    [配对成功] Ribo={ribo_str} <---> RNA={rna_str}")

    print("-" * 74)
    print(f"  共 {len(pairing)} 个有效配对")
    print("=" * 74)
    print()


def validate_and_align_columns(
    df_ribo: pd.DataFrame,
    df_rna: pd.DataFrame,
    pairing: dict,
) -> tuple:
    """
    根据配对字典严格对齐 ribo 和 rna 的列。

    - 仅保留在 pairing 中有映射、且在两个 DataFrame 中都存在的列
    - 输出的两个 DataFrame 列顺序严格一致 (第 i 列 ribo 对应第 i 列 rna)
    - 无法匹配的列被 DROP 并记录日志

    返回:
        (aligned_ribo, aligned_rna, aligned_info_rows)
    """
    ribo_cols = set(df_ribo.columns)
    rna_cols = set(df_rna.columns)

    # ── 处理 "all" 列名的特殊情况 ──────────────────────────────────────
    # 当使用 all.ribo (合并后的 .ribo 文件) 时, ribopy 提取的列名为 "all"
    # 此时如果只有单个配对, 直接将 "all" 映射到该配对
    ribo_is_all = (list(df_ribo.columns) == ["all"])
    rna_is_all = (list(df_rna.columns) == ["all"])

    if (ribo_is_all or rna_is_all) and len(pairing) == 1:
        p = list(pairing.values())[0]
        print(f"[对齐] 检测到 'all' 列名 (来自合并的 all.ribo), 单配对模式自动映射")
        if ribo_is_all:
            df_ribo = df_ribo.rename(columns={"all": p["ribo_gsm"]})
            ribo_cols = set(df_ribo.columns)
            print(f"  ribo_raw: 'all' -> '{p['ribo_gsm']}'")
        if rna_is_all:
            df_rna = df_rna.rename(columns={"all": p["rna_gsm"]})
            rna_cols = set(df_rna.columns)
            print(f"  rnaseq_raw: 'all' -> '{p['rna_gsm']}'")

    aligned_ribo_cols = []
    aligned_rna_cols = []
    aligned_info = []  # for infor_filter.csv
    dropped_pairs = []

    for key in sorted(pairing.keys()):
        p = pairing[key]
        # 每个 Ribo 样本可能有多个 SRR (多 run)
        # 在 count 矩阵中，列名可能是 SRR, GSM, 或 experiment_alias
        # 我们需要找到实际存在的列名
        ribo_found = None
        rna_found = None

        # 尝试匹配: SRR accession → GSM alias → SRX accession
        for candidate in p["ribo_srr"] + [p["ribo_gsm"], p["ribo_srx"]]:
            if candidate in ribo_cols:
                ribo_found = candidate
                break

        for candidate in p["rna_srr"] + [p["rna_gsm"], p["rna_srx"]]:
            if candidate in rna_cols:
                rna_found = candidate
                break

        if ribo_found and rna_found:
            aligned_ribo_cols.append(ribo_found)
            aligned_rna_cols.append(rna_found)
            aligned_info.append({
                "experiment_alias": ribo_found,
                "cell_line": p["cell_line"] or p["ribo_gsm"],
                "ribo_gsm": p["ribo_gsm"],
                "rna_gsm": p["rna_gsm"],
            })
        else:
            reason = ""
            if not ribo_found:
                reason += (f"Ribo 列 {p['ribo_srr']+[p['ribo_gsm']]} "
                           f"在 ribo_raw.csv 中未找到; ")
            if not rna_found:
                reason += (f"RNA 列 {p['rna_srr']+[p['rna_gsm']]} "
                           f"在 rnaseq_raw.csv 中未找到; ")
            dropped_pairs.append((p["ribo_gsm"], reason))

    if dropped_pairs:
        print(f"[对齐] 警告: {len(dropped_pairs)} 个配对因列名不匹配被 DROP:")
        for gsm, reason in dropped_pairs[:10]:
            print(f"  [DROP] {gsm}: {reason}")

    if not aligned_ribo_cols:
        raise RuntimeError(
            "严重错误: 没有任何 Ribo↔RNA 配对能在 count 矩阵中找到对应列!\n"
            f"  ribo_raw.csv 列名 (前10): {list(df_ribo.columns[:10])}\n"
            f"  rnaseq_raw.csv 列名 (前10): {list(df_rna.columns[:10])}\n"
            f"  配对字典中的 SRR 示例: "
            f"{list(pairing.values())[0] if pairing else '(空)'}")

    # 确认对齐后提取子集
    df_ribo_aligned = df_ribo[aligned_ribo_cols].copy()
    df_rna_aligned = df_rna[aligned_rna_cols].copy()

    # 统一列名为 ribo 侧列名 (保证两个 DF 列名完全一致)
    rename_map = dict(zip(aligned_rna_cols, aligned_ribo_cols))
    df_rna_aligned = df_rna_aligned.rename(columns=rename_map)

    # 打印对齐摘要
    orphan_ribo = ribo_cols - set(aligned_ribo_cols)
    orphan_rna = rna_cols - set(aligned_rna_cols)
    print(f"[对齐] 成功对齐 {len(aligned_ribo_cols)} 对样本列")
    if orphan_ribo:
        print(f"[对齐] ribo_raw.csv 中未配对列 ({len(orphan_ribo)}): "
              f"{sorted(orphan_ribo)[:5]}{'...' if len(orphan_ribo)>5 else ''}")
    if orphan_rna:
        print(f"[对齐] rnaseq_raw.csv 中未配对列 ({len(orphan_rna)}): "
              f"{sorted(orphan_rna)[:5]}{'...' if len(orphan_rna)>5 else ''}")

    return df_ribo_aligned, df_rna_aligned, aligned_info


# =============================================================================
# Stage 1: 数据预处理 — CPM 标准化、低表达基因过滤、dummy gene 合并
# =============================================================================
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  此段计算逻辑与过滤阈值 100% 严格复刻自 CenikLab 原文代码             │
# │  (src/ribobase_counts_processing.py)。                                 │
# │                                                                         │
# │  核心参数:                                                              │
# │    cpm_cutoff     = 1   (CPM 阈值，低于此值视为低表达)                  │
# │    overall_cutoff = 70  (百分比: 低表达样本占比超过此值的基因为 dummy)   │
# │                                                                         │
# │  过滤公式 (原文 ribobase_counts_processing.py:63-68):                   │
# │    row_cut_off = int(overall_cutoff / 100 * n_samples)                  │
# │    dummy_genes = genes where (CPM < cpm_cutoff).sum() > row_cut_off     │
# │                                                                         │
# │  Paired 模式: ribo dummy ∪ rna dummy (outer merge)                     │
# │  Dummy gene 合并: 所有 dummy 基因的 count 求和为一行 "dummy_gene"       │
# │  非 polyA 基因移除: 使用 nonpolyA_gene.csv (如提供)                    │
# │                                                                         │
# │  标准化策略: CPM (Counts Per Million) via bioinfokit.analys.norm().cpm() │
# └─────────────────────────────────────────────────────────────────────────┘
# =============================================================================

def CPM_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    CPM 标准化。
    100% 复刻 CenikLab ribobase_counts_processing.py:45-52
    """
    if HAS_BIOINFOKIT:
        nm = bioinfokit_norm()
        nm.cpm(df=df)
        return nm.cpm_norm
    else:
        # fallback: 手动 CPM (数学等价)
        return df.div(df.sum(axis=0), axis=1) * 1e6


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    分位数标准化。
    100% 复刻 CenikLab ribobase_counts_processing.py:30-43
    """
    df_sorted = pd.DataFrame(
        np.sort(df.values, axis=0),
        index=df.index,
        columns=df.columns,
    )
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn = df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return df_qn


def data_process(df_path_or_df):
    """
    数据处理: 读取 → groupby(index).mean() → strip 'dedup' → CPM → Quantile
    100% 复刻 CenikLab ribobase_counts_processing.py:54-60
    """
    if isinstance(df_path_or_df, (str, Path)):
        df_all = pd.read_csv(df_path_or_df, index_col=0)
    else:
        df_all = df_path_or_df.copy()

    df_count = df_all.groupby(df_all.index).mean()
    df_count.columns = df_count.columns.str.rstrip("dedup")
    df_cpm = CPM_normalize(df_count)
    df_quantile = quantile_normalize(df_count)
    return df_count, df_cpm, df_quantile


def dummy_gene_df(df_cpm: pd.DataFrame,
                  cpm_cutoff: int = 1,
                  overall_cutoff: int = 70) -> pd.Series:
    """
    识别低表达 dummy 基因。
    100% 复刻 CenikLab ribobase_counts_processing.py:63-69

    原文逻辑:
      row_cut_off = int(overall_cutoff / 100 * len(df.columns))
      dummy = genes where count(CPM < cpm_cutoff) > row_cut_off
    """
    row_cut_off = int(overall_cutoff / 100 * len(df_cpm.columns))
    df_dummy = df_cpm[(df_cpm < cpm_cutoff).sum(axis='columns') > row_cut_off]
    dummy_gene = df_dummy.index.to_series().rename("dummy_gene_series")
    return dummy_gene


def combine_dummy_gene(dummy_gene: pd.Series,
                       df: pd.DataFrame) -> pd.DataFrame:
    """
    将过滤掉的基因合并为单行 'dummy_gene'。
    100% 复刻 CenikLab ribobase_counts_processing.py:71-79
    """
    non_dummy = df[~df.index.isin(dummy_gene.index)]
    dummy_df = df[df.index.isin(dummy_gene.index)]
    dummy_result = pd.DataFrame(dummy_df.sum())
    dummy_result.rename(columns={0: 'dummy_gene'}, inplace=True)
    frames = [non_dummy, dummy_result.T]
    df_count_dummy = pd.concat(frames)
    return df_count_dummy


def stage1_preprocess(output_dir: str,
                      cpm_cutoff: int = DEFAULT_CPM_CUTOFF,
                      overall_cutoff: int = DEFAULT_OVERALL_CUTOFF,
                      nonpolya_csv: str = None,
                      pairing: dict = None) -> None:
    """
    Stage 1: paired 模式预处理。
    100% 复刻 CenikLab ribobase_counts_processing.py 的 'paired' 分支
    (lines 105-135)

    如果 pairing 不为 None，则在预处理前执行强约束配对对齐。
    """
    ribo_raw_path = os.path.join(output_dir, "ribo_raw.csv")
    rna_raw_path = os.path.join(output_dir, "rnaseq_raw.csv")

    if not os.path.exists(ribo_raw_path):
        raise FileNotFoundError(f"未找到 {ribo_raw_path}")
    if not os.path.exists(rna_raw_path):
        raise FileNotFoundError(
            f"未找到 {rna_raw_path}。Paired 模式需要 RNA-seq 数据。")

    # ── Step 0.5: 配对验证与对齐 (如果提供了 pairing 字典) ──────────────
    if pairing:
        print("[Stage 0.5] 执行强约束配对对齐...")
        df_ribo_raw = pd.read_csv(ribo_raw_path, index_col=0)
        df_rna_raw = pd.read_csv(rna_raw_path, index_col=0)

        print_pairing_report(pairing)

        df_ribo_aligned, df_rna_aligned, aligned_info = \
            validate_and_align_columns(df_ribo_raw, df_rna_raw, pairing)

        # 覆写对齐后的 CSV (后续 data_process 将读取这些文件)
        df_ribo_aligned.to_csv(ribo_raw_path)
        df_rna_aligned.to_csv(rna_raw_path)
        print(f"[Stage 0.5] 已覆写对齐后的 count 矩阵 "
              f"({len(df_ribo_aligned.columns)} 对配对样本)")

        # 覆写 infor_filter.csv
        info_list = []
        for i, info in enumerate(aligned_info, 1):
            info_list.append({
                "": i,
                "experiment_alias": info["experiment_alias"],
                "cell_line": info["cell_line"],
            })
        pd.DataFrame(info_list).to_csv(
            os.path.join(output_dir, "infor_filter.csv"), index=False)
    else:
        # ── 无 pairing 字典时的安全检查 ─────────────────────────────────
        df_ribo_check = pd.read_csv(ribo_raw_path, index_col=0, nrows=0)
        df_rna_check = pd.read_csv(rna_raw_path, index_col=0, nrows=0)
        ribo_cols = set(df_ribo_check.columns)
        rna_cols = set(df_rna_check.columns)
        if ribo_cols == rna_cols:
            print("[Stage 0.5] 列名完全一致 — 假定来自 .ribo 文件的内置配对 (安全)")
        else:
            common = ribo_cols & rna_cols
            only_ribo = ribo_cols - rna_cols
            only_rna = rna_cols - ribo_cols
            print("[Stage 0.5] ⚠ 警告: ribo_raw.csv 与 rnaseq_raw.csv 列名不一致!")
            print(f"  共同列: {len(common)}, 仅 Ribo: {len(only_ribo)}, 仅 RNA: {len(only_rna)}")
            if only_ribo:
                print(f"  仅 Ribo 列(前5): {sorted(only_ribo)[:5]}")
            if only_rna:
                print(f"  仅 RNA 列(前5): {sorted(only_rna)[:5]}")
            print("  提示: 建议提供 --metadata_csv 和 --mapping_csv 以启用强约束配对")
            print("  当前将按共同列名进行 INNER 匹配 (可能丢失数据)")
            if common:
                common_sorted = sorted(common)
                df_ribo_tmp = pd.read_csv(ribo_raw_path, index_col=0)[common_sorted]
                df_rna_tmp = pd.read_csv(rna_raw_path, index_col=0)[common_sorted]
                df_ribo_tmp.to_csv(ribo_raw_path)
                df_rna_tmp.to_csv(rna_raw_path)
                print(f"  已按 {len(common)} 个共同列对齐")
            else:
                raise RuntimeError(
                    "严重错误: ribo_raw.csv 和 rnaseq_raw.csv 没有任何共同列名!\n"
                    "  请提供 --metadata_csv 和 --mapping_csv 以建立配对映射。")

    # ── Step 1a: data_process ────────────────────────────────────────────
    print("[Stage 1] 预处理 Ribo-seq 原始计数...")
    ribo_count_temp, ribo_CPM_temp, ribo_Q_temp = data_process(ribo_raw_path)
    print("[Stage 1] 预处理 RNA-seq 原始计数...")
    rna_count_temp, rna_CPM_temp, rna_Q_temp = data_process(rna_raw_path)

    # ── Step 1b: identify dummy genes ────────────────────────────────────
    # 100% 复刻: 分别对 ribo 和 rna CPM 计算 dummy，然后 OUTER merge (union)
    print(f"[Stage 1] 识别 dummy 基因 (cpm_cutoff={cpm_cutoff}, overall_cutoff={overall_cutoff})...")
    ribo_dummy_gene = dummy_gene_df(ribo_CPM_temp, cpm_cutoff, overall_cutoff)
    rna_dummy_gene = dummy_gene_df(rna_CPM_temp, cpm_cutoff, overall_cutoff)
    all_dummy = pd.merge(
        ribo_dummy_gene, rna_dummy_gene,
        how="outer", left_index=True, right_index=True,
    )
    print(f"  Ribo dummy: {len(ribo_dummy_gene)}, RNA dummy: {len(rna_dummy_gene)}, "
          f"Union: {len(all_dummy)}")

    # ── Step 1c: combine dummy gene ──────────────────────────────────────
    print("[Stage 1] 合并 dummy 基因 (pre-final)...")
    ribo_count_dummy_wployA = combine_dummy_gene(all_dummy, ribo_count_temp)
    ribo_CPM_wployA = ribo_CPM_temp[~ribo_CPM_temp.index.isin(all_dummy.index)]
    ribo_Q_wployA = ribo_Q_temp[~ribo_Q_temp.index.isin(all_dummy.index)]
    rna_count_dummy_wployA = combine_dummy_gene(all_dummy, rna_count_temp)
    rna_CPM_wployA = rna_CPM_temp[~rna_CPM_temp.index.isin(all_dummy.index)]
    rna_Q_wployA = rna_Q_temp[~rna_Q_temp.index.isin(all_dummy.index)]

    # ── Step 1d: remove non-polyA genes ──────────────────────────────────
    # 100% 复刻 CenikLab ribobase_counts_processing.py:122-128
    if nonpolya_csv and os.path.exists(nonpolya_csv):
        print(f"[Stage 1] 移除非 polyA 基因 ({nonpolya_csv})...")
        polyA = pd.read_csv(nonpolya_csv, index_col=0)
        ribo_count_dummy = ribo_count_dummy_wployA[
            ~ribo_count_dummy_wployA.index.isin(polyA.index)]
        ribo_CPM = ribo_CPM_wployA[~ribo_CPM_wployA.index.isin(polyA.index)]
        ribo_Q = ribo_Q_wployA[~ribo_Q_wployA.index.isin(polyA.index)]
        rna_count_dummy = rna_count_dummy_wployA[
            ~rna_count_dummy_wployA.index.isin(polyA.index)]
        rna_CPM = rna_CPM_wployA[~rna_CPM_wployA.index.isin(polyA.index)]
        rna_Q = rna_Q_wployA[~rna_Q_wployA.index.isin(polyA.index)]
        print(f"  移除 {len(polyA)} 个非 polyA 基因")
    else:
        if nonpolya_csv:
            print(f"[Stage 1] 警告: nonpolyA 文件不存在 ({nonpolya_csv})，跳过该步骤")
        ribo_count_dummy = ribo_count_dummy_wployA
        ribo_CPM = ribo_CPM_wployA
        ribo_Q = ribo_Q_wployA
        rna_count_dummy = rna_count_dummy_wployA
        rna_CPM = rna_CPM_wployA
        rna_Q = rna_Q_wployA

    # ── Step 1e: save ────────────────────────────────────────────────────
    print("[Stage 1] 保存预处理结果...")
    ribo_count_dummy.to_csv(os.path.join(output_dir, "ribo_paired_count_dummy.csv"))
    ribo_CPM.to_csv(os.path.join(
        output_dir, f"ribo_paired_cpm_dummy_{overall_cutoff}.csv"))
    ribo_Q.to_csv(os.path.join(
        output_dir, f"ribo_paired_quantile_dummy_{overall_cutoff}.csv"))
    rna_count_dummy.to_csv(os.path.join(output_dir, "rna_paired_count_dummy.csv"))
    rna_CPM.to_csv(os.path.join(
        output_dir, f"rna_paired_cpm_dummy_{overall_cutoff}.csv"))
    rna_Q.to_csv(os.path.join(
        output_dir, f"rna_paired_quantile_dummy_{overall_cutoff}.csv"))

    print(f"  [OK] ribo_paired_count_dummy.csv  ({ribo_count_dummy.shape[0]} genes)")
    print(f"  [OK] rna_paired_count_dummy.csv   ({rna_count_dummy.shape[0]} genes)")


# =============================================================================
# Stage 2: CLR/ILR 成分回归 — 调用 TE.R
# =============================================================================
# ┌─────────────────────────────────────────────────────────────────────────┐
# │  此段调用的 TE.R 脚本 100% 严格复刻自 CenikLab 原文代码 (src/TE.R)。 │
# │                                                                         │
# │  TE 数学公式:                                                           │
# │    1. CLR 变换: propr(data, metric="rho", ivar="clr", alpha=NA, p=100) │
# │    2. CLR → ILR: compositions::clr2ilr()                               │
# │    3. 成分线性回归: lm(RIBO_ilr[,i] ~ RNA_ilr[,i]) 对每个基因          │
# │    4. TE = ilr2clr(residuals) — 回归残差转换回 CLR 空间                │
# │                                                                         │
# │  伪计数处理: propr 内部对零值执行乘性简单替换 (multiplicative simple    │
# │  replacement)，无需外部手动添加伪计数。                                 │
# └─────────────────────────────────────────────────────────────────────────┘
# =============================================================================

def stage2_run_te_r(output_dir: str) -> bool:
    """
    Stage 2: 调用 TE.R 执行 CLR/ILR 成分回归。
    R 脚本读取 {output_dir}/ribo_paired_count_dummy.csv 和
    rna_paired_count_dummy.csv，输出 human_TE_sample_level.rda
    """
    te_r = str(TE_R_PATH)
    if not os.path.exists(te_r):
        print(f"[Stage 2] 错误: 未找到 {te_r}")
        return False

    dummy_ribo = os.path.join(output_dir, "ribo_paired_count_dummy.csv")
    dummy_rna = os.path.join(output_dir, "rna_paired_count_dummy.csv")
    if not os.path.exists(dummy_ribo) or not os.path.exists(dummy_rna):
        print("[Stage 2] 错误: 缺少 Stage 1 输出的 count_dummy 文件")
        return False

    print(f"[Stage 2] 调用 TE.R (CLR/ILR 成分回归)...")
    print(f"  R 脚本: {te_r}")
    print(f"  工作目录: {output_dir}")

    cmd = ["Rscript", te_r, os.path.abspath(output_dir)]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200,
        )
        if result.stdout:
            for line in result.stdout.strip().split("\n"):
                print(f"  [R] {line}")
        if result.returncode != 0:
            print(f"[Stage 2] 错误: TE.R 返回码 {result.returncode}")
            if result.stderr:
                for line in result.stderr.strip().split("\n"):
                    print(f"  [R stderr] {line}")
            return False
        print("[Stage 2] TE.R 执行成功")
        return True
    except FileNotFoundError:
        print("[Stage 2] 错误: 未找到 Rscript，请确保 R 已安装")
        return False
    except subprocess.TimeoutExpired:
        print("[Stage 2] 错误: TE.R 执行超时 (>2h)")
        return False


# =============================================================================
# Stage 3: 后处理 — .rda → CSV, 按条件分组, 最终汇总
# =============================================================================

def stage3_postprocess(output_dir: str, info_csv: str = None) -> str:
    """
    Stage 3: 后处理 TE 结果。
    - 从 .rda 提取为 CSV
    - 按 cell_line 分组取均值 (复刻 TE.R:59-67)
    - 转置输出 (复刻 transpose_TE.py)
    - 生成 te_results_final.csv

    返回最终结果文件路径。
    """
    rda_file = os.path.join(output_dir, "human_TE_sample_level.rda")
    final_out = os.path.join(output_dir, "te_results_final.csv")

    if not os.path.exists(rda_file):
        print("[Stage 3] 警告: 未找到 human_TE_sample_level.rda")
        # 尝试直接读取已有的 CSV
        cellline_csv = os.path.join(output_dir, "human_TE_cellline_all.csv")
        if os.path.exists(cellline_csv):
            shutil.copy(cellline_csv, final_out)
            return final_out
        return None

    # 3a: .rda → sample-level CSV
    print("[Stage 3] 转换 .rda → CSV...")
    sample_csv = os.path.join(output_dir, "TE_sample_level.csv")
    r_cmd = (
        f'load("{rda_file}"); '
        f'write.csv(t(human_TE), "{sample_csv}")'
    )
    try:
        subprocess.run(
            ["Rscript", "-e", r_cmd],
            capture_output=True, check=True, timeout=300,
        )
        print(f"  [OK] {sample_csv}")
    except Exception as e:
        print(f"[Stage 3] 警告: .rda 转换失败: {e}")
        return None

    # 3b: 按 cell_line 分组 (复刻 TE.R:59-67)
    if info_csv and os.path.exists(info_csv):
        print("[Stage 3] 按 cell_line 分组取均值...")
        te_sample = pd.read_csv(sample_csv, index_col=0)
        infor = pd.read_csv(info_csv)
        te_sample_t = te_sample.T if te_sample.shape[0] < te_sample.shape[1] else te_sample

        if "experiment_alias" in infor.columns and "cell_line" in infor.columns:
            mapping = dict(zip(infor["experiment_alias"], infor["cell_line"]))
            te_sample_t.index = te_sample_t.index.map(
                lambda x: mapping.get(x, x))
            te_by_cellline = te_sample_t.groupby(te_sample_t.index).mean()
            cellline_csv = os.path.join(output_dir, "human_TE_cellline_all.csv")
            te_by_cellline.to_csv(cellline_csv)
            print(f"  [OK] {cellline_csv}")

    # 3c: 转置最终输出 (复刻 transpose_TE.py:9-12)
    print("[Stage 3] 生成最终结果...")
    if os.path.exists(os.path.join(output_dir, "human_TE_cellline_all.csv")):
        df_final = pd.read_csv(
            os.path.join(output_dir, "human_TE_cellline_all.csv"), index_col=0)
        transposed = df_final.T
        transposed.index = transposed.index.str.replace(r"\.(.*)", "", regex=True)
        transposed.to_csv(final_out)
    elif os.path.exists(sample_csv):
        # 无分组信息时直接输出 sample-level
        shutil.copy(sample_csv, final_out)
    else:
        return None

    print(f"  [OK] {final_out}")
    return final_out


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "翻译效率 (TE) 计算 — CenikLab CLR/ILR 成分回归方法\n"
            "核心算法 100% 复刻自 CenikLab TE_model 原文代码。"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--ribo_dir", "-i",
        type=str, default=None,
        help="包含 .ribo 文件的输入目录",
    )
    parser.add_argument(
        "--output_dir", "-o",
        type=str, default=str(DEFAULT_PROCESSED_DIR),
        help="输出目录 (默认: data/processed)",
    )
    parser.add_argument(
        "--cpm_cutoff", "-c",
        type=int, default=DEFAULT_CPM_CUTOFF,
        help=f"CPM 过滤阈值 (默认: {DEFAULT_CPM_CUTOFF})",
    )
    parser.add_argument(
        "--overall_cutoff", "-a",
        type=int, default=DEFAULT_OVERALL_CUTOFF,
        help=f"样本百分比阈值 (默认: {DEFAULT_OVERALL_CUTOFF})",
    )
    parser.add_argument(
        "--nonpolya_csv",
        type=str, default=None,
        help="非 polyA 基因列表 CSV 路径 (可选)",
    )
    parser.add_argument(
        "--info_csv",
        type=str, default=None,
        help="样本信息 CSV (含 experiment_alias, cell_line 列)",
    )
    parser.add_argument(
        "--metadata_csv",
        type=str, default=None,
        help="元数据 CSV (含 experiment_alias, matched_RNA-seq 等列, 用于强约束配对)",
    )
    parser.add_argument(
        "--mapping_csv",
        type=str, default=None,
        help="SRX→SRR 映射 CSV (含 Run, srx 列, 配合 --metadata_csv 使用)",
    )
    parser.add_argument(
        "--skip_extract",
        action="store_true",
        help="跳过 Stage 0 (假定 ribo_raw.csv 已存在)",
    )
    parser.add_argument(
        "--skip_r",
        action="store_true",
        help="跳过 Stage 2 R 脚本 (仅做预处理)",
    )
    parser.add_argument(
        "--dummy",
        action="store_true",
        help="生成含假 TE 值的测试 CSV (管线测试用)",
    )
    return parser.parse_args()


def generate_dummy_te(output_path: str) -> None:
    """生成假 TE 数据用于管线测试。"""
    np.random.seed(42)
    n_genes = 200
    n_samples = 6
    gene_ids = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
    sample_ids = [f"sample_{i}" for i in range(1, n_samples + 1)]

    te_values = np.random.randn(n_genes, n_samples) * 0.5
    df = pd.DataFrame(te_values, index=gene_ids, columns=sample_ids)
    df.index.name = "gene"

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.to_csv(output_path)
    print(f"[te_calculator] 已生成假 TE 数据 ({n_genes} genes × {n_samples} samples) → {output_path}")


def main() -> None:
    args = parse_args()
    output_dir = os.path.abspath(args.output_dir)

    # -- Dummy 模式 -------------------------------------------------------
    if args.dummy:
        final_out = os.path.join(output_dir, "te_results_final.csv")
        generate_dummy_te(final_out)
        return

    print("=" * 70)
    print("  TE Calculator — CenikLab CLR/ILR 成分回归方法")
    print("  核心算法 100% 严格复刻自 CenikLab TE_model 原文代码")
    print("=" * 70)
    print(f"  输出目录:      {output_dir}")
    print(f"  CPM 阈值:      {args.cpm_cutoff}")
    print(f"  Overall 阈值:  {args.overall_cutoff}%")
    print(f"  Skip extract:  {args.skip_extract}")
    print(f"  Skip R:        {args.skip_r}")
    print("-" * 70)

    os.makedirs(output_dir, exist_ok=True)

    # -- Stage 0: Extract -------------------------------------------------
    if not args.skip_extract:
        if not args.ribo_dir:
            print("错误: 需要 --ribo_dir 或 --skip_extract", file=sys.stderr)
            sys.exit(1)
        stage0_extract(args.ribo_dir, output_dir)
    else:
        print("[Stage 0] 跳过 (使用已有 ribo_raw.csv / rnaseq_raw.csv)")

    # -- Stage 0.5: Build pairing (if metadata provided) -------------------
    pairing = None
    if args.metadata_csv and args.mapping_csv:
        if not os.path.exists(args.metadata_csv):
            print(f"错误: --metadata_csv 文件不存在: {args.metadata_csv}",
                  file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.mapping_csv):
            print(f"错误: --mapping_csv 文件不存在: {args.mapping_csv}",
                  file=sys.stderr)
            sys.exit(1)
        pairing = build_sample_pairing(args.metadata_csv, args.mapping_csv)
    elif args.metadata_csv or args.mapping_csv:
        print("警告: --metadata_csv 和 --mapping_csv 必须同时提供，"
              "忽略配对验证", file=sys.stderr)

    # -- Stage 1: Preprocess ----------------------------------------------
    stage1_preprocess(
        output_dir,
        cpm_cutoff=args.cpm_cutoff,
        overall_cutoff=args.overall_cutoff,
        nonpolya_csv=args.nonpolya_csv,
        pairing=pairing,
    )

    # -- Stage 2: TE.R ----------------------------------------------------
    if not args.skip_r:
        success = stage2_run_te_r(output_dir)
        if not success:
            print("[警告] Stage 2 未成功，跳过 Stage 3")
            print("  预处理文件已保存，可手动运行:")
            print(f"  Rscript {TE_R_PATH} {output_dir}")
            sys.exit(1)
    else:
        print("[Stage 2] 跳过 R 脚本执行")

    # -- Stage 3: Postprocess ---------------------------------------------
    if not args.skip_r:
        info_csv = args.info_csv or os.path.join(output_dir, "infor_filter.csv")
        final_path = stage3_postprocess(output_dir, info_csv)
        if final_path:
            print(f"\n{'=' * 70}")
            print(f"  TE 计算完成 → {final_path}")
            print(f"{'=' * 70}")
        else:
            print("[警告] Stage 3 后处理未生成最终文件")
    else:
        print(f"\n{'=' * 70}")
        print(f"  预处理完成 (已跳过 R 脚本)")
        print(f"  中间文件位于: {output_dir}")
        print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
