import argparse
import os
import pandas as pd
import numpy as np
import yaml
from pprint import pprint

###############################################################################
#               GENERATE YAML FOR RIBOFLOW  (CSV-Based, CCDS 版)
###############################################################################
# 替代原 SQLite 驱动的 generate_yaml.py。
# 直接从本地 CSV 元数据表构建 RiboFlow 所需的 per-study project.yaml。
#
# 数据源:
#   1. srx_to_srr_mapping.csv  — SRX/SRR/GSM/organism/corrected_type/study_name
#   2. metadata.csv (header=1) — experiment_alias(GSM)/organism/matched_RNA-seq
#                                  /threep_adapter/fivep_adapter/corrected_type
#
# Sample Run:
#   python generate_yaml.py --study GSE100007 --template ../project.yaml \
#       --output ../../data/interim/project --mapping ../../data/external/srx_to_srr_mapping.csv \
#       --metadata ../../data/external/metadata/metadata.csv
###############################################################################

# ── 物种名归一化映射 ────────────────────────────────────────────────────────
# metadata.csv 中物种名不一致 (如 "Human", "Homo sapiens", "Mouse", "D. melanogaster")
# 需要统一映射到 references.yaml 的键 (全小写拉丁名)
ORGANISM_ALIAS = {
    # ── 英文俗名 / 缩写 → 拉丁名 ──────────────────────────────────────────
    "human":                     "homo sapiens",
    "mouse":                     "mus musculus",
    "rat":                       "rattus norvegicus",
    "arabidopsis":               "arabidopsis thaliana",
    "soybean":                   "zea mays",             # metadata 中 "Soybean" 暂映射到 zea_mays (需确认)
    "celegans":                  "caenorhabditis elegans",
    "d. melanogaster":           "drosophila melanogaster",
    # ── metadata.csv 中的 typo / 变体 → 正确拉丁名 ─────────────────────────
    "caenorhaboditis elegans":   "caenorhabditis elegans",
    "caenorhabitis brenneri":    "caenorhabditis brenneri",
    # ── 已在 references.yaml 中定义的全部拉丁名 (identity 映射, 确保不遗漏) ──
    "homo sapiens":              "homo sapiens",
    "mus musculus":              "mus musculus",
    "rattus norvegicus":         "rattus norvegicus",
    "arabidopsis thaliana":      "arabidopsis thaliana",
    "caenorhabditis elegans":    "caenorhabditis elegans",
    "caenorhabditis brenneri":   "caenorhabditis brenneri",
    "caenorhabditis briggsae":   "caenorhabditis briggsae",
    "caenorhabditis remanei":    "caenorhabditis remanei",
    "danio rerio":               "danio rerio",
    "drosophila melanogaster":   "drosophila melanogaster",
    "saccharomyces cerevisiae":  "saccharomyces cerevisiae",
    "saccharomyces paradoxus":   "saccharomyces paradoxus",
    "saccharomyces uvarum":      "saccharomyces uvarum",
    "schizosaccharomyces pombe": "schizosaccharomyces pombe",
    "escherichia coli":          "escherichia coli",
    "xenopus laevis":            "xenopus laevis",
    "zea mays":                  "zea mays",
    "trypanosoma brucei":        "trypanosoma brucei",
    "plasmodium falciparum":     "plasmodium falciparum",
    "leishmania donovani":       "leishmania donovani",
    "candida albicans":          "candida albicans",
    "staphylococcus aureus":     "staphylococcus aureus",
    "streptomyces clavuligerus": "streptomyces clavuligerus",
    # ── 业务目标物种 ──────────────────────────────────────────────────────
    "bombyx mori":               "bombyx mori",
    "spider_species":            "spider_species",
}


def normalize_organism(raw_name):
    """将 metadata 中不一致的物种名归一化为 references.yaml 键格式 (全小写拉丁名)。"""
    if pd.isna(raw_name) or raw_name.strip() == "":
        return None
    lowered = raw_name.strip().lower()
    # 处理复合条目如 "Saccharomyces cerevisiae* Saccharomyces paradoxus"
    if "*" in lowered:
        lowered = lowered.split("*")[0].strip()
    return ORGANISM_ALIAS.get(lowered, lowered)


# ── Adapter clip 序列生成 ────────────────────────────────────────────────────
# 默认 fallback adapter (Illumina universal)
FALLBACK_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def generate_clip_sequence(clipping_param_base, adapter_set, experiment_type):
    """从 adapter 集合构造 cutadapt 参数字符串。"""
    # 过滤空值
    adapter_set = {a for a in adapter_set if a and str(a) != "nan" and a.strip() != ""}

    if len(adapter_set) > 1:
        print("WARNING - Multiple adapters in {} experiments: {}".format(
            experiment_type, adapter_set))

    clip_sequence = clipping_param_base

    if len(adapter_set) >= 1:
        candidate_adapter = sorted(adapter_set)[0]  # 取字典序第一个，确保确定性
        overlap = 4
        # 计算 adapter 5' 端 N 数量以调整 overlap
        for ch in candidate_adapter:
            if ch == 'N':
                overlap += 1
            else:
                break

        if len(candidate_adapter) > 0:
            clip_sequence += " -a " + candidate_adapter + " --overlap=" + str(overlap)
            if experiment_type == 'Ribo-Seq':
                clip_sequence += " --trimmed-only"

    return clip_sequence


# ── 主函数: 从 CSV 生成 per-study YAML ─────────────────────────────────────

def generate_yaml(study, template, output, download_path,
                  mapping_csv, metadata_csv,
                  reference_file   = "scripts/references.yaml",
                  reference_folder = "reference"):
    """
    从本地 CSV 生成 RiboFlow project.yaml，取代原 SQLite 版本。

    Parameters
    ----------
    study : str
        Study 名称, 如 "GSE100007" 或 "GSE100007_dedup"
    template : str
        project.yaml 模板路径
    output : str
        输出目录 (study YAML 会写入 output/{gse_only}/{study}.yaml)
    download_path : str
        FASTQ 文件根目录 (CCDS: ../../data/raw/fastq)
    mapping_csv : str
        srx_to_srr_mapping.csv 路径
    metadata_csv : str
        metadata.csv 路径 (header=1)
    reference_file : str
        references.yaml 路径
    reference_folder : str
        reference 索引根目录
    """

    # ── 1. 解析 study 名称 ───────────────────────────────────────────────
    dedup_val = False
    gse_only  = study

    study_contents = study.split("_")
    if len(study_contents) > 1:
        gse_only = study_contents[0]
        if study_contents[1] == "dedup":
            dedup_val = True
        elif study_contents[1] != "test":
            raise RuntimeError(
                "Invalid run type '{}' for study '{}'".format(study_contents[1], study_contents[0]))

    # ── 2. 加载 CSV 数据 ────────────────────────────────────────────────
    mapping_df  = pd.read_csv(mapping_csv)
    metadata_df = pd.read_csv(metadata_csv, header=1)

    # 过滤 mapping 到当前 study
    study_mapping = mapping_df[mapping_df['study_name'] == gse_only].copy()
    if study_mapping.empty:
        raise Exception("Study '{}' not found in srx_to_srr_mapping.csv".format(gse_only))

    # 区分 Ribo-Seq / RNA-Seq
    ribo_mapping = study_mapping[study_mapping['corrected_type'] == 'Ribo-Seq'].copy()
    rna_mapping  = study_mapping[study_mapping['corrected_type'] == 'RNA-Seq'].copy()

    if ribo_mapping.empty:
        raise Exception("No Ribo-Seq experiments found for study '{}'".format(gse_only))

    # ── 3. 合并 metadata 获取 adapter / organism / matched RNA-Seq ───────
    # metadata 的 join key: experiment_alias (GSM) ↔ mapping 的 gsm
    meta_cols = ['experiment_alias', 'organism', 'matched_RNA-seq_experiment_alias',
                 'corrected_type', 'threep_adapter', 'fivep_adapter',
                 'experiment_accession']
    # 只保留 metadata 中存在的列
    meta_cols = [c for c in meta_cols if c in metadata_df.columns]
    meta_sub  = metadata_df[meta_cols].drop_duplicates(subset=['experiment_alias'])

    ribo_mapping = ribo_mapping.merge(
        meta_sub, left_on='gsm', right_on='experiment_alias', how='left', suffixes=('', '_meta'))
    rna_mapping = rna_mapping.merge(
        meta_sub, left_on='gsm', right_on='experiment_alias', how='left', suffixes=('', '_meta'))

    # ── 4. 物种检测 ─────────────────────────────────────────────────────
    organism_raw = ribo_mapping['organism'].dropna().unique()
    if len(organism_raw) == 0:
        # fallback: 尝试从 mapping 自带的 organism 列取
        organism_raw = ribo_mapping['organism'].dropna().unique()
    if len(organism_raw) == 0:
        raise Exception("No organism detected for study '{}'".format(gse_only))

    normalized_organisms = set(normalize_organism(o) for o in organism_raw if o)
    normalized_organisms.discard(None)

    if len(normalized_organisms) != 1:
        print("WARNING: Multiple organisms detected: {}".format(normalized_organisms))
    cur_organism = sorted(normalized_organisms)[0]

    # ── 5. 加载模板 YAML ────────────────────────────────────────────────
    with open(template) as f:
        ribo_yaml = yaml.load(f, Loader=yaml.FullLoader)

    yaml_name = study
    print("Currently Generating YAML File for: " + yaml_name)

    # ── 6. Adapter 序列处理 ──────────────────────────────────────────────
    ribo_adapters = set()
    if 'threep_adapter' in ribo_mapping.columns:
        ribo_adapters = set(
            ribo_mapping['threep_adapter'].dropna().unique()) - {"", "nan"}

    rna_adapters = set()
    if not rna_mapping.empty and 'threep_adapter' in rna_mapping.columns:
        rna_adapters = set(
            rna_mapping['threep_adapter'].dropna().unique()) - {"", "nan"}

    ribo_clip_base = '-u 1 --maximum-length=40 --minimum-length=15 --quality-cutoff=28'
    rna_clip_base  = '-u 5 -l 40 --quality-cutoff=28'

    # 如果 CSV 没有 adapter 信息, 使用 fallback
    if not ribo_adapters:
        print("  WARNING: No 3' adapter found for Ribo-Seq, using fallback: " + FALLBACK_ADAPTER)
        ribo_adapters = {FALLBACK_ADAPTER}

    ribo_yaml['clip_arguments'] = generate_clip_sequence(ribo_clip_base, ribo_adapters, 'Ribo-Seq')

    if 'rnaseq' in ribo_yaml:
        if not rna_adapters:
            # RNA-Seq adapter fallback: 不加 -a, 仅做质量修剪
            ribo_yaml['rnaseq']['clip_arguments'] = rna_clip_base
        else:
            ribo_yaml['rnaseq']['clip_arguments'] = generate_clip_sequence(
                rna_clip_base, rna_adapters, 'RNA-Seq')

    # ── 7. Reference 路径 ───────────────────────────────────────────────
    with open(reference_file, 'rt') as f:
        yaml_contents = yaml.load(f, Loader=yaml.FullLoader)

    if cur_organism not in yaml_contents:
        supported = sorted(yaml_contents.keys())
        raise Exception(
            "未找到物种 '{}' 的参考基因组配置！\n"
            "请先在 scripts/references.yaml 中添加该物种的索引路径，\n"
            "并将对应的 Bowtie2 索引文件放入 reference/ 目录。\n"
            "当前已支持的物种: {}".format(cur_organism, supported)
        )

    ref_contents = yaml_contents[cur_organism]
    for ref_key in ["filter", "regions", "transcript_lengths", "transcriptome"]:
        ref_path = ref_contents.get(ref_key, "")
        if not ref_path or str(ref_path).strip() == "":
            raise RuntimeError(
                "物种 '{}' 的参考路径 '{}' 为空!\n"
                "请先构建该物种的参考文件并更新 scripts/references.yaml。\n"
                "对于 C. elegans, 可运行 scripts/build_celegans_ref.sh 构建缺失的参考文件。".format(
                    cur_organism, ref_key))
        ribo_yaml['input']['reference'][ref_key] = os.path.join(
            reference_folder, ref_path)

    # ── 8. 构建 FASTQ 路径 (核心: 写死 CCDS 路径) ───────────────────────
    #  GSM → [SRR_1.fastq.gz, SRR_2.fastq.gz, ...]
    #  路径格式: {download_path}/{SRR}_1.fastq.gz  (扁平目录)

    # 获取 matched RNA-Seq GSM → Ribo GSM 的映射
    matched_rna_gsm_map = {}  # ribo_gsm -> rna_gsm
    if 'matched_RNA-seq_experiment_alias' in ribo_mapping.columns:
        for _, row in ribo_mapping.iterrows():
            ribo_gsm    = row['gsm']
            matched_rna = row.get('matched_RNA-seq_experiment_alias')
            if pd.notna(matched_rna) and str(matched_rna).strip() not in ("", "NA", "nan"):
                matched_rna_gsm_map[ribo_gsm] = str(matched_rna).strip()

    # 按 GSM 分组构建 Ribo-Seq FASTQ 路径
    ribo_fastq_dict = {}
    for gsm, group in ribo_mapping.groupby('gsm'):
        ribo_fastq_dict[gsm] = [
            os.path.join(download_path, row['Run'] + '_1.fastq.gz')
            for _, row in group.iterrows()
        ]

    # 按 GSM 分组构建 RNA-Seq FASTQ 路径
    rna_fastq_by_gsm = {}
    if not rna_mapping.empty:
        for gsm, group in rna_mapping.groupby('gsm'):
            rna_fastq_by_gsm[gsm] = [
                os.path.join(download_path, row['Run'] + '_1.fastq.gz')
                for _, row in group.iterrows()
            ]

    # 构建最终 per-ribo-GSM 的 RNA-Seq path dict (只保留已匹配的)
    rna_fastq_dict = {}
    for ribo_gsm in ribo_fastq_dict:
        rna_fastq_dict[ribo_gsm] = []
        matched_rna_gsm = matched_rna_gsm_map.get(ribo_gsm)
        if matched_rna_gsm and matched_rna_gsm in rna_fastq_by_gsm:
            rna_fastq_dict[ribo_gsm] = rna_fastq_by_gsm[matched_rna_gsm]

    # 校验: 至少有一个 Ribo-Seq 实验有 FASTQ 文件
    empty_experiments = [g for g, paths in ribo_fastq_dict.items() if not paths]
    if len(empty_experiments) == len(ribo_fastq_dict):
        raise Exception("All Ribo-Seq GSMs have no sequencing files for study '{}'".format(gse_only))
    if empty_experiments:
        print("  WARNING: {} Ribo-Seq GSMs have no SRR files: {}".format(
            len(empty_experiments), empty_experiments))

    # ── 9. 填充 YAML 结构 ───────────────────────────────────────────────
    ribo_yaml['input']['fastq_base'] = ""
    ribo_yaml['input']['fastq']      = ribo_fastq_dict

    # RNA-Seq 部分
    if 'rnaseq' in ribo_yaml:
        ribo_yaml['rnaseq']['fastq_base'] = ""
        ribo_yaml['rnaseq']['fastq']      = rna_fastq_dict
        ribo_yaml['rnaseq']['deduplicate'] = dedup_val

        # 清理无匹配 RNA-Seq 的条目
        delete_list = [k for k, v in ribo_yaml['rnaseq']['fastq'].items() if not v]
        for k in delete_list:
            del ribo_yaml['rnaseq']['fastq'][k]
        if not ribo_yaml['rnaseq']['fastq']:
            del ribo_yaml['rnaseq']

    ribo_yaml['do_rnaseq']                       = 'rnaseq' in ribo_yaml
    ribo_yaml['output']['output']['base']        = 'output/' + yaml_name
    ribo_yaml['output']['intermediates']['base'] = 'intermediates/' + yaml_name
    ribo_yaml['deduplicate']                     = dedup_val

    # ── 10. 写出 YAML 文件 ──────────────────────────────────────────────
    dir_base   = os.path.join(output, gse_only)
    final_name = os.path.join(dir_base, yaml_name + ".yaml")
    os.makedirs(os.path.dirname(final_name), exist_ok=True)

    with open(final_name, mode='w') as f:
        yaml.dump(ribo_yaml, f)

    print("  -> Written: " + final_name)


# ── CLI 入口 ─────────────────────────────────────────────────────────────

def get_parameters():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
Generate RiboFlow project.yaml from local CSV metadata (CCDS version).

Replaces the original SQLite-based generate_yaml.py.
Reads from srx_to_srr_mapping.csv and metadata.csv instead of db.sqlite3.

Sample Use:
    python generate_yaml.py --study GSE100007 --template ../project.yaml \\
        --output ../../data/interim/project \\
        --mapping ../../data/external/srx_to_srr_mapping.csv \\
        --metadata ../../data/external/metadata/metadata.csv
""")

    parser.add_argument("--template", type=str, required=True,
                        help="File Path of the Template Yaml File")
    parser.add_argument("--output", type=str, required=True,
                        help="Output Directory of the Generated Yaml")
    parser.add_argument("--study", type=str, required=False,
                        help="Study GSE of Interest (e.g., GSE100007 or GSE100007_dedup)")
    parser.add_argument("--text", type=str, required=False,
                        help="Text file containing list of studies (one per line)")
    parser.add_argument("--mapping", type=str, required=True,
                        help="Path to srx_to_srr_mapping.csv")
    parser.add_argument("--metadata", type=str, required=True,
                        help="Path to metadata.csv (header at row 1)")
    parser.add_argument("--download_path", type=str, required=False,
                        default='../../data/raw/fastq',
                        help="FASTQ root directory (default: ../../data/raw/fastq)")
    parser.add_argument("--reference_folder", type=str, required=False,
                        default='reference',
                        help="Folder containing reference indices")
    parser.add_argument("--reference_file", type=str, required=False,
                        default='scripts/references.yaml',
                        help="Yaml file containing reference paths")

    return parser.parse_args()


def main():
    params = get_parameters()

    if bool(params.study) == bool(params.text):
        print("Please specify one of either --study or --text.")
        exit(1)

    common_kwargs = dict(
        output           = params.output,
        template         = params.template,
        download_path    = params.download_path,
        mapping_csv      = params.mapping,
        metadata_csv     = params.metadata,
        reference_file   = params.reference_file,
        reference_folder = params.reference_folder,
    )

    if params.study:
        generate_yaml(study=params.study, **common_kwargs)
    else:
        with open(params.text) as f:
            for line in f:
                study = line.strip()
                if study:
                    generate_yaml(study=study, **common_kwargs)


if __name__ == "__main__":
    main()
