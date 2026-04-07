#!/usr/bin/env bash
# ===========================================================================
# build_celegans_ref.sh — 构建 C. elegans 转录组参考索引
#
# 从 WormBase/Ensembl 下载注释文件，使用 gffread 提取转录组序列，
# 构建 Bowtie2 索引，并生成 regions.bed 和 transcript_lengths.tsv。
#
# 前置依赖:
#   - gffread (https://github.com/gpertea/gffread)
#   - bowtie2-build
#   - samtools (用于 faidx)
#
# 用法:
#   cd workflow/snakescale
#   bash scripts/build_celegans_ref.sh
# ===========================================================================

set -euo pipefail

# ── 配置 ──────────────────────────────────────────────────────────────
SPECIES="celegans"
REF_BASE="reference/transcriptome/${SPECIES}"
GENOME_URL="https://ftp.ensembl.org/pub/release-112/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
GFF_URL="https://ftp.ensembl.org/pub/release-112/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.112.gff3.gz"
INDEX_PREFIX="${REF_BASE}/c_elegans_ref_selected"
REGIONS_BED="${REF_BASE}/c_elegans_ref_actual_regions.bed"
LENGTHS_TSV="${REF_BASE}/c_elegans_ref_transcript_lengths.tsv"

TMPDIR="${REF_BASE}/_build_tmp"

# ── 创建目录 ──────────────────────────────────────────────────────────
mkdir -p "${REF_BASE}" "${TMPDIR}"

echo "=== Step 1: 下载基因组和注释文件 ==="
if [ ! -f "${TMPDIR}/genome.fa" ]; then
    wget -q -O "${TMPDIR}/genome.fa.gz" "${GENOME_URL}"
    gunzip "${TMPDIR}/genome.fa.gz"
fi

if [ ! -f "${TMPDIR}/annotation.gff3" ]; then
    wget -q -O "${TMPDIR}/annotation.gff3.gz" "${GFF_URL}"
    gunzip "${TMPDIR}/annotation.gff3.gz"
fi

echo "=== Step 2: 提取转录组序列 (gffread) ==="
# 提取 mRNA 转录本序列
gffread "${TMPDIR}/annotation.gff3" \
    -g "${TMPDIR}/genome.fa" \
    -w "${TMPDIR}/transcripts.fa" \
    -M -K

echo "=== Step 3: 构建 Bowtie2 索引 ==="
bowtie2-build \
    --threads 8 \
    "${TMPDIR}/transcripts.fa" \
    "${INDEX_PREFIX}"

echo "=== Step 4: 生成 regions.bed (CDS/UTR 区域) ==="
# 从 GFF3 提取 CDS, five_prime_UTR, three_prime_UTR 区域
# 输出 BED 格式: transcript_id  start  end  region_type
awk -F'\t' '
$3 == "CDS" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" {
    # 提取 Parent transcript ID
    match($9, /Parent=transcript:([^;]+)/, arr)
    if (arr[1] != "") {
        # BED 格式 (0-based start)
        print arr[1] "\t" ($4 - 1) "\t" $5 "\t" $3
    }
}
' "${TMPDIR}/annotation.gff3" > "${REGIONS_BED}"

echo "=== Step 5: 生成 transcript_lengths.tsv ==="
# 使用 samtools faidx 获取每条转录本的长度
samtools faidx "${TMPDIR}/transcripts.fa"
awk -F'\t' '{print $1 "\t" $2}' "${TMPDIR}/transcripts.fa.fai" > "${LENGTHS_TSV}"

echo "=== Step 6: 清理临时文件 ==="
rm -rf "${TMPDIR}"

echo ""
echo "=== 完成! 请更新 scripts/references.yaml 中 caenorhabditis elegans 的路径: ==="
echo "  transcriptome: \"transcriptome/${SPECIES}/c_elegans_ref_selected*\""
echo "  regions: \"transcriptome/${SPECIES}/c_elegans_ref_actual_regions.bed\""
echo "  transcript_lengths: \"transcriptome/${SPECIES}/c_elegans_ref_transcript_lengths.tsv\""
