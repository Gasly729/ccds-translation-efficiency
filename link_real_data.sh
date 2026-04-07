#!/usr/bin/env bash
# ============================================================================
# link_real_data.sh — 将服务器公共归档目录中的真实数据软链接到项目 data/raw/
#
# 严禁使用 cp 复制数据，全部使用 ln -s 软链接映射。
# 执行前请根据实际情况修改下方 SOURCE_* 变量。
#
# 用法: bash link_real_data.sh
# ============================================================================
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ★ 请将以下路径替换为服务器上的真实归档路径 ★                          ║
# ╚══════════════════════════════════════════════════════════════════════════╝
SOURCE_FASTQ="/data/archive/silkworm_spider_riboseq"
SOURCE_GENOME="/data/archive/genomes"

# ── 项目内目标目录 ───────────────────────────────────────────────────────────
TARGET_FASTQ="${PROJECT_ROOT}/data/raw/fastq"
TARGET_GENOME="${PROJECT_ROOT}/data/raw/genome"

# ── 辅助函数 ─────────────────────────────────────────────────────────────────
safe_link() {
    # safe_link <源文件> <目标目录>
    # 如果目标已存在且是同一软链接则跳过，否则创建
    local src="$1"
    local dst_dir="$2"
    local dst="${dst_dir}/$(basename "${src}")"

    if [ ! -e "${src}" ]; then
        echo "  [警告] 源文件不存在，跳过: ${src}"
        return 0
    fi

    if [ -L "${dst}" ]; then
        local existing_target
        existing_target=$(readlink -f "${dst}" 2>/dev/null || true)
        local real_src
        real_src=$(readlink -f "${src}" 2>/dev/null || true)
        if [ "${existing_target}" = "${real_src}" ]; then
            echo "  [跳过] 已存在相同链接: $(basename "${src}")"
            return 0
        else
            echo "  [更新] 链接目标已变更，重新创建: $(basename "${src}")"
            rm -f "${dst}"
        fi
    elif [ -e "${dst}" ]; then
        echo "  [错误] 目标位置存在非软链接文件，拒绝覆盖: ${dst}"
        echo "         请手动检查并删除后重试。"
        return 1
    fi

    ln -s "${src}" "${dst}"
    echo "  [链接] $(basename "${src}") -> ${src}"
}

# ── 前置检查 ─────────────────────────────────────────────────────────────────
echo "============================================================"
echo " link_real_data.sh — 软链接真实数据到项目目录"
echo "============================================================"
echo ""

# 检查源目录是否存在
for dir_var in SOURCE_FASTQ SOURCE_GENOME; do
    dir_val="${!dir_var}"
    if [ ! -d "${dir_val}" ]; then
        echo "[警告] 数据源目录不存在: ${dir_val}"
        echo "       请修改脚本顶部的 ${dir_var} 变量为正确路径。"
        echo "       当前以 dry-run 模式继续，不会创建任何链接。"
    fi
done

# 确保目标目录存在
mkdir -p "${TARGET_FASTQ}"
mkdir -p "${TARGET_GENOME}"

# ── 1. 链接 FASTQ 文件 ──────────────────────────────────────────────────────
echo ""
echo ">>> 链接 FASTQ 测序数据..."
echo "    源目录: ${SOURCE_FASTQ}"
echo "    目标:   ${TARGET_FASTQ}"
echo ""

if [ -d "${SOURCE_FASTQ}" ]; then
    # 链接所有 .fastq.gz / .fq.gz 文件
    found_fastq=0
    for fq in "${SOURCE_FASTQ}"/*.fastq.gz "${SOURCE_FASTQ}"/*.fq.gz; do
        [ -e "${fq}" ] || continue
        safe_link "${fq}" "${TARGET_FASTQ}"
        found_fastq=1
    done

    # 也支持按子目录组织的数据 (例如按样本名建立子目录)
    for subdir in "${SOURCE_FASTQ}"/*/; do
        [ -d "${subdir}" ] || continue
        for fq in "${subdir}"*.fastq.gz "${subdir}"*.fq.gz; do
            [ -e "${fq}" ] || continue
            safe_link "${fq}" "${TARGET_FASTQ}"
            found_fastq=1
        done
    done

    if [ "${found_fastq}" -eq 0 ]; then
        echo "  [警告] 在 ${SOURCE_FASTQ} 下未找到 .fastq.gz / .fq.gz 文件"
    fi
else
    echo "  [跳过] 源目录不存在: ${SOURCE_FASTQ}"
fi

# ── 2. 链接参考基因组文件 ───────────────────────────────────────────────────
echo ""
echo ">>> 链接参考基因组文件..."
echo "    源目录: ${SOURCE_GENOME}"
echo "    目标:   ${TARGET_GENOME}"
echo ""

if [ -d "${SOURCE_GENOME}" ]; then
    # 链接常见基因组文件后缀: .fa, .fa.fai, .gtf, .gff, .gff3, .bed, .dict
    for ext in fa fa.fai fasta fasta.fai gtf gff gff3 bed dict; do
        for gf in "${SOURCE_GENOME}"/*.${ext}; do
            [ -e "${gf}" ] || continue
            safe_link "${gf}" "${TARGET_GENOME}"
        done
    done

    # 支持 STAR / HISAT2 预建索引目录
    for idx_dir in "${SOURCE_GENOME}"/star_index "${SOURCE_GENOME}"/hisat2_index; do
        if [ -d "${idx_dir}" ]; then
            local_name="${TARGET_GENOME}/$(basename "${idx_dir}")"
            if [ -L "${local_name}" ]; then
                echo "  [跳过] 索引目录链接已存在: $(basename "${idx_dir}")"
            else
                ln -s "${idx_dir}" "${local_name}"
                echo "  [链接] $(basename "${idx_dir}")/ -> ${idx_dir}"
            fi
        fi
    done
else
    echo "  [跳过] 源目录不存在: ${SOURCE_GENOME}"
fi

# ── 3. 验证摘要 ─────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " 链接完成。当前 data/raw/ 文件列表:"
echo "============================================================"
echo ""
echo "--- FASTQ ---"
ls -la "${TARGET_FASTQ}/" 2>/dev/null || echo "  (空)"
echo ""
echo "--- Genome ---"
ls -la "${TARGET_GENOME}/" 2>/dev/null || echo "  (空)"
echo ""
echo "提示: 如需更改数据源，请编辑脚本顶部的 SOURCE_FASTQ / SOURCE_GENOME 变量。"
