# ============================================================================
# Makefile — CCDS + Snakescale + TE_model 项目顶层调度 (服务器生产版)
#
# 目标:
#   make init            → 运行 setup_project.sh 创建目录骨架与占位文件
#   make link_data       → 软链接服务器归档数据到 data/raw/
#   make download_data   → 从 SRA 下载 FASTQ 数据 (SRX→SRR 映射 + sradownloader)
#   make run_snakescale  → 通过 SLURM 集群投递 Snakemake 流水线
#   make run_local       → 本地模式运行 Snakemake (仅限交互节点调试)
#   make calc_te         → 执行 Python TE 计算 (或生成假数据)
#   make plot            → 执行 R 脚本生成图表
#   make setup_slurm     → 初始化 Snakemake SLURM profile
#   make all             → 按顺序执行全部步骤 (集群模式)
#   make dry_run         → Snakemake 干运行，仅验证 DAG
#   make clean           → 清理中间与输出文件
# ============================================================================

.PHONY: all init link_data download_data download_prepare download_fetch \
        run_snakescale run_local calc_te plot \
        setup_slurm dry_run clean help check_not_login_node

# ── 路径配置 ─────────────────────────────────────────────────────────────────
PROJECT_ROOT   := $(shell pwd)
SNAKESCALE_DIR := workflow/snakescale
PROCESSED_DIR  := data/processed
INTERIM_DIR    := data/interim
FIGURES_DIR    := reports/figures
SLURM_PROFILE  := workflow/snakescale/slurm_profile
EXTERNAL_DIR   := data/external
RAW_FASTQ_DIR  := data/raw/fastq
LOG_DIR        := logs

# ── 工具配置 ─────────────────────────────────────────────────────────────────
PYTHON    := python3
RSCRIPT   := Rscript
SNAKEMAKE := snakemake

# ── SLURM 集群投递参数 ──────────────────────────────────────────────────────
# 分区名、账户等参数与 config.yaml 中 slurm 段保持一致
SLURM_PARTITION  ?= compute
SLURM_JOBS       ?= 50
SLURM_LATENCY    ?= 30

# Snakemake 通用参数
SNAKEMAKE_COMMON := --printshellcmds --use-conda --rerun-incomplete \
                    --latency-wait $(SLURM_LATENCY)

# 集群投递模式参数 (使用 --cluster sbatch)
SNAKEMAKE_CLUSTER := $(SNAKEMAKE_COMMON) \
    --jobs $(SLURM_JOBS) \
    --cluster "sbatch \
        --partition=$(SLURM_PARTITION) \
        --cpus-per-task={threads} \
        --mem={resources.mem_mb}M \
        --time={resources.runtime_min} \
        --output=logs/slurm/%j_%x.out \
        --error=logs/slurm/%j_%x.err \
        --job-name=snkscl_{rule}" \
    --default-resources mem_mb=8000 runtime_min=120 disk_mb=50000

# 本地调试模式参数
SNAKEMAKE_LOCAL := $(SNAKEMAKE_COMMON) --cores all

# ============================================================================
#  安全检查: 防止在登录节点直接运行计算密集型任务
# ============================================================================
# 通过检测 SLURM_JOB_ID (计算节点) 或 hostname 模式来判断。
# 如果你的登录节点命名有特定模式 (如 login*), 请按需修改。
HOSTNAME := $(shell hostname)

check_not_login_node:
	@if echo "$(HOSTNAME)" | grep -qiE '^login|^head|^master'; then \
		echo ""; \
		echo "╔══════════════════════════════════════════════════════════════╗"; \
		echo "║  ⚠  警告: 检测到当前处于登录节点 ($(HOSTNAME))              ║"; \
		echo "║                                                            ║"; \
		echo "║  请勿在登录节点直接运行计算密集型任务！                      ║"; \
		echo "║  推荐方式:                                                  ║"; \
		echo "║    1. make run_snakescale  (自动提交到 SLURM 计算节点)      ║"; \
		echo "║    2. srun --pty bash → 然后 make run_local (交互节点调试)  ║"; \
		echo "║                                                            ║"; \
		echo "║  如确需在此节点运行，请使用: make run_local FORCE=1          ║"; \
		echo "╚══════════════════════════════════════════════════════════════╝"; \
		echo ""; \
		if [ "$(FORCE)" != "1" ]; then exit 1; fi; \
	fi

# ── 默认目标 ─────────────────────────────────────────────────────────────────
all: init link_data download_data run_snakescale calc_te plot
	@echo "=========================================="
	@echo " 全部步骤执行完毕 (集群投递模式)"
	@echo "=========================================="

# ── 1. 初始化: 创建 CCDS 目录骨架与占位文件 ─────────────────────────────────
init:
	@echo ">>> [init] 创建 CCDS 目录骨架..."
	bash setup_project.sh
	@mkdir -p logs/slurm
	@echo ">>> [init] 完成"

# ── 2. 软链接真实数据 ───────────────────────────────────────────────────────
link_data:
	@echo ">>> [link_data] 软链接服务器归档数据..."
	bash link_real_data.sh
	@echo ">>> [link_data] 完成"

# ── 3. 下载 SRA 数据 (SRX→SRR 映射 + sradownloader 下载) ────────────────────
# download_data     = 完整流程 (prepare + fetch)
# download_prepare  = 仅生成映射文件，不下载
# download_fetch    = 仅下载 (假设映射文件已就绪)
download_data: download_prepare download_fetch

download_prepare:
	@echo ">>> [download_prepare] 生成 SRX→SRR 映射与 sradownloader 输入文件..."
	@mkdir -p $(EXTERNAL_DIR) $(LOG_DIR)/download
	$(PYTHON) -m src.data.download_sra --prepare
	@echo ">>> [download_prepare] 完成"
	@echo "    映射表:   $(EXTERNAL_DIR)/srx_to_srr_mapping.csv"
	@echo "    输入文件: $(EXTERNAL_DIR)/sradownloader_input.txt"

download_fetch:
	@echo ">>> [download_fetch] 调用 sradownloader 下载 FASTQ..."
	@mkdir -p $(RAW_FASTQ_DIR) $(LOG_DIR)/download
	$(PYTHON) -m src.data.download_sra --download --outdir $(RAW_FASTQ_DIR)
	@echo ">>> [download_fetch] 完成 → FASTQ 输出至 $(RAW_FASTQ_DIR)/"

# ── 4. 初始化 Snakemake SLURM profile ──────────────────────────────────────
# 创建一个本地 profile 目录，Snakemake 可通过 --profile 引用
setup_slurm:
	@echo ">>> [setup_slurm] 创建 SLURM profile..."
	@mkdir -p $(SLURM_PROFILE)
	@echo "cluster:" > $(SLURM_PROFILE)/config.yaml
	@echo "  mkdir -p logs/slurm && sbatch" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --partition=$(SLURM_PARTITION)" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --cpus-per-task={threads}" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --mem={resources.mem_mb}M" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --time={resources.runtime_min}" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --output=logs/slurm/%j_%x.out" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --error=logs/slurm/%j_%x.err" >> $(SLURM_PROFILE)/config.yaml
	@echo "    --job-name=snkscl_{rule}" >> $(SLURM_PROFILE)/config.yaml
	@echo "jobs: $(SLURM_JOBS)" >> $(SLURM_PROFILE)/config.yaml
	@echo "default-resources: [mem_mb=8000, runtime_min=120, disk_mb=50000]" >> $(SLURM_PROFILE)/config.yaml
	@echo "latency-wait: $(SLURM_LATENCY)" >> $(SLURM_PROFILE)/config.yaml
	@echo "printshellcmds: true" >> $(SLURM_PROFILE)/config.yaml
	@echo "rerun-incomplete: true" >> $(SLURM_PROFILE)/config.yaml
	@echo "use-conda: true" >> $(SLURM_PROFILE)/config.yaml
	@echo ">>> [setup_slurm] Profile 已写入 $(SLURM_PROFILE)/config.yaml"

# ── 5. 运行 Snakescale — SLURM 集群投递模式 (生产推荐) ─────────────────────
run_snakescale:
	@echo ">>> [run_snakescale] 以 SLURM 集群投递模式运行..."
	@echo "    分区: $(SLURM_PARTITION) | 最大并行任务: $(SLURM_JOBS)"
	@mkdir -p logs/slurm
	cd $(SNAKESCALE_DIR) && $(SNAKEMAKE) $(SNAKEMAKE_CLUSTER)
	@echo ">>> [run_snakescale] 所有任务已提交/完成"

# ── 6. 运行 Snakescale — 本地模式 (仅限交互节点调试) ───────────────────────
run_local: check_not_login_node
	@echo ">>> [run_local] 以本地模式运行 (调试用)..."
	cd $(SNAKESCALE_DIR) && $(SNAKEMAKE) $(SNAKEMAKE_LOCAL)
	@echo ">>> [run_local] 完成"

# ── 7. Dry Run — 仅验证 DAG，不执行任何任务 ─────────────────────────────────
dry_run:
	@echo ">>> [dry_run] 验证 Snakemake DAG..."
	cd $(SNAKESCALE_DIR) && $(SNAKEMAKE) -n --printshellcmds
	@echo ">>> [dry_run] DAG 验证完成"

# ── 8. 计算翻译效率 (TE) — CenikLab CLR/ILR 成分回归方法 ─────────────────
# 核心算法 100% 复刻自 CenikLab TE_model 原文代码。
# 模式 A: 有 .ribo 文件 → 全流程 (Stage 0-3)
# 模式 B: 有 ribo_raw.csv + rnaseq_raw.csv → 跳过提取 (Stage 1-3)
# 模式 C: 无数据 → --dummy 生成假数据用于管线测试
#
# 可覆盖变量:
#   RIBO_DIR=path/to/ribo_files  - .ribo 文件目录
#   CPM_CUTOFF=1                 - CPM 过滤阈值
#   OVERALL_CUTOFF=70            - 低表达样本百分比阈值
#   NONPOLYA_CSV=path            - 非 polyA 基因列表
#   SKIP_R=1                     - 仅做预处理，跳过 R 脚本
#   METADATA_CSV=path            - 元数据 CSV (启用强约束配对)
#   MAPPING_CSV=path             - SRX→SRR 映射 CSV (配合 METADATA_CSV)
RIBO_DIR         ?=
CPM_CUTOFF       ?= 1
OVERALL_CUTOFF   ?= 70
NONPOLYA_CSV     ?=
SKIP_R           ?=
METADATA_CSV     ?=
MAPPING_CSV      ?=

calc_te:
	@echo ">>> [calc_te] 计算翻译效率 (CenikLab CLR/ILR 方法)..."
	@if [ -n "$(RIBO_DIR)" ] && [ -d "$(RIBO_DIR)" ]; then \
		echo "    模式 A: 从 .ribo 文件执行全流程"; \
		$(PYTHON) -m src.te_calc.te_calculator \
			--ribo_dir $(RIBO_DIR) \
			--output_dir $(PROCESSED_DIR) \
			--cpm_cutoff $(CPM_CUTOFF) \
			--overall_cutoff $(OVERALL_CUTOFF) \
			$(if $(NONPOLYA_CSV),--nonpolya_csv $(NONPOLYA_CSV)) \
			$(if $(METADATA_CSV),--metadata_csv $(METADATA_CSV)) \
			$(if $(MAPPING_CSV),--mapping_csv $(MAPPING_CSV)) \
			$(if $(SKIP_R),--skip_r); \
	elif [ -f "$(PROCESSED_DIR)/ribo_raw.csv" ] && \
	     [ -f "$(PROCESSED_DIR)/rnaseq_raw.csv" ]; then \
		echo "    模式 B: 检测到已有 count 矩阵，跳过 .ribo 提取"; \
		$(PYTHON) -m src.te_calc.te_calculator \
			--skip_extract \
			--output_dir $(PROCESSED_DIR) \
			--cpm_cutoff $(CPM_CUTOFF) \
			--overall_cutoff $(OVERALL_CUTOFF) \
			$(if $(NONPOLYA_CSV),--nonpolya_csv $(NONPOLYA_CSV)) \
			$(if $(METADATA_CSV),--metadata_csv $(METADATA_CSV)) \
			$(if $(MAPPING_CSV),--mapping_csv $(MAPPING_CSV)) \
			$(if $(SKIP_R),--skip_r); \
	else \
		echo "    模式 C: 未检测到上游数据，使用 --dummy 生成测试数据"; \
		$(PYTHON) -m src.te_calc.te_calculator \
			--dummy \
			--output_dir $(PROCESSED_DIR); \
	fi
	@echo ">>> [calc_te] 完成"

# ── 9. 绘图: 执行 R 脚本生成 TE 可视化图表 ─────────────────────────────────
plot:
	@echo ">>> [plot] 生成 TE 特征可视化图表..."
	$(RSCRIPT) src/plots/plot_te_features.R \
		--input  $(PROCESSED_DIR) \
		--output $(FIGURES_DIR) \
		--format pdf
	@echo ">>> [plot] 完成 → 图表输出至 $(FIGURES_DIR)/"

# ── 清理 ────────────────────────────────────────────────────────────────────
clean:
	@echo ">>> [clean] 清理中间文件与输出..."
	rm -rf $(INTERIM_DIR)/*
	rm -rf $(PROCESSED_DIR)/*
	rm -rf $(FIGURES_DIR)/*
	rm -rf logs/slurm/*
	@echo ">>> [clean] 完成 (注意: data/raw/ 中的软链接未被清理)"

# ── 帮助信息 ─────────────────────────────────────────────────────────────────
help:
	@echo "═══════════════════════════════════════════════════════════════"
	@echo " Snakescale + TE Model — 服务器生产环境 Makefile"
	@echo "═══════════════════════════════════════════════════════════════"
	@echo ""
	@echo "  初始化与数据准备:"
	@echo "    make init            - 创建 CCDS 目录骨架与占位文件"
	@echo "    make link_data       - 软链接服务器归档数据到 data/raw/"
	@echo "    make download_data   - ★ SRA 下载完整流程 (映射+下载)"
	@echo "    make download_prepare - 仅生成 SRX→SRR 映射 (不下载)"
	@echo "    make download_fetch  - 仅执行 sradownloader 下载"
	@echo "    make setup_slurm     - 初始化 Snakemake SLURM profile"
	@echo ""
	@echo "  流水线执行:"
	@echo "    make run_snakescale  - ★ SLURM 集群投递模式 (生产推荐)"
	@echo "    make run_local       - 本地模式 (交互节点调试, 禁止登录节点)"
	@echo "    make dry_run         - 干运行, 仅验证 DAG"
	@echo ""
	@echo "  下游分析:"
	@echo "    make calc_te         - 执行 Python TE 计算"
	@echo "    make plot            - 执行 R 脚本生成图表"
	@echo ""
	@echo "  综合与维护:"
	@echo "    make all             - 按顺序执行全部步骤 (集群模式)"
	@echo "    make clean           - 清理中间文件与输出"
	@echo "    make help            - 显示此帮助信息"
	@echo ""
	@echo "  可覆盖变量:"
	@echo "    SLURM_PARTITION=xxx  - 覆盖 SLURM 分区 (默认: compute)"
	@echo "    SLURM_JOBS=N        - 覆盖最大并行任务数 (默认: 50)"
	@echo "    FORCE=1             - 强制允许在登录节点运行 run_local"
