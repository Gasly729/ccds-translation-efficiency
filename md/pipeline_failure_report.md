# Pipeline Failure Report — snakemake-ribo End-to-End Test

**Generated:** 2026-04-07 10:49 UTC+8  
**Pipeline:** `snakemake-ribo` (Snakescale + RiboFlow + TE Calculator)  
**Environment:** `snakemake-ribo` conda env, 208 cores, Linux  
**Log file:** `test_run.log` (1761 lines)  
**Exit code:** 2 (Snakemake), 2 (TE calc)

---

## 1. Test Configuration

| Species | Study | Ribo SRR | RNA SRR | Ribo GSM | RNA GSM |
|---|---|---|---|---|---|
| *Caenorhabditis elegans* | GSE48140 | SRR914339 | SRR914333 | GSM1169554 | GSM1169548 |
| *Arabidopsis thaliana* | GSE50597 | SRR966477 | SRR966478 | GSM1224478 | GSM1224479 |
| *Homo sapiens* | GSE49994 | SRR953773 | SRR953772 | GSM1211429 | GSM1211428 |
| *Saccharomyces cerevisiae* | GSE125038 | SRR8439615 | SRR8439643 | GSM3561545 | GSM3561573 |

All 8 FASTQ files confirmed present in `data/raw/fastq/` as symlinks to `/home/xrx/raw_data/`.

---

## 2. Overall Results Summary

| Study | Species | Pre-Run Check | RiboFlow | .ribo File | TE Calc |
|---|---|---|---|---|---|
| GSE125038 | S. cerevisiae | ✅ PASS | ✅ 33/33 processes | ✅ `all.ribo` (6.7 MB) | ❌ Column mismatch |
| GSE48140 | C. elegans | ✅ PASS | ❌ `create_ribo` collision | ❌ N/A | ❌ N/A |
| GSE50597 | A. thaliana | ✅ PASS | ❌ trace.txt race condition | ❌ N/A | ❌ N/A |
| GSE49994 | H. sapiens | ❌ FAIL | ⏭️ Skipped | ❌ N/A | ❌ N/A |

**Snakemake DAG:** 30 jobs planned, 27 completed (90%), 3 failed.  
**Net result:** 1 of 4 studies produced .ribo output; TE calculation failed for all.

---

## 3. Error Details by Stage

### 3.1 ❌ GSE49994 (Homo sapiens) — Pre-Run Check Failure

**Stage:** `classify_studies` (Snakemake rule)  
**Cause:** Adapter detection found low adapter presence for Ribo-Seq sample SRR953773, and the adapter guessing algorithm returned `None`.

**Raw log:**
```
Ribo-Seq files with low adapter presence:
SRR953773
Unable to guess adapter. Study invalid.
Guessed adapters: [None]
```

**Analysis:**  
The `check_adapter` rule ran `cutadapt --report=minimal` on a subsample of SRR953773 and found adapter trimming rate below the threshold (`adapter_threshold` in config.yaml, default 10%). The fallback `guess_adapter` routine (k-mer frequency analysis) could not identify a consensus 3' adapter. Because `has_all_uneven_lengths` was false and no adapter was guessed, the study was marked invalid per the classification logic in the Snakefile (lines 634-640).

**Root cause:** The selected Ribo-Seq sample (SRR953773, Hep3B RNC-mRNA) may use a non-standard library preparation without a typical 3' adapter, or the adapter sequence in metadata is incorrect/missing. The metadata CSV contains no `threep_adapter` for this sample.

**Fix options:**
1. Set `override: True` in `config.yaml` to force RiboFlow execution despite failed checks.
2. Select a different human Ribo-Seq sample with a known adapter.
3. Manually specify the correct adapter in `metadata.csv`.

---

### 3.2 ❌ GSE50597 (Arabidopsis thaliana) — Nextflow Trace File Race Condition

**Stage:** `run_riboflow` (Nextflow launch)  
**Cause:** Multiple Nextflow instances launched simultaneously and collided on the shared `nextflow_logs/trace.txt` file.

**Raw log:**
```
ERROR ~ No such file: /home/xrx/my_project/project/workflow/snakescale/nextflow_logs/trace.txt

 -- Check '.nextflow.log' file for details
[Tue Apr  7 10:19:23 2026]
Error in rule run_riboflow:
    jobid: 0
    input: ../../data/interim/modified_project/GSE50597/GSE50597_modified.yaml
    output: ../../data/processed/riboflow_status/GSE50597/riboflow_status.txt

RuleException:
CalledProcessError in file Snakefile, line 712:
Command 'nextflow riboflow/RiboFlow.groovy -params-file ...GSE50597_modified.yaml
  -profile stampede_local' returned non-zero exit status 1.
```

**Analysis:**  
Snakemake launched 3 `run_riboflow` jobs in parallel (GSE50597, GSE48140, GSE125038). All 3 Nextflow instances are configured to write trace/report/timeline to the **same** directory `./nextflow_logs/` (via `nextflow.config`). The first instance to start created the file; the second instance failed because Nextflow 19.04.1 does not handle concurrent file access gracefully when `overwrite = true`.

**Root cause:** `riboflow/nextflow.config` hardcodes shared output paths:
```groovy
trace {
    enabled = true
    file = "./nextflow_logs/trace.txt"
    overwrite = true
}
```

**Fix options:**
1. **Serialize RiboFlow calls** by setting `threads: 48` for `run_riboflow` rule (prevents parallel launches).
2. **Per-study trace paths:** Modify the Snakefile to set `NXF_WORK` and trace output to a per-study directory (e.g., `nextflow_logs/{study}/`).
3. Pass `-with-trace nextflow_logs/{study}_trace.txt` on the Nextflow command line in `run_riboflow`.

---

### 3.3 ❌ GSE48140 (Caenorhabditis elegans) — create_ribo Reference Collision

**Stage:** `run_riboflow` → Nextflow `create_ribo` process  
**Cause:** Nextflow detected multiple input files staged to the same name `reference` in the process work directory.

**Raw log:**
```
ERROR ~ Error executing process > 'create_ribo (1)'

Caused by:
  Process `create_ribo` input file name collision -- There are multiple input
  files for each of the following file names: reference
```

**Generated YAML reference paths for GSE48140:**
```yaml
input:
  reference:
    filter: reference/filter/celegans/celegans_rRNA_new*
    transcriptome: reference/transcriptome/celegans/celegans_rRNA_new*
    regions: reference/          # ← EMPTY → resolves to bare directory
    transcript_lengths: reference/   # ← EMPTY → resolves to bare directory
```

**Compared with working GSE125038 (S. cerevisiae):**
```yaml
input:
  reference:
    filter: reference/filter/saccharomyces_cerevisiae/saccharomyces_cerevisiae_rtRNA*
    transcriptome: reference/transcriptome/saccharomyces_cerevisiae/s_cerevisiae_ref_selected*
    regions: reference/transcriptome/saccharomyces_cerevisiae/s_cerevisiae_ref_actual_regions.bed
    transcript_lengths: reference/transcriptome/saccharomyces_cerevisiae/s_cerevisiae_ref_transcript_lengths.tsv
```

**Root cause:** `scripts/references.yaml` has **incomplete** C. elegans reference data:
```yaml
caenorhabditis elegans:
  folder_name: celegans
  filter: filter/celegans/celegans_rRNA_new*
  transcriptome: transcriptome/celegans/celegans_rRNA_new*   # ← Same as filter!
  regions: ""                                                  # ← MISSING
  transcript_lengths: ""                                       # ← MISSING
```

Additionally, the `transcriptome/celegans/` directory contains only rRNA filter files (`celegans_rRNA_new.*`), not actual transcriptome reference files. The C. elegans APPRIS/RefSeq transcriptome index has not been built.

**Fix:**
1. Build proper C. elegans transcriptome Bowtie2 indices (e.g., from WormBase canonical transcripts).
2. Create `regions.bed` and `transcript_lengths.tsv` files for C. elegans.
3. Update `scripts/references.yaml` with correct paths.

---

### 3.4 ❌ TE Calculation — Column Name Mismatch

**Stage:** `make calc_te` → `src.te_calc.te_calculator` Stage 0.5  
**Cause:** The .ribo file uses GSM-based experiment names internally, but the TE calculator's strong-pairing mode expects SRR-based column names.

**Raw log:**
```
RuntimeError: 严重错误: 没有任何 Ribo↔RNA 配对能在 count 矩阵中找到对应列!
  ribo_raw.csv 列名 (前10): ['all']
  rnaseq_raw.csv 列名 (前10): ['all']
  配对字典中的 SRR 示例: {'ribo_gsm': 'GSM1224478', ..., 'ribo_srr': ['SRR966477'], ...}
```

**Analysis:**  
- `ribopy` extracted from `all.ribo` produces count matrices with column name `"all"` (the merged .ribo filename stem).
- The strong-pairing logic in `validate_and_align_columns()` (line 444) searches for SRR IDs (e.g., `SRR8439615`) or GSM IDs in the column names. Neither `"all"` nor `"GSM3561545"` matches.
- Even if individual per-experiment .ribo files were used (with column `"GSM3561545"`), the TE calculator would still fail because it looks for SRR IDs first, then GSM IDs, but the rename logic may not handle the `all.ribo` → `GSM3561545` mapping correctly.
- Additionally, only 1 study produced valid output, meaning only 1 Ribo + 1 RNA sample are available. The CLR/ILR compositional regression requires ≥2 samples per component.

**Fix options:**
1. Use individual per-experiment .ribo files instead of the merged `all.ribo`.
2. Modify `te_calculator.py` `validate_and_align_columns()` to also match GSM-based column names.
3. Ensure ≥2 successful studies before attempting TE calculation.

---

## 4. Successful Components

Despite the failures, significant portions of the pipeline executed correctly:

| Component | Status | Details |
|---|---|---|
| YAML generation | ✅ | All 4 per-study YAMLs generated from CSV metadata |
| Adapter detection | ✅ | `check_adapter` completed for all 8 SRR files |
| Length detection | ✅ | `check_lengths` completed for all 4 studies |
| Adapter guessing | ✅ | `guess_adapter` completed for all 4 studies |
| Pre-run classification | ✅ | Correctly identified 3 valid + 1 invalid study |
| GSE125038 RiboFlow | ✅ | Full 33-process pipeline: clip → filter → align → create_ribo → inject_rnaseq → merge |
| .ribo output | ✅ | `GSM3561545.ribo` (S. cerevisiae, 6027 genes, Ribo+RNA) |
| Count extraction | ✅ | `ribo_raw.csv` and `rnaseq_raw.csv` produced (6027 genes) |

---

## 5. Recommended Remediation Priority

| Priority | Action | Impact | Effort |
|---|---|---|---|
| 🔴 P0 | Fix C. elegans references (build transcriptome index, regions, lengths) | Unblocks 1 species | High |
| 🔴 P0 | Fix Nextflow trace race condition (per-study output dirs) | Unblocks parallel runs | Low |
| 🟡 P1 | Fix TE calculator GSM↔SRR column mapping | Unblocks downstream | Medium |
| 🟡 P1 | Set `override: True` or fix H. sapiens adapter | Unblocks 1 species | Low |
| 🟢 P2 | Add `--keep-going` to Snakemake for resilient multi-study runs | Improves robustness | Low |

---

## 6. Appendix: Key File Locations

- **Full run log:** `test_run.log`
- **Snakemake log:** `workflow/snakescale/.snakemake/log/2026-04-07T101839.068095.snakemake.log`
- **Nextflow logs:** `workflow/snakescale/nextflow_logs/`
- **Pre-run status:** `data/interim/log/status.txt`
- **Valid studies:** `data/interim/log/valid_studies.txt`
- **Generated YAMLs:** `data/interim/project/{GSE}/{GSE}.yaml`
- **Modified YAMLs:** `data/interim/modified_project/{GSE}/{GSE}_modified.yaml`
- **Successful .ribo:** `workflow/snakescale/output/GSE125038/ribo/all.ribo`
- **Count matrices:** `data/processed/ribo_raw.csv`, `data/processed/rnaseq_raw.csv`
- **Config backup:** `workflow/snakescale/config.yaml.bak`
- **Test metadata:** `data/external/test_metadata.csv`, `data/external/test_srx_to_srr_mapping.csv`
