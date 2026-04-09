#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_pub_style_figures.py — 生成出版级 TE 分析图形
=====================================================================
读取 te_results_final.csv 和 sample_annotation.csv，生成 4 张公程式分析图形。

主要依赖:
    matplotlib, seaborn, sklearn, numpy, pandas

用法示例:
    python src/plots/generate_pub_style_figures.py \\
        --output-dir reports/figures/pub_style

注意: te_results_final.csv 是跨物种 outer join 结果，含大量 NaN。
         各图表处理策略见各 plot_fig* 函数内注释。
"""

from pathlib import Path
import argparse
import inspect
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import font_manager
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns

SPECIES_COLORS = {
    'human': '#E64B35',
    'mouse': '#4DBBD5',
    'arabidopsis': '#00A087',
    'yeast': '#F39B7F',
    'unknown': '#B0B0B0',
}

BIO_GROUP_COLORS = {
    'Human cell lines': '#E64B35',
    'Human immune cells': '#F39B7F',
    'Human tissues/other': '#C4705A',
    'Mouse neural': '#4DBBD5',
    'Mouse other tissues': '#91D1C2',
    'Plants': '#00A087',
    'Yeast': '#B5C7A3',
    'Unknown': '#B0B0B0',
}

BIO_GROUP_ORDER = [
    'Human cell lines',
    'Human immune cells',
    'Human tissues/other',
    'Mouse neural',
    'Mouse other tissues',
    'Plants',
    'Yeast',
    'Unknown',
]



def find_project_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here.parent] + list(here.parents):
        if (candidate / 'data/processed/te_results_final.csv').exists() and (candidate / 'data/processed/sample_annotation.csv').exists():
            return candidate
    raise FileNotFoundError('Could not locate project root from script location')



def configure_rcparams() -> str:
    available_fonts = {f.name for f in font_manager.fontManager.ttflist}
    font_choice = 'Arial' if 'Arial' in available_fonts else 'DejaVu Sans'
    matplotlib.rcParams.update({
        'font.family': font_choice,
        'font.size': 8,
        'axes.titlesize': 9,
        'axes.labelsize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'legend.frameon': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'figure.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })
    sns.set_style('white')
    return font_choice



def assign_bio_group(row: pd.Series) -> str:
    sp = str(row['species']).lower()
    ct = str(row['cell_type']).lower()
    tissue = str(row['tissue']).lower()
    if sp == 'human':
        if any(k in ct for k in ['hela', 'hek', 'hek293', 'hek293t', 'u2-os', 'a549', 'sw480', 'caco', 'bf', 'mcf', 'cfbe', 'nci', 'huh', 'bj', 'hl-60', 'k562', '293']):
            return 'Human cell lines'
        if any(k in tissue for k in ['macrophage', 'b cell', 'lymph', 'splenic']):
            return 'Human immune cells'
        return 'Human tissues/other'
    if sp == 'mouse':
        if any(k in tissue for k in ['brain', 'hippocampus', 'neuron', 'dentate']):
            return 'Mouse neural'
        return 'Mouse other tissues'
    if sp == 'arabidopsis':
        return 'Plants'
    if sp == 'yeast':
        return 'Yeast'
    return 'Unknown'



def load_data(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    te = pd.read_csv(project_root / 'data/processed/te_results_final.csv', index_col=0)
    anno = pd.read_csv(project_root / 'data/processed/sample_annotation.csv').set_index('sample_id')
    te = te[te.index != 'all']
    anno = anno[anno.index != 'all']
    common_idx = te.index.intersection(anno.index)
    te = te.loc[common_idx]
    anno = anno.loc[common_idx]
    anno['bio_group'] = anno.apply(assign_bio_group, axis=1)
    return te, anno



def add_confidence_ellipse(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    if len(x) < 4:
        return
    cov = np.cov(x, y)
    if not np.isfinite(cov).all() or np.linalg.matrix_rank(cov) < 2:
        return
    vals, vecs = np.linalg.eigh(cov)
    vals = np.clip(vals, a_min=0, a_max=None)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    scale = np.sqrt(3.2188758248682006)
    width, height = 2 * scale * np.sqrt(vals)
    if width == 0 or height == 0:
        return
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    ellipse = Ellipse(
        xy=(float(np.mean(x)), float(np.mean(y))),
        width=float(width),
        height=float(height),
        angle=float(angle),
        facecolor=color,
        edgecolor='none',
        alpha=0.12,
        zorder=1,
    )
    ax.add_patch(ellipse)



def spread_overlapping_points(x: np.ndarray, y: np.ndarray, labels: list[str]) -> tuple[np.ndarray, np.ndarray]:
    points = np.column_stack([x.astype(float), y.astype(float)])
    if len(points) <= 1:
        return points[:, 0], points[:, 1]
    span = max(float(np.ptp(points[:, 0])), float(np.ptp(points[:, 1])))
    cluster_tol = max(span * 0.0005, 0.05)
    offset_step = max(span * 0.015, 0.35)
    display = points.copy()
    visited = np.zeros(len(points), dtype=bool)
    for start in range(len(points)):
        if visited[start]:
            continue
        stack = [start]
        visited[start] = True
        cluster: list[int] = []
        while stack:
            idx = stack.pop()
            cluster.append(idx)
            distances = np.linalg.norm(points - points[idx], axis=1)
            neighbors = np.where((distances <= cluster_tol) & (~visited))[0]
            if len(neighbors) == 0:
                continue
            visited[neighbors] = True
            stack.extend(neighbors.tolist())
        if len(cluster) <= 1:
            continue
        cluster = sorted(cluster, key=lambda idx: labels[idx])
        center = points[cluster].mean(axis=0)
        for j, idx in enumerate(cluster):
            ring = j // 8
            pos = j % 8
            slots = min(len(cluster) - ring * 8, 8)
            angle = 2 * np.pi * pos / max(slots, 1)
            radius = offset_step * (1 + 0.45 * ring)
            display[idx, 0] = center[0] + radius * np.cos(angle)
            display[idx, 1] = center[1] + radius * np.sin(angle)
    return display[:, 0], display[:, 1]



def save_figure(fig: plt.Figure, stem: str, output_dir: Path) -> list[Path]:
    png_path = output_dir / f'{stem}.png'
    pdf_path = output_dir / f'{stem}.pdf'
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return [png_path, pdf_path]



def plot_fig1(te: pd.DataFrame, anno: pd.DataFrame, output_dir: Path) -> tuple[list[Path], np.ndarray, np.ndarray]:
    # te_results_final.csv 是跨物种 outer join，非本物种基因列为 NaN。
    # PCA 不接受 NaN，此处用 fillna(0) 将缺失基因位置补 0。
    # 这对 PCA 是合理的：跨物种基因本来不可比较，补 0 将其视为“不表达”而非虚假信号。
    x = StandardScaler().fit_transform(te.fillna(0).values)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(x)
    var = pca.explained_variance_ratio_ * 100
    plot_df = anno.copy()
    plot_df['PC1'] = coords[:, 0]
    plot_df['PC2'] = coords[:, 1]
    plot_df['PC1_display'], plot_df['PC2_display'] = spread_overlapping_points(
        plot_df['PC1'].to_numpy(dtype=float),
        plot_df['PC2'].to_numpy(dtype=float),
        plot_df.index.astype(str).tolist(),
    )
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    species_order = [k for k in SPECIES_COLORS if k in plot_df['species'].astype(str).str.lower().unique()]
    bio_group_order = [g for g in BIO_GROUP_ORDER if g in plot_df['bio_group'].unique()]
    panel_specs = [
        (axes[0], 'species', SPECIES_COLORS, species_order, 'By species', True),
        (axes[1], 'bio_group', BIO_GROUP_COLORS, bio_group_order, 'By biological group', False),
    ]
    for ax, column, palette, order, title, lowercase_match in panel_specs:
        for group in order:
            if lowercase_match:
                sub = plot_df[plot_df[column].astype(str).str.lower() == group]
            else:
                sub = plot_df[plot_df[column] == group]
            if sub.empty:
                continue
            ax.scatter(
                sub['PC1_display'],
                sub['PC2_display'],
                s=28,
                alpha=0.85,
                c=palette[group],
                edgecolors='white',
                linewidths=0.4,
                label=f'{group} (n={len(sub)})',
                zorder=3,
            )
            add_confidence_ellipse(ax, sub['PC1'].to_numpy(), sub['PC2'].to_numpy(), palette[group])
        ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
        ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
        ax.set_title(title, fontweight='bold')
        ax.legend(loc='best', handlelength=1.0, borderaxespad=0.2)
    all_x = plot_df['PC1_display'].to_numpy(dtype=float)
    all_y = plot_df['PC2_display'].to_numpy(dtype=float)
    x_span = float(all_x.max() - all_x.min())
    y_span = float(all_y.max() - all_y.min())
    x_pad = x_span * 0.15 if x_span > 0 else 1.0
    y_pad = y_span * 0.15 if y_span > 0 else 1.0
    shared_xlim = (float(all_x.min() - x_pad), float(all_x.max() + x_pad))
    shared_ylim = (float(all_y.min() - y_pad), float(all_y.max() + y_pad))
    for ax in axes:
        ax.set_xlim(shared_xlim)
        ax.set_ylim(shared_ylim)
    return save_figure(fig, 'Fig1_PCA', output_dir), var, coords



def plot_fig2(te: pd.DataFrame, anno: pd.DataFrame, output_dir: Path) -> list[Path]:
    sample_stats = pd.DataFrame({
        'sample_id': te.index,
        'median': te.median(axis=1),
        'q1': te.quantile(0.25, axis=1),
        'q3': te.quantile(0.75, axis=1),
    }).join(anno[['bio_group']])
    group_order = [g for g in BIO_GROUP_ORDER if g in sample_stats['bio_group'].unique()]
    ordered_frames = []
    separator_positions = []
    group_label_positions = {}
    y_cursor = 0.0
    gap = 0.9
    for group in group_order:
        sub = sample_stats[sample_stats['bio_group'] == group].sort_values('median', ascending=False).copy()
        if sub.empty:
            continue
        group_label_positions[group] = y_cursor - 0.55
        sub['y'] = np.arange(len(sub), dtype=float) + y_cursor
        ordered_frames.append(sub)
        y_cursor = float(sub['y'].max()) + 1.0
        separator_positions.append(y_cursor - 0.5 + gap / 2)
        y_cursor += gap
    ordered = pd.concat(ordered_frames, axis=0)
    ordered['label'] = ordered['sample_id'].str.replace('_all$', '', regex=True)
    fig, ax = plt.subplots(figsize=(4.5, 6.5))
    for _, row in ordered.iterrows():
        left_err = float(row['median'] - row['q1'])
        right_err = float(row['q3'] - row['median'])
        ax.errorbar(
            x=float(row['median']),
            y=float(row['y']),
            xerr=np.array([[left_err], [right_err]]),
            fmt='none',
            capsize=0,
            linewidth=0.8,
            alpha=0.6,
            color=BIO_GROUP_COLORS.get(row['bio_group'], '#B0B0B0'),
        )
        ax.scatter(
            float(row['median']),
            float(row['y']),
            s=22,
            color=BIO_GROUP_COLORS.get(row['bio_group'], '#B0B0B0'),
            zorder=3,
        )
    medians = ordered['median'].to_numpy(dtype=float)
    iqrs = np.maximum((ordered['median'] - ordered['q1']).to_numpy(dtype=float), (ordered['q3'] - ordered['median']).to_numpy(dtype=float))
    all_vals = list(medians + iqrs) + list(medians - iqrs)
    max_abs = max(abs(min(all_vals)), abs(max(all_vals))) * 1.15 if all_vals else 1.0
    ax.set_xlim(-max_abs, max_abs)
    ax.axvline(0, color='#999999', lw=0.6, ls='--')
    for sep in separator_positions[:-1]:
        ax.axhline(sep, color='#E0E0E0', lw=0.8)
    for group, ypos in group_label_positions.items():
        ax.text(-max_abs, ypos, group, color='#777777', fontsize=6.5, va='bottom')
    ax.set_yticks(ordered['y'].to_numpy())
    ax.set_yticklabels(ordered['label'].tolist())
    ax.set_xlabel('TE median (ILR-CLR)')
    ax.set_ylabel('')
    ax.invert_yaxis()
    return save_figure(fig, 'Fig2_TE_sample_median', output_dir)



def plot_fig3(te: pd.DataFrame, anno: pd.DataFrame, output_dir: Path) -> list[Path]:
    # 处理跨物种 NaN 的两步策略：
    # 1. dropna(axis=1, how='all')：删除全为 NaN 的基因列（即某物种在其他物种样本中无表达）
    # 2. melt 后 dropna(subset=['TE'])：删除跨物种造成的残余 NaN 行，避免 violin 分布计算失败
    plot_df = te.dropna(axis=1, how='all').copy()
    plot_df['sample_id'] = te.index
    plot_df = plot_df.merge(anno[['bio_group']], left_on='sample_id', right_index=True, how='left')
    plot_df = plot_df.melt(id_vars=['sample_id', 'bio_group'], var_name='gene', value_name='TE')
    plot_df = plot_df.dropna(subset=['TE'])
    order = [g for g in BIO_GROUP_ORDER if g in plot_df['bio_group'].unique() and g != 'Unknown']
    if 'Unknown' in plot_df['bio_group'].unique():
        order.append('Unknown')
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    violin_signature = inspect.signature(sns.violinplot)
    violin_kwargs = {
        'data': plot_df,
        'x': 'bio_group',
        'y': 'TE',
        'hue': 'bio_group',
        'order': order,
        'hue_order': order,
        'palette': {g: BIO_GROUP_COLORS[g] for g in order},
        'dodge': False,
        'inner': 'quartile',
        'linewidth': 0.5,
        'cut': 0.5,
        'width': 0.75,
        'ax': ax,
    }
    if 'legend' in violin_signature.parameters:
        violin_kwargs['legend'] = False
    if 'density_norm' in violin_signature.parameters:
        violin_kwargs['density_norm'] = 'width'
    else:
        violin_kwargs['scale'] = 'width'
    sns.violinplot(**violin_kwargs)
    if ax.legend_ is not None:
        ax.legend_.remove()
    for i, grp in enumerate(order):
        group_data = plot_df.loc[plot_df['bio_group'] == grp, 'TE']
        if float(group_data.std(ddof=0)) < 0.05:
            ax.scatter([i] * len(group_data), group_data, color=BIO_GROUP_COLORS.get(grp, '#999999'), s=12, alpha=0.7, zorder=5)
    sample_counts = anno['bio_group'].value_counts()
    ax.set_xlabel('')
    ax.set_ylabel('TE (ILR-CLR transformed)')
    ax.set_title('TE distribution by biological group', fontweight='bold')
    ax.axhline(0, color='#999999', lw=0.6, ls='--')
    ax.set_ylim(-5, 7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha='right')
    for i, grp in enumerate(order):
        ax.text(i, -5.3, f'n={int(sample_counts.get(grp, 0))}', ha='center', fontsize=6, clip_on=False)
    return save_figure(fig, 'Fig3_TE_violin', output_dir)



def plot_fig4(te: pd.DataFrame, anno: pd.DataFrame, output_dir: Path) -> list[Path]:
    counts = anno['bio_group'].value_counts().reindex(BIO_GROUP_ORDER).dropna().astype(int)
    order = [group for group in BIO_GROUP_ORDER if group in counts.index]
    total = int(sum(counts[group] for group in order))
    fig, ax = plt.subplots(figsize=(7, 1.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    fig.suptitle(f'Dataset composition (n={total})', fontsize=9, fontweight='bold', y=0.97)
    x = 0.0
    bar_y = 0.38
    bar_h = 0.35
    legend_handles = []
    for grp in order:
        count = int(counts[grp])
        width = count / total if total else 0.0
        rect = plt.Rectangle(
            (x, bar_y),
            width,
            bar_h,
            facecolor=BIO_GROUP_COLORS.get(grp, '#cccccc'),
            edgecolor='white',
            linewidth=0.6,
        )
        ax.add_patch(rect)
        if width > 0.10:
            ax.text(
                x + width / 2,
                bar_y + bar_h / 2,
                f'{grp}\n(n={count})',
                ha='center',
                va='center',
                fontsize=6,
                color='white',
                fontweight='bold',
            )
        else:
            legend_handles.append(
                plt.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=BIO_GROUP_COLORS.get(grp, '#cccccc'),
                    edgecolor='white',
                    linewidth=0.5,
                    label=f'{grp} (n={count})',
                )
            )
        x += width
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc='lower center',
            bbox_to_anchor=(0.5, -0.05),
            ncol=len(legend_handles),
            fontsize=6,
            frameon=False,
            handlelength=1.0,
            handleheight=0.8,
            columnspacing=0.8,
        )
    min_te = float(np.nanmin(te.to_numpy()))
    max_te = float(np.nanmax(te.to_numpy()))
    ax.text(
        0.0,
        0.02,
        f'TE matrix: {total} samples × {te.shape[1]:,} transcripts | Range: {min_te:.1f} to {max_te:.1f}',
        transform=ax.transAxes,
        fontsize=6.5,
        color='#555555',
    )
    return save_figure(fig, 'Fig4_composition', output_dir)



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', default='reports/figures/pub_style')
    parser.add_argument('--copy-dir', default='')
    args = parser.parse_args()
    project_root = find_project_root()
    font_choice = configure_rcparams()
    te, anno = load_data(project_root)
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = project_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = []
    fig1_paths, var, _ = plot_fig1(te, anno, output_dir)
    generated.extend(fig1_paths)
    generated.extend(plot_fig2(te, anno, output_dir))
    generated.extend(plot_fig3(te, anno, output_dir))
    generated.extend(plot_fig4(te, anno, output_dir))
    script_copy_path = output_dir / Path(__file__).name
    if Path(__file__).resolve() != script_copy_path.resolve():
        shutil.copy2(Path(__file__).resolve(), script_copy_path)
    copy_dir = Path(args.copy_dir) if args.copy_dir else None
    if copy_dir is not None:
        copy_dir.mkdir(parents=True, exist_ok=True)
        for path in generated:
            shutil.copy2(path, copy_dir / path.name)
    print('FONT', font_choice)
    print('PCA_VAR', f'{var[0]:.1f}', f'{var[1]:.1f}')
    for path in generated:
        print('WROTE', path)
        if copy_dir is not None:
            print('COPIED', copy_dir / path.name)
    print('SCRIPT', script_copy_path)


if __name__ == '__main__':
    main()
