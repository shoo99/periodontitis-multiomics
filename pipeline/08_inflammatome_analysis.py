"""
Inflammatome Analysis — 치주염 DEG 중 범염증 vs 특이 마커 분리
Cell Reports 2025/2026 (Díaz-Pinés Cort et al.)

분석:
  1. DEG vs Inflammatome overlap 계산
  2. 치주염 특이 마커 분리 (DEG - inflammatome)
  3. ssGSEA 염증 score 계산 → 환자별 중증도 상관
  4. Figure 생성
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

# ─────────────────────────────────────────────
# 0. 설정
# ─────────────────────────────────────────────
RESULT_DIR = Path("../results/inflammatome")
FIG_DIR    = Path("../figures")
REF_DIR    = Path("../references")
RESULT_DIR.mkdir(parents=True, exist_ok=True)

# Inflammatome 유전자 목록 (미리 다운로드된 파일)
INFLAM_TOP100_PATH  = REF_DIR / "inflammatome_top100.txt"
INFLAM_TOP2000_PATH = REF_DIR / "inflammatome_top2000.txt"
INFLAM_GOLD_PATH    = REF_DIR / "inflammation_goldstandard.txt"


# ─────────────────────────────────────────────
# 1. 유전자 목록 로드
# ─────────────────────────────────────────────
def load_inflammatome() -> dict:
    """Inflammatome 유전자 세트 로드"""
    infsets = {}

    if INFLAM_TOP100_PATH.exists():
        infsets['top100'] = set(INFLAM_TOP100_PATH.read_text().strip().split('\n'))
    if INFLAM_TOP2000_PATH.exists():
        infsets['top2000'] = set(INFLAM_TOP2000_PATH.read_text().strip().split('\n'))
    if INFLAM_GOLD_PATH.exists():
        infsets['gold_standard'] = set(INFLAM_GOLD_PATH.read_text().strip().split('\n'))

    print("Inflammatome 로드:")
    for k, v in infsets.items():
        print(f"  {k}: {len(v)}개 유전자")

    return infsets


# ─────────────────────────────────────────────
# 2. DEG vs Inflammatome Overlap 분석
# ─────────────────────────────────────────────
def analyze_deg_inflammatome_overlap(
    deg_df: pd.DataFrame,          # DESeq2 결과 (index=gene_name)
    inflam_sets: dict,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0
) -> dict:
    """
    DEG를 3가지로 분류:
      1. Concordant inflammatome: DEG ∩ inflammatome (범염증)
      2. Periodontitis-specific: DEG - inflammatome (치주염 특이 ⭐)
      3. Non-DEG inflammatome: inflammatome - DEG

    Returns:
        분류 결과 dict
    """
    # 유의 DEG
    sig_up   = deg_df[(deg_df['padj'] < padj_thresh) & (deg_df['log2FoldChange'] > lfc_thresh)].index
    sig_down = deg_df[(deg_df['padj'] < padj_thresh) & (deg_df['log2FoldChange'] < -lfc_thresh)].index
    sig_all  = set(sig_up) | set(sig_down)

    results = {}
    for set_name, inflam_genes in inflam_sets.items():
        overlap_up   = set(sig_up)   & inflam_genes
        overlap_down = set(sig_down) & inflam_genes
        specific_up  = set(sig_up)   - inflam_genes   # ⭐ 치주염 특이 상향
        specific_down= set(sig_down) - inflam_genes   # ⭐ 치주염 특이 하향

        results[set_name] = {
            'sig_up':        set(sig_up),
            'sig_down':      set(sig_down),
            'overlap_up':    overlap_up,
            'overlap_down':  overlap_down,
            'specific_up':   specific_up,
            'specific_down': specific_down,
            'overlap_pct':   len(overlap_up|overlap_down) / max(len(sig_all), 1) * 100
        }

        print(f"\n[{set_name}] Inflammatome 분석:")
        print(f"  전체 DEG: {len(sig_all)}개 (상향 {len(sig_up)}, 하향 {len(sig_down)})")
        print(f"  Inflammatome 겹침: {len(overlap_up|overlap_down)}개 ({results[set_name]['overlap_pct']:.1f}%)")
        print(f"  ⭐ 치주염 특이 상향: {len(specific_up)}개")
        print(f"  ⭐ 치주염 특이 하향: {len(specific_down)}개")

    return results


# ─────────────────────────────────────────────
# 3. ssGSEA — 샘플별 염증 Score
# ─────────────────────────────────────────────
def calculate_inflammation_score(
    expr_mat: pd.DataFrame,      # gene × sample
    inflam_genes: set,
    method: str = 'ssgsea'
) -> pd.Series:
    """
    ssGSEA 기반 샘플별 염증 score 계산
    (gsva 패키지 대안: Python 구현)

    method:
      'ssgsea': Barbie et al. 2009 방식
      'mean'  : 단순 평균 (빠른 대안)
    """
    # 교집합 유전자만
    common = list(set(expr_mat.index) & inflam_genes)
    if len(common) < 5:
        print(f"경고: 교집합 유전자 {len(common)}개 (너무 적음)")

    sub_mat = expr_mat.loc[common]

    if method == 'mean':
        return sub_mat.mean(axis=0)

    elif method == 'ssgsea':
        # 각 샘플별 rank-based enrichment score
        scores = {}
        for sample in expr_mat.columns:
            gene_ranks = expr_mat[sample].rank(ascending=True)
            n_genes    = len(gene_ranks)
            in_set     = gene_ranks[common]
            out_set    = gene_ranks.drop(common)

            # Running sum (simplified ssGSEA)
            es_in  = (in_set / n_genes).sum()
            es_out = (out_set / n_genes).sum()
            scores[sample] = es_in - es_out

        return pd.Series(scores)


# ─────────────────────────────────────────────
# 4. 염증 Score vs 임상 변수 상관
# ─────────────────────────────────────────────
def correlate_inflammation_score(
    inflam_score: pd.Series,
    meta_df: pd.DataFrame,           # 샘플 메타데이터
    clinical_vars: list = None
) -> pd.DataFrame:
    """
    샘플별 염증 score와 임상 변수 상관 분석
    예: 치주낭 깊이(PD), 임상부착수준(CAL), BOP, 방사선 골소실량
    """
    from scipy import stats

    if clinical_vars is None:
        clinical_vars = [c for c in meta_df.columns if c != 'group']

    results = []
    for var in clinical_vars:
        if var not in meta_df.columns:
            continue
        valid = meta_df[var].dropna()
        common_idx = inflam_score.index.intersection(valid.index)
        if len(common_idx) < 10:
            continue

        r, p = stats.spearmanr(inflam_score[common_idx], valid[common_idx])
        results.append({'variable': var, 'spearman_r': r, 'p_value': p})

    return pd.DataFrame(results).sort_values('p_value')


# ─────────────────────────────────────────────
# 5. Figure 생성
# ─────────────────────────────────────────────
def plot_venn_deg_inflammatome(overlap_results: dict,
                                set_name: str = 'top100',
                                save_path: str = None):
    """Figure: DEG vs Inflammatome Venn diagram style"""
    try:
        from matplotlib_venn import venn2
    except ImportError:
        print("matplotlib-venn 필요: pip install matplotlib-venn")
        return

    res = overlap_results[set_name]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 상향 Venn
    ax = axes[0]
    plt.sca(ax)
    venn = venn2(
        subsets=(
            len(res['specific_up']),
            len(set()),  # inflammatome only (상향 기준)
            len(res['overlap_up'])
        ),
        set_labels=('Periodontitis DEG (Up)', f'Inflammatome {set_name}'),
        set_colors=('#FF5722', '#2196F3'),
        alpha=0.6,
        ax=ax
    )
    ax.set_title(f'A. Upregulated DEG vs Inflammatome\n'
                 f'Specific={len(res["specific_up"])}, Shared={len(res["overlap_up"])}',
                 fontsize=12, fontweight='bold')

    # 하향 Venn
    ax = axes[1]
    plt.sca(ax)
    venn2(
        subsets=(
            len(res['specific_down']),
            len(set()),
            len(res['overlap_down'])
        ),
        set_labels=('Periodontitis DEG (Down)', f'Inflammatome {set_name}'),
        set_colors=('#4CAF50', '#9C27B0'),
        alpha=0.6,
        ax=ax
    )
    ax.set_title(f'B. Downregulated DEG vs Inflammatome\n'
                 f'Specific={len(res["specific_down"])}, Shared={len(res["overlap_down"])}',
                 fontsize=12, fontweight='bold')

    plt.suptitle('DEG Classification: Inflammatome vs Periodontitis-Specific',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


def plot_inflammation_score(inflam_score: pd.Series,
                              meta_df: pd.DataFrame,
                              save_path: str = None):
    """Figure: 정상 vs 환자 염증 score 분포"""
    from scipy import stats

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: violin plot
    ax = axes[0]
    score_df = pd.DataFrame({
        'score': inflam_score,
        'group': meta_df.loc[inflam_score.index, 'group']
    })

    palette = {'Control': '#4CAF50', 'Periodontitis': '#F44336'}
    sns.violinplot(data=score_df, x='group', y='score',
                   palette=palette, ax=ax, inner='box', alpha=0.8)
    sns.stripplot(data=score_df, x='group', y='score',
                  palette=palette, ax=ax, alpha=0.5, jitter=True, size=4)

    # Mann-Whitney U test
    ctrl  = score_df[score_df['group']=='Control']['score']
    perio = score_df[score_df['group']=='Periodontitis']['score']
    stat, p = stats.mannwhitneyu(ctrl, perio, alternative='two-sided')
    ax.set_title(f'A. Inflammation Score by Group\n(Mann-Whitney p={p:.2e})',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('ssGSEA Inflammation Score', fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    # Panel B: score 분포 히스토그램
    ax = axes[1]
    for group, color in palette.items():
        subset = score_df[score_df['group'] == group]['score']
        ax.hist(subset, bins=15, alpha=0.6, color=color,
                label=f'{group} (n={len(subset)})', edgecolor='white')

    ax.set_xlabel('Inflammation Score', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('B. Score Distribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


def plot_specific_vs_inflammatome_heatmap(
    deg_df: pd.DataFrame,
    expr_mat: pd.DataFrame,
    overlap_results: dict,
    meta_df: pd.DataFrame,
    n_top: int = 30,
    save_path: str = None
):
    """Figure: 치주염 특이 마커 vs 범염증 마커 비교 히트맵"""
    res = overlap_results.get('top100', list(overlap_results.values())[0])

    # 치주염 특이 상향 top genes (|log2FC| 기준)
    specific_up_ranked = deg_df.loc[
        list(res['specific_up']),
        ['log2FoldChange', 'padj']
    ].sort_values('log2FoldChange', ascending=False)

    # Inflammatome 공유 genes
    shared_ranked = deg_df.loc[
        list(res['overlap_up']),
        ['log2FoldChange', 'padj']
    ].sort_values('log2FoldChange', ascending=False)

    n_specific = min(n_top // 2, len(specific_up_ranked))
    n_shared   = min(n_top // 2, len(shared_ranked))

    genes_to_plot = (specific_up_ranked.head(n_specific).index.tolist() +
                     shared_ranked.head(n_shared).index.tolist())

    # 발현 행렬 준비
    sub_expr = expr_mat.loc[
        [g for g in genes_to_plot if g in expr_mat.index]
    ]

    # 샘플 그룹 정렬
    sample_order = (meta_df[meta_df['group']=='Control'].index.tolist() +
                    meta_df[meta_df['group']=='Periodontitis'].index.tolist())
    sub_expr = sub_expr[[s for s in sample_order if s in sub_expr.columns]]

    # row annotation: specific vs shared
    row_colors = ['#F44336' if g in res['specific_up'] else '#2196F3'
                  for g in sub_expr.index]

    # heatmap
    fig, ax = plt.subplots(figsize=(14, max(8, len(sub_expr) * 0.3)))

    # Z-score normalization per gene
    z_mat = sub_expr.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

    col_colors = ['#4CAF50' if meta_df.loc[s, 'group']=='Control'
                  else '#F44336'
                  for s in sub_expr.columns if s in meta_df.index]

    sns.heatmap(
        z_mat, cmap='RdBu_r', center=0,
        vmin=-3, vmax=3,
        yticklabels=True, xticklabels=False,
        ax=ax, cbar_kws={'label': 'Z-score'}
    )

    # 유전자 이름 색상
    for i, (tick, color) in enumerate(zip(ax.get_yticklabels(), row_colors)):
        tick.set_color(color)

    # 컬럼 구분선 (정상 vs 환자)
    n_ctrl = sum(1 for s in sub_expr.columns
                 if s in meta_df.index and meta_df.loc[s,'group']=='Control')
    ax.axvline(x=n_ctrl, color='black', linewidth=2)

    # 범례
    legend_elements = [
        mpatches.Patch(facecolor='#F44336', label='Periodontitis-specific DEG'),
        mpatches.Patch(facecolor='#2196F3', label='Shared with Inflammatome'),
        mpatches.Patch(facecolor='#4CAF50', label='Control samples'),
        mpatches.Patch(facecolor='#F44336', alpha=0.5, label='Periodontitis samples')
    ]
    ax.legend(handles=legend_elements, loc='upper left',
              bbox_to_anchor=(1.15, 1), fontsize=9)

    ax.set_title('Periodontitis-Specific vs Inflammatome-Shared DEGs',
                 fontsize=13, fontweight='bold')
    ax.set_ylabel('Genes', fontsize=11)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


# ─────────────────────────────────────────────
# 6. 메인 실행
# ─────────────────────────────────────────────
def run_inflammatome_pipeline(
    deg_path: str,
    mrna_expr_path: str,
    sample_meta_path: str
):
    print("=" * 60)
    print("Inflammatome Analysis Pipeline")
    print("=" * 60)

    # 데이터 로드
    deg_df  = pd.read_csv(deg_path, index_col=0)
    expr_df = pd.read_csv(mrna_expr_path, index_col=0)
    meta_df = pd.read_csv(sample_meta_path, index_col=0)

    # Inflammatome 로드
    inflam_sets = load_inflammatome()

    # 1. DEG vs Inflammatome overlap
    print("\n[1/4] DEG vs Inflammatome overlap 분석...")
    overlap_results = analyze_deg_inflammatome_overlap(deg_df, inflam_sets)

    # 결과 저장
    for set_name, res in overlap_results.items():
        # 치주염 특이 유전자 저장
        specific_df = deg_df.loc[
            list(res['specific_up'] | res['specific_down'])
        ].copy()
        specific_df['direction'] = specific_df['log2FoldChange'].apply(
            lambda x: 'Up' if x > 0 else 'Down'
        )
        specific_df['classification'] = 'Periodontitis-specific'
        specific_df.to_csv(RESULT_DIR / f'periodontitis_specific_DEG_{set_name}.csv')
        print(f"  저장: periodontitis_specific_DEG_{set_name}.csv")

    # 2. ssGSEA 염증 score
    print("\n[2/4] ssGSEA 염증 score 계산...")
    inflam_score = calculate_inflammation_score(
        expr_df,
        inflam_sets.get('top100', inflam_sets.get('gold_standard', set())),
        method='ssgsea'
    )

    # 3. Figure 생성
    print("\n[3/4] Figure 생성...")

    # Venn diagram
    if 'top100' in overlap_results:
        plot_venn_deg_inflammatome(
            overlap_results, 'top100',
            save_path=str(FIG_DIR / 'FigS_inflammatome_venn.png')
        )

    # 염증 score 분포
    plot_inflammation_score(
        inflam_score, meta_df,
        save_path=str(FIG_DIR / 'FigS_inflammation_score.png')
    )

    # 비교 히트맵
    plot_specific_vs_inflammatome_heatmap(
        deg_df, expr_df, overlap_results, meta_df,
        n_top=30,
        save_path=str(FIG_DIR / 'FigS_specific_vs_inflammatome_heatmap.png')
    )

    # 4. 임상 변수와 상관
    print("\n[4/4] 염증 score vs 임상 변수 상관...")
    # 임상 변수 있으면 활성화
    # corr_df = correlate_inflammation_score(inflam_score, meta_df)
    # corr_df.to_csv(RESULT_DIR / 'inflammation_score_clinical_corr.csv', index=False)

    # score 저장
    pd.DataFrame({
        'sample': inflam_score.index,
        'inflammation_score': inflam_score.values,
        'group': meta_df.loc[inflam_score.index, 'group']
    }).to_csv(RESULT_DIR / 'ssgsea_inflammation_scores.csv', index=False)

    print("\n✅ Inflammatome Analysis 완료!")
    return overlap_results, inflam_score


# 실행
if __name__ == "__main__":
    run_inflammatome_pipeline(
        deg_path        = "../results/mrna/DESeq2_results.csv",
        mrna_expr_path  = "../data/processed/mrna_vst.csv",
        sample_meta_path= "../data/sample_metadata.csv"
    )
