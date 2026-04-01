"""
PathIntegrate — Pathway-based Multi-omics Integration
PLoS Comput Biol 2024 (Wieder et al.)

분자 수준이 아닌 Pathway 수준에서 통합 → 더 해석 용이한 결과

설치: pip install PathIntegrate

분석:
  1. 분자 → 경로 점수 변환 (ssGSEA/kPCA)
  2. Single-view 모델 (PLS-DA on pathway scores)
  3. Multi-view 모델 (각 오믹스별 경로 공간 → 통합)
  4. 경로별 중요도 + 오믹스별 기여
  5. KEGG/Reactome pathway 시각화
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

RESULT_DIR = Path("../results/pathintegrate")
FIG_DIR    = Path("../figures")
RESULT_DIR.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────
# 1. PathIntegrate 설치 확인 + 실행
# ─────────────────────────────────────────────
def run_pathintegrate(
    mrna_expr: pd.DataFrame,    # gene × sample (VST or log2)
    prot_expr: pd.DataFrame,    # protein × sample
    metab_expr: pd.DataFrame,   # metabolite × sample
    y: pd.Series,               # 0=Control, 1=Periodontitis
    database: str = "KEGG"      # "KEGG" or "Reactome"
):
    """
    PathIntegrate Multi-View 모델 실행
    """
    try:
        import pathintegrate as pi
    except ImportError:
        print("설치 필요: pip install PathIntegrate sspa")
        return None

    # 1. pathway 데이터베이스 로드
    if database == "KEGG":
        pathways = pi.load_pathways(organism="hsa", database="KEGG")
    else:
        pathways = pi.load_pathways(organism="HSA", database="Reactome")

    print(f"경로 DB: {database}, 총 {len(pathways)}개 경로")

    # 2. 각 오믹스별 단일 뷰 모델
    print("\n[Single-View] mRNA...")
    sv_mrna = pi.SingleViewPathIntegrate(
        pathway_source=pathways,
        scoring_method="sspa_kpca",   # kPCA 기반 pathway score
        model=None                     # PLS-DA 자동 사용
    )
    sv_mrna.fit(mrna_expr.T, y)

    print("[Single-View] Proteomics...")
    sv_prot = pi.SingleViewPathIntegrate(
        pathway_source=pathways,
        scoring_method="sspa_kpca"
    )
    sv_prot.fit(prot_expr.T, y)

    print("[Single-View] Metabolomics...")
    sv_metab = pi.SingleViewPathIntegrate(
        pathway_source=pathways,
        scoring_method="sspa_kpca"
    )
    sv_metab.fit(metab_expr.T, y)

    # 3. Multi-View 모델 (3개 오믹스 동시)
    print("\n[Multi-View] 통합 모델...")
    mv_model = pi.MultiViewPathIntegrate(
        pathway_source=pathways,
        scoring_method="sspa_kpca"
    )
    mv_model.fit(
        [mrna_expr.T, prot_expr.T, metab_expr.T],
        y,
        omics_names=["mRNA", "Protein", "Metabolite"]
    )

    return {
        'sv_mrna': sv_mrna,
        'sv_prot': sv_prot,
        'sv_metab': sv_metab,
        'mv': mv_model
    }


# ─────────────────────────────────────────────
# 2. 수동 Pathway Score 계산 (PathIntegrate 없을 때 대안)
# ─────────────────────────────────────────────
def manual_pathway_score_analysis(
    mrna_expr: pd.DataFrame,
    prot_expr: pd.DataFrame,
    metab_expr: pd.DataFrame,
    y: pd.Series,
    top_n_pathway: int = 20
):
    """
    PathIntegrate 없을 때 수동 구현:
    gseapy ssGSEA + KEGG pathway gene sets
    """
    import gseapy as gp
    from scipy import stats

    print("gseapy 기반 수동 pathway 분석...")

    # KEGG pathway gene sets (gseapy 내장)
    kegg_sets = gp.get_library("KEGG_2021_Human")
    reactome_sets = gp.get_library("Reactome_2022")

    results = {}

    # mRNA pathway score (ssGSEA)
    print("  mRNA ssGSEA...")
    mrna_ssgsea = gp.ssgsea(
        data=mrna_expr,
        gene_sets=kegg_sets,
        min_size=5,
        max_size=500,
        no_plot=True
    )
    mrna_path_scores = mrna_ssgsea.res2d.pivot(
        index='Term', columns='Name', values='NES'
    ).fillna(0)
    results['mrna_pathway_scores'] = mrna_path_scores

    # 단백질은 gene symbol로 매핑 후 ssGSEA
    print("  Proteomics ssGSEA...")
    prot_ssgsea = gp.ssgsea(
        data=prot_expr,
        gene_sets=kegg_sets,
        min_size=3,  # 단백질은 커버리지 낮으므로 낮게
        max_size=500,
        no_plot=True
    )
    prot_path_scores = prot_ssgsea.res2d.pivot(
        index='Term', columns='Name', values='NES'
    ).fillna(0)
    results['prot_pathway_scores'] = prot_path_scores

    # 공통 경로 찾기
    common_paths = list(
        set(mrna_path_scores.index) & set(prot_path_scores.index)
    )
    print(f"  공통 경로: {len(common_paths)}개")

    # 각 경로별 통계 (Control vs Periodontitis)
    pathway_stats = []
    for path in common_paths:
        if path not in mrna_path_scores.index:
            continue
        scores = mrna_path_scores.loc[path]
        valid = scores.dropna()
        common_samples = valid.index.intersection(y.index)
        if len(common_samples) < 10:
            continue

        ctrl  = valid[common_samples][y[common_samples]==0]
        perio = valid[common_samples][y[common_samples]==1]

        if len(ctrl) < 3 or len(perio) < 3:
            continue

        stat, p = stats.mannwhitneyu(ctrl, perio, alternative='two-sided')
        fc = perio.mean() - ctrl.mean()

        pathway_stats.append({
            'pathway': path,
            'mean_diff': fc,
            'p_value': p,
            'ctrl_mean': ctrl.mean(),
            'perio_mean': perio.mean()
        })

    from statsmodels.stats.multitest import multipletests
    path_df = pd.DataFrame(pathway_stats)
    if len(path_df) > 0:
        path_df['padj'] = multipletests(path_df['p_value'], method='fdr_bh')[1]
        path_df = path_df.sort_values('padj')

    results['pathway_stats'] = path_df
    results['common_paths'] = common_paths

    return results


# ─────────────────────────────────────────────
# 3. Joint Pathway Enrichment — mRNA + Metabolomics
# ─────────────────────────────────────────────
def joint_pathway_enrichment(
    deg_list: list,            # 유의 유전자
    diff_metabolites: list,    # 유의 대사체 (HMDB ID or name)
    method: str = "metaboanalyst"
):
    """
    mRNA + 대사체 동시 Pathway 농축 분석
    (MetaboAnalyst 방식: Over-Representation Analysis)

    방법:
      1. 유전자 → KEGG gene ID 변환
      2. 대사체 → KEGG compound ID 변환
      3. 두 리스트 동시 KEGG pathway ORA
      4. 교집합이 있는 경로 → 진정한 multi-omics pathway

    주의:
      MetaboAnalyst 6.0의 Joint Pathway Analysis (웹) 또는
      MetaboAnalystR 패키지 사용 권장
    """
    print("Joint Pathway Enrichment:")
    print(f"  입력: 유전자 {len(deg_list)}개, 대사체 {len(diff_metabolites)}개")
    print()
    print("  권장 도구:")
    print("  1. MetaboAnalyst 6.0 웹 (https://www.metaboanalyst.ca/)")
    print("     → Joint Pathway Analysis 모듈")
    print("     → 입력: 유전자 목록 + 대사체 목록")
    print("     → 경로: KEGG / HumanCyc / SMPDB")
    print()
    print("  2. MetaboAnalystR (R 패키지):")
    print("     library(MetaboAnalystR)")
    print("     mSet <- InitDataObjects('conc', 'pathinteg', FALSE)")
    print("     mSet <- Read.TextData(mSet, 'genes.txt', 'rowu', 'disc')")
    print("     mSet <- SetKEGG.PathLib(mSet, 'hsa', 'current')")
    print("     mSet <- PerformIntegPathwayAnalysis(mSet, 'dc', 'hyper', 'global')")
    print()
    print("  3. gseapy (Python 근사):")

    # gseapy ORA 예시
    try:
        import gseapy as gp

        # 유전자 ORA
        gene_enr = gp.enrichr(
            gene_list=deg_list,
            gene_sets=["KEGG_2021_Human", "Reactome_2022",
                       "WikiPathway_2021_Human"],
            organism="Human",
            cutoff=0.05
        )

        print(f"\n  gseapy ORA 완료:")
        if gene_enr.res2d is not None:
            sig_paths = gene_enr.res2d[gene_enr.res2d['Adjusted P-value'] < 0.05]
            print(f"  유의 경로: {len(sig_paths)}개")

        return gene_enr

    except Exception as e:
        print(f"  gseapy 실행 오류: {e}")
        return None


# ─────────────────────────────────────────────
# 4. Pathway 통합 Figure
# ─────────────────────────────────────────────
def plot_pathway_comparison(
    mrna_path_df: pd.DataFrame,
    prot_path_df: pd.DataFrame,
    top_n: int = 15,
    save_path: str = None
):
    """
    Figure: mRNA vs Protein pathway score 비교 버블 차트
    (공통으로 유의한 경로 강조)
    """
    # 공통 유의 경로
    sig_mrna  = set(mrna_path_df[mrna_path_df['padj'] < 0.05]['pathway'])
    sig_prot  = set(prot_path_df[prot_path_df['padj'] < 0.05]['pathway'])
    common    = sig_mrna & sig_prot
    mrna_only = sig_mrna - sig_prot
    prot_only = sig_prot - sig_mrna

    print(f"\nPathway 요약:")
    print(f"  mRNA 유의: {len(sig_mrna)}개")
    print(f"  Protein 유의: {len(sig_prot)}개")
    print(f"  ⭐ 공통 유의: {len(common)}개")

    if not common:
        print("  공통 경로 없음 — threshold 조정 필요")
        return

    # 공통 경로 데이터 병합
    mrna_sub = mrna_path_df[mrna_path_df['pathway'].isin(common)].set_index('pathway')
    prot_sub  = prot_path_df[prot_path_df['pathway'].isin(common)].set_index('pathway')

    plot_df = pd.DataFrame({
        'pathway': list(common),
        'mrna_fc':   [mrna_sub.loc[p, 'mean_diff'] if p in mrna_sub.index else 0 for p in common],
        'prot_fc':   [prot_sub.loc[p, 'mean_diff'] if p in prot_sub.index else 0 for p in common],
        'mrna_padj': [mrna_sub.loc[p, 'padj'] if p in mrna_sub.index else 1 for p in common],
        'prot_padj': [prot_sub.loc[p, 'padj'] if p in prot_sub.index else 1 for p in common],
    })
    plot_df['-log10_p'] = (-np.log10(plot_df['mrna_padj'])).clip(0, 10)
    plot_df = plot_df.sort_values('-log10_p', ascending=False).head(top_n)

    # 경로명 축약 (50자 초과 시)
    plot_df['pathway_short'] = plot_df['pathway'].apply(
        lambda x: x[:47] + '...' if len(x) > 50 else x
    )

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.4)))

    scatter = ax.scatter(
        plot_df['mrna_fc'],
        range(len(plot_df)),
        c=plot_df['prot_fc'],
        s=plot_df['-log10_p'] * 30,
        cmap='RdBu_r', vmin=-2, vmax=2,
        alpha=0.8, edgecolors='black', linewidths=0.5
    )

    plt.colorbar(scatter, ax=ax, label='Protein FC (Periodontitis vs Control)')

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['pathway_short'], fontsize=9)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('mRNA Pathway Score Difference\n(Periodontitis - Control)', fontsize=11)
    ax.set_title('Shared Significant Pathways\n(mRNA + Protein, common)',
                 fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


# ─────────────────────────────────────────────
# 5. MetaboAnalyst 6.0 활용 가이드 출력
# ─────────────────────────────────────────────
def print_metaboanalyst_guide():
    guide = """
    ╔══════════════════════════════════════════════════════════════╗
    ║        MetaboAnalyst 6.0 Joint Pathway 분석 가이드           ║
    ╠══════════════════════════════════════════════════════════════╣
    ║                                                              ║
    ║  1. 접속: https://www.metaboanalyst.ca/                      ║
    ║  2. 메뉴: Multi-omics → Joint Pathway Analysis               ║
    ║                                                              ║
    ║  입력 A (Genes):                                             ║
    ║    - Differentially Expressed Genes (gene symbol)            ║
    ║    - DESeq2 padj < 0.05, |log2FC| > 1.0                      ║
    ║                                                              ║
    ║  입력 B (Metabolites):                                       ║
    ║    - Differential Metabolites (HMDB ID 권장)                  ║
    ║    - VIP > 1.0, p < 0.05                                      ║
    ║                                                              ║
    ║  설정:                                                       ║
    ║    - Organism: Homo sapiens                                   ║
    ║    - Pathway DB: KEGG (권장) or SMPDB                        ║
    ║    - Enrichment: Hypergeometric test                          ║
    ║    - Topology: Degree centrality                              ║
    ║                                                              ║
    ║  해석 기준:                                                   ║
    ║    - p.combine < 0.05 → 유의한 multi-omics pathway           ║
    ║    - Impact > 0.1 → 위상학적으로 중요한 경로                  ║
    ║    - 버블 크기 = Impact, 색 = -log10(p)                      ║
    ║                                                              ║
    ║  치주염에서 예상 TOP 경로:                                     ║
    ║    - Purine metabolism (퓨린 대사) ← 대사체+유전자 교집합      ║
    ║    - Arachidonic acid metabolism (염증 지질)                  ║
    ║    - Amino acid metabolism (아미노산)                         ║
    ║    - NF-kB signaling (면역 핵심)                              ║
    ║    - Cytokine-cytokine receptor interaction                   ║
    ║                                                              ║
    ╚══════════════════════════════════════════════════════════════╝
    """
    print(guide)


# 메인 실행
if __name__ == "__main__":
    print_metaboanalyst_guide()

    # 수동 분석 실행 (PathIntegrate 또는 gseapy 기반)
    # run_pathintegrate(mrna_expr, prot_expr, metab_expr, y)
    # manual_pathway_score_analysis(mrna_expr, prot_expr, metab_expr, y)
