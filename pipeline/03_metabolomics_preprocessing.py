"""
03_metabolomics_preprocessing.py
Untargeted Metabolomics 전처리 + OPLS-DA + 경로 분석
입력: metabolite × sample intensity table
출력: 유의 대사체, OPLS-DA 결과, 경로 분석
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

INDIR   = Path("../data/processed")
OUTDIR  = Path("../results/metabolomics")
FIGDIR  = Path("../figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("03. Metabolomics Preprocessing & Analysis")
print("=" * 60)


# ══════════════════════════════════════════════════════════════
# STEP 1: 데이터 로드 + QC
# ══════════════════════════════════════════════════════════════
print("\n[1/5] 데이터 로드...")

metab = pd.read_csv(INDIR / "metabolomics_intensity.csv", index_col=0)
meta  = pd.read_csv(INDIR / "sample_metadata.csv",  index_col=0)
meta  = meta.loc[metab.columns]
y     = (meta['group'] == 'Periodontitis').astype(int)

print(f"  대사체: {metab.shape[0]}개, 샘플: {metab.shape[1]}개")

# 0/음수 → NaN
metab[metab <= 0] = np.nan

# 결측치 필터 (그룹별 50% 이상 결측 → 제거)
ctrl_cols  = meta[meta['group']=='Control'].index.tolist()
perio_cols = meta[meta['group']=='Periodontitis'].index.tolist()

keep_ctrl  = metab[ctrl_cols].isna().mean(axis=1) < 0.5
keep_perio = metab[perio_cols].isna().mean(axis=1) < 0.5
metab_filt = metab[keep_ctrl | keep_perio].copy()
print(f"  필터 후 대사체: {metab_filt.shape[0]}개")


# ══════════════════════════════════════════════════════════════
# STEP 2: 정규화 + 스케일링
# ══════════════════════════════════════════════════════════════
print("\n[2/5] 정규화 + 스케일링...")

# KNN imputation
from sklearn.impute import KNNImputer
imputer   = KNNImputer(n_neighbors=5)
metab_imp = pd.DataFrame(
    imputer.fit_transform(metab_filt.T).T,
    index=metab_filt.index,
    columns=metab_filt.columns
)

# Log2 변환
metab_log2 = np.log2(metab_imp + 1)

# Pareto scaling (PCA/PLS-DA 전 권장)
def pareto_scale(df):
    mean_ = df.mean(axis=1)
    std_  = df.std(axis=1)
    sqrt_std = np.sqrt(std_)
    return df.subtract(mean_, axis=0).divide(sqrt_std.replace(0, 1), axis=0)

metab_pareto = pareto_scale(metab_log2)
metab_pareto.to_csv(OUTDIR / "metabolomics_log2_pareto.csv")
metab_log2.to_csv(OUTDIR / "metabolomics_log2.csv")

print(f"  Log2 + Pareto scaling 완료")

# PCA — 배치/이상치 확인
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pca    = PCA(n_components=2)
X_pca  = pca.fit_transform(metab_pareto.T)
pca_df = pd.DataFrame({
    'PC1': X_pca[:,0], 'PC2': X_pca[:,1],
    'group': meta['group'], 'sample': meta.index
})
pca_var = pca.explained_variance_ratio_ * 100

fig, ax = plt.subplots(figsize=(8,6))
for grp, color in [('Control','#4CAF50'),('Periodontitis','#E91E63')]:
    sub = pca_df[pca_df['group']==grp]
    ax.scatter(sub['PC1'], sub['PC2'], c=color, s=80,
               alpha=0.8, label=grp, edgecolors='white')
for _, row in pca_df.iterrows():
    ax.annotate(row['sample'], (row['PC1'], row['PC2']),
                fontsize=6, alpha=0.7)
ax.set_xlabel(f"PC1 ({pca_var[0]:.1f}%)", fontsize=12)
ax.set_ylabel(f"PC2 ({pca_var[1]:.1f}%)", fontsize=12)
ax.set_title("Metabolomics PCA", fontsize=13, fontweight='bold')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(FIGDIR / "Fig_Metab_PCA.pdf", dpi=300, bbox_inches='tight')
plt.close()


# ══════════════════════════════════════════════════════════════
# STEP 3: OPLS-DA + Permutation Test
# ══════════════════════════════════════════════════════════════
print("\n[3/5] OPLS-DA + VIP 계산...")

try:
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.model_selection import cross_val_score, LeaveOneOut

    X_plsda = metab_pareto.T.values
    y_vals  = y.values

    # PLS-DA (OPLS-DA 근사 — n=60 소규모)
    pls = PLSRegression(n_components=2, scale=True, max_iter=1000)
    pls.fit(X_plsda, y_vals)

    # LOOCV R²/Q² 계산
    loo = LeaveOneOut()
    scores_q2 = cross_val_score(pls, X_plsda, y_vals,
                                  cv=loo, scoring='r2')
    Q2 = scores_q2.mean()
    R2 = pls.score(X_plsda, y_vals)

    print(f"  R² = {R2:.3f}, Q² = {Q2:.3f}")
    if Q2 <= 0:
        print("  ⚠️ Q² ≤ 0: 과적합 가능성. 피처 수 줄이기 권장")

    # Permutation test (200회)
    print("  Permutation test 200회...")
    perm_q2_scores = []
    for _ in range(200):
        y_perm  = np.random.permutation(y_vals)
        pls_p   = PLSRegression(n_components=2, scale=True, max_iter=500)
        perm_q2 = cross_val_score(pls_p, X_plsda, y_perm,
                                   cv=5, scoring='r2').mean()
        perm_q2_scores.append(perm_q2)

    p_perm = sum(np.array(perm_q2_scores) >= Q2) / len(perm_q2_scores)
    print(f"  Permutation p = {p_perm:.3f} (< 0.05 필요)")

    # VIP score 계산
    def calc_vip(pls_model, X):
        T = pls_model.x_scores_         # scores (n × components)
        W = pls_model.x_weights_         # weights (features × components)
        Q = pls_model.y_loadings_        # y loadings

        p    = W.shape[0]   # n_features
        nc   = W.shape[1]   # n_components
        SSY  = np.sum(T**2 * Q**2, axis=0)  # SS explained per component

        VIP = np.zeros(p)
        for j in range(p):
            WW = (W[j,:] / np.linalg.norm(W[:,k]) for k in range(nc))
            VIP[j] = np.sqrt(p * np.sum(
                [SSY[k] * (W[j,k] / np.linalg.norm(W[:,k]))**2
                 for k in range(nc)]
            ) / np.sum(SSY))
        return VIP

    vip_scores = calc_vip(pls, X_plsda)
    vip_df     = pd.DataFrame({
        'metabolite': metab_pareto.index,
        'VIP':        vip_scores
    }).sort_values('VIP', ascending=False)

    print(f"  VIP > 1.0: {(vip_df['VIP'] > 1.0).sum()}개")

    # PLS-DA score plot
    T = pls.x_scores_
    score_df = pd.DataFrame({
        'LV1': T[:,0], 'LV2': T[:,1],
        'group': meta['group'].values
    })

    fig, axes = plt.subplots(1,2, figsize=(14,6))

    # Score plot
    ax = axes[0]
    for grp, color in [('Control','#4CAF50'),('Periodontitis','#E91E63')]:
        sub = score_df[score_df['group']==grp]
        ax.scatter(sub['LV1'], sub['LV2'], c=color, s=80,
                   alpha=0.8, label=grp, edgecolors='white')
    ax.set_xlabel("LV1", fontsize=12); ax.set_ylabel("LV2", fontsize=12)
    ax.set_title(f"PLS-DA Scores\nR²={R2:.3f}, Q²={Q2:.3f}", fontsize=12)
    ax.legend(); ax.grid(True, alpha=0.3)
    ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)

    # VIP top 20
    ax = axes[1]
    top_vip = vip_df.head(20)
    ax.barh(range(len(top_vip)), top_vip['VIP'].values[::-1],
            color=['#F44336' if v>1 else '#9E9E9E' for v in top_vip['VIP'].values[::-1]],
            alpha=0.8)
    ax.set_yticks(range(len(top_vip)))
    ax.set_yticklabels(top_vip['metabolite'].values[::-1], fontsize=8)
    ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel("VIP Score", fontsize=11)
    ax.set_title("Top 20 VIP (PLS-DA)", fontsize=12)
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(FIGDIR / "Fig_Metab_PLSDA.pdf", dpi=300, bbox_inches='tight')
    plt.close()

except Exception as e:
    print(f"  PLS-DA 오류: {e}")
    vip_df = pd.DataFrame({'metabolite': metab_pareto.index, 'VIP': np.ones(len(metab_pareto))})
    Q2 = 0
    p_perm = 1.0


# ══════════════════════════════════════════════════════════════
# STEP 4: 단변량 통계 + 최종 유의 대사체 선별
# ══════════════════════════════════════════════════════════════
print("\n[4/5] 단변량 분석 + 유의 대사체 선별...")

# t-test + Fold Change
stat_results = []
for metab_name in metab_log2.index:
    ctrl_vals  = metab_log2.loc[metab_name, ctrl_cols].dropna()
    perio_vals = metab_log2.loc[metab_name, perio_cols].dropna()

    if len(ctrl_vals) < 3 or len(perio_vals) < 3:
        continue

    # Levene test (등분산 확인)
    _, p_lev = stats.levene(ctrl_vals, perio_vals)
    equal_var = p_lev > 0.05

    t_stat, p_val = stats.ttest_ind(perio_vals, ctrl_vals,
                                     equal_var=equal_var)
    fc = perio_vals.mean() - ctrl_vals.mean()  # log2 scale
    fc_ratio = 2**fc  # actual fold change

    stat_results.append({
        'metabolite': metab_name,
        'mean_ctrl':  ctrl_vals.mean(),
        'mean_perio': perio_vals.mean(),
        'log2FC':     fc,
        'FC':         fc_ratio,
        'pvalue':     p_val,
        'equal_var':  equal_var
    })

stat_df = pd.DataFrame(stat_results)
stat_df['padj'] = multipletests(stat_df['pvalue'], method='fdr_bh')[1]
stat_df = stat_df.merge(vip_df, on='metabolite', how='left')
stat_df['VIP'] = stat_df['VIP'].fillna(0)
stat_df = stat_df.sort_values('pvalue')

# 3중 필터: VIP > 1.0 AND padj < 0.05 AND |FC| > 1.5
sig_metab = stat_df[
    (stat_df['VIP']  > 1.0) &
    (stat_df['padj'] < 0.05) &
    (stat_df['FC'].abs() > 1.5)
].copy()
sig_metab['direction'] = sig_metab['log2FC'].apply(
    lambda x: 'Up' if x > 0 else 'Down'
)

print(f"  유의 대사체 (3중 필터): {len(sig_metab)}개")
print(f"    상향: {(sig_metab['direction']=='Up').sum()}, "
      f"하향: {(sig_metab['direction']=='Down').sum()}")

stat_df.to_csv(OUTDIR / "metabolomics_all_results.csv", index=False)
sig_metab.to_csv(OUTDIR / "metabolomics_significant.csv", index=False)

# Volcano
fig, ax = plt.subplots(figsize=(8,7))
colors  = stat_df['padj'].apply(
    lambda p: '#F44336' if p < 0.05 else '#9E9E9E')
ax.scatter(stat_df['log2FC'], -np.log10(stat_df['padj']+1e-10),
           c=colors, alpha=0.5, s=20)

# VIP > 1 표시
vip_sig = stat_df[(stat_df['VIP']>1) & (stat_df['padj']<0.05)]
ax.scatter(vip_sig['log2FC'], -np.log10(vip_sig['padj']+1e-10),
           c='#9C27B0', alpha=0.8, s=50, label='VIP>1 & sig')

# Top 15 레이블
for _, row in sig_metab.head(15).iterrows():
    ax.annotate(row['metabolite'][:25],
                (row['log2FC'], -np.log10(row['padj']+1e-10)),
                fontsize=6, alpha=0.8,
                xytext=(3,3), textcoords='offset points')

ax.axvline(x=0, color='gray', lw=0.5, linestyle='--')
ax.axhline(y=-np.log10(0.05), color='red', lw=0.8, linestyle='--')
ax.set_xlabel("log2 Fold Change", fontsize=12)
ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
ax.set_title("Volcano — Metabolomics", fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(FIGDIR / "Fig_Metab_Volcano.pdf", dpi=300, bbox_inches='tight')
plt.close()

# Boxplot (Top 12 대사체)
top_metabs = sig_metab.head(12)['metabolite'].tolist()
if top_metabs:
    fig, axes = plt.subplots(3, 4, figsize=(16, 12))
    axes = axes.flatten()
    for i, m in enumerate(top_metabs):
        ax = axes[i]
        plot_df = pd.DataFrame({
            'value': metab_log2.loc[m].values,
            'group': meta['group'].values
        })
        sns.boxplot(data=plot_df, x='group', y='value',
                    palette={'Control':'#4CAF50','Periodontitis':'#E91E63'},
                    ax=ax, width=0.5)
        sns.stripplot(data=plot_df, x='group', y='value',
                      color='black', ax=ax, alpha=0.4, size=4)
        ax.set_title(m[:30], fontsize=9, fontweight='bold')
        ax.set_xlabel(""); ax.set_ylabel("log2 intensity", fontsize=8)
        row = sig_metab[sig_metab['metabolite']==m].iloc[0]
        ax.text(0.5, 0.97, f"FC={row['FC']:.2f}, p={row['padj']:.3f}",
                transform=ax.transAxes, ha='center', va='top', fontsize=7)
    for j in range(len(top_metabs), len(axes)):
        axes[j].axis('off')
    plt.suptitle("Top Significant Metabolites", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGDIR / "Fig_Metab_Boxplots.pdf", dpi=300, bbox_inches='tight')
    plt.close()


# ══════════════════════════════════════════════════════════════
# STEP 5: gseapy 경로 분석 (대사체 → HMDB → KEGG)
# ══════════════════════════════════════════════════════════════
print("\n[5/5] 대사체 경로 분석...")

try:
    import gseapy as gp

    # 상향 대사체 경로 분석
    up_metabs = sig_metab[sig_metab['direction']=='Up']['metabolite'].tolist()

    if len(up_metabs) >= 5:
        enr = gp.enrichr(
            gene_list=up_metabs,
            gene_sets=["HMDB_Metabolites", "KEGG_2021_Human"],
            organism="Human",
            cutoff=0.05,
            outdir=str(OUTDIR / "enrichr_metabolites")
        )
        if enr.res2d is not None:
            sig_path = enr.res2d[enr.res2d['Adjusted P-value'] < 0.05]
            sig_path.to_csv(OUTDIR / "metabolite_pathway_enrichment.csv", index=False)
            print(f"  유의 경로: {len(sig_path)}개")

except Exception as e:
    print(f"  gseapy 경로 분석 오류: {e}")

# 요약 저장
summary = {
    'n_metabolites_input':    metab.shape[0],
    'n_metabolites_filtered': metab_filt.shape[0],
    'n_metabolites_sig':      len(sig_metab),
    'n_sig_up':               (sig_metab['direction']=='Up').sum(),
    'n_sig_down':             (sig_metab['direction']=='Down').sum(),
    'PLS_R2':                 round(R2, 3),
    'PLS_Q2':                 round(Q2, 3),
    'permutation_p':          round(p_perm, 3)
}
pd.DataFrame([summary]).to_csv(OUTDIR / "metabolomics_summary.csv", index=False)

print("\n✅ Metabolomics Analysis 완료!")
print(f"  유의 대사체: {len(sig_metab)}개 | R²={R2:.3f} | Q²={Q2:.3f}")
