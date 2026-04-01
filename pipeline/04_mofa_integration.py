"""
04_mofa_integration.py
MOFA+ 비지도 멀티오믹스 통합
3개 오믹스 레이어 → Latent Factor → 치주염 연관 Factor 발굴
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

OUTDIR = Path("../results/mofa")
FIGDIR = Path("../figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("04. MOFA+ Multi-Omics Integration")
print("=" * 60)


# ══════════════════════════════════════════════════════════════
# STEP 1: 데이터 준비
# ══════════════════════════════════════════════════════════════
print("\n[1/4] 데이터 준비...")

meta       = pd.read_csv("../data/processed/sample_metadata.csv", index_col=0)
mrna_vst   = pd.read_csv("../results/mrna/mrna_vst_matrix.csv", index_col=0)
prot_log2  = pd.read_csv("../results/proteomics/proteomics_log2_normalized.csv", index_col=0)
metab_log2 = pd.read_csv("../results/metabolomics/metabolomics_log2.csv", index_col=0)

# 공통 샘플 확인
common_samples = list(
    set(mrna_vst.columns) & set(prot_log2.columns) & set(metab_log2.columns)
)
common_samples = [s for s in meta.index if s in common_samples]
print(f"  공통 샘플: {len(common_samples)}")

# 각 오믹스 상위 변동 피처만 (MOFA 속도 최적화)
def top_variable_features(df, n=2000):
    mad = df.mad(axis=1) if hasattr(df, 'mad') else df.apply(lambda x: np.median(np.abs(x - x.median())), axis=1)
    return df.loc[mad.nlargest(n).index]

mrna_sub  = top_variable_features(mrna_vst[common_samples], n=2000)
prot_sub  = top_variable_features(prot_log2[common_samples], n=1000)
metab_sub = metab_log2[common_samples]

print(f"  mRNA: {mrna_sub.shape[0]}개")
print(f"  Protein: {prot_sub.shape[0]}개")
print(f"  Metabolite: {metab_sub.shape[0]}개")

y = (meta.loc[common_samples, 'group'] == 'Periodontitis').astype(int)


# ══════════════════════════════════════════════════════════════
# STEP 2: MOFA+ 훈련
# ══════════════════════════════════════════════════════════════
print("\n[2/4] MOFA+ 훈련...")

try:
    from mofapy2.run.entry_point import entry_point

    ent = entry_point()

    # 데이터 설정 (feature × sample → MOFA는 sample × feature로 변환)
    ent.set_data_options(
        scale_groups=False,
        scale_views=False      # 각 view 독립 스케일 유지
    )

    # 데이터 입력 (리스트 형식)
    ent.set_data_matrix(
        data=[
            [mrna_sub.T.values],    # 1 group × 1 view (mRNA)
            [prot_sub.T.values],    # 1 group × 1 view (Protein)
            [metab_sub.T.values]    # 1 group × 1 view (Metabolite)
        ],
        likelihoods=["gaussian", "gaussian", "gaussian"],
        views_names=["mRNA", "Protein", "Metabolite"],
        groups_names=["Periodontitis_Study"],
        samples_names=[common_samples],
        features_names=[
            mrna_sub.index.tolist(),
            prot_sub.index.tolist(),
            metab_sub.index.tolist()
        ]
    )

    # 모델 옵션
    ent.set_model_options(
        factors           = 15,      # 초기 넉넉하게 → auto-pruning
        spikeslab_weights = True,    # sparse weights
        ard_factors       = True,
        ard_weights       = True
    )

    # 훈련 옵션
    ent.set_train_options(
        iter              = 1000,
        convergence_mode  = "medium",
        startELBO         = 1,
        freqELBO          = 5,
        seed              = 42,
        verbose           = False
    )

    # 훈련 실행
    ent.build()
    ent.run()

    # 모델 저장
    mofa_outfile = str(OUTDIR / "mofa_model.hdf5")
    ent.save(mofa_outfile)
    print(f"  모델 저장: {mofa_outfile}")

    # MOFA2 R 패키지로 분석하는 경우를 위해 hdf5 저장
    # R에서: library(MOFA2); model <- load_model("mofa_model.hdf5")

    # Factor 정보 추출
    factors    = ent.model.nodes["Z"].getExpectation()   # (n_factors, n_samples)
    r2_df      = pd.DataFrame(ent.model.calculate_variance_explained())

    print(f"  남은 Factor 수: {factors.shape[0]}")
    factors_df = pd.DataFrame(
        factors.T,
        index=common_samples,
        columns=[f"Factor{i+1}" for i in range(factors.shape[0])]
    )
    factors_df.to_csv(OUTDIR / "mofa_factor_scores.csv")

except ImportError:
    print("  mofapy2 미설치 → pip install mofapy2")
    print("  R 대안: library(MOFA2) 사용 권장")

    # R 코드 출력
    r_code = """
# R에서 MOFA2 실행:
library(MOFA2)

# 데이터 준비 (각 오믹스: feature × sample)
data_list <- list(
  mRNA       = as.matrix(read.csv("mrna_vst_matrix.csv", row.names=1)),
  Protein    = as.matrix(read.csv("proteomics_log2_normalized.csv", row.names=1)),
  Metabolite = as.matrix(read.csv("metabolomics_log2.csv", row.names=1))
)

# MOFA 객체 생성
mofa <- create_mofa(data_list)

# 옵션 설정
data_opts  <- get_default_data_options(mofa)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors    <- 15
model_opts$spikeslab_weights <- TRUE

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "medium"
train_opts$seed             <- 42

mofa <- prepare_mofa(mofa,
  data_options=data_opts,
  model_options=model_opts,
  training_options=train_opts
)

# 훈련
mofa <- run_mofa(mofa, outfile="mofa_model.hdf5")

# Factor-표현형 상관
factors <- get_factors(mofa)$Periodontitis_Study
group_label <- ifelse(meta$group=="Periodontitis", 1, 0)
cor_results <- cor(factors, group_label, use="pairwise.complete.obs")

# 분산 설명 시각화
plot_variance_explained(mofa, max_r2=15)
plot_factor(mofa, factor=1:5, color_by="group")
"""
    print(r_code)

    # 대체 분석: PCA 기반 간이 통합
    print("\n  대체: sklearn PCA 기반 간이 통합...")
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    # 스케일 통일 후 concat
    def scale_df(df):
        sc = StandardScaler()
        return pd.DataFrame(sc.fit_transform(df.T).T,
                             index=df.index, columns=df.columns)

    combined = pd.concat([
        scale_df(mrna_sub).add_prefix("mRNA_"),
        scale_df(prot_sub).add_prefix("prot_"),
        scale_df(metab_sub).add_prefix("metab_")
    ])
    combined = combined[common_samples]

    pca_joint = PCA(n_components=10)
    factors_arr = pca_joint.fit_transform(combined.T)
    factors_df  = pd.DataFrame(
        factors_arr,
        index=common_samples,
        columns=[f"Factor{i+1}" for i in range(10)]
    )
    factors_df.to_csv(OUTDIR / "mofa_factor_scores_pca_approx.csv")
    print(f"  PCA 근사 완료 (Factor 10개)")


# ══════════════════════════════════════════════════════════════
# STEP 3: Factor-표현형 상관 분석
# ══════════════════════════════════════════════════════════════
print("\n[3/4] Factor-표현형 상관 분석...")

factor_files = list(OUTDIR.glob("mofa_factor_scores*.csv"))
if factor_files:
    factors_df = pd.read_csv(factor_files[0], index_col=0)

    corr_results = []
    for factor in factors_df.columns:
        r, p = stats.spearmanr(factors_df[factor], y)
        corr_results.append({
            'factor': factor,
            'spearman_r': r,
            'p_value':    p
        })

    corr_df = pd.DataFrame(corr_results)
    corr_df['padj'] = pd.Series(
        __import__('statsmodels').stats.multitest.multipletests(
            corr_df['p_value'], method='fdr_bh')[1]
    )
    corr_df = corr_df.sort_values('p_value')

    # 유의 Factor
    sig_factors = corr_df[corr_df['padj'] < 0.05]['factor'].tolist()
    print(f"  치주염 연관 유의 Factor: {sig_factors}")
    corr_df.to_csv(OUTDIR / "mofa_factor_phenotype_corr.csv", index=False)

    # Factor score plot (Factor1 vs Factor2)
    if factors_df.shape[1] >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))

        # Factor score scatter
        ax = axes[0]
        colors = ['#4CAF50' if g=='Control' else '#F44336'
                  for g in meta.loc[factors_df.index, 'group']]
        ax.scatter(factors_df.iloc[:,0], factors_df.iloc[:,1],
                   c=colors, s=80, alpha=0.8, edgecolors='white')
        ax.set_xlabel(f"Factor 1 (r={corr_df.iloc[0]['spearman_r']:.2f})", fontsize=11)
        ax.set_ylabel("Factor 2", fontsize=11)
        ax.set_title("MOFA+ Factor Scores", fontsize=12, fontweight='bold')
        from matplotlib.patches import Patch
        legend_el = [Patch(facecolor='#4CAF50', label='Control'),
                     Patch(facecolor='#F44336', label='Periodontitis')]
        ax.legend(handles=legend_el)
        ax.grid(True, alpha=0.3)

        # Factor-Trait 상관 막대
        ax = axes[1]
        colors_bar = ['#F44336' if p < 0.05 else '#9E9E9E'
                      for p in corr_df['padj']]
        bars = ax.bar(corr_df['factor'], corr_df['spearman_r'].abs(),
                      color=colors_bar, alpha=0.8, edgecolor='white')
        ax.axhline(y=0.3, color='red', linestyle='--', alpha=0.7, label='|r|=0.3')
        ax.set_xlabel("Factor"); ax.set_ylabel("|Spearman r|")
        ax.set_title("Factor–Periodontitis Correlation", fontsize=12, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(); ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(FIGDIR / "Fig_MOFA_factors.pdf", dpi=300, bbox_inches='tight')
        plt.close()


# ══════════════════════════════════════════════════════════════
# STEP 4: Factor Loading — Top 피처 추출
# ══════════════════════════════════════════════════════════════
print("\n[4/4] Factor Loading 피처 추출...")
# (MOFA2 hdf5에서 loading 추출은 R MOFA2 패키지로 수행 권장)
# R 코드:
print("""
# R에서 Loading 추출:
# loadings <- get_weights(mofa, views="all", factors=1:5)
# top_mrna_f1 <- head(sort(abs(loadings$mRNA[,1]), decreasing=T), 30)
# top_prot_f1 <- head(sort(abs(loadings$Protein[,1]), decreasing=T), 20)
# top_metab_f1 <- head(sort(abs(loadings$Metabolite[,1]), decreasing=T), 20)
""")

print("\n✅ MOFA+ Analysis 완료!")
