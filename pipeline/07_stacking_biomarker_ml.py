"""
Multi-omics Biomarker Discovery: Stacking Ensemble + SHAP
치주염 정상(30) vs 환자(30) — mRNA + Proteomics + Metabolomics

파이프라인:
  1. LASSO (feature pre-selection, 1-SE rule)
  2. Random Forest (SHAP 기반)
  3. XGBoost (SHAP 기반)
  4. Stacking Ensemble (LR meta-learner)
  5. LOOCV 검증
  6. Bootstrap AUC CI
  7. 단독 vs 통합 AUC 비교
  8. 최소 임상 패널 선정
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ── sklearn / ML
from sklearn.linear_model import LogisticRegression, LassoCV
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
from sklearn.utils import resample
from statsmodels.stats.multitest import multipletests

from xgboost import XGBClassifier
import shap

# ─────────────────────────────────────────────
# 0. 설정
# ─────────────────────────────────────────────
SEED = 42
np.random.seed(SEED)

RESULT_DIR = Path("../results/ml")
FIG_DIR    = Path("../figures")
RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────
# 1. 데이터 로드 (단일 오믹스 분석 결과 입력)
# ─────────────────────────────────────────────
def load_integrated_features(
    mrna_deg_path: str,       # DESeq2 결과 (padj<0.05, |log2FC|>1)
    prot_dep_path: str,       # limma 결과 (adj.p<0.05, |log2FC|>0.58)
    metab_path: str,          # OPLS-DA VIP>1 + t-test p<0.05 결과
    mrna_expr_path: str,      # 전체 mRNA 발현 행렬 (gene × sample)
    prot_expr_path: str,      # 전체 단백질 발현 행렬 (protein × sample)
    metab_expr_path: str,     # 전체 대사체 발현 행렬 (metabolite × sample)
    sample_meta_path: str     # 샘플 메타데이터 (sample, group: Control/Periodontitis)
) -> tuple:
    """
    Returns:
        X_dict: {'mRNA': df, 'Protein': df, 'Metabolite': df}
        X_integrated: pd.DataFrame (모든 오믹스 통합)
        y: pd.Series (0=Control, 1=Periodontitis)
    """
    # 샘플 메타데이터
    meta = pd.read_csv(sample_meta_path, index_col=0)
    y = (meta['group'] == 'Periodontitis').astype(int)
    samples = meta.index.tolist()

    # 유의 피처 목록 로드
    deg  = pd.read_csv(mrna_deg_path, index_col=0)
    dep  = pd.read_csv(prot_dep_path, index_col=0)
    dmet = pd.read_csv(metab_path, index_col=0)

    # 유의 피처 필터링
    sig_genes   = deg[deg['padj'] < 0.05].index.tolist()
    sig_prots   = dep[dep['adj.P.Val'] < 0.05].index.tolist()
    sig_metabs  = dmet[dmet['pvalue'] < 0.05].index.tolist()

    # 발현 행렬 로드 및 유의 피처만 선택
    mrna_mat  = pd.read_csv(mrna_expr_path, index_col=0).loc[sig_genes, samples].T
    prot_mat  = pd.read_csv(prot_expr_path, index_col=0).loc[sig_prots, samples].T
    metab_mat = pd.read_csv(metab_expr_path, index_col=0).loc[sig_metabs, samples].T

    # 오믹스별 prefix 추가
    mrna_mat.columns  = ['mRNA_' + c for c in mrna_mat.columns]
    prot_mat.columns  = ['prot_' + c for c in prot_mat.columns]
    metab_mat.columns = ['metab_' + c for c in metab_mat.columns]

    X_dict = {
        'mRNA':       mrna_mat,
        'Protein':    prot_mat,
        'Metabolite': metab_mat
    }
    X_integrated = pd.concat([mrna_mat, prot_mat, metab_mat], axis=1)

    print(f"피처 수: mRNA={len(sig_genes)}, Protein={len(sig_prots)}, Metabolite={len(sig_metabs)}")
    print(f"통합 피처: {X_integrated.shape[1]}, 샘플: {X_integrated.shape[0]}")

    return X_dict, X_integrated, y


# ─────────────────────────────────────────────
# 2. STEP 1: LASSO Feature Pre-selection (1-SE rule)
# ─────────────────────────────────────────────
def lasso_feature_selection(X: pd.DataFrame, y: pd.Series, cv: int = 10) -> list:
    """
    LASSO with 1-SE rule (보수적 feature 선택 → 작은 임상 패널)
    
    1-SE rule: lambda_min 대신 lambda_min + 1 × SE 사용
    → 더 희소한 모델 → 임상 적용 용이
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    lasso_cv = LassoCV(
        cv=cv,
        random_state=SEED,
        max_iter=10000,
        n_alphas=100
    ).fit(X_scaled, y)

    # 1-SE rule 적용
    # LassoCV는 mean_mse_path_ 제공
    mean_mse  = lasso_cv.mse_path_.mean(axis=-1)
    std_mse   = lasso_cv.mse_path_.std(axis=-1)
    min_idx   = np.argmin(mean_mse)
    threshold = mean_mse[min_idx] + std_mse[min_idx]   # 1-SE 기준
    alphas    = lasso_cv.alphas_

    # 1-SE rule에 해당하는 가장 큰 alpha (희소 모델)
    alpha_1se_idx = np.where(mean_mse <= threshold)[0][-1]   # 가장 희소
    alpha_1se = alphas[alpha_1se_idx]

    from sklearn.linear_model import Lasso
    lasso_1se = Lasso(alpha=alpha_1se, max_iter=10000, random_state=SEED)
    lasso_1se.fit(X_scaled, y)

    selected = X.columns[lasso_1se.coef_ != 0].tolist()
    print(f"\nLASSO 1-SE rule: {len(selected)}개 피처 선택 (전체 {X.shape[1]}개 중)")
    print(f"  alpha_min={lasso_cv.alpha_:.4f}, alpha_1SE={alpha_1se:.4f}")

    return selected, lasso_1se.coef_[lasso_1se.coef_ != 0], X.columns[lasso_1se.coef_ != 0]


# ─────────────────────────────────────────────
# 3. STEP 2: LOOCV Stacking Ensemble
# ─────────────────────────────────────────────
def build_stacking_ensemble():
    """
    Base learners: LASSO(LR) + RF + XGBoost
    Meta-learner: Logistic Regression (최종 통합)
    """
    base_learners = [
        ('lasso_lr', Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(
                penalty='l1', solver='liblinear',
                C=0.1, max_iter=1000, random_state=SEED
            ))
        ])),
        ('rf', RandomForestClassifier(
            n_estimators=500,
            max_depth=3,          # 소규모 데이터 → 얕은 트리
            min_samples_leaf=3,
            random_state=SEED,
            n_jobs=-1
        )),
        ('xgb', XGBClassifier(
            n_estimators=200,
            max_depth=3,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            reg_alpha=0.1,        # L1 regularization
            reg_lambda=1.0,       # L2 regularization
            use_label_encoder=False,
            eval_metric='logloss',
            random_state=SEED,
            verbosity=0
        ))
    ]

    meta_learner = LogisticRegression(
        C=1.0, max_iter=1000, random_state=SEED
    )

    stack = StackingClassifier(
        estimators=base_learners,
        final_estimator=meta_learner,
        cv=5,                     # ⚠️ n=60에서 stacking 내부 CV
        stack_method='predict_proba',
        passthrough=False,        # 원본 피처는 meta에 전달 안 함
        n_jobs=-1
    )
    return stack


def loocv_evaluation(X: pd.DataFrame, y: pd.Series, model) -> dict:
    """
    LOOCV로 각 모델 평가
    n=60에서 가장 적합한 검증 방법
    """
    loo = LeaveOneOut()
    oof_probs  = np.zeros(len(y))
    oof_labels = np.zeros(len(y))

    for i, (train_idx, test_idx) in enumerate(loo.split(X)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        model_clone = model  # sklearn clone 필요시 추가
        try:
            import copy
            m = copy.deepcopy(model)
            m.fit(X_train, y_train)
            oof_probs[test_idx]  = m.predict_proba(X_test)[:, 1]
            oof_labels[test_idx] = y_test.values
        except Exception as e:
            print(f"  Fold {i} error: {e}")
            oof_probs[test_idx] = 0.5

    auc = roc_auc_score(oof_labels, oof_probs)
    fpr, tpr, thresholds = roc_curve(oof_labels, oof_probs)

    # 최적 임계값 (Youden's J)
    j_scores = tpr - fpr
    opt_thresh = thresholds[np.argmax(j_scores)]
    y_pred = (oof_probs >= opt_thresh).astype(int)

    return {
        'auc': auc,
        'fpr': fpr,
        'tpr': tpr,
        'probs': oof_probs,
        'labels': oof_labels,
        'opt_thresh': opt_thresh,
        'y_pred': y_pred
    }


# ─────────────────────────────────────────────
# 4. Bootstrap AUC 신뢰구간
# ─────────────────────────────────────────────
def bootstrap_auc_ci(y_true, y_prob, n_boot=1000, ci=0.95) -> dict:
    """
    Bootstrap으로 AUC 95% CI 계산
    n=60 소규모 데이터에서 필수
    """
    auc_scores = []
    for _ in range(n_boot):
        idx = resample(np.arange(len(y_true)), random_state=None)
        if len(np.unique(y_true[idx])) < 2:
            continue
        auc_scores.append(roc_auc_score(y_true[idx], y_prob[idx]))

    lower = np.percentile(auc_scores, (1 - ci) / 2 * 100)
    upper = np.percentile(auc_scores, (1 + ci) / 2 * 100)
    return {'mean': np.mean(auc_scores), 'lower': lower, 'upper': upper}


# ─────────────────────────────────────────────
# 5. 단독 vs 통합 AUC 비교 (핵심 Figure)
# ─────────────────────────────────────────────
def compare_omics_auc(X_dict: dict, X_integrated: pd.DataFrame,
                       y: pd.Series, selected_features: list) -> pd.DataFrame:
    """
    각 오믹스 단독 vs 통합 패널 AUC 비교
    → 멀티오믹스 통합의 당위성 수치 증명
    """
    results = []

    # 각 오믹스 단독
    for omics_name, X_omics in X_dict.items():
        # 해당 오믹스의 선택된 피처만
        prefix = {'mRNA': 'mRNA_', 'Protein': 'prot_', 'Metabolite': 'metab_'}[omics_name]
        sel = [f for f in selected_features if f.startswith(prefix)]
        if len(sel) < 2:
            continue

        model = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(C=0.1, random_state=SEED, max_iter=1000))
        ])
        res = loocv_evaluation(X_integrated[sel], y, model)
        ci  = bootstrap_auc_ci(res['labels'], res['probs'])

        results.append({
            'Model': f'{omics_name} only',
            'AUC': res['auc'],
            'CI_lower': ci['lower'],
            'CI_upper': ci['upper'],
            'n_features': len(sel)
        })
        print(f"  {omics_name} only: AUC={res['auc']:.3f} [{ci['lower']:.3f}-{ci['upper']:.3f}]")

    # 통합 패널
    stack = build_stacking_ensemble()
    res_int = loocv_evaluation(X_integrated[selected_features], y, stack)
    ci_int  = bootstrap_auc_ci(res_int['labels'], res_int['probs'])

    results.append({
        'Model': 'Integrated Panel',
        'AUC': res_int['auc'],
        'CI_lower': ci_int['lower'],
        'CI_upper': ci_int['upper'],
        'n_features': len(selected_features)
    })
    print(f"  Integrated: AUC={res_int['auc']:.3f} [{ci_int['lower']:.3f}-{ci_int['upper']:.3f}]")

    return pd.DataFrame(results), res_int


# ─────────────────────────────────────────────
# 6. SHAP 분석
# ─────────────────────────────────────────────
def shap_analysis(X: pd.DataFrame, y: pd.Series,
                  selected_features: list) -> dict:
    """
    XGBoost + SHAP beeswarm + 오믹스별 기여도
    """
    X_sel = X[selected_features].copy()

    # 최종 XGBoost 학습 (전체 데이터)
    xgb_final = XGBClassifier(
        n_estimators=200, max_depth=3, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        use_label_encoder=False, eval_metric='logloss',
        random_state=SEED, verbosity=0
    )
    xgb_final.fit(X_sel, y)

    # SHAP 계산
    explainer   = shap.TreeExplainer(xgb_final)
    shap_values = explainer.shap_values(X_sel)  # (n_samples, n_features)

    # 오믹스별 SHAP 기여도
    omics_contrib = {}
    for prefix, name in [('mRNA_', 'mRNA'), ('prot_', 'Protein'), ('metab_', 'Metabolite')]:
        idx = [i for i, f in enumerate(selected_features) if f.startswith(prefix)]
        if idx:
            omics_contrib[name] = np.abs(shap_values[:, idx]).mean()

    print("\nSHAP 오믹스별 기여도:")
    total = sum(omics_contrib.values())
    for k, v in omics_contrib.items():
        print(f"  {k}: {v:.4f} ({v/total*100:.1f}%)")

    return {
        'shap_values': shap_values,
        'explainer': explainer,
        'model': xgb_final,
        'X': X_sel,
        'omics_contrib': omics_contrib,
        'features': selected_features
    }


# ─────────────────────────────────────────────
# 7. Figure 생성
# ─────────────────────────────────────────────
def plot_roc_comparison(auc_df: pd.DataFrame, res_integrated: dict,
                         save_path: str = None):
    """Figure: 단독 vs 통합 ROC 곡선 + AUC 막대 (2-panel)"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: ROC curves
    ax = axes[0]
    colors = {'mRNA only': '#2196F3', 'Protein only': '#FF5722',
              'Metabolite only': '#4CAF50', 'Integrated Panel': '#9C27B0'}

    # 통합 ROC
    fpr, tpr = res_integrated['fpr'], res_integrated['tpr']
    auc_val = res_integrated['auc']
    ax.plot(fpr, tpr, color='#9C27B0', lw=2.5,
            label=f'Integrated ({auc_val:.3f})', zorder=5)
    ax.fill_between(fpr, tpr, alpha=0.1, color='#9C27B0')

    ax.plot([0,1],[0,1], 'k--', lw=1, alpha=0.5)
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('A. ROC Curves (LOOCV)', fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.05])
    ax.grid(True, alpha=0.3)

    # Panel B: AUC bar + CI
    ax = axes[1]
    models = auc_df['Model'].tolist()
    aucs   = auc_df['AUC'].tolist()
    lowers = auc_df['CI_lower'].tolist()
    uppers = auc_df['CI_upper'].tolist()
    bar_colors = [colors.get(m, '#607D8B') for m in models]

    bars = ax.bar(range(len(models)), aucs, color=bar_colors,
                  alpha=0.85, edgecolor='white', linewidth=1.5)
    # CI error bars
    yerr_low  = [a - l for a, l in zip(aucs, lowers)]
    yerr_high = [u - a for a, u in zip(aucs, uppers)]
    ax.errorbar(range(len(models)), aucs,
                yerr=[yerr_low, yerr_high],
                fmt='none', color='black', capsize=5, linewidth=1.5)

    # AUC 수치 표시
    for i, (bar, auc) in enumerate(zip(bars, aucs)):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{auc:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=20, ha='right', fontsize=10)
    ax.set_ylabel('AUC (LOOCV)', fontsize=12)
    ax.set_title('B. Single vs Integrated Omics AUC', fontsize=13, fontweight='bold')
    ax.set_ylim([0.4, 1.05])
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random')
    ax.axhline(y=0.8, color='green', linestyle=':', alpha=0.5, label='AUC=0.8')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


def plot_shap_summary(shap_data: dict, save_path: str = None):
    """Figure: SHAP beeswarm (top20) + 오믹스 기여도 파이차트 (2-panel)"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: SHAP beeswarm
    ax = axes[0]
    plt.sca(ax)
    shap.summary_plot(
        shap_data['shap_values'],
        shap_data['X'],
        max_display=20,
        show=False,
        plot_size=None
    )
    ax.set_title('A. SHAP Feature Importance (Top 20)', fontsize=13, fontweight='bold')

    # Panel B: 오믹스별 기여도 파이차트
    ax = axes[1]
    contrib = shap_data['omics_contrib']
    colors_pie = ['#2196F3', '#FF5722', '#4CAF50']
    wedges, texts, autotexts = ax.pie(
        list(contrib.values()),
        labels=list(contrib.keys()),
        colors=colors_pie[:len(contrib)],
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops={'edgecolor': 'white', 'linewidth': 2}
    )
    for text in autotexts:
        text.set_fontsize(12)
        text.set_fontweight('bold')
    ax.set_title('B. Omics Layer SHAP Contribution', fontsize=13, fontweight='bold')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


def plot_minimal_panel_roc(X: pd.DataFrame, y: pd.Series,
                            final_panel: list, save_path: str = None):
    """Figure: 최소 임상 패널 ROC + 상세 성능 지표"""
    from sklearn.metrics import confusion_matrix, roc_auc_score

    stack = build_stacking_ensemble()
    res   = loocv_evaluation(X[final_panel], y, stack)
    ci    = bootstrap_auc_ci(res['labels'], res['probs'])

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A: ROC
    ax = axes[0]
    ax.plot(res['fpr'], res['tpr'], color='#9C27B0', lw=2.5,
            label=f'AUC={res["auc"]:.3f} [{ci["lower"]:.3f}-{ci["upper"]:.3f}]')
    ax.fill_between(res['fpr'], res['tpr'], alpha=0.15, color='#9C27B0')
    ax.plot([0,1],[0,1],'k--',lw=1,alpha=0.5)
    ax.scatter([res['opt_thresh']], [res['tpr'][np.argmin(np.abs(res['fpr']-res['opt_thresh']))]],
               color='red', s=100, zorder=5, label='Optimal threshold')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'A. Minimal Biomarker Panel ROC\n(n={len(final_panel)} features)',
                 fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel B: Confusion matrix
    ax = axes[1]
    cm = confusion_matrix(res['labels'], res['y_pred'])
    sns.heatmap(cm, annot=True, fmt='d', cmap='Purples',
                xticklabels=['Control', 'Periodontitis'],
                yticklabels=['Control', 'Periodontitis'],
                ax=ax, cbar=False, annot_kws={'size': 14})
    tn, fp, fn, tp = cm.ravel()
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    ppv  = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv  = tn / (tn + fn) if (tn + fn) > 0 else 0
    ax.set_title(
        f'B. Confusion Matrix\n'
        f'Sensitivity={sens:.3f}, Specificity={spec:.3f}\n'
        f'PPV={ppv:.3f}, NPV={npv:.3f}',
        fontsize=11, fontweight='bold'
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


# ─────────────────────────────────────────────
# 8. 메인 실행
# ─────────────────────────────────────────────
def run_ml_pipeline(
    mrna_deg_path, prot_dep_path, metab_path,
    mrna_expr_path, prot_expr_path, metab_expr_path,
    sample_meta_path
):
    print("=" * 60)
    print("Multi-omics Biomarker ML Pipeline — Periodontitis")
    print("=" * 60)

    # 1. 데이터 로드
    print("\n[1/6] 데이터 로드 중...")
    X_dict, X_integrated, y = load_integrated_features(
        mrna_deg_path, prot_dep_path, metab_path,
        mrna_expr_path, prot_expr_path, metab_expr_path,
        sample_meta_path
    )

    # 2. LASSO feature selection
    print("\n[2/6] LASSO 1-SE rule feature selection...")
    selected, coefs, feat_names = lasso_feature_selection(X_integrated, y, cv=10)

    # LASSO 선택 피처 저장
    lasso_df = pd.DataFrame({
        'feature': feat_names,
        'lasso_coef': coefs,
        'omics': ['mRNA' if f.startswith('mRNA_') else
                  'Protein' if f.startswith('prot_') else
                  'Metabolite' for f in feat_names]
    })
    lasso_df.to_csv(RESULT_DIR / 'lasso_selected_features.csv')
    print(f"  omic breakdown: {lasso_df.groupby('omics').size().to_dict()}")

    # 3. 오믹스별 단독 vs 통합 AUC 비교
    print("\n[3/6] 오믹스별 단독 vs 통합 AUC 비교...")
    auc_df, res_integrated = compare_omics_auc(X_dict, X_integrated, y, selected)
    auc_df.to_csv(RESULT_DIR / 'omics_auc_comparison.csv', index=False)

    # Figure 9: ROC + AUC 비교
    plot_roc_comparison(auc_df, res_integrated,
                         save_path=str(FIG_DIR / 'Fig9_ROC_comparison.png'))

    # 4. SHAP 분석
    print("\n[4/6] SHAP 분석...")
    shap_data = shap_analysis(X_integrated, y, selected)

    # Figure 10: SHAP beeswarm + 오믹스 기여도
    plot_shap_summary(shap_data,
                       save_path=str(FIG_DIR / 'Fig10_SHAP_summary.png'))

    # 5. 최소 임상 패널 선정 (SHAP top 10 기준)
    print("\n[5/6] 최소 임상 패널 선정...")
    mean_shap = np.abs(shap_data['shap_values']).mean(axis=0)
    top_idx   = np.argsort(mean_shap)[::-1][:10]
    final_panel = [selected[i] for i in top_idx]

    print(f"최소 패널 ({len(final_panel)}개):")
    for i, f in enumerate(final_panel):
        omics = 'mRNA' if f.startswith('mRNA_') else 'Protein' if f.startswith('prot_') else 'Metabolite'
        print(f"  {i+1}. [{omics}] {f.split('_',1)[1]} (SHAP={mean_shap[top_idx[i]]:.4f})")

    # Figure 11: 최소 패널 ROC
    plot_minimal_panel_roc(X_integrated, y, final_panel,
                            save_path=str(FIG_DIR / 'Fig11_minimal_panel_ROC.png'))

    # 6. 결과 저장
    print("\n[6/6] 결과 저장...")
    auc_df.to_csv(RESULT_DIR / 'final_auc_summary.csv', index=False)

    final_panel_df = pd.DataFrame({
        'rank': range(1, len(final_panel)+1),
        'feature': final_panel,
        'gene_name': [f.split('_',1)[1] for f in final_panel],
        'omics': ['mRNA' if f.startswith('mRNA_') else
                  'Protein' if f.startswith('prot_') else
                  'Metabolite' for f in final_panel],
        'mean_shap': [mean_shap[top_idx[i]] for i in range(len(final_panel))]
    })
    final_panel_df.to_csv(RESULT_DIR / 'final_biomarker_panel.csv', index=False)

    print("\n" + "=" * 60)
    print("✅ ML Pipeline 완료!")
    print(f"  최종 통합 AUC: {res_integrated['auc']:.3f}")
    print(f"  최소 패널: {len(final_panel)}개 바이오마커")
    print("=" * 60)

    return {
        'selected_features': selected,
        'final_panel': final_panel,
        'auc_summary': auc_df,
        'shap_data': shap_data,
        'res_integrated': res_integrated
    }


# ─────────────────────────────────────────────
# 실행 예시 (데이터 준비 후 경로 수정)
# ─────────────────────────────────────────────
if __name__ == "__main__":
    results = run_ml_pipeline(
        mrna_deg_path    = "../results/mrna/DESeq2_results.csv",
        prot_dep_path    = "../results/proteomics/limma_results.csv",
        metab_path       = "../results/metabolomics/oplsda_results.csv",
        mrna_expr_path   = "../data/processed/mrna_vst.csv",
        prot_expr_path   = "../data/processed/proteomics_log2.csv",
        metab_expr_path  = "../data/processed/metabolomics_log2.csv",
        sample_meta_path = "../data/sample_metadata.csv"
    )
