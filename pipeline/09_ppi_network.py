"""
PPI Network Analysis — STRING API + Hub Gene 발굴
DEG + DEP 통합 → PPI → Hub 유전자 → 바이오마커 후보 연계

분석:
  1. STRING API로 PPI 네트워크 구축
  2. NetworkX로 네트워크 분석
  3. Hub gene: degree + betweenness + closeness 중심성
  4. 공통 Hub: DEG ∩ DEP ∩ PPI_hub → 최강 바이오마커
  5. 시각화
"""

import requests
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

RESULT_DIR = Path("../results/ppi")
FIG_DIR    = Path("../figures")
RESULT_DIR.mkdir(parents=True, exist_ok=True)

STRING_API = "https://string-db.org/api"
SPECIES    = 9606   # Human


# ─────────────────────────────────────────────
# 1. STRING API — PPI 네트워크 구축
# ─────────────────────────────────────────────
def get_string_network(
    gene_list: list,
    score_threshold: int = 400,    # 0.4 = medium confidence (표준)
    add_nodes: int = 0             # 외부 노드 확장 (0=없음)
) -> pd.DataFrame:
    """
    STRING API로 PPI 상호작용 데이터 가져오기
    score_threshold:
      150 = low confidence
      400 = medium confidence ← 권장
      700 = high confidence
      900 = highest confidence
    """
    # 유전자 목록을 STRING ID로 변환
    params_map = {
        "identifiers": "\r".join(gene_list),
        "species": SPECIES,
        "format": "json"
    }
    response = requests.post(
        f"{STRING_API}/json/get_string_ids",
        data=params_map, timeout=60
    )
    string_ids = response.json()
    id_map = {item['queryItem']: item['stringId']
              for item in string_ids if 'stringId' in item}

    print(f"STRING ID 매핑: {len(id_map)}/{len(gene_list)}개 성공")

    # 상호작용 네트워크 요청
    mapped_ids = list(id_map.values())
    params_net = {
        "identifiers": "\r".join(mapped_ids),
        "species": SPECIES,
        "required_score": score_threshold,
        "add_nodes": add_nodes,
        "caller_identity": "multiomics_periodontitis"
    }
    response = requests.post(
        f"{STRING_API}/tsv/network",
        data=params_net, timeout=120
    )

    lines = response.text.strip().split('\n')
    if len(lines) < 2:
        print("⚠️ 상호작용 데이터 없음. threshold를 낮춰보세요.")
        return pd.DataFrame()

    # TSV 파싱
    headers = lines[0].split('\t')
    rows = [line.split('\t') for line in lines[1:]]
    df = pd.DataFrame(rows, columns=headers)

    # STRING ID → gene symbol 역매핑
    rev_map = {v: k for k, v in id_map.items()}
    if 'stringId_A' in df.columns:
        df['gene_A'] = df['stringId_A'].map(rev_map).fillna(df['stringId_A'])
        df['gene_B'] = df['stringId_B'].map(rev_map).fillna(df['stringId_B'])
        df['score']  = df['score'].astype(float)

    print(f"PPI 상호작용: {len(df)}개 (threshold={score_threshold})")
    return df


# ─────────────────────────────────────────────
# 2. NetworkX 네트워크 구축 + Hub Gene 분석
# ─────────────────────────────────────────────
def build_network_and_find_hubs(
    ppi_df: pd.DataFrame,
    top_n_hub: int = 20
) -> tuple:
    """
    NetworkX 그래프 구축 + 중심성 지표로 Hub gene 선별

    Hub gene 기준:
      - Degree centrality (연결 수) ← 가장 직관적
      - Betweenness centrality (정보 흐름의 병목)
      - Closeness centrality (네트워크 중심 근접성)
    """
    G = nx.Graph()

    for _, row in ppi_df.iterrows():
        G.add_edge(
            row['gene_A'], row['gene_B'],
            weight=float(row.get('score', 1.0))
        )

    print(f"\n네트워크 통계:")
    print(f"  노드 수: {G.number_of_nodes()}")
    print(f"  엣지 수: {G.number_of_edges()}")
    print(f"  평균 degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")

    # 중심성 계산
    degree_cent      = nx.degree_centrality(G)
    betweenness_cent = nx.betweenness_centrality(G, normalized=True)
    closeness_cent   = nx.closeness_centrality(G)

    # 통합 점수 (가중 합산)
    nodes = list(G.nodes())
    cent_df = pd.DataFrame({
        'gene':        nodes,
        'degree':      [G.degree(n) for n in nodes],
        'degree_cent': [degree_cent[n] for n in nodes],
        'betweenness': [betweenness_cent[n] for n in nodes],
        'closeness':   [closeness_cent[n] for n in nodes]
    })

    # 각 중심성 정규화 후 합산
    for col in ['degree_cent', 'betweenness', 'closeness']:
        min_v, max_v = cent_df[col].min(), cent_df[col].max()
        if max_v > min_v:
            cent_df[f'{col}_norm'] = (cent_df[col] - min_v) / (max_v - min_v)
        else:
            cent_df[f'{col}_norm'] = 0.0

    cent_df['hub_score'] = (
        cent_df['degree_cent_norm'] * 0.4 +
        cent_df['betweenness_norm'] * 0.4 +
        cent_df['closeness_norm']   * 0.2
    )

    cent_df = cent_df.sort_values('hub_score', ascending=False)
    hub_genes = cent_df.head(top_n_hub)['gene'].tolist()

    print(f"\nTop {top_n_hub} Hub Genes:")
    for i, row in cent_df.head(top_n_hub).iterrows():
        print(f"  {row['gene']}: score={row['hub_score']:.3f}, "
              f"degree={row['degree']}, "
              f"betweenness={row['betweenness']:.3f}")

    return G, cent_df, hub_genes


# ─────────────────────────────────────────────
# 3. 공통 Hub 유전자 발굴 (DEG ∩ DEP ∩ PPI_hub)
# ─────────────────────────────────────────────
def find_cross_omics_hubs(
    hub_genes_ppi: list,
    deg_genes: list,
    dep_proteins: list,
    wgcna_hubs_mrna: list = None,
    wgcna_hubs_prot: list = None
) -> dict:
    """
    다중 증거 교집합으로 최강 바이오마커 후보 선별
    증거 수준:
      ⭐⭐⭐⭐⭐ DEG ∩ DEP ∩ PPI_hub (+ WGCNA hub)
      ⭐⭐⭐⭐   DEG ∩ DEP ∩ PPI_hub
      ⭐⭐⭐    DEG ∩ PPI_hub 또는 DEP ∩ PPI_hub
    """
    hub_set = set(hub_genes_ppi)
    deg_set = set(deg_genes)
    dep_set = set(dep_proteins)

    result = {
        'level5': set(),  # 모든 증거
        'level4': set(),  # DEG + DEP + PPI hub
        'level3': set(),  # DEG + PPI hub 또는 DEP + PPI hub
    }

    base = hub_set & deg_set & dep_set
    result['level4'] = base

    if wgcna_hubs_mrna and wgcna_hubs_prot:
        wgcna_set = set(wgcna_hubs_mrna) | set(wgcna_hubs_prot)
        result['level5'] = base & wgcna_set

    result['level3'] = (hub_set & deg_set) | (hub_set & dep_set)

    print("\n공통 Hub Gene 요약:")
    print(f"  ⭐⭐⭐⭐⭐ (DEG∩DEP∩PPI∩WGCNA): {len(result['level5'])}개 → {result['level5']}")
    print(f"  ⭐⭐⭐⭐  (DEG∩DEP∩PPI)       : {len(result['level4'])}개 → {result['level4']}")
    print(f"  ⭐⭐⭐   (DEG∩PPI or DEP∩PPI): {len(result['level3'])}개")

    return result


# ─────────────────────────────────────────────
# 4. PPI 네트워크 시각화
# ─────────────────────────────────────────────
def plot_ppi_network(
    G: nx.Graph,
    cent_df: pd.DataFrame,
    hub_genes: list,
    deg_set: set = None,
    dep_set: set = None,
    save_path: str = None
):
    """
    PPI 네트워크 시각화
    - 노드 크기 = degree 중심성
    - 노드 색 = 증거 유형 (DEG only / DEP only / Both / Hub)
    - 엣지 두께 = interaction score
    """
    fig, ax = plt.subplots(figsize=(14, 12))

    # 서브그래프 (Hub + 1-hop 이웃)
    hub_nodes = set(hub_genes)
    neighbors = set()
    for h in hub_genes:
        if h in G:
            neighbors.update(G.neighbors(h))
    subgraph_nodes = hub_nodes | (neighbors & hub_nodes)  # hub들만
    # 더 넓게: hub + 직접 연결된 유의 유전자
    all_sig = (deg_set or set()) | (dep_set or set())
    subgraph_nodes = hub_nodes | (neighbors & all_sig)
    H = G.subgraph(subgraph_nodes)

    # 레이아웃
    pos = nx.spring_layout(H, k=2.0, seed=42, weight='weight')

    # 노드 색 결정
    def get_node_color(n):
        in_deg = n in (deg_set or set())
        in_dep = n in (dep_set or set())
        in_hub = n in hub_genes
        if in_hub and in_deg and in_dep: return '#9C27B0'  # 보라: 모든 증거
        elif in_hub and in_deg:           return '#F44336'  # 빨강: DEG+hub
        elif in_hub and in_dep:           return '#FF9800'  # 주황: DEP+hub
        elif in_hub:                      return '#2196F3'  # 파랑: hub only
        elif in_deg and in_dep:           return '#E91E63'  # 분홍: DEG+DEP
        elif in_deg:                      return '#FF5722'  # 연빨: DEG only
        elif in_dep:                      return '#FFC107'  # 노랑: DEP only
        else:                             return '#9E9E9E'  # 회색: 연결만

    node_colors = [get_node_color(n) for n in H.nodes()]
    node_sizes  = [cent_df[cent_df['gene']==n]['degree'].values[0] * 80
                   if n in cent_df['gene'].values else 100
                   for n in H.nodes()]

    # 엣지 그리기
    edges = H.edges(data=True)
    edge_weights = [e[2].get('weight', 0.5) for e in edges]
    nx.draw_networkx_edges(H, pos, ax=ax, width=[w*2 for w in edge_weights],
                           alpha=0.3, edge_color='gray')

    # 노드 그리기
    nx.draw_networkx_nodes(H, pos, ax=ax,
                           node_color=node_colors,
                           node_size=node_sizes, alpha=0.9)

    # Hub gene 레이블만 표시
    hub_labels = {n: n for n in H.nodes() if n in hub_genes}
    nx.draw_networkx_labels(H, pos, hub_labels, ax=ax,
                            font_size=8, font_weight='bold')

    # 범례
    from matplotlib.patches import Patch
    legend = [
        Patch(facecolor='#9C27B0', label='Hub + DEG + DEP (all evidence)'),
        Patch(facecolor='#F44336', label='Hub + DEG'),
        Patch(facecolor='#FF9800', label='Hub + DEP'),
        Patch(facecolor='#2196F3', label='Hub only'),
        Patch(facecolor='#E91E63', label='DEG + DEP'),
        Patch(facecolor='#9E9E9E', label='Connected node'),
    ]
    ax.legend(handles=legend, loc='upper left', fontsize=9,
              bbox_to_anchor=(1, 1))

    ax.set_title('PPI Network — Hub Genes (Periodontitis)\n'
                 'Node size ∝ Degree, Edge width ∝ STRING score',
                 fontsize=13, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"저장: {save_path}")
    plt.show()


# ─────────────────────────────────────────────
# 5. 메인 실행
# ─────────────────────────────────────────────
def run_ppi_pipeline(
    deg_path: str,
    dep_path: str,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 1.0
):
    print("=" * 60)
    print("PPI Network Analysis — Periodontitis")
    print("=" * 60)

    # 유의 유전자/단백 로드
    deg_df = pd.read_csv(deg_path, index_col=0)
    dep_df = pd.read_csv(dep_path, index_col=0)

    sig_genes = deg_df[
        (deg_df['padj'] < padj_thresh) &
        (deg_df['log2FoldChange'].abs() > lfc_thresh)
    ].index.tolist()

    sig_prots = dep_df[
        (dep_df['adj.P.Val'] < padj_thresh) &
        (dep_df['logFC'].abs() > 0.58)
    ].index.tolist()

    print(f"\n입력: DEG={len(sig_genes)}, DEP={len(sig_prots)}")

    # 통합 유전자 목록 (mRNA + Protein 합집합)
    all_genes = list(set(sig_genes + sig_prots))

    # 1. STRING 네트워크 구축
    print("\n[1/4] STRING API 호출...")
    ppi_df = get_string_network(all_genes, score_threshold=400)

    if ppi_df.empty:
        print("PPI 데이터 없음. threshold를 낮춰 재시도합니다.")
        ppi_df = get_string_network(all_genes, score_threshold=200)

    ppi_df.to_csv(RESULT_DIR / 'ppi_interactions.csv', index=False)

    # 2. 네트워크 분석 + Hub gene
    print("\n[2/4] Hub gene 분석...")
    G, cent_df, hub_genes = build_network_and_find_hubs(ppi_df, top_n_hub=20)
    cent_df.to_csv(RESULT_DIR / 'hub_gene_centrality.csv', index=False)

    # 3. 공통 Hub 교집합
    print("\n[3/4] 공통 Hub 발굴...")
    cross_hubs = find_cross_omics_hubs(
        hub_genes, sig_genes, sig_prots
    )

    # 저장
    for level, genes in cross_hubs.items():
        if genes:
            pd.Series(list(genes)).to_csv(
                RESULT_DIR / f'cross_hub_{level}.csv', index=False, header=['gene']
            )

    # 4. 시각화
    print("\n[4/4] Figure 생성...")
    plot_ppi_network(
        G, cent_df, hub_genes,
        deg_set=set(sig_genes),
        dep_set=set(sig_prots),
        save_path=str(FIG_DIR / 'Fig_PPI_network.png')
    )

    print("\n✅ PPI Analysis 완료!")
    return G, cent_df, hub_genes, cross_hubs


if __name__ == "__main__":
    run_ppi_pipeline(
        deg_path="../results/mrna/DESeq2_results.csv",
        dep_path="../results/proteomics/limma_results.csv"
    )
