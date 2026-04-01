"""
11_llm_analysis_agent.py
LLM 기반 멀티오믹스 결과 해석 + 바이오마커 인사이트 발굴 에이전트

기능:
  1. 분석 결과 자동 수집 (DEG, DEP, 대사체, MOFA, SHAP, PPI)
  2. PubMed 자동 검색 (바이오마커 후보별 문헌 근거)
  3. Claude API 추론:
     - 생물학적 의미 해석
     - Cross-omics 패턴 → 새 가설
     - Discussion 초안 생성
     - 결론 도출
  4. Markdown 보고서 자동 생성

사용법:
  ANTHROPIC_API_KEY 환경변수 설정 후:
  python 11_llm_analysis_agent.py

또는 OpenAI API 사용 시:
  OPENAI_API_KEY 설정 후 --model gpt-4o 옵션
"""

import os
import json
import time
import requests
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Optional
import xml.etree.ElementTree as ET

# ─────────────────────────────────────────────
# 0. 설정
# ─────────────────────────────────────────────
RESULT_DIR  = Path("../results")
OUTDIR      = Path("../results/llm_insights")
OUTDIR.mkdir(parents=True, exist_ok=True)

TIMESTAMP   = datetime.now().strftime("%Y%m%d_%H%M%S")
REPORT_PATH = OUTDIR / f"biomarker_insights_{TIMESTAMP}.md"

# LLM 모델 선택 (Claude 권장)
DEFAULT_MODEL = "claude-opus-4-5"  # 또는 "claude-sonnet-4-5"


# ═══════════════════════════════════════════════════════════════
# MODULE 1: 결과 자동 수집
# ═══════════════════════════════════════════════════════════════
class ResultsCollector:
    """모든 분석 결과를 읽어 구조화된 컨텍스트로 변환"""

    def __init__(self, result_dir: Path):
        self.result_dir = result_dir

    def load_all(self) -> dict:
        ctx = {}
        ctx['mrna']       = self._load_mrna()
        ctx['proteomics'] = self._load_proteomics()
        ctx['metabolomics']= self._load_metabolomics()
        ctx['mofa']       = self._load_mofa()
        ctx['diablo']     = self._load_diablo()
        ctx['ml']         = self._load_ml()
        ctx['ppi']        = self._load_ppi()
        ctx['immune']     = self._load_immune()
        ctx['inflam']     = self._load_inflammatome()
        return ctx

    def _safe_read(self, path, top_n=20):
        try:
            df = pd.read_csv(path)
            return df.head(top_n).to_dict('records')
        except Exception:
            return []

    def _load_mrna(self) -> dict:
        base = self.result_dir / "mrna"
        return {
            "top_deg": self._safe_read(base / "DESeq2_significant_DEG.csv", 30),
            "gsea_kegg": self._safe_read(base / "GSEA_KEGG_results.csv", 15),
            "wgcna_hubs": self._safe_read(base / "WGCNA_hub_genes.csv", 20),
            "summary": self._safe_read(base / "mrna_analysis_summary.csv", 1)
        }

    def _load_proteomics(self) -> dict:
        base = self.result_dir / "proteomics"
        return {
            "top_dep":     self._safe_read(base / "limma_significant_DEP.csv", 30),
            "wgcna_hubs":  self._safe_read(base / "WGCNA_hub_proteins.csv", 20),
            "concordance": self._safe_read(base / "mRNA_Protein_concordance.csv", 50)
        }

    def _load_metabolomics(self) -> dict:
        base = self.result_dir / "metabolomics"
        return {
            "top_metabolites": self._safe_read(base / "metabolomics_significant.csv", 30),
            "pathway":         self._safe_read(base / "metabolite_pathway_enrichment.csv", 15),
            "summary":         self._safe_read(base / "metabolomics_summary.csv", 1)
        }

    def _load_mofa(self) -> dict:
        base = self.result_dir / "mofa"
        return {
            "factor_corr": self._safe_read(base / "mofa_factor_phenotype_corr.csv", 15)
        }

    def _load_diablo(self) -> dict:
        base = self.result_dir / "diablo"
        return {
            "selected_features": self._safe_read(base / "diablo_all_selected_features.csv", 30),
            "ber":               self._safe_read(base / "diablo_loocv_BER.csv", 5)
        }

    def _load_ml(self) -> dict:
        base = self.result_dir / "ml"
        return {
            "final_panel":  self._safe_read(base / "final_biomarker_panel.csv", 15),
            "auc_summary":  self._safe_read(base / "final_auc_summary.csv", 10),
            "lasso_features": self._safe_read(base / "lasso_selected_features.csv", 20)
        }

    def _load_ppi(self) -> dict:
        base = self.result_dir / "ppi"
        return {
            "hub_genes":    self._safe_read(base / "hub_gene_centrality.csv", 20),
            "cross_hubs_l4": self._safe_read(base / "cross_hub_level4.csv", 15),
            "cross_hubs_l5": self._safe_read(base / "cross_hub_level5.csv", 10)
        }

    def _load_immune(self) -> dict:
        base = self.result_dir / "immune"
        return {
            "cell_stats": self._safe_read(base / "CIBERSORT_stats.csv", 15),
            "hub_corr":   self._safe_read(base / "hub_gene_immune_correlation.csv", 10)
        }

    def _load_inflammatome(self) -> dict:
        base = self.result_dir / "inflammatome"
        return {
            "specific_deg": self._safe_read(
                base / "periodontitis_specific_DEG_top100.csv", 20),
            "scores":       self._safe_read(base / "ssgsea_inflammation_scores.csv", 5)
        }


# ═══════════════════════════════════════════════════════════════
# MODULE 2: PubMed 자동 검색
# ═══════════════════════════════════════════════════════════════
class PubMedSearcher:
    """바이오마커 후보에 대한 PubMed 문헌 검색"""

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search_gene_periodontitis(self, gene: str, max_results: int = 3) -> list:
        """특정 유전자 + 치주염 관련 논문 검색"""
        query = f"{gene}[Title/Abstract] AND periodontitis[Title/Abstract]"
        return self._search(query, max_results)

    def search_pathway(self, pathway: str, max_results: int = 3) -> list:
        """경로 + 치주염 검색"""
        query = f"{pathway}[Title/Abstract] AND periodontitis[Title/Abstract]"
        return self._search(query, max_results)

    def _search(self, query: str, max_results: int) -> list:
        try:
            # ESearch
            search_url = f"{self.BASE_URL}/esearch.fcgi"
            search_params = {
                "db": "pubmed", "term": query,
                "retmax": max_results, "retmode": "json",
                "sort": "relevance"
            }
            r = requests.get(search_url, params=search_params, timeout=10)
            ids = r.json().get("esearchresult", {}).get("idlist", [])

            if not ids:
                return []

            # EFetch — 초록 가져오기
            fetch_url = f"{self.BASE_URL}/efetch.fcgi"
            fetch_params = {
                "db": "pubmed", "id": ",".join(ids),
                "rettype": "abstract", "retmode": "xml"
            }
            rf = requests.get(fetch_url, params=fetch_params, timeout=15)
            root = ET.fromstring(rf.content)

            results = []
            for article in root.findall(".//PubmedArticle"):
                title = article.findtext(".//ArticleTitle", "")
                year  = article.findtext(".//PubDate/Year", "")
                pmid  = article.findtext(".//PMID", "")
                abstract = article.findtext(".//AbstractText", "")[:300] if article.findtext(".//AbstractText") else ""
                results.append({
                    "pmid": pmid, "title": title,
                    "year": year, "abstract": abstract[:300]
                })

            time.sleep(0.35)  # NCBI rate limit
            return results

        except Exception as e:
            return [{"error": str(e)}]


# ═══════════════════════════════════════════════════════════════
# MODULE 3: LLM 분석 에이전트
# ═══════════════════════════════════════════════════════════════
class LLMAnalysisAgent:
    """Claude API 기반 멀티오믹스 결과 해석 에이전트"""

    def __init__(self, model: str = DEFAULT_MODEL):
        self.model   = model
        self.api_key = os.environ.get("ANTHROPIC_API_KEY") or \
                       os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError(
                "API key 없음. 환경변수 설정:\n"
                "  export ANTHROPIC_API_KEY='sk-ant-...'\n"
                "  또는 export OPENAI_API_KEY='sk-...'"
            )
        self.provider = "anthropic" if "ANTHROPIC_API_KEY" in os.environ else "openai"
        self.history  = []   # 대화 이력 유지

    def chat(self, user_message: str, system: str = None) -> str:
        """LLM에 메시지 전송 + 응답 반환"""
        self.history.append({"role": "user", "content": user_message})

        if self.provider == "anthropic":
            return self._call_anthropic(system)
        else:
            return self._call_openai(system)

    def _call_anthropic(self, system: str = None) -> str:
        headers = {
            "x-api-key":         self.api_key,
            "anthropic-version": "2023-06-01",
            "content-type":      "application/json"
        }
        payload = {
            "model":      self.model,
            "max_tokens": 4096,
            "messages":   self.history
        }
        if system:
            payload["system"] = system

        r = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers=headers, json=payload, timeout=120
        )
        r.raise_for_status()
        response = r.json()["content"][0]["text"]
        self.history.append({"role": "assistant", "content": response})
        return response

    def _call_openai(self, system: str = None) -> str:
        from openai import OpenAI
        client = OpenAI(api_key=self.api_key)
        messages = self.history.copy()
        if system:
            messages.insert(0, {"role": "system", "content": system})

        r = client.chat.completions.create(
            model=self.model, messages=messages,
            max_tokens=4096
        )
        response = r.choices[0].message.content
        self.history.append({"role": "assistant", "content": response})
        return response

    def reset_history(self):
        self.history = []


# ═══════════════════════════════════════════════════════════════
# MODULE 4: 분석 워크플로우
# ═══════════════════════════════════════════════════════════════
class MultiOmicsInsightEngine:
    """전체 분석 → LLM 해석 → 보고서 생성 통합 엔진"""

    SYSTEM_PROMPT = """당신은 멀티오믹스 생물정보학 전문가이자 치주의학 연구자입니다.

역할:
- 치주염(Periodontitis) 멀티오믹스(mRNA + DDA Proteomics + Untargeted Metabolomics) 분석 결과를 해석
- 바이오마커 후보의 생물학적 의미와 임상적 가치를 평가
- 기존 문헌과 비교하여 새로운 인사이트를 발굴
- 논문 Discussion 초안 작성 지원

분석 맥락:
- 샘플: 정상 30명 vs 치주염 30명 (n=60, 동일 조직에서 3개 오믹스 동시 측정)
- 치주염: 만성 염증성 치주 조직 파괴 질환
- 목표: 진단 바이오마커 패널 발굴

답변 원칙:
- 한국어로 답변 (기술 용어는 영어 병기)
- 구체적 수치와 문헌 근거 포함
- 단순 요약이 아닌 "왜 중요한가" 깊이 있는 해석
- 실제 임상 적용 가능성 평가
- 불확실한 부분은 솔직하게 언급"""

    def __init__(self, model: str = DEFAULT_MODEL):
        self.collector = ResultsCollector(RESULT_DIR)
        self.pubmed    = PubMedSearcher()
        self.agent     = LLMAnalysisAgent(model)
        self.ctx       = {}
        self.report    = []

    def _add_section(self, title: str, content: str):
        self.report.append(f"\n## {title}\n")
        self.report.append(content)
        print(f"\n{'='*60}\n{title}\n{'='*60}")
        print(content[:500] + "..." if len(content) > 500 else content)

    def _format_ctx(self, data: dict, max_items: int = 10) -> str:
        """dict → 간결한 텍스트 변환"""
        lines = []
        for k, v in data.items():
            if isinstance(v, list) and v:
                items = v[:max_items]
                lines.append(f"**{k}** ({len(v)}건):")
                for item in items[:5]:
                    if isinstance(item, dict):
                        # 주요 컬럼만 추출
                        summary = {kk: vv for kk, vv in item.items()
                                   if kk in ['gene','protein','metabolite',
                                              'log2FoldChange','logFC','padj',
                                              'adj.P.Val','VIP','FC',
                                              'direction','module','hub_score',
                                              'Description','NES','mean_shap']}
                        lines.append(f"  - {summary}")
        return "\n".join(lines) if lines else "데이터 없음"


    # ── 분석 1: 단일 오믹스 결과 해석 ──────────────────────
    def analyze_single_omics(self):
        print("\n[1/7] 단일 오믹스 결과 해석 중...")
        self.agent.reset_history()

        prompt = f"""
아래는 치주염(n=30) vs 정상(n=30) 비교 멀티오믹스 분석 결과입니다.

=== mRNA 결과 ===
{self._format_ctx(self.ctx['mrna'])}

=== Proteomics 결과 ===
{self._format_ctx(self.ctx['proteomics'])}

=== Metabolomics 결과 ===
{self._format_ctx(self.ctx['metabolomics'])}

위 결과를 바탕으로 다음을 분석해주세요:

1. **핵심 발견 요약** (오믹스별 top 5 마커와 그 의미)
2. **생물학적 해석** — 이 마커들이 치주염 병인에서 하는 역할은?
3. **예상치 못한 발견** — 기존 문헌과 다르거나 새로운 것은?
4. **각 오믹스의 강점** — 어느 오믹스가 가장 강한 신호를 보이는가?
"""
        response = self.agent.chat(prompt, self.SYSTEM_PROMPT)
        self._add_section("1. 단일 오믹스 결과 해석", response)


    # ── 분석 2: mRNA-Protein Concordance 해석 ──────────────
    def analyze_concordance(self):
        print("\n[2/7] mRNA-Protein Concordance 해석 중...")

        conc_data = self.ctx['proteomics'].get('concordance', [])
        concordant = [x for x in conc_data if x.get('type') == 'Concordant']
        discordant = [x for x in conc_data if x.get('type') == 'Discordant']

        prompt = f"""
mRNA-Protein Concordance 분석 결과:

**Concordant 유전자** (mRNA↑ & Protein↑ 또는 양방향 일치, {len(concordant)}개):
{json.dumps(concordant[:10], indent=2, ensure_ascii=False)}

**Discordant 유전자** (mRNA와 Protein 방향 불일치, {len(discordant)}개):
{json.dumps(discordant[:10], indent=2, ensure_ascii=False)}

분석 요청:
1. **Concordant 마커의 임상적 가치** — 왜 양방향 일치가 가장 강력한 바이오마커인가?
2. **Discordant 마커 해석** — 불일치의 생물학적 의미는?
   (Post-translational modification, protein stability, ubiquitination 등)
3. **치주염에서 특히 주목할 Discordant 유전자** — 어떤 단백질 조절 메커니즘이 관여할까?
4. **바이오마커 선택 전략** — Concordant vs Discordant 중 어느 것을 우선시해야 하나?
"""
        response = self.agent.chat(prompt)
        self._add_section("2. mRNA-Protein Concordance 심층 해석", response)


    # ── 분석 3: Cross-omics 패턴 → 새 가설 발굴 ───────────
    def discover_cross_omics_patterns(self):
        print("\n[3/7] Cross-omics 패턴 탐색 중...")

        # 각 오믹스 상위 마커 통합
        top_genes   = [x.get('gene','') for x in self.ctx['mrna'].get('top_deg', [])[:15]]
        top_prots   = [x.get('protein','') for x in self.ctx['proteomics'].get('top_dep', [])[:15]]
        top_metabs  = [x.get('metabolite','') for x in self.ctx['metabolomics'].get('top_metabolites', [])[:15]]
        hub_genes   = [x.get('gene','') for x in self.ctx['ppi'].get('cross_hubs_l4', [])[:10]]
        shap_panel  = [x.get('gene_name','') for x in self.ctx['ml'].get('final_panel', [])[:10]]
        mofa_sig    = [x.get('factor','') for x in self.ctx['mofa'].get('factor_corr', [])
                       if x.get('padj', 1) < 0.05]

        prompt = f"""
멀티오믹스 통합 분석 결과 요약:

**mRNA Top DEG:** {', '.join(filter(None, top_genes))}
**Protein Top DEP:** {', '.join(filter(None, top_prots))}
**Metabolite Top 마커:** {', '.join(filter(None, top_metabs))}
**PPI Hub (다중 증거):** {', '.join(filter(None, hub_genes))}
**ML 최종 패널:** {', '.join(filter(None, shap_panel))}
**MOFA 유의 Factor:** {', '.join(map(str, mofa_sig))}

면역세포 변화: {self._format_ctx(self.ctx['immune'], 5)}
Inflammatome 특이 마커: {self._format_ctx(self.ctx['inflam'], 5)}

이 정보를 종합하여 분석해주세요:

1. **가장 강력한 통합 바이오마커 후보** (3-5개)
   - 여러 오믹스에서 동시에 나타나는 것
   - 근거: [오믹스 레이어], [증거 강도], [예상 AUC 기여]

2. **예상치 못한 Cross-omics 패턴**
   - 기존에 알려지지 않은 mRNA-대사체 연결
   - 면역세포-분자마커 연관성
   - 새로운 생물학적 가설

3. **치주염 분자 서브타입 가능성**
   - 데이터가 시사하는 환자 이질성
   - 어떤 분자 기반으로 서브타입 구분 가능한가?

4. **시스템 생물학적 해석**
   - MOFA factor가 나타내는 질병 축은 무엇인가?
   - Gene-Protein-Metabolite 조절 축 제안
"""
        response = self.agent.chat(prompt)
        self._add_section("3. Cross-omics 패턴 & 새로운 가설", response)


    # ── 분석 4: 바이오마커별 PubMed 근거 검색 + LLM 해석 ──
    def research_biomarkers_with_literature(self):
        print("\n[4/7] 바이오마커 문헌 검색 + 해석 중...")

        # SHAP top 패널 마커 대상
        panel = self.ctx['ml'].get('final_panel', [])[:8]
        if not panel:
            panel = self.ctx['ppi'].get('cross_hubs_l4', [])[:8]

        lit_results = {}
        for item in panel:
            gene = item.get('gene_name') or item.get('gene') or item.get('feature', '')
            gene = str(gene).replace('mRNA_','').replace('prot_','').replace('metab_','')
            if gene and gene != 'nan':
                print(f"  PubMed 검색: {gene} + periodontitis")
                lit_results[gene] = self.pubmed.search_gene_periodontitis(gene, max_results=3)

        # 문헌 요약 → LLM
        lit_summary = json.dumps(lit_results, indent=2, ensure_ascii=False)[:6000]

        prompt = f"""
최종 바이오마커 패널과 PubMed 문헌 검색 결과:

**ML 최종 패널:**
{json.dumps(panel, indent=2, ensure_ascii=False)}

**문헌 검색 결과 (각 마커별 관련 논문):**
{lit_summary}

분석 요청:
1. **문헌 근거 평가** — 각 마커에 대한 기존 연구 현황
   - 치주염에서 이미 알려진 것 vs 새로 발견된 것
   - 다른 염증 질환에서의 역할 (근거 이전 가능성)

2. **새로운 바이오마커 인사이트**
   - 기존에 치주염 바이오마커로 제안된 적 없는 분자는?
   - 왜 이번 연구에서 처음 나타났을 가능성이 높은가?
   (3-omics 매칭 설계의 독특성?)

3. **임상 적용 가능성 평가**
   - 어느 마커가 진단 키트로 가장 적합한가?
   - 샘플 종류 (혈청/타액/치은열구액) 전환 가능성?

4. **추가 검증 실험 제안**
   - 어떤 마커를 우선적으로 in vitro / in vivo 검증해야 하나?
   - 가장 효율적인 검증 전략은?
"""
        response = self.agent.chat(prompt)
        self._add_section("4. 바이오마커 문헌 근거 + 임상 해석", response)


    # ── 분석 5: Discussion 초안 생성 ───────────────────────
    def generate_discussion_draft(self):
        print("\n[5/7] Discussion 초안 생성 중...")
        self.agent.reset_history()  # 새 대화 시작

        # 지금까지의 핵심 발견을 요약해서 전달
        key_findings = "\n".join(self.report[-2:]) if len(self.report) >= 2 else "앞 분석 참고"

        prompt = f"""
치주염 멀티오믹스 연구의 논문 Discussion 섹션 초안을 작성해주세요.

**연구 설계:**
- 정상 30 vs 치주염 30 (동일 조직 3개 오믹스)
- mRNA (DESeq2), DDA Proteomics (limma), Untargeted Metabolomics (PLS-DA)
- MOFA+, DIABLO 통합, Stacking Ensemble ML + SHAP
- Inflammatome 분류, PPI 네트워크

**핵심 발견 (이전 분석에서 도출):**
{key_findings[:3000]}

**참고 선행 연구:**
- Luo et al. 2023 (Arch Oral Biol): transcriptomics+metabolomics → 퓨린대사, ABC transporter
- Chu et al. 2024 (J Proteome Res): gingival metabolomics → 포도당/퓨린/아미노산 대사 이상
- Front. Med. 2025: transcriptomics+methylation → 9개 진단 유전자 (AUC 높음)

Discussion 구성 (각 소절 3-5문장):

### 4.1 멀티오믹스 통합의 가치
### 4.2 핵심 바이오마커의 생물학적 의미
### 4.3 면역 microenvironment 해석
### 4.4 Inflammatome 적용의 의의
### 4.5 임상 적용 가능성과 한계
### 4.6 향후 연구 방향

각 소절은 영어로 작성 (논문 투고용), 주요 용어는 한글 해설 병기.
"""
        response = self.agent.chat(prompt, self.SYSTEM_PROMPT)
        self._add_section("5. Discussion 초안 (영문)", response)


    # ── 분석 6: 새로운 바이오마커 가설 심층 발굴 ──────────
    def deep_biomarker_hypothesis(self):
        print("\n[6/7] 새로운 바이오마커 가설 심층 발굴 중...")

        # Inflammatome 특이 마커 (완전히 새로운 후보)
        specific_deg = self.ctx['inflam'].get('specific_deg', [])[:15]
        discordant   = [x for x in self.ctx['proteomics'].get('concordance', [])
                        if x.get('type') == 'Discordant'][:10]

        prompt = f"""
치주염 특이 바이오마커 가설 발굴을 위한 심층 분석:

**Inflammatome 제외 후 치주염 특이 DEG** (가장 독창적 후보):
{json.dumps(specific_deg, indent=2, ensure_ascii=False)}

**mRNA-Protein Discordant 유전자** (PTM 조절 가능성):
{json.dumps(discordant, indent=2, ensure_ascii=False)}

**PPI 다중증거 Hub**:
{self._format_ctx({'hubs': self.ctx['ppi'].get('cross_hubs_l4', [])[:10]})}

다음에 대해 깊이 있는 가설을 제시해주세요:

1. **가장 독창적인 바이오마커 후보 TOP 3**
   각 후보에 대해:
   - 왜 기존 연구에서 놓쳤을 가능성이 있는가?
   - 치주염 병인에서의 역할 메커니즘 가설
   - 어떤 실험으로 가설을 검증할 수 있나?

2. **Post-translational regulation 가설**
   Discordant 마커 중 가장 흥미로운 것:
   - 단백질이 mRNA와 반대 방향인 이유
   - 관련 E3 ubiquitin ligase / kinase 예측
   - 임상적 의미 (단백질 레벨 측정의 중요성)

3. **치주염 분자 서브타입 가설**
   - 데이터가 시사하는 2가지 이상의 분자 서브타입
   - 각 서브타입의 특징적 마커
   - 서브타입별 치료 전략 차별화 가능성

4. **구강-전신 연관성 바이오마커**
   - 치주염 마커 중 심혈관/당뇨와 연관 가능성 있는 것
   - 치주염이 전신 염증에 미치는 영향을 설명하는 분자
"""
        response = self.agent.chat(prompt)
        self._add_section("6. 새로운 바이오마커 가설 심층 발굴", response)


    # ── 분석 7: 결론 + 다음 단계 ───────────────────────────
    def generate_conclusion_and_next_steps(self):
        print("\n[7/7] 결론 + 다음 단계 도출 중...")

        prompt = f"""
지금까지의 분석을 종합하여:

1. **핵심 결론 3줄 요약**
   - 이 연구가 치주염 분야에 기여하는 핵심 내용
   - 기존 연구와의 차별점

2. **최종 추천 바이오마커 패널**
   형식:
   | 순위 | 마커명 | 오믹스 | 증거강도 | 임상의미 | 우선검증 |
   |------|--------|--------|---------|---------|--------|
   최소 6개, 최대 10개

3. **즉시 할 수 있는 추가 분석** (데이터 재활용)
   - 현재 데이터로 추가로 끌어낼 수 있는 분석
   - 예: 치주염 중증도별 용량-반응 분석, 상관 네트워크 시각화 등

4. **논문 accept 가능성 높이는 전략**
   - 추가하면 좋을 분석 (외부 검증 등)
   - 어느 저널에 투고할 때 어떤 점을 강조해야 하나

5. **향후 연구 제안** (후속 논문 씨앗)
   - 이 연구 결과로 설계할 수 있는 후속 연구 2-3가지
"""
        response = self.agent.chat(prompt)
        self._add_section("7. 결론 + 다음 단계 제안", response)


    # ── 전체 실행 ────────────────────────────────────────────
    def run(self):
        print("=" * 60)
        print("LLM Multi-omics Insight Engine 시작")
        print("=" * 60)

        # 결과 로드
        print("\n[0/7] 분석 결과 로드 중...")
        self.ctx = self.collector.load_all()

        # 보고서 헤더
        self.report.append(f"# 치주염 멀티오믹스 LLM 인사이트 보고서\n")
        self.report.append(f"**생성일시:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.report.append(f"**모델:** {self.agent.model}\n")
        self.report.append("---\n")

        # 7개 분석 순차 실행
        self.analyze_single_omics()
        self.analyze_concordance()
        self.discover_cross_omics_patterns()
        self.research_biomarkers_with_literature()
        self.generate_discussion_draft()
        self.deep_biomarker_hypothesis()
        self.generate_conclusion_and_next_steps()

        # 보고서 저장
        report_text = "\n".join(self.report)
        with open(REPORT_PATH, 'w', encoding='utf-8') as f:
            f.write(report_text)

        print(f"\n{'='*60}")
        print(f"✅ 완료! 보고서: {REPORT_PATH}")
        print(f"{'='*60}")
        return report_text


# ═══════════════════════════════════════════════════════════════
# MODULE 5: 인터랙티브 대화 모드
# ═══════════════════════════════════════════════════════════════
class InteractiveChatMode:
    """연구자와 LLM이 결과에 대해 자유롭게 대화"""

    def __init__(self, model: str = DEFAULT_MODEL):
        self.agent     = LLMAnalysisAgent(model)
        self.collector = ResultsCollector(RESULT_DIR)
        self.ctx       = {}

    def start(self):
        print("=" * 60)
        print("인터랙티브 분석 모드 (종료: 'quit' 또는 'exit')")
        print("=" * 60)
        print("\n분석 결과 로드 중...")
        self.ctx = self.collector.load_all()

        # 컨텍스트 요약을 시스템에 주입
        ctx_summary = self._build_context_summary()

        system = MultiOmicsInsightEngine.SYSTEM_PROMPT + f"""

=== 현재 분석 결과 요약 ===
{ctx_summary}

이 데이터를 바탕으로 연구자의 질문에 답하세요.
대화 이력을 유지하며 이전 답변을 참고하세요.
"""
        print("\n✅ 준비 완료. 분석 결과에 대해 자유롭게 질문하세요.\n")
        print("예시 질문:")
        print("  - 'IL6이 왜 중요한가?'")
        print("  - '퓨린 대사 이상과 치주염의 관계를 설명해줘'")
        print("  - '바이오마커 패널 중 타액에서 측정 가능한 것은?'")
        print("  - 'Discussion 4.2 소절 더 자세히 써줘'")
        print("  - '이 데이터로 쓸 수 있는 후속 논문 아이디어는?'\n")

        while True:
            try:
                user_input = input("🔬 질문: ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\n대화 종료.")
                break

            if user_input.lower() in ('quit', 'exit', '종료', 'q'):
                print("대화를 종료합니다.")
                break
            if not user_input:
                continue

            print("\n🤖 분석 중...\n")
            response = self.agent.chat(user_input, system if not self.agent.history else None)
            print(f"답변:\n{response}\n")
            print("-" * 60)

            # 대화 저장
            chat_log_path = OUTDIR / f"chat_{datetime.now().strftime('%Y%m%d')}.md"
            with open(chat_log_path, 'a', encoding='utf-8') as f:
                f.write(f"\n**Q:** {user_input}\n\n**A:** {response}\n\n---\n")

    def _build_context_summary(self) -> str:
        lines = []
        # DEG top 10
        degs = self.ctx['mrna'].get('top_deg', [])[:10]
        if degs:
            genes = [f"{x.get('gene','')}({x.get('direction','')})" for x in degs]
            lines.append(f"Top DEG: {', '.join(genes)}")
        # DEP top 10
        deps = self.ctx['proteomics'].get('top_dep', [])[:10]
        if deps:
            prots = [f"{x.get('protein','')}({x.get('direction','')})" for x in deps]
            lines.append(f"Top DEP: {', '.join(prots)}")
        # 대사체 top 10
        mets = self.ctx['metabolomics'].get('top_metabolites', [])[:10]
        if mets:
            ms = [f"{x.get('metabolite','')}" for x in mets]
            lines.append(f"Top 대사체: {', '.join(ms)}")
        # ML 패널
        panel = self.ctx['ml'].get('final_panel', [])[:8]
        if panel:
            ps = [x.get('gene_name','') for x in panel]
            lines.append(f"ML 바이오마커 패널: {', '.join(filter(None, ps))}")
        return "\n".join(lines) if lines else "결과 파일 없음 (분석 먼저 실행)"


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description="LLM 기반 멀티오믹스 인사이트 분석 에이전트"
    )
    parser.add_argument(
        "--mode", choices=["auto", "chat"], default="auto",
        help="auto: 자동 보고서 생성 | chat: 인터랙티브 대화"
    )
    parser.add_argument(
        "--model", default=DEFAULT_MODEL,
        help=f"LLM 모델 (기본: {DEFAULT_MODEL})"
    )
    args = parser.parse_args()

    if args.mode == "auto":
        engine = MultiOmicsInsightEngine(model=args.model)
        engine.run()
    else:
        chat = InteractiveChatMode(model=args.model)
        chat.start()


if __name__ == "__main__":
    main()
