# 🤖 LLM 분석 에이전트 사용 가이드

## 개요

`11_llm_analysis_agent.py`는 멀티오믹스 분석 결과를 Claude/GPT-4에게 전달하여  
**생물학적 해석 → 바이오마커 가설 발굴 → Discussion 초안 → 결론**을  
자동으로 생성하는 AI 연구 보조 에이전트입니다.

---

## 🔄 작동 방식

```
[01~10 파이프라인 결과]
  results/mrna/     → DEG, GSEA, WGCNA hub
  results/prot/     → DEP, concordance
  results/metab/    → 유의 대사체
  results/mofa/     → Factor 상관
  results/ml/       → 최종 패널, SHAP, AUC
  results/ppi/      → Hub gene
  results/immune/   → 면역세포 비율
  results/inflam/   → 특이 DEG
         ↓
[자동 수집 + 구조화]
         ↓
[PubMed 자동 검색]
  각 바이오마커 후보 → 관련 논문 검색 (NCBI Entrez API)
         ↓
[Claude / GPT-4 추론]
  7단계 순차 분석:
  ① 단일 오믹스 해석
  ② mRNA-Protein Concordance
  ③ Cross-omics 패턴 & 새 가설
  ④ 문헌 근거 + 임상 해석
  ⑤ Discussion 초안 (영문)
  ⑥ 새로운 바이오마커 가설 심층
  ⑦ 결론 + 다음 단계
         ↓
[Markdown 보고서 자동 저장]
  results/llm_insights/biomarker_insights_YYYYMMDD_HHMMSS.md
```

---

## 🚀 실행 방법

### 1. API Key 설정

```bash
# Claude (Anthropic) — 권장
export ANTHROPIC_API_KEY="sk-ant-api03-..."

# 또는 OpenAI
export OPENAI_API_KEY="sk-..."
```

### 2. 실행 모드

#### 모드 A: 자동 보고서 생성 (전체 7단계)
```bash
cd pipeline
python 11_llm_analysis_agent.py --mode auto
# 소요 시간: 약 5~15분
# 출력: results/llm_insights/biomarker_insights_YYYYMMDD.md
```

#### 모드 B: 인터랙티브 대화 (자유 질문)
```bash
python 11_llm_analysis_agent.py --mode chat
```

### 3. 모델 선택

```bash
# Claude Opus (가장 깊은 추론, 권장)
python 11_llm_analysis_agent.py --model claude-opus-4-5

# Claude Sonnet (빠른 속도, 비용 절감)
python 11_llm_analysis_agent.py --model claude-sonnet-4-5

# GPT-4o (OpenAI 사용 시)
python 11_llm_analysis_agent.py --model gpt-4o
```

---

## 💬 인터랙티브 모드 예시 질문

### 바이오마커 해석
```
🔬 "IL6가 왜 가장 중요한 바이오마커인가?"
🔬 "MMP8과 MMP13의 역할 차이는?"
🔬 "퓨린 대사 이상이 치주염에서 어떤 의미인가?"
🔬 "deoxyinosine이 바이오마커로 적합한 이유는?"
```

### 새로운 가설 발굴
```
🔬 "mRNA↑ 단백질↓인 유전자 중 가장 흥미로운 것은?"
🔬 "치주염 면역 microenvironment에서 Treg 감소의 임상적 의미는?"
🔬 "MOFA Factor1이 나타내는 생물학적 의미는 무엇인가?"
🔬 "inflammatome에 없는 치주염 특이 마커 중 새로운 것을 분석해줘"
```

### 논문 작성 지원
```
🔬 "Discussion 4.2 소절을 더 자세히 써줘"
🔬 "이 결과를 Journal of Clinical Periodontology에 맞게 강조점을 잡아줘"
🔬 "Abstract를 250단어 이내로 써줘"
🔬 "Reviewer 예상 질문과 답변을 준비해줘"
```

### 후속 연구
```
🔬 "이 데이터로 쓸 수 있는 후속 논문 아이디어 3가지는?"
🔬 "바이오마커 검증을 위한 최소 비용 실험 설계는?"
🔬 "타액 기반 진단 키트로 전환 가능한 마커는?"
🔬 "치주염-당뇨 연관성을 이 데이터로 분석할 수 있는가?"
```

---

## 📊 자동 생성 보고서 구조

```markdown
# 치주염 멀티오믹스 LLM 인사이트 보고서

## 1. 단일 오믹스 결과 해석
   - 핵심 발견 요약 (오믹스별 top 5 마커)
   - 생물학적 해석
   - 예상치 못한 발견

## 2. mRNA-Protein Concordance 심층 해석
   - Concordant 마커의 임상적 가치
   - Discordant 마커 해석 (PTM 조절)
   - 바이오마커 선택 전략

## 3. Cross-omics 패턴 & 새로운 가설
   - 가장 강력한 통합 바이오마커 후보
   - 예상치 못한 Cross-omics 패턴
   - 치주염 분자 서브타입 가능성

## 4. 바이오마커 문헌 근거 + 임상 해석
   - 마커별 PubMed 논문 근거
   - 새로운 바이오마커 vs 기존 알려진 것
   - 임상 적용 가능성 평가

## 5. Discussion 초안 (영문)
   - 4.1 멀티오믹스 통합의 가치
   - 4.2 핵심 바이오마커의 생물학적 의미
   - 4.3 면역 microenvironment 해석
   - 4.4 Inflammatome 적용의 의의
   - 4.5 임상 적용 가능성과 한계
   - 4.6 향후 연구 방향

## 6. 새로운 바이오마커 가설 심층 발굴
   - 독창적 바이오마커 후보 TOP 3
   - Post-translational regulation 가설
   - 치주염 분자 서브타입 가설
   - 구강-전신 연관성 바이오마커

## 7. 결론 + 다음 단계 제안
   - 핵심 결론 3줄 요약
   - 최종 추천 바이오마커 패널 (표)
   - 즉시 할 수 있는 추가 분석
   - 논문 accept 전략
   - 향후 연구 제안
```

---

## ⚙️ 고급 사용법

### 특정 섹션만 실행
```python
from pipeline.llm_analysis_agent import MultiOmicsInsightEngine

engine = MultiOmicsInsightEngine()
engine.ctx = engine.collector.load_all()

# 특정 분석만 실행
engine.analyze_concordance()          # mRNA-Protein 분석
engine.deep_biomarker_hypothesis()    # 가설 발굴만
engine.generate_discussion_draft()    # Discussion 초안만
```

### 대화 이력 기반 심화 질문
```python
from pipeline.llm_analysis_agent import LLMAnalysisAgent

agent = LLMAnalysisAgent()
# 1차 질문
r1 = agent.chat("MMP8의 역할을 설명해줘")
# 대화 이력 유지된 상태로 심화 질문
r2 = agent.chat("그렇다면 MMP8을 억제하는 치료 전략은?")
r3 = agent.chat("현재 임상에서 MMP8 억제제가 사용되고 있나?")
```

---

## 📌 주의사항

- API 비용: Claude Opus 기준 전체 auto 실행 시 약 $0.5~2 소요
- 인터넷 연결 필요 (PubMed API + LLM API)
- 분석 결과 파일 (`results/` 디렉토리)이 있어야 정확한 해석 가능
- 결과 파일 없이도 실행 가능하나 LLM이 일반적 답변 제공
