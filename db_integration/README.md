# Global Omics DB Integration Platform

치주염 멀티오믹스 바이오마커 발굴을 위한 글로벌 공개 DB 통합 플랫폼

## 아키텍처 요약
- 수집: GEO, PRIDE, MetaboLights, TCGA, TOPMed, STRING, HMDB 등
- 큐레이션: 표준화, 품질관리, 배치 보정
- 통합 DB: PostgreSQL + Parquet
- AI 엔진: 전이학습 기반 바이오마커 예측
- 검증: 다층 교차 검증 프레임워크
