-- ================================================================
-- 02_db_schema.sql
-- 통합 오믹스 DB 스키마 — PostgreSQL
-- 자체 연구 + 공개 DB 통합
-- ================================================================

-- 확장
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS pg_trgm;    -- 텍스트 검색

-- ================================================================
-- 1. 데이터 소스 등록
-- ================================================================
CREATE TABLE data_sources (
    id          SERIAL PRIMARY KEY,
    source_name VARCHAR(50)  NOT NULL,  -- GEO, PRIDE, MetaboLights, OWN
    source_type VARCHAR(30),            -- public, proprietary
    accession   VARCHAR(50),
    title       TEXT,
    url         TEXT,
    organism    VARCHAR(100) DEFAULT 'Homo sapiens',
    tissue      VARCHAR(100),
    disease     VARCHAR(200),
    n_samples   INTEGER,
    omics_type  VARCHAR(50),            -- mRNA, Proteomics, Metabolomics
    platform    VARCHAR(100),
    year        INTEGER,
    doi         VARCHAR(200),
    created_at  TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_sources_accession ON data_sources(accession);
CREATE INDEX idx_sources_omics ON data_sources(omics_type);
CREATE INDEX idx_sources_disease ON data_sources USING gin(disease gin_trgm_ops);

-- ================================================================
-- 2. 샘플 메타데이터
-- ================================================================
CREATE TABLE samples (
    id              UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    source_id       INTEGER REFERENCES data_sources(id),
    sample_id       VARCHAR(100) NOT NULL,  -- GSM*, 자체 ID 등
    group_label     VARCHAR(50),            -- Control, Periodontitis, Disease
    disease_stage   VARCHAR(50),            -- Stage I/II/III/IV
    age             INTEGER,
    sex             VARCHAR(10),
    smoking         BOOLEAN,
    diabetes        BOOLEAN,
    probing_depth   FLOAT,                  -- 임상 지표 (mm)
    cal             FLOAT,                  -- Clinical Attachment Level
    bop             FLOAT,                  -- Bleeding on Probing (%)
    tissue_source   VARCHAR(100),           -- gingival, saliva, GCF, blood
    metadata_json   JSONB,                  -- 기타 메타데이터
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_samples_source ON samples(source_id);
CREATE INDEX idx_samples_group ON samples(group_label);
CREATE INDEX idx_samples_disease ON samples(disease_stage);

-- ================================================================
-- 3. 유전자/단백질/대사체 참조 테이블
-- ================================================================
CREATE TABLE genes (
    id          SERIAL PRIMARY KEY,
    hgnc_symbol VARCHAR(50) UNIQUE NOT NULL,
    ensembl_id  VARCHAR(20),
    entrez_id   INTEGER,
    full_name   TEXT,
    chromosome  VARCHAR(10),
    gene_type   VARCHAR(50),           -- protein_coding, lncRNA 등
    uniprot_id  VARCHAR(20),
    updated_at  TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_genes_symbol ON genes(hgnc_symbol);
CREATE INDEX idx_genes_ensembl ON genes(ensembl_id);

CREATE TABLE proteins (
    id              SERIAL PRIMARY KEY,
    uniprot_id      VARCHAR(20) UNIQUE NOT NULL,
    gene_id         INTEGER REFERENCES genes(id),
    protein_name    TEXT,
    molecular_weight FLOAT,
    subcellular_loc TEXT[],
    go_bp           TEXT[],
    go_mf           TEXT[],
    is_secreted     BOOLEAN,           -- 분비 단백질 여부 (바이오마커 가치)
    updated_at      TIMESTAMP DEFAULT NOW()
);

CREATE TABLE metabolites (
    id              SERIAL PRIMARY KEY,
    hmdb_id         VARCHAR(20) UNIQUE,
    inchikey        VARCHAR(30),
    name            VARCHAR(200) NOT NULL,
    formula         VARCHAR(100),
    molecular_weight FLOAT,
    kegg_id         VARCHAR(20),
    chebi_id        VARCHAR(20),
    metabolite_class VARCHAR(100),     -- Amino acid, Lipid, Nucleoside 등
    pathways        TEXT[],
    biofluid        TEXT[],            -- blood, saliva, urine 등
    updated_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_metabolites_hmdb ON metabolites(hmdb_id);
CREATE INDEX idx_metabolites_name ON metabolites USING gin(name gin_trgm_ops);

-- ================================================================
-- 4. 발현 데이터 (정규화된 값)
-- ================================================================
-- 주의: 대용량 행렬은 Parquet 파일로 저장, DB에는 요약값만

CREATE TABLE mrna_expression (
    id          BIGSERIAL PRIMARY KEY,
    sample_id   UUID REFERENCES samples(id),
    gene_id     INTEGER REFERENCES genes(id),
    value       FLOAT NOT NULL,        -- VST / log2CPM normalized
    raw_count   INTEGER,
    created_at  TIMESTAMP DEFAULT NOW()
);

-- 파티셔닝 (대용량 최적화)
CREATE INDEX idx_mrna_sample ON mrna_expression(sample_id);
CREATE INDEX idx_mrna_gene ON mrna_expression(gene_id);

CREATE TABLE protein_expression (
    id          BIGSERIAL PRIMARY KEY,
    sample_id   UUID REFERENCES samples(id),
    protein_id  INTEGER REFERENCES proteins(id),
    value       FLOAT NOT NULL,        -- log2 normalized abundance
    is_imputed  BOOLEAN DEFAULT FALSE,
    created_at  TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_prot_sample ON protein_expression(sample_id);
CREATE INDEX idx_prot_protein ON protein_expression(protein_id);

CREATE TABLE metabolite_expression (
    id              BIGSERIAL PRIMARY KEY,
    sample_id       UUID REFERENCES samples(id),
    metabolite_id   INTEGER REFERENCES metabolites(id),
    value           FLOAT NOT NULL,    -- log2 normalized intensity
    vip_score       FLOAT,
    is_imputed      BOOLEAN DEFAULT FALSE,
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_metab_sample ON metabolite_expression(sample_id);
CREATE INDEX idx_metab_metabolite ON metabolite_expression(metabolite_id);

-- ================================================================
-- 5. 통계 분석 결과 저장
-- ================================================================
CREATE TABLE differential_analysis (
    id              BIGSERIAL PRIMARY KEY,
    source_id       INTEGER REFERENCES data_sources(id),
    omics_type      VARCHAR(20) NOT NULL,  -- mRNA, Protein, Metabolite
    feature_id      INTEGER,               -- gene_id / protein_id / metabolite_id
    feature_name    VARCHAR(200),
    log2fc          FLOAT,
    pvalue          FLOAT,
    padj            FLOAT,
    direction       VARCHAR(10),           -- Up, Down
    method          VARCHAR(30),           -- DESeq2, limma, t-test
    comparison      VARCHAR(100),          -- "Periodontitis_vs_Control"
    n_ctrl          INTEGER,
    n_case          INTEGER,
    analysis_date   DATE DEFAULT CURRENT_DATE
);

CREATE INDEX idx_diff_source ON differential_analysis(source_id);
CREATE INDEX idx_diff_omics ON differential_analysis(omics_type);
CREATE INDEX idx_diff_feature ON differential_analysis(feature_name);
CREATE INDEX idx_diff_padj ON differential_analysis(padj);

-- ================================================================
-- 6. 바이오마커 후보 테이블 (핵심!)
-- ================================================================
CREATE TABLE biomarker_candidates (
    id              SERIAL PRIMARY KEY,
    feature_name    VARCHAR(200) NOT NULL,
    omics_type      VARCHAR(20) NOT NULL,
    gene_id         INTEGER REFERENCES genes(id),
    protein_id      INTEGER REFERENCES proteins(id),
    metabolite_id   INTEGER REFERENCES metabolites(id),

    -- 증거 점수
    n_studies_sig   INTEGER DEFAULT 0,     -- 유의한 연구 수
    n_omics_layers  INTEGER DEFAULT 0,     -- 오믹스 레이어 수 (최대 3)
    is_ppi_hub      BOOLEAN DEFAULT FALSE,
    is_wgcna_hub    BOOLEAN DEFAULT FALSE,
    lasso_selected  BOOLEAN DEFAULT FALSE,
    shap_rank       INTEGER,               -- SHAP 랭킹

    -- 성능 지표
    meta_log2fc     FLOAT,                 -- 메타분석 평균 FC
    meta_pvalue     FLOAT,
    pooled_auc      FLOAT,                 -- 단독 AUC (메타분석)

    -- 통합 랭킹 점수 (0-100)
    biomarker_score FLOAT,
    confidence      VARCHAR(20),           -- High / Medium / Low

    -- 임상 가치
    is_secreted     BOOLEAN,
    detectable_in_saliva BOOLEAN,
    detectable_in_gcf    BOOLEAN,
    existing_assay  BOOLEAN,               -- 기존 측정법 존재 여부

    -- 문헌
    pubmed_ids      TEXT[],
    first_reported  INTEGER,               -- 최초 보고 연도

    created_at      TIMESTAMP DEFAULT NOW(),
    updated_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX idx_bm_feature ON biomarker_candidates(feature_name);
CREATE INDEX idx_bm_score ON biomarker_candidates(biomarker_score DESC);
CREATE INDEX idx_bm_omics ON biomarker_candidates(omics_type);
CREATE INDEX idx_bm_confidence ON biomarker_candidates(confidence);

-- ================================================================
-- 7. 경로 분석 결과
-- ================================================================
CREATE TABLE pathway_analysis (
    id              SERIAL PRIMARY KEY,
    source_id       INTEGER REFERENCES data_sources(id),
    omics_type      VARCHAR(20),
    pathway_id      VARCHAR(50),           -- KEGG: hsa05110, Reactome: R-HSA-xxx
    pathway_db      VARCHAR(20),           -- KEGG, Reactome, GO
    pathway_name    TEXT,
    pvalue          FLOAT,
    padj            FLOAT,
    nes             FLOAT,                 -- Normalized Enrichment Score (GSEA)
    n_genes_overlap INTEGER,
    direction       VARCHAR(20),           -- activated, suppressed
    method          VARCHAR(20)            -- GSEA, ORA
);

CREATE INDEX idx_pathway_name ON pathway_analysis USING gin(pathway_name gin_trgm_ops);

-- ================================================================
-- 8. PPI 네트워크 (STRING)
-- ================================================================
CREATE TABLE ppi_interactions (
    id              BIGSERIAL PRIMARY KEY,
    protein_a       INTEGER REFERENCES proteins(id),
    protein_b       INTEGER REFERENCES proteins(id),
    combined_score  INTEGER,               -- STRING score (0-1000)
    experimental    INTEGER,
    coexpression    INTEGER,
    textmining      INTEGER,
    database        VARCHAR(20) DEFAULT 'STRING',
    version         VARCHAR(10) DEFAULT 'v12.0'
);

CREATE INDEX idx_ppi_a ON ppi_interactions(protein_a);
CREATE INDEX idx_ppi_b ON ppi_interactions(protein_b);
CREATE INDEX idx_ppi_score ON ppi_interactions(combined_score DESC);

-- ================================================================
-- 9. 메타분석 결과
-- ================================================================
CREATE TABLE meta_analysis (
    id              SERIAL PRIMARY KEY,
    feature_name    VARCHAR(200) NOT NULL,
    omics_type      VARCHAR(20),
    n_studies       INTEGER,
    pooled_log2fc   FLOAT,
    pooled_se       FLOAT,
    pooled_pvalue   FLOAT,
    i_squared       FLOAT,                 -- 이질성 (%)
    method          VARCHAR(30),           -- random_effects, fixed_effects
    analysis_date   DATE DEFAULT CURRENT_DATE,
    studies_json    JSONB                  -- 포함된 연구 목록
);

CREATE INDEX idx_meta_feature ON meta_analysis(feature_name);
CREATE INDEX idx_meta_pvalue ON meta_analysis(pooled_pvalue);

-- ================================================================
-- 10. AI 예측 결과
-- ================================================================
CREATE TABLE ai_predictions (
    id              SERIAL PRIMARY KEY,
    feature_name    VARCHAR(200) NOT NULL,
    model_name      VARCHAR(50),           -- TransferLearning, GNN, Ensemble
    model_version   VARCHAR(20),
    predicted_auc   FLOAT,
    predicted_fc    FLOAT,
    confidence_score FLOAT,
    training_data   TEXT,                  -- 학습에 사용된 데이터셋
    prediction_date TIMESTAMP DEFAULT NOW()
);

-- ================================================================
-- 유용한 뷰 (View)
-- ================================================================

-- 바이오마커 종합 랭킹 뷰
CREATE VIEW v_biomarker_ranking AS
SELECT
    bc.feature_name,
    bc.omics_type,
    bc.biomarker_score,
    bc.confidence,
    bc.n_studies_sig,
    bc.n_omics_layers,
    bc.meta_log2fc,
    bc.pooled_auc,
    bc.is_ppi_hub,
    bc.is_wgcna_hub,
    bc.detectable_in_saliva,
    bc.detectable_in_gcf,
    g.full_name AS gene_full_name,
    p.protein_name,
    m.metabolite_class,
    bc.pubmed_ids
FROM biomarker_candidates bc
LEFT JOIN genes g ON bc.gene_id = g.id
LEFT JOIN proteins p ON bc.protein_id = p.id
LEFT JOIN metabolites m ON bc.metabolite_id = m.id
ORDER BY bc.biomarker_score DESC;

-- 다중 오믹스 공통 마커 뷰 (가장 강력한 후보)
CREATE VIEW v_multi_omics_consensus AS
SELECT
    feature_name,
    COUNT(DISTINCT omics_type) AS n_omics,
    MAX(biomarker_score)       AS max_score,
    ARRAY_AGG(DISTINCT omics_type) AS omics_types
FROM biomarker_candidates
WHERE padj_threshold < 0.05
GROUP BY feature_name
HAVING COUNT(DISTINCT omics_type) >= 2
ORDER BY n_omics DESC, max_score DESC;

-- 데이터셋 현황 뷰
CREATE VIEW v_dataset_summary AS
SELECT
    ds.source_name,
    ds.omics_type,
    COUNT(DISTINCT s.id) AS n_samples,
    SUM(CASE WHEN s.group_label='Periodontitis' THEN 1 ELSE 0 END) AS n_cases,
    SUM(CASE WHEN s.group_label='Control' THEN 1 ELSE 0 END) AS n_controls,
    ds.tissue,
    ds.year
FROM data_sources ds
LEFT JOIN samples s ON s.source_id = ds.id
GROUP BY ds.id, ds.source_name, ds.omics_type, ds.tissue, ds.year
ORDER BY ds.year DESC;
