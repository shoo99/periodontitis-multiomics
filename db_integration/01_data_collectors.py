"""
01_data_collectors.py
글로벌 공개 DB 데이터 수집기
GEO / PRIDE / MetaboLights / TCGA / STRING / HMDB / GTEx
"""

import os
import re
import time
import json
import requests
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional
from datetime import datetime
import xml.etree.ElementTree as ET

RAW_DIR   = Path("../data/public_db")
META_DIR  = RAW_DIR / "metadata"
RAW_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
# 1. GEO Collector — mRNA / methylation / microarray
# ═══════════════════════════════════════════════════════════════
class GEOCollector:
    """
    NCBI GEO에서 치주염 관련 오믹스 데이터 자동 수집
    API: NCBI Entrez E-utilities
    """
    BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search_datasets(self,
                        query: str = "periodontitis gingival RNA-seq",
                        min_samples: int = 10,
                        max_results: int = 100) -> list:
        """GEO 데이터셋 검색"""
        search_url = f"{self.BASE}/esearch.fcgi"
        params = {
            "db": "gds", "term": query,
            "retmax": max_results, "retmode": "json"
        }
        r = requests.get(search_url, params=params, timeout=30)
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        print(f"  GEO 검색 결과: {len(ids)}개")

        datasets = []
        for uid in ids:
            meta = self._fetch_metadata(uid)
            if meta and meta.get("n_samples", 0) >= min_samples:
                datasets.append(meta)
            time.sleep(0.34)  # NCBI rate limit

        return datasets

    def _fetch_metadata(self, uid: str) -> Optional[dict]:
        fetch_url = f"{self.BASE}/esummary.fcgi"
        params = {"db": "gds", "id": uid, "retmode": "json"}
        try:
            r = requests.get(fetch_url, params=params, timeout=15)
            result = r.json().get("result", {}).get(uid, {})
            if not result:
                return None

            return {
                "source":    "GEO",
                "accession": result.get("accession", ""),
                "title":     result.get("title", ""),
                "n_samples": int(result.get("n_samples", 0)),
                "gpl":       result.get("gpl", ""),
                "pdat":      result.get("pdat", ""),
                "summary":   result.get("summary", "")[:300],
                "organism":  result.get("organism", ""),
                "uid":       uid
            }
        except Exception:
            return None

    def download_matrix(self, accession: str, outdir: Path = None) -> Optional[Path]:
        """GEO soft matrix 파일 다운로드"""
        if outdir is None:
            outdir = RAW_DIR / "geo" / accession
        outdir.mkdir(parents=True, exist_ok=True)

        # GEO FTP URL
        prefix = accession[:7] if len(accession) > 7 else accession
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/"
               f"{accession}/matrix/{accession}_series_matrix.txt.gz")

        outfile = outdir / f"{accession}_matrix.txt.gz"
        try:
            r = requests.get(url, stream=True, timeout=120)
            if r.status_code == 200:
                with open(outfile, 'wb') as f:
                    for chunk in r.iter_content(1024*1024):
                        f.write(chunk)
                print(f"  ✅ {accession} 다운로드 완료: {outfile}")
                return outfile
            else:
                print(f"  ⚠️  {accession} FTP 없음 (status={r.status_code})")
                return None
        except Exception as e:
            print(f"  ❌ {accession} 오류: {e}")
            return None

    def parse_matrix(self, matrix_file: Path) -> Optional[pd.DataFrame]:
        """GEO matrix 파일 → expression DataFrame"""
        import gzip
        try:
            open_fn = gzip.open if str(matrix_file).endswith('.gz') else open
            rows, meta_lines = [], []
            with open_fn(matrix_file, 'rt', errors='replace') as f:
                header_found = False
                for line in f:
                    if line.startswith("!"):
                        meta_lines.append(line.strip())
                        continue
                    if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                        header_found = True
                        cols = line.strip().split("\t")
                        continue
                    if header_found:
                        rows.append(line.strip().split("\t"))

            if not rows:
                return None
            df = pd.DataFrame(rows, columns=cols)
            df = df.set_index("ID_REF")
            df = df.apply(pd.to_numeric, errors='coerce')
            return df
        except Exception as e:
            print(f"  파싱 오류: {e}")
            return None

    def collect_periodontitis_datasets(self) -> pd.DataFrame:
        """치주염 GEO 데이터셋 전체 수집 + 메타데이터 저장"""
        queries = [
            "periodontitis gingival RNA-seq homo sapiens",
            "periodontal disease transcriptomics human gingival tissue",
            "periodontitis gene expression microarray human"
        ]
        all_meta = []
        for q in queries:
            print(f"\n검색: {q}")
            results = self.search_datasets(q, min_samples=10)
            all_meta.extend(results)

        # 중복 제거
        seen = set()
        unique = []
        for m in all_meta:
            if m['accession'] not in seen:
                seen.add(m['accession'])
                unique.append(m)

        df_meta = pd.DataFrame(unique)
        df_meta.to_csv(META_DIR / "geo_periodontitis_datasets.csv", index=False)
        print(f"\n✅ GEO 데이터셋: {len(unique)}개 확인")
        print(f"   알려진 주요 셋: GSE16134, GSE173078, GSE152042 포함 여부 확인")
        return df_meta


# ═══════════════════════════════════════════════════════════════
# 2. PRIDE Collector — DDA/DIA Proteomics
# ═══════════════════════════════════════════════════════════════
class PRIDECollector:
    """
    PRIDE Archive (EBI) — 단백체 데이터 수집
    REST API: https://www.ebi.ac.uk/pride/ws/archive/v2
    """
    BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"

    def search_projects(self,
                        keyword: str = "periodontitis",
                        page_size: int = 100) -> list:
        """PRIDE 프로젝트 검색"""
        url = f"{self.BASE}/search/projects"
        params = {
            "keyword":  keyword,
            "pageSize": page_size,
            "page":     0,
            "sortDirection": "DESC",
            "sortCondition": "submission_date"
        }
        try:
            r = requests.get(url, params=params, timeout=30)
            data = r.json()
            projects = data.get("_embedded", {}).get("projects", [])
            print(f"  PRIDE 검색 결과: {len(projects)}개 프로젝트")
            return projects
        except Exception as e:
            print(f"  PRIDE 검색 오류: {e}")
            return []

    def get_project_files(self, accession: str) -> list:
        """특정 프로젝트의 파일 목록"""
        url = f"{self.BASE}/files/byProject"
        params = {"accession": accession, "pageSize": 200}
        try:
            r = requests.get(url, params=params, timeout=30)
            files = r.json().get("_embedded", {}).get("files", [])
            # txt/tsv/csv 결과 파일만 필터
            result_files = [
                f for f in files
                if any(ext in f.get("fileName","").lower()
                       for ext in ['.txt', '.tsv', '.csv', 'result', 'protein'])
            ]
            return result_files
        except Exception:
            return []

    def download_result_file(self,
                              accession: str,
                              file_url: str,
                              filename: str) -> Optional[Path]:
        """PRIDE 결과 파일 다운로드"""
        outdir = RAW_DIR / "pride" / accession
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / filename

        try:
            r = requests.get(file_url, stream=True, timeout=300)
            if r.status_code == 200:
                with open(outfile, 'wb') as f:
                    for chunk in r.iter_content(1024*1024):
                        f.write(chunk)
                return outfile
        except Exception as e:
            print(f"  다운로드 오류: {e}")
        return None

    def collect_periodontitis_proteomics(self) -> pd.DataFrame:
        """치주염 관련 PRIDE 프로젝트 메타데이터 수집"""
        keywords = ["periodontitis", "periodontal", "gingival proteomics"]
        all_projects = []

        for kw in keywords:
            print(f"\nPRIDE 검색: {kw}")
            projects = self.search_projects(kw)
            for p in projects:
                meta = {
                    "source":      "PRIDE",
                    "accession":   p.get("accession", ""),
                    "title":       p.get("title", ""),
                    "description": p.get("projectDescription", "")[:300],
                    "n_samples":   p.get("numSamples", 0),
                    "organism":    p.get("organisms", [{}])[0].get("name", "") if p.get("organisms") else "",
                    "submit_date": p.get("submissionDate", ""),
                    "keywords":    str(p.get("keywords", []))
                }
                if meta["accession"]:
                    all_projects.append(meta)

        # 중복 제거
        df = pd.DataFrame(all_projects).drop_duplicates("accession")
        df.to_csv(META_DIR / "pride_periodontitis_projects.csv", index=False)
        print(f"\n✅ PRIDE 프로젝트: {len(df)}개")
        return df


# ═══════════════════════════════════════════════════════════════
# 3. MetaboLights Collector — Untargeted Metabolomics
# ═══════════════════════════════════════════════════════════════
class MetaboLightsCollector:
    """
    MetaboLights (EBI) — 대사체 데이터 수집
    REST API: https://www.ebi.ac.uk/metabolights/ws
    """
    BASE = "https://www.ebi.ac.uk/metabolights/ws"

    def search_studies(self, keyword: str = "periodontitis") -> list:
        """MetaboLights 연구 검색"""
        url = f"{self.BASE}/studies/search"
        params = {"query": keyword}
        try:
            r = requests.get(url, params=params, timeout=30)
            studies = r.json().get("content", [])
            print(f"  MetaboLights 검색: {len(studies)}개")
            return studies
        except Exception as e:
            print(f"  MetaboLights 오류: {e}")
            return []

    def get_study_metadata(self, accession: str) -> dict:
        """연구 메타데이터"""
        url = f"{self.BASE}/studies/{accession}/summary"
        try:
            r = requests.get(url, timeout=15)
            return r.json()
        except Exception:
            return {}

    def get_data_files(self, accession: str) -> list:
        """데이터 파일 목록"""
        url = f"{self.BASE}/studies/{accession}/files"
        try:
            r = requests.get(url, timeout=15)
            files = r.json().get("study", [])
            return [f for f in files
                    if any(ext in f.get("file","").lower()
                           for ext in ['.tsv', '.txt', 'm_', 'assay'])]
        except Exception:
            return []

    def download_maf_file(self, accession: str) -> Optional[Path]:
        """MAF (Metabolite Assignment File) 다운로드"""
        outdir = RAW_DIR / "metabolights" / accession
        outdir.mkdir(parents=True, exist_ok=True)

        files = self.get_data_files(accession)
        maf_files = [f for f in files if f.get("file","").startswith("m_")]

        for maf in maf_files[:1]:  # 첫 번째 MAF만
            file_url = (f"https://www.ebi.ac.uk/metabolights/ws/studies/"
                        f"{accession}/download?file={maf['file']}")
            outfile = outdir / maf['file']
            try:
                r = requests.get(file_url, stream=True, timeout=120)
                if r.status_code == 200:
                    with open(outfile, 'wb') as f:
                        for chunk in r.iter_content(1024*1024):
                            f.write(chunk)
                    return outfile
            except Exception:
                pass
        return None

    def collect_periodontitis_metabolomics(self) -> pd.DataFrame:
        """치주염 대사체 데이터 수집"""
        keywords = ["periodontitis", "periodontal", "gingival metabolomics"]
        all_studies = []

        for kw in keywords:
            print(f"\nMetaboLights 검색: {kw}")
            studies = self.search_studies(kw)
            for s in studies:
                acc = s.get("accession", "")
                if acc:
                    all_studies.append({
                        "source":    "MetaboLights",
                        "accession": acc,
                        "title":     s.get("title", ""),
                        "organism":  s.get("organism", ""),
                        "status":    s.get("studyStatus", ""),
                        "n_assays":  len(s.get("assays", []))
                    })

        df = pd.DataFrame(all_studies).drop_duplicates("accession")
        # 알려진 치주염 대사체 연구 수동 추가
        known = pd.DataFrame([
            {"source":"MetaboLights","accession":"MTBLS8357",
             "title":"Gingival tissue metabolomics - severe periodontitis (Chu 2024)",
             "organism":"Homo sapiens","status":"PUBLIC","n_assays":2}
        ])
        df = pd.concat([df, known]).drop_duplicates("accession")
        df.to_csv(META_DIR / "metabolights_periodontitis_studies.csv", index=False)
        print(f"\n✅ MetaboLights 연구: {len(df)}개")
        return df


# ═══════════════════════════════════════════════════════════════
# 4. STRING Collector — PPI Network
# ═══════════════════════════════════════════════════════════════
class STRINGCollector:
    """
    STRING v12.0 — 단백질 상호작용 네트워크 수집
    """
    BASE = "https://string-db.org/api"

    def get_functional_enrichment(self, gene_list: list,
                                   species: int = 9606) -> pd.DataFrame:
        """유전자 목록의 기능적 농축 분석"""
        params = {
            "identifiers": "%0d".join(gene_list),
            "species":     species,
            "caller_identity": "multiomics_periodontitis"
        }
        r = requests.post(
            f"{self.BASE}/tsv/enrichment",
            data=params, timeout=60
        )
        lines = r.text.strip().split('\n')
        if len(lines) < 2:
            return pd.DataFrame()

        rows = [l.split('\t') for l in lines]
        return pd.DataFrame(rows[1:], columns=rows[0])

    def get_network_image(self, gene_list: list,
                           outfile: Path) -> bool:
        """네트워크 이미지 다운로드"""
        params = {
            "identifiers": "%0d".join(gene_list),
            "species":     9606,
            "network_type": "functional",
            "caller_identity": "multiomics_periodontitis"
        }
        r = requests.post(
            f"{self.BASE}/image/network",
            data=params, timeout=60
        )
        if r.status_code == 200:
            with open(outfile, 'wb') as f:
                f.write(r.content)
            return True
        return False

    def download_species_network(self, species_id: int = 9606,
                                  min_score: int = 700) -> Optional[Path]:
        """
        전체 인간 PPI 네트워크 다운로드 (고신뢰도)
        - 파일: 9606.protein.links.v12.0.txt.gz (~150MB)
        - high confidence만: score >= 700
        """
        outdir = RAW_DIR / "string"
        outdir.mkdir(parents=True, exist_ok=True)
        url = (f"https://stringdb-downloads.org/download/"
               f"protein.links.v12.0/{species_id}.protein.links.v12.0.txt.gz")
        outfile = outdir / f"{species_id}.protein.links.v12.0.txt.gz"

        if outfile.exists():
            print(f"  STRING 파일 이미 존재: {outfile}")
            return outfile

        print(f"  STRING 전체 네트워크 다운로드 중 (~150MB)...")
        try:
            r = requests.get(url, stream=True, timeout=600)
            total = int(r.headers.get('content-length', 0))
            downloaded = 0
            with open(outfile, 'wb') as f:
                for chunk in r.iter_content(1024*1024):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        print(f"\r  {pct:.1f}%", end='', flush=True)
            print(f"\n  ✅ 완료: {outfile}")
            return outfile
        except Exception as e:
            print(f"  ❌ 오류: {e}")
            return None


# ═══════════════════════════════════════════════════════════════
# 5. HMDB Collector — 대사체 어노테이션
# ═══════════════════════════════════════════════════════════════
class HMDBCollector:
    """
    HMDB (Human Metabolome Database) — 대사체 정보
    """
    BASE = "https://hmdb.ca"

    def search_metabolite(self, name: str) -> list:
        """대사체 이름으로 HMDB ID + 구조 검색"""
        url = f"{self.BASE}/metabolites/search"
        params = {"query": name}
        try:
            r = requests.get(url, params=params, timeout=15)
            # HTML 파싱 필요 (API가 XML 반환)
            return []  # 실제 구현에서 BeautifulSoup 사용
        except Exception:
            return []

    def get_metabolite_info(self, hmdb_id: str) -> dict:
        """HMDB ID로 대사체 상세 정보"""
        url = f"{self.BASE}/metabolites/{hmdb_id}.xml"
        try:
            r = requests.get(url, timeout=15)
            root = ET.fromstring(r.content)

            def get_text(tag):
                el = root.find(f".//{tag}")
                return el.text if el is not None else ""

            return {
                "hmdb_id":   hmdb_id,
                "name":      get_text("name"),
                "formula":   get_text("chemical_formula"),
                "inchikey":  get_text("inchikey"),
                "kegg_id":   get_text("kegg_id"),
                "pathways":  [p.findtext("name","") for p in root.findall(".//pathway")[:10]],
                "diseases":  [d.findtext("name","") for d in root.findall(".//disease")[:10]],
                "biofluid":  [b.text for b in root.findall(".//biofluid_locations//biofluid")]
            }
        except Exception:
            return {"hmdb_id": hmdb_id}

    def build_metabolite_annotation_db(self,
                                        metabolite_names: list) -> pd.DataFrame:
        """대사체 목록 → HMDB 어노테이션 DB 구축"""
        records = []
        for name in metabolite_names:
            print(f"  HMDB 검색: {name}")
            # 이름 → HMDB ID 변환 (간단 매핑)
            info = {"name": name, "hmdb_id": "unknown",
                    "formula": "", "kegg_id": ""}
            records.append(info)
            time.sleep(0.2)

        df = pd.DataFrame(records)
        df.to_csv(META_DIR / "hmdb_metabolite_annotations.csv", index=False)
        return df


# ═══════════════════════════════════════════════════════════════
# 6. GTEx Collector — 정상 조직 발현 참고값
# ═══════════════════════════════════════════════════════════════
class GTExCollector:
    """
    GTEx Portal — 조직별 정상 유전자 발현
    (치주염 마커의 정상 기준값 확보용)
    """
    BASE = "https://gtexportal.org/api/v2"

    def get_gene_expression(self, gene_id: str,
                             tissue: str = "Minor_Salivary_Gland") -> dict:
        """특정 조직의 유전자 발현값"""
        url = f"{self.BASE}/expression/geneExpression"
        params = {
            "gencodeId": gene_id,
            "tissueSiteDetailId": tissue
        }
        try:
            r = requests.get(url, params=params, timeout=15)
            return r.json()
        except Exception:
            return {}

    def get_available_tissues(self) -> list:
        """사용 가능한 조직 목록"""
        url = f"{self.BASE}/dataset/tissueSiteDetail"
        try:
            r = requests.get(url, timeout=15)
            return r.json().get("tissueSiteDetail", [])
        except Exception:
            return []

    def get_oral_tissue_baselines(self, gene_list: list) -> pd.DataFrame:
        """
        구강 관련 조직의 정상 발현값 수집
        tissues: Minor_Salivary_Gland, Skin 등
        """
        oral_tissues = [
            "Minor_Salivary_Gland",
            "Skin_Sun_Exposed_Lower_leg",
            "Esophagus_Mucosa"
        ]
        records = []
        for gene in gene_list[:50]:  # 상위 50개만
            for tissue in oral_tissues:
                expr = self.get_gene_expression(gene, tissue)
                if expr:
                    records.append({
                        "gene": gene,
                        "tissue": tissue,
                        "median_tpm": expr.get("median", np.nan),
                        "n_samples": expr.get("numSamples", 0)
                    })
                time.sleep(0.1)

        df = pd.DataFrame(records)
        df.to_csv(META_DIR / "gtex_oral_baselines.csv", index=False)
        return df


# ═══════════════════════════════════════════════════════════════
# MAIN: 전체 수집 실행
# ═══════════════════════════════════════════════════════════════
def run_all_collectors():
    print("=" * 60)
    print("Global Omics DB Data Collection")
    print("=" * 60)

    summary = {}

    # GEO
    print("\n[1/5] GEO 수집...")
    geo = GEOCollector()
    geo_df = geo.collect_periodontitis_datasets()
    summary["geo"] = len(geo_df)

    # PRIDE
    print("\n[2/5] PRIDE 수집...")
    pride = PRIDECollector()
    pride_df = pride.collect_periodontitis_proteomics()
    summary["pride"] = len(pride_df)

    # MetaboLights
    print("\n[3/5] MetaboLights 수집...")
    ml_coll = MetaboLightsCollector()
    ml_df = ml_coll.collect_periodontitis_metabolomics()
    summary["metabolights"] = len(ml_df)

    # STRING (전체 네트워크는 선택적)
    print("\n[4/5] STRING 메타데이터...")
    # string = STRINGCollector()
    # string.download_species_network()  # 대용량 → 별도 실행
    summary["string"] = "manual (see 02_db_schema.py)"

    # 결과 요약
    print("\n" + "=" * 60)
    print("✅ 수집 완료 요약:")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print("=" * 60)

    # 통합 메타데이터 저장
    all_meta = pd.concat([
        geo_df.assign(omics="mRNA"),
        pride_df.assign(omics="Proteomics"),
        ml_df.assign(omics="Metabolomics")
    ], ignore_index=True)
    all_meta.to_csv(META_DIR / "all_public_datasets.csv", index=False)
    print(f"\n통합 메타데이터: {len(all_meta)}개 데이터셋")
    return all_meta


if __name__ == "__main__":
    run_all_collectors()
