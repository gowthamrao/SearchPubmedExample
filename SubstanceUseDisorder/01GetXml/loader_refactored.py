# loader_refactored.py  (trimmed to the important bits)

import os, pandas as pd
from datetime import datetime, timezone
from typing import List, Dict, Any
from Bio import Entrez
from searchpubmed.pubmed import (
    get_pmid_from_pubmed,
    map_pmids_to_pmcids,
    get_pmc_full_xml,
    get_pubmed_metadata_pmid
)
from searchpubmed.query_builder import build_query, STRATEGY2_OPTS

import xml.etree.ElementTree as ET
from rapidfuzz import fuzz





# 🔹  NEW:  import the full term lists
from terms import ALL_TERMS_LIST

# ------------------------------------------------------------------
# Config
# ------------------------------------------------------------------
PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")
RETURN_MAX  = int(os.getenv("RETURN_MAX", "100"))
DEFAULT_START = 2025
DEFAULT_END   = datetime.now(timezone.utc).year

Entrez.email   = os.getenv("ENTREZ_EMAIL", "you@example.com")
Entrez.api_key = PM_KEY

# RWD / real-world-data block (unchanged)
RWD_TERMS = build_query(STRATEGY2_OPTS)

# ------------------------------------------------------------------
# Build the exact PubMed query you used in Databricks
# ------------------------------------------------------------------
def build_pubmed_query() -> str:
    query = '((("Substance-Related Disorders"[Mesh] OR "Behavior, Addictive"[Mesh] OR "Substance Abuse, Intravenous"[Mesh] OR "Drug Users"[Mesh] OR "Addiction Medicine"[Mesh] OR "Substance Withdrawal Syndrome"[Mesh] OR "Substance Abuse, Oral"[Mesh] OR substance use disorder*[tiab] OR drug use disorder*[tiab] OR substance misuse[tiab] OR drug misuse[tiab] OR substance depend*[tiab] OR drug depend*[tiab] OR substance abuse[tiab] OR drug abuse[tiab] OR addictive behavio?r[tiab] OR addiction[tiab] OR addicted person*[tiab] OR "nonmedical use"[tiab] OR "non-medical use"[tiab] OR illicit drug use[tiab] OR recreational drug use[tiab] OR problematic use[tiab] OR polysubstance use[tiab] OR poly-substance use[tiab] OR chemical depend*[tiab] OR compulsive drug use[tiab] OR drug seeking[tiab] OR SUD[tiab] OR "Tobacco Use Disorder"[Mesh] OR "Smoking"[Mesh] OR "Vaping"[Mesh] OR "Electronic Nicotine Delivery Systems"[Mesh] OR cigarette smoking[tiab] OR cigarette depend*[tiab] OR tobacco use[tiab] OR tobacco misuse[tiab] OR tobacco depend*[tiab] OR nicotine depend*[tiab] OR nicotine addiction[tiab] OR smoker*[tiab] OR vaping[tiab] OR e-cigarette*[tiab] OR e cigarette*[tiab] OR ENDS[tiab] OR "heated tobacco product*"[tiab] OR HTP[tiab] OR electronic cigarette*[tiab] OR "Marijuana Abuse"[Mesh] OR "Marijuana Smoking"[Mesh] OR cannabis use disorder*[tiab] OR marijuana use disorder*[tiab] OR cannabis misuse[tiab] OR marijuana misuse[tiab] OR cannabis depend*[tiab] OR marijuana depend*[tiab] OR cannabis abuse[tiab] OR marijuana abuse[tiab] OR THC use[tiab] OR cannabinoid use[tiab] OR THC misuse[tiab] OR tetrahydrocannabinol misuse[tiab] OR weed use[tiab] OR weed misuse[tiab] OR cannabis smoking[tiab] OR CBD misuse[tiab] OR "Opioid-Related Disorders"[Mesh] OR opioid use disorder*[tiab] OR opiate use disorder*[tiab] OR opioid misuse[tiab] OR opiate misuse[tiab] OR opioid depend*[tiab] OR opiate depend*[tiab] OR opioid abuse[tiab] OR opiate abuse[tiab] OR heroin depend*[tiab] OR heroin misuse[tiab] OR fentanyl misuse[tiab] OR "nonmedical prescription opioid use"[tiab] OR OUD[tiab] OR MAT[tiab] OR opioid addiction[tiab] OR synthetic opioid*[tiab] OR "medication-assisted treatment"[tiab] OR "medication assisted treatment"[tiab] OR "Alcohol-Related Disorders"[Mesh] OR "Alcoholism"[Mesh] OR "Binge Drinking"[Mesh] OR alcohol use disorder*[tiab] OR AUD[tiab] OR alcohol misuse[tiab] OR hazardous drink*[tiab] OR problem drink*[tiab] OR alcohol depend*[tiab] OR alcohol addiction[tiab] OR heavy drink*[tiab] OR binge drink*[tiab] OR harmful alcohol use[tiab] OR excessive alcohol[tiab] OR "alcohol overuse"[tiab] OR alcohol consumption disorder*[tiab]) AND (((Electronic Health Records[MeSH] OR Medical Record Systems, Computerized[MeSH] OR "routinely collected health data"[MeSH] OR EHR[TIAB] OR EMR[TIAB] OR "electronic health record"[TIAB] OR "electronic medical record"[TIAB] OR Insurance Claim Review[MeSH] OR Insurance Claim Reporting[MeSH] OR "claims data"[TIAB] OR "administrative data"[TIAB] OR "insurance claims"[TIAB] OR Databases, Factual[MeSH] OR "Real-World Data"[TIAB] OR "Real-World Evidence"[TIAB] OR "real-world data"[TIAB] OR "real-world evidence"[TIAB] OR "SEER"[TIAB] OR "NHANES"[TIAB] OR "CPRD"[TIAB] OR "MarketScan"[TIAB] OR "Optum"[TIAB] OR "Truven"[TIAB] OR "IQVIA"[TIAB] OR "PharMetrics"[TIAB] OR "Symphony Health"[TIAB] OR "Premier Healthcare"[TIAB] OR "Medicare"[TIAB] OR "Medicaid"[TIAB] OR "All-Payer"[TIAB] OR "All Payer"[TIAB] OR "TriNetX"[TIAB] OR "Cerner"[TIAB] OR "Komodo"[TIAB] OR "Kaiser"[TIAB] OR "Explorys"[TIAB] OR "The Health Improvement Network"[TIAB] OR "Vizient"[TIAB] OR "HealthVerity"[TIAB] OR "Datavant"[TIAB] OR "Merative"[TIAB]) AND (Observational Study[PT] OR Observational Studies as Topic[MeSH] OR observational[TIAB] OR "observational study"[TIAB] OR observational stud*[TIAB] OR Retrospective Studies[MeSH] OR retrospective[TIAB] OR "retrospective study"[TIAB] OR Secondary Data Analysis[MeSH] OR "secondary analysis"[TIAB] OR "secondary data analysis"[TIAB] OR Health Services Research[MeSH] OR Outcome Assessment, Health Care[MeSH] OR Comparative Effectiveness Research[MeSH] OR Cohort Studies[MeSH] OR cohort[TIAB] OR "cohort study"[TIAB] OR cohort stud*[TIAB] OR Longitudinal Studies[MeSH] OR "longitudinal study"[TIAB])) english[lang] ("2010"[dp] : "3000"[dp]) NOT (Clinical Trials as Topic[MeSH] OR Controlled Clinical Trials as Topic[MeSH] OR Randomized Controlled Trial[PT] OR Clinical Trial[PT]))) AND (("2016"[PDAT] : "2025"[PDAT])))'
    return f"query"


# ------------------------------------------------------------------
# Download XMLs (PMC if possible, PubMed fallback otherwise)
# ------------------------------------------------------------------
def download_pubmed_files(start_year: int,
                          end_year:   int,
                          out_dir:    str) -> List[str]:

    query  = build_pubmed_query()
    print(query)

    
##    pmids  = get_pmid_from_pubmed(query=query,
##                                  retmax=RETURN_MAX,
##                                  api_key=PM_KEY)
    all_pmids = get_pmid_from_pubmed(query=query,
                                     retmax=RETURN_MAX,
                                     api_key=PM_KEY)
    pmids = all_pmids[:RETURN_MAX]      # enforce at most RETURN_MAX IDs
    raw_mapping = map_pmids_to_pmcids(pmids, api_key=PM_KEY)

    # Ensure we always work with list[dict]
    if isinstance(raw_mapping, pd.DataFrame):
        mappings = raw_mapping.to_dict("records")
    else:
        mappings = list(raw_mapping)

    os.makedirs(out_dir, exist_ok=True)
    paths: List[str] = []

    for rec in mappings:
        pmid, pmcid = rec["pmid"], rec.get("pmcid")
        ident = pmcid or pmid

        filepath = os.path.join(out_dir, f"{ident}.xml")
        # Skip if already downloaded
        if os.path.exists(filepath) and os.path.getsize(filepath) > 100:
            paths.append(filepath)
            continue

        try:
            # 1) Fetch raw bytes
            if pmcid:
                handle = Entrez.efetch(db="pmc", id=pmcid, retmode="xml")
            else:
                handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml")
            raw_xml = handle.read()
            handle.close()

            # Decode if needed
            if isinstance(raw_xml, (bytes, bytearray)):
                xml = raw_xml.decode("utf-8", errors="ignore")
            else:
                xml = raw_xml

            # Skip if the response is empty
            if not xml.strip():
                print(f"[WARN] Empty XML for {ident}")
                continue

            #filepath = os.path.join(out_dir, f"{ident}.xml")
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(xml)
               
            paths.append(filepath)

        except Exception as e:
            print(f"[ERROR] Fetch failed for {ident}: {e}")
            continue

    return paths


# ------------------------------------------------------------------
# Loader summary identical to the Spark notebook but with pandas
# ------------------------------------------------------------------
def run_pubmed_loader(start_year: int,
                      end_year:   int) -> Dict[str, Any]:
    query = build_pubmed_query()
    pmids = get_pmid_from_pubmed(query=query,
                                 retmax=RETURN_MAX,
                                 api_key=PM_KEY)
    raw = map_pmids_to_pmcids(pmids, api_key=PM_KEY)
    mappings = raw.to_dict("records") if isinstance(raw, pd.DataFrame) else list(raw)

    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for rec in mappings:
        rec["retrieved"] = now

    df = pd.DataFrame(mappings)
    unique_pmid = int(df["pmid"].nunique())
    unique_pmcid = int(df["pmcid"].dropna().nunique())
    pmid_with_pmc = int(df[df["pmcid"].notna()]["pmid"].nunique())
    pmid_without_pmc = unique_pmid - pmid_with_pmc

    summary = {
        "unique_pmid":        unique_pmid,
        "unique_pmcid":       unique_pmcid,
        "pmid_with_pmcid":    pmid_with_pmc,
        "pmid_without_pmcid": pmid_without_pmc,
    }

    pmcid_dist = (
        df.dropna(subset=["pmcid"])  
          .groupby("pmid")["pmcid"]
          .nunique()
          .to_dict()
    )

    return {
        "query":               query,
        "summary":             summary,
        "pmcid_distribution":  pmcid_dist,
        "mappings":            mappings,
    }
# ------------------------------------------------------------------
# PubMed metadata download & QC
# ------------------------------------------------------------------
def download_pubmed_metadata(start_year: int,
                             end_year:   int,
                             out_dir:    str) -> pd.DataFrame:
    query = build_pubmed_query()
    pmids = get_pmid_from_pubmed(query=query,
                                 retmax=RETURN_MAX,
                                 api_key=PM_KEY)
    meta_df = get_pubmed_metadata_pmid(pmids=pmids, api_key=PM_KEY)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"pmid_metadata_{start_year}_{end_year}.csv")
    meta_df.to_csv(csv_path, index=False)
    print(f"📥 Wrote PubMed metadata to {csv_path}")
    return meta_df

def qc_pubmed_metadata(meta_df: pd.DataFrame, out_dir: str = "outputs"):
    """
    Generates QC summaries:
      1) Counts missing titles/abstracts (printed + CSV)
      2) Histogram of publication years (printed + CSV)
      3) Top 10 journals by count (printed + CSV)
      4) Top 10 MeSH terms & keywords (printed + CSV)
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) Missing fields
    missing = {
        "missing_title": int(meta_df["title"].isna().sum()),
        "missing_abstract": int(meta_df["abstract"].isna().sum()),
        "total_records": len(meta_df)
    }
    pd.DataFrame([missing]).to_csv(f"{out_dir}/qc_missing_counts.csv", index=False)
    print("\n📊 Missing‐Field Counts:")
    for k, v in missing.items():
        print(f"  • {k}: {v}")

    # 2) Year histogram
    years = pd.to_datetime(meta_df["publicationDate"], errors="coerce").dt.year
    year_counts = years.value_counts().sort_index()
    hist_df = year_counts.rename_axis("year").reset_index(name="count")
    hist_df.to_csv(f"{out_dir}/qc_year_histogram.csv", index=False)
    print("\n📈 Publication Year Distribution:")
    print(hist_df.to_string(index=False))

    # 3) Top journals
    top_journals = meta_df["journal"].value_counts().head(10)\
                          .rename_axis("journal").reset_index(name="count")
    top_journals.to_csv(f"{out_dir}/qc_top_journals.csv", index=False)
    print("\n🏛️ Top 10 Journals:")
    print(top_journals.to_string(index=False))

    # 4) MeSH & keyword frequencies
    def explode_and_report(col_name):
        exploded = (
            meta_df[col_name]
            .dropna()
            .explode()
            .value_counts()
            .head(10)
            .rename_axis(col_name)
            .reset_index(name="count")
        )
        exploded.to_csv(f"{out_dir}/qc_top_{col_name}.csv", index=False)
        print(f"\n🔑 Top 10 {col_name.replace('_', ' ').title()}:")
        print(exploded.to_string(index=False))

    if "mesh_terms" in meta_df.columns:
        explode_and_report("mesh_terms")
    if "keywords" in meta_df.columns:
        explode_and_report("keywords")

# ------------------------------------------------------------------
# PMC Full-text download & reconciliation stubs
# ------------------------------------------------------------------

def download_pmc_fulltext(pmcids: List[str], out_dir: str) -> List[str]:
    os.makedirs(out_dir, exist_ok=True)
    paths: List[str] = []
    for pmcid in pmcids:
        filename = f"{pmcid}.xml"
        filepath = os.path.join(out_dir, filename)
        if os.path.exists(filepath) and os.path.getsize(filepath) > 100:
            paths.append(filepath)
            continue
        try:
            xml = get_pmc_full_xml(pmc_id=pmcid, api_key=PM_KEY)
            if not xml.strip(): print(f"[WARN] Empty PMC XML for {pmcid}"); continue
            with open(filepath, "w", encoding="utf-8") as f: f.write(xml)
            paths.append(filepath)
        except Exception as e:
            print(f"[ERROR] PMC fetch failed for {pmcid}: {e}")
    return paths


def fetch_pmc_fulltexts(start_year: int, end_year: int, out_dir: str) -> List[str]:
    res = run_pubmed_loader(start_year, end_year)
    pmcids = [rec.get("pmcid") for rec in res["mappings"] if rec.get("pmcid")]
    return download_pmc_fulltext(pmcids, out_dir)


# ------------------------------------------------------------------
# PMC reconciliation
# ------------------------------------------------------------------
def qc_pmc_reconciliation(pubmed_df: pd.DataFrame,
                          xml_dir:      str,
                          out_dir:      str = "outputs"):
    """
    Reconcile PubMed vs PMC titles using fuzzy matching.
    Writes a CSV report and prints summary counts.
    """
    # Parse PMC XMLs into DataFrame
    records = []
    for fname in os.listdir(xml_dir):
        if not fname.endswith('.xml'): continue
        path = os.path.join(xml_dir, fname)
        tree = ET.parse(path)
        root = tree.getroot()
        # Extract IDs and title
        meta = {'pmcid': None, 'pmid': None, 'title_pmc': None}
        for el in root.findall('.//article-id'):
            if el.attrib.get('pub-id-type') == 'pmcid': meta['pmcid'] = el.text
            if el.attrib.get('pub-id-type') == 'pmid':  meta['pmid'] = el.text
        title_el = root.find('.//article-title')
        meta['title_pmc'] = title_el.text.strip() if title_el is not None else ''
        records.append(meta)
    pmc_df = pd.DataFrame(records)

    # Merge with pubmed_df on pmid
    merged = pubmed_df.merge(pmc_df, on='pmid', how='inner')
    # Compute similarity
    merged['title_pubmed'] = merged['title']
    merged['sim_score'] = merged.apply(
        lambda row: fuzz.token_sort_ratio(row['title_pubmed'], row['title_pmc']), axis=1)

    # Classify
    total = len(merged)
    exact = (merged['sim_score'] == 100).sum()
    near = ((merged['sim_score'] >= 90) & (merged['sim_score'] < 100)).sum()
    mismatches = (merged['sim_score'] < 90).sum()

    # Save full report
    os.makedirs(out_dir, exist_ok=True)
    report_path = os.path.join(out_dir, 'qc_pmc_reconciliation.csv')
    merged.to_csv(report_path, index=False)

    # Print summary
    print(f"\n🔍 PMC Reconciliation Report ({total} records):")
    print(f"  • Exact matches: {exact}")
    print(f"  • Near matches (>=90): {near}")
    print(f"  • Mismatches (<90): {mismatches}")
    print(f"Report saved to {report_path}")

# ------------------------------------------------------------------
# XML QA
# ------------------------------------------------------------------
def qc_pmc_xml(xml_dir: str, out_dir: str = "outputs"):
    """
    QA of PMC XML files: check for <body>, count <sec> and <p> tags.
    Outputs CSV summary and prints overview.
    """
    summary = []
    for fname in os.listdir(xml_dir):
        if not fname.endswith('.xml'): continue
        path = os.path.join(xml_dir, fname)
        tree = ET.parse(path)
        root = tree.getroot()
        body = root.find('.//body')
        secs = root.findall('.//sec')
        pars = root.findall('.//p')
        summary.append({
            'file': fname,
            'has_body': bool(body),
            'sec_count': len(secs),
            'p_count': len(pars)
        })
    qa_df = pd.DataFrame(summary)
    os.makedirs(out_dir, exist_ok=True)
    qa_path = os.path.join(out_dir, 'qc_pmc_xml_summary.csv')
    qa_df.to_csv(qa_path, index=False)

    # Print overview
    total = len(qa_df)
    with_body = qa_df['has_body'].sum()
    avg_secs = qa_df['sec_count'].mean()
    avg_ps = qa_df['p_count'].mean()
    print(f"\n📄 PMC XML QA ({total} files):")
    print(f"  • Files with <body>: {with_body}/{total}")
    print(f"  • Avg <sec> tags per file: {avg_secs:.1f}")
    print(f"  • Avg <p> tags per file: {avg_ps:.1f}")
    print(f"Summary saved to {qa_path}")
