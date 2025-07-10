
# Folder 1 - get xml
# set up - with configuration details 
# script to get pmids - save them (array - sorted)  # commit
# read pmid - get pmid metadata - save them (csv, sorted by pmid)  # commit
# read pmid - get pmcid  - save them (csv - pmid, pmcid, sorted by pmid pmcid)  # commit
# read pmcid and get metadata - save them (csv - sorted pmcid)  # commit
# match pmid - pmcid using metadata - so that there is only one pmcid for each pmid. left join. there maybe pmid without pmcid - save them (csv - pmid, pmcid, merged metadata no duplicated columns pmid, pmcid s0rted
# for each pmcid with pmid - get xml. save xml

# Folder 2 - chunk xml using xml tree - create one csv per pmcid with paragraph chunks. Greedy, do not loose text. If chunk size is less than 100 character, merge to next paragraph iteratively.

# Folder 3 - for each csv, loop thru each row, run pyregularexpression and create indicator 
## output is pmid, pmcid, paragraph chunk, indicator variables (see databricks)


############################################
### set up/config/connection

import os, pandas as pd
from datetime import datetime, timezone
from typing import List, Dict, Any
from Bio import Entrez
from searchpubmed.query_builder import build_query, STRATEGY2_OPTS
from searchpubmed.pubmed import (
    get_pmid_from_pubmed, 
    map_pmids_to_pmcids,
    get_pmc_full_xml,
    get_pubmed_metadata_pmid,
    get_pubmed_metadata_pmcid
)
import pytz
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")

# Create outputs/ and xmls/ under this script’s folder
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))


# Current UTC timestamp (timezone-aware)
now_utc = datetime.now(timezone.utc)

# base constants
BASE_NAME        = "pubmedsearch"
tz               = pytz.timezone("Asia/Kolkata")
today_suffix     = datetime.now(tz).strftime("%Y%m%d")

TABLE_NAME = f"{BASE_NAME}_{today_suffix}"
PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")
RETURN_MAX  = int(os.getenv("RETURN_MAX", "80"))
CATALOG      = "odysseus"        # catalog
SCHEMA       = "ods_pd_0160"     # schema / database

##Entrez.email   = os.getenv("ENTREZ_EMAIL", "you@example.com")
Entrez.api_key = PM_KEY


# -------------------------------------------------------------------
# 0 Create search criteria pubmed query constructor
# -------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Search-specific configuration
# ---------------------------------------------------------------------------
RWD_TERMS     = build_query(STRATEGY2_OPTS)
# -------------------------------------------------------------------
# Addiction, SUD – MeSH + free-text
# -------------------------------------------------------------------
ADDICTION_TERMS: list[str] = [
    # Core MeSH headings (explode by default) – maximise recall
    '"Substance-Related Disorders"[Mesh]',
    '"Behavior, Addictive"[Mesh]',
    '"Substance Abuse, Intravenous"[Mesh]',
    '"Drug Users"[Mesh]',
    '"Addiction Medicine"[Mesh]',
    # Extra MeSH to catch detox / route-specific misuse
    '"Substance Withdrawal Syndrome"[Mesh]',
    '"Substance Abuse, Oral"[Mesh]',
    # Key free-text variants (DSM-5 wording, UK/US spelling, wildcards)
    'substance use disorder*[tiab]',
    'drug use disorder*[tiab]',
    'substance misuse[tiab]',
    'drug misuse[tiab]',
    'substance depend*[tiab]',
    'drug depend*[tiab]',
    'substance abuse[tiab]',
    'drug abuse[tiab]',
    'addictive behavio?r[tiab]',      # handles UK “behaviour”
    'addiction[tiab]',
    'addicted person*[tiab]',
    '"nonmedical use"[tiab]',          # prescription diversion
    '"non-medical use"[tiab]',
    'illicit drug use[tiab]',
    'recreational drug use[tiab]',
    'problematic use[tiab]',
    # Emerging / specialist phrases
    'polysubstance use[tiab]',
    'poly-substance use[tiab]',
    'chemical depend*[tiab]',
    'compulsive drug use[tiab]',
    'drug seeking[tiab]',
    'SUD[tiab]',                       # low-noise acronym in biomedical corpus
]

# -------------------------------------------------------------------
# Tobacco / Nicotine – add e-cigs & heated products
# -------------------------------------------------------------------
TOBACCO_TERMS: list[str] = [
    '"Tobacco Use Disorder"[Mesh]',
    '"Smoking"[Mesh]',
    '"Vaping"[Mesh]',                               # 2017 MeSH
    '"Electronic Nicotine Delivery Systems"[Mesh]', # “ENDS”, 2018 MeSH
    'cigarette smoking[tiab]',
    'cigarette depend*[tiab]',
    'tobacco use[tiab]',
    'tobacco misuse[tiab]',
    'tobacco depend*[tiab]',
    'nicotine depend*[tiab]',
    'nicotine addiction[tiab]',
    'smoker*[tiab]',
    'vaping[tiab]',
    'e-cigarette*[tiab] OR e cigarette*[tiab]',
    'ENDS[tiab]',
    '"heated tobacco product*"[tiab] OR HTP[tiab]',
    'electronic cigarette*[tiab]',
]

# -------------------------------------------------------------------
# Cannabis – DSM-5 term + street synonyms
# -------------------------------------------------------------------
CANNABIS_TERMS: list[str] = [
    '"Marijuana Abuse"[Mesh]',
    '"Marijuana Smoking"[Mesh]',
    'cannabis use disorder*[tiab]',
    'marijuana use disorder*[tiab]',
    'cannabis misuse[tiab]',
    'marijuana misuse[tiab]',
    'cannabis depend*[tiab]',
    'marijuana depend*[tiab]',
    'cannabis abuse[tiab]',
    'marijuana abuse[tiab]',
    'THC use[tiab]',
    'cannabinoid use[tiab]',
    # New / colloquial turns of phrase
    'THC misuse[tiab] OR tetrahydrocannabinol misuse[tiab]',
    'weed use[tiab] OR weed misuse[tiab]',
    'cannabis smoking[tiab]',
    'CBD misuse[tiab]',
]

# -------------------------------------------------------------------
# Opioid – prescription & illicit, treatment context
# -------------------------------------------------------------------
OPIOID_TERMS: list[str] = [
    '"Opioid-Related Disorders"[Mesh]',
    'opioid use disorder*[tiab]',
    'opiate use disorder*[tiab]',
    'opioid misuse[tiab]',
    'opiate misuse[tiab]',
    'opioid depend*[tiab]',
    'opiate depend*[tiab]',
    'opioid abuse[tiab]',
    'opiate abuse[tiab]',
    'heroin depend*[tiab]',
    'heroin misuse[tiab]',
    'fentanyl misuse[tiab]',
    '"nonmedical prescription opioid use"[tiab]',
    'OUD[tiab]',                          # clinically common abbreviation
    'MAT[tiab]',                          # medication-assisted treatment
    'opioid addiction[tiab]',
    'synthetic opioid*[tiab]',
    '"medication-assisted treatment"[tiab] OR "medication assisted treatment"[tiab]'
]

# -------------------------------------------------------------------
# Alcohol – hazardous/binge wording, ICD/DSM acronyms
# -------------------------------------------------------------------
ALCOHOL_TERMS: list[str] = [
    '"Alcohol-Related Disorders"[Mesh]',
    '"Alcoholism"[Mesh]',
    '"Binge Drinking"[Mesh]',
    'alcohol use disorder*[tiab]',
    'AUD[tiab]',                          # “alcohol-use disorder”
    'alcohol misuse[tiab]',
    'hazardous drink*[tiab]',
    'problem drink*[tiab]',
    'alcohol depend*[tiab]',
    'alcohol addiction[tiab]',
    'heavy drink*[tiab]',
    'binge drink*[tiab]',
    'harmful alcohol use[tiab]',
    'excessive alcohol[tiab]',
    '"alcohol overuse"[tiab]',
    'alcohol consumption disorder*[tiab]'
]

ALL_TERMS_LIST = (
    ADDICTION_TERMS
    + TOBACCO_TERMS
    + CANNABIS_TERMS
    + OPIOID_TERMS
    + ALCOHOL_TERMS
)

year_start = 2016
year_end   = 2025

# Construct query
base = " OR ".join(ALL_TERMS_LIST)
q = f"({base})"
q = f"({q} AND ({RWD_TERMS}))"
date_block = (
    f'("{2016}"[PDAT] : "{2025}"[PDAT])'
)
pubmed_query = f"({q} AND ({date_block}))"
                         
                                 
# ─────────────────────────────────────────────────────────────────────────────
# 0 Query PubMed & report unique PMIDs
# script to get pmids - save them (array - sorted)  # commit
# ─────────────────────────────────────────────────────────────────────────────

# Paths & dirs
# add outputs and qc directories in the future

## keep relative path within substance use disorder (i.e. keep results within project - project is substance use disorder.)
pmid_csv     = os.path.join(BASE_DIR, "pmids.csv")

pmids = get_pmid_from_pubmed(
    query=pubmed_query,
    retmax=RETURN_MAX,
    api_key=PM_KEY,
)

n_unique = len(set(pmids))
print(f"ESearch returned {n_unique} unique PMIDs")
table_full_name = f"{CATALOG}.{SCHEMA}.{TABLE_NAME}_01_pmid_pmcid"

pmid_list = sorted(set(pmids))
pd.Series(pmid_list, name="pmid").to_csv(pmid_csv, index=False)
print(f"✔ Saved {len(pmid_list)} PMIDs to {pmid_csv}")

# ─────────────────────────────────────────────────────────────────────────────
# 1  Read PMIDs from CSV → Fetch PMID metadata → Save sorted CSV
# ─────────────────────────────────────────────────────────────────────────────

##mapping_csv   = os.path.join(OUTPUT_DIR, f"{TABLE_NAME}_01_pmid_pmcid.csv")
##dist_csv      = os.path.join(OUTPUT_DIR, f"{TABLE_NAME}_01_pmcid_distribution.csv")

# 1.1 Load the saved PMIDs
pmid_df      = pd.read_csv(pmid_csv, dtype=str)
pmids        = pmid_df["pmid"].dropna().unique().tolist()
print(f"▶ Loaded {len(pmids)} PMIDs from {pmid_csv}")

# 1.2 Fetch metadata for each PMID
meta_df_pmid = get_pubmed_metadata_pmid(pmids=pmids, api_key=PM_KEY)

# 1.3 Sort by PMID and save
meta_df_pmid = meta_df_pmid.sort_values("pmid")
pmid_meta_csv = f"{TABLE_NAME}_pmid_METADATA.csv"
meta_df_pmid.to_csv(pmid_meta_csv, index=False)
print(f"✔ Wrote PubMed metadata ({len(meta_df_pmid)} rows) to {pmid_meta_csv}")


# ─────────────────────────────────────────────────────────────────────────────
# 2  Read PMIDs → Map to PMCIDs → Save mapping CSV (sorted)   # commit
# ─────────────────────────────────────────────────────────────────────────────
# (re-use the same pmids list & PM_KEY from above)

# 2.1 Reload the saved PMIDs
##print(f"▶ Loaded {len(pmids)} PMIDs from {pmid_csv}")

# 2.2 Map PMIDs → PMCIDs using your helper
mapping_records = map_pmids_to_pmcids(pmids, api_key=PM_KEY)
mapping_df      = pd.DataFrame(mapping_records)

# 2.3 Sort and save
mapping_df = mapping_df.sort_values(["pmid", "pmcid"])
mapping_csv = f"{TABLE_NAME}_01_pmid_pmcid.csv"
mapping_df.to_csv(mapping_csv, index=False)
print(f"✔ Wrote PMID→PMCID mapping ({len(mapping_df)} rows) to {mapping_csv}")
      

## read from csv in step 1 - and get crosswak csv read pmid - get pmcid  - save them (csv - pmid, pmcid, sorted by pmid pmcid)  # commit
## Completed - next stats

# ─────────────────────────────────────────────────────────────────────────────
# 3  Basic cardinalities (unique counts) and save QC CSV
# ─────────────────────────────────────────────────────────────────────────────
unique_pmid        = mapping_df["pmid"].nunique()
unique_pmcid       = mapping_df["pmcid"].dropna().nunique()
pmid_with_pmcid    = mapping_df.dropna(subset=["pmcid"])["pmid"].nunique()
pmid_without_pmcid = unique_pmid - pmid_with_pmcid

qc_summary = {
    "unique_pmid":        unique_pmid,
    "unique_pmcid":       unique_pmcid,
    "pmid_with_pmcid":    pmid_with_pmcid,
    "pmid_without_pmcid": pmid_without_pmcid,
}
summary_df = pd.DataFrame([qc_summary])

# Print to console
print("\n▶ PMID→PMCID Mapping Summary:")
print(summary_df.to_string(index=False))

# Save QC CSV
qc_csv = f"{TABLE_NAME}_qc_pmid_pmcid_SUMMARY.csv"
summary_df.to_csv(qc_csv, index=False)
print(f"✔ Wrote mapping QC summary to {qc_csv}")


# ─────────────────────────────────────────────────────────────────────────────
# 4  Load PubMed metadata → run QC
# ─────────────────────────────────────────────────────────────────────────────
# Paths - add if need separate dir

# 4.1 Load the PubMed‐level metadata you saved in Step 1
meta_df_pmid  = pd.read_csv(pmid_meta_csv, dtype=str)
print(f"▶ Loaded PubMed metadata ({len(meta_df_pmid)} rows) from {pmid_meta_csv}")
  
# 4.2 Headline QC: missing fields & duplicates
cols = ["title","abstract","journal","publicationDate","doi","firstAuthor","lastAuthor","meshTags","keywords"]
qc = {f"missing_{c}": int(meta_df_pmid[c].isna().sum()) for c in cols}
qc.update({
    "rows_in_table":   len(meta_df_pmid),
    "unique_pmids":    meta_df_pmid["pmid"].nunique(),
})
qc["duplicate_pmids"] = qc["rows_in_table"] - qc["unique_pmids"]
qc_df = pd.DataFrame([qc])
qc_filename = "qc_PMID_headline_counts.csv"
qc_df.to_csv(qc_filename, index=False)
print(f"✔ Saved QC data to {qc_filename}")

# 4.3 Year distribution

# 4.4 Top journals

# 4.5 MeSH & keyword frequencies

# 4.6 Publication‐year bar chart

# ─────────────────────────────────────────────────────────────────────────────
# 5  Distribution of #PMCIDs per PMID (only where ≥1 PMCID) and save CSV
# ─────────────────────────────────────────────────────────────────────────────
##dist_df = (
##    mapping_df
##    .dropna(subset=["pmcid"])
##    .groupby("pmid")["pmcid"]
##    .nunique()
##    .reset_index(name="pmcid_count")
##    .sort_values("pmcid_count")
##)
dist_df = (
    mapping_df
    .loc[mapping_df["pmcid"].notna()]
    .groupby("pmid", as_index=False)
    .agg(pmcid_count=("pmcid", "nunique"))
    .sort_values("pmcid_count")
)

dist_csv = f"{TABLE_NAME}_01_pmcid_distribution.csv"
dist_df.to_csv(dist_csv, index=False)
print(f"✔ Wrote PMCIDs‐per‐PMID distribution to {dist_csv}")

# ─────────────────────────────────────────────────────────────────────────────
# 6  Load PMCID list → Fetch PMC metadata → Save raw CSV & run QC
# ─────────────────────────────────────────────────────────────────────────────

# 6.1  Load mapping CSV and extract unique PMCIDs
map_df  = pd.read_csv(mapping_csv, dtype=str)
pmcids  = sorted(map_df["pmcid"].dropna().unique())
print(f"▶ Found {len(pmcids)} unique PMCIDs for metadata fetch")


# 6.2  Fetch PMC metadata in retry‐chunked batches
def chunked_get_pmc_metadata(ids, chunk_size=50, retries=3, backoff=5):
    frames = []
    for i in range(0, len(ids), chunk_size):
        batch = ids[i:i+chunk_size]
        for attempt in range(1, retries+1):
            try:
                df = get_pubmed_metadata_pmcid(pmcids=batch, api_key=PM_KEY)
                frames.append(pd.DataFrame(df))
                break
            except Exception as e:
                if attempt < retries:
                    print(f"⚠️  Batch {i//chunk_size+1} failed ({e}), retrying in {backoff}s…")
                    time.sleep(backoff)
                else:
                    raise
    return pd.concat(frames, ignore_index=True)

pmc_meta_df = chunked_get_pmc_metadata(pmcids, chunk_size=100)

# 6.4  Save raw PMC metadata (sorted by pmcid)
pmc_meta_df.sort_values("pmcid", inplace=True)
raw_csv = f"{TABLE_NAME}_pmc_METADATA.csv"
pmc_meta_df.to_csv(raw_csv, index=False)
print(f"✔ Wrote PMC metadata ({len(pmc_meta_df)} rows) to {raw_csv}")
    

# 6.5  QC Step 1: headline counts (rows, unique PMCIDs, unique PMIDs)
  
# 6.6  QC Step 2: distribution – PMCIDs per PMID & PMIDs per PMCID

# 6.7  QC Step 3: top 20 journals
    
# 6.9  QC Step 5: MeSH tags & keyword freqs (top 50)

# 6.10  Generate PMC publication‐year bar chart (PNG)



xxx


# 1) Download all article XMLs (PMC full-text or PubMed fallback)
print("📄 Downloading article XMLs (PMC + PubMed fallback)...")
xml_paths = download_pubmed_files(year_start, year_end, out_dir='xmls')
print(f"Downloaded {len(xml_paths)} XML files to 'xmls/'")

# 2) PubMed loader summary
print("🔍 Running PubMed loader summary...")
loader_result = run_pubmed_loader(year_start, year_end)
print(pd.DataFrame([loader_result['summary']]))

# 3) Download PubMed metadata & QC
print("📥 Downloading PubMed metadata...")
meta_df = download_pubmed_metadata(year_start, year_end, out_dir='outputs')
print(f"{len(meta_df)} PubMed records fetched")

print("📊 Running PubMed QC...")
qc_pubmed_metadata(meta_df, out_dir='outputs')

# 4) PMC reconciliation
if xml_paths:
    print("🔍 Running PMC reconciliation...")
    qc_pmc_reconciliation(meta_df, 'xmls', out_dir='outputs')
else:
    print("⚠️ No XML files found, skipping PMC reconciliation.")

# 5) PMC XML QA
if xml_paths:
    print("✅ Running PMC XML QA...")
    qc_pmc_xml('xmls', out_dir='outputs')
else:
    print("⚠️ No XML files found, skipping PMC XML QA.")

print("🎉 Pipeline complete! Check 'outputs/' and 'xmls/' for results.")
