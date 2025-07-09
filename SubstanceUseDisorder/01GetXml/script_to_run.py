
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
import json

PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")

# Ensure output directories exist
os.makedirs('outputs', exist_ok=True)
os.makedirs('xmls', exist_ok=True)

# Current UTC timestamp (timezone-aware)
now_utc = datetime.now(timezone.utc)

# base constants
BASE_NAME        = "pubmedsearch"
tz               = pytz.timezone("Asia/Kolkata")
today_suffix     = datetime.now(tz).strftime("%Y%m%d")

TABLE_NAME = f"{BASE_NAME}_{today_suffix}"
PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")
RETURN_MAX  = int(os.getenv("RETURN_MAX", "150"))
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
                         
                                 
# step 1: return pmids by boolean search
# Query PubMed & report unique PMIDs
# script to get pmids - save them (array - sorted)  # commit
# ─────────────────────────────────────────────────────────────────────────────
# Paths 
## keep relative path within substance use disorder (i.e. keep results within project - project is substance use disorder.)
pmid_csv     = "outputs/pmids.csv"

pmids = get_pmid_from_pubmed(
    query=pubmed_query,
    retmax=RETURN_MAX,
    api_key=PM_KEY,
)

n_unique = len(set(pmids))
print(f"ESearch returned {n_unique} unique PMIDs")
table_full_name = f"{CATALOG}.{SCHEMA}.{TABLE_NAME}_01_pmid_pmcid"

###!!! directly save pmids -- save file

def remove_if_exists(path):
    if os.path.isfile(path):
        os.remove(path)
        print(f"⚠️  Removed stale cache: {path}")


pmid_list = sorted(set(pmids))
os.makedirs("outputs", exist_ok=True)

###!!! save pmids as csv . save pmid_csv. overwrite prior files.

# ─────────────────────────────────────────────────────────────────────────────
# 1  fetch pmid metadata
# ─────────────────────────────────────────────────────────────────────────────

mapping_csv  = f"outputs/{TABLE_NAME}_01_pmid_pmcid.csv"
dist_csv     = f"outputs/{TABLE_NAME}_01_pmcid_distribution.csv"
### !!! read pmid (read from csv) - get pmid metadata - save them (csv, sorted by pmid)  # commit 

## write the csv table with metadata.

if os.path.isfile(pmid_json):
    with open(pmid_json) as f:
        loaded_pmids = json.load(f)
    print(f"✔ Loaded {len(loaded_pmids)} PMIDs from {pmid_json}")

    # check true hit count
    handle = Entrez.esearch(db="pubmed", term=pubmed_query, retmax=0, api_key=PM_KEY)
    true_count = int(Entrez.read(handle)["Count"])

    if len(loaded_pmids) == true_count:
        pmids = loaded_pmids
        print("✔ PMID list is up to date; skipping ESearch")
    else:
        print(f"ℹ️ PMID list ({len(loaded_pmids)}) is stale vs true count ({true_count}); re-running ESearch")
        # remove old downstream caches
        remove_if_exists(mapping_csv)
        remove_if_exists(dist_csv)

        pmids = get_pmid_from_pubmed(query=pubmed_query, retmax=RETURN_MAX, api_key=PM_KEY)
        pmids = sorted(set(pmids))
        # overwrite PMID cache
        with open(pmid_json, "w") as f:
            json.dump(pmids, f, indent=2)
        pd.Series(pmids, name="pmid").to_csv(pmid_csv, index=False)
        print(f"✔ Saved updated PMID list ({len(pmids)}) to {pmid_json} & {pmid_csv}")
else:
    # first run ever
    pmids = get_pmid_from_pubmed(query=pubmed_query, retmax=RETURN_MAX, api_key=PM_KEY)
    pmids = sorted(set(pmids))
    with open(pmid_json, "w") as f:
        json.dump(pmids, f, indent=2)
    pd.Series(pmids, name="pmid").to_csv(pmid_csv, index=False)
    print(f"✔ Saved PMID list ({len(pmids)}) to {pmid_json} & {pmid_csv}")

    
# ─────────────────────────────────────────────────────────────────────────────
# 1  Map PMIDs → PMCIDs (pandas DataFrame)
# ─────────────────────────────────────────────────────────────────────────────
# (re-use the same pmids list & PM_KEY from above)


## read from csv in step 1 - and get crosswak csv read pmid - get pmcid  - save them (csv - pmid, pmcid, sorted by pmid pmcid)  # commit
mapping_csv = f"outputs/{TABLE_NAME}_01_pmid_pmcid.csv"
if os.path.isfile(mapping_csv):
    # Load existing mapping
    mapping_df = pd.read_csv(mapping_csv)
    print(f"✔ Loaded existing mapping CSV; skipping PMID→PMCID download: {mapping_csv}")
else:
    # Do the ELink batch download
    mapping_df = map_pmids_to_pmcids(pmids, api_key=PM_KEY)
    if isinstance(mapping_df, list):
        mapping_df = pd.DataFrame(mapping_df)
    mapping_df["retrieved"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    mapping_df.to_csv(mapping_csv, index=False)
    print(f"✔ Wrote PMID→PMCID mapping to {mapping_csv}")

# ─────────────────────────────────────────────────────────────────────────────
# 2  Basic cardinalities (unique counts)
# ─────────────────────────────────────────────────────────────────────────────
unique_pmid         = mapping_df["pmid"].nunique()
unique_pmcid        = mapping_df["pmcid"].dropna().nunique()
pmid_with_pmcid     = mapping_df.dropna(subset=["pmcid"])["pmid"].nunique()
pmid_without_pmcid  = unique_pmid - pmid_with_pmcid

summary = {
    "unique_pmid":        unique_pmid,
    "unique_pmcid":       unique_pmcid,
    "pmid_with_pmcid":    pmid_with_pmcid,
    "pmid_without_pmcid": pmid_without_pmcid,
}
summary_df = pd.DataFrame([summary])
print("\n▶ PMCID Mapping Summary:")
print(summary_df.to_string(index=False))

# ─────────────────────────────────────────────────────────────────────────────
# 3  Distribution of #PMCIDs per PMID (only where ≥1 PMCID)
# ─────────────────────────────────────────────────────────────────────────────
dist_csv = f"outputs/{TABLE_NAME}_01_pmcid_distribution.csv"
if os.path.isfile(dist_csv):
    dist_df = pd.read_csv(dist_csv)
    print(f"✔ Loaded existing distribution CSV; skipping aggregation: {dist_csv}")
else:
    dist_df = (
        mapping_df
        .dropna(subset=["pmcid"])
        .groupby("pmid")["pmcid"]
        .nunique()
        .reset_index(name="pmcid_count")
        .sort_values("pmcid_count")
    )
    dist_df.to_csv(dist_csv, index=False)
    print(f"✔ Wrote PMCIDs‐per‐PMID distribution to {dist_csv}")

# (optional) save the distribution
csv_dist = f"outputs/{TABLE_NAME}_01_pmcid_distribution.csv"
if not os.path.isfile(csv_dist):
    dist_df.to_csv(csv_dist, index=False)
    print(f"Wrote distribution CSV: {csv_dist}")
else:
    print(f"Distribution CSV already exists—skipping: {csv_dist}")



# ─────────────────────────────────────────────────────────────────────────────
# 4  Load PMIDs → Fetch PubMed metadata → Save raw metadata & run QC
# ─────────────────────────────────────────────────────────────────────────────
import json

# Paths and directories
pmid_json = "outputs/pmids.json"
meta_dir  = "outputs/meta"
qc_dir    = os.path.join(meta_dir, "qc")
os.makedirs(qc_dir, exist_ok=True)

# 4.1 Load PMID list
with open(pmid_json) as f:
    pmids = json.load(f)
print(f"▶ Loaded {len(pmids)} PMIDs for metadata fetch")

# 4.2 Fetch PubMed metadata
meta_df = get_pubmed_metadata_pmid(pmids=pmids, api_key=PM_KEY)

# 4.3 Save raw metadata
meta_csv = os.path.join(meta_dir, f"pmid_metadata_{TABLE_NAME}.csv")
meta_df.to_csv(meta_csv, index=False)
print(f"✔ Wrote PubMed metadata ({len(meta_df)} rows) to {meta_csv}")

# 4.4 QC Step 1: Headline missing‐field counts & duplicates
cols = ["title","abstract","journal","publicationDate",
        "doi","firstAuthor","lastAuthor","meshTags","keywords"]
missing_counts = {f"missing_{c}": int(meta_df[c].isna().sum()) for c in cols}
total_rows     = len(meta_df)
unique_pmids   = meta_df["pmid"].nunique()
duplicate_pmids = total_rows - unique_pmids

headline_qc = {
    **missing_counts,
    "rows_in_table":    total_rows,
    "unique_pmids":     unique_pmids,
    "duplicate_pmids":  duplicate_pmids,
}
pd.DataFrame([headline_qc])\
  .to_csv(os.path.join(qc_dir, "qc_headline_counts.csv"), index=False)
print("▶ Headline QC saved to qc_headline_counts.csv")

# 4.5 QC Step 2: Publication‐year distribution
years = (
    meta_df["publicationDate"]
    .str.extract(r"(\d{4})", expand=False)
    .dropna()
    .astype(int)
)
year_counts = (
    years.value_counts()
         .sort_index()
         .rename_axis("year")
         .reset_index(name="count")
)
year_counts.to_csv(os.path.join(qc_dir, "qc_year_distribution.csv"), index=False)
print("▶ Year distribution saved to qc_year_distribution.csv")

# 4.6 QC Step 3: Top‐20 journals by article count
top_journals = (
    meta_df["journal"]
    .dropna()
    .value_counts()
    .head(20)
    .rename_axis("journal")
    .reset_index(name="count")
)
top_journals.to_csv(os.path.join(qc_dir, "qc_top_journals.csv"), index=False)
print("▶ Top journals saved to qc_top_journals.csv")

# 4.7 QC Step 4: MeSH‐tag & keyword frequencies (top-30 each)
def explode_and_count(df, col, sep=r";\s*", top_n=30):
    return (
        df[col]
        .dropna()
        .str.split(sep)
        .explode()
        .str.strip()
        .value_counts()
        .head(top_n)
        .rename_axis(col)
        .reset_index(name="count")
    )

mesh_counts = explode_and_count(meta_df, "meshTags", top_n=30)
kw_counts   = explode_and_count(meta_df, "keywords", top_n=30)

mesh_counts.to_csv(os.path.join(qc_dir, "qc_top_mesh_tags.csv"), index=False)
kw_counts.to_csv(os.path.join(qc_dir, "qc_top_keywords.csv"), index=False)
print("▶ MeSH & keyword counts saved to qc_top_mesh_tags.csv & qc_top_keywords.csv")


# ─────────────────────────────────────────────────────────────────────────────
# 4.8 Generate Publication‐Year Bar Chart (PNG)
# ─────────────────────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Use the year_counts DataFrame from Step 4.5
years = year_counts["year"].tolist()
counts = year_counts["count"].tolist()

plt.figure(figsize=(10, 5))
plt.bar(years, counts)  
plt.title("Articles per Publication Year")
plt.xlabel("Year")
plt.ylabel("Number of Articles")

# ─── Force integer ticks ─────────────────────────────────────────────────────
ax = plt.gca()
ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
plt.xticks(rotation=90)

plt.tight_layout()

# Save as before
chart_path = os.path.join(qc_dir, "qc_year_distribution.png")
plt.savefig(chart_path, bbox_inches="tight")
print(f"▶ Saved year distribution chart to {chart_path}")


# (Optionally) plt.show() if you want it displayed in an interactive session



# ─────────────────────────────────────────────────────────────────────────────
# 5  Load PMCID list → Fetch PMC metadata → Save raw CSV & run QC
# ─────────────────────────────────────────────────────────────────────────────
import pandas as pd
import json
import time

# Paths & dirs
map_csv   = f"outputs/{TABLE_NAME}_01_pmid_pmcid.csv"
pmc_meta_dir = "outputs/meta/pmc"
qc_pmc_dir   = os.path.join(pmc_meta_dir, "qc")
os.makedirs(qc_pmc_dir, exist_ok=True)

# 5.1 Load mapping CSV and get unique PMCIDs
map_df = pd.read_csv(map_csv, dtype=str)
pmcids = sorted(
    map_df["pmcid"].dropna().unique().tolist()
)
print(f"▶ Found {len(pmcids)} unique PMCIDs for metadata fetch")

# 5.2 Optional: strip 'PMC' prefix if your helper expects digits only
# pmcids = [x.lstrip("PMC") for x in pmcids]

# 5.3 Fetch PMC metadata in chunks (retry on failure)
def chunked_get_pmc_metadata(ids, chunk_size=50, max_retries=3, backoff=5):
    all_dfs = []
    for i in range(0, len(ids), chunk_size):
        batch = ids[i:i+chunk_size]
        for attempt in range(1, max_retries+1):
            try:
                df = get_pubmed_metadata_pmcid(pmcids=batch, api_key=PM_KEY)
                if isinstance(df, list):
                    df = pd.DataFrame(df)
                all_dfs.append(df)
                break
            except Exception as e:
                if attempt < max_retries:
                    print(f"⚠️  Batch {i//chunk_size+1} failed (attempt {attempt}): {e!r}, retrying in {backoff}s")
                    time.sleep(backoff)
                else:
                    raise
    return pd.concat(all_dfs, ignore_index=True)

pmc_meta_df = chunked_get_pmc_metadata(pmcids, chunk_size=100)

# 5.4 Save raw PMC metadata (sorted by pmcid)
pmc_meta_df = pmc_meta_df.sort_values("pmcid")
raw_csv = os.path.join(pmc_meta_dir, f"pmc_metadata_{TABLE_NAME}.csv")
pmc_meta_df.to_csv(raw_csv, index=False)
print(f"✔ Wrote PMC metadata ({len(pmc_meta_df)} rows) to {raw_csv}")

# 5.5 QC Step 1: Dedupe counts
total_rows    = len(pmc_meta_df)
unique_pmcids = pmc_meta_df["pmcid"].nunique()
unique_pmids  = pmc_meta_df["pmid"].nunique()
qc1 = {
    "rows_in_table": total_rows,
    "unique_pmcids": unique_pmcids,
    "unique_pmids":  unique_pmids
}
pd.DataFrame([qc1])\
  .to_csv(os.path.join(qc_pmc_dir, "qc_pmc_headline_counts.csv"), index=False)
print("▶ PMC headline QC saved to qc_pmc_headline_counts.csv")

# 5.6 QC Step 2: Distribution of PMCIDs per PMID & PMIDs per PMCID
dist_pmcid_per_pmid = (
    pmc_meta_df.groupby("pmid")["pmcid"]
               .nunique()
               .rename("pmcid_count")
               .reset_index()
               .pmcid_count.value_counts()
               .sort_index()
               .rename_axis("pmcid_per_pmid")
               .reset_index(name="freq")
)
dist_pmid_per_pmcid = (
    pmc_meta_df.groupby("pmcid")["pmid"]
               .nunique()
               .rename("pmid_count")
               .reset_index()
               .pmid_count.value_counts()
               .sort_index()
               .rename_axis("pmid_per_pmcid")
               .reset_index(name="freq")
)
dist_pmcid_per_pmid.to_csv(os.path.join(qc_pmc_dir, "qc_pmcid_per_pmid.csv"), index=False)
dist_pmid_per_pmcid.to_csv(os.path.join(qc_pmc_dir, "qc_pmid_per_pmcid.csv"), index=False)
print("▶ PMC/PMID distributions saved to qc_pmcid_per_pmid.csv & qc_pmid_per_pmcid.csv")

# 5.7 QC Step 3: Top‐20 journals
top_journals_pmc = (
    pmc_meta_df["journal"].dropna()
    .value_counts()
    .head(20)
    .rename_axis("journal")
    .reset_index(name="count")
)
top_journals_pmc.to_csv(os.path.join(qc_pmc_dir, "qc_top_journals_pmc.csv"), index=False)
print("▶ PMC top journals saved to qc_top_journals_pmc.csv")

# 5.8 QC Step 4: Year histogram
years_pmc = (
    pmc_meta_df["publicationDate"]
    .str.extract(r"(\d{4})", expand=False)
    .dropna().astype(int)
)
year_counts_pmc = (
    years_pmc.value_counts()
             .sort_index()
             .rename_axis("year")
             .reset_index(name="count")
)
year_counts_pmc.to_csv(os.path.join(qc_pmc_dir, "qc_year_histogram_pmc.csv"), index=False)
print("▶ PMC year histogram saved to qc_year_histogram_pmc.csv")

# 5.9 QC Step 5: MeSH tags & keywords
def explode_and_count(df, col, sep=r";\s*", top_n=30):
    return (
        df[col].dropna()
             .str.split(sep).explode().str.strip()
             .value_counts().head(top_n)
             .rename_axis(col)
             .reset_index(name="count")
    )

mesh_pmc = explode_and_count(pmc_meta_df, "meshTags", top_n=50)
kw_pmc   = explode_and_count(pmc_meta_df, "keywords", top_n=50)
mesh_pmc.to_csv(os.path.join(qc_pmc_dir, "qc_top_mesh_tags_pmc.csv"), index=False)
kw_pmc.to_csv(os.path.join(qc_pmc_dir, "qc_top_keywords_pmc.csv"), index=False)
print("▶ PMC MeSH & keyword counts saved to qc_top_mesh_tags_pmc.csv & qc_top_keywords_pmc.csv")


# ─────────────────────────────────────────────
# 5.10  Generate PMC Year‐of‐Publication Histogram (PNG)
# ─────────────────────────────────────────────
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Reuse the year_counts_pmc DataFrame from Step 5.8
years = year_counts_pmc["year"].tolist()
counts = year_counts_pmc["count"].tolist()

plt.figure(figsize=(12, 4))
plt.bar(years, counts)
plt.xlabel("Publication Year")
plt.ylabel("Number of Articles")
plt.title("PMC Articles by Publication Year")

# Force integer ticks on the x‐axis
ax = plt.gca()
ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
plt.xticks(rotation=90)
plt.tight_layout()

# Save the figure
hist_png = os.path.join(qc_pmc_dir, "qc_year_histogram_pmc.png")
plt.savefig(hist_png, bbox_inches="tight")
print(f"▶ Saved PMC year histogram chart to {hist_png}")




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
