
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





# set up

import os
import pandas as pd

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

PM_KEY      = os.getenv("NCBI_API_KEY", "7ace9dd51ab7d522ad634bee5a1f4c46d409")

# Ensure output directories exist
os.makedirs('outputs', exist_ok=True)
os.makedirs('xmls', exist_ok=True)


# step 1: return pmids by boolean search
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

# Configuration
year_start = 2016
year_end   = 2025

query = '((("Substance-Related Disorders"[Mesh] OR "Behavior, Addictive"[Mesh] OR "Substance Abuse, Intravenous"[Mesh] OR "Drug Users"[Mesh] OR "Addiction Medicine"[Mesh] OR "Substance Withdrawal Syndrome"[Mesh] OR "Substance Abuse, Oral"[Mesh] OR substance use disorder*[tiab] OR drug use disorder*[tiab] OR substance misuse[tiab] OR drug misuse[tiab] OR substance depend*[tiab] OR drug depend*[tiab] OR substance abuse[tiab] OR drug abuse[tiab] OR addictive behavio?r[tiab] OR addiction[tiab] OR addicted person*[tiab] OR "nonmedical use"[tiab] OR "non-medical use"[tiab] OR illicit drug use[tiab] OR recreational drug use[tiab] OR problematic use[tiab] OR polysubstance use[tiab] OR poly-substance use[tiab] OR chemical depend*[tiab] OR compulsive drug use[tiab] OR drug seeking[tiab] OR SUD[tiab] OR "Tobacco Use Disorder"[Mesh] OR "Smoking"[Mesh] OR "Vaping"[Mesh] OR "Electronic Nicotine Delivery Systems"[Mesh] OR cigarette smoking[tiab] OR cigarette depend*[tiab] OR tobacco use[tiab] OR tobacco misuse[tiab] OR tobacco depend*[tiab] OR nicotine depend*[tiab] OR nicotine addiction[tiab] OR smoker*[tiab] OR vaping[tiab] OR e-cigarette*[tiab] OR e cigarette*[tiab] OR ENDS[tiab] OR "heated tobacco product*"[tiab] OR HTP[tiab] OR electronic cigarette*[tiab] OR "Marijuana Abuse"[Mesh] OR "Marijuana Smoking"[Mesh] OR cannabis use disorder*[tiab] OR marijuana use disorder*[tiab] OR cannabis misuse[tiab] OR marijuana misuse[tiab] OR cannabis depend*[tiab] OR marijuana depend*[tiab] OR cannabis abuse[tiab] OR marijuana abuse[tiab] OR THC use[tiab] OR cannabinoid use[tiab] OR THC misuse[tiab] OR tetrahydrocannabinol misuse[tiab] OR weed use[tiab] OR weed misuse[tiab] OR cannabis smoking[tiab] OR CBD misuse[tiab] OR "Opioid-Related Disorders"[Mesh] OR opioid use disorder*[tiab] OR opiate use disorder*[tiab] OR opioid misuse[tiab] OR opiate misuse[tiab] OR opioid depend*[tiab] OR opiate depend*[tiab] OR opioid abuse[tiab] OR opiate abuse[tiab] OR heroin depend*[tiab] OR heroin misuse[tiab] OR fentanyl misuse[tiab] OR "nonmedical prescription opioid use"[tiab] OR OUD[tiab] OR MAT[tiab] OR opioid addiction[tiab] OR synthetic opioid*[tiab] OR "medication-assisted treatment"[tiab] OR "medication assisted treatment"[tiab] OR "Alcohol-Related Disorders"[Mesh] OR "Alcoholism"[Mesh] OR "Binge Drinking"[Mesh] OR alcohol use disorder*[tiab] OR AUD[tiab] OR alcohol misuse[tiab] OR hazardous drink*[tiab] OR problem drink*[tiab] OR alcohol depend*[tiab] OR alcohol addiction[tiab] OR heavy drink*[tiab] OR binge drink*[tiab] OR harmful alcohol use[tiab] OR excessive alcohol[tiab] OR "alcohol overuse"[tiab] OR alcohol consumption disorder*[tiab]) AND (((Electronic Health Records[MeSH] OR Medical Record Systems, Computerized[MeSH] OR "routinely collected health data"[MeSH] OR EHR[TIAB] OR EMR[TIAB] OR "electronic health record"[TIAB] OR "electronic medical record"[TIAB] OR Insurance Claim Review[MeSH] OR Insurance Claim Reporting[MeSH] OR "claims data"[TIAB] OR "administrative data"[TIAB] OR "insurance claims"[TIAB] OR Databases, Factual[MeSH] OR "Real-World Data"[TIAB] OR "Real-World Evidence"[TIAB] OR "real-world data"[TIAB] OR "real-world evidence"[TIAB] OR "SEER"[TIAB] OR "NHANES"[TIAB] OR "CPRD"[TIAB] OR "MarketScan"[TIAB] OR "Optum"[TIAB] OR "Truven"[TIAB] OR "IQVIA"[TIAB] OR "PharMetrics"[TIAB] OR "Symphony Health"[TIAB] OR "Premier Healthcare"[TIAB] OR "Medicare"[TIAB] OR "Medicaid"[TIAB] OR "All-Payer"[TIAB] OR "All Payer"[TIAB] OR "TriNetX"[TIAB] OR "Cerner"[TIAB] OR "Komodo"[TIAB] OR "Kaiser"[TIAB] OR "Explorys"[TIAB] OR "The Health Improvement Network"[TIAB] OR "Vizient"[TIAB] OR "HealthVerity"[TIAB] OR "Datavant"[TIAB] OR "Merative"[TIAB]) AND (Observational Study[PT] OR Observational Studies as Topic[MeSH] OR observational[TIAB] OR "observational study"[TIAB] OR observational stud*[TIAB] OR Retrospective Studies[MeSH] OR retrospective[TIAB] OR "retrospective study"[TIAB] OR Secondary Data Analysis[MeSH] OR "secondary analysis"[TIAB] OR "secondary data analysis"[TIAB] OR Health Services Research[MeSH] OR Outcome Assessment, Health Care[MeSH] OR Comparative Effectiveness Research[MeSH] OR Cohort Studies[MeSH] OR cohort[TIAB] OR "cohort study"[TIAB] OR cohort stud*[TIAB] OR Longitudinal Studies[MeSH] OR "longitudinal study"[TIAB])) english[lang] ("2010"[dp] : "3000"[dp]) NOT (Clinical Trials as Topic[MeSH] OR Controlled Clinical Trials as Topic[MeSH] OR Randomized Controlled Trial[PT] OR Clinical Trial[PT]))) AND (("2016"[PDAT] : "2025"[PDAT])))'

all_pmids = get_pmid_from_pubmed(query=query,
                                 retmax=99999999,
                                 api_key=PM_KEY)
                                 
                                 


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
