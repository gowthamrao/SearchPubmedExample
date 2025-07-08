import os
# ─── Test‐mode override ──────────────────────────────────────────────────────
# os.environ["RETURN_MAX"] = "15"   # ← cap at 5 records for all subsequent loader calls
# ──────────────────────────────────────────────────────────────────────────

import pandas as pd
from loader_refactored import (
    download_pubmed_files,
    run_pubmed_loader,
    download_pubmed_metadata,
    qc_pubmed_metadata,
    qc_pmc_reconciliation,
    qc_pmc_xml
)


# Configuration
year_start = 2016
year_end   = 2025

# Ensure output directories exist
os.makedirs('outputs', exist_ok=True)
os.makedirs('xmls', exist_ok=True)

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
