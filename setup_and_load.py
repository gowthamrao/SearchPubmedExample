import subprocess
import sys

PACKAGES = [
    "pandas",
    "biopython",
    "rapidfuzz",
    "searchpubmed",
]

def install_packages():
    for pkg in PACKAGES:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])


def load_dependencies():
    import pandas as pd  # noqa: F401
    from Bio import Entrez  # noqa: F401
    from searchpubmed.pubmed import (  # noqa: F401
        get_pmid_from_pubmed,
        map_pmids_to_pmcids,
        get_pmc_full_xml,
        get_pubmed_metadata_pmid,
    )


def main():
    install_packages()
    load_dependencies()
    print("Dependencies installed and loaded.")


if __name__ == "__main__":
    main()
