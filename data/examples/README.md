# Example Datasets

This directory contains example datasets for demonstrating MucOneUp analysis features.

## Contents

### vntr_database.tsv

**Description:** VNTR structures from published MUC1 research

**Source:** PMID 29520014 - MUC1 VNTR variation in families with gastric and breast cancer

**Format:** Tab-separated values (TSV)

**Columns:**
- `publication`: PubMed ID reference
- `id`: Family/individual identifier (e.g., F1_III-10)
- `allele`: Allele number (1 or 2)
- `vntr`: VNTR structure (dash-separated repeat units)

**Statistics:**
- 44 alleles from published family studies
- 36 unique VNTR structures (after deduplication)
- Repeat lengths: 42-85 units
- Mean repeat length: 63.3 units
- Median repeat length: 70 units

**Usage:**

```bash
# Analyze transition probabilities
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv --header --structure-column vntr

# Save results to file
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv --header -o probabilities.json

# Extract specific metrics with jq
muconeup --config config.json analyze vntr-stats \
  data/examples/vntr_database.tsv --header | jq '.mean_repeats'
```

**Citation:**

If you use this dataset in your research, please cite:

> PMID 29520014 - [Insert full citation when available]

## Adding New Example Datasets

Example datasets should be:

- **Small** (<100 KB) to avoid repository bloat
- **Well-documented** with clear descriptions and metadata
- **Publicly shareable** without privacy/licensing concerns
- **Properly cited** if derived from published research
- **Representative** of real-world use cases

Place new datasets in this directory with a descriptive filename and update this README with documentation.
