# MucOneUp

[![Test & Quality status](https://github.com/vntyper-org/MucOneUp/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/vntyper-org/MucOneUp/actions/workflows/test.yml)
[![Python 3.10 to 3.12](https://img.shields.io/badge/python-3.10%E2%80%933.12-blue)](https://github.com/vntyper-org/MucOneUp/blob/main/pyproject.toml)
[![Documentation](https://img.shields.io/badge/docs-vntyper--org.github.io-blue)](https://vntyper-org.github.io/MucOneUp/)
[![Container image on ghcr.io](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/users/berntpopp/packages/container/package/muconeup%2Fmuconeup)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://vntyper-org.github.io/MucOneUp/about/license/)
[![Zenodo DOI 10.5281/zenodo.19740405](https://zenodo.org/badge/DOI/10.5281/zenodo.19740405.svg)](https://doi.org/10.5281/zenodo.19740405)

MucOneUp simulates diploid MUC1 VNTR haplotypes with known mutations and generates
sequencing reads from them, so that tools which call variants in the MUC1 VNTR can be
tested against a known truth. It is research software, for research use only.

It does three things:

- **Haplotypes.** Builds two VNTR alleles from configurable repeat-length distributions and
  repeat-to-repeat transition probabilities, then inserts, deletes or replaces bases at
  chosen repeats (for example the `dupC` frameshift). Random SNPs and a fixed `--seed` are
  supported, and every run writes a JSON truth record.
- **Reads.** Simulates Illumina paired-end reads (ReSeq2), Oxford Nanopore reads (NanoSim),
  PacBio HiFi reads (pbsim3) and PacBio/ONT amplicon reads. Built-in read profiles for ONT
  R10, calibrated on public PRJEB92208 data, add strand-specific errors, homopolymer and
  whole-unit PCR stutter, PCR artefacts and a per-read truth manifest
  (`muconeup reads profiles`).
- **Analysis.** Predicts open reading frames and scores the frameshift protein produced by
  ADTKD-MUC1 variants (`analyze orfs`), derives repeat-transition probabilities from a VNTR
  database (`analyze vntr-stats`) and simulates the SNaPshot assay in silico
  (`analyze snapshot-validate`).

## Install

MucOneUp is not on PyPI. Install it from a clone; the clone also gives you the default
`config.json`.

```bash
git clone https://github.com/vntyper-org/MucOneUp.git
cd MucOneUp
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

Haplotype simulation and ORF analysis need only Python 3.10 or newer. Read simulation also
needs external tools (BWA, samtools, ReSeq2, NanoSim, pbsim3, minimap2); the
[installation guide](https://vntyper-org.github.io/MucOneUp/getting-started/installation/)
sets them up with conda.

The container image bundles MucOneUp with all read simulators:

```bash
docker run --rm -v "$(pwd)/output:/data" ghcr.io/berntpopp/muconeup/muconeup:latest \
  --config /app/config.json simulate --out-dir /data --out-base sample
```

## Quick start

From the repository root, simulate a diploid sample with 60 repeats per allele and a `dupC`
frameshift in repeat 25 of haplotype 1, then predict its ORFs:

```bash
muconeup --config config.json simulate --fixed-lengths 60 \
  --mutation-name dupC --mutation-targets 1,25 --seed 42 \
  --out-dir output --out-base sample
muconeup --config config.json analyze orfs output/sample.001.simulated.fa \
  --out-dir output --out-base sample
```

`output/` then holds `sample.001.simulated.fa` (both haplotypes),
`sample.001.simulation_stats.json` (the truth record), `sample.001.mutated_unit.fa`,
`sample.orfs.fa` and `sample.orf_stats.json`.

To simulate reads from the FASTA, for example 100x Illumina coverage once the external
tools are installed:

```bash
muconeup --config config.json reads illumina output/sample.001.simulated.fa --coverage 100
```

The [quick start tutorial](https://vntyper-org.github.io/MucOneUp/getting-started/quickstart/)
continues from here, and the [CLI reference](https://vntyper-org.github.io/MucOneUp/reference/cli/)
lists every command and option.

## Documentation

The documentation at <https://vntyper-org.github.io/MucOneUp/> covers core concepts, the
simulation and read-simulation guides, the configuration reference and the published test
datasets. The sibling project [VNtyper](https://github.com/vntyper-org/VNtyper) genotypes
the MUC1 VNTR from sequencing data.

## Development

```bash
make init       # install uv, dev dependencies and pre-commit hooks
make ci-check   # the lint, format and type checks that CI runs
make test       # the test suite with coverage
```

See [CONTRIBUTING.md](CONTRIBUTING.md) before opening a pull request. Report bugs and ask
questions in the [issue tracker](https://github.com/vntyper-org/MucOneUp/issues).

## Citing

MucOneUp is pre-release software under active development; a manuscript is in preparation.
Until then, cite the Zenodo archive (concept DOI
[10.5281/zenodo.19740405](https://doi.org/10.5281/zenodo.19740405)) or use GitHub's
"Cite this repository" button, which reads [CITATION.cff](CITATION.cff). The
[citation guide](https://vntyper-org.github.io/MucOneUp/about/citation/) has BibTeX for
specific versions and citations for the bundled simulators.

## License

MIT. See the [license page](https://vntyper-org.github.io/MucOneUp/about/license/).
Maintained by [Bernt Popp](https://github.com/berntpopp).
