# Citation Guide

How to cite MucOneUp in your research publications and acknowledge contributors.

---

## Software Citation

If you use MucOneUp in your research, please cite:

### BibTeX Format

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/berntpopp/MucOneUp},
  note = {Software version available at https://github.com/berntpopp/MucOneUp/releases}
}
```

### Text Format

> Popp, B. (2025). *MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit*. Available at: https://github.com/berntpopp/MucOneUp

### Version-Specific Citation

Include the specific version used in your research:

```bibtex
@software{muconeup2025,
  author = {Popp, Bernt},
  title = {MucOneUp: MUC1 VNTR Simulation and Analysis Toolkit},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/berntpopp/MucOneUp},
  version = {0.16.0},
  note = {Version 0.16.0, https://github.com/berntpopp/MucOneUp/releases/tag/v0.16.0}
}
```

**Find your version:**

```bash
muconeup --version
```

---

## Manuscript Citation

A manuscript describing MucOneUp is **in preparation**. This page will be updated when the manuscript is published.

**For now, cite using the software reference above.**

If you publish results using MucOneUp before the manuscript is available, please:

1. Use the software citation (BibTeX above)
2. Reference the GitHub repository URL
3. Specify the version number used
4. [Open an issue](https://github.com/berntpopp/MucOneUp/issues) to let us know about your publication

---

## Methods Section Template

Use this template when describing MucOneUp in your methods:

### Synthetic Data Generation

> Synthetic MUC1 VNTR sequences were generated using MucOneUp v0.16.0 (Popp, 2025; https://github.com/berntpopp/MucOneUp). Diploid haplotypes were simulated with [specify: fixed/random] VNTR repeat counts ([specify range/distribution]) using probability-based repeat transitions defined in the configuration file (provided in supplementary materials). [If applicable: Mutations were applied at positions [specify] using the [mutation name] mutation definition.] [If applicable: SNPs were integrated with a density of [X] per kb.] Sequencing reads were simulated using [Illumina/Oxford Nanopore/PacBio] parameters with [X]× coverage. All simulations used seed [number] for reproducibility.

### Example (Filled)

> Synthetic MUC1 VNTR sequences were generated using MucOneUp v0.16.0 (Popp, 2025; https://github.com/berntpopp/MucOneUp). Diploid haplotypes were simulated with fixed VNTR repeat counts of 60 repeats per haplotype using probability-based repeat transitions defined in the configuration file (provided in supplementary materials). The dupC mutation was applied at haplotype 1, position 25, to generate a mutated reference for variant caller benchmarking. Illumina paired-end reads (150 bp) were simulated with 100× coverage using the w-Wessim2 pipeline integrated in MucOneUp. All simulations used seed 42 for reproducibility.

---

## Example VNTR Database Citation

If you use the included example dataset (`data/examples/vntr_database.tsv`), cite the appropriate source for your research data.

**Note:** Verify the original source of any VNTR database you use and cite accordingly. The example database in this repository is provided for demonstration purposes.

---

## Acknowledging Tool Dependencies

MucOneUp integrates several external tools. If you use read simulation features, acknowledge the relevant tools:

### Illumina Read Simulation (w-Wessim2)

> Illumina read simulation was performed using the w-Wessim2 pipeline integrated in MucOneUp, which utilizes ReSeq (Schmeing & Robinson, 2021) for error modeling and BWA (Li & Durbin, 2009) for read alignment.

**ReSeq Citation:**

```bibtex
@article{schmeing2021reseq,
  title={ReSeq simulates realistic Illumina high-throughput sequencing data},
  author={Schmeing, Stephan and Robinson, Mark D},
  journal={Genome Biology},
  volume={22},
  number={1},
  pages={1--15},
  year={2021},
  publisher={BioMed Central}
}
```

**BWA Citation:**

```bibtex
@article{li2009bwa,
  title={Fast and accurate short read alignment with Burrows-Wheeler transform},
  author={Li, Heng and Durbin, Richard},
  journal={Bioinformatics},
  volume={25},
  number={14},
  pages={1754--1760},
  year={2009},
  publisher={Oxford University Press}
}
```

---

### Oxford Nanopore Read Simulation (NanoSim)

> Oxford Nanopore read simulation was performed using NanoSim (Yang et al., 2017) integrated in MucOneUp with diploid split-simulation to ensure balanced allelic coverage.

**NanoSim Citation:**

```bibtex
@article{yang2017nanosim,
  title={NanoSim: nanopore sequence read simulator based on statistical characterization},
  author={Yang, Chen and Chu, Justin and Warren, Ren{\'e} L and Birol, Inan{\c{c}}},
  journal={GigaScience},
  volume={6},
  number={4},
  pages={1--6},
  year={2017},
  publisher={Oxford University Press}
}
```

---

### PacBio HiFi Read Simulation (pbsim3)

> PacBio HiFi read simulation was performed using pbsim3 (Ono et al., 2021) with CCS consensus calling integrated in MucOneUp.

**pbsim3 Citation:**

```bibtex
@article{ono2021pbsim3,
  title={PBSIM3: a simulator for all types of PacBio and ONT long reads},
  author={Ono, Yukiteru and Hamada, Michiaki and Asai, Kiyoshi},
  journal={NAR Genomics and Bioinformatics},
  volume={4},
  number={2},
  pages={lqac092},
  year={2021},
  publisher={Oxford University Press}
}
```

---

### ORF Prediction (orfipy)

> Open reading frames were predicted using orfipy (Singh & Wurtele, 2021) integrated in MucOneUp.

**orfipy Citation:**

```bibtex
@article{singh2021orfipy,
  title={Orfipy: a fast and flexible tool for extracting ORFs},
  author={Singh, Urminder and Wurtele, Eve Syrkin},
  journal={Bioinformatics},
  volume={37},
  number={21},
  pages={3970--3972},
  year={2021},
  publisher={Oxford University Press}
}
```

---

## Supplementary Materials

When publishing results using MucOneUp, consider providing:

### Configuration File

Include your `config.json` in supplementary materials for full reproducibility.

**Example caption:**

> **Supplementary File 1:** MucOneUp configuration file (config.json) used for VNTR simulation, including repeat definitions, probability transitions, and mutation specifications.

---

### Simulation Parameters

Document all parameters used:

```
Simulation Parameters:
- MucOneUp version: 0.16.0
- Seed: 42
- VNTR length: Fixed at 60 repeats
- Reference assembly: hg38
- Mutation: dupC at haplotype 1, position 25
- Read simulator: Illumina (w-Wessim2)
- Coverage: 100×
- Read length: 150 bp paired-end
- Fragment size: 350 ± 50 bp
```

---

### Ground Truth Data

If benchmarking variant callers, include ground truth mutation positions:

```
Ground Truth Variants:
- Haplotype: 1
- Position: 25 (VNTR repeat index)
- Mutation type: dupC insertion
- Affected sequence: GCCCCACCCCTCCTCCCGCCGCGCCG
```

**Extract from simulation_stats.json:**

```bash
jq '.mutations' sample.simulation_stats.json > ground_truth.json
```

---

## Community Contributions

If MucOneUp enabled your research, we'd love to hear about it!

**Please:**

1. [Open an issue](https://github.com/berntpopp/MucOneUp/issues/new) with:
   - Publication title and journal
   - DOI or preprint link
   - Brief description of how you used MucOneUp

2. We'll showcase your work in our documentation and README

**Benefits:**

- Increase visibility for your research
- Help other users discover relevant applications
- Contribute to the MucOneUp user community

---

## License

MucOneUp is licensed under the **MIT License**:

```
MIT License

Copyright (c) 2025 Bernt Popp

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

**In short:**

- ✅ Commercial use allowed
- ✅ Modification allowed
- ✅ Distribution allowed
- ✅ Private use allowed
- ⚠️ No warranty provided
- ⚠️ No liability accepted

---

## Contact

**Maintainer:** Bernt Popp

**GitHub:** [@berntpopp](https://github.com/berntpopp)

**Repository:** [berntpopp/MucOneUp](https://github.com/berntpopp/MucOneUp)

**Issues:** [GitHub Issue Tracker](https://github.com/berntpopp/MucOneUp/issues)

**Discussions:** [GitHub Discussions](https://github.com/berntpopp/MucOneUp/discussions)

---

## Citation Checklist

When citing MucOneUp, ensure you include:

- [ ] Software citation (BibTeX or text format)
- [ ] Specific version number used (e.g., v0.16.0)
- [ ] GitHub repository URL
- [ ] Seed value used (for reproducibility)
- [ ] Configuration file in supplementary materials
- [ ] Relevant tool citations (ReSeq, NanoSim, etc.)
- [ ] Ground truth data (if benchmarking)

---

## Frequently Asked Questions

### Should I cite the GitHub repository or wait for the manuscript?

**For now, cite the GitHub repository using the software citation format.** The manuscript is in preparation and this page will be updated when published.

### How do I cite a specific version?

Include the version number in your citation:

```bibtex
@software{muconeup2025,
  ...
  version = {0.16.0},
  note = {Version 0.16.0, https://github.com/berntpopp/MucOneUp/releases/tag/v0.16.0}
}
```

### Do I need to cite external tools (ReSeq, NanoSim, etc.)?

**Yes, if you use them.** If you simulate Illumina reads, cite ReSeq and BWA. If you simulate ONT reads, cite NanoSim. See [Acknowledging Tool Dependencies](#acknowledging-tool-dependencies).

### Can I use MucOneUp in commercial research?

**Yes.** MucOneUp is licensed under MIT, which allows commercial use. However, please cite appropriately in any publications.

### How do I acknowledge MucOneUp in presentations?

Include:

- MucOneUp logo (if available in repository)
- GitHub repository URL
- Version number used
- Brief description: "Synthetic MUC1 VNTR data generated with MucOneUp v0.16.0"

---

## Update History

| Date | Update |
|------|--------|
| 2025-10-20 | Initial citation guide created |
| TBD | Will be updated when manuscript is published |

---

**Last Updated:** 2025-10-20
