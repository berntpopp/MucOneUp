# MucOneUp

**MucOneUp** is a Python tool for simulating **MUC1 VNTR diploid references**. It builds customized references that:

1. **Generate** haplotypes containing a variable‐length VNTR region with specified probabilities.  
2. **Force** a canonical terminal block (`6` or `6p` → `7 → 8 → 9`) before the right‐hand constant.  
3. Optionally **introduce mutations** (inserts, deletes, or replaces) in selected repeats.

## Table of Contents

- [Installation](#installation)  
- [Quick Start](#quick-start)  
- [Usage](#usage)  
  - [Command-Line Arguments](#command-line-arguments)  
  - [Example Commands](#example-commands)  
- [Project Structure and Logic](#project-structure-and-logic)  
  - [Modules](#modules)  
  - [Config File Layout](#config-file-layout)  
- [License](#license)

---

## Installation

1. **Clone or download** this repository.  
2. **Install** in a Python 3.7+ environment:
   ```bash
   pip install .
   ```
   This will install the `muc_one_up` Python package locally.

3. (Optional) If you want to make edits and test them, install in **editable** mode:
   ```bash
   pip install -e .
   ```

Once installed, you’ll have a command-line program called **`muconeup`** available.

---

## Quick Start

1. **Create** or update a **JSON config** (see [Config File Layout](#config-file-layout)) describing your repeats, probabilities, mutations, etc.  
2. **Run** the tool, specifying your config and any other parameters:
   ```bash
   muconeup --config config.json --output muc1_simulated.fa --structure-output muc1_struct.txt
   ```
3. Inspect the resulting:
   - **`muc1_simulated.fa`**: multi-FASTA file of haplotype sequences.  
   - **`muc1_struct.txt`**: textual representation of each haplotype’s chain of repeats.

---

## Usage

Below are the available **command-line arguments**. Use `--help` for more details:

```bash
muconeup --help
```

### Command-Line Arguments

| Argument                 | Description                                                                                                                                                                                                                          |
|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--config <path>`       | **Required**. Path to the JSON config file containing repeats, probabilities, constants, length model, mutations, etc.                                                                                                               |
| `--output <filename>`   | Filename for the **simulated FASTA** output. Defaults to `muc1_simulated.fa`.                                                                                                                                                       |
| `--structure-output <filename>` | (Optional) Text file listing the chain of repeats for each haplotype.                                                                                                                                |
| `--num-haplotypes N`     | Number of haplotypes to simulate. Typically `2` for diploid. Defaults to 2.                                                                                                                                                         |
| `--fixed-lengths <ints>` | One or more integers specifying the VNTR repeat length for each haplotype. If only one integer is provided, it applies to **all** haplotypes. If multiple are provided, the number must match `--num-haplotypes`. Otherwise, random. |
| `--seed <int>`          | Random seed for reproducible simulations (affects both VNTR building and random mutation target selection, if used).                                                                                                                  |
| `--mutation-name <str>` | (Optional) Name of a mutation from the config to apply. If omitted, no mutation is applied.                                                                                                                                           |
| `--mutation-targets <pairs>` | (Optional) One or more `haplotype_index,repeat_index` pairs (1-based). E.g., `1,5 2,7`. If provided, each pair indicates which haplotype and which repeat to mutate. If **omitted**, the mutation is randomly applied once to a valid repeat. |

### Example Commands

1. **Generate two diploid haplotypes** with random VNTR lengths based on the config’s distribution:
   ```bash
   muconeup \
     --config config.json \
     --output muc1_sim.fa \
     --structure-output muc1_struct.txt
   ```

2. **Force a fixed length of 60 repeats** for each haplotype:
   ```bash
   muconeup \
     --config config.json \
     --output muc1_fixed.fa \
     --structure-output muc1_fixed.txt \
     --num-haplotypes 2 \
     --fixed-lengths 60
   ```

3. **Apply a known mutation** (`dupC`) at a specific repeat (haplotype #1, repeat #5):
   ```bash
   muconeup \
     --config config.json \
     --output muc1_mutated.fa \
     --structure-output muc1_mutated.txt \
     --mutation-name dupC \
     --mutation-targets 1,5
   ```

4. **Apply a known mutation** (`snpA`) to a **random allowed** repeat:
   ```bash
   muconeup \
     --config config.json \
     --output muc1_random.fa \
     --mutation-name snpA
   ```

---

## Project Structure and Logic

The **muc_one_up** Python package is organized into modules. Here is a brief summary:

```
muc_one_up/
├── cli.py           # Defines the main CLI logic and argument parsing
├── config.py        # Loads and validates the JSON config
├── distribution.py  # Samples the target VNTR length from a specified distribution
├── fasta_writer.py  # Helper (optional) for writing FASTA (unused in the final code, but example)
├── mutate.py        # Logic to apply specified mutations to haplotypes
├── probabilities.py # Weighted random picks for next repeats, etc.
├── simulate.py      # Core simulation code to build haplotypes:
│                   #    - chain repeats with given probabilities
│                   #    - forcibly insert final block 6->7->8->9
└── __init__.py      # Package initialization, version info
```

### High-Level Logic

1. **CLI** (`cli.py`):  
   - Parses user arguments (config path, output file, lengths, mutation name/targets, etc.).  
   - Loads the config.  
   - Calls **`simulate_diploid`** to build the haplotypes.  
   - (Optional) calls **`apply_mutations`** if user requested a specific mutation.  
   - Writes the resulting haplotypes to a multi-FASTA and optionally a textual structure file.

2. **Simulation** (`simulate.py`):  
   - For each haplotype, picks a **target_length** (either from `--fixed-lengths` or from a probability distribution in `config["length_model"]`).  
   - Builds the VNTR chain by sampling repeats from `config["probabilities"]`, **excluding** end repeats (`6`, `6p`, `9`) until the final 4 repeats.  
   - Forces the block `6/6p -> 7 -> 8 -> 9` at `(target_length - 4)`.  
   - Appends **left** and **right** constants from config.  

3. **Mutations** (`mutate.py`):  
   - If a user requests `--mutation-name`, either:  
     - Apply at user-supplied `--mutation-targets`, or  
     - Randomly pick **one** haplotype + repeat that is “allowed” to mutate.  
   - A “mutation” can **insert**, **delete**, or **replace** bases within that repeat.  
   - If the current repeat symbol is **not** in the mutation’s `allowed_repeats`, the code **forces** a new symbol (randomly picked from `allowed_repeats`).  
   - Rebuilds the haplotype sequence to keep it consistent.  
   - Adds an “m” suffix (e.g. `Xm`) to the mutated repeat symbol in the structure output.

4. **Config** (`config.py`):  
   - Loads a JSON file with keys:
     - **`repeats`**: A dictionary mapping repeat symbols to sequences.  
     - **`constants`**: left and right flanks of MUC1.  
     - **`probabilities`**: For each repeat symbol, the probability distribution of valid “next” repeats.  
     - **`length_model`**: Info for random length sampling (e.g., normal distribution, min/max).  
     - **`mutations`**: Named mutations, each with `allowed_repeats` and a list of changes.  

---

## Config File Layout

A simplified example:

```jsonc
{
  "repeats": {
    "1": "CACAGCATTCTTCTC...", 
    "2": "CTGAGTGGTGGAGGA...", 
    // ...
    "9": "TGAGCCTGATGCAGA..."
  },
  "constants": {
    "left": "ACGTACGTACGTACGT",
    "right": "TGCAAGCTTTGCAAGC"
  },
  "probabilities": {
    "1": { "2": 1.0 },
    "2": { "3": 1.0 },
    "3": { "4": 1.0 },
    "4": { "5": 1.0 },
    "5": { "X": 0.2, "C": 0.8 },
    // ...
    "9": { "END": 1.0 }
  },
  "length_model": {
    "distribution": "normal",
    "min_repeats": 20,
    "max_repeats": 130,
    "mean_repeats": 70,
    "median_repeats": 65
  },
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X", "C", "B", "A"], 
      "changes": [
        {
          "type": "insert",
          "start": 2,
          "end": 3,
          "sequence": "G"
        }
      ]
    }
    // Additional mutations...
  }
}
```

---

## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

