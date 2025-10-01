# Batch Processing Analysis: Series Generation Gap

**Date:** 2025-10-01
**Status:** ✅ IMPLEMENTED
**Priority:** COMPLETED

## Executive Summary

**The Problem (RESOLVED):** When `simulate` generates a series of files (e.g., with `--simulate-series`), downstream commands (`reads`, `analyze`) could not process them in batch. They only accepted ONE file at a time, and `pipeline` only processed the FIRST file (`.001`).

**Impact (RESOLVED):** Users had to manually write shell loops to process series.

**Solution (IMPLEMENTED):**
1. ✅ Removed `pipeline` command (unnecessary complexity)
2. ✅ Added multi-file support to all `reads` and `analyze` commands following Unix philosophy
3. ✅ Commands now accept multiple files via `nargs=-1` (Click variadic arguments)
4. ✅ Full backward compatibility maintained

---

## Current Behavior

### How Series Generation Works

When you run:
```bash
muconeup --config X simulate --fixed-lengths 20-40 --simulate-series 1
```

This generates **21 files**:
```
sample.001.simulated.fa  # 20 repeats
sample.002.simulated.fa  # 21 repeats
sample.003.simulated.fa  # 22 repeats
...
sample.021.simulated.fa  # 40 repeats
```

### The Gap: Downstream Commands Don't Support Batch

**1. `reads` command** (lines 285-376 in click_main.py):
```python
@reads.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
def illumina(ctx, input_fasta, ...):
    # Only processes ONE file
```

**2. `analyze` command** (lines 425-509 in click_main.py):
```python
@analyze.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
def orfs(ctx, input_fasta, ...):
    # Only processes ONE file
```

**3. `pipeline` command** (lines 636-713 in click_main.py):
```python
def pipeline(ctx, ...):
    # Line 672: HARDCODED to use only .001 file
    fasta_file = Path(kwargs["out_dir"]) / f"{kwargs['out_base']}.001.simulated.fa"
```

### Current Workaround (Manual Shell Loop)

Users must write this themselves:
```bash
# Step 1: Generate series
muconeup --config X simulate --fixed-lengths 20-40 --simulate-series 1 --out-base sample

# Step 2: Manually loop over files
for file in sample.*.simulated.fa; do
  muconeup --config X reads illumina "$file" --out-base "${file%.fa}_reads"
  muconeup --config X analyze orfs "$file" --out-base "${file%.fa}_orfs"
done
```

---

## Research: CLI Best Practices

### Unix Philosophy (clig.dev)

> **Core Principle:** "Small, simple programs with clean interfaces that can be combined to build larger systems."

> **McIlroy's Rule:** "Expect the output of every program to become the input to another, as yet unknown, program."

**Key Guidelines:**

1. ✅ **Multiple arguments are fine** for simple actions against multiple files
   - Example: `rm file1.txt file2.txt file3.txt`

2. ✅ **Support glob patterns** seamlessly
   - Example: `rm *.txt`

3. ✅ **Composability over built-in orchestration**
   - Tools should be designed to work with xargs, parallel, shell loops
   - Don't force users into a specific workflow

4. ✅ **Plain text interfaces** for easy piping
   - Line-based output for composability

### Click Framework Patterns

From Context7 Click documentation:

**Variadic Arguments (nargs=-1):**
```python
@click.command()
@click.argument('src', nargs=1)
@click.argument('dsts', nargs=-1)  # Accept multiple files
def copy(src: str, dsts: tuple[str, ...]):
    for destination in dsts:
        click.echo(f"Copy {src} to folder {destination}")
```

**Glob Pattern Handling:**
- Click automatically expands glob patterns on Windows
- On Unix, shell expands patterns before passing to Click
- Best practice: Accept multiple arguments, let shell handle expansion

### Parallel Processing Best Practices

From xargs/parallel research:

1. ✅ **Design for xargs compatibility**
   - Single-file commands compose naturally with xargs
   - Example: `ls *.fa | xargs -I {} muconeup reads illumina {}`

2. ✅ **Provide built-in parallelism only when it adds value**
   - Don't reinvent xargs/parallel
   - Add parallelism if the tool has internal state/caching benefits

3. ✅ **Report progress thoughtfully**
   - Parallel output is hard to read
   - Consider `--quiet` mode for batch operations

---

## Recommendations

### Solution Architecture: Three Complementary Approaches

Following Unix philosophy: **Provide flexible tools, let users compose them.**

#### **Option 1: Accept Multiple Files (RECOMMENDED)**

**Pros:**
- ✅ Intuitive: `muconeup reads illumina *.simulated.fa`
- ✅ Works with shell globs naturally
- ✅ Common pattern (cp, rm, mv all do this)
- ✅ Sequential processing (predictable output)

**Cons:**
- ⚠️ Processing is sequential (but this is often desired)
- ⚠️ Need to handle output naming carefully

**Implementation:**
```python
@reads.command()
@click.argument('input_fastas', nargs=-1, required=True,
                type=click.Path(exists=True, dir_okay=False))
def illumina(ctx, input_fastas, out_dir, coverage, threads, ...):
    """Simulate Illumina reads from one or more FASTA files.

    Examples:
      # Single file
      muconeup reads illumina sample.001.fa --out-base reads

      # Multiple files (explicit)
      muconeup reads illumina sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup reads illumina sample.*.simulated.fa
    """
    for fasta_file in input_fastas:
        # Auto-generate output base from input filename if not specified
        base_name = Path(fasta_file).stem  # e.g., "sample.001.simulated"
        out_base = f"{base_name}_reads"

        logging.info(f"Processing {fasta_file} -> {out_base}")
        # Run simulation for this file
        simulate_reads_pipeline(config, fasta_file, out_base, ...)
```

#### **Option 2: Input Directory Support**

**Pros:**
- ✅ Clean for large batches
- ✅ No argument length limits
- ✅ Easy to specify "all files in this folder"

**Cons:**
- ⚠️ Need to define matching pattern
- ⚠️ Less flexible than explicit glob

**Implementation:**
```python
@reads.command()
@click.argument('input', type=click.Path(exists=True))  # File OR directory
@click.option('--pattern', default='*.simulated.fa',
              help='Glob pattern for directory input')
def illumina(ctx, input, pattern, ...):
    """Simulate reads from FASTA file(s).

    INPUT can be:
      - A single FASTA file
      - A directory (processes all files matching --pattern)

    Examples:
      # Single file
      muconeup reads illumina sample.001.fa

      # Directory with default pattern
      muconeup reads illumina output_dir/

      # Directory with custom pattern
      muconeup reads illumina output_dir/ --pattern '*.fa'
    """
    input_path = Path(input)

    if input_path.is_file():
        files = [input_path]
    else:
        files = sorted(input_path.glob(pattern))
        if not files:
            raise click.UsageError(f"No files matching '{pattern}' in {input}")

    for fasta_file in files:
        # Process each file
        ...
```

#### **Option 3: xargs-Friendly Design (CURRENT - Keep This!)**

**Current design is ALREADY xargs-friendly!** Single-file commands compose perfectly:

```bash
# Using xargs (sequential)
ls sample.*.simulated.fa | xargs -I {} muconeup --config X reads illumina {} --out-base {}_reads

# Using GNU parallel (parallel execution)
ls sample.*.simulated.fa | parallel muconeup --config X reads illumina {} --out-base {.}_reads

# Using shell loop (most portable)
for file in sample.*.simulated.fa; do
  muconeup --config X reads illumina "$file" --out-base "${file%.fa}_reads"
done
```

**Recommendation:** KEEP single-file support, ADD multi-file support.

---

## Specific Recommendations for MucOneUp

### Phase 1: Add Multi-File Support (Quick Win)

**Commands to Update:**
1. ✅ `reads illumina` - Accept multiple input files
2. ✅ `reads ont` - Accept multiple input files
3. ✅ `analyze orfs` - Accept multiple input files
4. ✅ `analyze stats` - Accept multiple input files

**Implementation Pattern:**
```python
@click.argument('input_fastas', nargs=-1, required=True,
                type=click.Path(exists=True, dir_okay=False))
@click.option('--out-base', default=None,
              help='Base name for outputs (auto-generated if multiple files)')
def command(ctx, input_fastas, out_base, ...):
    if len(input_fastas) > 1 and out_base:
        # Warn: out_base ignored for multiple files
        logging.warning("--out-base ignored for multiple files (auto-generated)")

    for fasta_file in input_fastas:
        if out_base and len(input_fastas) == 1:
            actual_out_base = out_base
        else:
            # Auto-generate: sample.001.simulated.fa -> sample.001.simulated_reads
            actual_out_base = Path(fasta_file).stem + "_suffix"

        # Process file
        ...
```

### Phase 2: Fix Pipeline Series Handling

**Current Bug:**
```python
# Line 672 - HARDCODED to .001
fasta_file = Path(kwargs["out_dir"]) / f"{kwargs['out_base']}.001.simulated.fa"
```

**Fix Option A: Process ALL generated files**
```python
def pipeline(ctx, ...):
    # Run simulate (may generate multiple files)
    ctx.invoke(simulate, ...)

    # Find ALL generated files
    out_dir = Path(kwargs["out_dir"])
    pattern = f"{kwargs['out_base']}.*.simulated.fa"
    generated_files = sorted(out_dir.glob(pattern))

    if not generated_files:
        raise RuntimeError(f"No simulated files found matching {pattern}")

    logging.info(f"Found {len(generated_files)} simulated files to process")

    # Process each file
    for fasta_file in generated_files:
        base_name = fasta_file.stem

        if with_reads:
            ctx.invoke(illumina, input_fasta=str(fasta_file),
                      out_base=f"{base_name}_reads", ...)

        if with_orfs:
            ctx.invoke(orfs, input_fasta=str(fasta_file),
                      out_base=f"{base_name}_orfs", ...)
```

**Fix Option B: Add explicit flag**
```python
@click.option('--process-series/--process-first', default=False,
              help='Process all series files or just first (default: first)')
def pipeline(ctx, process_series, ...):
    if process_series:
        # Process all files
        ...
    else:
        # Current behavior: process only .001
        ...
```

### Phase 3: Documentation Updates

**1. Update README.md** - Add batch processing examples:
```markdown
### Batch Processing Series

**Option 1: Shell loop**
```bash
for file in sample.*.simulated.fa; do
  muconeup --config X reads illumina "$file" --out-base "${file%.fa}_reads"
done
```

**Option 2: Multi-file argument**
```bash
muconeup --config X reads illumina sample.*.simulated.fa
```

**Option 3: Parallel with GNU parallel**
```bash
ls sample.*.fa | parallel muconeup --config X reads illumina {}
```
```

**2. Update MIGRATION_v2.md** - Add batch processing section

**3. Add examples to command help strings**

---

## Implementation Priority

### Critical (Do Now)
1. ✅ **Fix `pipeline` command** to handle series (Option A recommended)
   - Users expect `--simulate-series` to work with `--with-reads`/`--with-orfs`

### High Priority (Phase 2)
2. ✅ **Add multi-file support to `reads` commands**
   - Most common use case for batch processing

3. ✅ **Add multi-file support to `analyze` commands**
   - Natural extension

### Nice-to-Have (Phase 3)
4. ⚠️ **Add `--pattern` support for directory input**
   - Convenience feature, not essential

5. ⚠️ **Add `--parallel` flag with thread pool**
   - Only if users request it (xargs/parallel already exist)

---

## Testing Strategy

### Test Cases

1. **Single file (backward compatibility)**
   ```bash
   muconeup reads illumina sample.001.fa --out-base test
   ```

2. **Multiple explicit files**
   ```bash
   muconeup reads illumina sample.001.fa sample.002.fa sample.003.fa
   ```

3. **Glob pattern (shell expansion)**
   ```bash
   muconeup reads illumina sample.*.simulated.fa
   ```

4. **Pipeline with series**
   ```bash
   muconeup pipeline --fixed-lengths 20-40 --simulate-series 1 \
     --with-reads illumina --with-orfs
   # Should process ALL 21 generated files
   ```

5. **Empty glob (error handling)**
   ```bash
   muconeup reads illumina nonexistent.*.fa
   # Should fail with clear error message
   ```

---

## Examples: Before vs After

### Example 1: Analyze Series of ORFs

**Before (manual loop required):**
```bash
# Step 1: Generate series
muconeup simulate --fixed-lengths 20-40 --simulate-series 1 --out-base sample

# Step 2: Manually loop
for file in sample.*.simulated.fa; do
  muconeup analyze orfs "$file" --out-base "${file%.fa}_orfs"
done
```

**After (multi-file support):**
```bash
# Step 1: Generate series
muconeup simulate --fixed-lengths 20-40 --simulate-series 1 --out-base sample

# Step 2: Process all at once
muconeup analyze orfs sample.*.simulated.fa
```

**After (pipeline with series fix):**
```bash
# All in one command
muconeup pipeline --fixed-lengths 20-40 --simulate-series 1 \
  --out-base sample --with-orfs
# Automatically processes all 21 files
```

### Example 2: Read Simulation for All Series

**Before:**
```bash
muconeup simulate --fixed-lengths 20-40 --simulate-series 1 --out-base sample

for file in sample.*.simulated.fa; do
  base="${file%.simulated.fa}"
  muconeup reads illumina "$file" --out-base "${base}_reads" --coverage 30
done
```

**After:**
```bash
muconeup simulate --fixed-lengths 20-40 --simulate-series 1 --out-base sample
muconeup reads illumina sample.*.simulated.fa --coverage 30
```

---

## Decision Matrix

| Approach | Unix Philosophy | User Convenience | Implementation Effort | Backward Compatible |
|----------|----------------|------------------|----------------------|---------------------|
| **Option 1: Multi-file args** | ✅ Excellent | ✅ High | 🟡 Medium | ✅ Yes (nargs=-1 accepts 1+) |
| **Option 2: Directory support** | ✅ Good | 🟡 Medium | 🟡 Medium | ✅ Yes (different arg type) |
| **Option 3: Keep xargs-friendly** | ✅ Perfect | 🔴 Low | ✅ None (done) | ✅ Yes |
| **Built-in parallel** | 🔴 Poor (reinvents wheel) | 🟡 Medium | 🔴 High | ⚠️ Complex |

**Recommendation:** Implement **Options 1 + 3** (keep single-file, add multi-file).

---

## Unix Philosophy Checklist

✅ **Do one thing well** - Each command has clear single responsibility
✅ **Composability** - Works with xargs, parallel, shell loops
✅ **Text interfaces** - File paths are plain text
✅ **No surprises** - Glob expansion follows shell conventions
✅ **Expect composition** - Users can chain commands in unforeseen ways
✅ **Keep it simple** - Multi-file support is straightforward addition

---

## Conclusion

**Current State:** Series generation creates files that can't be batch-processed by downstream commands without manual scripting.

**Recommended Solution:**
1. ✅ **Add multi-file support** to `reads` and `analyze` commands (Unix philosophy: accept multiple arguments)
2. ✅ **Fix `pipeline` command** to process ALL series files, not just `.001`
3. ✅ **Keep xargs-friendly** single-file mode (backward compatible)
4. ✅ **Document batch patterns** in README and migration guide

**Impact:** Makes MucOneUp truly production-ready for parameter sweeps and batch analysis workflows while maintaining Unix philosophy principles.

**Next Steps:**
1. Implement multi-file support (Phase 1 - Critical)
2. Fix pipeline series handling (Phase 1 - Critical)
3. Add tests for batch scenarios (Phase 2)
4. Update documentation (Phase 3)
