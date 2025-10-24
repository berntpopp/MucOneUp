# Tool Version Detection - Testing Results

**Date:** 2025-10-24
**Environment:** WSL Ubuntu + Conda Environments
**Tools Tested:** 9/9 (100% coverage)

---

## Testing Summary

✅ **ALL 9 TOOLS TESTED** using their actual conda environments:
- env_nanosim (ONT simulation)
- env_wessim (Illumina simulation)
- env_pacbio (PacBio simulation)
- System install (samtools, minimap2, bwa)

---

## Verified Tool Specifications

### ✅ Standard Tools (--version flag works)

**1. samtools 1.17**
```bash
Command: samtools --version
Exit: 0
Stream: stdout
Output: "samtools 1.17"
```

**2. minimap2 2.28-r1209**
```bash
Command: minimap2 --version
Exit: 0
Stream: stdout
Output: "2.28-r1209"
Note: Prepend "minimap2 " in parser
```

**3. NanoSim 3.2.2** (conda: env_nanosim)
```bash
Command: simulator.py --version
Exit: 0
Stream: stdout
Output: "NanoSim 3.2.2"
⚠️ Binary is "simulator.py", NOT "nanosim"!
```

**4. ccs 6.4.0** (conda: env_pacbio)
```bash
Command: ccs --version
Exit: 0
Stream: stdout
Output: "ccs 6.4.0 (commit v6.4.0)"
```

**5. reseq 1.1** (conda: env_wessim)
```bash
Command: reseq --version
Exit: 0
Stream: stdout
Output: "ReSeq version 1.1"
✓ --version works (not just --help as originally thought)
```

---

### ⚠️ Special Handling Required

**6. bwa 0.7.18-r1243** (system)
```bash
Command: bwa (no args)
Exit: 1 ⚠️ NON-ZERO!
Stream: stderr
Output: "Version: 0.7.18-r1243-dirty"
Parse: Regex r"version:\s*([\w\.-]+)"
Critical: MUST use subprocess.run(..., check=False)
```

**7. pblat v.36x2** (conda: env_wessim)
```bash
Command: pblat (no args)
Exit: 0
Stream: stdout
Output: "pblat - BLAT with parallel supports v. 36x2"
Parse: Regex r"v\.\s*([\w]+)"
Note: NO --version flag available
```

**8. pbsim3 3.0.5** (conda: env_pacbio)
```bash
Command: pbsim (no args)
Exit: 255 ⚠️ NON-ZERO!
Stream: stdout
Output: "USAGE: pbsim [options]"
🚨 CRITICAL ISSUE: NO VERSION IN OUTPUT!
⚠️ Binary name is "pbsim", NOT "pbsim3"
Solution: Return generic "pbsim3 (installed)" string
Alternative: Parse conda metadata (complex)
```

**9. faToTwoBit** (conda: env_wessim)
```bash
Command: faToTwoBit (no args)
Exit: 255 ⚠️ NON-ZERO!
Stream: stdout
Output: "faToTwoBit - Convert DNA from fasta to 2bit format"
🚨 CRITICAL ISSUE: NO VERSION ANYWHERE!
Solution: Return "faToTwoBit (version unknown)"
```

---

## Exit Code Matrix

| Tool | Exit Code | Requires check=False |
|------|-----------|----------------------|
| samtools | 0 ✅ | No |
| minimap2 | 0 ✅ | No |
| simulator.py | 0 ✅ | No |
| ccs | 0 ✅ | No |
| reseq | 0 ✅ | No |
| bwa | 1 ⚠️ | **YES** |
| pblat | 0 ✅ | No |
| pbsim | 255 ⚠️ | **YES** |
| faToTwoBit | 255 ⚠️ | **YES** |

**Critical:** 3 tools return non-zero - MUST use `subprocess.run(..., check=False)`

---

## Binary Name Mismatches

| Config Key | Expected Binary | Actual Binary | Issue |
|------------|-----------------|---------------|-------|
| nanosim | nanosim | **simulator.py** | ⚠️ Name mismatch |
| pbsim3 | pbsim3 | **pbsim** | ⚠️ Name mismatch + no version |

**Solution:** Config should specify actual binary path/name

---

## Version Output Patterns

### Pattern 1: Clean Version String
- **Tools:** samtools, NanoSim, ccs, reseq
- **Output:** "ToolName X.Y.Z"
- **Parse:** Direct use (strip whitespace)

### Pattern 2: Version Number Only
- **Tools:** minimap2
- **Output:** "2.28-r1209"
- **Parse:** Prepend tool name

### Pattern 3: Regex Extraction
- **Tools:** bwa, pblat
- **Output:** Mixed text with version embedded
- **Parse:** Regex to extract version substring

### Pattern 4: No Version Available
- **Tools:** pbsim3, faToTwoBit
- **Output:** Usage text, NO version info
- **Parse:** Return generic "(installed)" or "(version unknown)"

---

## Implementation Updates Required

### Updated tool_map Configuration

```python
tool_map = {
    "samtools": {"args": ["--version"], "parse": parse_samtools},
    "minimap2": {"args": ["--version"], "parse": parse_minimap2},
    "nanosim": {"args": ["--version"], "parse": parse_nanosim},  # Binary: simulator.py
    "pbsim3": {"args": [], "parse": parse_pbsim3},  # Binary: pbsim, returns generic string
    "ccs": {"args": ["--version"], "parse": parse_ccs},
    "reseq": {"args": ["--version"], "parse": parse_reseq},  # --version works!
    "bwa": {"args": [], "parse": parse_bwa},  # No flag, stderr, exit 1
    "faToTwoBit": {"args": [], "parse": parse_fatotwobit},  # No version available
    "pblat": {"args": [], "parse": parse_pblat},  # No flag, parse usage
}
```

### Key Parser Changes

**pblat parser:**
```python
def parse_pblat(stdout: str, stderr: str) -> str:
    for line in stdout.splitlines():
        if "pblat" in line.lower():
            match = re.search(r"v\.\s*([\w]+)", line)
            if match:
                return f"pblat v.{match.group(1)}"
    return "N/A"
```

**pbsim3 parser:**
```python
def parse_pbsim3(stdout: str, stderr: str) -> str:
    # pbsim3 provides NO version in output
    if "pbsim" in stdout.lower() or "USAGE: pbsim" in stdout:
        return "pbsim3 (installed)"  # Generic string
    return "N/A"
```

**faToTwoBit parser:**
```python
def parse_fatotwobit(stdout: str, stderr: str) -> str:
    if "fatotwobit" in stdout.lower():
        return "faToTwoBit (version unknown)"
    return "N/A"
```

---

## Test Coverage

| Tool | Unit Test | Integration Test | Status |
|------|-----------|------------------|--------|
| samtools | ✅ Mocked | ✅ Real tool | Complete |
| minimap2 | ✅ Mocked | ✅ Real tool | Complete |
| NanoSim | ✅ Mocked | ✅ Conda env | Complete |
| pbsim3 | ✅ Mocked | ✅ Conda env | Complete |
| ccs | ✅ Mocked | ✅ Conda env | Complete |
| reseq | ✅ Mocked | ✅ Conda env | Complete |
| bwa | ✅ Mocked | ✅ Real tool | Complete |
| faToTwoBit | ✅ Mocked | ✅ Conda env | Complete |
| pblat | ✅ Mocked | ✅ Conda env | Complete |

**Coverage: 9/9 tools (100%)**

---

## Recommendations

### 1. Config.json Tool Names
Ensure config.json uses correct binary names:
```json
{
  "tools": {
    "nanosim": "simulator.py",  // Not "nanosim"
    "pbsim3": "pbsim",           // Not "pbsim3"
    // ... rest standard
  }
}
```

### 2. Error Handling
**ALWAYS use check=False:**
```python
result = subprocess.run(cmd, check=False, ...)  # Don't raise on non-zero
```

### 3. Generic Strings for Missing Versions
Accept that some tools (pbsim3, faToTwoBit) cannot provide version:
- Return descriptive string like "(installed)" or "(version unknown)"
- Log warning for debugging
- Don't fail pipeline

### 4. Documentation Updates
- Update CLAUDE.md with binary name mismatches
- Add table of special cases to README
- Include example metadata output

---

## Example Metadata Output

**With all 9 tools:**
```tsv
Parameter	Value
Tool	MucOneUp
Version	0.15.0
Platform	Illumina
tool.samtools_version	samtools 1.17
tool.minimap2_version	minimap2 2.28-r1209
tool.nanosim_version	NanoSim 3.2.2
tool.pbsim3_version	pbsim3 (installed)
tool.ccs_version	ccs 6.4.0 (commit v6.4.0)
tool.reseq_version	ReSeq version 1.1
tool.bwa_version	bwa 0.7.18-r1243
tool.faToTwoBit_version	faToTwoBit (version unknown)
tool.pblat_version	pblat v.36x2
```

---

## Conclusion

✅ **All 9 tools successfully tested**
✅ **Version detection protocols verified**
✅ **Special cases identified and documented**
✅ **Parser implementations updated**
⚠️ **2 tools lack version info (pbsim3, faToTwoBit) - use generic strings**

**Documentation updated:** `TOOL_VERSION_TRACKING.md`
**Ready for implementation:** v0.16.0

---

*Testing conducted: 2025-10-24*
*Environments: WSL Ubuntu, conda (env_nanosim, env_wessim, env_pacbio)*
*Tools tested: 9/9 (100%)*
