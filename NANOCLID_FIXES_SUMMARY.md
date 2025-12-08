# NanoCliD Pipeline Fixes Summary

**Date:** December 5, 2025  
**Author:** GitHub Copilot Session  
**Branch:** my-local-branch

---

## 1. Container Path Issue (CRITICAL FIX)

### Problem
The pipeline was failing with `dorado.sif: no such file or directory` because container paths were not being resolved correctly.

### Root Cause
- Multiple `.snk` files had inconsistent `CONTAINERS_PATH` definitions
- Some used `workflow.singularity_prefix`, others used `config.get("singularity-prefix", "")`
- The `-C/--containersFolder` argument was optional with a wrong default

### Solution

#### 1.1 Made containersFolder a required argument
**File:** `nanoclid_FO.py`
```python
# Changed from optional to required
run_parser.add_argument("-C", "--containersFolder", required=True, 
    help="Path to singularity containers folder (required).")
```

#### 1.2 Simplified CONTAINERS_PATH in all Snakemake files

**File:** `workflow/Snakefile` (line ~21)
```python
# CONTAINERS_PATH is set from the required -C/--containersFolder argument
CONTAINERS_PATH = config["containers_path"]
```

**File:** `workflow/rules/common.snk` (line ~14)
```python
# CONTAINERS_PATH is set from the required -C/--containersFolder argument
CONTAINERS_PATH = config["containers_path"]
```

**File:** `workflow/rules/demultiplexing.snk` (line ~20)
```python
# CONTAINERS_PATH is set from the required -C/--containersFolder argument
CONTAINERS_PATH = config["containers_path"]
```

#### 1.3 Fixed profile update to set singularity-prefix
**File:** `nanoclid_FO.py` - `__updateProfile()` method
```python
# Set singularity-prefix from the containersFolder argument
config["singularity-prefix"] = self.containersFolder
```

#### 1.4 Created missing symlink
```bash
cd /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/singularity
ln -s dorado_v0.7.2.sif dorado.sif
```

---

## 2. Test Command Fix

### Problem
The `test` subcommand was creating a `NanoClid` instance without passing `containersFolder`.

### Solution
**File:** `nanoclid_FO.py` (line ~897)
```python
# Before:
NanoClid(profile = "standalone")._runTest(...)

# After:
NanoClid(profile = "standalone", containersFolder=args.containersFolder)._runTest(...)
```

---

## 3. hg19.genome File Fix

### Problem
The `hg19.genome` file was incorrectly formatted (contained region coordinates instead of chromosome sizes), causing `bedtools complement` to fail.

### Solution
Downloaded proper chromosome sizes:
```bash
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes \
  -o /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/nanoclid_ref/hg19.genome
```

**Expected format:**
```
chr1    249250621
chr2    243199373
chr3    198022430
...
```

---

## 4. Dorado Model Configuration

### Problem
Dorado basecaller fails with: `toml::parse: file open error -> dna_r9.4.1_e8_hac@v3.3/config.toml`

### Cause
Dorado needs either:
1. Pre-downloaded model with full path, or
2. Ability to auto-download (requires internet access in container)

### Solution
Download the model and update config to use full path:

```bash
# Create models directory
mkdir -p /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/models

# Download model
singularity exec \
  -B /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD \
  /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/singularity/dorado.sif \
  dorado download --model dna_r9.4.1_e8_hac@v3.3 \
  --directory /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/models
```

Then update config template to use full path:
```yaml
dorado:
  model: /path/to/models/dna_r9.4.1_e8_hac@v3.3
```

---

## 5. Command Reference

### Running the Test
```bash
cd /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD
source venv/bin/activate
python nanoclid_FO.py test \
  -p standalone \
  -R /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/nanoclid_ref \
  -C /hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD/singularity
```

### Running the Pipeline
```bash
python nanoclid_FO.py run \
  -g hg19 \
  -I /path/to/input/data \
  -r RUN_ID \
  -p standalone \
  -R /path/to/reference \
  -S /path/to/samplesheet.csv \
  -C /path/to/singularity/containers  # REQUIRED
```

### Testing Containers Manually
```bash
# Test bedtools
singularity exec /path/to/singularity/bedtools.sif bedtools --version

# Test dorado
singularity exec /path/to/singularity/dorado.sif dorado --version

# Test blue-crab
singularity exec /path/to/singularity/bluecrab.sif blue-crab --version

# Test with actual command and bindings
singularity exec \
  -B /path/to/data,/path/to/ref \
  /path/to/singularity/bedtools.sif \
  bedtools complement -i input.bed -g genome.genome
```

---

## 6. Files Modified

| File | Changes |
|------|---------|
| `nanoclid_FO.py` | Made `-C` required, fixed test command, snakemakeBin auto-detection |
| `workflow/Snakefile` | Simplified CONTAINERS_PATH |
| `workflow/rules/common.snk` | Simplified CONTAINERS_PATH |
| `workflow/rules/demultiplexing.snk` | Simplified CONTAINERS_PATH |
| `profiles/externe/*/config.yaml` | Updated singularity-prefix handling |
| `singularity/dorado.sif` | Created symlink to dorado_v0.7.2.sif |
| `nanoclid_ref/hg19.genome` | Replaced with proper chromosome sizes |

---

## 7. Remaining Issues

1. **Dorado model path**: Need to download models and update config template
2. **GPU requirement**: Dorado basecalling requires GPU access for production use
3. **SyntaxWarnings**: Minor escape sequence warnings in preprocessing.snk (cosmetic)

---

## 8. Quick Troubleshooting

| Error | Cause | Solution |
|-------|-------|----------|
| `dorado.sif: no such file` | Missing container or wrong path | Check `-C` argument points to singularity folder |
| `Exit status 255` | Container execution failed | Test container manually with `singularity exec` |
| `chromosome chr1 does not exist` | Wrong genome file format | Use proper `.genome` file with chr sizes |
| `toml::parse: file open error` | Dorado model not found | Download model or use full path |

---

*Generated from NanoCliD debugging session - December 5, 2025*
