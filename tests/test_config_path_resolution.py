"""Tests for config path resolution (#80)."""

import json

from muc_one_up.config import load_config_raw


def test_resolves_relative_paths_from_config_dir(tmp_path, monkeypatch):
    """Relative paths in config are resolved relative to config file directory."""
    ref_dir = tmp_path / "reference"
    ref_dir.mkdir()
    model = ref_dir / "test.model"
    model.write_text("model")

    config = {
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": "reference/test.model",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    # Load from a different working directory
    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["pacbio_params"]["model_file"] == str(model)


def test_preserves_absolute_paths(tmp_path):
    """Absolute paths in config are not modified."""
    config = {
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": "/absolute/path/to/model.file",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    loaded = load_config_raw(str(config_path))
    assert loaded["pacbio_params"]["model_file"] == "/absolute/path/to/model.file"


def test_preserves_paths_that_exist_from_cwd(tmp_path, monkeypatch):
    """Paths that already resolve from cwd are not changed."""
    config = {
        "read_simulation": {
            "human_reference": "config.json",
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir(tmp_path)
    loaded = load_config_raw(str(config_path))
    assert loaded["read_simulation"]["human_reference"] == "config.json"


def test_resolves_reference_genomes_fasta_path(tmp_path, monkeypatch):
    """reference_genomes.*.fasta_path is resolved relative to config directory."""
    ref_dir = tmp_path / "reference" / "hg38"
    ref_dir.mkdir(parents=True)
    fa = ref_dir / "hg38.fa"
    fa.write_text(">hg38\nACGT\n")

    config = {
        "reference_genomes": {
            "hg38": {
                "fasta_path": "reference/hg38/hg38.fa",
                "vntr_region": "chr1:1-100",
            }
        }
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["reference_genomes"]["hg38"]["fasta_path"] == str(fa.resolve())


def test_resolves_nanosim_training_data_path_with_model_profile(tmp_path, monkeypatch):
    """nanosim_params.training_data_path prefix is resolved when <prefix>_model_profile exists."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()
    profile = models_dir / "ecoli_model_profile"
    profile.write_text("profile data")

    config = {
        "nanosim_params": {
            "training_data_path": "models/ecoli",
            "coverage": 30,
        }
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["nanosim_params"]["training_data_path"] == str((models_dir / "ecoli").resolve())


def test_resolves_capture_bed(tmp_path, monkeypatch):
    """read_simulation.vntr_capture_efficiency.capture_bed is resolved relative to config dir."""
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t1\t100\n")

    config = {
        "read_simulation": {
            "vntr_capture_efficiency": {
                "capture_bed": "targets.bed",
            }
        }
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["read_simulation"]["vntr_capture_efficiency"]["capture_bed"] == str(bed.resolve())


def test_resolves_flank_fasta(tmp_path, monkeypatch):
    """ont_fragment_params.flank_fasta is resolved relative to config dir."""
    flank = tmp_path / "flanks.fa"
    flank.write_text(">left\nACGT\n>right\nTGCA\n")

    config = {
        "ont_fragment_params": {
            "flank_fasta": "flanks.fa",
        }
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    monkeypatch.chdir("/tmp")
    loaded = load_config_raw(str(config_path))
    assert loaded["ont_fragment_params"]["flank_fasta"] == str(flank.resolve())
