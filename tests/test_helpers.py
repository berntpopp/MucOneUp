"""Unit tests for reference download and installation helper scripts.

Tests cover security fixes and enhancements for Issue #155:
1. Command injection prevention (shell=False, argument lists) in create_minimap2_index
2. Preset-specific minimap2 index generation (map-ont, map-hifi) and NanoSim training data path
3. TarSlip protection (filter="data"), HTTPS URL usage, and MD5 checksum verification
4. Safe whitespace splitting in execute_index_command
"""

import gzip
import hashlib
import json
import sys
import tarfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure helpers directory is on sys.path for direct importing
REPO_ROOT = Path(__file__).resolve().parent.parent
HELPERS_DIR = REPO_ROOT / "helpers"
if str(HELPERS_DIR) not in sys.path:
    sys.path.insert(0, str(HELPERS_DIR))

import download_references  # noqa: E402
import install_references  # noqa: E402

# -----------------------------------------------------------------------------
# Problem 1 & 2: Minimap2 index creation, shell=False, presets, and whitespace
# -----------------------------------------------------------------------------


def test_create_minimap2_index_no_shell_and_preset(tmp_path):
    """create_minimap2_index must use shell=False and pass preset arguments safely."""
    ref_file = tmp_path / "reference.fasta"
    ref_file.write_text(">chr1\nACGT\n")

    with patch("subprocess.run") as mock_run:
        result = download_references.create_minimap2_index(
            reference_path=str(ref_file),
            output_dir=str(tmp_path),
            minimap2_cmd="minimap2",
            preset="map-ont",
        )

        expected_index = tmp_path / "reference.fasta.map-ont.mmi"
        assert result == expected_index

        mock_run.assert_called_once()
        args, kwargs = mock_run.call_args
        cmd_args = args[0]

        # Verify command argument list and shell=False
        assert isinstance(cmd_args, list)
        assert kwargs.get("shell") is False
        assert kwargs.get("check") is True
        assert cmd_args == [
            "minimap2",
            "-x",
            "map-ont",
            "-d",
            str(expected_index),
            str(ref_file),
        ]


def test_create_minimap2_index_handles_spaces_and_special_chars(tmp_path):
    """Paths containing whitespace or shell metacharacters must not execute shell commands."""
    space_dir = tmp_path / "path with spaces"
    space_dir.mkdir()
    ref_file = space_dir / "ref; injected && ls.fasta"
    ref_file.write_text(">chr1\nACGT\n")

    with patch("subprocess.run") as mock_run:
        result = download_references.create_minimap2_index(
            reference_path=str(ref_file),
            output_dir=str(space_dir),
            minimap2_cmd="mamba run -n env_nanosim minimap2",
            preset="map-hifi",
        )

        expected_index = space_dir / f"{ref_file.name}.map-hifi.mmi"
        assert result == expected_index

        mock_run.assert_called_once()
        cmd_args, kwargs = mock_run.call_args[0], mock_run.call_args[1]
        cmd = cmd_args[0]

        assert kwargs.get("shell") is False
        # The arguments must be preserved as discrete strings without shell splitting
        assert cmd[:5] == ["mamba", "run", "-n", "env_nanosim", "minimap2"]
        assert "-x" in cmd and "map-hifi" in cmd
        assert str(expected_index) in cmd
        assert str(ref_file) in cmd


def test_create_minimap2_index_skips_when_exists(tmp_path):
    """create_minimap2_index must skip execution if the index file already exists."""
    ref_file = tmp_path / "ref.fasta"
    ref_file.write_text(">chr1\nACGT\n")
    existing_index = tmp_path / "ref.fasta.map-ont.mmi"
    existing_index.write_bytes(b"EXISTING_INDEX")

    with patch("subprocess.run") as mock_run:
        result = download_references.create_minimap2_index(
            reference_path=str(ref_file),
            output_dir=str(tmp_path),
            preset="map-ont",
        )
        assert result == existing_index
        mock_run.assert_not_called()


def test_create_minimap2_index_missing_reference(tmp_path):
    """create_minimap2_index returns None if reference file does not exist."""
    missing_ref = tmp_path / "nonexistent.fa"
    result = download_references.create_minimap2_index(
        reference_path=str(missing_ref),
        output_dir=str(tmp_path),
    )
    assert result is None


# -----------------------------------------------------------------------------
# Problem 2: NanoSim training_data_path and dual preset minimap2 indexing
# -----------------------------------------------------------------------------


def test_download_references_main_training_data_path_and_presets(tmp_path):
    """main() must write training_data_path with '/training' suffix and build both presets."""
    dummy_ref = tmp_path / "hg38.fa"
    dummy_ref.write_text(">chr1\nACGT\n")

    config_data = {
        "constants": {"hg19": {}, "hg38": {}},
        "read_simulation": {"human_reference": str(dummy_ref)},
        "tools": {"minimap2": "minimap2"},
        "nanosim_params": {"training_data_path": "old_path"},
    }
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(config_data))

    dummy_model_dir = tmp_path / "reference" / "nanosim" / "model_dir"
    dummy_model_dir.mkdir(parents=True)

    with (
        patch("download_references.download_nanosim_model", return_value=dummy_model_dir),
        patch("download_references.create_minimap2_index") as mock_create_idx,
        patch(
            "sys.argv",
            [
                "download_references.py",
                "--config",
                str(config_file),
                "--references-dir",
                str(tmp_path / "reference"),
                "--setup-nanosim",
                "--create-minimap2-index",
            ],
        ),
    ):
        download_references.main()

        # Check updated config.json
        updated_config = json.loads(config_file.read_text())
        expected_rel_path = str(dummy_model_dir.relative_to(config_file.parent) / "training")
        assert updated_config["nanosim_params"]["training_data_path"] == expected_rel_path

        # Check that both map-ont and map-hifi presets were built
        assert mock_create_idx.call_count == 2
        called_presets = [call.kwargs.get("preset") for call in mock_create_idx.call_args_list]
        assert "map-ont" in called_presets
        assert "map-hifi" in called_presets


# -----------------------------------------------------------------------------
# Problem 3: TarSlip (filter="data"), HTTPS URLs, and MD5 Checksum Verification
# -----------------------------------------------------------------------------


def test_download_nanosim_model_uses_tar_data_filter(tmp_path):
    """download_nanosim_model must extract tar archives with filter='data'."""
    mock_response = MagicMock()
    mock_response.iter_content.return_value = [b"chunk"]

    with (
        patch("download_references.requests") as mock_requests,
        patch("download_references.NANOSIM_MODEL_MD5", hashlib.md5(b"chunk").hexdigest()),
        patch("tarfile.open") as mock_tar_open,
    ):
        mock_requests.get.return_value = mock_response
        mock_tar = MagicMock()
        mock_tar_open.return_value.__enter__.return_value = mock_tar

        result = download_references.download_nanosim_model(output_dir=tmp_path, force=True)
        assert result is not None

        mock_tar.extractall.assert_called_once()
        assert mock_tar.extractall.call_args.kwargs.get("filter") == "data"


def test_download_nanosim_model_flattens_nested_directory(tmp_path):
    """download_nanosim_model flattens top-level model_name/ folder from upstream tarball."""
    import io

    model_name = "human_giab_hg002_sub1M_kitv14_dorado_v3.2.1"
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w:gz") as tar:
        data = b"profile_data"
        info = tarfile.TarInfo(name=f"{model_name}/training_model_profile")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    tar_bytes = tar_buf.getvalue()

    mock_response = MagicMock()
    mock_response.iter_content.return_value = [tar_bytes]

    with (
        patch("download_references.requests") as mock_requests,
        patch("download_references.NANOSIM_MODEL_MD5", hashlib.md5(tar_bytes).hexdigest()),
    ):
        mock_requests.get.return_value = mock_response
        extract_dir = download_references.download_nanosim_model(output_dir=tmp_path, force=True)
        assert extract_dir is not None
        assert (extract_dir / "training_model_profile").is_file()
        assert not (extract_dir / model_name).exists()


def test_install_references_extract_tar_uses_data_filter(tmp_path):
    """process_reference in install_references.py must extract tar files with filter='data'."""
    tar_file = tmp_path / "test_model.tar.gz"
    tar_file.write_bytes(b"dummy")

    ref_info = {
        "url": "https://example.com/test_model.tar.gz",
        "target_path": "test_model.tar.gz",
        "extract_tar": True,
        "extract_path": "extracted_dir",
    }

    with (
        patch("install_references.download_file"),
        patch("install_references.calculate_md5", return_value="abc123"),
        patch("tarfile.open") as mock_tar_open,
    ):
        mock_tar = MagicMock()
        mock_tar_open.return_value.__enter__.return_value = mock_tar

        install_references.process_reference(
            ref_name="test_ref",
            ref_info=ref_info,
            output_dir=tmp_path,
            skip_indexing=True,
            bwa_path="bwa",
            gatk_path="gatk",
            minimap2_path="minimap2",
            force=False,
        )

        mock_tar.extractall.assert_called_once()
        assert mock_tar.extractall.call_args.kwargs.get("filter") == "data"


def test_ucsc_das_url_uses_https():
    """get_ucsc_region_das must fetch sequences via https:// protocol."""
    with patch("download_references.requests") as mock_requests:
        mock_response = MagicMock()
        mock_response.content = b"<DASDNA><SEQUENCE><DNA>ACGT</DNA></SEQUENCE></DASDNA>"
        mock_requests.get.return_value = mock_response

        seq = download_references.get_ucsc_region_das("hg38", "chr1", 100, 200)
        assert seq == "ACGT"

        called_url = mock_requests.get.call_args[0][0]
        assert called_url.startswith("https://genome.ucsc.edu/")


def test_install_references_config_urls_and_checksums():
    """install_references_config.json must have https:// URLs and expected_md5 defined."""
    config_path = HELPERS_DIR / "install_references_config.json"
    assert config_path.exists()

    config = json.loads(config_path.read_text())
    refs = config.get("references", {})
    assert len(refs) > 0

    for name, info in refs.items():
        url = info.get("url", "")
        assert url.startswith("https://"), f"{name} URL is not HTTPS: {url}"
        assert "expected_md5" in info, f"{name} is missing expected_md5 checksum"
        assert len(info["expected_md5"]) == 32, f"{name} has invalid MD5 length"


def test_install_references_checksum_verification_success(tmp_path):
    """process_reference succeeds when MD5 checksum matches expected_md5."""
    target = tmp_path / "ref.txt"
    target.write_text("content")

    ref_info = {
        "url": "https://example.com/ref.txt",
        "target_path": "ref.txt",
        "expected_md5": "d41d8cd98f00b204e9800998ecf8427e",
    }

    with (
        patch("install_references.download_file"),
        patch("install_references.calculate_md5", return_value="D41D8CD98F00B204E9800998ECF8427E"),
    ):
        result = install_references.process_reference(
            ref_name="test_ref",
            ref_info=ref_info,
            output_dir=tmp_path,
            skip_indexing=True,
            bwa_path="bwa",
            gatk_path="gatk",
            minimap2_path="minimap2",
            force=False,
        )
        assert Path(result) == target.resolve()


def test_install_references_checksum_verification_failure(tmp_path):
    """process_reference exits with code 1 when MD5 checksum does not match expected_md5."""
    target = tmp_path / "ref.txt"
    target.write_text("content")

    ref_info = {
        "url": "https://example.com/ref.txt",
        "target_path": "ref.txt",
        "expected_md5": "expected_hash_value",
    }

    with (
        patch("install_references.download_file"),
        patch("install_references.calculate_md5", return_value="different_hash_value"),
        pytest.raises(SystemExit) as exc_info,
    ):
        install_references.process_reference(
            ref_name="test_ref",
            ref_info=ref_info,
            output_dir=tmp_path,
            skip_indexing=True,
            bwa_path="bwa",
            gatk_path="gatk",
            minimap2_path="minimap2",
            force=False,
        )
    assert exc_info.value.code == 1


# -----------------------------------------------------------------------------
# Problem 4: Safe Whitespace Splitting in execute_index_command
# -----------------------------------------------------------------------------


def test_execute_index_command_with_spaces_in_path(tmp_path):
    """execute_index_command must safely handle paths containing spaces without breaking."""
    space_dir = tmp_path / "dir with spaces"
    space_dir.mkdir()
    fasta_file = space_dir / "my reference.fa"
    fasta_file.write_text(">seq\nATCG\n")

    template = "bwa index {path}"

    with patch("subprocess.run") as mock_run:
        install_references.execute_index_command(template, fasta_file)

        mock_run.assert_called_once()
        cmd_args = mock_run.call_args[0][0]

        # bwa index must receive the path with spaces as a single argument element
        assert cmd_args == ["bwa", "index", str(fasta_file)]


def test_execute_index_command_failure_exits(tmp_path):
    """execute_index_command must exit with code 1 when subprocess fails."""
    import subprocess

    fasta_file = tmp_path / "ref.fa"
    fasta_file.write_text(">seq\nATCG\n")

    with (
        patch(
            "subprocess.run",
            side_effect=subprocess.CalledProcessError(1, ["cmd"], stderr=b"command error"),
        ),
        pytest.raises(SystemExit) as exc_info,
    ):
        install_references.execute_index_command("bwa index {path}", fasta_file)

    assert exc_info.value.code == 1


# -----------------------------------------------------------------------------
# Problem 2 (Install References): Preset-specific indexing in install_references.py
# -----------------------------------------------------------------------------


def test_install_references_preset_minimap2_indexing(tmp_path):
    """process_reference must build preset-specific indexes (map-ont, map-hifi)."""
    fasta_file = tmp_path / "reference.fna"
    fasta_file.write_text(">seq\nATCG\n")

    ref_info = {
        "url": "https://example.com/reference.fna",
        "target_path": "reference.fna",
        "minimap2_index_command": "minimap2 -x {preset} -d {path}.{preset}.mmi {path}",
    }

    with (
        patch("install_references.download_file"),
        patch("install_references.calculate_md5", return_value="hash"),
        patch("install_references.execute_index_command") as mock_exec,
    ):
        install_references.process_reference(
            ref_name="test_ref",
            ref_info=ref_info,
            output_dir=tmp_path,
            skip_indexing=False,
            bwa_path="bwa",
            gatk_path="gatk",
            minimap2_path="minimap2",
            force=False,
        )

        # Both map-ont and map-hifi indexing commands must have been executed
        assert mock_exec.call_count == 2
        commands_executed = [call.args[0] for call in mock_exec.call_args_list]

        assert any("-x map-ont -d {path}.map-ont.mmi {path}" in cmd for cmd in commands_executed)
        assert any("-x map-hifi -d {path}.map-hifi.mmi {path}" in cmd for cmd in commands_executed)


# -----------------------------------------------------------------------------
# Audit #7 / k / l / m: robust installs
# -----------------------------------------------------------------------------


def _nanosim_tar_bytes() -> bytes:
    import io

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        data = b"profile"
        info = tarfile.TarInfo(name="training_model_profile")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def test_download_nanosim_model_corrupt_download_leaves_no_folder(tmp_path):
    """A corrupt tarball fails MD5, leaves no model folder, and the next run re-downloads."""
    good = _nanosim_tar_bytes()
    corrupt = MagicMock()
    corrupt.iter_content.return_value = [b"not a tarball"]
    ok = MagicMock()
    ok.iter_content.return_value = [good]

    with (
        patch("download_references.requests") as mock_requests,
        patch("download_references.NANOSIM_MODEL_MD5", hashlib.md5(good).hexdigest()),
    ):
        mock_requests.get.side_effect = [corrupt, ok]
        assert download_references.download_nanosim_model(output_dir=tmp_path) is None
        model_dirs = [p for p in (tmp_path / "nanosim").iterdir() if p.is_dir()]
        assert model_dirs == []

        extract_dir = download_references.download_nanosim_model(output_dir=tmp_path)
        assert (extract_dir / "training_model_profile").is_file()
        assert mock_requests.get.call_count == 2
        assert mock_requests.get.call_args.kwargs["timeout"] == 60


def test_download_nanosim_model_failed_extract_removes_folder(tmp_path):
    """If extraction fails after MD5 passes, the partial model folder is removed."""
    data = b"garbage"
    resp = MagicMock()
    resp.iter_content.return_value = [data]
    with (
        patch("download_references.requests") as mock_requests,
        patch("download_references.NANOSIM_MODEL_MD5", hashlib.md5(data).hexdigest()),
    ):
        mock_requests.get.return_value = resp
        assert download_references.download_nanosim_model(output_dir=tmp_path) is None
    assert [p for p in (tmp_path / "nanosim").iterdir() if p.is_dir()] == []


def test_extract_gzip_failure_leaves_no_final_file(tmp_path):
    """An interrupted gzip extract must not leave the final (or .part) file behind."""
    gz = tmp_path / "ref.fna.gz"
    gz.write_bytes(gzip.compress(b">chr1\nACGT\n" * 1000))

    with (
        patch("install_references.shutil.copyfileobj", side_effect=OSError("disk full")),
        pytest.raises(SystemExit),
    ):
        install_references.extract_gzip(gz)
    assert not (tmp_path / "ref.fna").exists()
    assert not (tmp_path / "ref.fna.part").exists()

    assert install_references.extract_gzip(gz).read_bytes() == b">chr1\nACGT\n" * 1000


def test_process_reference_skips_download_when_md5_matches(tmp_path):
    """An existing file with matching MD5 is not downloaded again (unless --force)."""
    target = tmp_path / "ref.txt"
    target.write_text("content")
    ref_info = {
        "url": "https://example.com/ref.txt",
        "target_path": "ref.txt",
        "expected_md5": hashlib.md5(b"content").hexdigest(),
    }
    kwargs = {
        "ref_name": "r",
        "ref_info": ref_info,
        "output_dir": tmp_path,
        "skip_indexing": True,
        "bwa_path": "bwa",
        "gatk_path": "gatk",
        "minimap2_path": "minimap2",
    }

    with patch("install_references.download_file") as mock_dl:
        install_references.process_reference(force=False, **kwargs)
    mock_dl.assert_not_called()

    def fake_download(url, dest):
        dest.write_text("content")

    with patch("install_references.download_file", side_effect=fake_download) as mock_dl:
        install_references.process_reference(force=True, **kwargs)
    mock_dl.assert_called_once()


def test_process_reference_uses_htsjdk_dict_name(tmp_path):
    """An existing ref.dict (htsjdk naming) skips CreateSequenceDictionary."""
    (tmp_path / "ref.fna").write_text(">chr1\nACGT\n")
    (tmp_path / "ref.dict").write_text("@HD\n")
    ref_info = {
        "url": "https://example.com/ref.fna",
        "target_path": "ref.fna",
        "seq_dict_command": "gatk CreateSequenceDictionary -R {path}",
    }
    with (
        patch("install_references.download_file"),
        patch("install_references.execute_index_command") as mock_exec,
    ):
        install_references.process_reference(
            ref_name="r",
            ref_info=ref_info,
            output_dir=tmp_path,
            skip_indexing=True,
            bwa_path="bwa",
            gatk_path="gatk",
            minimap2_path="minimap2",
            force=False,
        )
    mock_exec.assert_not_called()


def test_create_tarball_records_relative_tarball_path(tmp_path, monkeypatch):
    """create_tarball must not crash on its own relative tarball path."""
    import generate_test_data

    monkeypatch.chdir(tmp_path)
    dataset = tmp_path / "work" / "testdata_40-70"
    dataset.mkdir(parents=True)
    (dataset / "a.txt").write_text("x")

    tarball = generate_test_data.create_tarball(dataset, "v9.9.9")

    info = json.loads(Path("data/testdata_releases/_latest/VERSION.json").read_text())
    assert info["source_tarball"] == str(tarball)
    assert info["generated_utc"].endswith("Z")
