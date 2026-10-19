"""Unit tests for the Illumina read simulation pipeline."""

from __future__ import annotations

from pathlib import Path

from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline
from muc_one_up.read_simulator.stages import AlignmentResult, FragmentResult


def _minimal_illumina_config(tmp_path: Path) -> dict:
    ref_file = tmp_path / "hg38.fa"
    ref_file.write_text(">chr1\nACGT\n")
    for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
        (tmp_path / f"hg38.fa{ext}").touch()

    return {
        "tools": {
            "reseq": "reseq",
            "faToTwoBit": "faToTwoBit",
            "samtools": "samtools",
            "pblat": "pblat",
            "bwa": "bwa",
        },
        "read_simulation": {
            "reseq_model": str(tmp_path / "model"),
            "sample_bam": str(tmp_path / "sample.bam"),
            "output_dir": str(tmp_path),
            "vntr_capture_efficiency": {"enabled": False},
        },
        "reference_assembly": "hg38",
        "reference_genomes": {
            "hg38": {
                "fasta_path": str(ref_file),
                "vntr_region": "chr1:155188487-155192239",
            }
        },
    }


def test_illumina_pipeline_respects_skip_alignment(tmp_path, mocker):
    """When skip_alignment=True, alignment and reference resolution are skipped."""
    config = _minimal_illumina_config(tmp_path)
    config["read_simulation"]["skip_alignment"] = True

    input_fa = tmp_path / "input.fa"
    input_fa.write_text(">test\nACGTACGT\n")

    r1_path = str(tmp_path / "out_R1.fastq.gz")
    r2_path = str(tmp_path / "out_R2.fastq.gz")

    mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
    mocker.patch("muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={})
    mocker.patch("muc_one_up.read_simulator.pipeline.log_tool_versions")
    mocker.patch("muc_one_up.read_simulator.pipeline.cleanup_intermediates")
    mock_meta = mocker.patch("muc_one_up.read_simulator.pipeline.create_pipeline_metadata")
    mock_manifest = mocker.patch("muc_one_up.read_simulator.pipeline.generate_read_manifest")
    mocker.patch(
        "muc_one_up.read_simulator.pipeline.prepare_fragments",
        return_value=FragmentResult(r1_fastq=r1_path, r2_fastq=r2_path),
    )
    mock_resolve_ref = mocker.patch("muc_one_up.read_simulator.pipeline.resolve_human_reference")
    mock_align = mocker.patch("muc_one_up.read_simulator.pipeline.align_and_refine")

    result = simulate_reads_pipeline(config, str(input_fa))

    assert result == r1_path
    mock_resolve_ref.assert_not_called()
    mock_align.assert_not_called()
    mock_manifest.assert_not_called()
    assert mock_meta.called
    assert "bwa" not in mock_meta.call_args[1]["tools_used"]


def test_illumina_pipeline_aligns_by_default(tmp_path, mocker):
    """When skip_alignment is False/omitted, alignment is performed and BAM is returned."""
    config = _minimal_illumina_config(tmp_path)

    input_fa = tmp_path / "input.fa"
    input_fa.write_text(">test\nACGTACGT\n")

    r1_path = str(tmp_path / "out_R1.fastq.gz")
    r2_path = str(tmp_path / "out_R2.fastq.gz")
    final_bam = str(tmp_path / "out.bam")

    mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
    mocker.patch("muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={})
    mocker.patch("muc_one_up.read_simulator.pipeline.log_tool_versions")
    mocker.patch("muc_one_up.read_simulator.pipeline.cleanup_intermediates")
    mock_meta = mocker.patch("muc_one_up.read_simulator.pipeline.create_pipeline_metadata")
    mock_manifest = mocker.patch("muc_one_up.read_simulator.pipeline.generate_read_manifest")
    mocker.patch(
        "muc_one_up.read_simulator.pipeline.prepare_fragments",
        return_value=FragmentResult(r1_fastq=r1_path, r2_fastq=r2_path),
    )
    mock_resolve_ref = mocker.patch(
        "muc_one_up.read_simulator.pipeline.resolve_human_reference", return_value="/ref/hg38.fa"
    )
    mock_align = mocker.patch(
        "muc_one_up.read_simulator.pipeline.align_and_refine",
        return_value=AlignmentResult(final_bam=final_bam),
    )

    result = simulate_reads_pipeline(config, str(input_fa))

    assert result == final_bam
    mock_resolve_ref.assert_called_once()
    mock_align.assert_called_once()
    mock_manifest.assert_called_once()
    assert "bwa" in mock_meta.call_args[1]["tools_used"]
