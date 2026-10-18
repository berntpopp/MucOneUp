"""Tests for muc_one_up.cli.outputs module.

Tests cover:
- FASTA output writing (single and dual mutation mode)
- Structure file writing
- Mutated VNTR unit writing
- Haplotype-specific comment generation
- SNP integration in outputs
- Error handling
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from muc_one_up.cli.outputs import (
    _translate_snps_to_mutated,
    write_fasta_outputs,
    write_mutated_units,
    write_structure_files,
)
from muc_one_up.exceptions import FileOperationError
from muc_one_up.type_defs import HaplotypeResult, MutationTarget
from tests.utils import (
    assert_valid_fasta,
    count_fasta_records,
)


def _hr(seq: str, chain: list[str]) -> HaplotypeResult:
    """Helper to build HaplotypeResult from legacy-style arguments."""
    return HaplotypeResult.from_tuple((seq, chain))


@pytest.fixture
def mock_args():
    """Create mock args object for testing."""
    args = Mock()
    args.mutation_name = None
    args.mutation_targets = None
    args.output_structure = False
    args.random_snps = False
    args.snp_file = None
    args.snp_input_file = None  # Fix: properly set to None
    args.random_snp_density = None
    args.reference_assembly = "hg38"
    return args


@pytest.fixture
def sample_results(minimal_config: dict) -> list[HaplotypeResult]:
    """Sample simulation results as HaplotypeResult objects."""
    left = minimal_config["constants"]["hg38"]["left"]
    right = minimal_config["constants"]["hg38"]["right"]

    seq1 = left + "ATCGATCGATCG" + right
    chain1 = ["1", "2", "X", "B", "6", "7", "8", "9"]

    seq2 = left + "GCTAGCTAGCTA" + right
    chain2 = ["1", "2", "A", "B", "6p", "7", "8", "9"]

    return [_hr(seq1, chain1), _hr(seq2, chain2)]


@pytest.mark.unit
class TestWriteFastaOutputs:
    """Tests for write_fasta_outputs function."""

    def test_write_fasta_single_mode_no_mutation(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA in single mode without mutation."""
        mock_args.output_structure = False

        out_results, mutated_results, snp_info_normal, snp_info_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        # Check output file exists
        output_file = tmp_path / "test.001.simulated.fa"
        assert output_file.exists()
        assert_valid_fasta(output_file)
        assert count_fasta_records(output_file) == 2

        # Verify results are returned
        assert len(out_results) == len(sample_results)
        assert mutated_results is None
        assert snp_info_normal == {}  # No SNPs applied
        assert snp_info_mut == {}

    def test_write_fasta_single_mode_with_mutation(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA in single mode with mutation."""
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,25", "2,30"]

        _out_results, _mutated_results, _snp_info_normal, _snp_info_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.simulated.fa"
        assert output_file.exists()

        # Check that mutation comments are in headers
        with output_file.open() as f:
            content = f.read()

        assert "Mutation Applied: dupC" in content

    def test_write_fasta_single_mode_with_random_mutation_targets(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Random mutation targets (from mutation_positions) are written to FASTA comments."""
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = None

        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=[MutationTarget(1, 2)],
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.simulated.fa"
        assert output_file.exists()
        content = output_file.read_text()
        assert "Mutation Applied: dupC (Targets: [(1, 2)])" in content

    def test_write_fasta_dual_mutation_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA in dual mutation mode."""
        mutated_results = [
            _hr(sample_results[0].sequence, ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            _hr(sample_results[1].sequence, ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
        ]

        _out_results, _mut_results, _snp_info_normal, _snp_info_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[(1, 25), (2, 30)],
            structure_mutation_info=None,
        )

        # Check both output files exist
        normal_file = tmp_path / "test.001.normal.simulated.fa"
        mut_file = tmp_path / "test.001.mut.simulated.fa"

        assert normal_file.exists()
        assert mut_file.exists()
        assert_valid_fasta(normal_file)
        assert_valid_fasta(mut_file)

        # Check content of normal file
        with normal_file.open() as f:
            content = f.read()
        assert "Normal sequence (no mutations applied)" in content

        # Check content of mutated file
        with mut_file.open() as f:
            content = f.read()
        assert "Mutation Applied: dupC" in content
        assert "[(1, 25), (2, 30)]" in content

    def test_write_fasta_with_structure_mutation_info(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA with mutation info from structure file."""
        structure_mutation_info = {
            "name": "dupC",
            "targets": [(1, 25)],
        }

        _out_results, _, _, _ = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=structure_mutation_info,
        )

        output_file = tmp_path / "test.001.simulated.fa"
        with output_file.open() as f:
            content = f.read()

        assert "Mutation Applied: dupC" in content

    def test_write_fasta_haplotype_specific_comments(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test that mutation comments are haplotype-specific."""
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,25"]  # Only haplotype 1

        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.simulated.fa"
        with output_file.open() as f:
            lines = [line.strip() for line in f if line.startswith(">")]

        # First haplotype should have mutation comment
        assert "Mutation Applied" in lines[0]
        # Second haplotype should not
        assert len(lines[1].split()) == 1  # Only header, no comment

    def test_write_fasta_tuple_mutation_targets(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test handling mutation targets as tuples instead of strings."""
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = [(1, 25), (2, 30)]  # Tuples not strings

        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.simulated.fa"
        assert output_file.exists()

    def test_write_fasta_with_nested_directory(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA to nested directory (must exist)."""
        out_dir = tmp_path / "subdir" / "nested"
        out_dir.mkdir(parents=True, exist_ok=True)  # Create directory first

        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(out_dir),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = out_dir / "test.001.simulated.fa"
        assert output_file.exists()

    def test_write_fasta_error_handling_dual_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test error handling in dual mutation mode."""
        # Use invalid path to trigger error
        with pytest.raises(FileOperationError, match="Writing dual FASTA outputs failed"):
            write_fasta_outputs(
                args=mock_args,
                config=minimal_config,
                out_dir="/nonexistent/path/that/cannot/exist",
                out_base="test",
                sim_index=1,
                results=sample_results,
                mutated_results=sample_results,
                dual_mutation_mode=True,
                mutation_pair=["normal", "dupC"],
                mutation_positions=[(1, 25)],
                structure_mutation_info=None,
            )

    def test_write_fasta_error_handling_single_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test error handling in single mode."""
        with pytest.raises(FileOperationError, match="Writing FASTA output failed"):
            write_fasta_outputs(
                args=mock_args,
                config=minimal_config,
                out_dir="/nonexistent/path/that/cannot/exist",
                out_base="test",
                sim_index=1,
                results=sample_results,
                mutated_results=None,
                dual_mutation_mode=False,
                mutation_pair=None,
                mutation_positions=None,
                structure_mutation_info=None,
            )


@pytest.mark.unit
class TestWriteMutatedUnits:
    """Tests for write_mutated_units function."""

    def test_write_mutated_units_single_mode(self, tmp_path: Path, mock_args, minimal_config: dict):
        """Test writing mutated VNTR units in single mode."""
        mock_args.mutation_name = "dupC"

        mutated_units = {
            1: [(25, "ATCGATCGATCGDUP")],
            2: [(30, "GCTAGCTAGCTADUP")],
        }

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=mutated_units,
            dual_mutation_mode=False,
        )

        output_file = tmp_path / "test.001.mutated_unit.fa"
        assert output_file.exists()

        # Check content
        with output_file.open() as f:
            content = f.read()

        assert ">haplotype_1_repeat_25" in content
        assert "ATCGATCGATCGDUP" in content
        assert ">haplotype_2_repeat_30" in content
        assert "GCTAGCTAGCTADUP" in content

    def test_write_mutated_units_dual_mode(self, tmp_path: Path, mock_args, minimal_config: dict):
        """Test writing mutated VNTR units in dual mutation mode."""
        mock_args.mutation_name = "dupC"

        mutated_units = {
            1: [(25, "ATCGATCGATCGDUP")],
        }

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=mutated_units,
            dual_mutation_mode=True,
        )

        # In dual mode, should have .mut variant suffix
        output_file = tmp_path / "test.001.mut.mutated_unit.fa"
        assert output_file.exists()

    def test_write_mutated_units_no_mutation(self, tmp_path: Path, mock_args, minimal_config: dict):
        """Test that no file is written when no mutation is applied."""
        mock_args.mutation_name = None

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=None,
            dual_mutation_mode=False,
        )

        # No output files should be created
        assert not list(tmp_path.glob("*.mutated_unit.fa"))

    def test_write_mutated_units_multiple_per_haplotype(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """Test writing multiple mutated units for same haplotype."""
        mock_args.mutation_name = "dupC"

        mutated_units = {
            1: [(10, "SEQ1"), (20, "SEQ2"), (30, "SEQ3")],
        }

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=mutated_units,
            dual_mutation_mode=False,
        )

        output_file = tmp_path / "test.001.mutated_unit.fa"
        with output_file.open() as f:
            content = f.read()

        assert content.count(">haplotype_1_repeat_") == 3
        assert "SEQ1" in content
        assert "SEQ2" in content
        assert "SEQ3" in content

    def test_write_mutated_units_error_handling(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """Test error handling when writing mutated units."""
        mock_args.mutation_name = "dupC"
        mutated_units = {1: [(25, "ATCG")]}

        with pytest.raises(FileOperationError, match="Writing mutated VNTR unit FASTA failed"):
            write_mutated_units(
                args=mock_args,
                out_dir="/nonexistent/path",
                out_base="test",
                sim_index=1,
                mutated_units=mutated_units,
                dual_mutation_mode=False,
            )


@pytest.mark.unit
class TestWriteStructureFiles:
    """Tests for write_structure_files function."""

    def test_write_structure_file_disabled(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test that no structure file is written when output_structure is False."""
        mock_args.output_structure = False

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        # No structure files should be created
        assert not list(tmp_path.glob("*.vntr_structure.txt"))

    def test_write_structure_file_single_mode_no_mutation(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing structure file in single mode without mutation."""
        mock_args.output_structure = True

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        assert output_file.exists()

        with output_file.open() as f:
            lines = [line.strip() for line in f]

        # Should have 2 haplotype lines
        assert "haplotype_1\t1-2-X-B-6-7-8-9" in lines
        assert "haplotype_2\t1-2-A-B-6p-7-8-9" in lines

    def test_write_structure_file_with_mutation(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing structure file with mutation."""
        mock_args.output_structure = True
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,25", "2,30"]

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        with output_file.open() as f:
            content = f.read()

        assert "# Mutation Applied: dupC" in content
        assert "Targets:" in content

    def test_write_structure_file_with_random_target(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing structure file with random mutation target."""
        mock_args.output_structure = True
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = None  # Random target

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        with output_file.open() as f:
            content = f.read()

        assert "# Mutation Applied: dupC (Target: random)" in content

    def test_write_structure_file_with_random_target_resolved_positions(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """When random targets are resolved into mutation_positions, write specific targets."""
        mock_args.output_structure = True
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = None

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=[MutationTarget(1, 2)],
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        content = output_file.read_text()
        assert "# Mutation Applied: dupC (Targets: [(1, 2)])" in content

    def test_write_structure_file_with_structure_mutation_info(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing structure file with mutation info from structure file."""
        mock_args.output_structure = True

        structure_mutation_info = {
            "name": "dupC",
            "targets": [(1, 25), (2, 30)],
        }

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=structure_mutation_info,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        with output_file.open() as f:
            content = f.read()

        assert "# Mutation Applied: dupC" in content
        assert "[(1, 25), (2, 30)]" in content

    def test_write_structure_file_dual_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing structure files in dual mutation mode."""
        mock_args.output_structure = True

        mutated_results = [
            _hr(sample_results[0].sequence, ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            _hr(sample_results[1].sequence, ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
        ]

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[(1, 25), (2, 30)],
            structure_mutation_info=None,
        )

        # Check both files exist
        normal_file = tmp_path / "test.001.normal.vntr_structure.txt"
        mut_file = tmp_path / "test.001.mut.vntr_structure.txt"

        assert normal_file.exists()
        assert mut_file.exists()

        # Check normal file content
        with normal_file.open() as f:
            content = f.read()
        assert "# Normal sequence (no mutations applied)" in content
        assert "haplotype_1\t1-2-X-B-6-7-8-9" in content

        # Check mutated file content
        with mut_file.open() as f:
            content = f.read()
        assert "# Mutation Applied: dupC" in content
        assert "haplotype_1\t1-2-Xm-B-6-7-8-9" in content
        assert "haplotype_2\t1-2-Am-B-6p-7-8-9" in content

    def test_write_structure_file_preserves_mutation_markers(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """Test that mutation markers (m suffix) are preserved in structure file."""
        mock_args.output_structure = True

        results = [
            _hr("SEQUENCE1", ["1", "2m", "X", "B", "6", "7", "8", "9"]),
            _hr("SEQUENCE2", ["1", "2", "Xm", "B", "6p", "7", "8", "9"]),
        ]

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        output_file = tmp_path / "test.001.vntr_structure.txt"
        with output_file.open() as f:
            content = f.read()

        assert "1-2m-X-B-6-7-8-9" in content
        assert "1-2-Xm-B-6p-7-8-9" in content

    def test_write_structure_file_error_handling_dual_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test error handling in dual mutation mode."""
        mock_args.output_structure = True

        with pytest.raises(FileOperationError, match="Writing dual structure files failed"):
            write_structure_files(
                args=mock_args,
                out_dir="/nonexistent/path",
                out_base="test",
                sim_index=1,
                results=sample_results,
                mutated_results=sample_results,
                dual_mutation_mode=True,
                mutation_pair=["normal", "dupC"],
                mutation_positions=[(1, 25)],
                structure_mutation_info=None,
            )

    def test_write_structure_file_error_handling_single_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test error handling in single mode."""
        mock_args.output_structure = True

        with pytest.raises(FileOperationError, match="Writing structure file failed"):
            write_structure_files(
                args=mock_args,
                out_dir="/nonexistent/path",
                out_base="test",
                sim_index=1,
                results=sample_results,
                mutated_results=None,
                dual_mutation_mode=False,
                mutation_pair=None,
                mutation_positions=None,
                structure_mutation_info=None,
            )


@pytest.mark.integration
class TestOutputsIntegration:
    """Integration tests for cli.outputs module."""

    def test_full_output_workflow_single_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test complete output workflow in single mode."""
        mock_args.output_structure = True
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,25"]

        mutated_units = {1: [(25, "MUTATEDSEQ")]}

        # Write all outputs
        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=mutated_units,
            dual_mutation_mode=False,
        )

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=None,
            structure_mutation_info=None,
        )

        # Verify all expected files exist
        assert (tmp_path / "test.001.simulated.fa").exists()
        assert (tmp_path / "test.001.mutated_unit.fa").exists()
        assert (tmp_path / "test.001.vntr_structure.txt").exists()

    def test_full_output_workflow_dual_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test complete output workflow in dual mutation mode."""
        mock_args.output_structure = True
        mock_args.mutation_name = "dupC"  # Required for write_mutated_units

        mutated_results = [
            _hr(sample_results[0].sequence, ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            _hr(sample_results[1].sequence, ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
        ]

        mutated_units = {1: [(25, "MUTATEDSEQ")]}

        # Write all outputs
        write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[(1, 25)],
            structure_mutation_info=None,
        )

        write_mutated_units(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            mutated_units=mutated_units,
            dual_mutation_mode=True,
        )

        write_structure_files(
            args=mock_args,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[(1, 25)],
            structure_mutation_info=None,
        )

        # Verify all expected files exist with correct naming
        assert (tmp_path / "test.001.normal.simulated.fa").exists()
        assert (tmp_path / "test.001.mut.simulated.fa").exists()
        assert (tmp_path / "test.001.mut.mutated_unit.fa").exists()
        assert (tmp_path / "test.001.normal.vntr_structure.txt").exists()
        assert (tmp_path / "test.001.mut.vntr_structure.txt").exists()

    def test_dual_mode_random_snps_applied_and_shifted_once(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Dual mode applies random SNPs to normal, translates them to mutant with indel shift,
        and writes random_snp_output_file only once.
        """
        import random

        mock_args.output_structure = False
        mock_args.random_snps = True
        mock_args.random_snp_density = 50.0
        snp_out = tmp_path / "random_snps.tsv"
        mock_args.random_snp_output_file = str(snp_out)
        mock_args.random_snp_region = "all"
        mock_args.random_snp_haplotypes = "all"
        mock_args.seed = 42

        # Create mutated sequence with 1-bp insertion (e.g. dupC at position after left flank)
        left = minimal_config["constants"]["hg38"]["left"]
        normal_seq = sample_results[0].sequence
        # Insert a nucleotide at position len(left) + 2
        mut_pos = len(left) + 2
        mut_seq = normal_seq[:mut_pos] + "C" + normal_seq[mut_pos:]

        mutated_results = [
            HaplotypeResult(sequence=mut_seq, chain=sample_results[0].chain),
            HaplotypeResult(sequence=sample_results[1].sequence, chain=sample_results[1].chain),
        ]

        rng = random.Random(42)
        _out_res, _out_mut_res, snp_normal, snp_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test",
            sim_index=1,
            results=sample_results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[MutationTarget(1, 1)],
            structure_mutation_info=None,
            rng=rng,
        )

        assert snp_out.exists()
        # Ensure random_snps.tsv has entries
        lines = [
            line.strip()
            for line in snp_out.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
        assert len(lines) > 0

        # Normal and mutant results have SNPs
        assert len(snp_normal.get(0, [])) > 0
        assert len(snp_mut.get(0, [])) > 0

        # For any SNP downstream of mut_pos, mutant coordinate must be shifted by +1
        normal_snps_hap0 = snp_normal[0]
        mut_snps_hap0 = snp_mut[0]
        assert len(normal_snps_hap0) == len(mut_snps_hap0)
        for n_snp, m_snp in zip(normal_snps_hap0, mut_snps_hap0, strict=True):
            if n_snp["position"] >= mut_pos:
                assert m_snp["position"] == n_snp["position"] + 1

    def test_dual_mode_snp_inside_mutated_repeat_upstream_of_insertion(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """(a) SNP inside mutated repeat upstream of mid-repeat insertion remains unshifted."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]
        rep3 = minimal_config["repeats"]["X"]

        normal_seq = left + norm_rep1 + rep2 + rep3 + right
        results = [_hr(normal_seq, ["1", "2", "X"])]

        # Mutate repeat 1 with insertion at offset 6
        mut_rep1 = norm_rep1[:6] + "C" + norm_rep1[6:]
        mut_seq = left + mut_rep1 + rep2 + rep3 + right
        mut_hr = _hr(mut_seq, ["1m", "2", "X"])
        mutated_results = [mut_hr]

        # SNP at offset 2 (upstream of insertion at 6)
        snp_pos = len(left) + 2
        ref_base = normal_seq[snp_pos]
        alt_base = "T" if ref_base != "T" else "A"

        snp_tsv = tmp_path / "snps_a.tsv"
        snp_tsv.write_text(f"1\t{snp_pos}\t{ref_base}\t{alt_base}\n")
        mock_args.snp_input_file = str(snp_tsv)
        mock_args.random_snps = False

        _, _, snp_normal, snp_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test_a",
            sim_index=1,
            results=results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[MutationTarget(1, 1)],
            structure_mutation_info=None,
        )

        assert len(snp_normal[0]) == 1
        assert len(snp_mut[0]) == 1
        assert snp_mut[0][0]["position"] == snp_pos
        assert snp_mut[0][0]["ref_base"] == ref_base

    def test_dual_mode_snp_between_two_repeat_insertions(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """(b) Insertions in repeat 1 and repeat 3; SNP in repeat 2 is shifted only by repeat 1's insertion."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]
        norm_rep3 = minimal_config["repeats"]["X"]

        normal_seq = left + norm_rep1 + rep2 + norm_rep3 + right
        results = [_hr(normal_seq, ["1", "2", "X"])]

        # Repeat 1 has +1 insertion, Repeat 3 has +2 insertion
        mut_rep1 = norm_rep1[:4] + "A" + norm_rep1[4:]
        mut_rep2 = rep2
        mut_rep3 = norm_rep3[:2] + "CC" + norm_rep3[2:]
        mut_seq = left + mut_rep1 + mut_rep2 + mut_rep3 + right
        mut_hr = _hr(mut_seq, ["1m", "2", "Xm"])
        mut_hr._repeat_lengths = [len(mut_rep1), len(mut_rep2), len(mut_rep3)]
        mutated_results = [mut_hr]

        # SNP inside repeat 2 at offset 3
        snp_pos = len(left) + len(norm_rep1) + 3
        ref_base = normal_seq[snp_pos]
        alt_base = "G" if ref_base != "G" else "T"

        snp_tsv = tmp_path / "snps_b.tsv"
        snp_tsv.write_text(f"1\t{snp_pos}\t{ref_base}\t{alt_base}\n")
        mock_args.snp_input_file = str(snp_tsv)
        mock_args.random_snps = False

        _, _, snp_normal, snp_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test_b",
            sim_index=1,
            results=results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[MutationTarget(1, 1), MutationTarget(1, 3)],
            structure_mutation_info=None,
        )

        assert len(snp_normal[0]) == 1
        assert len(snp_mut[0]) == 1
        # Shifted by +1 from repeat 1's insertion; unaffected by repeat 3's +2 insertion
        assert snp_mut[0][0]["position"] == snp_pos + 1
        assert snp_mut[0][0]["ref_base"] == ref_base

    def test_dual_mode_snp_inside_deleted_span_omitted(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """(c) SNP falling inside a deleted span in a mutated repeat is omitted from mutant SNPs."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]
        rep3 = minimal_config["repeats"]["X"]

        normal_seq = left + norm_rep1 + rep2 + rep3 + right
        results = [_hr(normal_seq, ["1", "2", "X"])]

        # Repeat 1 has bases 2..4 deleted (-2 bp)
        mut_rep1 = norm_rep1[:2] + norm_rep1[4:]
        mut_seq = left + mut_rep1 + rep2 + rep3 + right
        mut_hr = _hr(mut_seq, ["1m", "2", "X"])
        mutated_results = [mut_hr]

        # SNP inside the deleted span (offset 3 in repeat 1)
        snp_pos = len(left) + 3
        ref_base = normal_seq[snp_pos]
        alt_base = "A" if ref_base != "A" else "C"

        snp_tsv = tmp_path / "snps_c.tsv"
        snp_tsv.write_text(f"1\t{snp_pos}\t{ref_base}\t{alt_base}\n")
        mock_args.snp_input_file = str(snp_tsv)
        mock_args.random_snps = False

        _, _, snp_normal, snp_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test_c",
            sim_index=1,
            results=results,
            mutated_results=mutated_results,
            dual_mutation_mode=True,
            mutation_pair=["normal", "dupC"],
            mutation_positions=[MutationTarget(1, 1)],
            structure_mutation_info=None,
        )

        assert len(snp_normal[0]) == 1
        # Mutated haplotype should omit this deleted SNP
        assert len(snp_mut.get(0, [])) == 0

    def test_translate_snps_to_mutated_fallback_multi_repeat_indels(
        self, mock_args, minimal_config: dict
    ):
        """Test fallback repeat-length inference when _repeat_lengths is None with multi-repeat indels."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]
        norm_rep3 = minimal_config["repeats"]["X"]

        normal_seq = left + norm_rep1 + rep2 + norm_rep3 + right
        normal_hr = _hr(normal_seq, ["1", "2", "X"])

        # Repeat 1 has +1 insertion, Repeat 3 has +1 insertion (total delta = +2)
        mut_rep1 = norm_rep1[:4] + "A" + norm_rep1[4:]
        mut_rep2 = rep2
        mut_rep3 = norm_rep3[:2] + "C" + norm_rep3[2:]
        mut_seq = left + mut_rep1 + mut_rep2 + mut_rep3 + right

        # Notice: _repeat_lengths is intentionally NOT set (None) to exercise fallback
        mut_hr = _hr(mut_seq, ["1m", "2", "Xm"])
        assert mut_hr._repeat_lengths is None

        # SNP inside repeat 2 at offset 3
        snp_pos = len(left) + len(norm_rep1) + 3
        ref_base = normal_seq[snp_pos]
        applied_snps = {
            0: [{"haplotype": 1, "position": snp_pos, "ref_base": ref_base, "alt_base": "T"}]
        }

        translated = _translate_snps_to_mutated(
            args=mock_args,
            config=minimal_config,
            normal_results=[normal_hr],
            mutated_results=[mut_hr],
            applied_snp_info_normal=applied_snps,
            mutation_positions=[MutationTarget(1, 1), MutationTarget(1, 3)],
        )

        assert len(translated) == 1
        # Shifted by +1 from repeat 1's insertion
        assert translated[0]["position"] == snp_pos + 1
        assert translated[0]["ref_base"] == ref_base

    def test_single_mode_snp_coordinate_translation_across_insertion(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """Test that single-mutation mode translates --snp-input-file coordinates across indels."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]
        norm_rep3 = minimal_config["repeats"]["X"]

        normal_seq = left + norm_rep1 + rep2 + norm_rep3 + right
        normal_hr = _hr(normal_seq, ["1", "2", "X"])

        # Mutate repeat 1 with +1 insertion
        mut_rep1 = norm_rep1[:4] + "A" + norm_rep1[4:]
        mut_seq = left + mut_rep1 + rep2 + norm_rep3 + right
        mut_hr = _hr(mut_seq, ["1m", "2", "X"])
        mut_hr._repeat_lengths = [len(mut_rep1), len(rep2), len(norm_rep3)]

        # SNP in repeat 2 (downstream of insertion)
        snp_pos = len(left) + len(norm_rep1) + 2
        ref_base = normal_seq[snp_pos]
        alt_base = "G" if ref_base != "G" else "T"

        snp_file = tmp_path / "snps_single.tsv"
        snp_file.write_text(f"1\t{snp_pos}\t{ref_base}\t{alt_base}\n")

        mock_args.snp_input_file = str(snp_file)
        mock_args.random_snps = False
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,1"]

        out_results, mutated_res, applied_normal, _applied_mut = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test_single_trans",
            sim_index=1,
            results=[mut_hr],
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=[MutationTarget(1, 1)],
            structure_mutation_info=None,
            normal_results=[normal_hr],
        )

        assert mutated_res is None
        assert 0 in applied_normal
        assert len(applied_normal[0]) == 1
        # Position should be translated from snp_pos to snp_pos + 1
        assert applied_normal[0][0]["position"] == snp_pos + 1

        # Check sequence has the mutation applied at snp_pos + 1
        final_seq = out_results[0].sequence
        assert final_seq[snp_pos + 1] == alt_base

    def test_single_mode_snp_inside_deleted_span_omitted(
        self, tmp_path: Path, mock_args, minimal_config: dict
    ):
        """Test that single-mutation mode omits SNPs falling inside a deleted span."""
        left = minimal_config["constants"]["hg38"]["left"]
        right = minimal_config["constants"]["hg38"]["right"]
        norm_rep1 = minimal_config["repeats"]["1"]
        rep2 = minimal_config["repeats"]["2"]

        normal_seq = left + norm_rep1 + rep2 + right
        normal_hr = _hr(normal_seq, ["1", "2"])

        # Repeat 1 has bases 2..4 deleted
        mut_rep1 = norm_rep1[:2] + norm_rep1[4:]
        mut_seq = left + mut_rep1 + rep2 + right
        mut_hr = _hr(mut_seq, ["1m", "2"])

        # SNP inside deleted span (offset 3)
        snp_pos = len(left) + 3
        ref_base = normal_seq[snp_pos]
        alt_base = "A" if ref_base != "A" else "C"

        snp_file = tmp_path / "snps_del.tsv"
        snp_file.write_text(f"1\t{snp_pos}\t{ref_base}\t{alt_base}\n")

        mock_args.snp_input_file = str(snp_file)
        mock_args.random_snps = False
        mock_args.mutation_name = "del"
        mock_args.mutation_targets = ["1,1"]

        _out_results, _, applied_normal, _ = write_fasta_outputs(
            args=mock_args,
            config=minimal_config,
            out_dir=str(tmp_path),
            out_base="test_single_del",
            sim_index=1,
            results=[mut_hr],
            mutated_results=None,
            dual_mutation_mode=False,
            mutation_pair=None,
            mutation_positions=[MutationTarget(1, 1)],
            structure_mutation_info=None,
            normal_results=[normal_hr],
        )

        # SNP inside deleted span is omitted
        assert len(applied_normal.get(0, [])) == 0
