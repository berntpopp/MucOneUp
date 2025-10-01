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
    write_fasta_outputs,
    write_mutated_units,
    write_structure_files,
)
from muc_one_up.exceptions import FileOperationError
from tests.utils import (
    assert_valid_fasta,
    count_fasta_records,
)


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
    return args


@pytest.fixture
def sample_results(minimal_config: dict) -> list[tuple[str, list[str]]]:
    """Sample simulation results (sequence, chain) tuples."""
    left = minimal_config["constants"]["hg38"]["left"]
    right = minimal_config["constants"]["hg38"]["right"]

    seq1 = left + "ATCGATCGATCG" + right
    chain1 = ["1", "2", "X", "B", "6", "7", "8", "9"]

    seq2 = left + "GCTAGCTAGCTA" + right
    chain2 = ["1", "2", "A", "B", "6p", "7", "8", "9"]

    return [(seq1, chain1), (seq2, chain2)]


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
        assert out_results == sample_results
        assert mutated_results is None
        assert snp_info_normal == {}  # No SNPs applied
        assert snp_info_mut == {}

    def test_write_fasta_single_mode_with_mutation(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA in single mode with mutation."""
        mock_args.mutation_name = "dupC"
        mock_args.mutation_targets = ["1,25", "2,30"]

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

        output_file = tmp_path / "test.001.simulated.fa"
        assert output_file.exists()

        # Check that mutation comments are in headers
        with output_file.open() as f:
            content = f.read()

        assert "Mutation Applied: dupC" in content

    def test_write_fasta_dual_mutation_mode(
        self, tmp_path: Path, mock_args, minimal_config: dict, sample_results: list
    ):
        """Test writing FASTA in dual mutation mode."""
        mutated_results = [
            (sample_results[0][0], ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            (sample_results[1][0], ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
        ]

        out_results, mut_results, snp_info_normal, snp_info_mut = write_fasta_outputs(
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

        out_results, _, _, _ = write_fasta_outputs(
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
            (sample_results[0][0], ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            (sample_results[1][0], ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
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
            ("SEQUENCE1", ["1", "2m", "X", "B", "6", "7", "8", "9"]),
            ("SEQUENCE2", ["1", "2", "Xm", "B", "6p", "7", "8", "9"]),
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
            (sample_results[0][0], ["1", "2", "Xm", "B", "6", "7", "8", "9"]),
            (sample_results[1][0], ["1", "2", "Am", "B", "6p", "7", "8", "9"]),
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
