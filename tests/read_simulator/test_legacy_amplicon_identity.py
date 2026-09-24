"""Legacy amplicon byte-identity guard (0.45 compatibility promise, #132 review).

Simulations without a read profile must keep producing the same reads. With
pbsim3/ccs/samtools mocked at the ``run_command`` level, a seeded legacy run
must hand pbsim3 the same template FASTAs (pinned md5) with the same argument
lists. Any change to PCR-bias sampling, template writing or tool arguments of
the legacy path fails here, without needing the external tools.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli

FWD = "GGAGAAAAGGAGACTTCGGCTACCCAG"
REV_RC = "TCAGCTTCTACTCTGGTGCACAACGGC"  # reverse complement of the config's reverse primer
UNIT = "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA"
HAPLOTYPES = {
    "haplotype_1": "ACGTTGCA" * 10 + FWD + UNIT * 20 + REV_RC + "TTGCAACG" * 10,
    "haplotype_2": "ACGTTGCA" * 10 + FWD + UNIT * 32 + REV_RC + "TTGCAACG" * 10,
}
SEED = 42
COVERAGE = 100
FAKE_BAM = b"BAM\x01" + b"\x00" * 4096  # above MIN_VALID_CCS_OUTPUT_SIZE

WRAPPER_MODULES = (
    "muc_one_up.read_simulator.wrappers.pbsim3_wrapper",
    "muc_one_up.read_simulator.wrappers.ccs_wrapper",
    "muc_one_up.read_simulator.wrappers.samtools_convert",
    "muc_one_up.read_simulator.wrappers.samtools_core",
)

# Pinned for the seeded legacy runs below (identical on origin/main eabb721, i.e.
# MucOneUp 0.45.0). The template md5s depend on amplicon extraction, PCR-bias
# sampling (seed 42) and template writing; the calls on config.json defaults.
# Re-pin only for an intentional, documented legacy output change.
EXPECTED: dict[str, dict[str, Any]] = {
    "ont": {
        "templates": {
            "template_hap1.fa": "876285d754365e28ddb8090a393e9303",
            "template_hap2.fa": "0b8cc660c00c40eb9ba01df8d868c19f",
        },
        "calls": [
            [
                "pbsim",
                "--strategy",
                "templ",
                "--method",
                "qshmm",
                "--qshmm",
                "<dir>/QSHMM-ONT-HQ.model",
                "--template",
                "<dir>/template_hap1.fa",
                "--pass-num",
                "1",
                "--accuracy-mean",
                "0.95",
                "--prefix",
                "<dir>/ont_hap1",
                "--seed",
                "43",
            ],
            [
                "pbsim",
                "--strategy",
                "templ",
                "--method",
                "qshmm",
                "--qshmm",
                "<dir>/QSHMM-ONT-HQ.model",
                "--template",
                "<dir>/template_hap2.fa",
                "--pass-num",
                "1",
                "--accuracy-mean",
                "0.95",
                "--prefix",
                "<dir>/ont_hap2",
                "--seed",
                "44",
            ],
        ],
    },
    "pacbio": {
        "templates": {
            "template_hap1.fa": "876285d754365e28ddb8090a393e9303",
            "template_hap2.fa": "0b8cc660c00c40eb9ba01df8d868c19f",
        },
        "calls": [
            [
                "pbsim",
                "--strategy",
                "templ",
                "--method",
                "errhmm",
                "--errhmm",
                "<dir>/ERRHMM-SEQUEL.model",
                "--template",
                "<dir>/template_hap1.fa",
                "--pass-num",
                "10",
                "--accuracy-mean",
                "0.85",
                "--prefix",
                "<dir>/clr_hap1",
                "--seed",
                "43",
            ],
            [
                "pbsim",
                "--strategy",
                "templ",
                "--method",
                "errhmm",
                "--errhmm",
                "<dir>/ERRHMM-SEQUEL.model",
                "--template",
                "<dir>/template_hap2.fa",
                "--pass-num",
                "10",
                "--accuracy-mean",
                "0.85",
                "--prefix",
                "<dir>/clr_hap2",
                "--seed",
                "44",
            ],
            [
                "ccs",
                "<dir>/clr_hap1_0001.bam",
                "<dir>/hifi_hap1_0001.bam",
                "--min-passes",
                "3",
                "--min-rq",
                "0.99",
                "--num-threads",
                "8",
            ],
            [
                "ccs",
                "<dir>/clr_hap2_0001.bam",
                "<dir>/hifi_hap2_0001.bam",
                "--min-passes",
                "3",
                "--min-rq",
                "0.99",
                "--num-threads",
                "8",
            ],
            [
                "samtools",
                "merge",
                "-f",
                "-@",
                "8",
                "<dir>/legacy_amplicon_hifi.bam",
                "<dir>/hifi_hap1_0001.bam",
                "<dir>/hifi_hap2_0001.bam",
            ],
            [
                "samtools",
                "fastq",
                "-@",
                "8",
                "-0",
                "<dir>/legacy_amplicon_hifi.fastq",
                "<dir>/legacy_amplicon_hifi.bam",
            ],
        ],
    },
}


def _norm(arg: str) -> str:
    """Replace temporary absolute paths by their file name."""
    return f"<dir>/{Path(arg).name}" if "/" in arg else arg


class _FakeTools:
    """Records tool calls and creates the files the wrappers expect."""

    def __init__(self) -> None:
        self.calls: list[list[str]] = []
        self.templates: dict[str, str] = {}

    def __call__(self, cmd: list[str], **_: Any) -> None:
        args = [str(a) for a in cmd]
        self.calls.append([_norm(a) for a in args])
        if "--template" in args:
            template = Path(args[args.index("--template") + 1])
            self.templates[template.name] = hashlib.md5(template.read_bytes()).hexdigest()
            prefix = args[args.index("--prefix") + 1]
            if args[args.index("--pass-num") + 1] == "1":
                with gzip.open(f"{prefix}.fq.gz", "wt") as handle:
                    handle.write("@r1\nACGT\n+\nIIII\n")
            else:
                Path(f"{prefix}_0001.bam").write_bytes(FAKE_BAM)
            return
        for arg in args[1:]:
            path = Path(arg)
            wanted = path.suffix in (".bam", ".fastq", ".fq")
            if wanted and not path.exists() and path.parent.exists():
                path.write_bytes(FAKE_BAM if path.suffix == ".bam" else b"@r\nA\n+\nI\n")


def _config(tmp_path: Path) -> Path:
    for name in ("QSHMM-ONT-HQ.model", "ERRHMM-SEQUEL.model"):
        (tmp_path / name).write_text("model")
    config = json.loads((Path(__file__).parents[2] / "config.json").read_text())
    config["tools"].update(pbsim3="pbsim", ccs="ccs", samtools="samtools")
    config["ont_amplicon_params"]["model_file"] = str(tmp_path / "QSHMM-ONT-HQ.model")
    config["pacbio_params"]["model_file"] = str(tmp_path / "ERRHMM-SEQUEL.model")
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return path


def _run_legacy(tmp_path: Path, platform: str, monkeypatch: pytest.MonkeyPatch) -> _FakeTools:
    fake = _FakeTools()
    for module in WRAPPER_MODULES:
        monkeypatch.setattr(f"{module}.run_command", fake)
    monkeypatch.setattr(
        "muc_one_up.read_simulator.utils.metadata_writer.capture_tool_versions",
        lambda tools: dict.fromkeys(tools, "mocked"),
    )
    fasta = tmp_path / "legacy.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in HAPLOTYPES.items()))
    result = CliRunner().invoke(
        cli,
        [
            "--config",
            str(_config(tmp_path)),
            "--log-level",
            "WARNING",
            "reads",
            "amplicon",
            str(fasta),
            "--platform",
            platform,
            "--coverage",
            str(COVERAGE),
            "--seed",
            str(SEED),
            "--no-align",
            "--out-dir",
            str(tmp_path / "out"),
            "--out-base",
            "legacy",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    return fake


@pytest.mark.parametrize("platform", ["ont", "pacbio"])
def test_seeded_legacy_amplicon_run_is_pinned(
    tmp_path: Path, platform: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = _run_legacy(tmp_path, platform, monkeypatch)
    assert fake.templates == EXPECTED[platform]["templates"]
    assert fake.calls == EXPECTED[platform]["calls"]
